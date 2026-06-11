#!/usr/bin/env python3
"""
MPact m6A Scoring Pipeline - Strand-aware variant effect prediction.

Scores variants using the MPact model with strand-aware logic:
- Scans +/- N nt around each SNV for candidate A-centered windows
- Uses both reference and alternate allele sequence contexts
- Detects m6A disruptions, creations, and contextual changes
- Optional annotation against A-to-I editing sites (no filtering)

Input:
    VCF(.vcf/.vcf.gz) or TSV with columns: #Chromosome, Position, Reference, Alteration
  (optional: strand, transcript_id columns for strand inference)

Output:
  - Scored TSV: detailed per-candidate predictions with m6A deltas
  - Optional PNG histogram of delta distribution

Model:
  - window_501: 501nt sequence context with 6-channel positional encoding
  - Input: one-hot + sin/cos positional embeddings
  - Output: sigmoid probability of m6A modification
"""

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import gzip

import numpy as np
import pandas as pd
import tensorflow as tf

try:
    import pysam
except ImportError:
    pysam = None

# Ensure local modules can be imported from this directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Try to import AtoIFilter from local ato_i_filter.py
try:
    from ato_i_filter import AtoIFilter
except ImportError:
    AtoIFilter = None


HALF = 250


CODING_CONSEQUENCES = {
    "coding_sequence_variant",
    "frameshift_variant",
    "inframe_deletion",
    "inframe_insertion",
    "missense_variant",
    "protein_altering_variant",
    "start_lost",
    "start_retained_variant",
    "stop_gained",
    "stop_lost",
    "stop_retained_variant",
    "synonymous_variant",
    "incomplete_terminal_codon_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_region_variant",
}


class ReduceSumAxis1(tf.keras.layers.Layer):
    """Custom layer for proper serialization of reduce_sum operations."""
    def call(self, inputs):
        return tf.reduce_sum(inputs, axis=1)

    def compute_output_shape(self, input_shape):
        return (input_shape[0], input_shape[2])


_COMP = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(seq):
    return seq.translate(_COMP)[::-1]


def complement_base(base):
    return reverse_complement(str(base).upper())


def open_text_maybe_gzip(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


# ── Genome / annotation version helpers ──────────────────────────────────────

_GENOME_BUILD_ALIASES = {
    # Ensembl-style -> canonical label
    "grch38": "GRCh38", "hg38": "GRCh38",
    "grch37": "GRCh37", "hg19": "GRCh37", "hg37": "GRCh37",
    "mm10":   "mm10",   "grcm38": "mm10",
    "mm39":   "mm39",   "grcm39": "mm39",
}

def _canonical_build(token):
    return _GENOME_BUILD_ALIASES.get(token.lower())


def parse_genome_version(path):
    """Extract genome build and, for GTF, the annotation release from a filepath.
    Returns (build_label, release_str) where release_str is None for FASTA.
    Examples:
      Homo_sapiens.GRCh38.110.gtf.gz  -> ('GRCh38', '110')
      hg38.fa                          -> ('GRCh38', None)
    """
    import re as _re
    stem = os.path.basename(str(path)).lower()
    build = None
    release = None
    for token in _re.split(r'[._\-]+', stem):
        if build is None:
            b = _canonical_build(token)
            if b:
                build = b
        elif _re.fullmatch(r'\d+', token) and release is None:
            release = token  # first bare integer after build = release number
    return build, release


def check_fasta_gtf_version_consistency(fasta_path, gtf_path):
    """Warn if FASTA and GTF appear to use different genome builds."""
    if not fasta_path or not gtf_path:
        return
    fb, _ = parse_genome_version(fasta_path)
    gb, gr = parse_genome_version(gtf_path)
    if fb and gb and fb != gb:
        print(
            f"WARNING: Genome build mismatch detected!\n"
            f"  FASTA  -> {os.path.basename(fasta_path)}: build={fb}\n"
            f"  GTF    -> {os.path.basename(gtf_path)}: build={gb}\n"
            f"  Strand annotation may be incorrect if reference coordinates differ."
        )
    else:
        build_label = gb or fb or "unknown"
        release_label = f" release={gr}" if gr else ""
        print(f"Genome version check: build={build_label}{release_label}  [FASTA and GTF consistent]")


def gtf_version_label(gtf_path):
    """Return a short human-readable label for the GTF, e.g. 'GRCh38.110'."""
    build, release = parse_genome_version(gtf_path)
    parts = [p for p in (build, release) if p]
    return ".".join(parts) if parts else os.path.basename(str(gtf_path))


def parse_gtf_attributes(attr_text):
    attrs = {}
    for part in attr_text.strip().split(";"):
        part = part.strip()
        if not part or " " not in part:
            continue
        key, val = part.split(" ", 1)
        attrs[key] = val.strip().strip('"')
    return attrs


def load_transcript_strand_map(gtf_path):
    """Load transcript -> strand mapping from GTF file."""
    tx2strand = {}
    if not gtf_path:
        return tx2strand
    if not os.path.exists(gtf_path):
        print(f"WARNING: GTF not found: {gtf_path}")
        return tx2strand

    with open_text_maybe_gzip(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            strand = f[6]
            if strand not in {"+", "-"}:
                continue
            attrs = parse_gtf_attributes(f[8])
            tx = attrs.get("transcript_id")
            if not tx:
                continue
            tx2strand[tx] = strand
            tx2strand[tx.split(".")[0]] = strand

    print(f"Loaded transcript strands for {len(tx2strand)} transcript keys from {gtf_path}")
    return tx2strand


def extract_enst_from_sample(sample_value):
    if sample_value is None:
        return None
    m = re.search(r"(ENST\d+(?:\.\d+)?)", str(sample_value))
    if m:
        return m.group(1)
    return None


def normalize_strand_value(v):
    s = str(v).strip()
    return s if s in {"+", "-"} else None


class FastaFetcher:
    """FASTA accessor with fast in-process backend and samtools fallback."""

    def __init__(self, fasta_path):
        self.fasta_path = fasta_path
        self.samtools_path = shutil.which("samtools") or "/apps/samtools/1.22/bin/samtools"
        self.backend = "samtools"
        self._fa = None

        if pysam is not None:
            try:
                self._fa = pysam.FastaFile(fasta_path)
                self.backend = "pysam"
            except Exception:
                self._fa = None

    def _chrom_candidates(self, chrom):
        c = str(chrom)
        if c.startswith("chr"):
            return [c, c[3:]]
        return [c, "chr" + c]

    def fetch(self, chrom, start1, end1):
        """Fetch sequence for 1-based inclusive [start1, end1]."""
        if start1 < 1 or end1 < start1:
            return ""

        expected = end1 - start1 + 1

        if self.backend == "pysam" and self._fa is not None:
            for cand in self._chrom_candidates(chrom):
                try:
                    seq = self._fa.fetch(cand, start1 - 1, end1).upper()
                except Exception:
                    continue
                if len(seq) == expected:
                    return seq
            return ""

        region = f"{chrom}:{start1}-{end1}"
        r = subprocess.run(
            [self.samtools_path, "faidx", self.fasta_path, region],
            capture_output=True,
            text=True,
            check=False,
        )
        if r.returncode != 0:
            return ""
        lines = r.stdout.strip().splitlines()
        if len(lines) < 2:
            return ""
        seq = "".join(lines[1:]).upper()
        return seq if len(seq) == expected else ""


def encode_with_position(seqs, center_index=HALF):
    """Encode sequences with one-hot + positional embeddings (sin/cos)."""
    mapping = {
        "A": [1, 0, 0, 0],
        "C": [0, 1, 0, 0],
        "G": [0, 0, 1, 0],
        "T": [0, 0, 0, 1],
        "N": [0, 0, 0, 0],
    }
    n = len(seqs)
    X = np.zeros((n, 501, 6), dtype=np.float32)
    for i, seq in enumerate(seqs):
        for j, nt in enumerate(seq):
            X[i, j, :4] = mapping.get(nt, [0, 0, 0, 0])
            rel = j - center_index
            X[i, j, 4] = np.sin(rel / 10.0)
            X[i, j, 5] = np.cos(rel / 10.0)
    return X


def normalize_chrom(c):
    c = str(c).strip()
    if not c.startswith("chr"):
        return "chr" + c
    return c


def parse_vcf_info(info_text):
    info = {}
    for item in str(info_text).split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info


def parse_gene_symbols_from_geneinfo(geneinfo_text):
    """Extract gene symbols from ClinVar GENEINFO field."""
    symbols = []
    raw = str(geneinfo_text or "").strip()
    if not raw or raw in {".", "NA", "N/A"}:
        return symbols

    for entry in raw.split("|"):
        entry = entry.strip()
        if not entry:
            continue
        symbol = entry.split(":", 1)[0].strip()
        if symbol:
            symbols.append(symbol)
    return symbols


def load_gtf_strand_interval_index(gtf_path, include_features=None, bin_size=1_000_000):
    """Load a lightweight genomic interval index for strand inference.

    The index is keyed by chromosome and coarse bins; each entry stores
    (start1, end1, strand, gene_name_upper, gene_id_upper).
    """
    if include_features is None:
        include_features = {"gene"}

    if not gtf_path:
        return {}
    if not os.path.exists(gtf_path):
        print(f"WARNING: GTF not found for interval strand inference: {gtf_path}")
        return {}

    include_features = {str(v).strip().lower() for v in include_features}
    idx = {}
    n_added = 0

    with open_text_maybe_gzip(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue

            feature = f[2].strip().lower()
            if feature not in include_features:
                continue

            chrom = normalize_chrom(f[0].strip())
            strand = f[6].strip()
            if strand not in {"+", "-"}:
                continue

            try:
                start1 = int(f[3])
                end1 = int(f[4])
            except Exception:
                continue
            if end1 < start1:
                continue

            attrs = parse_gtf_attributes(f[8])
            gene_name = str(attrs.get("gene_name", "") or "").strip().upper()
            gene_id = str(attrs.get("gene_id", "") or "").strip().upper()

            b0 = max(0, start1 // bin_size)
            b1 = max(0, end1 // bin_size)
            chr_bins = idx.setdefault(chrom, {})
            rec = (start1, end1, strand, gene_name, gene_id)
            for b in range(b0, b1 + 1):
                chr_bins.setdefault(b, []).append(rec)
            n_added += 1

    print(
        f"Loaded GTF interval strand index: {n_added} features across {len(idx)} chromosomes "
        f"(bin_size={bin_size})"
    )
    return idx


def infer_strand_from_gtf_intervals(chrom, pos1, interval_index, gene_symbols=None, bin_size=1_000_000):
    """Infer strand at a genomic position from overlapping GTF intervals.

    If gene_symbols is provided, hits are constrained to matching gene_name/gene_id
    when possible. Returns '+', '-', or None if ambiguous/unresolved.
    """
    if not interval_index:
        return None

    chr_bins = interval_index.get(normalize_chrom(chrom))
    if not chr_bins:
        return None

    b = max(0, int(pos1) // bin_size)
    hits = []
    for rec in chr_bins.get(b, []):
        start1, end1, strand, gene_name, gene_id = rec
        if start1 <= pos1 <= end1:
            hits.append(rec)
    if not hits:
        return None

    if gene_symbols:
        gset = {str(x).strip().upper() for x in gene_symbols if str(x).strip()}
        if gset:
            constrained = [
                rec for rec in hits
                if (rec[3] and rec[3] in gset) or (rec[4] and rec[4] in gset)
            ]
            if constrained:
                hits = constrained

    strands = {rec[2] for rec in hits}
    if len(strands) == 1:
        return next(iter(strands))
    return None


def is_genic_from_clinvar_info(info_text):
    """Infer genic status from ClinVar-style INFO tags (MC/GENEINFO)."""
    info = parse_vcf_info(info_text)

    # Prefer molecular consequence if present.
    mc_val = str(info.get("MC", "") or "")
    if mc_val:
        tokens = []
        for item in mc_val.split(","):
            item = item.strip()
            if not item:
                continue
            # Common ClinVar form: SO:0001583|missense_variant
            if "|" in item:
                item = item.split("|", 1)[1]
            tokens.append(item.strip().lower())

        if "intron_variant" in tokens:
            return True
        if any("utr_variant" in t for t in tokens):
            return True
        if any(t in CODING_CONSEQUENCES for t in tokens):
            return True

    # Fallback: gene annotation present implies non-intergenic.
    geneinfo = str(info.get("GENEINFO", "") or "").strip()
    if geneinfo and geneinfo not in {".", "NA", "N/A"}:
        return True

    return False


def parse_vcf_annotation_fields(header_line, info_id):
    marker = f'ID={info_id}'
    if marker not in header_line:
        return None
    m = re.search(r'Format: ([^\"]+)', header_line)
    if not m:
        return None
    return [field.strip() for field in m.group(1).split("|")]


def first_annotation_value(annotation_fields, annotation_values, candidates):
    for candidate in candidates:
        if candidate not in annotation_fields:
            continue
        index = annotation_fields.index(candidate)
        if index >= len(annotation_values):
            continue
        value = annotation_values[index].strip()
        if value:
            return value
    return None


def extract_transcript_id_from_info(info_text, csq_fields=None, ann_fields=None):
    info = parse_vcf_info(info_text)

    if csq_fields and info.get("CSQ"):
        for entry in str(info["CSQ"]).split(","):
            values = entry.split("|")
            feature_type = first_annotation_value(csq_fields, values, ["Feature_type", "BIOTYPE"])
            transcript_id = first_annotation_value(
                csq_fields,
                values,
                ["Feature", "Transcript", "Transcript_ID", "transcript_id"],
            )
            if transcript_id and (feature_type in {None, "Transcript", "transcript", "mRNA"}):
                return transcript_id

    if ann_fields and info.get("ANN"):
        for entry in str(info["ANN"]).split(","):
            values = entry.split("|")
            feature_type = first_annotation_value(ann_fields, values, ["Feature_Type"])
            transcript_id = first_annotation_value(
                ann_fields,
                values,
                ["Feature_ID", "Transcript_ID", "transcript_id"],
            )
            if transcript_id and (feature_type in {None, "transcript", "Transcript", "mRNA"}):
                return transcript_id

    return None


def load_variants_input(input_path, genic_only=False):
    """Load variants from VCF(.gz) or TSV into a dataframe."""
    lp = str(input_path).lower()

    if lp.endswith(".vcf") or lp.endswith(".vcf.gz"):
        rows = []
        total_records = 0
        total_alt_alleles = 0
        kept_snv_alleles = 0
        csq_fields = None
        ann_fields = None
        with open_text_maybe_gzip(input_path) as fin:
            for line in fin:
                if not line:
                    continue
                if line.startswith("##INFO=<ID=CSQ"):
                    csq_fields = parse_vcf_annotation_fields(line, "CSQ")
                    continue
                if line.startswith("##INFO=<ID=ANN"):
                    ann_fields = parse_vcf_annotation_fields(line, "ANN")
                    continue
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    continue
                total_records += 1
                chrom, pos, vid, ref, alt_field, _, _, info_text = fields[:8]
                ref = str(ref).upper()
                alts = [a.strip().upper() for a in str(alt_field).split(",")]
                total_alt_alleles += len(alts)
                if len(ref) != 1 or ref not in {"A", "C", "G", "T"}:
                    continue
                info_map = parse_vcf_info(info_text)
                transcript_id = extract_transcript_id_from_info(info_text, csq_fields=csq_fields, ann_fields=ann_fields)
                is_genic = is_genic_from_clinvar_info(info_text)
                gene_symbols = parse_gene_symbols_from_geneinfo(info_map.get("GENEINFO", ""))
                for alt in alts:
                    if len(alt) != 1 or alt not in {"A", "C", "G", "T"}:
                        continue
                    kept_snv_alleles += 1
                    if genic_only and not is_genic:
                        continue
                    rows.append(
                        {
                            "#Chromosome": chrom,
                            "Position": pos,
                            "Reference": ref,
                            "Alteration": alt,
                            "VariantID": vid,
                            "transcript_id": transcript_id,
                            "gene_symbols": "|".join(gene_symbols),
                        }
                    )
        dropped_non_snv = total_alt_alleles - kept_snv_alleles
        n_before_genic = kept_snv_alleles
        if genic_only:
            dropped_non_genic = n_before_genic - len(rows)
            print(
                f"VCF genic filter: kept={len(rows)}, dropped_non_genic={dropped_non_genic}"
            )

        df = pd.DataFrame(rows)

        print(
            f"VCF parse summary: records={total_records}, alt_alleles={total_alt_alleles}, "
            f"kept_snv={kept_snv_alleles}, dropped_non_snv={dropped_non_snv}"
        )
        return df, list(df.columns), "vcf"

    df = pd.read_csv(input_path, sep="\t", dtype=str)
    return df, list(df.columns), "tsv"


def iter_input_chunks(input_path, genic_only=False, chunk_size=10000):
    """Yield input dataframe chunks for bounded-memory processing."""
    lp = str(input_path).lower()
    if lp.endswith(".vcf") or lp.endswith(".vcf.gz"):
        df, input_cols, input_kind = load_variants_input(input_path, genic_only=genic_only)
        yield df, input_cols, input_kind
        return

    for chunk in pd.read_csv(input_path, sep="\t", dtype=str, chunksize=int(chunk_size)):
        yield chunk, list(chunk.columns), "tsv"


def resolve_required_columns(df):
    """Resolve key variant columns from common aliases and standardize names."""

    def first_present(cands):
        for c in cands:
            if c in df.columns:
                return c
        return None

    chrom_col = first_present(["#Chromosome", "Chromosome", "CHROM", "chrom", "chr"])
    pos_col = first_present(["Position", "POS", "pos", "Start", "start"])
    ref_col = first_present(["Reference", "REF", "ref"])
    alt_col = first_present(["Alteration", "ALT", "alt", "Alternate"])

    missing = []
    if chrom_col is None:
        missing.append("#Chromosome/CHROM")
    if pos_col is None:
        missing.append("Position/POS")
    if ref_col is None:
        missing.append("Reference/REF")
    if alt_col is None:
        missing.append("Alteration/ALT")
    if missing:
        raise ValueError(
            "Input is missing required variant columns: "
            + ", ".join(missing)
            + f". Available columns: {list(df.columns)}"
        )

    if "#Chromosome" not in df.columns:
        df["#Chromosome"] = df[chrom_col]
    if "Position" not in df.columns:
        df["Position"] = df[pos_col]
    if "Reference" not in df.columns:
        df["Reference"] = df[ref_col]
    if "Alteration" not in df.columns:
        df["Alteration"] = df[alt_col]

    return df


def flatten_predict_output(pred):
    """Handle Keras predict outputs across model output styles (ndarray/dict/list)."""
    if isinstance(pred, dict):
        pred = next(iter(pred.values()))
    elif isinstance(pred, (list, tuple)):
        pred = pred[0]
    return np.asarray(pred).reshape(-1)


def add_delta_stats(df):
    """Add z-score and p-value columns to delta scores."""
    delta = df["mpact_delta_alt_minus_ref"].astype(float).to_numpy()
    mu = float(delta.mean())
    sd = float(delta.std(ddof=0))
    if sd == 0.0:
        z = np.zeros_like(delta)
    else:
        z = (delta - mu) / sd

    root2 = math.sqrt(2.0)
    pvals = np.array([math.erfc(abs(float(v)) / root2) for v in z], dtype=float)

    df["delta_zscore"] = z
    df["delta_p_two_sided"] = pvals
    return mu, sd


def compute_delta_stats_from_tsv(tsv_path, chunk_size=250000):
    """Compute mean/std for mpact deltas from a large TSV without loading it all."""
    n = 0
    s = 0.0
    ss = 0.0
    for part in pd.read_csv(
        tsv_path,
        sep="\t",
        usecols=["mpact_delta_alt_minus_ref"],
        chunksize=int(chunk_size),
    ):
        arr = part["mpact_delta_alt_minus_ref"].astype(float).to_numpy()
        if arr.size == 0:
            continue
        n += int(arr.size)
        s += float(arr.sum())
        ss += float((arr * arr).sum())

    if n == 0:
        return 0.0, 0.0

    mu = s / n
    var = max(0.0, (ss / n) - (mu * mu))
    sd = math.sqrt(var)
    return mu, sd


def predict_scores_in_batches(model, seqs, batch_size):
    """Predict MPact scores for sequence list in bounded-memory batches."""
    n = len(seqs)
    out = np.zeros(n, dtype=np.float32)
    if n == 0:
        return out

    for i in range(0, n, batch_size):
        j = min(i + batch_size, n)
        X = encode_with_position(seqs[i:j])
        out[i:j] = flatten_predict_output(model.predict(X, batch_size=batch_size, verbose=0))
    return out


def z_and_p_from_delta(delta, mu, sd):
    """Compute z-scores and two-sided p-values from delta vector and global stats."""
    delta = np.asarray(delta, dtype=float)
    if sd == 0.0:
        z = np.zeros_like(delta)
    else:
        z = (delta - mu) / sd
    root2 = math.sqrt(2.0)
    pvals = np.array([math.erfc(abs(float(v)) / root2) for v in z], dtype=float)
    return z, pvals


def plot_delta_hist(delta_values, out_png):
    """Generate histogram of m6A delta distribution."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        if isinstance(delta_values, pd.DataFrame):
            delta = delta_values["mpact_delta_alt_minus_ref"].astype(float).to_numpy()
        else:
            delta = np.asarray(delta_values, dtype=float)

        if delta.size == 0:
            print("WARNING: no delta values available for histogram")
            return

        plt.figure(figsize=(8, 5))
        plt.hist(delta, bins=80, color="#2f6db3", alpha=0.9)
        plt.axvline(0.0, color="black", linestyle="--", linewidth=1)
        plt.xlabel("m6A-delta (ALT - REF)")
        plt.ylabel("Count")
        plt.title("MPact m6A score delta distribution")
        plt.tight_layout()
        plt.savefig(out_png, dpi=200)
        plt.close()
        print(f"Delta histogram saved: {out_png}")
    except Exception as e:
        print(f"WARNING: could not render plot: {e}")


def main():
    p = argparse.ArgumentParser(
        description="Score variants with MPact m6A model using strand-aware scanning"
    )
    p.add_argument(
        "--input",
        required=True,
        help="Input variants file (.tsv, .vcf, or .vcf.gz)"
    )
    p.add_argument(
        "--output-tsv",
        required=True,
        help="Output TSV with MPact predictions"
    )
    p.add_argument(
        "--fasta",
        required=True,
        help="Reference FASTA (indexed with samtools faidx)"
    )
    p.add_argument(
        "--model-path",
        required=True,
        help="Path to MPact model (HDF5 format)"
    )
    p.add_argument(
        "--output-plot",
        default="",
        help="Optional: output PNG histogram of delta scores"
    )
    p.add_argument(
        "--gtf",
        default="Homo_sapiens.GRCh38.110.gtf.gz",
        required=False,
        help="GTF/GTF.gz for Gencode strand inference (REQUIRED for strand certainty; default: bundled Ensembl GRCh38.110)"
    )
    p.add_argument(
        "--gtf-interval-features",
        default="gene",
        help="Comma-separated GTF feature types to index for interval strand fallback (default: gene; genic-locus based)",
    )
    gtf_geneinfo_group = p.add_mutually_exclusive_group()
    gtf_geneinfo_group.add_argument(
        "--gtf-interval-use-geneinfo",
        dest="gtf_interval_use_geneinfo",
        action="store_true",
        default=True,
        help="When falling back to GTF interval strand inference, constrain overlaps by GENEINFO symbols when present (default: on)",
    )
    gtf_geneinfo_group.add_argument(
        "--no-gtf-interval-use-geneinfo",
        dest="gtf_interval_use_geneinfo",
        action="store_false",
        help="Disable GENEINFO-constrained interval strand inference",
    )
    p.add_argument(
        "--rediportal-gz",
        default="TABLE1_hg38_v3.txt.gz",
        help="REDIportal A-to-I sites file for annotation (required; bundled local symlink by default)"
    )
    p.add_argument(
        "--scan-radius",
        type=int,
        default=5,
        help="Scanning radius around SNV for candidate A sites (default: 5 nt)"
    )
    p.add_argument(
        "--batch-size",
        type=int,
        default=1024,
        help="Batch size for model prediction (default: 1024 in optimized version)"
    )
    p.add_argument(
        "--input-chunk-size",
        type=int,
        default=50000,
        help="Rows per input chunk for streaming mode (default: 50000 in optimized version)"
    )
    p.add_argument(
        "--resume",
        action="store_true",
        help="Resume from a previous interrupted run using checkpoint + temporary scored TSV"
    )
    p.add_argument(
        "--resume-from-row",
        type=int,
        default=0,
        help="0-based input row index to resume from; rows before this are skipped"
    )
    p.add_argument(
        "--checkpoint-path",
        default="",
        help="Optional checkpoint JSON path (default: <output-tsv>.checkpoint.json)"
    )
    p.add_argument(
        "--keep-raw-temp",
        action="store_true",
        help="Keep intermediate <output-tsv>.raw_scored.tsv after successful completion"
    )
    genic_group = p.add_mutually_exclusive_group()
    genic_group.add_argument(
        "--genic-only",
        dest="genic_only",
        action="store_true",
        default=True,
        help="For VCF input, keep only genic variants inferred from INFO tags (MC/GENEINFO) (default: on)",
    )
    genic_group.add_argument(
        "--allow-nongenic",
        dest="genic_only",
        action="store_false",
        help="For VCF input, disable genic-only filtering and keep intergenic variants too",
    )
    args = p.parse_args()

    out_dir = os.path.dirname(args.output_tsv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # ── Version consistency check and GTF validation ──────────────────────────
    if not args.gtf or not os.path.exists(args.gtf):
        raise FileNotFoundError(
            f"GTF file is REQUIRED for strand certainty (provided or default not found): {args.gtf}"
        )
    check_fasta_gtf_version_consistency(args.fasta, args.gtf)
    _gtf_label = gtf_version_label(args.gtf)
    print(f"Gencode strand source: {args.gtf}  (label: {_gtf_label})")

    # Load required A-to-I annotation index.
    if not args.rediportal_gz or not os.path.exists(args.rediportal_gz):
        raise FileNotFoundError(
            f"REDIportal A-to-I file is REQUIRED for scoring and was not found: {args.rediportal_gz}"
        )
    if AtoIFilter is None:
        raise RuntimeError("ato_i_filter module is required for A-to-I annotation but is not available")
    print(f"Loading A-to-I sites from {args.rediportal_gz}...")
    ato_i = AtoIFilter(args.rediportal_gz)
    print("A-to-I sites indexed")

    # Load model once and score in bounded-memory chunks.
    print(f"Loading model: {args.model_path}")
    model = tf.keras.models.load_model(
        args.model_path,
        compile=False,
        custom_objects={"_ReduceSumAxis1": ReduceSumAxis1, "ReduceSumAxis1": ReduceSumAxis1},
    )
    print("Model loaded.")

    fasta_fetcher = FastaFetcher(args.fasta)
    print(f"FASTA backend: {fasta_fetcher.backend}")

    tx_strand_map = load_transcript_strand_map(args.gtf)
    gtf_interval_bin_size = 1_000_000
    features = {
        v.strip().lower()
        for v in str(args.gtf_interval_features or "").split(",")
        if v.strip()
    }
    if not features:
        features = {"gene"}
    gtf_interval_index = load_gtf_strand_interval_index(
        args.gtf,
        include_features=features,
        bin_size=gtf_interval_bin_size,
    )

    checkpoint_path = args.checkpoint_path or f"{args.output_tsv}.checkpoint.json"
    raw_tmp_path = f"{args.output_tsv}.raw_scored.tsv"
    resume_from_row = int(args.resume_from_row or 0)

    if args.resume:
        if resume_from_row <= 0 and os.path.exists(checkpoint_path):
            try:
                with open(checkpoint_path, "r") as f:
                    ckpt = json.load(f)
                resume_from_row = int(ckpt.get("last_input_row", -1)) + 1
                print(f"Resuming from checkpoint row: {resume_from_row}")
            except Exception as e:
                print(f"WARNING: failed to read checkpoint {checkpoint_path}: {e}")
        if not os.path.exists(raw_tmp_path):
            if resume_from_row > 0:
                print(
                    f"Resume mode with explicit start row ({resume_from_row}) and no raw temp TSV; "
                    "starting a fresh temp output from that row"
                )
            else:
                print("WARNING: --resume set but raw temp TSV is missing; starting a fresh temp output")
    else:
        if os.path.exists(raw_tmp_path):
            os.remove(raw_tmp_path)
        if os.path.exists(checkpoint_path):
            os.remove(checkpoint_path)

    # strand_user_provided is already in input cols if user provided a strand column.
    # strand_gencode is always the Gencode-inferred strand used for scoring.
    # strand_gencode_source records the exact GTF version used.
    # strand_gencode_inference records how strand was inferred.
    score_cols = [
        "strand_gencode",
        "strand_gencode_source",
        "strand_gencode_inference",
        "score_type",
        "a_genomic_pos1",
        "snp_to_a_mRNA_offset",
        "ref_center_is_A",
        "alt_center_is_A",
        "overlaps_AtoI_exact",
        "near_AtoI_5nt",
        "near_AtoI_10nt",
        "mpact_ref_score",
        "mpact_alt_score",
        "mpact_ref_stoichiometry_pct",
        "mpact_alt_stoichiometry_pct",
        "mpact_delta_stoichiometry_pct",
        "alt_center_A_destroyed",
        "mpact_delta_alt_minus_ref",
        "mpact_abs_delta",
        "delta_zscore",
        "delta_p_two_sided",
        "ref10",
        "alt10",
    ]

    raw_cols = [c for c in score_cols if c not in {"delta_zscore", "delta_p_two_sided"}]

    effective_input_cols = None
    strand_col_name = None
    transcript_col_name = None
    scan_radius = int(args.scan_radius)
    total_input_rows_seen = 0
    total_input_rows_used = 0
    total_candidates = 0

    write_mode = "a" if (args.resume and os.path.exists(raw_tmp_path) and os.path.getsize(raw_tmp_path) > 0) else "w"
    wrote_header = write_mode == "a"

    for chunk_idx, (df_in, input_cols, input_kind) in enumerate(
        iter_input_chunks(args.input, genic_only=args.genic_only, chunk_size=args.input_chunk_size),
        start=1,
    ):
        chunk_n = len(df_in)
        if chunk_n == 0:
            continue

        chunk_start = total_input_rows_seen
        chunk_end = total_input_rows_seen + chunk_n
        total_input_rows_seen = chunk_end

        if chunk_end <= resume_from_row:
            continue

        if chunk_start < resume_from_row:
            skip = resume_from_row - chunk_start
            df_in = df_in.iloc[skip:].copy()
            chunk_start = resume_from_row

        if effective_input_cols is None:
            strand_col_name = next((c for c in ["strand", "Strand"] if c in input_cols), None)
            transcript_col_name = next((c for c in ["transcript_id", "Transcript", "ENST", "Sample"] if c in input_cols), None)
            effective_input_cols = [
                "strand_user_provided" if c == strand_col_name else c for c in input_cols
            ]

            if input_kind == "vcf":
                n_with_tx = 0
                if transcript_col_name is not None and transcript_col_name in df_in.columns:
                    n_with_tx = int(df_in[transcript_col_name].fillna("").astype(str).str.strip().ne("").sum())
                if n_with_tx == 0:
                    print(
                        "WARNING: VCF input has no transcript annotations usable for strand inference. "
                        "Strand will fall back to genome orientation unless transcript IDs can be inferred."
                    )

            print(
                f"Initialized streaming input: kind={input_kind}, chunk_size={args.input_chunk_size}, "
                f"resume_from_row={resume_from_row}"
            )

        df_in = resolve_required_columns(df_in.copy())

        # Normalize and validate key fields
        df_in["_chrom"] = df_in["#Chromosome"].map(normalize_chrom)
        df_in["_pos1"] = pd.to_numeric(df_in["Position"], errors="coerce").astype("Int64")
        df_in["_ref"] = df_in["Reference"].astype(str).str.upper()
        df_in["_alt"] = df_in["Alteration"].astype(str).str.upper()
        if "gene_symbols" not in df_in.columns:
            df_in["gene_symbols"] = ""

        n_before_filter = len(df_in)
        valid = (
            df_in["_pos1"].notna()
            & (df_in["_ref"].str.len() == 1)
            & (df_in["_alt"].str.len() == 1)
            & df_in["_ref"].str.match(r"^[ACGT]$")
            & df_in["_alt"].str.match(r"^[ACGT]$")
        )
        df_in = df_in.loc[valid].copy()
        df_in["_pos1"] = df_in["_pos1"].astype(int)
        n_dropped_filter = n_before_filter - len(df_in)
        if n_dropped_filter:
            print(f"Chunk {chunk_idx}: dropped {n_dropped_filter} non-SNV/invalid rows")

        # Preserve user strand as metadata only; scoring always uses Gencode strand inference.
        if strand_col_name is not None and strand_col_name in df_in.columns:
            df_in["strand_user_provided"] = df_in[strand_col_name].map(normalize_strand_value)
            df_in = df_in.drop(columns=[strand_col_name])

        df_in["_strand"] = None
        df_in["_strand_inference"] = None
        if tx_strand_map and transcript_col_name is not None and transcript_col_name in df_in.columns:
            def _strand_from_tx(v):
                tx = extract_enst_from_sample(v) if transcript_col_name == "Sample" else (None if v is None else str(v))
                if not tx:
                    return None
                return tx_strand_map.get(tx) or tx_strand_map.get(tx.split(".")[0])

            df_in["_strand"] = df_in[transcript_col_name].map(_strand_from_tx)
            df_in.loc[df_in["_strand"].notna(), "_strand_inference"] = "transcript"

        missing_after_tx = df_in["_strand"].isna()
        if bool(missing_after_tx.any()):
            sub = df_in.loc[missing_after_tx, ["_chrom", "_pos1", "gene_symbols"]].copy()
            gene_symbol_cache = {}
            inferred = []
            inferred_mode = []
            for chrom, pos1, gene_symbols_raw in sub.itertuples(index=False, name=None):
                genes = None
                if args.gtf_interval_use_geneinfo:
                    key = str(gene_symbols_raw or "")
                    genes = gene_symbol_cache.get(key)
                    if genes is None:
                        genes = parse_gene_symbols_from_geneinfo(key)
                        gene_symbol_cache[key] = genes

                # First pass: constrained by GENEINFO if enabled.
                strand = infer_strand_from_gtf_intervals(
                    chrom,
                    int(pos1),
                    gtf_interval_index,
                    gene_symbols=genes,
                    bin_size=gtf_interval_bin_size,
                )
                mode = "gtf_geneinfo" if (strand is not None and args.gtf_interval_use_geneinfo) else None

                # Second pass fallback: relaxed interval lookup only when constrained pass fails.
                if strand is None and args.gtf_interval_use_geneinfo:
                    strand = infer_strand_from_gtf_intervals(
                        chrom,
                        int(pos1),
                        gtf_interval_index,
                        gene_symbols=None,
                        bin_size=gtf_interval_bin_size,
                    )
                    if strand is not None:
                        mode = "gtf_relaxed"

                if strand is not None and mode is None:
                    mode = "gtf_relaxed"

                inferred.append(strand)
                inferred_mode.append(mode)
            df_in.loc[missing_after_tx, "_strand"] = inferred
            df_in.loc[missing_after_tx, "_strand_inference"] = inferred_mode

        records = []
        wide_half = scan_radius + HALF
        expected_len = 2 * wide_half + 1

        for _, row in df_in.iterrows():
            chrom = row["_chrom"]
            snp_pos1 = int(row["_pos1"])
            ref = row["_ref"]
            alt = row["_alt"]
            base = {c: row.get(c, "") for c in effective_input_cols}

            # Fetch once per SNV instead of per-candidate
            wide_seq_plus = fasta_fetcher.fetch(chrom, snp_pos1 - wide_half, snp_pos1 + wide_half)
            if len(wide_seq_plus) != expected_len:
                continue

            strand = normalize_strand_value(row.get("_strand", None))
            if strand is None:
                print(f"WARNING: Could not infer strand for {chrom}:{snp_pos1} {ref}->{alt}; skipping row.")
                continue

            ref_oriented = ref if strand == "+" else complement_base(ref)
            alt_oriented = alt if strand == "+" else complement_base(alt)

            # Generate all candidate positions at once
            for delta in range(-scan_radius, scan_radius + 1):
                a_pos1 = snp_pos1 + delta
                wide_idx = wide_half + delta
                sub_start = wide_idx - HALF
                sub_end = wide_idx + HALF + 1
                if sub_start < 0 or sub_end > len(wide_seq_plus):
                    continue

                seq_plus_501 = wide_seq_plus[sub_start:sub_end]
                oriented = seq_plus_501 if strand == "+" else reverse_complement(seq_plus_501)

                snp_idx = HALF + (snp_pos1 - a_pos1) if strand == "+" else HALF + (a_pos1 - snp_pos1)
                snp_to_a = snp_pos1 - a_pos1 if strand == "+" else a_pos1 - snp_pos1

                if not (0 <= snp_idx < 501):
                    continue
                if oriented[snp_idx] != ref_oriented:
                    continue

                alt_seq = oriented[:snp_idx] + alt_oriented + oriented[snp_idx + 1:]
                ref_center_is_A = oriented[HALF] == "A"
                alt_center_is_A = alt_seq[HALF] == "A"
                if not (ref_center_is_A or alt_center_is_A):
                    continue

                score_type = ("disruption" if ref_center_is_A and not alt_center_is_A
                             else "creation" if (not ref_center_is_A) and alt_center_is_A
                             else "context")

                overlaps_ato_i_exact = False
                near_ato_i_5nt = False
                near_ato_i_10nt = False
                if ato_i is not None:
                    overlaps_ato_i_exact = ato_i.overlaps_exact(chrom, a_pos1, strand)
                    near_dist = ato_i.nearest_distance(chrom, a_pos1, strand)
                    if near_dist is not None:
                        near_ato_i_5nt = near_dist <= 5
                        near_ato_i_10nt = near_dist <= 10

                records.append(
                    {
                        **base,
                        "strand_gencode": strand,
                        "strand_gencode_source": _gtf_label,
                        "strand_gencode_inference": row.get("_strand_inference", ""),
                        "score_type": score_type,
                        "a_genomic_pos1": a_pos1,
                        "snp_to_a_mRNA_offset": snp_to_a,
                        "ref501": oriented,
                        "alt501": alt_seq,
                        "ref_center_is_A": ref_center_is_A,
                        "alt_center_is_A": alt_center_is_A,
                        "overlaps_AtoI_exact": overlaps_ato_i_exact,
                        "near_AtoI_5nt": near_ato_i_5nt,
                        "near_AtoI_10nt": near_ato_i_10nt,
                    }
                )

        total_input_rows_used += len(df_in)

        if not records:
            with open(checkpoint_path, "w") as f:
                json.dump({"last_input_row": chunk_end - 1, "candidates_written": total_candidates}, f)
            print(f"Chunk {chunk_idx}: no scoreable candidates")
            continue

        df = pd.DataFrame(records)
        ref_scores = predict_scores_in_batches(model, df["ref501"].tolist(), args.batch_size)
        alt_center_is_a = (df["alt501"].str[HALF] == "A").to_numpy()
        alt_scores = np.zeros(len(df), dtype=np.float32)
        if alt_center_is_a.any():
            alt_scores_sub = predict_scores_in_batches(
                model,
                df.loc[alt_center_is_a, "alt501"].tolist(),
                args.batch_size,
            )
            alt_scores[alt_center_is_a] = alt_scores_sub

        df["mpact_ref_score"] = ref_scores
        df["mpact_alt_score"] = alt_scores
        df["alt_center_A_destroyed"] = ~alt_center_is_a
        df["mpact_delta_alt_minus_ref"] = df["mpact_alt_score"] - df["mpact_ref_score"]
        df["mpact_abs_delta"] = df["mpact_delta_alt_minus_ref"].abs()
        df["mpact_ref_stoichiometry_pct"] = df["mpact_ref_score"] * 100.0
        df["mpact_alt_stoichiometry_pct"] = df["mpact_alt_score"] * 100.0
        df["mpact_delta_stoichiometry_pct"] = df["mpact_delta_alt_minus_ref"] * 100.0
        df["ref10"] = df["ref501"].str[HALF - 5 : HALF + 5]
        df["alt10"] = df["alt501"].str[HALF - 5 : HALF + 5]

        out_part = df[effective_input_cols + raw_cols].copy()
        out_part.to_csv(raw_tmp_path, sep="\t", index=False, mode=write_mode, header=(not wrote_header))
        write_mode = "a"
        wrote_header = True

        total_candidates += len(out_part)
        with open(checkpoint_path, "w") as f:
            json.dump({"last_input_row": chunk_end - 1, "candidates_written": total_candidates}, f)

        print(
            f"Chunk {chunk_idx}: input_rows={len(df_in)}, candidates={len(out_part)}, "
            f"total_candidates={total_candidates}"
        )

    if not os.path.exists(raw_tmp_path) or os.path.getsize(raw_tmp_path) == 0:
        cols = (effective_input_cols or []) + ["strand_gencode", "strand_gencode_source"]
        pd.DataFrame(columns=cols).to_csv(args.output_tsv, sep="\t", index=False)
        if os.path.exists(checkpoint_path):
            os.remove(checkpoint_path)
        print(f"No scoreable candidates. Wrote empty: {args.output_tsv}")
        return

    print(
        f"Scoring pass complete: processed_input_rows={total_input_rows_used}, "
        f"candidate_records={total_candidates}"
    )

    mu, sd = compute_delta_stats_from_tsv(raw_tmp_path)
    print(f"Global delta stats from temp output: mean={mu:.6g}, std={sd:.6g}")

    out_cols = (effective_input_cols or []) + score_cols
    plot_sample = []
    first_chunk = True
    for part in pd.read_csv(raw_tmp_path, sep="\t", chunksize=250000):
        delta = part["mpact_delta_alt_minus_ref"].astype(float).to_numpy()
        z, pvals = z_and_p_from_delta(delta, mu, sd)
        part["delta_zscore"] = z
        part["delta_p_two_sided"] = pvals

        if args.output_plot and len(plot_sample) < 200000:
            room = 200000 - len(plot_sample)
            plot_sample.extend(delta[:room].tolist())

        part[out_cols].to_csv(
            args.output_tsv,
            sep="\t",
            index=False,
            mode="w" if first_chunk else "a",
            header=first_chunk,
        )
        first_chunk = False

    if args.output_plot:
        plot_delta_hist(np.asarray(plot_sample, dtype=float), args.output_plot)

    if not args.keep_raw_temp and os.path.exists(raw_tmp_path):
        os.remove(raw_tmp_path)
    if os.path.exists(checkpoint_path):
        os.remove(checkpoint_path)

    print(f"\nWrote {total_candidates} scored candidates to: {args.output_tsv}")
    print(f"Delta mean={mu:.6g}, std={sd:.6g}")


if __name__ == "__main__":
    main()
