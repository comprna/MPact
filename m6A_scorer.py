#!/usr/bin/env python3
"""Score known or candidate m6A sites with MPact.

This script evaluates MPact directly on m6A-centered genomic sites rather than
on SNV effects. Input rows are oriented to transcript sense using a strand
column: plus-strand sites are scored as genomic sequence and minus-strand sites
as reverse-complement sequence. Scoreable sites have an A at the center after
orientation.
"""

import argparse
import json
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import tensorflow as tf


SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from annotation_features import (  # noqa: E402
    ACCESSIBILITY_LENGTHS,
    AccessibilityCalculator,
    ConservationTrack,
)
from score_mpact import (  # noqa: E402
    FastaFetcher,
    ReduceSumAxis1,
    class_predictions,
    encode_with_position,
    normalize_chrom,
    reverse_complement,
)


SUPPORTED_WINDOW_SIZES = (101, 201, 501, 801, 1001)
THRESHOLDS = (0.5, 0.7, 0.8, 0.9, 0.95)


def first_present(columns, candidates):
    lookup = {str(c).lower(): c for c in columns}
    for candidate in candidates:
        if candidate in columns:
            return candidate
        found = lookup.get(candidate.lower())
        if found is not None:
            return found
    return None


def percentile(values, q):
    if not values:
        return None
    return float(np.percentile(np.asarray(values, dtype=float), q))


def empty_summary():
    d = {
        "n": 0,
        "score_sum": 0.0,
        "score_min": None,
        "score_max": None,
    }
    for threshold in THRESHOLDS:
        d[f"score_ge_{threshold:g}"] = 0
    return d


def update_score_summary(summary, scores):
    scores = np.asarray(scores, dtype=float)
    if scores.size == 0:
        return
    summary["n"] += int(scores.size)
    summary["score_sum"] += float(scores.sum())
    summary["score_min"] = (
        float(scores.min())
        if summary["score_min"] is None
        else min(summary["score_min"], float(scores.min()))
    )
    summary["score_max"] = (
        float(scores.max())
        if summary["score_max"] is None
        else max(summary["score_max"], float(scores.max()))
    )
    for threshold in THRESHOLDS:
        summary[f"score_ge_{threshold:g}"] += int((scores >= threshold).sum())


def finalize_summary(summary, score_sample):
    n = summary["n"]
    out = {"n": n}
    if n == 0:
        return out
    out.update(
        {
            "score_mean": summary["score_sum"] / n,
            "score_min": summary["score_min"],
            "score_q05": percentile(score_sample, 5),
            "score_q25": percentile(score_sample, 25),
            "score_median": percentile(score_sample, 50),
            "score_q75": percentile(score_sample, 75),
            "score_q95": percentile(score_sample, 95),
            "score_max": summary["score_max"],
        }
    )
    for threshold in THRESHOLDS:
        count = summary[f"score_ge_{threshold:g}"]
        out[f"n_score_ge_{threshold:g}"] = count
        out[f"frac_score_ge_{threshold:g}"] = count / n
    return out


def sample_extend(sample, values, max_items):
    if len(sample) >= max_items:
        return
    remaining = max_items - len(sample)
    sample.extend(np.asarray(values[:remaining], dtype=float).tolist())


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="Input m6A site TSV")
    p.add_argument("--output-tsv", required=True, help="Per-row scored output TSV")
    p.add_argument("--summary-json", required=True, help="Summary metrics JSON")
    p.add_argument("--summary-tsv", required=True, help="Summary metrics TSV")
    p.add_argument("--fasta", required=True, help="Reference FASTA")
    p.add_argument("--model-path", required=True, help="MPact Keras model")
    p.add_argument("--chrom-col", default=None, help="Chromosome column override")
    p.add_argument("--pos-col", default=None, help="1-based site position column override")
    p.add_argument("--end-col", default=None, help="Optional end column override")
    p.add_argument("--strand-col", default=None, help="Strand column override")
    p.add_argument("--group-col", default=None, help="Optional grouping column for summaries")
    p.add_argument(
        "--window-size",
        type=int,
        choices=SUPPORTED_WINDOW_SIZES,
        default=501,
        help="Input length expected by the selected model (default: 501)",
    )
    p.add_argument(
        "--conservation-bigwig",
        required=True,
        help="Required phyloP/phastCons bigWig used to add the site's genomic score",
    )
    p.add_argument(
        "--conservation-label",
        default=None,
        help="Output label for --conservation-bigwig (default: bigWig basename)",
    )
    p.add_argument("--chunksize", type=int, default=100000)
    p.add_argument("--batch-size", type=int, default=2048)
    p.add_argument("--max-percentile-sample", type=int, default=500000)
    p.add_argument("--max-rows", type=int, default=None, help="Optional smoke-test row limit")
    return p.parse_args()


def resolve_columns(args, columns):
    chrom_col = args.chrom_col or first_present(
        columns, ["chr", "chrom", "Chromosome", "#Chromosome", "CHROM"]
    )
    pos_col = args.pos_col or first_present(
        columns, ["start", "Start", "Position", "POS", "pos", "site", "Site"]
    )
    end_col = args.end_col or first_present(columns, ["end", "End"])
    strand_col = args.strand_col or first_present(columns, ["strand", "Strand"])

    missing = []
    if chrom_col is None:
        missing.append("chromosome")
    if pos_col is None:
        missing.append("position")
    if strand_col is None:
        missing.append("strand")
    if missing:
        raise SystemExit(
            "Missing required input columns: "
            + ", ".join(missing)
            + ". Use --chrom-col, --pos-col, and --strand-col to override."
        )

    group_col = args.group_col
    if group_col is None and "dataset" in columns:
        group_col = "dataset"
    return chrom_col, pos_col, end_col, strand_col, group_col


def main():
    args = parse_args()
    window_size = int(args.window_size)
    half = window_size // 2
    for path in [args.output_tsv, args.summary_json, args.summary_tsv]:
        out_dir = os.path.dirname(path)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
    accessibility_calculator = AccessibilityCalculator()
    conservation_track = ConservationTrack(args.conservation_bigwig, args.conservation_label)

    header = pd.read_csv(args.input, sep="\t", dtype=str, nrows=0)
    chrom_col, pos_col, end_col, strand_col, group_col = resolve_columns(args, header.columns)

    fetcher = FastaFetcher(args.fasta)
    model = tf.keras.models.load_model(
        args.model_path,
        custom_objects={
            "ReduceSumAxis1": ReduceSumAxis1,
            "_ReduceSumAxis1": ReduceSumAxis1,
            "MPact>ReduceSumAxis1": ReduceSumAxis1,
        },
        compile=False,
    )
    output_names = list(getattr(model, "output_names", []) or [])

    total_rows = 0
    invalid_rows = 0
    non_a_center_rows = 0
    model_window_size = int(model.input_shape[1])
    if model_window_size != window_size:
        raise ValueError(
            f"--window-size {window_size} does not match model input length "
            f"{model_window_size}: {args.model_path}"
        )
    scoreable_rows = 0
    all_summary = empty_summary()
    per_group = defaultdict(empty_summary)
    score_sample = []
    seen_sites = {}

    first_write = True
    for chunk in pd.read_csv(args.input, sep="\t", dtype=str, chunksize=args.chunksize):
        if args.max_rows is not None:
            remaining = int(args.max_rows) - total_rows
            if remaining <= 0:
                break
            chunk = chunk.head(remaining)
        total_rows += len(chunk)

        records = []
        seqs = []
        row_meta = []

        for rowd in chunk.to_dict("records"):
            rec = dict(rowd)
            chrom = normalize_chrom(rowd.get(chrom_col))
            strand = str(rowd.get(strand_col, "")).strip()
            try:
                pos1 = int(rowd.get(pos_col))
                end1 = int(rowd.get(end_col)) if end_col else pos1
            except Exception:
                invalid_rows += 1
                continue

            if pos1 != end1 or strand not in {"+", "-"}:
                invalid_rows += 1
                continue

            seq = fetcher.fetch(chrom, pos1 - half, pos1 + half)
            if len(seq) != window_size or "N" in seq:
                invalid_rows += 1
                continue

            oriented = seq if strand == "+" else reverse_complement(seq)
            center_base = oriented[half]
            scoreable = center_base == "A"
            if not scoreable:
                non_a_center_rows += 1

            rec.update(
                {
                    "mpact_chr": chrom,
                    "mpact_pos1": pos1,
                    "mpact_strand": strand,
                    "center_base_transcript": center_base,
                    "scoreable_center_A": scoreable,
                    "mpact_score": np.nan,
                    "accessibility_source": "ViennaRNA_RNAplfold_W200_L150",
                    "conservation_source": conservation_track.label,
                    "site_conservation_score": conservation_track.score(chrom, pos1),
                }
            )
            rec["mpact_ref_unpaired_probability_1nt"] = np.nan
            for length in ACCESSIBILITY_LENGTHS:
                rec[f"mpact_ref_accessibility_{length}nt"] = np.nan
            for suffix, value in accessibility_calculator.score(oriented).items():
                rec[f"mpact_ref_{suffix}"] = value
            records.append(rec)

            if scoreable:
                seqs.append(oriented)
                row_meta.append(len(records) - 1)

        if seqs:
            encoded = encode_with_position(seqs, window_size=window_size)
            pred = model.predict(encoded, batch_size=args.batch_size, verbose=0)
            scores = class_predictions(pred)
            scoreable_rows += len(scores)
            update_score_summary(all_summary, scores)
            sample_extend(score_sample, scores, args.max_percentile_sample)

            for local_i, score in zip(row_meta, scores):
                rec = records[local_i]
                score = float(score)
                rec["mpact_score"] = score

                group = rec.get(group_col, "all") if group_col else "all"
                update_score_summary(per_group[str(group)], [score])

                key = (rec["mpact_chr"], rec["mpact_pos1"], rec["mpact_strand"])
                previous = seen_sites.get(key)
                if previous is None or score > previous:
                    seen_sites[key] = score
        out = pd.DataFrame.from_records(records)
        out.to_csv(
            args.output_tsv,
            sep="\t",
            index=False,
            mode="w" if first_write else "a",
            header=first_write,
        )
        first_write = False
        print(
            f"Processed {total_rows} rows; scoreable_center_A={scoreable_rows}; "
            f"invalid={invalid_rows}; non_A_center={non_a_center_rows}",
            flush=True,
        )
    unique_summary = empty_summary()
    if seen_sites:
        unique_scores = np.asarray(list(seen_sites.values()), dtype=float)
        update_score_summary(unique_summary, unique_scores)
        unique_score_sample = unique_scores[: args.max_percentile_sample].tolist()
    else:
        unique_score_sample = []

    summary = {
        "input": args.input,
        "output_tsv": args.output_tsv,
        "model_path": args.model_path,
        "model_outputs": output_names,
        "window_size": window_size,
        "accessibility_source": "ViennaRNA_RNAplfold_W200_L150",
        "conservation_source": conservation_track.label,
        "columns": {
            "chrom": chrom_col,
            "position": pos_col,
            "end": end_col,
            "strand": strand_col,
            "group": group_col,
        },
        "total_input_rows": total_rows,
        "invalid_or_unfetchable_rows": invalid_rows,
        "non_A_center_rows": non_a_center_rows,
        "scoreable_center_A_rows": scoreable_rows,
        "unique_scoreable_sites": len(seen_sites),
        "all_scoreable_rows": finalize_summary(all_summary, score_sample),
        "unique_sites_best_score_per_chr_pos_strand": finalize_summary(
            unique_summary, unique_score_sample
        ),
        "groups": {
            group: finalize_summary(stats, [])
            for group, stats in sorted(per_group.items())
        },
    }

    with open(args.summary_json, "w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    rows = []
    for label, stats in [
        ("all_scoreable_rows", summary["all_scoreable_rows"]),
        ("unique_sites_best_score_per_chr_pos_strand", summary["unique_sites_best_score_per_chr_pos_strand"]),
    ]:
        row = {"group": label}
        row.update(stats)
        rows.append(row)
    pd.DataFrame(rows).to_csv(args.summary_tsv, sep="\t", index=False)

    print(f"Wrote: {args.output_tsv}")
    print(f"Wrote: {args.summary_json}")

    if conservation_track is not None:
        conservation_track.close()
    print(f"Wrote: {args.summary_tsv}")


if __name__ == "__main__":
    main()
