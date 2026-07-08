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

from score_mpact import (  # noqa: E402
    FastaFetcher,
    ReduceSumAxis1,
    encode_with_position,
    normalize_chrom,
    reverse_complement,
    split_predict_outputs,
)


HALF = 250
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
        "stoich_sum": 0.0,
        "stoich_min": None,
        "stoich_max": None,
    }
    for threshold in THRESHOLDS:
        d[f"score_ge_{threshold:g}"] = 0
    return d


def update_score_summary(summary, scores, stoich):
    scores = np.asarray(scores, dtype=float)
    stoich = np.asarray(stoich, dtype=float) if stoich is not None else None
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
    if stoich is not None:
        summary["stoich_sum"] += float(stoich.sum())
        summary["stoich_min"] = (
            float(stoich.min())
            if summary["stoich_min"] is None
            else min(summary["stoich_min"], float(stoich.min()))
        )
        summary["stoich_max"] = (
            float(stoich.max())
            if summary["stoich_max"] is None
            else max(summary["stoich_max"], float(stoich.max()))
        )


def finalize_summary(summary, score_sample, stoich_sample):
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
    if stoich_sample:
        out.update(
            {
                "stoich_mean_pct": summary["stoich_sum"] / n,
                "stoich_min_pct": summary["stoich_min"],
                "stoich_q05_pct": percentile(stoich_sample, 5),
                "stoich_median_pct": percentile(stoich_sample, 50),
                "stoich_q95_pct": percentile(stoich_sample, 95),
                "stoich_max_pct": summary["stoich_max"],
            }
        )
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
    for path in [args.output_tsv, args.summary_json, args.summary_tsv]:
        out_dir = os.path.dirname(path)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

    header = pd.read_csv(args.input, sep="\t", dtype=str, nrows=0)
    chrom_col, pos_col, end_col, strand_col, group_col = resolve_columns(args, header.columns)

    fetcher = FastaFetcher(args.fasta)
    model = tf.keras.models.load_model(
        args.model_path,
        custom_objects={"ReduceSumAxis1": ReduceSumAxis1, "_ReduceSumAxis1": ReduceSumAxis1},
        compile=False,
    )
    output_names = list(getattr(model, "output_names", []) or [])
    has_stoich_head = (len(getattr(model, "outputs", []) or []) > 1) or any(
        "stoich" in str(name).lower() for name in output_names
    )

    total_rows = 0
    invalid_rows = 0
    non_a_center_rows = 0
    scoreable_rows = 0
    all_summary = empty_summary()
    per_group = defaultdict(empty_summary)
    score_sample = []
    stoich_sample = []
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

            seq = fetcher.fetch(chrom, pos1 - HALF, pos1 + HALF)
            if len(seq) != 501 or "N" in seq:
                invalid_rows += 1
                continue

            oriented = seq if strand == "+" else reverse_complement(seq)
            center_base = oriented[HALF]
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
                    "mpact_stoich_source": "",
                    "mpact_stoichiometry_pct": np.nan,
                }
            )
            records.append(rec)

            if scoreable:
                seqs.append(oriented)
                row_meta.append(len(records) - 1)

        if seqs:
            encoded = encode_with_position(seqs)
            pred = model.predict(encoded, batch_size=args.batch_size, verbose=0)
            scores, stoich = split_predict_outputs(pred)
            if stoich is None:
                stoich = scores * 100.0
                stoich_source = "class_output_scaled_x100"
            else:
                stoich_source = "stoich_output"

            scoreable_rows += len(scores)
            update_score_summary(all_summary, scores, stoich)
            sample_extend(score_sample, scores, args.max_percentile_sample)
            sample_extend(stoich_sample, stoich, args.max_percentile_sample)

            for local_i, score, stoich_pct in zip(row_meta, scores, stoich):
                rec = records[local_i]
                score = float(score)
                stoich_pct = float(stoich_pct)
                rec["mpact_score"] = score
                rec["mpact_stoich_source"] = stoich_source
                rec["mpact_stoichiometry_pct"] = stoich_pct

                group = rec.get(group_col, "all") if group_col else "all"
                update_score_summary(per_group[str(group)], [score], [stoich_pct])

                key = (rec["mpact_chr"], rec["mpact_pos1"], rec["mpact_strand"])
                previous = seen_sites.get(key)
                if previous is None or score > previous[0]:
                    seen_sites[key] = (score, stoich_pct)

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
        unique_scores = np.asarray([v[0] for v in seen_sites.values()], dtype=float)
        unique_stoich = np.asarray([v[1] for v in seen_sites.values()], dtype=float)
        update_score_summary(unique_summary, unique_scores, unique_stoich)
        unique_score_sample = unique_scores[: args.max_percentile_sample].tolist()
        unique_stoich_sample = unique_stoich[: args.max_percentile_sample].tolist()
    else:
        unique_score_sample = []
        unique_stoich_sample = []

    summary = {
        "input": args.input,
        "output_tsv": args.output_tsv,
        "model_path": args.model_path,
        "model_outputs": output_names,
        "has_stoich_head": has_stoich_head,
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
        "all_scoreable_rows": finalize_summary(all_summary, score_sample, stoich_sample),
        "unique_sites_best_score_per_chr_pos_strand": finalize_summary(
            unique_summary, unique_score_sample, unique_stoich_sample
        ),
        "groups": {
            group: finalize_summary(stats, [], [])
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
    print(f"Wrote: {args.summary_tsv}")


if __name__ == "__main__":
    main()
