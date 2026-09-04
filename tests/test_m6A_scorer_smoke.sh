#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

outdir="${TMPDIR:-/tmp}/mpact_m6a_scorer_smoke"
python_bin="${PYTHON:-python}"
conservation_bw="$outdir/test_conservation.bw"
mkdir -p "$outdir"

"$python_bin" - "$conservation_bw" <<'PY'
import pyBigWig
import sys

track = pyBigWig.open(sys.argv[1], "w")
track.addHeader([
    ("chr1", 248956422),
    ("chr3", 198295559),
    ("chr5", 181538259),
    ("chr21", 46709983),
])
track.addEntries(
    ["chr1", "chr3", "chr5", "chr21"],
    [42682718, 13321818, 141628575, 33517364],
    ends=[42682719, 13321819, 141628576, 33517365],
    values=[0.1, 0.2, 0.3, 0.4],
)
track.close()
PY

"$python_bin" m6A_scorer.py \
  --input mini_m6A_sites.tsv \
  --output-tsv "$outdir/mini_m6A_sites.scored.tsv" \
  --summary-json "$outdir/mini_m6A_sites.summary.json" \
  --summary-tsv "$outdir/mini_m6A_sites.summary.tsv" \
  --fasta hg38.fa \
  --model-path models/mpact_dtm6a_501nt_seed42.keras \
  --window-size 501 \
  --conservation-bigwig "$conservation_bw" \
  --batch-size 4

"$python_bin" - "$outdir/mini_m6A_sites.summary.json" "$outdir/mini_m6A_sites.scored.tsv" <<'PY'
import csv
import json
import sys

with open(sys.argv[1]) as handle:
    summary = json.load(handle)

assert summary["total_input_rows"] == 4, summary
assert summary["scoreable_center_A_rows"] == 4, summary
assert summary["invalid_or_unfetchable_rows"] == 0, summary
assert summary["all_scoreable_rows"]["n"] == 4, summary
assert summary["window_size"] == 501, summary

with open(sys.argv[2], newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    columns = set(reader.fieldnames)
    rows = list(reader)
required = {
    "mpact_score",
    "mpact_ref_unpaired_probability_1nt",
    "mpact_ref_accessibility_20nt",
    "site_conservation_score",
}
assert required <= columns, required - columns
assert not any("stoich" in column.lower() for column in columns)
assert all(row["site_conservation_score"] for row in rows)
PY

echo "m6A_scorer smoke test passed: $outdir"
