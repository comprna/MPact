#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

outdir="${TMPDIR:-/tmp}/mpact_m6a_scorer_smoke"
python_bin="${PYTHON:-python}"
mkdir -p "$outdir"

"$python_bin" m6A_scorer.py \
  --input mini_m6A_sites.tsv \
  --output-tsv "$outdir/mini_m6A_sites.scored.tsv" \
  --summary-json "$outdir/mini_m6A_sites.summary.json" \
  --summary-tsv "$outdir/mini_m6A_sites.summary.tsv" \
  --fasta hg38.fa \
  --model-path model_window_501.h5 \
  --batch-size 4

"$python_bin" - "$outdir/mini_m6A_sites.summary.json" <<'PY'
import json
import sys

with open(sys.argv[1]) as handle:
    summary = json.load(handle)

assert summary["total_input_rows"] == 4, summary
assert summary["scoreable_center_A_rows"] == 4, summary
assert summary["invalid_or_unfetchable_rows"] == 0, summary
assert summary["all_scoreable_rows"]["n"] == 4, summary
PY

echo "m6A_scorer smoke test passed: $outdir"
