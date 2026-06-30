#!/usr/bin/env bash
set -euo pipefail

# Example run for MPact_Scoring_Pipeline using bundled sample input.
# You can override defaults via environment variables, for example:
#   PYTHON_BIN=/path/to/python OUTPUT_DIR=/tmp/mpact_out bash example_run.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-${ROOT_DIR}/.venv/bin/python}"
INPUT_TSV="${INPUT_TSV:-${SCRIPT_DIR}/mini_scoreable.tsv}"
OUTPUT_DIR="${OUTPUT_DIR:-${SCRIPT_DIR}/example_output}"
FASTA="${FASTA:-${SCRIPT_DIR}/hg38.fa}"
GTF="${GTF:-${SCRIPT_DIR}/Homo_sapiens.GRCh38.110.gtf.gz}"
REDIPORTAL_GZ="${REDIPORTAL_GZ:-${SCRIPT_DIR}/TABLE1_hg38_v3.txt.gz}"
MODEL="${MODEL:-${SCRIPT_DIR}/model_window_501.h5}"
OUTPUT_TSV="${OUTPUT_DIR}/sample_predictions.tsv"
OUTPUT_PLOT="${OUTPUT_DIR}/sample_delta_hist.png"

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "ERROR: Python executable not found or not executable: ${PYTHON_BIN}" >&2
  exit 1
fi
if [[ ! -f "${INPUT_TSV}" ]]; then
  echo "ERROR: Input TSV not found: ${INPUT_TSV}" >&2
  exit 1
fi
if [[ ! -f "${FASTA}" ]]; then
  echo "ERROR: FASTA not found: ${FASTA}" >&2
  exit 1
fi
if [[ ! -f "${FASTA}.fai" ]]; then
  echo "ERROR: FASTA index missing: ${FASTA}.fai" >&2
  echo "Run: samtools faidx ${FASTA}" >&2
  exit 1
fi
if [[ ! -f "${MODEL}" ]]; then
  echo "ERROR: Model file not found: ${MODEL}" >&2
  exit 1
fi
if [[ ! -f "${GTF}" ]]; then
  echo "ERROR: GTF not found: ${GTF}" >&2
  exit 1
fi
if [[ ! -f "${REDIPORTAL_GZ}" ]]; then
  echo "ERROR: REDIportal A-to-I file not found: ${REDIPORTAL_GZ}" >&2
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"

echo "Running MPact sample scoring..."
"${PYTHON_BIN}" "${SCRIPT_DIR}/score_mpact.py" \
  --input "${INPUT_TSV}" \
  --output-tsv "${OUTPUT_TSV}" \
  --fasta "${FASTA}" \
  --model-path "${MODEL}" \
  --gtf "${GTF}" \
  --rediportal-gz "${REDIPORTAL_GZ}" \
  --output-plot "${OUTPUT_PLOT}" \
  --scan-radius 5 \
  --batch-size 256

echo
echo "Done. Outputs:"
echo "  ${OUTPUT_TSV}"
echo "  ${OUTPUT_PLOT}"
echo
echo "Preview:"
head -n 5 "${OUTPUT_TSV}" | column -t -s $'\t'
