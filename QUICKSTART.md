# MPact Scoring Pipeline - Quick Start

This file is the minimal runbook. Full details are in README.md.

## 1. Setup

```bash
git clone https://github.com/comprna/MPact.git
cd MPact

# First time only
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## 2. Bundled paths (all local to this directory)

- Script: score_mpact.py
- Default model: model_window_501.h5 (501 nt window)
- Alternative trained models: model_window_101.h5 and model_window_201.h5 (use with --window-size 101 or --window-size 201)
- FASTA: hg38.fa
- GTF: Homo_sapiens.GRCh38.110.gtf.gz

## 3. Run on bundled sample (recommended smoke test)

```bash
cd MPact
python score_mpact.py \
  --input mini_scoreable.tsv \
  --output-tsv sample_predictions.tsv \
  --fasta hg38.fa \
  --model-path model_window_501.h5 \
  --output-plot sample_delta_hist.png
```

## 4. Run on your own file

### TSV input

```bash
cd MPact
python score_mpact.py \
  --input /path/to/variants.tsv \
  --output-tsv /path/to/predictions.tsv \
  --fasta hg38.fa \
  --model-path model_window_501.h5 \
  --output-plot /path/to/delta_histogram.png
```

### VCF input

```bash
cd MPact
python score_mpact.py \
  --input /path/to/variants.vcf.gz \
  --output-tsv /path/to/predictions.tsv \
  --fasta hg38.fa \
  --model-path model_window_501.h5 \
  --gtf Homo_sapiens.GRCh38.110.gtf.gz
```

## 5. Required and optional flags

**Required:**
- `--gtf Homo_sapiens.GRCh38.110.gtf.gz` (or custom) — Strand inference from GTF is mandatory for correctness. Bundled default is recommended.

**Optional:**
- Keep intergenic VCF variants: `--allow-nongenic`
- Smaller memory footprint: `--batch-size 256`
- Larger neighborhood scan: `--scan-radius 20`
- Add A-to-I annotation: `--rediportal-gz /path/to/TABLE1_hg38_v3.txt.gz`

## 6. Input requirements

- Supported input: TSV, VCF, VCF.GZ
- TSV required columns: #Chromosome, Position, Reference, Alteration
- SNV input expected (non-SNV rows are dropped)

## 7. PBS submission

Edit submit_mpact_scoring.pbs, then run:

```bash
qsub submit_mpact_scoring.pbs
```
