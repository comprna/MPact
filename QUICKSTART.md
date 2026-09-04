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
- Recommended model: models/mpact_dtm6a_501nt_seed42.keras
- Current models: models/mpact_dtm6a_{101,201,501,801,1001}nt_seed42.keras
- Model metadata and SHA-256 checksums: models/model_manifest.json
- Required conservation input: a GRCh38 phyloP or phastCons bigWig
- Accessibility is always calculated with ViennaRNA RNAplfold
- The public pipeline does not report stoichiometry predictions
- FASTA: hg38.fa
- GTF: Homo_sapiens.GRCh38.110.gtf.gz

## 3. Run on bundled sample (recommended smoke test)

```bash
cd MPact
python score_mpact.py \
  --input mini_scoreable.tsv \
  --output-tsv sample_predictions.tsv \
  --fasta hg38.fa \
  --model-path models/mpact_dtm6a_501nt_seed42.keras \
  --window-size 501 \
  --conservation-bigwig /path/to/hg38.phyloP_or_phastCons.bw \
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
  --model-path models/mpact_dtm6a_501nt_seed42.keras \
  --window-size 501 \
  --conservation-bigwig /path/to/hg38.phyloP_or_phastCons.bw \
  --output-plot /path/to/delta_histogram.png
```

### VCF input

```bash
cd MPact
python score_mpact.py \
  --input /path/to/variants.vcf.gz \
  --output-tsv /path/to/predictions.tsv \
  --fasta hg38.fa \
  --model-path models/mpact_dtm6a_501nt_seed42.keras \
  --window-size 501 \
  --conservation-bigwig /path/to/hg38.phyloP_or_phastCons.bw \
  --gtf Homo_sapiens.GRCh38.110.gtf.gz
```

## 5. Required and optional flags

**Required:**
- `--gtf Homo_sapiens.GRCh38.110.gtf.gz` (or custom) — Strand inference from GTF is mandatory for correctness. Bundled default is recommended.

- `--rediportal-gz /path/to/TABLE1_hg38_v3.txt.gz` — A-to-I annotations are mandatory.
- `--conservation-bigwig /path/to/track.bw` — The track must use the same assembly as the FASTA.
**Optional:**
- Keep intergenic VCF variants: `--allow-nongenic`
- Smaller memory footprint: `--batch-size 256`
- Larger neighborhood scan: `--scan-radius 20`

## 6. Input requirements

- Supported input: TSV, VCF, VCF.GZ
- TSV required columns: #Chromosome, Position, Reference, Alteration
- SNV input expected (non-SNV rows are dropped)

## 7. PBS submission

Edit submit_mpact_scoring.pbs, then run:

```bash
qsub submit_mpact_scoring.pbs
```
