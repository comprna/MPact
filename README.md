# MPact m6A Variant Scoring Pipeline

Standalone pipeline to score SNVs with the MPact model and estimate m6A effect size.


## Quick Start (Copy-Paste)

```bash
cd /g/data/qq78/akanksha/m6A-snp/MPact_Scoring_Pipeline

# Create env (first time only)
python -m venv /g/data/qq78/akanksha/m6A-snp/.venv
source /g/data/qq78/akanksha/m6A-snp/.venv/bin/activate
pip install -r requirements.txt

# Run on bundled sample TSV
/g/data/qq78/akanksha/m6A-snp/.venv/bin/python score_mpact.py \
  --input mini_mixed.tsv \
  --output-tsv smoke_test_results.tsv \
  --fasta hg38.fa \
  --model-path model_window_501.h5 \
  --rediportal-gz TABLE1_hg38_v3.txt.gz \
  --scan-radius 5 \
  --batch-size 256
```

## Bundled Paths (All In This Directory)

- Script: `score_mpact.py`
- Model: `model_window_501.h5`
- FASTA: `hg38.fa`
- FASTA index: `hg38.fa.fai`
- GTF: `Homo_sapiens.GRCh38.110.gtf.gz`
- A-to-I reference: `TABLE1_hg38_v3.txt.gz`
- Sample input TSV: `mini_mixed.tsv`
- Sample input VCF: `mini_unannotated_test.vcf`

Required annotation files:
- `Homo_sapiens.GRCh38.110.gtf.gz` for strand inference
- `TABLE1_hg38_v3.txt.gz` for A-to-I annotation

## Inputs

### Supported file types
- TSV
- VCF
- VCF.GZ

### Required variant fields

TSV requires these columns (tab-separated):

| #Chromosome | Position | Reference | Alteration |
|-------------|----------|-----------|------------|
| chr1        | 100000   | A         | G          |

VCF uses standard `CHROM POS ID REF ALT` columns.

### Optional fields
- `transcript_id` / `ENST` to improve transcript-aware strand mapping.
- `strand` (if user-provided).
- Extra columns are preserved into output.

## Core Run Commands

### 1. Run on your own TSV

```bash
cd /g/data/qq78/akanksha/m6A-snp/MPact_Scoring_Pipeline
/g/data/qq78/akanksha/m6A-snp/.venv/bin/python score_mpact.py \
  --input /path/to/variants.tsv \
  --output-tsv /path/to/predictions.tsv \
  --fasta hg38.fa \
  --model-path model_window_501.h5 \
  --gtf Homo_sapiens.GRCh38.110.gtf.gz \
  --rediportal-gz TABLE1_hg38_v3.txt.gz \
  --output-plot /path/to/delta_histogram.png \
  --scan-radius 5 \
  --batch-size 512
```

### 2. Run on your own VCF

```bash
cd /g/data/qq78/akanksha/m6A-snp/MPact_Scoring_Pipeline
/g/data/qq78/akanksha/m6A-snp/.venv/bin/python score_mpact.py \
  --input /path/to/variants.vcf.gz \
  --output-tsv /path/to/predictions.tsv \
  --fasta hg38.fa \
  --model-path model_window_501.h5 \
  --gtf Homo_sapiens.GRCh38.110.gtf.gz \
  --rediportal-gz TABLE1_hg38_v3.txt.gz \
  --scan-radius 5 \
  --batch-size 512
```

### 3. Run on included sample (recommended first check)

```bash
cd /g/data/qq78/akanksha/m6A-snp/MPact_Scoring_Pipeline
/g/data/qq78/akanksha/m6A-snp/.venv/bin/python score_mpact.py \
  --input mini_mixed.tsv \
  --output-tsv sample_predictions.tsv \
  --fasta hg38.fa \
  --model-path model_window_501.h5 \
  --gtf Homo_sapiens.GRCh38.110.gtf.gz \
  --rediportal-gz TABLE1_hg38_v3.txt.gz \
  --output-plot sample_delta_hist.png
```

## Important Runtime Defaults

From current `score_mpact.py`:
- **`--gtf` is required** for strand inference and uses `Homo_sapiens.GRCh38.110.gtf.gz` by default.
- **`--rediportal-gz` is required** for A-to-I annotation and uses `TABLE1_hg38_v3.txt.gz` by default.
- Default `--scan-radius`: `5`
- Default `--batch-size`: `512`
- Default `--input-chunk-size`: `10000`
- VCF mode defaults to genic-only filtering (`--genic-only` is on).

### Resume an interrupted run

If a run stops partway through, rerun the same command with `--resume`.
You can also set `--resume-from-row` to continue from a specific input row and `--checkpoint-path` to point at the saved checkpoint JSON.

Example:

```bash
/g/data/qq78/akanksha/m6A-snp/.venv/bin/python score_mpact.py \
  --input mini_mixed.tsv \
  --output-tsv smoke_test_results.tsv \
  --fasta hg38.fa \
  --model-path model_window_501.h5 \
  --gtf Homo_sapiens.GRCh38.110.gtf.gz \
  --rediportal-gz TABLE1_hg38_v3.txt.gz \
  --resume
```

If you want to keep intergenic VCF rows too, add:

```bash
--allow-nongenic
```

## Output

Main output is a TSV with one row per scored candidate A-site.

### Header dictionary (all output columns)

The output contains:
- All normalized input variant columns (for example `#Chromosome`, `Position`, `Reference`, `Alteration`, and any additional input metadata columns)
- Plus the MPact-generated columns below

| Header | Meaning |
|---|---|
| strand_gencode | Final strand used for scoring (`+` or `-`), inferred from GTF transcript/interval logic. |
| strand_gencode_source | Label of the GTF reference used (for example `GRCh38.110`). |
| score_type | Candidate class: `disruption`, `creation`, or `context`. |
| a_genomic_pos1 | 1-based genomic coordinate of the candidate A center that was scored. |
| snp_to_a_mRNA_offset | Signed offset from SNP to candidate A in oriented transcript coordinates. |
| ref_center_is_A | Whether reference 501-nt oriented sequence has `A` at center position. |
| alt_center_is_A | Whether alternate 501-nt oriented sequence has `A` at center position. |
| overlaps_AtoI_exact | `True` if candidate A exactly overlaps an indexed A-to-I site; otherwise `False`. |
| near_AtoI_5nt | `True` if nearest A-to-I site is within 5 nt; otherwise `False`. |
| near_AtoI_10nt | `True` if nearest A-to-I site is within 10 nt; otherwise `False`. |
| mpact_ref_score | MPact model score on reference-oriented 501-nt sequence. |
| mpact_alt_score | MPact model score on alternate-oriented 501-nt sequence. |
| mpact_ref_stoichiometry_pct | Reference score scaled to percent (`mpact_ref_score * 100`). |
| mpact_alt_stoichiometry_pct | Alternate score scaled to percent (`mpact_alt_score * 100`). |
| mpact_delta_stoichiometry_pct | Percent delta (`(alt - ref) * 100`). |
| alt_center_A_destroyed | `True` when ALT no longer has center A (`not alt_center_is_A`). |
| mpact_delta_alt_minus_ref | Raw effect size: `mpact_alt_score - mpact_ref_score`. |
| mpact_abs_delta | Absolute effect size: `abs(mpact_delta_alt_minus_ref)`. |
| delta_zscore | Z-score of raw delta against global delta distribution for the run. |
| delta_p_two_sided | Two-sided p-value computed from `delta_zscore`. |
| ref10 | 10-nt centered reference context slice from scored window. |
| alt10 | 10-nt centered alternate context slice from scored window. |

Note:
- `delta_zscore` and `delta_p_two_sided` are computed after chunk scoring using global delta mean/std over the temporary scored output.

### Detailed interpretation guide

#### 1) Input passthrough columns

All input columns are preserved first in the output. This means original variant annotations (for example ClinVar labels, INFO-derived fields, cohort tags) stay attached to each scored candidate row.

Important behavior:
- One input SNV can produce multiple output rows, because MPact scans candidate A centers around the SNV (`--scan-radius`).
- The same input variant may therefore appear many times, each with a different `a_genomic_pos1` and possibly different `score_type`.

#### 2) Strand and coordinate columns

- `strand_gencode`: The strand actually used for sequence orientation and scoring. This is the critical strand field to trust for downstream analysis.
- `strand_gencode_source`: The GTF build label used to infer strand. Useful for reproducibility and cross-run auditing.
- `a_genomic_pos1`: Candidate m6A-centered genomic coordinate (1-based).
- `snp_to_a_mRNA_offset`: SNP-to-center distance in transcript orientation.

Offset sign meaning:
- Positive value: SNP is downstream of candidate A in the oriented transcript frame.
- Negative value: SNP is upstream of candidate A in the oriented transcript frame.
- Zero: SNP overlaps the candidate center position.

#### 3) Event typing columns

- `ref_center_is_A` and `alt_center_is_A` describe whether the center base is A before and after applying ALT.
- `score_type` is derived from these booleans:
  - `disruption`: REF center is A and ALT center is not A
  - `creation`: REF center is not A and ALT center is A
  - `context`: REF center is A and ALT center is A
- `alt_center_A_destroyed` is simply `not alt_center_is_A` and is most informative for disruptive events.

#### 4) A-to-I proximity columns

These annotate known editing context around the candidate center:
- `overlaps_AtoI_exact`: direct overlap
- `near_AtoI_5nt`: local neighborhood overlap within 5 nt
- `near_AtoI_10nt`: broader neighborhood overlap within 10 nt

Practical use:
- Use `overlaps_AtoI_exact == True` for strict known-site overlap analyses.
- Use `near_AtoI_5nt` or `near_AtoI_10nt` when testing local editing-environment enrichment.

#### 5) Score and delta columns

Core numeric fields:
- `mpact_ref_score`: model score on REF sequence
- `mpact_alt_score`: model score on ALT sequence
- `mpact_delta_alt_minus_ref`: ALT minus REF (main direction-aware effect size)
- `mpact_abs_delta`: absolute magnitude of effect

Percent fields are linear transforms of score fields:
- `mpact_ref_stoichiometry_pct = mpact_ref_score * 100`
- `mpact_alt_stoichiometry_pct = mpact_alt_score * 100`
- `mpact_delta_stoichiometry_pct = mpact_delta_alt_minus_ref * 100`

Direction interpretation for `mpact_delta_alt_minus_ref`:
- Negative: ALT decreases predicted m6A signal relative to REF
- Positive: ALT increases predicted m6A signal relative to REF
- Near zero: limited predicted effect for that candidate center

#### 6) Statistical normalization columns

- `delta_zscore`: standardized delta across all candidate rows in the run
- `delta_p_two_sided`: two-sided p-value from that z-score

These are run-level normalized values, so they depend on the cohort/distribution in that specific run. If you compare different runs, compare carefully because z-score baselines can shift.

#### 7) Local context columns

- `ref10` and `alt10` are short centered context strings from the full 501-nt windows.
- They are useful for quick motif-level sanity checks and visual inspection.
- They are not a replacement for full-window model inputs, but they are helpful for debugging and reporting.

#### 8) Row-level interpretation example

If a row has:
- `score_type = disruption`
- `mpact_delta_alt_minus_ref = -0.42`
- `mpact_abs_delta = 0.42`
- `delta_p_two_sided = 0.001`

Then that candidate center is predicted to show a strong ALT-driven loss of m6A signal, with relatively extreme effect size within the run distribution.

Optional histogram PNG is written when `--output-plot` is set.

## Test command


```bash
cd /g/data/qq78/akanksha/m6A-snp/MPact_Scoring_Pipeline
/g/data/qq78/akanksha/m6A-snp/.venv/bin/python score_mpact.py \
  --input mini_mixed.tsv \
  --output-tsv smoke_test_results.tsv \
  --fasta hg38.fa \
  --model-path model_window_501.h5 \
  --scan-radius 5 \
  --batch-size 256
```

Observed result:
- Exit code `0`
- Output TSV created successfully
- 5 lines total (header + 4 scored candidates)

## HPC PBS Usage

Template script: `submit_mpact_scoring.pbs`

Before `qsub`, update:
- `#PBS -P`, queue, walltime, ncpus, mem, storage
- `INPUT_PATH`, `OUTPUT_DIR`, `FASTA`, `MODEL`
- optional `REDIPORTAL_GZ`, `GTF`

Submit with:

```bash
qsub submit_mpact_scoring.pbs
```

## Resume Note

The run is resumable with the updated PBS wrapper:

```bash
qsub /g/data/qq78/akanksha/m6A-snp/Scripts/run_gnomad_genic_snv_dtm6a501_scan20.pbs
```

The wrapper calls `score_mpact.py` with `--resume`, so already-scored chunk outputs in [external_validation/gnomad/chrom_chunks](../external_validation/gnomad/chrom_chunks) are skipped and only missing chunks are processed.

The launcher now passes absolute paths for:
- FASTA: `MPact_Scoring_Pipeline/hg38.fa`
- GTF: `MPact_Scoring_Pipeline/Homo_sapiens.GRCh38.110.gtf.gz`
- REDIportal: `MPact_Scoring_Pipeline/TABLE1_hg38_v3.txt.gz`

## Troubleshooting

### Missing FASTA index

```bash
samtools faidx /path/to/hg38.fa
```

### `samtools` not found

```bash
module load samtools
which samtools
```

### Slow startup or no GPU on login node

TensorFlow may print CUDA warnings on CPU-only nodes. This is expected if no GPU is present, and CPU inference still runs.

### Memory pressure

Use a smaller batch size:

```bash
--batch-size 256
```
