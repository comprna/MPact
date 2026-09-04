# MPact m6A Variant Scoring Pipeline

Standalone pipeline to score SNVs with the MPact model and estimate m6A effect size.


## Quick Start (Copy-Paste)

```bash
git clone https://github.com/comprna/MPact.git
cd MPact

# Create env (first time only)
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Run on bundled sample TSV
python score_mpact.py \
  --input mini_scoreable.tsv \
  --output-tsv smoke_test_results.tsv \
  --fasta hg38.fa \
  --model-path models/mpact_dtm6a_501nt_seed42.keras \
  --window-size 501 \
  --conservation-bigwig /path/to/hg38.phyloP_or_phastCons.bw \
  --rediportal-gz TABLE1_hg38_v3.txt.gz \
  --scan-radius 5 \
  --batch-size 256
```

## Bundled Paths (All In This Directory)

- Script: `score_mpact.py` for SNV effect scoring
- Script: `m6A_scorer.py` for direct scoring of known or candidate m6A-centered sites
- Recommended model: `models/mpact_dtm6a_501nt_seed42.keras`
- Current DTM6A models for every evaluated window:
  - `models/mpact_dtm6a_101nt_seed42.keras`
  - `models/mpact_dtm6a_201nt_seed42.keras`
  - `models/mpact_dtm6a_501nt_seed42.keras`
  - `models/mpact_dtm6a_801nt_seed42.keras`
  - `models/mpact_dtm6a_1001nt_seed42.keras`
- Checkpoint hashes, validation metrics, and thresholds: `models/model_manifest.json`
- FASTA: `hg38.fa`
- FASTA index: `hg38.fa.fai`
- GTF: `Homo_sapiens.GRCh38.110.gtf.gz`
- A GRCh38 phyloP or phastCons bigWig supplied with `--conservation-bigwig`
- A-to-I reference: `TABLE1_hg38_v3.txt.gz`
- Sample scoreable TSV: `mini_scoreable.tsv`
- Sample m6A-site TSV: `mini_m6A_sites.tsv`
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
cd MPact
python score_mpact.py \
  --input /path/to/variants.tsv \
  --output-tsv /path/to/predictions.tsv \
  --fasta hg38.fa \
  --model-path models/mpact_dtm6a_501nt_seed42.keras \
  --window-size 501 \
  --conservation-bigwig /path/to/hg38.phyloP_or_phastCons.bw \
  --gtf Homo_sapiens.GRCh38.110.gtf.gz \
  --rediportal-gz TABLE1_hg38_v3.txt.gz \
  --output-plot /path/to/delta_histogram.png \
  --scan-radius 5 \
  --batch-size 512
```

### 2. Run on your own VCF

```bash
cd MPact
python score_mpact.py \
  --input /path/to/variants.vcf.gz \
  --output-tsv /path/to/predictions.tsv \
  --fasta hg38.fa \
  --model-path models/mpact_dtm6a_501nt_seed42.keras \
  --window-size 501 \
  --conservation-bigwig /path/to/hg38.phyloP_or_phastCons.bw \
  --gtf Homo_sapiens.GRCh38.110.gtf.gz \
  --rediportal-gz TABLE1_hg38_v3.txt.gz \
  --scan-radius 5 \
  --batch-size 512
```

### 3. Score known m6A sites directly

Use `m6A_scorer.py` when the input is a table of known or candidate m6A-centered genomic sites, not SNVs. The script fetches the selected sequence window around each site, orients it using the strand column, checks that the transcript-oriented center base is `A`, and reports the MPact site score, accessibility, and conservation.

The input TSV must contain chromosome, 1-based position, and strand columns. Recognized defaults include `chr`, `chrom`, `#Chromosome`, or `CHROM`; `start`, `Position`, `POS`, or `pos`; and `strand` or `Strand`. Use `--chrom-col`, `--pos-col`, and `--strand-col` for other column names. If an `end` column is present, start and end must match.

```bash
cd MPact
python m6A_scorer.py \
  --input mini_m6A_sites.tsv \
  --output-tsv mini_m6A_sites.scored.tsv \
  --summary-json mini_m6A_sites.summary.json \
  --summary-tsv mini_m6A_sites.summary.tsv \
  --fasta hg38.fa \
  --model-path models/mpact_dtm6a_501nt_seed42.keras \
  --window-size 501 \
  --conservation-bigwig /path/to/hg38.phyloP_or_phastCons.bw \
  --batch-size 4
```

Run the bundled smoke test with:

```bash
bash tests/test_m6A_scorer_smoke.sh
```

### 4. Run on included variant sample (recommended first check)

```bash
cd MPact
python score_mpact.py \
  --input mini_scoreable.tsv \
  --output-tsv sample_predictions.tsv \
  --fasta hg38.fa \
  --model-path models/mpact_dtm6a_501nt_seed42.keras \
  --window-size 501 \
  --conservation-bigwig /path/to/hg38.phyloP_or_phastCons.bw \
  --gtf Homo_sapiens.GRCh38.110.gtf.gz \
  --rediportal-gz TABLE1_hg38_v3.txt.gz \
  --output-plot sample_delta_hist.png
```

## Important Runtime Defaults

From current `score_mpact.py`:
- **`--conservation-bigwig` is required** and must match the reference assembly.
- Accessibility is always calculated with ViennaRNA RNAplfold.
- **`--gtf` is required** for strand inference and uses `Homo_sapiens.GRCh38.110.gtf.gz` by default.
- **`--rediportal-gz` is required** for A-to-I annotation and uses `TABLE1_hg38_v3.txt.gz` by default.
- Default `--scan-radius`: `5`
- Default `--batch-size`: `1024`
- Default `--input-chunk-size`: `50000`
- VCF mode defaults to genic-only filtering (`--genic-only` is on).

### Resume an interrupted run

If a run stops partway through, rerun the same command with `--resume`.
You can also set `--resume-from-row` to continue from a specific input row and `--checkpoint-path` to point at the saved checkpoint JSON.

Example:

```bash
python score_mpact.py \
  --input mini_scoreable.tsv \
  --output-tsv smoke_test_results.tsv \
  --fasta hg38.fa \
  --model-path models/mpact_dtm6a_501nt_seed42.keras \
  --window-size 501 \
  --conservation-bigwig /path/to/hg38.phyloP_or_phastCons.bw \
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
| ref_center_is_A | Whether the reference-oriented sequence has `A` at its center. |
| alt_center_is_A | Whether the alternate-oriented sequence has `A` at its center. |
| overlaps_AtoI_exact | `True` if candidate A exactly overlaps an indexed A-to-I site; otherwise `False`. |
| near_AtoI_5nt | `True` if nearest A-to-I site is within 5 nt; otherwise `False`. |
| near_AtoI_10nt | `True` if nearest A-to-I site is within 10 nt; otherwise `False`. |
| mpact_ref_score | MPact model score on the reference-oriented sequence. |
| mpact_alt_score | MPact model score on the alternate-oriented sequence. |
| alt_center_A_destroyed | `True` when ALT no longer has center A (`not alt_center_is_A`). |
| mpact_delta_alt_minus_ref | Raw effect size: `mpact_alt_score - mpact_ref_score`. |
| mpact_abs_delta | Absolute effect size: `abs(mpact_delta_alt_minus_ref)`. |
| delta_zscore | Z-score of raw delta against global delta distribution for the run. |
| delta_p_two_sided | Two-sided p-value computed from `delta_zscore`. |
| ref_scan_seq | Centered reference context slice from the scored window; length follows `--scan-radius`, capped between 5 and 21 nt. |
| accessibility_source | ViennaRNA RNAplfold parameter label. |
| mpact_ref/alt/delta_unpaired_probability_1nt | Center-base unpaired probabilities for REF, ALT, and ALT minus REF. |
| mpact_ref/alt/delta_accessibility_5nt | Centered 5-nt unpaired probabilities for REF, ALT, and their difference. |
| mpact_ref/alt/delta_accessibility_10nt | Centered 10-nt unpaired probabilities for REF, ALT, and their difference. |
| mpact_ref/alt/delta_accessibility_20nt | Centered 20-nt unpaired probabilities for REF, ALT, and their difference. |
| conservation_source | User-supplied phyloP or phastCons track label. |
| variant_conservation_score | Conservation score at the input variant coordinate. |
| a_site_conservation_score | Conservation score at the candidate adenosine coordinate. |
| alt_scan_seq | Centered alternate context slice from the scored window; length follows `--scan-radius`, capped between 5 and 21 nt. |

Note:
- `delta_zscore` and `delta_p_two_sided` are computed after chunk scoring using global delta mean/std over the temporary scored output.
- Accessibility and conservation are reported annotations; they are not inputs to the released sequence-only classifiers.

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

The released models report classification scores and variant deltas only. Stoichiometry is intentionally not predicted because the experimental model did not generalize reliably.

Direction interpretation for `mpact_delta_alt_minus_ref`:
- Negative: ALT decreases predicted m6A signal relative to REF
- Positive: ALT increases predicted m6A signal relative to REF
- Near zero: limited predicted effect for that candidate center

#### 6) Statistical normalization columns

- `delta_zscore`: standardized delta across all candidate rows in the run
- `delta_p_two_sided`: two-sided p-value from that z-score

These are run-level normalized values, so they depend on the cohort/distribution in that specific run. If you compare different runs, compare carefully because z-score baselines can shift.

#### 7) Local context columns

- `ref_scan_seq` and `alt_scan_seq` are short centered context strings from the full 501-nt windows.
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
cd MPact
python score_mpact.py \
  --input mini_scoreable.tsv \
  --output-tsv smoke_test_results.tsv \
  --fasta hg38.fa \
  --model-path models/mpact_dtm6a_501nt_seed42.keras \
  --window-size 501 \
  --conservation-bigwig /path/to/hg38.phyloP_or_phastCons.bw \
  --scan-radius 5 \
  --batch-size 256
```

Observed result:
- Exit code `0`
- Output TSV created successfully
- 6 lines total (header + 5 scored candidates)

## HPC PBS Usage

Template script: `submit_mpact_scoring.pbs`

Before `qsub`, update:
- `#PBS -P`, queue, walltime, ncpus, mem, storage
- `INPUT_PATH`, `OUTPUT_DIR`, `FASTA`, `MODEL`, `WINDOW_SIZE`, and `CONSERVATION_BIGWIG`
- `REDIPORTAL_GZ` and `GTF`

Submit with:

```bash
qsub submit_mpact_scoring.pbs
```

## Resume Note

Runs are resumable with `--resume`. Rerun the same command and keep the same `--output-tsv`; the scorer will use the checkpoint and temporary scored TSV to continue interrupted work.

```bash
python score_mpact.py \
  --input mini_scoreable.tsv \
  --output-tsv smoke_test_results.tsv \
  --fasta hg38.fa \
  --model-path models/mpact_dtm6a_501nt_seed42.keras \
  --window-size 501 \
  --conservation-bigwig /path/to/hg38.phyloP_or_phastCons.bw \
  --resume
```

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
