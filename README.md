# MPact: m6A site and variant scoring

MPact predicts m6A potential at candidate adenosines and how SNVs change that
potential. The bundled classification-only CNN–BiLSTM–attention models use
strand-oriented genomic sequence and positional encoding. They do not predict
stoichiometry. Accessibility and conservation are optional output annotations,
not model inputs.

## Installation

```bash
git clone https://github.com/comprna/MPact.git
cd MPact
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Only for optional accessibility / conservation:
pip install -r requirements-annotations.txt
```

Provide an indexed GRCh38 FASTA (`samtools faidx hg38.fa`). Variant scoring also
requires a matching GTF and REDIportal A-to-I table. Large reference resources
are not distributed in GitHub. Default local filenames are
`Homo_sapiens.GRCh38.110.gtf.gz` and `TABLE1_hg38_v3.txt.gz`; override their paths
as needed. Models expect TensorFlow 2.15; CPU inference is supported.

## Models

The current Models are:

| Window | Validation AUPRC | Validation AUROC |
|---|---:|---:|
| 101 | 0.9947 | 0.9962 |
| 201 | 0.9966 | 0.9975 |
| **501 (default)** | 0.9975 | 0.9980 |
| 801 | 0.9979 | 0.9983 |
| 1001 | 0.9983 | 0.9986 |


## Score SNVs

```bash
python score_mpact.py \
  --input variants.tsv \
  --output-tsv predictions.tsv \
  --fasta /path/to/hg38.fa \
  --gtf /path/to/annotation.gtf.gz \
  --rediportal-gz /path/to/REDIportal.txt.gz \
  --batch-size 64
```

TSV input requires `#Chromosome`, `Position` (1-based), `Reference` and
`Alteration`; VCF and VCF.GZ are also supported. Extra TSV input metadata, including
ClinVar `CLNSIG`, are retained. Transcript IDs and gene symbols improve
GTF-based strand inference. Input strand is retained as metadata; inferred
`strand_gencode` is used for scoring. Unresolved strand, inconsistent reference
alleles or insufficient context are not assigned valid predictions.

VCF input defaults to genic-only filtering. `--allow-nongenic` disables that
filter but does not bypass strand or context requirements.

The default scan is ±5 nt around each SNV (`--scan-radius 5`). A row represents
a variant/candidate-adenosine pair, so one SNV can produce multiple rows.
Context is strand-oriented genomic sequence, not spliced transcript sequence,
matching training.

```bash
# Select the bundled 1001-nt model:
python score_mpact.py --input variants.tsv --output-tsv predictions_1001.tsv \
  --fasta /path/to/hg38.fa --window-size 1001 \
  --gtf /path/to/annotation.gtf.gz --rediportal-gz /path/to/REDIportal.txt.gz

# Optional annotations; either flag can also be used alone:
python score_mpact.py --input variants.tsv --output-tsv annotated.tsv \
  --fasta /path/to/hg38.fa \
  --gtf /path/to/annotation.gtf.gz --rediportal-gz /path/to/REDIportal.txt.gz \
  --accessibility --conservation-bigwig /path/to/hg38.phyloP100way.bw
```

## Score candidate m6A sites

`m6A_scorer.py` takes a site TSV containing chromosome, 1-based position and
strand. It recognizes `chr`/`chrom`/`#Chromosome`, `start`/`Position`/`POS`, and
`strand`/`Strand`; override with `--chrom-col`, `--pos-col`, `--strand-col`.
If an end column is present, it must equal the start. Convert BED coordinates
to this convention before scoring.

```bash
python m6A_scorer.py \
  --input mini_m6A_sites.tsv \
  --output-tsv sites.scored.tsv \
  --summary-json sites.summary.json \
  --summary-tsv sites.summary.tsv \
  --fasta /path/to/hg38.fa \
  --batch-size 64
```

The same model/window and optional-annotation flags are supported.
`--accessibility-processes N` retains parallel folding support, but only when
`--accessibility` is enabled. Non-A centers receive a missing site score, not a
negative classification. Invalid/unfetchable rows are counted in the summary.

## Compact output

Input columns come first. The variant scorer adds:

| Columns | Meaning |
|---|---|
| `strand_gencode` | Strand used for scoring |
| `score_type` | Center-base event: creation, disruption or context |
| `a_genomic_pos1`, `snp_to_a_mRNA_offset` | Candidate coordinate and oriented variant offset |
| `ref_center_is_A`, `alt_center_is_A` | Whether each allele retains a candidate adenosine |
| `mpact_ref_score`, `mpact_alt_score` | Reference and alternative classification scores |
| `mpact_delta_alt_minus_ref`, `mpact_abs_delta` | Signed and absolute score differences |
| `ref_motif5`, `alt_motif5` | Centered 5-nt contexts |
| `ref_is_DRACH`, `alt_is_DRACH`, `ref_is_RAC`, `alt_is_RAC` | Motif membership |
| `overlaps_AtoI_exact`, `near_AtoI_5nt`, `near_AtoI_10nt` | A-to-I flags; sites are not filtered |

The site scorer adds `mpact_score`, `center_base_transcript`,
`scoreable_center_A`, `motif5`, `is_DRACH` and `is_RAC`.


`score_type` describes the nucleotide event.


Optional annotation columns appear only when enabled:

- `--accessibility` adds reference/alternative/delta values for center-base
  unpaired probability and centered 5-, 10- and 20-nt accessibility. Site
  scoring adds only reference values. ViennaRNA parameters are W=200 and
  maximum span=150, capped to sequence length.
- `--conservation-bigwig PATH` adds `variant_conservation_score` and
  `a_site_conservation_score`; site scoring adds `site_conservation_score`.
  Missing track coverage remains missing, not zero. The track must match the
  reference assembly. `--conservation-label` sets a provenance label.
- `--full-output` adds strand/context.
  


## Resume and provenance

Rerun interrupted variant commands with `--resume`, the same output filename
and settings. Retain `.raw_scored.tsv`, `.checkpoint.json` and `.run.json`
sidecars. 
`--keep-raw-temp` retains the raw TSV after successful completion.
`--output-plot plot.png` adds a delta histogram.

Reduce `--batch-size` or `--input-chunk-size` for lower memory use. The PBS
template `submit_mpact_scoring.pbs` must be adapted to your project, queue,
environment and reference paths.

## Tests

```bash
python -m unittest discover -s tests -p 'test_*.py'
# Requires annotation dependencies and local hg38.fa:
bash tests/test_m6A_scorer_smoke.sh
```
