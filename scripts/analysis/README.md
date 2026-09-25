# Chromosome 15 cis-methylation analysis

This directory separates inference from figure rendering. The primary unit is a
participant, not a CpG or genomic window.

## Evidence rules

- PWS-DEL combined methylation is maternal-retained only inside that
  participant's independently confirmed, buffered CN=1 interval.
- AS-DEL combined methylation is paternal-retained only inside that
  participant's independently confirmed, buffered CN=1 interval.
- Direct PWS-versus-AS parental contrasts use only the intersection of all
  buffered deletion intervals.
- H1 and H2 in unaffected and DiGeorge samples remain numerical, unoriented
  haplotypes. Their absolute difference is used as phase-invariant ASM.
- The control/DiGeorge combined signal estimates the shared diploid scaffold.
- Outside validated CN=1 sequence, no observed maternal or paternal label is
  assigned without parental DNA.

## Run

Each script contains its paths and analysis constants at the top. Generate
Figure 1 so its deletion provenance table exists, verify the GTF path near the
top of `00_build_prespecified_regions.py`, and then run everything without
command-line arguments:

```bash
python3 scripts/analysis/run_analysis.py
```

The individual stages can also be run separately:

```bash
python3 scripts/analysis/00_build_prespecified_regions.py
python3 scripts/analysis/01_build_chr15_evidence_matrix.py
python3 scripts/analysis/02_depth_missingness_sensitivity.py
python3 scripts/analysis/03_reciprocal_cis_architecture.py
```

## Script responsibilities

| Script | Scientific purpose | Principal outputs |
|---|---|---|
| `00_build_prespecified_regions.py` | Builds the locus catalog from the T2T GTF before testing effects, preventing selection of genes after seeing the methylation results. | `prespecified_regions.tsv` |
| `01_build_chr15_evidence_matrix.py` | Converts pb-CpG-tools tracks into a validated participant-by-window evidence matrix and applies the CN=1 rules that determine where parental direction is directly observable. | Evidence matrix, track inventory, validated CN=1 intervals, common deletion core and depth QC |
| `02_depth_missingness_sensitivity.py` | Tests whether sequencing depth, CpG recovery, filtering threshold or read weighting changes the reciprocal PWS/AS contrast. | Per-sample QC, per-group QC, depth-recovery associations and estimator concordance |
| `03_reciprocal_cis_architecture.py` | Estimates the shared diploid scaffold, phase-invariant ASM and the direct maternal-retained minus paternal-retained contrast with participant bootstrap intervals. | Window architecture, participant ASM and regional effects |
| `run_analysis.py` | Runs scripts `00` through `03` in dependency order. | All analysis outputs |

`cis_analysis/` contains reusable parsers, evidence rules, window construction
and bootstrap utilities. It contains no server paths or cohort-specific
configuration.

## Outputs

`01_evidence_matrix` contains one row per participant, track, and 1-kb window,
including the evidence class, CpG count, depth, site-weighted beta, read-weighted
beta, depth-capped sensitivity estimates, and CN=1 eligibility.

`02_depth_sensitivity` quantifies CpG recovery, depth dependence, and the
stability of the reciprocal-deletion contrast across CpG thresholds and beta
estimators. The site-weighted beta is primary, so a deeply sequenced CpG does
not receive greater biological weight than a lower-depth CpG.

`03_cis_architecture` produces three distinct signals for downstream figures:

1. shared diploid scaffold from combined unaffected and DiGeorge tracks;
2. phase-invariant absolute ASM from H1/H2 in those diploid references;
3. directly observed maternal-retained minus paternal-retained contrast in the
   common reciprocal CN=1 interval, with participant-bootstrap intervals.

## Interpretation limits

- H1 and H2 are never converted into global maternal and paternal labels.
- A parental direction outside the common reciprocal CN=1 core is not reported
  as observed.
- Depth-capped estimators limit the influence of high-depth CpGs but cannot
  reconstruct CpGs that pb-CpG-tools did not emit.
- Genomic windows are technical measurement units; participant resampling is
  used for biological uncertainty.
