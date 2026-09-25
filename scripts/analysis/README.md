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
- Participants are the biological replicates. Group means weight every
  participant equally, and all intervals resample participants, never windows
  or CpGs.
- Methylation beta values are never divided by or weighted towards sequencing
  depth. The 15x/30x difference is handled with common-depth read thinning,
  capped effective coverage, minimum effective observations, a shared-CpG
  sensitivity analysis and explicit missingness.
- Smoothing never crosses a non-evaluable 1-kb window.
- DiGeorge participants are an independent diploid disease-control cohort, not
  parental references.
- Outside validated CN=1 sequence, no observed maternal or paternal label is
  assigned without parental DNA.

## Run

Each script contains its paths and analysis constants at the top. Generate
Figure 1 so its deletion provenance table exists, verify the GTF path near the
top of `00_build_prespecified_regions.py` (the first existing entry of
`GTF_CANDIDATES` is used), and then run everything without command-line
arguments:

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
| `01_build_chr15_evidence_matrix.py` | Converts pb-CpG-tools tracks into a validated participant-by-window evidence matrix, applies the CN=1 rules that determine where parental direction is directly observable, and adds full-depth, common-depth downsampled, capped-coverage and shared-CpG estimates. | Evidence matrix, track inventory, validated CN=1 intervals, common deletion core and depth QC |
| `02_depth_missingness_sensitivity.py` | Tests whether sequencing depth, CpG recovery, filtering threshold or read weighting changes the reciprocal PWS/AS contrast. | Per-sample QC, per-group QC, depth-recovery associations and estimator concordance |
| `03_reciprocal_cis_architecture.py` | Estimates the PWS maternal-retained minus AS paternal-retained contrast in the common CN=1 interval (windows, focal intervals and prespecified regions) and regional phase-invariant ASM in controls and DiGeorge participants. | All Figure 2 tables |
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

`03_cis_architecture` contains the Figure 2 inputs:

| File | Content |
|---|---|
| `parent_associated_windows.tsv.gz` | Per 1-kb window: equally weighted PWS and AS retained-copy means, `delta_beta` (PWS - AS), participant counts, participant-bootstrap CI, genomic support, effective observations, missingness, common CN=1 flag, contiguous segment, gap-safe 21-kb median with bootstrap ribbon, and the same contrast for every sensitivity estimator. |
| `parent_window_participant_values.tsv.gz` | Participant-by-window retained-copy beta inside the common CN=1 interval. |
| `focal_intervals.tsv` | Candidate focal intervals and whether each meets the prespecified window, participant, effect-size, leave-one-participant-out and sensitivity criteria. |
| `regional_parent_contrasts.tsv` | One row per prespecified region: group means, `delta_beta`, participant-bootstrap and Welch 95% CIs, participant counts, evaluable windows and CpGs, status (`maternal_retained_higher`, `paternal_retained_higher`, `inconclusive`, `not_evaluable`), and coverage and shared-CpG sensitivity. |
| `regional_parent_contrasts_by_estimator.tsv` | The regional contrast for each estimator. |
| `regional_participant_values.tsv.gz` | One regional value per participant, region and estimator for both analyses. |
| `regional_phase_invariant_asm.tsv` | Regional mean \|beta_H1 - beta_H2\| for controls and DiGeorge participants separately; intervals only for n >= 3 and labelled descriptive. |
| `participant_missingness.tsv` | Fraction of common CN=1 windows evaluable per participant and estimator, with retained-copy depth. |
| `figure2_analysis_report.tsv` | Input validation, interval restriction, bootstrap unit, labelling, gap-safety checks and the main estimates. |

The run stops with an error if any contrast is estimated outside the common
CN=1 interval, if a control or DiGeorge row carries a parental label, or if a
smoothed value spans a non-evaluable window.

## Interpretation limits

- H1 and H2 are never converted into maternal and paternal labels, including by
  orienting phase blocks at the imprinting centre.
- Absolute ASM is positive even without allelic methylation because of
  sampling noise; the OCA2 downstream control provides the empirical baseline.
- Percentile bootstrap intervals with three AS participants or two controls
  are descriptive.
- A parental direction outside the common reciprocal CN=1 core is not reported
  as observed.
- Depth-capped estimators limit the influence of high-depth CpGs but cannot
  reconstruct CpGs that pb-CpG-tools did not emit.
- Genomic windows are technical measurement units; participant resampling is
  used for biological uncertainty.
