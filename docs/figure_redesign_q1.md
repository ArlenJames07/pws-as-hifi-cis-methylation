# Figure redesign for the 15q11–q13 resubmission

**Working figure map.** This is a design and evidence contract for the current
17-person cohort (5 PWS-DEL, 3 AS-DEL, 1 PWS-mUPD, 2 unaffected controls, 6
22q11.2-deletion/DiGeorge controls). It does not claim that unpublished cohort
outputs were inspected: participant-level result tables and rendered images are
local, not tracked in this repository. Use T2T-CHM13v2.0 coordinates throughout.

The argument has five steps: establish the molecular configuration; locate the
parent-associated methylation divergence; map the focal transition; inspect
multi-CpG patterns on reads; assess how structural context and repeats relate to
the methylation architecture. Keep the biological participant as the inference
unit wherever a group contrast is made.

| Main figure | Question and suggested panels | Current source | Manuscript statement that the panel can support |
|---|---|---|---|
| **1. Molecular identity** | (a) locus and parent-of-origin schematic; (b) chr15 dosage/confirmed deletions; (c) IC methylation and chromosome-15 configuration per participant. Show both unaffected and all six DiGeorge controls. | `FIGURE_1.py`, including deletion provenance and IC validation. | Reciprocal deletions leave opposite parental chromosomes; mUPD retains two maternal-like copies. This is the sole diagnostic classification display. |
| **2. Direct cis architecture** | (a) signed PWS-DEL maternal-retained minus AS-DEL paternal-retained Δβ in the shared, buffered CN=1 core with unevaluable gaps; (b) prespecified regional effects with participant-bootstrap intervals and actual group n; (c) control and DiGeorge **unoriented** \|H1−H2\|. | `analysis/01–03`; `FIGURE_2.py`. | Parental divergence is focal against a substantial shared methylation background. Only (a,b) identify parental direction, and only within the common CN=1 interval. |
| **3. Focal transition** | (a) one locus-level zoom showing raw CpGs plus window support for the IC/SNHG14 region; (b) participant-specific transition intervals and a clearly distinguished consensus; (c) change-point/bootstrap/threshold robustness; (d) annotation and distance to BP/SD landmarks. | `FIGURE_3.py` Phase 3 tables. | A local methylation transition may recur across configurations. Its position and uncertainty must be stated separately from structural deletion breakpoints. |
| **4. Single-molecule coordination** | (a) small key defining read/CpG states; (b) matched SNORD116 reads with missing calls visibly missing; (c) per-participant methylation fraction or within-read heterogeneity, showing every donor; (d) participant-level AS-DEL versus PWS-DEL effect. | `FIGURE_4.py` and its nine precomputed inputs. | HiFi reads exhibit local multi-CpG patterns; reads are nested within people. Molecule counts alone cannot support a population-level p-value or domain-wide spreading. |
| **5. Duplicon context and repeat methylation** | (a) per-carrier breakpoint evidence matrix distinguishing SD geometry, crossover excess, homologous assembly, split reads and final verdict; (b) signed *raw* repeat-element Δβ by 250-kb bin inside the common CN=1 core; (c) PWS-DEL, AS-DEL and the single mUPD compared descriptively with the biparental combined-track median. | `duplicons/05–08`; new `FIGURE_5_DUPLICON_REPEAT.py`. | Structural evidence has differing resolution across carriers; parent-associated methylation can be assessed in testable repeats as well as genes. A compatible SD geometry alone is not a confirmed NAHR event. |

## Suggested visual rules

- Use a stable orange for maternal-retained, blue for paternal-retained, purple
  for mUPD; use marker shape or text so direction remains legible in grayscale.
- Put the sample count in each panel or regional row. Leave missing CpGs,
  unevaluable windows, nonassembled contigs and untested junctions visibly
  missing. Avoid zero filling and smoothing across missing intervals.
- A parent-of-origin **observation** is restricted to the independently
  confirmed CN=1 intersection in PWS-DEL and AS-DEL. Numerical H1/H2 in
  Control/DiGeorge are unoriented. IC- or assembly-anchored orientation outside
  this interval is a separate, explicitly labeled inference.
- Sign of Δβ means PWS-DEL maternal retained minus AS-DEL paternal retained;
  positive is maternal-higher. Use the raw difference in Figure 5b; its
  additional centered difference is a sensitivity view and may change the
  apparent direction of a small effect.
- Keep the 1-kb methylation analysis distinct from the 250-kb repeat summary.
  Figure 5 repeat bins are medians of eligible reference-annotated repeat
  elements, not direct measurements of every base or all repeats in a bin.
- Mark nominal 250-kb permutation results as **nominal**. With eight deletion
  carriers there are only 56 possible 5/3 label assignments; across many bins,
  a per-bin p ≤ 0.05 is not a familywise or FDR result. Do not describe the
  derived domains as independent discoveries without a multiplicity strategy.

## Supplementary placement

| Source | Where to place it | Reason |
|---|---|---|
| Figure 2 participant heatmap, retained-copy profiles, diploid combined scaffold, depth/downsampling/shared-CpG checks and full regional estimates | Supplementary Fig. S2 and tables | These validate panel 2 while keeping its direct contrast readable. |
| Figure 3 complete boundary intervals, single-CpG plots, all threshold sweeps and annotation provenance | Supplementary Fig. S3 and tables | One main transition zoom suffices; full uncertainty stays accessible. |
| Figure 4 extended read stacks, CpG recovery, all read-level statistics | Supplementary Fig. S4 and tables | Show the nesting and missingness rather than treating each read as a donor. |
| Legacy `FIGURE_5.py` genome-wide CNV/SV burdens and breakpoint-flank distance decay; SUNK counts, all junction-like reads, all repeat-element scatterplots, repeat-family summaries, assembled-contig coverage | Supplementary Fig. S5 and tables | Genome-wide burdens do not by themselves test the cis-methylation mechanism; junction-like reads occur in some diploid genomes. |

## Figure 5 source and gating rules

Run `analysis/01` before `duplicons/03`, `05` and `07`. Then run `07`, `05b`
if needed, and `08`; review `followup_report.md` and the status table. The new
renderer reads only the following outputs and does not read participant BAMs:

| Panel | Required file | Interpretation gate |
|---|---|---|
| 5a | `results/08_duplicons/followup/breakpoint_status.tsv` | Confirm both CN=1 edges and inspect carrier-specific `evidence` and `missing`; an assembled deletion junction with no homologous SD switch is **not** NAHR confirmation. Split-read support requires at least two deletion-type reads and no such reads in the comparator panel. |
| 5b | `results/08_duplicons/followup/methylation_regions.tsv` plus `results/analysis/01_evidence_matrix/common_reciprocal_cn1_core.tsv` | `raw_median_delta` is defined only where deletions expose both parental states. Show nominal permutation markers and coverage/evaluable-element numbers in supplementary tables. |
| 5c | Same `methylation_regions.tsv` | Differences versus the combined biparental tracks are descriptive; the lone mUPD is a holdout direction check, not a replicated cohort. Do not equate this with a phased maternal-versus-paternal contrast outside CN=1. |

The existing Figure 5 renderer still produces `Figure5_v7` and uses a different
table contract. The new renderer exports `Figure5_duplicon_repeat.{png,pdf,svg}`
under `results/07_figures/figure_5/`. Once real cohort tables are available,
inspect full-size PDF/PNG for label collision and inspect the analysis report
for unexpected absences before selecting the final manuscript export.

The Nextflow `MAKE_FIGURE_2` process still invokes `FIGURE_2.py` with options
that the current renderer does not accept. Run the upstream Nextflow stages
with figure generation disabled, then run the standalone analysis and figure
scripts. The new Figure 5 renderer is not wired into Nextflow.

## Checks before a Q1 resubmission

1. Confirm that all 17 rows are represented where intended. Report the actual
   evaluable n in each contrast; the eight deletion carriers define the direct
   contrast, while DiGeorge helps validate the intact-chr15 background.
2. Ensure all data-dependent quantities in the Results and captions come from
   the regenerated output tables. The repository does not contain those tables.
3. For Figure 3, verify a boundary interval and its bootstrap/threshold
   stability before using a coordinate such as a 3.7-kb core in a headline.
4. For Figure 4, recompute uncertainty at participant level or with an
   appropriate nested model. A read-level test across thousands of reads is
   descriptive unless it respects donor nesting.
5. For Figure 5, check the corrected NAHR rule, mixed-depth/missing-CpG
   sensitivity and all missing inputs. The retained chromosome in a deletion
   carrier is not the deleted transmitting chromosome; no association between
   its structural haplotype and deletion predisposition can be inferred without
   parental genomes.
6. Audit methylation BED score scales in historical scripts before reusing
   their older plotted results. The `analysis/cis_analysis/methylation.py`
   reader in this revision converts every pb-CpG-tools score from 0–100 to β.

The journal's permitted width, panel count and image requirements should be
applied to the final render; no specific journal has been selected here.
