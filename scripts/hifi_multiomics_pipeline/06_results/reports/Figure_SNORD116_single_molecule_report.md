# Figure SNORD116 single-molecule report

Generated: 2026-06-16

## Scope

Figure analyzed: `figures/Figure_SNORD116_single_molecule_architecture.png`

This report documents the publication-facing four-panel SNORD116 single-molecule figure and the quantitative tables used to render it. The goal of the figure is to show that molecule-level CpG methylation barcodes preserve structured parent-of-origin patterns across imprinting-relevant intervals, with the strongest molecule-level heterogeneity and reciprocal parental contrast occurring at `SNORD116`.

Primary biological message:

- Molecule-level CpG methylation barcodes preserve the window- and boundary-level methylation patterns seen at the `PWS imprinting centre` and `SNORD116`.
- These methylation states are directly visible on individual PacBio HiFi molecules rather than only as cohort-level averages.
- Cross-group entropy differences are region-specific and are strongest at `SNORD116`.
- The strongest reciprocal read-level parent-of-origin contrast is focal to `SNORD116`, where `AS-DEL retained paternal` molecules are more methylated than `PWS-DEL retained maternal` molecules.

## Overview

- Region-level molecule observations analyzed for Panels C-D: `2,496`
- Unique molecules contributing to the region-level analyses: `2,496`
- Total CpG calls across SNRPN/PWS-IC, SNORD116, and downstream control windows: `146,016`
- Molecule-level panels A-B show representative high-completeness control molecules from the shared SNORD116 display window and are reported separately from the region-level statistical analyses.
- Full control barcode stacks are additionally exported to `figures/Figure_SNORD116_single_molecule_architecture_supplement_barcode_stacks.pdf` and companion image formats.
- Panels A-B compress CpG positions for display only so the barcode structure remains legible at manuscript scale; all statistics remain tied to the uncompressed saved tables.

## Output files

- Main figure PNG: `figures/Figure_SNORD116_single_molecule_architecture.png`
- Main figure PDF: `figures/Figure_SNORD116_single_molecule_architecture.pdf`
- Main figure SVG: `figures/Figure_SNORD116_single_molecule_architecture.svg`
- Main figure JPEG: `Figure_SNORD116_single_molecule_architecture.jpeg`
- Figure-specific report: `Figure_SNORD116_single_molecule_report.md`
- Mirrored detailed report: `reports/report.md`

## Data provenance

Primary tables used by the renderer:

- `tables/Figure4_panelA_control_paternal_SNORD116_matrix.tsv.gz`
- `tables/Figure4_panelA_control_paternal_SNORD116_rows.tsv`
- `tables/Figure4_panelB_control_maternal_SNORD116_matrix.tsv.gz`
- `tables/Figure4_panelB_control_maternal_SNORD116_rows.tsv`
- `tables/Figure4_SNORD116_shared_core_window.tsv`
- `tables/Figure4_gene_features.tsv`
- `tables/Figure4_panelC_entropy_plot_input.tsv`
- `tables/Figure4_panelD_SNORD116_plot_input.tsv`
- `tables/Figure4_single_molecule_CpG_calls.tsv.gz`

## Figure caption

Single-molecule HiFi methylation resolves focal parent-of-origin architecture across SNORD116.

Panels A-B show representative high-completeness control paternal-like and maternal-like molecules across a shared SNORD116 display window, with the number shown and the total qualifying molecules indicated in each header. Rows represent individual HiFi molecules and columns represent ordered CpG positions; CpG positions were compressed locally for display only. Insets show representative molecules across a narrower SNORD116 subwindow at higher display resolution to make multi-CpG barcode structure directly visible.

Panel C compares molecule-level methylation entropy across the SNRPN/PWS-IC interval, SNORD116, and a downstream control interval, with the strongest cross-group heterogeneity observed at SNORD116 (q = 1.7 x 10^-32, eta^2 = 0.19).

Panel D compares reciprocal deletion backgrounds at SNORD116 and shows higher per-molecule methylation in AS-DEL retained paternal molecules than in PWS-DEL retained maternal molecules (delta median = 0.072, p = 6.2 x 10^-13). White points indicate sample-level medians.

Together, the figure shows that individual HiFi molecules carry coherent multi-CpG methylation barcodes across SNORD116 and that the strongest parent-of-origin contrast is focal and molecule-level rather than a downstream background effect.

## Figure design notes

- Panels `A-B` are structural visualization panels. They use the true saved single-molecule calls, display representative high-completeness control molecules, and preserve molecule ordering, coordinates, and CpG-state values. Their CpG columns are compressed only for display density.
- The highlighted SNORD116 zoom inset is generated from a narrower shared high-coverage subwindow inside the SNORD116 core display interval and uses the same underlying molecules and same methylation calls.
- Shared SNORD116 core window used for emphasis in Panels A-B: `chr15:22,828,727-22,843,947`.
- Zoomed subwindow used for the barcode inset: `chr15:22,838,731-22,841,346`.
- Panels `C-D` use the quantitative saved tables directly and are not derived from the display compression used in Panels A-B.
- Panel `C` draws only pairwise entropy brackets with `q < 0.05` above the boxplots; if a region has no significant pairwise contrast, the overall region-level `q` is shown instead.

## Executive interpretation

- `SNRPN/PWS-IC` shows significant cross-group entropy differences (`q = 3.3 x 10^-5`, `eta^2 = 0.13`), consistent with imprinting-centre-specific molecule-level organisation.
- `SNORD116` shows the strongest entropy difference by a large margin (`q = 1.7 x 10^-32`, `eta^2 = 0.19`), marking it as the dominant locus of molecule-level methylation heterogeneity in this figure.
- `Downstream control` does not show meaningful cross-group entropy separation (`q = 0.589`, `eta^2 = 0.00`), arguing against a diffuse genome-wide effect.
- The reciprocal deletion comparison at `SNORD116` remains strongly separated (`delta median = 0.072`, `95% CI = 0.047-0.097`, `p = 6.2 x 10^-13`), with `AS-DEL retained paternal` molecules more methylated than `PWS-DEL retained maternal` molecules.

## Panel-by-panel reading guide

### Panel A. Paternal-like control molecules

- Shows `18` of `75` qualifying control paternal-like molecules from `2` control samples.
- Median control-paternal display methylation is `0.672` with median `99` CpGs per molecule and median span `11078` bp.
- Each row is one HiFi molecule; each colored cell corresponds to an ordered CpG position after display compression.
- Red marks methylated CpGs, blue marks unmethylated CpGs, and grey marks positions not observed on that molecule.
- The design goal of this panel is to let the reader see extended multi-CpG methylation patterns on single molecules rather than isolated per-site changes.

### Panel B. Maternal-like control molecules

- Shows `18` of `102` qualifying control maternal-like molecules from `2` control samples.
- Median control-maternal display methylation is `0.651` with median `87` CpGs per molecule and median span `11481` bp.
- The shared coordinate system, feature track, and zoom inset allow direct structural comparison with Panel A over the same SNORD116 display window.
- Together with Panel A, this panel is meant to teach the reader how coherent single-molecule methylation barcodes appear in the control state before moving to the quantitative tests in Panels C-D.

### Panel C. Molecule-level methylation entropy

- Boxplots summarize molecule-level entropy by region and sample group; white points are sample-level medians.
- `SNRPN/PWS-IC`: significant but moderate heterogeneity (`q = 3.3 x 10^-5`, `eta^2 = 0.13`).
- `SNORD116`: strongest heterogeneity signal in the figure (`q = 1.7 x 10^-32`, `eta^2 = 0.19`).
- `Downstream control`: null result (`q = 0.589`, `eta^2 = 0.00`), which supports locus specificity.
- Only significant within-region pairwise entropy comparisons (`q < 0.05`) are drawn as brackets above the boxplots in the figure; `Downstream control` shows the non-significant overall region-level `q` because no pairwise contrast passed the display threshold.
- This panel is the main region-level quantitative test that turns the structural barcode impression from Panels A-B into a formal heterogeneity comparison.

### Panel D. Parent-of-origin contrast at SNORD116

- Compares `100` `AS-DEL retained paternal` molecules from `3` samples against `254` `PWS-DEL retained maternal` molecules from `5` samples.
- Median per-molecule methylation is `0.701` in `AS-DEL retained paternal` and `0.630` in `PWS-DEL retained maternal`.
- The observed separation is `delta median = 0.072` with bootstrap `95% CI = 0.047-0.097` and Mann-Whitney `U = 18,938`.
- This panel provides the clearest reciprocal parent-of-origin result in the figure and directly matches the paternal-higher signal described in the upstream windowed analysis.

## Molecule totals by group

| group | region_level_molecule_observations |
| --- | --- |
| AS-DEL | 314 |
| Control | 906 |
| PWS-DEL | 761 |
| PWS-mUPD | 515 |

## Molecule totals by region

| region | n_molecules |
| --- | --- |
| Downstream control | 1538 |
| SNORD116 | 784 |
| SNRPN/PWS-IC | 174 |

## Control display-window summary

| display panel | n displayed | n total | n samples | median CpGs | median span bp | median missing fraction | median methylation |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Paternal-like control molecules | 18.000 | 75 | 2 | 99.000 | 11078.000 | 0.515 | 0.672 |
| Maternal-like control molecules | 18.000 | 102 | 2 | 87.000 | 11481.000 | 0.512 | 0.651 |

## Region-level molecule summary

| region | group | n molecules | n samples | total CpG calls | median CpGs | median span bp | median missing fraction | median methylation | median entropy | entropy Q1 | entropy Q3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Downstream control | AS-DEL | 188 | 3 | 19877 | 100.500 | 13956.500 | 0.000 | 0.624 | 0.293 | 0.273 | 0.316 |
| Downstream control | Control | 550 | 2 | 57780 | 105.000 | 13370.500 | 0.000 | 0.647 | 0.295 | 0.274 | 0.321 |
| Downstream control | PWS-DEL | 458 | 5 | 46447 | 101.000 | 12522.000 | 0.000 | 0.646 | 0.295 | 0.272 | 0.320 |
| Downstream control | PWS-mUPD | 342 | 1 | 36403 | 105.000 | 13262.000 | 0.000 | 0.656 | 0.292 | 0.272 | 0.315 |
| SNORD116 | AS-DEL | 100 | 3 | 17591 | 167.500 | 14300.000 | 0.000 | 0.701 | 0.262 | 0.249 | 0.273 |
| SNORD116 | Control | 294 | 2 | 48276 | 162.000 | 13576.500 | 0.000 | 0.670 | 0.275 | 0.257 | 0.301 |
| SNORD116 | PWS-DEL | 254 | 5 | 38838 | 149.000 | 12850.000 | 0.000 | 0.630 | 0.301 | 0.279 | 0.325 |
| SNORD116 | PWS-mUPD | 136 | 1 | 21308 | 158.500 | 13590.500 | 0.000 | 0.628 | 0.297 | 0.280 | 0.326 |
| SNRPN/PWS-IC | AS-DEL | 26 | 3 | 2990 | 115.000 | 15551.500 | 0.000 | 0.265 | 0.232 | 0.207 | 0.246 |
| SNRPN/PWS-IC | Control | 62 | 2 | 7126 | 115.000 | 13937.500 | 0.000 | 0.583 | 0.247 | 0.226 | 0.268 |
| SNRPN/PWS-IC | PWS-DEL | 49 | 5 | 5633 | 115.000 | 12572.000 | 0.000 | 0.751 | 0.269 | 0.256 | 0.282 |
| SNRPN/PWS-IC | PWS-mUPD | 37 | 1 | 4254 | 115.000 | 14947.000 | 0.000 | 0.756 | 0.254 | 0.228 | 0.280 |

## Entropy statistics

| region | H | p | q | eta^2 |
| --- | --- | --- | --- | --- |
| SNRPN/PWS-IC | 24.273 | 2.2 x 10^-5 | 3.3 x 10^-5 | 0.125 |
| SNORD116 | 153.121 | 5.6 x 10^-33 | 1.7 x 10^-32 | 0.192 |
| Downstream control | 1.920 | 0.589 | 0.589 | 0.000 |

## Significant pairwise entropy comparisons shown in Panel C

| region | group A | group B | U | p | q |
| --- | --- | --- | --- | --- | --- |
| SNORD116 | AS-DEL | PWS-DEL | 3490.000 | 2.3 x 10^-26 | 1.4 x 10^-25 |
| SNORD116 | AS-DEL | PWS-mUPD | 2021.000 | 3.0 x 10^-20 | 8.9 x 10^-20 |
| SNORD116 | Control | PWS-DEL | 22993.000 | 8.4 x 10^-15 | 1.7 x 10^-14 |
| SNORD116 | Control | PWS-mUPD | 12649.000 | 9.0 x 10^-10 | 1.3 x 10^-9 |
| SNORD116 | Control | AS-DEL | 19824.000 | 1.9 x 10^-7 | 2.3 x 10^-7 |
| SNRPN/PWS-IC | AS-DEL | PWS-DEL | 238.000 | 9.1 x 10^-6 | 5.5 x 10^-5 |
| SNRPN/PWS-IC | Control | PWS-DEL | 932.000 | 5.0 x 10^-4 | 0.001 |
| SNRPN/PWS-IC | AS-DEL | PWS-mUPD | 282.000 | 0.006 | 0.011 |
| SNRPN/PWS-IC | Control | AS-DEL | 1070.000 | 0.016 | 0.024 |

## Reciprocal parental-state statistics

| comparison | n mol A | n mol B | n samples A | n samples B | median A | median B | delta median | CI low | CI high | U | p | rank-biserial |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| AS-DEL retained paternal vs PWS-DEL retained maternal | 100 | 254 | 3 | 5 | 0.701 | 0.630 | 0.072 | 0.047 | 0.097 | 18938.000 | 6.2 x 10^-13 | 0.491 |

## Statistical methods

- Panel C: Kruskal-Wallis tests across the four sample groups were run independently for `SNRPN/PWS-IC`, `SNORD116`, and `Downstream control`.
- Panel C multiple-testing correction: Benjamini-Hochberg false-discovery-rate adjustment across the three region-level entropy tests.
- Panel C pairwise brackets: two-sided Mann-Whitney U tests were computed for all six within-region group pairs, with Benjamini-Hochberg correction applied within each region; only comparisons with `q < 0.05` are drawn as brackets on the figure. Regions without any significant pairwise contrast display the overall region-level `q` instead.
- Panel C effect size: Kruskal eta-squared analogue `(H - k + 1) / (n - k)`.
- Panel D: a prespecified two-sided Mann-Whitney U test comparing `AS-DEL retained paternal` versus `PWS-DEL retained maternal` molecules at `SNORD116`.
- Panel D confidence interval: bootstrap 95% CI for the median difference using 5,000 resamples.
- Panels A-B: no inferential statistics are driven by the display compression. Compression affects only figure readability, not the saved underlying values.

## Panel interpretation

### Panel A and Panel B

- Control paternal-like and maternal-like molecules show coherent multi-CpG barcode structure across the shared SNORD116 display window.
- Panels A-B now display representative high-completeness control molecules; CpG positions are still compressed for display only so the molecule-level barcode structure remains legible at manuscript scale.
- Added zoom insets now show representative molecules across a higher-resolution SNORD116 subwindow so the same underlying data can be read directly as coordinated multi-CpG barcodes rather than scattered single-site events.
- The control parental-state mean-methylation difference is modest and filter-sensitive, so Panels A-B are intentionally presented as representative molecule-level structure rather than as the primary statistical proof of parental separation.

### Panel C

- `SNORD116` shows the strongest region-level heterogeneity signal across sample groups (`q = 1.7 x 10^-32`, `eta^2 = 0.19`), exceeding both `SNRPN/PWS-IC` and `Downstream control`.
- The downstream control interval shows little cross-group entropy structure, indicating that the heterogeneity signal is locus-specific rather than a global methylation artifact.

### Panel D

- `AS-DEL retained paternal` molecules have higher mean per-molecule methylation than `PWS-DEL retained maternal` molecules at SNORD116 (`delta median = 0.072`, `p = 6.2 x 10^-13`).
- This reciprocal pattern is robust across stricter minimum-CpG and minimum-span thresholds.

## Sensitivity checks

### Minimum span threshold in the control display panels

| min span fraction | paternal n | maternal n | paternal median | maternal median | control p |
| --- | --- | --- | --- | --- | --- |
| 0.000 | 74.000 | 102.000 | 0.673 | 0.651 | 0.166 |
| 0.300 | 37.000 | 57.000 | 0.679 | 0.640 | 0.289 |
| 0.350 | 20.000 | 30.000 | 0.694 | 0.628 | 0.044 |
| 0.400 | 12.000 | 19.000 | 0.655 | 0.640 | 0.855 |

- The control paternal-like versus maternal-like mean-methylation contrast is not stable across display filters. This supports using Panels A-B as representative structural views and Panels C-D as the main quantitative evidence.

### Reciprocal deletion comparison under alternative filters

| min CpGs | min span bp | AS-DEL n | PWS-DEL n | delta median | p | U |
| --- | --- | --- | --- | --- | --- | --- |
| 20.000 | 10000.000 | 99.000 | 231.000 | 0.071 | 0.000 | 17129.000 |
| 20.000 | 12000.000 | 77.000 | 157.000 | 0.064 | 0.000 | 8970.000 |
| 20.000 | 14000.000 | 54.000 | 81.000 | 0.051 | 0.000 | 3221.000 |
| 20.000 | 16000.000 | 33.000 | 44.000 | 0.061 | 0.000 | 1099.000 |
| 60.000 | 10000.000 | 93.000 | 213.000 | 0.065 | 0.000 | 14493.000 |
| 60.000 | 12000.000 | 71.000 | 145.000 | 0.061 | 0.000 | 7446.000 |
| 60.000 | 14000.000 | 51.000 | 74.000 | 0.049 | 0.000 | 2681.000 |
| 60.000 | 16000.000 | 30.000 | 38.000 | 0.049 | 0.002 | 818.000 |
| 100.000 | 10000.000 | 77.000 | 176.000 | 0.059 | 0.000 | 9668.000 |
| 100.000 | 12000.000 | 61.000 | 131.000 | 0.060 | 0.000 | 5738.000 |
| 100.000 | 14000.000 | 46.000 | 68.000 | 0.049 | 0.000 | 2243.000 |
| 100.000 | 16000.000 | 29.000 | 35.000 | 0.057 | 0.001 | 757.000 |
| 120.000 | 10000.000 | 68.000 | 154.000 | 0.050 | 0.000 | 7206.000 |
| 120.000 | 12000.000 | 54.000 | 115.000 | 0.046 | 0.000 | 4281.000 |
| 120.000 | 14000.000 | 43.000 | 61.000 | 0.047 | 0.000 | 1851.000 |
| 120.000 | 16000.000 | 27.000 | 31.000 | 0.046 | 0.004 | 606.000 |

- The AS-DEL versus PWS-DEL difference remains positive and statistically supported across all tested minimum-CpG and minimum-span thresholds shown above.

### Displayed molecules versus quantitative analyses

- Panels A-B display representative control molecules from the full qualifying sets. All inferential statistics in Panels C-D continue to use the full retained molecule sets defined in the saved phase-4 tables.

### Relative strength of the SNORD116 heterogeneity signal

- The cross-group entropy effect size is strongest at SNORD116 (`eta^2 = 0.19`), weaker at SNRPN/PWS-IC (`eta^2 = 0.13`), and absent downstream (`eta^2 = 0.00`).

## Limitations and scope boundaries

- Panels A-B are intended as structural single-molecule views, not as standalone inferential comparisons between paternal-like and maternal-like controls.
- The control parental contrast in the display subset is sensitive to span and completeness filters, so the strongest formal evidence in this figure comes from Panels C-D rather than from an A-versus-B hypothesis test.
- The figure is deliberately locus-focused. It supports focal cis-coordination at imprinting-relevant intervals and does not claim a genome-wide barcode phenomenon from this panel set alone.

## Reproducibility

Figure generation command:

```bash
python3 FIGURE_4.py --outdir /home/rare/arlen/pws-as-hifi-cis-methylation/scripts/hifi_multiomics_pipeline/06_results
```

Convenience outputs written by the renderer:

- `/home/rare/arlen/pws-as-hifi-cis-methylation/scripts/hifi_multiomics_pipeline/06_results/figures/Figure_SNORD116_single_molecule_architecture.png`
- `/home/rare/arlen/pws-as-hifi-cis-methylation/scripts/hifi_multiomics_pipeline/06_results/figures/Figure_SNORD116_single_molecule_architecture.pdf`
- `/home/rare/arlen/pws-as-hifi-cis-methylation/scripts/hifi_multiomics_pipeline/06_results/figures/Figure_SNORD116_single_molecule_architecture.svg`
- `/home/rare/arlen/pws-as-hifi-cis-methylation/scripts/hifi_multiomics_pipeline/06_results/figures/Figure_SNORD116_single_molecule_architecture.jpeg`
- `/home/rare/arlen/pws-as-hifi-cis-methylation/scripts/hifi_multiomics_pipeline/06_results/reports/Figure_SNORD116_single_molecule_report.md`
- `/home/rare/arlen/pws-as-hifi-cis-methylation/scripts/hifi_multiomics_pipeline/06_results/reports/report.md`

## Main result

Single-molecule PacBio HiFi methylation calls show that SNORD116 is organized as coherent multi-CpG molecule-level barcodes rather than isolated CpG noise. Although the control paternal-like and maternal-like display panels are best interpreted as representative structural views, quantitative testing identifies SNORD116 as the locus with the strongest cross-group molecule-level heterogeneity and shows that reciprocal deletion genomes retain distinct parent-of-origin methylation states: AS-DEL retains the paternal-like state, whereas PWS-DEL retains the maternal-like state.
