# Figure 2 Report: Reciprocal cis-methylation architecture across the 15q11-q13 imprinted domain

Generated: 2026-06-16 12:48:46 -0500

## 1. Purpose of Figure 2
Figure 2 tests whether reciprocal deletion genomes reconstruct the parental cis-methylation architecture across the T2T chr15 imprinted domain: controls provide the biparental reference, PWS-DEL exposes the retained maternal-only profile, AS-DEL exposes the retained paternal-only profile, and PWS-mUPD confirms duplicated maternal identity.

## 2. Input files
- `phase2_run_summary.tsv`
- `tables/Phase2_control_reference_architecture.tsv`
- `tables/Phase2_retained_haplotype_profiles.tsv`
- `tables/Phase2_all_samples_haplotype_1kb_methylation.tsv.gz`
- `tables/Phase2_retained_and_mUPD_correlations.tsv`
- `tables/Phase2_reciprocal_delta_profile.tsv`
- `tables/Phase2_reciprocal_delta_boundary_candidates.tsv`
- `tables/Phase2_gene_track_features.tsv`
- `/home/rare/arlen/reference/chm13v22.sorted.gtf`
- `/home/rare/arlen/reference/ICR_t2t.bed`
- `scripts/paper_vf/create_figure2_reciprocal_cis_architecture_improved.py`

## 3. Coordinate system
- Reference: T2T-CHM13v2.0
- Analysis domain: `chr15:17,600,000-28,000,000`
- Display domain in the improved figure: `chr15:18,000,000-28,000,000`
- Major landmarks retained in the main figure: BP1, BP2, BP3, and PWS-IC/SNRPN

## 4. Window size and smoothing parameters
- Base window size: 1 kb
- Smoothed display tracks: centered rolling median, 31 windows
- Approximate smoothing span: 31 kb

## 5. Number of informative windows per panel
- Shared annotation track: structural and gene context only, no methylation windows
- Panel A (Controls): 9,658 informative windows with both control maternal and control paternal means
- Panel B (PWS-DEL): 30,845 within-footprint sample-windows across 5 deletion genomes
- Panel C (AS-DEL): 14,550 within-footprint sample-windows across 3 deletion genomes
- Panel D (PWS-mUPD): 8,518 windows for haplotype 1 and 8,488 windows for haplotype 2
- Panel E (Delta architecture): 5,313 shared windows with both PWS-DEL retained maternal and AS-DEL retained paternal means

## 6. Correlation statistics
- Controls: `r = 0.772` across 9,658 informative windows
- PWS-DEL retained maternal versus control maternal: mean `r = 0.746`, range `0.700-0.792`
- AS-DEL retained paternal versus control paternal: mean `r = 0.692`, range `0.677-0.702`
- PWS-mUPD haplotype 1 versus control maternal: `r = 0.775`
- PWS-mUPD haplotype 2 versus control maternal: `r = 0.772`
- Reciprocal overlay (PWS-DEL retained maternal versus AS-DEL retained paternal): `r = 0.691` across 5,313 shared windows

## 7. Interpretation of each panel
- Shared annotation track: provides T2T structural and gene context for BP1, BP2, BP3, the PWS imprinting center, and broader 15q11-q13 genes.
- Panel A: controls provide the biparental reference and show the reciprocal maternal and paternal methylation architecture used throughout the rest of the figure.
- Panel B: PWS-DEL retained methylation closely follows the control maternal reference, supporting exposure of a maternal-only architecture by paternal deletion.
- Panel C: AS-DEL retained methylation closely follows the control paternal reference, supporting exposure of a paternal-only architecture by maternal deletion.
- Panel D: both 004P haplotypes track the control maternal reference, consistent with copy-neutral duplicated maternal architecture rather than biparental inheritance.
- Panel E: the delta track is near zero across broad intervals, indicating a shared methylation scaffold, but shows localized parental divergence centered on the SNRPN/SNHG14 and SNORD116 interval.

## 8. BP1/BP2/BP3 are structural landmarks, not methylation-boundary calls
BP1, BP2, and BP3 are shown only as T2T structural breakpoint-cluster landmarks. They are not interpreted here as methylation boundaries. The localized methylation divergences highlighted in Panel E are data-driven parental differences within the imprinted gene domain, not calls anchored to BP1/BP2/BP3.

## 9. Main conclusion
Reciprocal PWS and AS deletion genomes reconstruct the parental cis-methylation architecture across the chr15 imprinted domain. PWS-DEL reveals the retained maternal scaffold, AS-DEL reveals the retained paternal scaffold, PWS-mUPD independently confirms duplicated maternal identity, and controls supply the biparental reference needed to interpret localized parental divergence near SNRPN/SNHG14 and SNORD116.

## 10. Output file list
- `figures/Figure2_reciprocal_cis_architecture_improved.png`
- `figures/Figure2_reciprocal_cis_architecture_improved.pdf`
- `figures/Figure2_reciprocal_cis_architecture_improved.svg`
- `reports/Figure2_reciprocal_cis_architecture_report.md`
- `tables/Figure2_shared_annotation_track.tsv`
- `tables/Figure2B_control_reference_plot.tsv`
- `tables/Figure2C_pwsdel_retained_maternal_plot.tsv`
- `tables/Figure2D_asdel_retained_paternal_plot.tsv`
- `tables/Figure2E_pwsmupd_maternal_plot.tsv`
- `tables/Figure2F_reciprocal_delta_architecture_plot.tsv`
