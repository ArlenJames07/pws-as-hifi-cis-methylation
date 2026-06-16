# Figure 1 Report: Diagnostic validation at the PWS-AS imprinting centre

## 1. Purpose
Figure 1 tests whether haplotype-resolved PacBio HiFi methylation recovers the expected diagnostic parent-of-origin states at the PWS-AS imprinting centre across PWS-DEL, PWS-mUPD, AS-DEL, and control samples.

## 2. Input data
- Allele-level methylation matrix: `tables/Figure1A_allele_methylation_matrix.tsv`
- Per-CpG contrast table: `tables/Figure1B_per_CpG_contrast.tsv`
- Diagnostic state summary: `tables/Figure1C_diagnostic_state_summary.tsv`
- Coverage/phasing support table: `tables/Figure1D_coverage_phasing_support.tsv`
- Metadata / parameters: `/home/rare/arlen/outputs/methylation/metadata/metadata_methylation.csv` and `/home/rare/arlen/pws-as-hifi-cis-methylation/scripts/hifi_multiomics_pipeline/06_results/figures/phase1_run_parameters.json`
- Script: `FIGURE_1.py`

## 3. Coordinate system
- Reference: T2T-CHM13v2.0
- PWS-AS IC core interval: `chr15:22,691,258-22,693,494`
- No hg38/GRCh38 coordinates were used.

## 4. Panel A interpretation
- PWS-DEL samples retain one methylated maternal-pattern allele and lack the paternal allele.
- PWS-mUPD retains two methylated maternal-pattern haplotypes.
- AS-DEL samples retain one unmethylated paternal-pattern allele and lack the maternal allele.
- Controls retain one maternal-pattern and one paternal-pattern allele.

| Group | Expected state | Observed state | Interpretation |
| --- | --- | --- | --- |
| PWS-DEL | M / absent | M / absent | paternal deletion |
| PWS-mUPD | M / M | M / M | maternal UPD (duplicated maternal state) |
| AS-DEL | absent / P | absent / P | maternal deletion |
| Control | M / P | M / P | canonical biparental state |

## 5. Panel B interpretation
- PWS-DEL shows a positive maternal-pattern contrast signal across the IC.
- AS-DEL shows a negative paternal-pattern contrast signal.
- Controls show the canonical maternal-minus-paternal contrast.
- PWS-mUPD shows near-zero contrast because both retained haplotypes are maternal-pattern.
- Near-zero contrast in PWS-mUPD does not indicate absence of methylation signal.

## 6. Panel C interpretation
- PWS-DEL = M / absent
- PWS-mUPD = M / M
- AS-DEL = absent / P
- Control = M / P

## 7. Panel D interpretation
- Panel D summarizes total IC depth, minimum allele depth, CpGs per allele, and phased span across the IC.
- Filled circles indicate calls passing the nominal IC support threshold.
- Open circles indicate calls below the nominal IC support threshold.
- Deletion samples can show lower apparent allele support because the affected interval is biologically hemizygous.

## 8. Main conclusion
Together, these results validate that allele-resolved long-read methylation, interpreted with copy-number and phasing support, recovers the expected molecular configurations at the PWS-AS imprinting centre. This diagnostic validation supports downstream reconstruction of parental cis-methylation architecture across 15q11-q13.

## 9. Figure caption draft
Figure 1. Diagnostic validation at the PWS-AS imprinting centre. (A) Mean methylation across the canonical PWS-AS imprinting centre core (`chr15:22,691,258-22,693,494`, T2T-CHM13v2.0) shown for two physical alleles per sample. PWS-DEL samples retain one maternal-pattern methylated allele and lack the paternal allele, PWS-mUPD retains two maternal-pattern physical haplotypes, AS-DEL retains one paternal-pattern unmethylated allele and lacks the maternal allele, and controls show the canonical maternal/paternal biparental state. Hatched cells indicate absent/deleted alleles. (B) Per-CpG parent-of-origin methylation contrast across the IC, shown as group-level median profiles with interquartile ribbons and faint individual-sample traces. PWS-DEL retains positive maternal-pattern signal, AS-DEL retains negative paternal-pattern signal, controls show the canonical maternal-minus-paternal profile, and PWS-mUPD remains near zero because both physical haplotypes are maternal-pattern. (C) Compact diagnostic summary of expected and observed allele configurations by mechanism. (D) Coverage and phasing support at the IC, summarizing total IC depth, minimum allele depth, CpGs per allele, and phased IC span. M, maternal-pattern methylation; P, paternal-pattern methylation; absent, deleted/absent allele.

## 10. Output files
- `figures/Figure1_improved.png`
- `figures/Figure1_improved.pdf`
- `figures/Figure1_improved.svg`
- `reports/Figure1_report.md`
- `tables/Figure1A_allele_methylation_matrix.tsv`
- `tables/Figure1B_per_CpG_contrast.tsv`
- `tables/Figure1C_diagnostic_state_summary.tsv`
- `tables/Figure1D_coverage_phasing_support.tsv`

## 11. Quality-control checks
- The PWS-mUPD sample is shown as two maternal-pattern physical haplotypes and is not collapsed into one maternal cell.
- Absent alleles are represented by hatching, not by grey fill alone.
- Panel B uses an external legend so group labels do not overlap the main traces.
- Panels A and D use the same anonymized sample order: `PW-1` to `PW-5`, `UPD-1`, `AS-1` to `AS-3`, `CTRL-1` to `CTRL-2`.
- All coordinates shown are T2T-CHM13v2.0 coordinates.
- This Markdown report was generated automatically.
