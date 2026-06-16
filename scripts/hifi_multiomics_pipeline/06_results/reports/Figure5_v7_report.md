# Figure 5 v7 report

## Revised plotting strategy

- Build the main figure around a single message: chr15q11-q13 deletion architecture is clear, whereas genome-wide CNV/SV burden and breakpoint-aligned methylation effects are modest.
- Move dense full-sample chr15 coverage tracks out of the main figure and into a supplementary companion panel.
- Keep only sample-level burden summaries that are interpretable at a glance: chr15 deletion size, non-chr15 large CNV burden, total SV count, and total SV span.
- Treat PWS-UPD as a breakpoint-coordinate-matched descriptive reference rather than as a true breakpoint-flanking deletion analysis.

## Figure layout

- Panel A: chr15 HiFi coverage and deletion classes using the v1-style coverage track layout.
- Panel B: single Manhattan-style CNV panel with chr15 highlighted and compact non-chr15 burden statistics.
- Panel C: total SV count and total SV span per sample, with p-values displayed in a separate header band above the plots.
- Panel D: breakpoint-aligned methylation difference versus controls across 0-10, 10-25, 25-50, and 50-100 kb bins.
- Panel E: forest-style near-versus-far breakpoint methylation effect summary.
- Supplementary coverage figure: full chr15 HiFi coverage tracks for all deletion carriers.

## Panel-specific recommendations implemented

- Panel A retains the coverage-track view but uses cleaner left-side labels and spacing.
- Panel B now uses one condensed Manhattan-style panel with group-colored chr15 deletion points, separated sample callouts, and the chr15 label above the plotting area.
- Panel C now keeps only total count and total span, and moves p-values into a dedicated header band rather than inside the plotting panels.
- Panels D-E use conservative language and emphasize effect sizes, confidence intervals, and null-crossing intervals rather than implying strong distance-decay.
- Multiple-testing correction is applied across SV metrics and across the four formal methylation near-versus-far tests.

## Revised script structure

1. Load existing Figure 5 input tables and methylation file inventory.
2. Recompute sample-level non-chr15 CNV burden and SV burden statistics.
3. Re-bin the existing breakpoint-coordinate-aligned methylation table into distance-decay intervals.
4. Build distance-bin summaries, exact permutation/sign-flip statistics, and bootstrap confidence intervals.
5. Render a simplified main figure plus a supplementary full-coverage figure and write analysis tables and narrative report.

## Suggested panel titles

- A. chr15 HiFi coverage and deletion classes
- B. chr15 deletion dominates large CNV signal
- C. Global SV burden does not clearly separate diagnostic groups
- D. Limited evidence for breakpoint-proximal methylation decay
- E. Near-versus-far breakpoint methylation effects are small

## Key quantitative observations

- Non-chr15 CNV count burden remains weakly separated across groups (`exact permutation p=0.337`, `q=0.438`, `epsilon^2=0.07`).
- Non-chr15 CNV total span is similarly non-dominant (`exact permutation p=0.438`).
- Total SV count does not separate groups strongly (`p=0.292`, `q=0.340`, `epsilon^2=0.12`), and the same is true for total SV span (`p=0.267`).
- `BND` burden is the only nominal SV signal (`p=0.013`), but it does not remain significant after FDR (`q=0.090`).
- The deletion architecture remains dominated by recurrent BP1/BP2-to-BP3/BP4-like events with `007P` as the only clearly atypical extended deletion.
- Methylation interpretation is intentionally conservative because parental-haplotype assignments are incomplete for most deletion carriers; the v7 analysis uses retained haplotype where labeled and combined methylation otherwise.

## Suggested caption text

Figure 5. Structural deletion architecture, genome-wide structural burden, and breakpoint-associated methylation. (A) chr15 deletion intervals are shown schematically per sample with canonical BP1-BP5 guides, separating recurrent BP1/BP2-to-BP3/BP4-like classes from the single atypical extended deletion. (B) chr15 deletion size is shown per sample together with non-chr15 large autosomal CNV burden (`>=2 Mb`); canonical chr15 deletions are the dominant CNV events, whereas non-chr15 burden overlaps across groups. (C) Global SV burden is summarized as total SV count, total SV span, and a compact type-specific statistics table; no SV burden metric survives FDR correction. (D) Breakpoint-aligned methylation is summarized as signed delta methylation relative to matched controls across distance bins from `0-10 kb` to `50-100 kb`. (E) Near-versus-far breakpoint methylation effects are summarized per group and breakpoint side. PWS-UPD is shown as a breakpoint-coordinate-matched descriptive reference rather than a true breakpoint-flanking deletion analysis. A supplementary panel provides the full chr15 HiFi coverage tracks for all deletion carriers. Across panels, the figure supports a recurrent chr15 structural mechanism with limited evidence that global genome-wide CNV/SV burden or broad breakpoint-flanking methylation change is the primary discriminating signal.

## Results-ready interpretation template

Deletion carriers showed a predominantly recurrent chr15 architecture, with most samples mapping to BP1/BP2-to-BP3/BP4-like classes and a single atypical extended deletion. Outside the canonical chr15 event, large autosomal CNV burden did not separate groups strongly (`exact permutation p=0.337` for count burden), and genome-wide SV burden showed similarly shallow differences (`total SV count p=0.292`; all SV burden `q>=0.05`). Breakpoint-coordinate-aligned methylation differences relative to controls remained small overall and were most consistent with weak, local deviations rather than a broad or uniform epigenetic bleed effect. The UPD sample was analyzed only at canonical breakpoint-matched coordinates and is therefore interpreted descriptively rather than as evidence for true breakpoint-flanking methylation change.

## Methylation effect summary

- PWS-DEL 5': near `|Δ|=0.027`, far `|Δ|=0.008`, near-minus-far `+0.022`, p=0.250, q=0.500; Formal exact sign-flip test across samples.
- PWS-DEL 3': near `|Δ|=0.026`, far `|Δ|=0.016`, near-minus-far `+0.010`, p=0.625, q=0.833; Formal exact sign-flip test across samples.
- AS-DEL 5': near `|Δ|=0.008`, far `|Δ|=0.009`, near-minus-far `-0.001`, p=1.000, q=1.000; Formal exact sign-flip test across samples.
- AS-DEL 3': near `|Δ|=0.031`, far `|Δ|=0.012`, near-minus-far `+0.019`, p=0.250, q=0.500; Formal exact sign-flip test across samples.
- PWS-UPD BP-matched 5': near `|Δ|=0.023`, far `|Δ|=0.020`, near-minus-far `+0.003`, p=n/a, q=n/a; Descriptive only; BP-coordinate-matched regions in one PWS-UPD sample.
- PWS-UPD BP-matched 3': near `|Δ|=0.009`, far `|Δ|=0.004`, near-minus-far `+0.005`, p=n/a, q=n/a; Descriptive only; BP-coordinate-matched regions in one PWS-UPD sample.
