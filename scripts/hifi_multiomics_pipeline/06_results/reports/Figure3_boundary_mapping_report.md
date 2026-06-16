# Figure 3 boundary mapping report

## Input files used
| Sample | Layer | Path |
| --- | --- | --- |
| 001P | combined | /home/rare/arlen/outputs/methylation/genomes_2/08_1_A01_bc2043_001P.combined.bed |
| 002P | combined | /home/rare/arlen/outputs/methylation/genomes_2/08_1_A01_bc2044_002P.combined.bed |
| 005P | combined | /home/rare/arlen/outputs/methylation/genomes_2/08_1_A01_bc2047_005P.combined.bed |
| 006P | combined | /home/rare/arlen/outputs/methylation/genomes_2/08_1_B01_bc2048_006P.combined.bed |
| 007P | combined | /home/rare/arlen/outputs/methylation/genomes_2/08_1_C01_bc2049_007P.combined.bed |
| 013A | combined | /home/rare/arlen/outputs/methylation/genomes_2/08_1_A01_bc2055_013A.combined.bed |
| 014A | combined | /home/rare/arlen/outputs/methylation/genomes_2/08_1_A01_bc2056_014A.combined.bed |
| 016A | combined | /home/rare/arlen/outputs/methylation/genomes_2/08_1_B01_bc2058_016A.combined.bed |
| 017C | combined | /home/rare/arlen/outputs/methylation/genomes_2/08_1_C01_bc2059_017C.combined.bed |
| 017C | hap1 | /home/rare/arlen/outputs/methylation/genomes_2/08_1_C01_bc2059_017C.hap1.bed |
| 017C | hap2 | /home/rare/arlen/outputs/methylation/genomes_2/08_1_C01_bc2059_017C.hap2.bed |
| 018C | combined | /home/rare/arlen/outputs/methylation/genomes_2/08_1_D01_bc2060_018C.combined.bed |
| 018C | hap1 | /home/rare/arlen/outputs/methylation/genomes_2/08_1_D01_bc2060_018C.hap1.bed |
| 018C | hap2 | /home/rare/arlen/outputs/methylation/genomes_2/08_1_D01_bc2060_018C.hap2.bed |

## Coordinate system used
- Coordinate system: `T2T-CHM13 chr15 coordinate (Mb)`
- Control signal mode: `absolute`
- Control contrast formula: `|maternal - paternal|`
- PWS-DEL contrast formula: `retained maternal - control biallelic baseline`
- AS-DEL contrast formula: `control biallelic baseline - retained paternal`

## Interval definitions
| Interval | Coordinates | Width (bp) | Width (kb) | Shared core fraction |
| --- | --- | --- | --- | --- |
| Controls | chr15:22,690,600-22,696,850 | 6,250 | 6.25 | 59.2% |
| PWS-DEL | chr15:22,690,900-22,694,700 | 3,800 | 3.80 | 97.4% |
| AS-DEL | chr15:22,691,000-22,694,800 | 3,800 | 3.80 | 97.4% |
| Shared core | chr15:22,691,000-22,694,700 | 3,700 | 3.70 | 100% |

## Pairwise overlaps
| Comparison | Overlap (bp) | Overlap (kb) | Fraction of first interval | Fraction of second interval | Jaccard |
| --- | --- | --- | --- | --- | --- |
| Controls ∩ PWS-DEL | 3,800 | 3.80 | 60.8% | 100% | 0.608 |
| Controls ∩ AS-DEL | 3,800 | 3.80 | 60.8% | 100% | 0.608 |
| PWS-DEL ∩ AS-DEL | 3,700 | 3.70 | 97.4% | 97.4% | 0.949 |
| Controls ∩ PWS-DEL ∩ AS-DEL | 3,700 | 3.70 | 59.2% | 97.4% |  |

## Shared core fractions
- Shared core width = `3,700 bp` = `3.70 kb`
- Controls width = `6,250 bp` = `6.25 kb`
- PWS-DEL width = `3,800 bp` = `3.80 kb`
- AS-DEL width = `3,800 bp` = `3.80 kb`
- Shared core / Controls = `59.2%`
- Shared core / PWS-DEL = `97.4%`
- Shared core / AS-DEL = `97.4%`

## Breakpoint distances
| Breakpoint | Interval | Distance to shared start (bp) | Distance to shared midpoint (bp) | Distance to shared end (bp) |
| --- | --- | --- | --- | --- |
| BP1 | chr15:17,691,439-20,454,275 | 2,236,725 | 2,238,575 | 2,240,425 |
| BP2 | chr15:20,753,698-21,183,655 | 1,507,345 | 1,509,195 | 1,511,045 |
| BP3 | chr15:25,875,912-26,632,507 | 3,184,912 | 3,183,062 | 3,181,212 |

## Sensitivity analysis
| Method | Called interval | Width (bp) | Width (kb) | Shared-core overlap (bp) | Shared-core overlap (%) | Called-interval overlap (%) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 250-bp windows | chr15:22,691,400-22,694,350 | 2,950 | 2.95 | 2,950 | 79.7% | 100% | partial overlap |
| 500-bp windows | chr15:22,691,300-22,694,400 | 3,100 | 3.10 | 3,100 | 83.8% | 100% | partial overlap |
| 1-kb windows | chr15:22,691,000-22,694,700 | 3,700 | 3.70 | 3,700 | 100% | 100% | matches shared core |
| CpG-adaptive windows | chr15:22,691,609-22,695,739 | 4,130 | 4.13 | 3,091 | 83.5% | 74.8% | partial overlap |
| Alternative threshold | chr15:22,691,000-22,694,600 | 3,600 | 3.60 | 3,600 | 97.3% | 100% | partial overlap |
| Change-point detection | chr15:22,691,200-22,693,900 | 2,700 | 2.70 | 2,700 | 73% | 100% | partial overlap |
| Bootstrap resampling | chr15:22,692,300-22,693,848 | 1,548 | 1.55 | 1,548 | 41.8% | 100% | partial overlap |

## Warnings
- At least one primary sliding window had zero CpGs.
- LOESS smoothing was not added because the current workflow has no dedicated LOESS-based boundary caller; omitted rather than fabricated.

## Short interpretation
The plotted Figure 3 intervals support a focal regulatory boundary rather than a deletion-wide methylation effect. Controls, PWS-DEL, and AS-DEL converge on the same shared 3.7-kb core at `chr15:22,691,000-22,694,700`, while the structural BP1/BP2/BP3 deletion architecture remains distinct from that focal methylation-transition interval. Across available windowing and boundary-calling sensitivity analyses, the called intervals remain centred on the shared core, supporting positional convergence rather than amplitude-based interpretation.
