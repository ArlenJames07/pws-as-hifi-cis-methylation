# 15q11–q13 duplicons, breakpoints and repeats

Python scripts that add one new analysis to the repository without touching the
existing workflow. Everything lives in `scripts/duplicons/`; the only other new
files are `tests/test_duplicons.py` and your local `assemblies.local.csv` (and,
if needed, `bams.local.csv`).

## Where it fits

```
pws-as-hifi-cis-methylation/
├── main.nf, modules/, envs/ ...        unchanged — stages 01–07 as today
├── scripts/
│   ├── analysis/                        unchanged — 00–03 and cis_analysis (reused, not modified)
│   ├── figures/                         unchanged
│   └── duplicons/                       NEW
│       ├── config.py                    every path, region and threshold (reads params.local.yml)
│       ├── 00_check_inputs.py           checks BAMs, assemblies, tracks, deletions before 01–07
│       ├── 01_prepare_reference.py      CHM13 window, SD pairs, SUNK set, RepeatMasker on the window
│       ├── 02_count_sunks.py            read k-mer counts at the SUNKs, per participant
│       ├── 03_place_contigs.py          hifiasm contigs → CHM13, 15q extraction, haplotype RepeatMasker
│       ├── 04_assembly_methylation.py   optional: IC-anchored label for the contig carrying the IC
│       ├── 05_duplicon_breakpoints.py   breakpoints inside the duplicons: copy number + junction reads
│       ├── 05b_breakpoint_summary.py    rebuilds 05's per-carrier tables in seconds
│       ├── 06_haplotype_repeats.py      repeats per assembled haplotype, repeat content of indels
│       ├── 07_repeat_methylation.py     parent-of-origin methylation of repeat elements
│       ├── 08_followup.py               split reads, assembled fusions, mUPD check, genes; one status table
│       ├── run_duplicons.py             runs 01–07 in order
│       ├── assemblies.example.csv
│       ├── bams.example.csv
│       └── duplicon_analysis/           shared code (PAF projection, k-mers, readers, tools)
├── tests/test_duplicons.py              NEW
├── assemblies.local.csv                 NEW, local only (add it to .gitignore)
├── bams.local.csv                       NEW, optional, local only
└── results/08_duplicons/                NEW outputs: reference/ sunk_counts/ assembly/
                                         breakpoints/ haplotype_repeats/ repeat_methylation/
```

It reads what the workflow already produces:

| Input | Used by |
|---|---|
| `params.local.yml` (reference, cpg_model, *_bin paths) | all |
| `assets/metadata.csv` | all (via `cis_analysis.load_cohort`) |
| HiFi BAMs (+ .bai) aligned to the `reference` in params.local.yml; MM/ML tags for 04 — see "BAM files" | 02, 04, 05 |
| `results/04_phasing/<sample>.blocks.tsv` | 07 (optional IC-anchored table) |
| `results/06_methylation/<sample>/<sample>.cpg.*.bed` | 07 |
| `results/analysis/01_evidence_matrix/chr15_structural_evidence.tsv` | 03, 05, 07 (via `load_deletion_map`) |
| hifiasm contigs listed in `assemblies.local.csv` | 03, 04 |

## Setup

1. Copy `scripts/duplicons/` and `tests/test_duplicons.py` into the repository.
2. `cp scripts/duplicons/assemblies.example.csv assemblies.local.csv` and fill in the
   paths to the hifiasm contigs (GFA or FASTA, gzipped or not; contigs, not RagTag
   scaffolds). Columns `sample,hap1,hap2,primary`; leave a cell empty if you don't
   use that assembly. Add `assemblies.local.csv` to `.gitignore`.

   | Participants | Columns to fill | Why |
   |---|---|---|
   | Control, DiGeorge, PWS-mUPD | `hap1`, `hap2` (`*.bp.hap1/hap2.p_ctg.gfa`) | two chromosomes 15 to compare |
   | PWS-DEL, AS-DEL | `primary` (`*.bp.p_ctg.gfa`), plus `hap1`/`hap2` if you want to compare | inside the deletion there is one chromosome; hifiasm may drop, fragment or duplicate it in hap1/hap2, but it is always in p_ctg |

   After step 03, `<sample>.window_cover.tsv` reports for every assembly how much of
   the HiFiCNV CN=1 interval it carries (`cn1_deletion_interval`) and whether it is
   present more than once (`copies_in_assembly`). Use it to confirm which assembly
   holds the retained chromosome. Inside that interval the primary contigs are the
   single retained chromosome (maternal in PWS-DEL, paternal in AS-DEL); outside it,
   p_ctg mixes both parental chromosomes, so compare flanks with hap1/hap2 instead.
3. BAM files: see "BAM files" below.
4. Tools: `minimap2`, `samtools`, `meryl`, `RepeatMasker` (with TRF), and
   `pb-CpG-tools` for step 04 only. Put their paths in `params.local.yml` if they are
   not on `PATH`: `minimap2_bin`, `meryl_bin`, `repeatmasker_bin`
   (`samtools_bin` and `pbcpgtools_bin` are already there). Python needs numpy,
   pandas and matplotlib (the `envs/figures.yml` environment has them).
   RepeatMasker 4.2 (bioconda) ships without a repeat library: run `download_dfam.py`
   once and choose `1,2` (Dfam root + curated consensus), then check
   `famdb.py check 'Homo sapiens'` reports Curated Consensus as `present`. If
   `~/.local` holds another pandas/numpy, run
   `conda env config vars set PYTHONNOUSERSITE=1` in the environment.
5. Check the regions in `config.py` (see "Checks" below). Threads and meryl memory:
   `DUPLICON_THREADS`, `DUPLICON_MERYL_GB` environment variables.

## BAM files

The scripts look for each participant's BAM in `results/01_alignment/` (another
folder: `DUPLICON_ALIGNMENT_DIR=/path`), in this order:

1. `bams.local.csv` at the repository root (`sample,bam`; several BAMs of one
   participant separated by `;`) — always wins;
2. `<sample>.aligned.bam` (the workflow's name);
3. any `*.bam` whose name contains the sample ID as a separate token, e.g. the
   SMRT Link name `08_1_A01_bc2043_001P.bam` (`x_1005P.bam` does not match `005P`).
   Links are followed. Two matches (e.g. `…_bc2044_002P.bam` and
   `…_bc2044v2_002P.bam`) stop the run until you choose in `bams.local.csv`:

```
sample,bam
002P,/path/to/results/01_alignment/08_1_A01_bc2044v2_002P.bam
```

List two BAMs for one participant only if they hold different reads (two
sequencing runs); the same reads twice would double its k-mer counts.

Step 02 reads every primary read whatever the reference. Step 05 looks reads up by
CHM13 coordinates, so each BAM must be aligned to the same reference as
`params.local.yml`: the scripts compare the chr15 length in the BAM header with the
reference `.fai` (CHM13 99,753,195 bp; GRCh38 101,991,189 bp) and skip junction
reads for a BAM that differs. Step 04 re-maps the reads itself but needs MM/ML tags.
`00_check_inputs.py` reports all of this per participant.

Methylation tracks are found as `<sample>.cpg.<kind>.bed[.gz]` or, failing that, one
pb-CpG-tools file under `results/06_methylation/` named `*<sample>*.<kind>.bed[.gz]`.

## Run

```bash
python3 scripts/duplicons/00_check_inputs.py          # inputs only; writes results/08_duplicons/input_check.tsv
python3 scripts/duplicons/run_duplicons.py            # 00 (stops on input problems), 01, 02, 03, 05, 06, 07
python3 scripts/duplicons/run_duplicons.py --with-04  # also the optional contig labels
python3 scripts/duplicons/run_duplicons.py --from 05  # resume (skips 00)
python3 scripts/duplicons/02_count_sunks.py --samples 001P,002P   # any step alone
python3 tests/test_duplicons.py -v
```

Steps 01–04 skip outputs that already exist (`--force` recomputes). Run analysis
`01_build_chr15_evidence_matrix.py` first: 03, 05 and 07 read its CN=1 intervals.

### Step 08: one table with what is confirmed and what is missing

`08_followup.py` reads 03, 05/05b and 07 and adds split reads (SA tags linking the two
edges, counted in the carrier and in every genome with two copies of chr15 at the same
positions), assembled fusions (with the nearest confidently placed contig blocks on each
side) and genes from the GTF in params.local.yml. Verdicts in `breakpoint_status.tsv`:

| status | meaning |
|---|---|
| junction resolved at bp level | >= 2 deletion-type split reads in the carrier, none in the panel |
| confirmed NAHR | edges at homologous positions of one direct pair plus junction-read excess at the carrier's crossover; or an assembled contig that switches between the two copies of a direct pair within 20 kb of homologous positions. A contig crossing the deletion without that homologous switch supports the deletion, not NAHR. |
| confirmed by assembly | a contig spans the deletion, but its switch is not at homologous positions of a direct pair |
| compatible with NAHR | edges homologous, no independent confirmation yet |
| edges near one SD pair, not homologous | same pair, > 20 kb from homologous positions |
| edge outside WINDOW | CN 1 reaches the end of the window |
| unresolved | no pair links the edges |

Assembled junctions are placed more precisely than SUNK edges (the contig switch is within
the last/first confidently placed block, usually < 1 kb); `crossover_estimate` gives the best
available position: split reads (bp) > assembly > SUNK edges (several kb).

`methylation_domains.tsv` reports runs of 250-kb bins with permutation p <= 0.05 and a raw
difference of at least 0.02 (smaller shifts are of the size of the centring offset). Each
bin and each domain also get three contrasts against the biparental genomes (combined
methylation): PWS-DEL (maternal copy only), AS-DEL (paternal copy only) and PWS-mUPD (two
maternal copies). If the paternal chromosome is more methylated by d inside a domain, PWS-DEL
and PWS-mUPD drop by about d/2 relative to the bins outside it and AS-DEL rises by d/2
(`*_inside_minus_outside`). PWS-mUPD takes no part in finding the domains, so it is the
independent check. It needs 07 from v6 on, which stores the PWS-mUPD combined values
(`scaffold_combined` in `element_by_participant.tsv.gz`); rerun 07 (minutes), then 08.

### What needs recomputing after a change

| Change | Rerun |
|---|---|
| `DUPLICON_DOMAIN_START/END`, thresholds in 07 | `07_repeat_methylation.py` only (minutes), then 08 |
| a new version of 08 or its GTF | `08_followup.py` only (`--skip-split-reads` if they were already counted: seconds) |
| parameters of 05 (bin size, `P_SWITCH`), panel | `05_duplicon_breakpoints.py` only (minutes) |
| a new or corrected BAM for one participant | `02 --samples X --force`, then 05 |
| `WINDOW` | 01–03 with `--force`; 01–04 skip existing outputs, so without `--force` the old window is silently kept |

To look at a region outside `WINDOW` (e.g. the proximal edge of one deletion) without
redoing the main run, analyse a small window in its own folder and count only the reads
aligned there (`--window-reads`: minutes per genome instead of hours):

```bash
export DUPLICON_OUT=results/08_duplicons_extra DUPLICON_WINDOW=chr15:16000000-18500000
python3 scripts/duplicons/01_prepare_reference.py
python3 scripts/duplicons/02_count_sunks.py --window-reads --samples <panel genomes>,<carrier>
python3 scripts/duplicons/05_duplicon_breakpoints.py
```

Use `--window-reads` for every genome of that run (panel included): counts from reads
aligned to the window and full-read counts must not be mixed in one run.

## What each step answers, and what to read first

| Step | Question | Read first |
|---|---|---|
| 05 | Where inside the BP duplicons did each deletion break? Is it NAHR? | `junction_crossovers.tsv`, then `sunk15q.overview.png` |
| 03 | How much of BP1–BP5 did each haplotype assemble? Inversions? | `assembly/<sample>/<sample>.window_cover.tsv`, `*.junctions.tsv` |
| 06 | Which repeats do the haplotypes carry, and how do they differ from CHM13? | `haplotype_repeats/cohort_sv_repeats.tsv` |
| 07 | Is the methylation of repeat elements parent-of-origin specific? | `repeat_methylation/validation.tsv` |
| 08 | For each deletion: confirmed, compatible or unresolved, and what is missing? Which regions differ by parent, and which genes are there? | `followup/breakpoint_status.tsv`, `followup/followup_report.md` |

**What the real cohort showed (read this before 05's outputs).** Reads with copy-A SUNKs
followed by copy-B SUNKs ("junction-like" reads) occur in the genomes with two copies of
chr15 at similar rates to the carriers: paralog-specific variants of the 15q duplicons
are polymorphic and gene conversion moves them between copies, so a normal chromosome
can read as a "fusion" relative to CHM13. The HMM likewise switches copy number inside
the duplicons of every genome. Therefore:
- `breakpoints_vs_hificnv.tsv` takes the deletion as the CN=1 run overlapping the HiFiCNV
  interval (runs separated only by bins inside that interval are joined) and reports its
  two outer edges; it does not pick "the nearest 2->1 switch".
- `nahr_test.tsv` tests only those two edges, and counts junction-like reads of that pair
  whose crossover lies at the carrier's own edge, against the panel at the same position.
- `junction_enrichment.tsv` reports, per carrier and pair, junction-like reads observed
  vs expected from the panel rate; only an excess at a pair spanning the carrier's edges
  is evidence. `junction_crossovers.tsv` (raw counts) is not evidence by itself.
- `05b_breakpoint_summary.py` rebuilds these three tables from existing 05 outputs in
  seconds. Confirm candidate crossovers independently: split-read alignments (SA tag)
  when one side is unique sequence, and `assembly/<s>/<s>.junctions.tsv` (a contig of the
  deleted chromosome jumping across the deletion).

**Breakpoints (05) have two readouts.** *Junction reads* are single HiFi reads that
carry copy-A SUNKs followed by copy-B SUNKs of one direct SD pair — the NAHR fusion
itself. The crossover lies between the last A-SUNK and the first B-SUNK, so its
precision is the spacing of paralog-specific variants. Control, DiGeorge and
PWS-mUPD genomes should give none. *Copy number at SUNKs* scans the whole window
and finds CN changes in either copy, but read-depth noise limits it to several kb
and it can produce false segments at low depth. Report the junction-read crossover
as the breakpoint; use copy number for events without junction reads.

## Interpretation rules

- Repeat **sequence** has no parent-of-origin label; its **methylation** can.
  Differences between the two haplotypes of a biparental genome (06) are structural
  polymorphism between homologs.
- The retained chromosome in PWS-DEL/AS-DEL comes from the parent in whom the
  deletion did not arise: it samples the population, not the predisposing allele.
  Predisposition (e.g. a BP2–BP3 inversion) needs the transmitting parent's DNA.
- 07 calls a direct parent-of-origin element when the PWS-DEL and AS-DEL values do
  not overlap and differ by >= 0.20. With 5 + 3 carriers, non-overlap arises by chance
  for ~3.6 % of elements (2 of the 56 labellings), so 07 repeats the call for all 56
  labellings of the carriers: `validation.tsv` reports the calls expected under
  relabelling and a false-discovery estimate, `element_summary.tsv` the fraction of
  labellings calling each element, and `direct_delta_by_region.tsv` 250-kb medians of
  the difference against their permutation range (regional shifts).
- 07 follows `scripts/analysis/README.md`: parental direction only from PWS-DEL minus
  AS-DEL inside the common CN=1 core; Control and DiGeorge contribute unoriented
  |H1 − H2|. IC-anchored orientations (`--ic-anchored`, `--contig-anchored`) are
  sensitivity analyses written to `*_inferred.tsv.gz`, never reported as observed.
- Elements inside SDs or within 5 kb of the IC are reported but not tested.

## Checks before interpreting

- `WINDOW` (default chr15:17.5–33.0 Mb) must contain BP1–BP5 and every HiFiCNV
  deletion interval, including the 14.25-Mb event of 007P (17,592,000–31,842,000).
  `00_check_inputs.py` stops the run if one does not fit. After changing it, rerun
  01–03 with `--force` (the SUNK set and the extracted pieces depend on it).
- `CONTROL_REGION` (default chr15:60–64 Mb) must be CN 2 and SD-free in all 17
  genomes.
- In `scripts/analysis/01`, BP4 = 26.46 Mb and BP5 = 31.84 Mb fall outside the 15q13.3
  interval given for CHM13 by Höps et al. 2026 (27.5–30.8 Mb). Verify them before
  naming breakpoint classes.
- A BP block that reads ~2× the median depth when reads are mapped back to the
  assembly is collapsed; do not interpret its structure.

## Validation

- Unit tests (`tests/test_duplicons.py`): PAF projection on both strands; an NAHR
  deletion between 99.5%-identical direct duplicons recovered at homologous positions;
  a fusion read called with a crossover interval under 1 kb in both orientations,
  and a normal read not called; parent-of-origin repeat methylation recovered with
  the correct direction from randomly oriented biparental haplotypes; BAM and
  methylation-track lookup for workflow and SMRT Link names (ambiguous names and
  broken links refused, `bams.local.csv` override, several BAMs per participant).
- The eight-genome run below repeated with SMRT Link BAM names, 001P split into two
  BAMs and one BAM with a different chr15 length: identical SUNK counts and junction
  crossovers; the mismatched BAM was reported by 00 and skipped by 05.
- End-to-end run on a simulated eight-genome cohort (real minimap2, samtools, meryl;
  12-kb reads at 10× per haplotype; three NAHR deletions):
  - SD pair found exactly by 01;
  - junction reads: 7–10 per deletion carrier and 0 in the five diploid genomes;
    every crossover interval (41–350 bp) contained the true crossover;
  - copy number: all six switches were placed in the correct SD copy, but only three
    support intervals contained the true crossover (the other three missed it by
    0.1–9 kb), and one control showed a false CN 3 segment — hence the rule above;
  - 07: both testable parent-specific elements recovered, 0 false; mUPD |H1 − H2| 0.003.
- Not run on real data. RepeatMasker and pb-CpG-tools were replaced by stand-ins in
  the end-to-end test.

## Verify methylation score scales

The current `scripts/analysis/cis_analysis/methylation.py::read_track` divides
every pb-CpG-tools `mod_score` by 100. The duplicon analysis reader also selects
the score scale once per file. Other historical figure paths may still use a
per-value conversion: audit any figure that re-reads methylation BEDs directly
before using it for a manuscript. To estimate the number of low-score sites:

```bash
awk '!/^#/ && $4 > 0 && $4 <= 1' results/06_methylation/013A/013A.cpg.combined.bed | wc -l
```

If legacy outputs were made with per-value conversion, rerun the affected figure
analysis with the corrected reader before interpreting methylation differences.

## References

- Höps W, Porubsky D, … Eichler EE, Gilissen C. Distinct mechanisms of CNV formation
  at the human 15q13.3 locus. bioRxiv 2026. doi:10.64898/2026.03.03.709017
- Sudmant PH, et al. Diversity of human copy number variation and multicopy genes.
  Science 2010;330:641–646. doi:10.1126/science.1197005
- Gimelli G, et al. Genomic inversions of human chromosome 15q11–q13 in mothers of
  Angelman syndrome patients with class II (BP2/3) deletions. Hum Mol Genet
  2003;12:849–858. doi:10.1093/hmg/ddg101
