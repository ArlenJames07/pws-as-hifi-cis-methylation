# pws-as-hifi-cis-methylation

Reproducible Nextflow DSL2 workflow for PacBio HiFi variant calling, phasing,
copy-number analysis, native CpG methylation, and manuscript figures in the
Prader–Willi/Angelman syndrome 15q11–q13 imprinted domain.

## What the workflow runs

```text
HiFi BAM
   │
   ├── pbmm2 alignment ── DeepVariant ── PASS SNVs/indels ─┐
   │                                                        ├── HiPhase
   └── pbsv discover ── pbsv call ── PASS SVs ─────────────┘     │
                                                                 ├── HiFiCNV
                                                                 ├── pb-CpG-tools
                                                                 └── optional Figures 1–3
```

Each sample is processed independently. Nextflow records the commands, resource
use, task hashes, execution report, timeline, trace, and DAG. Failed or
interrupted runs can be resumed without repeating completed tasks.

The older Python launchers under `scripts/01_variant_calling` through
`scripts/04_haplotype_methylation` are retained as analysis provenance. They
contain machine-specific paths and are not called by the Nextflow workflow.

## Repository structure

```text
pws-as-hifi-cis-methylation/
├── main.nf                    # workflow entry point and channel wiring
├── run_step.sh                # run one cumulative stage at a time
├── import_existing_results.sh # link completed legacy outputs; run no tools
├── nextflow.config            # defaults, resources, and execution profiles
├── nextflow_schema.json       # documented parameter interface
├── conf/
│   └── params.example.yml     # editable run configuration
├── assets/
│   ├── samplesheet.csv        # sample-sheet template (fake paths only)
│   └── metadata.csv           # non-identifiable figure metadata template
├── envs/                      # pinned Conda environments
├── modules/local/             # one Nextflow process per analysis stage
├── tests/data/                # tiny files for a safe stub test
└── scripts/
    ├── 01_variant_calling/    # legacy provenance scripts
    ├── 02_phasing/
    ├── 03_hifi_cnvs/
    ├── 04_haplotype_methylation/
    └── figures/               # manuscript figure programs
```

## Requirements

- Linux and Java 17 or newer
- Nextflow 24.10 or newer
- Mamba/Conda for `-profile conda`
- Docker for the optional `-profile docker` DeepVariant container
- `pb-CpG-tools` 2.3.2 installed locally, because it is distributed separately
  from the Conda environments in this repository

The Conda profile pins pbmm2, pbsv, bcftools, samtools, HiPhase, HiFiCNV, and
the Python figure libraries. DeepVariant is pinned to 1.8.0 by the Docker
profile. If all programs are already installed, use the `standard` profile.

## Configure a run

1. Copy the templates to files that will contain your local paths:

   ```bash
   cp assets/samplesheet.csv samplesheet.local.csv
   cp conf/params.example.yml params.local.yml
   ```

2. Edit `samplesheet.local.csv`. It must have this header and one existing HiFi
   BAM per row:

   ```csv
   sample,bam
   001P,/data/private/001P.hifi.bam
   002P,/data/private/002P.hifi.bam
   ```

   The tracked template lists the complete 17-sample cohort, including the six
   DiGeorge samples `008D`, `009D`, `010D`, `011D`, `012D`, and `015D`.
   Replace every placeholder BAM path before running. Every pre-figure stage
   reads this same sheet, so these samples are included automatically in
   alignment, small-variant calling, structural-variant calling, phasing, CNV,
   and methylation analysis.

3. Edit `params.local.yml` with absolute paths to the T2T-CHM13v2.0 reference
   resources, FASTA index, and the pb-CpG-tools model. Set
   `input: samplesheet.local.csv`.
   Local parameter files are ignored by Git when named `*.local.yml`.

Required parameters:

| Parameter | Purpose |
|---|---|
| `input` | CSV containing `sample,bam` |
| `reference` | T2T-CHM13v2.0 FASTA |
| `reference_fai` | FASTA index; defaults to `<reference>.fai` |
| `tandem_repeats` | tandem-repeat BED for `pbsv discover` |
| `cnv_exclude` | excluded-region BED for HiFiCNV |
| `cpg_model` | pb-CpG-tools `.tflite` pileup model |

If `aligned_bam_to_cpg_scores` is not on `PATH`, add its executable path to the
parameter file:

```yaml
pbcpgtools_bin: /opt/pb-CpG-tools/bin/aligned_bam_to_cpg_scores
```

The same pattern is supported for `pbmm2_bin`, `samtools_bin`,
`bcftools_bin`, `pbsv_bin`, `hiphase_bin`, and `hificnv_bin`. This allows a
local run to reproduce the exact executable locations recorded by the legacy
scripts without editing workflow modules.

## Run

### Reuse the outputs produced by the legacy scripts

If the programs in `scripts/01_variant_calling/` through
`scripts/04_haplotype_methylation/` have already been run, do not start the
Nextflow analysis again merely to recreate the numbered results layout. The
repository includes a non-computing importer that creates symbolic links to
the existing files:

```bash
# Preview every link without changing results/
./import_existing_results.sh --dry-run

# Create the links
./import_existing_results.sh
```

The importer uses the verified legacy locations and creates this layout:

```text
results/
├── 01_alignment/                         # /mnt/diskrare/.../aligned_reads/t2t
├── 02_small_variants/{raw,filtered}/
├── 03_structural_variants/{signatures,calls}/
├── 04_phasing/
├── 05_cnv/<sample>/
└── 06_methylation/<sample>/
```

It links all primary files found at the top level of the legacy folders,
including legacy reference or older samples, preserves the original
filenames, and never replaces an existing file or link. Symbolic links avoid
duplicating the large BAM, VCF, BigWig, and BED files. Keep the original legacy
directories available because the links depend on them. Existing downstream
PCA images inside the HiFiCNV directory are intentionally not imported because
this command stops before figure products; it imports the primary top-level
HiFiCNV outputs only.

This operation organizes existing outputs but does not add them to the
Nextflow task cache. Do not subsequently run `run_step.sh` for stages that you
intend to reuse only through these links; Nextflow cache reuse is based on
`work/` and `.nextflow/`, not on files present under `results/`.

### Run one stage at a time

The `--stage` option is cumulative. Each command stops after the requested
stage. Always keep the same `work_dir` and use `-resume`: previously completed
processes are then loaded from the Nextflow cache, and only the new stage runs.
The local `workstation` profile limits memory-heavy processes to one sample at
a time, which is appropriate for the configured 32-core host.

Run the following commands in order:

```bash
# Short form using the included launcher:
./run_step.sh alignment
./run_step.sh small_variants
./run_step.sh structural_variants
./run_step.sh phasing
./run_step.sh cnv
./run_step.sh methylation
```

The equivalent complete commands are:

```bash
# 1. Alignment -> results/01_alignment/
nextflow run main.nf -profile workstation,docker \
  -params-file params.local.yml --stage alignment -resume

# 2. DeepVariant and PASS filtering -> results/02_small_variants/
nextflow run main.nf -profile workstation,docker \
  -params-file params.local.yml --stage small_variants -resume

# 3. pbsv discovery/calling -> results/03_structural_variants/
nextflow run main.nf -profile workstation,docker \
  -params-file params.local.yml --stage structural_variants -resume

# 4. HiPhase -> results/04_phasing/
nextflow run main.nf -profile workstation,docker \
  -params-file params.local.yml --stage phasing -resume

# 5. HiFiCNV -> results/05_cnv/<sample>/
nextflow run main.nf -profile workstation,docker \
  -params-file params.local.yml --stage cnv -resume

# 6. pb-CpG-tools -> results/06_methylation/<sample>/
nextflow run main.nf -profile workstation,docker \
  -params-file params.local.yml --stage methylation -resume
```

For example, the structural-variant invocation reconstructs the dependency
graph through structural variants, but alignment and small-variant tasks are
reported as `Cached process`; they are not executed again. Do not remove
`work/` until all six stages are complete.

### Run all pre-figure stages together

The launcher accepts `all` and explicitly disables figure generation. It runs
alignment through methylation and stops after `results/06_methylation/`:

```bash
./run_step.sh all
```

To leave the complete pre-figure workflow running after closing the terminal:

```bash
mkdir -p logs
nohup ./run_step.sh all > logs/prefigure_pipeline.log 2>&1 &
echo $! > logs/prefigure_pipeline.pid
```

Monitor it with:

```bash
tail -f logs/prefigure_pipeline.log
```

The equivalent direct Nextflow command for the configured workstation uses
the local resource limits plus Docker for DeepVariant:

```bash
nextflow run main.nf \
  -profile workstation,docker \
  -params-file params.local.yml \
  --stage all \
  --run_figures false \
  -resume
```

Because `-resume` is enabled, successfully completed tasks are loaded from the
cache and only missing or interrupted work is submitted. Keep the same
`work_dir`, `.nextflow/` cache, input paths, reference, and program parameters.
The published files in `results/` alone are not sufficient for cache reuse.

On a machine where every dependency is already on `PATH`:

```bash
nextflow run main.nf -profile standard -params-file params.local.yml -resume
```

For SLURM compute nodes on which Docker is permitted, combine the scheduler
profile with the dependency profiles:

```bash
nextflow run main.nf -profile slurm,conda,docker -params-file params.local.yml -resume
```

View the command-line summary with `nextflow run main.nf --help`.


## Output structure

```text
results/
├── 01_alignment/
├── 02_small_variants/{raw,filtered}/
├── 03_structural_variants/{signatures,calls}/
├── 04_phasing/
├── 05_cnv/<sample>/
├── 06_methylation/<sample>/
├── 07_figures/
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.txt
    └── pipeline_dag.html
```

The `work/` directory is Nextflow's cache. Keep it while a run may need
`-resume`; remove it only after results and provenance reports have been safely
archived.

## Test the workflow structure

The stub test executes every core process without genomics software or real
human data:

```bash
nextflow run main.nf -profile test -stub-run
```

This checks the DSL2 graph, sample parsing, channel joins, declared outputs,
and publish paths. It does not validate scientific results.

## What belongs in Git

Commit workflow source, configuration templates, environment definitions,
documentation, and tiny test fixtures. Do not commit participant BAM/VCF/BED
files, reference genomes, Nextflow cache, or generated results. The included
`.gitignore` enforces those boundaries.

From the `main` branch, review and stage this reproducibility change with:

```bash
git status --short
git add \
  .gitignore README.md main.nf nextflow.config nextflow_schema.json \
  assets conf envs modules tests scripts/figures/README.md
git diff --cached --stat
git diff --cached
git commit -m "Add reproducible Nextflow HiFi analysis pipeline"
git push origin main
```

Always inspect `git diff --cached` before the commit. In particular, confirm
that no private genomic data or local absolute paths were staged.

## Data availability and privacy

Raw human genomic data are not included because of participant privacy,
consent restrictions, and the sensitive nature of long-read genome data. Do
not add FASTQ, BAM, CRAM, participant VCF, full methylation BED, or identifiable
clinical metadata to this repository.

The study data are being deposited under NCBI BioProject **PRJNA1469122**:
“Implementation and evaluation of a predictive genomic association model for
Rare Diseases based on DNA repeat configurations and structural variants.”

## Citation

Mosquera-Ruiz A, Tobar-Tosse F, Londoño Velasco E, Lores J, Losada-Casallas KD,
Ortega JG, Riccio-Rengifo C, Jaramillo-Botero A, Sharma A. (2026).
*Reciprocal imprinting disorder deletions reveal coordinated cis-methylation
architecture of the human 15q11–q13 domain*. Manuscript in preparation / under
review.

## Author and license

Arlen James Mosquera-Ruiz, Pontificia Universidad Javeriana Cali

[ORCID 0009-0008-0796-9099](https://orcid.org/0009-0008-0796-9099) ·
[arlen22@javerianacali.edu.co](mailto:arlen22@javerianacali.edu.co)

MIT License. See [LICENSE](LICENSE).
