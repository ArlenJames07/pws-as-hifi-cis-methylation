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
    └── 05_make_figures/       # manuscript figure programs
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

## Run

Recommended local execution uses Conda for the open Bioconda tools and Docker
for DeepVariant:

```bash
nextflow run main.nf \
  -profile conda,docker \
  -params-file params.local.yml \
  -resume
```

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

## Manuscript figures

Figures 1–3 can run after the core workflow by setting `run_figures: true` and
providing `gtf`, `metadata`, `icr_bed`, `segdup_bed`, `imprintome_bed`, and
`repeats_bed` in the parameter file.

Figures 4 and 5 are render-only in this repository. Their scripts require
precomputed `Figure4_*.tsv` or `Figure5*.tsv` tables that the current repository
does not generate. Supply those directories with `figure4_tables` and
`figure5_tables`; otherwise those two processes are skipped. This distinction
prevents the workflow from claiming that missing upstream analyses are
reproducible.

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
  assets conf envs modules tests scripts/05_make_figures/README.md
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
