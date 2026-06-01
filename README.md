# pws-as-hifi-cis-methylation

**PacBio HiFi workflows for haplotype- and molecule-resolved methylation analysis of the Prader–Willi/Angelman syndrome 15q11–q13 imprinted domain.**

---

## Repository information

| Field               | Description                                                         |
| ------------------- | ------------------------------------------------------------------- |
| **Pipeline author** | Arlen James Mosquera-Ruiz                                           |
| **Institution**     | Pontificia Universidad Javeriana Cali                               |
| **ORCID**           | [0009-0008-0796-9099](https://orcid.org/0009-0008-0796-9099)        |
| **Contact**         | [arlen22@javerianacali.edu.co](mailto:arlen22@javerianacali.edu.co) |

---

## Overview

This repository contains reproducible analysis scripts, configuration files, and non-identifiable processed summary tables supporting the manuscript:

**“Reciprocal imprinting disorder deletions reveal coordinated cis-methylation architecture of the human 15q11–q13 domain.”**

The workflow uses **PacBio HiFi long-read sequencing** aligned to **T2T-CHM13v2.0** to integrate:

* structural variant analysis
* copy-number profiling
* haplotype phasing
* native CpG methylation calling
* parent-of-origin methylation reconstruction
* SNHG14/ICR-proximal boundary mapping
* single-molecule methylation analysis
* breakpoint-flanking methylation assessment

The repository is designed to support transparent reproduction of the main analyses and manuscript figures for the PWS/AS 15q11–q13 imprinted domain.

---

## Scientific scope

This pipeline was developed to study how reciprocal deletion mechanisms in **Prader–Willi syndrome** and **Angelman syndrome** can be used as a natural hemizygous system to reconstruct parental methylation architecture across chromosome 15q11–q13.

The main analyses include:

1. validation of parent-of-origin methylation at the PWS/AS imprinting centre
2. reconstruction of maternal and paternal cis-methylation profiles
3. identification of the SNHG14/ICR-proximal methylation-transition boundary
4. molecule-level methylation analysis across the SNORD116 cluster
5. integration of methylation, copy number, structural variation and haplotype phase

---

## Repository contents

```text
pws-as-hifi-cis-methylation/
├── README.md
├── LICENSE
├── CITATION.cff
├── config/
│   ├── sample_metadata_template.tsv
│   ├── regions_15q11q13_chm13.bed
│   └── figure_config.yaml
├── scripts/
│   ├── 01_alignment_qc/
│   ├── 02_variant_cnv_phasing/
│   ├── 03_methylation_calling/
│   ├── 04_parental_methylation_architecture/
│   ├── 05_boundary_mapping/
│   ├── 06_single_molecule_analysis/
│   ├── 07_structural_context/
│   └── 08_figures/
├── data/
│   ├── example/
│   └── processed_summary_tables/
├── results/
│   ├── figures/
│   └── tables/
└── docs/
    ├── workflow_overview.md
    ├── input_requirements.md
    └── reproduce_figures.md
```

---

## Data availability and privacy

Raw human genomic data are **not included** in this repository because of participant privacy, consent restrictions, and the sensitive nature of long-read human genome sequencing data.

This repository may include:

* analysis scripts
* configuration files
* example input templates
* non-identifiable processed summary tables
* figure-generation code

This repository does **not** include:

* raw FASTQ files
* BAM/CRAM files
* full VCF files from individual genomes
* full per-sample methylation BED files
* identifiable clinical or genomic metadata

---

## Citation

If you use this pipeline, code, or processed summary outputs, please cite:

Mosquera-Ruiz A, Tobar-Tosse F, Londoño Velasco E, Lores J, Losada-Casallas KD, Ortega JG, Riccio-Rengifo C, Jaramillo-Botero A, Sharma A. (2026). *Reciprocal imprinting disorder deletions reveal coordinated cis-methylation architecture of the human 15q11–q13 domain*. Manuscript in preparation / under review.

---

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

---

## Contact

**Arlen James Mosquera-Ruiz**
Doctorate in Engineering and Applied Sciences
Pontificia Universidad Javeriana Cali
Cali, Colombia

Email: [arlen22@javerianacali.edu.co](mailto:arlen22@javerianacali.edu.co)
ORCID: [0009-0008-0796-9099](https://orcid.org/0009-0008-0796-9099)
