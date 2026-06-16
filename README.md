# pws-as-hifi-cis-methylation

**PacBio HiFi workflows for haplotype- and molecule-resolved methylation, structural-variant, copy-number and phasing analysis of the Prader–Willi/Angelman syndrome 15q11–q13 imprinted domain.**

---

## Repository information

| Field | Description |
|---|---|
| **Pipeline author** | Arlen James Mosquera-Ruiz |
| **Institution** | Pontificia Universidad Javeriana Cali |
| **ORCID** | [0009-0008-0796-9099](https://orcid.org/0009-0008-0796-9099) |
| **Contact** | [arlen22@javerianacali.edu.co](mailto:arlen22@javerianacali.edu.co) |

---

## Overview

This repository contains the analysis scripts supporting the manuscript:

**“Reciprocal imprinting disorder deletions reveal coordinated cis-methylation architecture of the human 15q11–q13 domain.”**

The workflow was developed for **PacBio HiFi long-read sequencing** data aligned to **T2T-CHM13v2.0**, integrating:

- structural variant discovery  
- haplotype phasing  
- copy-number analysis  
- native CpG methylation processing  
- haplotype-resolved methylation profiling  
- parent-of-origin methylation reconstruction  
- SNHG14/ICR-proximal boundary mapping  
- SNORD116 single-molecule methylation analysis  
- manuscript figure generation  

The repository is intended to support transparent reuse of the custom scripts used to analyse the PWS/AS 15q11–q13 imprinted domain. Raw human genomic data are not stored in this GitHub repository.

---

## Scientific scope

This pipeline uses reciprocal molecular classes of **Prader–Willi syndrome** and **Angelman syndrome** as a natural hemizygous system to reconstruct parental methylation architecture across chromosome 15q11–q13.

The main analyses include:

1. validation of parent-of-origin methylation at the PWS/AS imprinting centre  
2. reconstruction of maternal and paternal cis-methylation profiles  
3. identification of a focal SNHG14/ICR-proximal methylation-transition boundary  
4. molecule-level methylation analysis across the SNORD116 cluster  
5. integration of methylation, copy number, structural variation and haplotype phase  

---

## Repository structure

Current repository layout:

```text
pws-as-hifi-cis-methylation/
├── README.md
└── scripts/
    └── hifi_multiomics_pipeline/
        ├── README.md
        ├── 01_structural_variants/
        │   └── pbsv_discover.py
        ├── 02_phasing/
        │   └── hiphase.py
        ├── 03_hifi_cnvs/
        │   └── CNV.py
        ├── 04_haplotype_methylation/
        │   └── pbcpgtools.py
        └── 05_figures/
            ├── README.md
            ├── FIGURE_1.py
            ├── FIGURE_2.py
            ├── FIGURE_3.py
            ├── FIGURE_4.py
            ├── FIGURE_5.py
            ├── q1_pws_as_master_pipeline.py
            ├── phase1_figure1_v2.py
            ├── paper_vf_phase2_reciprocal_cis_architecture.py
            ├── create_figure2_reciprocal_cis_architecture_improved.py
            ├── phase3_boundary_mapping.py
            ├── phase4_per_molecule_cis_architecture.py
            ├── figure_snord116_single_molecule_architecture.py
            ├── update_figure4_manuscript_layout.py
            └── make_figure5_structural_context_v2.py
```

---

## Pipeline modules

| Module | Script | Purpose |
|---|---|---|
| `01_structural_variants` | `pbsv_discover.py` | Structural-variant discovery and processing using PacBio HiFi data |
| `02_phasing` | `hiphase.py` | Haplotype phasing of long-read variant calls |
| `03_hifi_cnvs` | `CNV.py` | Copy-number analysis across chromosome 15 and the PWS/AS critical region |
| `04_haplotype_methylation` | `pbcpgtools.py` | Processing of PacBio native CpG methylation calls and haplotype-resolved methylation information |
| `05_figures` | figure scripts | Generation of manuscript figures and downstream summary visualizations |

## Data availability and privacy

Raw human genomic data are **not included** in this GitHub repository because of participant privacy, consent restrictions and the sensitive nature of long-read human genome sequencing data.

This repository may include:

- analysis scripts  
- figure-generation code  
- configuration examples  
- non-identifiable processed summaries, when applicable  

This repository does **not** include:

- raw FASTQ files  
- BAM, BAI, CRAM or CRAI files  
- full individual-level VCF files  
- full individual-level methylation BED files  
- identifiable clinical metadata  
- participant-level private genomic information  

Raw PacBio HiFi sequencing data and aligned BAM files used in the study are being deposited under NCBI BioProject accession **PRJNA1469122**:

**“Implementation and evaluation of a predictive genomic association model for Rare Diseases based on DNA repeat configurations and structural variants.”**

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
