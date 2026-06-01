# HiFi Multiomics Pipeline

This folder groups the current long-read analysis scripts into a simple repo-style layout.

## Stages

1. `01_structural_variants/pbsv_discover.py`
   Structural variant discovery and calling with `pbmm2` and `pbsv`.
2. `02_phasing/hiphase.py`
   Phasing with `HiPhase` using aligned BAMs, SNVs, and SVs.
3. `03_hifi_cnvs/CNV.py`
   HiFi CNV calling from phased BAMs and phased variant calls.
4. `04_haplotype_methylation/pbcpgtools.py`
   Haplotype-aware CpG methylation calling from phased BAMs.

## Notes

- The files in this folder are wrappers that run the existing source scripts.
- Original script locations remain unchanged:
  - `/home/rare/arlen/scripts/SV_calling/pbsv_discover.py`
  - `/home/rare/arlen/scripts/Hiphase/hiphase.py`
  - `/home/rare/arlen/scripts/SV_calling/CNV.py`
  - `/home/rare/arlen/scripts/methylation_genomes/pbcpgtools.py`
