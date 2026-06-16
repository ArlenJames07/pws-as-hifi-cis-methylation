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
5. `05_figures/phase1_figure1_v2.py`
   Figure generation for the manuscript, including `Figure1_improved.{png,pdf,svg}`, `Figure2_reciprocal_cis_architecture_improved.{png,pdf,svg}`, `Figure3_boundary_mapping_improved.{png,pdf,svg}`, `Figure_SNORD116_single_molecule_architecture.{png,pdf,svg}`, and `Figure5_v7.{png,pdf}`.

## Notes

- The files in this folder are wrappers that run the existing source scripts.
- Original script locations remain unchanged:
  - `/home/rare/arlen/scripts/SV_calling/pbsv_discover.py`
  - `/home/rare/arlen/scripts/Hiphase/hiphase.py`
  - `/home/rare/arlen/scripts/SV_calling/CNV.py`
  - `/home/rare/arlen/scripts/methylation_genomes/pbcpgtools.py`
  - `/home/rare/arlen/scripts/paper_vf/phase1_figure1_v2.py`
  - `/home/rare/arlen/scripts/paper_vf/q1_pws_as_master_pipeline.py`
  - `/home/rare/arlen/scripts/paper_vf/paper_vf_phase2_reciprocal_cis_architecture.py`
  - `/home/rare/arlen/scripts/paper_vf/create_figure2_reciprocal_cis_architecture_improved.py`
  - `/home/rare/arlen/scripts/paper_vf/phase3_boundary_mapping.py`
  - `/home/rare/arlen/scripts/paper_vf/phase4_per_molecule_cis_architecture.py`
  - `/home/rare/arlen/scripts/paper_vf/update_figure4_manuscript_layout.py`
  - `/home/rare/arlen/scripts/paper_vf/make_figure5_structural_context_v2.py`
