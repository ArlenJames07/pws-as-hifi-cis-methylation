# Figure Scripts

This folder now exposes one self-contained canonical script per manuscript
figure:

- `FIGURE_1.py`
- `FIGURE_2.py`
- `FIGURE_3.py`
- `FIGURE_4.py`
- `FIGURE_5.py`

Each `FIGURE_n.py` file now contains the figure-making code directly instead of
delegating to another script at runtime.

## Figure 1

- Canonical script: `FIGURE_1.py`
- Legacy aliases:
  - `phase1_figure1_v2.py`
  - `q1_pws_as_master_pipeline.py`
- Vendored source:
  - `/home/rare/arlen/scripts/paper_vf/phase1_figure1_v2.py`
- Primary outputs:
  - `Figure1_improved.png`
  - `Figure1_improved.pdf`
  - `Figure1_improved.svg`
  - `Figure1_ABC.png`
  - `Figure1_ABCD.png`

## Figure 2

- Canonical script: `FIGURE_2.py`
- Vendored sources:
  - `/home/rare/arlen/scripts/paper_vf/paper_vf_phase2_reciprocal_cis_architecture.py`
  - `/home/rare/arlen/scripts/paper_vf/create_figure2_reciprocal_cis_architecture_improved.py`
- Internal flow:
  1. run the reciprocal cis-architecture analysis
  2. render `Figure2_reciprocal_cis_architecture_improved`
- Primary outputs:
  - `Figure2_reciprocal_cis_architecture_improved.png`
  - `Figure2_reciprocal_cis_architecture_improved.pdf`
  - `Figure2_reciprocal_cis_architecture_improved.svg`

## Figure 3

- Canonical script: `FIGURE_3.py`
- Vendored sources:
  - `/home/rare/arlen/scripts/paper_vf/paper_vf_q1_pipeline.py`
  - `/home/rare/arlen/scripts/paper_vf/paper_vf_phase2_reciprocal_cis_architecture.py`
  - `/home/rare/arlen/scripts/paper_vf/create_figure2_reciprocal_cis_architecture_improved.py`
  - `/home/rare/arlen/scripts/paper_vf/phase3_boundary_mapping.py`
- Primary outputs:
  - `Figure3_boundary_mapping_improved.png`
  - `Figure3_boundary_mapping_improved.pdf`
  - `Figure3_boundary_mapping_improved.svg`

## Figure 4

- Canonical script: `FIGURE_4.py`
- Legacy alias:
  - `figure_snord116_single_molecule_architecture.py`
- Vendored source:
  - `/home/rare/arlen/scripts/paper_vf/update_figure4_manuscript_layout.py`
- Primary outputs:
  - `Figure_SNORD116_single_molecule_architecture.png`
  - `Figure_SNORD116_single_molecule_architecture.pdf`
  - `Figure_SNORD116_single_molecule_architecture.svg`
  - `Figure4_per_molecule_cis_architecture.png`
  - `Figure4_per_molecule_cis_architecture.pdf`

## Figure 5

- Canonical script: `FIGURE_5.py`
- Vendored source:
  - `/home/rare/arlen/scripts/paper_vf/make_figure5_structural_context_v2.py`
- Primary outputs:
  - `Figure5_v7.png`
  - `Figure5_v7.pdf`
