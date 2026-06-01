# Figure Scripts

This stage exposes the manuscript Figure 1 generator inside the
`hifi_multiomics_pipeline` layout.

## Figure 1

- Script: `phase1_figure1_v2.py`
- Launcher: `q1_pws_as_master_pipeline.py`
- Primary outputs:
  - `Figure1_improved.png`
  - `Figure1_improved.pdf`
  - `Figure1_improved.svg`

These wrappers execute the current source scripts from:

- `/home/rare/arlen/scripts/paper_vf/phase1_figure1_v2.py`
- `/home/rare/arlen/scripts/paper_vf/q1_pws_as_master_pipeline.py`

## Figure 2

- Phase 2 generator: `paper_vf_phase2_reciprocal_cis_architecture.py`
- Improved renderer: `create_figure2_reciprocal_cis_architecture_improved.py`
- Primary outputs:
  - `Figure2_reciprocal_cis_architecture_improved.png`
  - `Figure2_reciprocal_cis_architecture_improved.pdf`
  - `Figure2_reciprocal_cis_architecture_improved.svg`

Run order:

1. `paper_vf_phase2_reciprocal_cis_architecture.py`
2. `create_figure2_reciprocal_cis_architecture_improved.py`

These wrappers execute the current source scripts from:

- `/home/rare/arlen/scripts/paper_vf/paper_vf_phase2_reciprocal_cis_architecture.py`
- `/home/rare/arlen/scripts/paper_vf/create_figure2_reciprocal_cis_architecture_improved.py`

## Figure 3

- Boundary mapper: `phase3_boundary_mapping.py`
- Primary outputs:
  - `Figure3_boundary_mapping_improved.png`
  - `Figure3_boundary_mapping_improved.pdf`
  - `Figure3_boundary_mapping_improved.svg`

This wrapper executes the current source script from:

- `/home/rare/arlen/scripts/paper_vf/phase3_boundary_mapping.py`

## Figure 4 / SNORD116 Single-Molecule Architecture

- Phase 4 generator: `phase4_per_molecule_cis_architecture.py`
- Publication renderer: `update_figure4_manuscript_layout.py`
- Primary outputs:
  - `Figure_SNORD116_single_molecule_architecture.png`
  - `Figure_SNORD116_single_molecule_architecture.pdf`
  - `Figure_SNORD116_single_molecule_architecture.svg`

Run order:

1. `phase4_per_molecule_cis_architecture.py`
2. `update_figure4_manuscript_layout.py`

These wrappers execute the current source scripts from:

- `/home/rare/arlen/scripts/paper_vf/phase4_per_molecule_cis_architecture.py`
- `/home/rare/arlen/scripts/paper_vf/update_figure4_manuscript_layout.py`

## Figure 5

- Structural-context renderer: `make_figure5_structural_context_v2.py`
- Primary outputs:
  - `Figure5_v7.png`
  - `Figure5_v7.pdf`

This wrapper executes the current source script from:

- `/home/rare/arlen/scripts/paper_vf/make_figure5_structural_context_v2.py`
