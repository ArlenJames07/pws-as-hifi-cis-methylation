#!/usr/bin/env python3
"""Wrapper for the publication-facing SNORD116 single-molecule renderer."""

from __future__ import annotations

import runpy
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[4]
SOURCE = PROJECT_ROOT / "scripts" / "paper_vf" / "update_figure4_manuscript_layout.py"


if __name__ == "__main__":
    if not SOURCE.exists():
        raise FileNotFoundError(f"Source script not found: {SOURCE}")
    runpy.run_path(str(SOURCE), run_name="__main__")
