#!/usr/bin/env python3
"""Wrapper for the manuscript Figure 5 v7 structural-context renderer."""

from __future__ import annotations

import runpy
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[4]
SOURCE = PROJECT_ROOT / "scripts" / "paper_vf" / "make_figure5_structural_context_v2.py"


if __name__ == "__main__":
    if not SOURCE.exists():
        raise FileNotFoundError(f"Source script not found: {SOURCE}")
    runpy.run_path(str(SOURCE), run_name="__main__")
