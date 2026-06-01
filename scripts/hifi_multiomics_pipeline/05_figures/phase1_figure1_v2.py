#!/usr/bin/env python3
"""Wrapper for the manuscript Figure 1 generator."""

from __future__ import annotations

import runpy
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[4]
SOURCE = PROJECT_ROOT / "scripts" / "paper_vf" / "phase1_figure1_v2.py"


if __name__ == "__main__":
    if not SOURCE.exists():
        raise FileNotFoundError(f"Source script not found: {SOURCE}")
    runpy.run_path(str(SOURCE), run_name="__main__")
