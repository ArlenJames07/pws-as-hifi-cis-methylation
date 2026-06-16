#!/usr/bin/env python3
"""Wrapper for the manuscript Figure 3 boundary-mapping generator."""

from __future__ import annotations

import runpy
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[4]
SOURCE_DIR = PROJECT_ROOT / "scripts" / "paper_vf"
SOURCE = SOURCE_DIR / "phase3_boundary_mapping.py"


if __name__ == "__main__":
    if not SOURCE.exists():
        raise FileNotFoundError(f"Source script not found: {SOURCE}")
    sys.path.insert(0, str(SOURCE_DIR))
    runpy.run_path(str(SOURCE), run_name="__main__")
