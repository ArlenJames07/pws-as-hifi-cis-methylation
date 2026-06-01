#!/usr/bin/env python3
"""Convenience wrapper for the SNORD116 single-molecule manuscript figure."""

from __future__ import annotations

import runpy
from pathlib import Path


if __name__ == "__main__":
    script = Path(__file__).with_name("update_figure4_manuscript_layout.py")
    runpy.run_path(str(script), run_name="__main__")
