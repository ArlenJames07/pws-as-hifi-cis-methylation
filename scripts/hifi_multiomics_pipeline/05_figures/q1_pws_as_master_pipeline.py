#!/usr/bin/env python3
"""Convenience wrapper for the Phase 1 PWS/AS Figure 1 generator."""

from __future__ import annotations

import runpy
from pathlib import Path


if __name__ == "__main__":
    script = Path(__file__).with_name("phase1_figure1_v2.py")
    runpy.run_path(str(script), run_name="__main__")
