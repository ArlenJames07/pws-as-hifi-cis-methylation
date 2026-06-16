#!/usr/bin/env python3
"""Legacy alias for the canonical Figure 1 wrapper."""

from __future__ import annotations

import runpy
from pathlib import Path


if __name__ == "__main__":
    script = Path(__file__).with_name("FIGURE_1.py")
    runpy.run_path(str(script), run_name="__main__")
