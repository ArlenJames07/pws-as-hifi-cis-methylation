#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ANALYSIS_SCRIPTS = (
    "00_build_prespecified_regions.py",
    "01_build_chr15_evidence_matrix.py",
    "02_depth_missingness_sensitivity.py",
    "03_reciprocal_cis_architecture.py",
)


def execute(script_name: str) -> None:
    command = [sys.executable, str(HERE / script_name)]
    print(" ".join(command), flush=True)
    subprocess.run(command, check=True)


def main() -> None:
    for script_name in ANALYSIS_SCRIPTS:
        execute(script_name)


if __name__ == "__main__":
    main()
