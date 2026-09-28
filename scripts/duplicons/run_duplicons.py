#!/usr/bin/env python3
"""
Run the duplicon, breakpoint and repeat analyses in order.

    python3 scripts/duplicons/run_duplicons.py                 # 00 (input check), 01, 02, 03, 05, 06, 07
    python3 scripts/duplicons/run_duplicons.py --from 05        # resume at a step (skips 00)
    python3 scripts/duplicons/run_duplicons.py --with-04        # include the optional
                                                                # contig-label step
Step 00 stops the run if an input problem would break or distort a later step.
Each step skips outputs that already exist; pass --force to recompute them.
Steps 00 and 02-04 accept --samples 001P,002P to run a subset.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEPS = (
    ("00", "00_check_inputs.py"),
    ("01", "01_prepare_reference.py"),
    ("02", "02_count_sunks.py"),
    ("03", "03_place_contigs.py"),
    ("04", "04_assembly_methylation.py"),
    ("05", "05_duplicon_breakpoints.py"),
    ("06", "06_haplotype_repeats.py"),
    ("07", "07_repeat_methylation.py"),
)
PER_SAMPLE = {"00", "02", "03", "04"}
CACHED = {"01", "02", "03", "04"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from", dest="start", default="00")
    ap.add_argument("--to", dest="stop", default="07")
    ap.add_argument("--with-04", action="store_true", help="run the optional contig-label step")
    ap.add_argument("--samples")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    for number, script in STEPS:
        if not (a.start <= number <= a.stop) or (number == "04" and not a.with_04):
            continue
        command = [sys.executable, str(HERE / script)]
        if a.samples and number in PER_SAMPLE:
            command += ["--samples", a.samples]
        if a.force and number in CACHED:
            command.append("--force")
        if number == "07" and a.with_04:
            command.append("--contig-anchored")
        print(" ".join(command), flush=True)
        subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
