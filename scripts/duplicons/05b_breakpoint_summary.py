#!/usr/bin/env python3
"""
05b -- Rebuild the per-carrier breakpoint tables from the existing outputs of 05, in
seconds (no reads are re-read):

  breakpoints_vs_hificnv.tsv, nahr_test.tsv, junction_enrichment.tsv

Step 05 calls the same code at its end; run this after updating the scripts, or after
changing BP_CLUSTERS in config.py, instead of rerunning 05.
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
from duplicon_analysis.tools import log  # noqa: E402


def main(argv=None) -> None:
    argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter).parse_args(argv)
    out = C.BREAKPOINT_DIR
    for name in ("sunk15q.bins.tsv", "sunk15q.transitions.tsv"):
        if not (out / name).is_file():
            raise SystemExit(f"{out / name} missing: run 05_duplicon_breakpoints.py first")
    importlib.import_module("05_duplicon_breakpoints").summarise(out)
    log(f"written: {out}/breakpoints_vs_hificnv.tsv, nahr_test.tsv, junction_enrichment.tsv")


if __name__ == "__main__":
    main()
