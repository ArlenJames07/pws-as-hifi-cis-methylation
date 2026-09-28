#!/usr/bin/env python3
"""
Run the manuscript figures in order, after checking that their inputs exist.

  python3 scripts/figures/run_figures.py                 # Figures 2-7
  python3 scripts/figures/run_figures.py --check         # only list what is present / missing
  python3 scripts/figures/run_figures.py 5 7             # selected figures
  python3 scripts/figures/run_figures.py --render-only   # redraw 3, 4, 5, 7 from their tables
  python3 scripts/figures/run_figures.py --with-figure-1 # Figure 1 first (BAM + samtools; always
                                                         # reads results/ of this repository)

Order of the upstream steps (each writes under results/):
  scripts/analysis/run_analysis.py      -> results/analysis/        (00-03)
  scripts/duplicons/run_duplicons.py    -> results/08_duplicons/    (01-08)
  scripts/figures/FIGURE_1.py           -> results/07_figures/figure_1/ (parental assignment used by 3, 4)
Options --results, --params, --metadata and --gtf are passed to every figure.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from figlib import paths as paths_lib  # noqa: E402

RENDER_ONLY = {3, 4, 5, 7}


def requirements(p) -> dict[int, list[tuple[Path, bool, str]]]:
    """figure -> [(path, required, what it gives)]."""
    a, d = p.analysis, p.duplicons
    return {
        1: [(p.phasing, True, "HiPhase BAM/VCF"), (p.cnv, True, "HiFiCNV"), (p.methylation, True, "pb-CpG-tools")],
        2: [(a / "01_evidence_matrix" / "common_reciprocal_cn1_core.tsv", True, "analysis 01"),
            (a / "01_evidence_matrix" / "chr15_window_evidence_matrix.tsv.gz", True, "analysis 01"),
            (a / "03_cis_architecture" / "parent_associated_windows.tsv.gz", True, "analysis 03"),
            (a / "annotation_genes.tsv", True, "analysis 00"),
            (d / "repeat_methylation" / "direct_delta_by_region.tsv", False, "duplicons 07 (repeat bins, b-c)"),
            (d / "followup" / "methylation_domains.tsv", False, "duplicons 08 (domain, b-c)")],
        3: [(p.methylation, True, "pb-CpG-tools hap1/hap2/combined"),
            (p.figures / "figure_1" / "tables" / "Figure1_parental_like_assignment.tsv", False,
             "Figure 1 (control parent assignment; else IC methylation)"),
            (d / "repeat_methylation" / "element_summary.tsv", False, "duplicons 07 (repeats at the core, d)"),
            (d / "reference" / "window.sd.bed", False, "duplicons 01 (duplications)")],
        4: [(p.phasing, True, "HiPhase BAMs with MM/ML (pysam)"),
            (a / "prespecified_regions.tsv", True, "analysis 00"),
            (p.figures / "figure_1" / "tables" / "Figure1_parental_like_assignment.tsv", False,
             "Figure 1 (control parent assignment)")],
        5: [(a / "01_evidence_matrix" / "chr15_structural_evidence.tsv", True, "analysis 01"),
            (p.cnv, True, "HiFiCNV"), (p.sv_calls, False, "pbsv (d)"), (p.methylation, True, "pb-CpG-tools (e-f)"),
            (d / "breakpoints" / "sunk15q.bins.tsv", False, "duplicons 05 (SUNK copy number; else HiFiCNV)"),
            (d / "followup" / "breakpoint_status.tsv", False, "duplicons 08 (edges, class)")],
        6: [(d / "repeat_methylation" / "element_summary.tsv", True, "duplicons 07"),
            (d / "repeat_methylation" / "element_by_participant.tsv.gz", True, "duplicons 07"),
            (d / "repeat_methylation" / "validation.tsv", True, "duplicons 07"),
            (d / "followup" / "methylation_domains.tsv", False, "duplicons 08 (domain)"),
            (d / "haplotype_repeats" / "cohort_sv_repeats.tsv", False, "duplicons 06 (d)")],
        7: [(a / "01_evidence_matrix" / "chr15_structural_evidence.tsv", True, "analysis 01"),
            (d / "breakpoints" / "nahr_test.tsv", False, "duplicons 05 (a-c)"),
            (d / "breakpoints" / "sunk15q.bins.tsv", False, "duplicons 05 (b, e)"),
            (d / "followup" / "breakpoint_status.tsv", False, "duplicons 08 (c)"),
            (d / "followup" / "assembled_fusions.tsv", False, "duplicons 08 (d)"),
            (d / "assembly", False, "duplicons 03 (contig alignments, d-e)"),
            (d / "reference" / "window.sd_pairs.bedpe", False, "duplicons 01 (pairs, a-b)"),
            (p.phasing, False, "BAMs (split reads, e)")],
    }


def present(path: Path) -> bool:
    return path.is_dir() and any(path.iterdir()) if path.is_dir() else path.is_file() and path.stat().st_size > 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    paths_lib.add_common_arguments(ap)
    ap.add_argument("figures", nargs="*", type=int, help="figure numbers (default 2-7)")
    ap.add_argument("--with-figure-1", action="store_true", help="run FIGURE_1.py first")
    ap.add_argument("--render-only", action="store_true", help="Figures 3, 4, 5, 7: redraw from their tables")
    ap.add_argument("--check", action="store_true", help="list inputs and stop")
    ap.add_argument("--keep-going", action="store_true", help="run the remaining figures after a failure")
    a = ap.parse_args(argv)
    p = paths_lib.resolve(a)
    figs = a.figures or [2, 3, 4, 5, 6, 7]
    if a.with_figure_1 and 1 not in figs:
        figs = [1] + figs
    req = requirements(p)
    print(f"results: {p.results}\nduplicons: {p.duplicons}\n")
    blocked = set()
    for n in figs:
        for path, required, what in req.get(n, []):
            ok = present(path)
            if required and not ok:
                blocked.add(n)
            flag = "ok     " if ok else ("MISSING" if required else "absent ")
            print(f"Figure {n}  {flag}  {what:<55} {path}")
        print()
    if a.check:
        return 1 if blocked else 0
    passthrough = []
    for opt in ("results", "params", "metadata", "gtf"):
        v = getattr(a, opt)
        if v:
            passthrough += [f"--{opt}", v]
    failed = []
    for n in figs:
        if n in blocked:
            print(f"== Figure {n}: skipped (required input missing)")
            failed.append(n)
            continue
        # FIGURE_1.py takes no options: it always reads results/ of this repository
        cmd = [sys.executable, str(HERE / f"FIGURE_{n}.py"), *(passthrough if n != 1 else [])]
        if a.render_only and n in RENDER_ONLY:
            cmd.append("--render-only")
        print(f"== Figure {n}: {' '.join(cmd)}", flush=True)
        t0 = time.time()
        r = subprocess.run(cmd)
        print(f"   exit {r.returncode} after {time.time() - t0:.0f} s", flush=True)
        if r.returncode:
            failed.append(n)
            if not a.keep_going:
                break
    print("\nfailed or skipped: " + (", ".join(map(str, failed)) if failed else "none"))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
