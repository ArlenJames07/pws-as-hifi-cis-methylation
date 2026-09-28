#!/usr/bin/env python3
"""
01 -- Reference preparation for the 15q11-q13 duplicon analyses (run once).

  window.fa / control.fa     T2T-CHM13 BP1-BP5 window and a diploid control region
  window.self.paf            window aligned to itself (minimap2 -x asm20 -DP)
  window.sd_pairs.bedpe      paralog pairs (>= 10 kb, >= 90% identity), chr15 coordinates
  window.sd.bed              merged SD blocks
  sunks.k31.meryl / .txt     singly-unique k-mers: k-mers of window+control that occur
                             once in the whole of CHM13 (Sudmant et al. 2010)
  window.fa.out              RepeatMasker on the window (skipped if not installed)

Output: results/08_duplicons/reference/
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
from duplicon_analysis.assembly import window_offset  # noqa: E402
from duplicon_analysis.tools import log, require, run, up_to_date  # noqa: E402

SD_MIN_LENGTH = 10_000
SD_MIN_IDENTITY = 0.90
SUNK_K = 31


def sd_pairs_from_self_alignment(paf: Path, bedpe: Path, bed: Path) -> int:
    """Self-alignment PAF (window-relative) -> SD pairs and merged SD blocks (chr15)."""
    pairs, sides = [], []
    with open(paf) as handle:
        for n, line in enumerate(handle, 1):
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            length, matches = int(f[10]), int(f[9])
            if length < SD_MIN_LENGTH or matches / length < SD_MIN_IDENTITY:
                continue
            if f[2] == f[7] and f[3] == f[8]:          # the trivial diagonal
                continue
            chrom, off = window_offset(f[0])
            a0, a1 = int(f[2]) + off, int(f[3]) + off
            b0, b1 = int(f[7]) + off, int(f[8]) + off
            pairs.append((chrom, a0, a1, chrom, b0, b1, f"SDpair_{n}", f"{matches / length:.4f}", "+", f[4]))
            sides += [(chrom, a0, a1), (chrom, b0, b1)]
    with open(bedpe, "w") as out:
        out.writelines("\t".join(map(str, p)) + "\n" for p in pairs)
    merged = []
    for chrom, s, e in sorted(sides):
        if merged and merged[-1][0] == chrom and s <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], e)
        else:
            merged.append([chrom, s, e])
    with open(bed, "w") as out:
        out.writelines(f"{c}\t{s}\t{e}\tSD_{c}:{s}-{e}\n" for c, s, e in merged)
    return len(pairs)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--threads", type=int, default=C.THREADS)
    ap.add_argument("--force", action="store_true", help="recompute outputs that already exist")
    a = ap.parse_args(argv)
    require(C.SAMTOOLS, C.MINIMAP2, C.MERYL)
    if not C.REFERENCE.is_file():
        raise SystemExit(f"Reference not found: {C.REFERENCE} (set 'reference' in {C.PARAMS_FILE})")
    out = C.REFERENCE_DIR
    out.mkdir(parents=True, exist_ok=True)
    window, control = out / "window.fa", out / "control.fa"

    if not up_to_date([window, control], a.force):
        run([C.SAMTOOLS, "faidx", C.REFERENCE, C.WINDOW], stdout=window)
        run([C.SAMTOOLS, "faidx", C.REFERENCE, C.CONTROL_REGION], stdout=control)

    paf, bedpe, bed = out / "window.self.paf", out / "window.sd_pairs.bedpe", out / "window.sd.bed"
    if not up_to_date([paf], a.force):
        # -D ignores the diagonal, -P keeps every chain (parameters of Hoeps et al. 2026)
        run([C.MINIMAP2, "-t", a.threads, "-x", "asm20", "-c", "--eqx", "-DP", "-m200", window, window],
            stdout=paf)
    if not up_to_date([bedpe, bed], a.force):
        n = sd_pairs_from_self_alignment(paf, bedpe, bed)
        log(f"{n} SD pairs >= {SD_MIN_LENGTH} bp and >= {SD_MIN_IDENTITY:.0%} identity")

    sunk_db, sunk_txt = out / f"sunks.k{SUNK_K}.meryl", out / f"sunks.k{SUNK_K}.txt"
    if not up_to_date([sunk_db, sunk_txt], a.force):
        tmp = out / "meryl_tmp"
        tmp.mkdir(exist_ok=True)
        run([C.MERYL, "count", f"k={SUNK_K}", f"threads={a.threads}", f"memory={C.MERYL_MEMORY_GB}",
             C.REFERENCE, "output", tmp / "chm13.meryl"])
        run([C.MERYL, "equal-to", "1", tmp / "chm13.meryl", "output", tmp / "chm13.unique.meryl"])
        run([C.MERYL, "count", f"k={SUNK_K}", f"threads={a.threads}", window, control,
             "output", tmp / "targets.meryl"])
        if sunk_db.exists():
            shutil.rmtree(sunk_db)
        run([C.MERYL, "intersect", tmp / "targets.meryl", tmp / "chm13.unique.meryl", "output", sunk_db])
        run([C.MERYL, "print", sunk_db], stdout=sunk_txt)
        shutil.rmtree(tmp)
        log(f"{sum(1 for _ in open(sunk_txt)):,} SUNKs in window + control")

    rm_out = out / "window.fa.out"
    if shutil.which(C.REPEATMASKER) is None:
        log(f"{C.REPEATMASKER} not found: skipping RepeatMasker (needed by 06 and 07)")
    elif not up_to_date([rm_out], a.force):
        run([C.REPEATMASKER, "-species", "human", "-xsmall", "-pa", a.threads, "-dir", out, window])
    log(f"done: {out}")


if __name__ == "__main__":
    main()
