#!/usr/bin/env python3
"""
02 -- Read k-mer counts at the SUNKs, per participant.

For each participant: primary reads from its BAM(s) in results/01_alignment/
(<sample>.aligned.bam, SMRT Link names such as 08_1_A01_bc2043_001P.bam, or the
paths in bams.local.csv) are counted with meryl, intersected with the SUNK set
from 01, and printed as "KMER<TAB>COUNT". The reads are used whatever reference
they were aligned to. The large per-sample meryl database is deleted afterwards.

--window-reads counts only reads aligned to WINDOW (+/- 1 Mb) and CONTROL_REGION, plus
unmapped reads: minutes instead of hours per genome. Reads from the window that the
aligner placed on another chromosome are missed; use it for a supplementary window
(DUPLICON_OUT=...) and for the panel genomes of that run, not mixed with full counts.

Output: results/08_duplicons/sunk_counts/<sample>.sunk.k31.txt
Usage:  python3 scripts/duplicons/02_count_sunks.py [--samples 001P,002P] [--force] [--window-reads]
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
from duplicon_analysis import bams as bam_lookup  # noqa: E402
from duplicon_analysis import cis_analysis  # noqa: E402
from duplicon_analysis.assembly import parse_region  # noqa: E402
from duplicon_analysis.tools import log, parse_samples_arg, require, run, run_pipeline, up_to_date  # noqa: E402

SUNK_K = 31


def region_reads_command(bams: list[Path], threads: int) -> str:
    """Primary reads overlapping WINDOW +/- 1 Mb and CONTROL_REGION, plus unmapped reads."""
    chrom, w0, w1 = parse_region(C.WINDOW)
    intervals = sorted([(chrom, max(1, w0 - 1_000_000), w1 + 1_000_000), parse_region(C.CONTROL_REGION)])
    merged = [list(intervals[0])]
    for c, a, b in intervals[1:]:          # overlapping regions would output the same read twice
        if c == merged[-1][0] and a <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], b)
        else:
            merged.append([c, a, b])
    regions = " ".join(f"{c}:{a}-{b}" for c, a, b in merged)
    parts = []
    for b in bams:
        parts.append(f"{C.SAMTOOLS} view -u -F 0x900 {b} {regions} | {C.SAMTOOLS} fastq -@ {threads} -")
        parts.append(f"{C.SAMTOOLS} view -u -f 4 {b} '*' | {C.SAMTOOLS} fastq -@ {threads} -")
    return "{ " + "; ".join(parts) + "; }"


def count_sample(sample: str, bams: list[Path], threads: int, force: bool, window_reads: bool = False) -> None:
    out = C.SUNK_DIR / f"{sample}.sunk.k{SUNK_K}.txt"
    if up_to_date([out], force):
        log(f"{sample}: up to date")
        return
    if not bams:
        log(f"{sample}: no BAM in {C.ALIGNMENT_DIR} or {C.BAMS_PATH.name}; skipped")
        return
    log(f"{sample}: " + ", ".join(str(b) for b in bams))
    tmp = C.SUNK_DIR / f"tmp_{sample}"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    try:
        # samtools fastq drops secondary and supplementary records by default
        reads = (region_reads_command(bams, threads) if window_reads
                 else bam_lookup.fastq_command(C.SAMTOOLS, bams, threads))
        run_pipeline(f"{reads} | gzip -1 > {tmp / 'reads.fq.gz'}")
        run([C.MERYL, "count", f"k={SUNK_K}", f"threads={threads}", f"memory={C.MERYL_MEMORY_GB}",
             tmp / "reads.fq.gz", "output", tmp / "reads.meryl"])
        # intersect keeps the value of the FIRST database: the read count
        run([C.MERYL, "intersect", tmp / "reads.meryl", C.REFERENCE_DIR / f"sunks.k{SUNK_K}.meryl",
             "output", tmp / "sunks.meryl"])
        run([C.MERYL, "print", tmp / "sunks.meryl"], stdout=out)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    log(f"{sample}: {out}")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", help="comma-separated subset (default: every sample in metadata)")
    ap.add_argument("--threads", type=int, default=C.THREADS)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--window-reads", action="store_true",
                    help="count only reads aligned to WINDOW +/- 1 Mb, CONTROL_REGION and unmapped reads")
    a = ap.parse_args(argv)
    require(C.SAMTOOLS, C.MERYL)
    if not (C.REFERENCE_DIR / f"sunks.k{SUNK_K}.meryl").exists():
        raise SystemExit("SUNK set missing: run 01_prepare_reference.py first")
    C.SUNK_DIR.mkdir(parents=True, exist_ok=True)
    cohort = cis_analysis().load_cohort(C.METADATA_PATH)
    wanted = parse_samples_arg(a.samples)
    table = bam_lookup.read_bam_table(C.BAMS_PATH)
    problems = []
    for sample in cohort.samples:
        if wanted is not None and sample not in wanted:
            continue
        try:
            bams = bam_lookup.find_bams(sample, C.ALIGNMENT_DIR, table)
        except bam_lookup.BamLookupError as error:
            log(f"SKIPPED {error}")
            problems.append(str(error))
            continue
        count_sample(sample, bams, a.threads, a.force, a.window_reads)
    if problems:
        raise SystemExit("Fix these BAMs and rerun (finished samples are kept):\n  " + "\n  ".join(problems))


if __name__ == "__main__":
    main()
