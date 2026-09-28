#!/usr/bin/env python3
"""
05 -- Deletion breakpoints inside the BP segmental duplications.

HiFiCNV loses read depth inside the BP1-BP5 duplicon blocks, so the CN=1
intervals from analysis 01 stop where unique sequence stops. Two independent
readouts locate the breakpoints inside the duplicons:

  copy number   paralog-specific copy number at the SUNKs finds the 2->1 and 1->2
                switches and tests whether they fall at homologous positions of one
                paralog pair (NAHR). Resolution is limited by read-depth noise
                (typically several kb at 20-30x). Control, DiGeorge and PWS-mUPD
                genomes (two copies of chr15) form the per-SUNK reference panel.
  junction reads single HiFi reads carrying copy-A SUNKs followed by copy-B SUNKs
                of one direct SD pair: the NAHR fusion read directly. The
                crossover lies between the last A-SUNK and the first B-SUNK, so
                the precision is the spacing of paralog-specific variants. Panel
                genomes should give none (negative control).

Inputs:  results/08_duplicons/reference/, results/08_duplicons/sunk_counts/,
         results/analysis/01_evidence_matrix/chr15_structural_evidence.tsv,
         the participants' BAMs (junction reads only; they must be aligned to the
         same reference as params.local.yml -- checked from the chr15 length)
Outputs: results/08_duplicons/breakpoints/
  sunk15q.transitions.tsv      switches with SUNK-level position and support interval
  sunk15q.nahr_pairs.tsv       homologous-position test per deletion participant
  sunk15q.bins.tsv             copy number per 5-kb bin, every participant
  sunk15q.overview.png         look at this first: panel genomes must read CN 2
  sunk15q.zoom_<sample>.png
  breakpoints_vs_hificnv.tsv   per carrier and edge: outer edge of the CN=1 run that carries
                               the deletion vs the HiFiCNV edge, SD block, BP cluster
  nahr_test.tsv                per carrier: are both edges at homologous positions of one
                               direct SD pair?
  junction_enrichment.tsv      junction-like reads in carriers vs the rate in genomes with
                               two copies of chr15 (they occur in normal chromosomes too)
  junction_reads.tsv           every read with SUNKs of both copies of a pair
  junction_crossovers.tsv      per participant and SD pair: junction reads and the
                               crossover interval in copy A and copy B coordinates
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
from duplicon_analysis import ANALYSIS_DIR, cis_analysis, junction_reads, sunk  # noqa: E402
from duplicon_analysis import bams as bam_lookup  # noqa: E402
from duplicon_analysis import breakpoint_summary  # noqa: E402
from duplicon_analysis.tools import log  # noqa: E402

SUNK_K = 31
BIN_BP = 5_000
MIN_SUNKS = 8
P_SWITCH = 1e-3


def breakpoint_landmarks() -> dict[str, int]:
    """BP1-BP5 positions as defined once in scripts/analysis/01 (not copied here)."""
    cis_analysis()
    path = ANALYSIS_DIR / "01_build_chr15_evidence_matrix.py"
    spec = importlib.util.spec_from_file_location("evidence_matrix_constants", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return dict(module.BREAKPOINT_LANDMARKS)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args(argv)
    ref = C.REFERENCE_DIR
    for name in ("window.fa", "control.fa", f"sunks.k{SUNK_K}.txt", "window.sd.bed", "window.sd_pairs.bedpe"):
        if not (ref / name).is_file():
            raise SystemExit(f"{ref / name} missing: run 01_prepare_reference.py")
    out = C.BREAKPOINT_DIR
    out.mkdir(parents=True, exist_ok=True)
    cohort = cis_analysis().load_cohort(C.METADATA_PATH)
    rows = []
    for sample in cohort.samples:
        counts = C.SUNK_DIR / f"{sample}.sunk.k{SUNK_K}.txt"
        if counts.is_file():
            rows.append({"sample": sample, "group": cohort.mechanism(sample), "counts": str(counts)})
        else:
            log(f"{sample}: no SUNK counts (run 02_count_sunks.py)")
    if not rows:
        raise SystemExit("no SUNK count files")
    sheet = out / "sunk_samples.tsv"
    pd.DataFrame(rows).to_csv(sheet, sep="\t", index=False)
    landmarks = out / "bp_landmarks.bed"
    landmarks.write_text("".join(f"{C.CHROM}\t{p}\t{p + 1}\t{n}\n" for n, p in breakpoint_landmarks().items()))
    sunk.main([
        "--k", str(SUNK_K), "--window-fasta", str(ref / "window.fa"), "--control-fasta", str(ref / "control.fa"),
        "--sunks", str(ref / f"sunks.k{SUNK_K}.txt"), "--samples", str(sheet),
        "--panel-groups", ",".join(C.MECHANISMS_DIPLOID_CHR15),
        "--sd-bed", str(ref / "window.sd.bed"), "--sd-pairs", str(ref / "window.sd_pairs.bedpe"),
        "--landmarks-bed", str(landmarks), "--bin-bp", str(BIN_BP), "--min-sunks", str(MIN_SUNKS),
        "--p-switch", str(P_SWITCH), "--outdir", str(out), "--prefix", "sunk15q",
    ])

    pairs, index = junction_reads.prepare(ref / "window.fa", ref / f"sunks.k{SUNK_K}.txt",
                                          ref / "window.sd_pairs.bedpe", SUNK_K)
    per_read, summary = [], []
    table = bam_lookup.read_bam_table(C.BAMS_PATH)
    for sample in cohort.samples:
        try:
            bams = bam_lookup.find_bams(sample, C.ALIGNMENT_DIR, table)
        except bam_lookup.BamLookupError as error:
            log(f"junction reads SKIPPED {error}")
            continue
        if not bams:
            log(f"{sample}: no BAM found; junction reads skipped")
            continue
        usable = []
        for bam in bams:
            ok, message = bam_lookup.reference_check(C.SAMTOOLS, bam, C.REFERENCE, C.CHROM)
            if not ok:
                log(f"{sample}: {bam.name}: {message}; junction reads skipped for this BAM")
            elif bam_lookup.bam_index(bam) is None:
                log(f"{sample}: {bam.name} has no index (samtools index); junction reads skipped for this BAM")
            else:
                usable.append(bam)
        if not usable:
            continue
        reads, crossovers = junction_reads.junction_reads(sample, usable, C.SAMTOOLS, pairs, index, out, SUNK_K)
        per_read.append(reads)
        if len(crossovers):
            crossovers.insert(1, "mechanism", cohort.mechanism(sample))
        summary.append(crossovers)
        n = int(crossovers["junction_reads"].sum()) if len(crossovers) else 0
        log(f"{sample} ({cohort.mechanism(sample)}): {n} junction read(s)")
    frames = [f for f in per_read if len(f)]
    (pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()).to_csv(
        out / "junction_reads.tsv", sep="\t", index=False)
    frames = [f for f in summary if len(f)]
    (pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()).to_csv(
        out / "junction_crossovers.tsv", sep="\t", index=False, na_rep="NA")
    summarise(out)
    log(f"done: {out}")


def summarise(out: Path) -> None:
    """breakpoints_vs_hificnv.tsv, nahr_test.tsv, junction_enrichment.tsv (also 05b)."""
    if not C.STRUCTURAL_EVIDENCE_PATH.is_file():
        log(f"{C.STRUCTURAL_EVIDENCE_PATH} not found; breakpoint summary skipped")
        return
    breakpoint_summary.summarise(out, C.STRUCTURAL_EVIDENCE_PATH, C.REFERENCE_DIR / "window.sd.bed",
                                 C.REFERENCE_DIR / "window.sd_pairs.bedpe", C.MECHANISMS_DELETION,
                                 C.MECHANISMS_DIPLOID_CHR15, C.BP_CLUSTERS, BIN_BP)


if __name__ == "__main__":
    main()
