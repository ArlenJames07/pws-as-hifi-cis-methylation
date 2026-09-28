#!/usr/bin/env python3
"""
06 -- Repeats predicted on the assembled 15q11-q13 haplotypes.

RepeatMasker hits on each haplotype's extracted pieces (03) are placed on
T2T-CHM13 through the piece alignments and compared with CHM13 over the SAME
aligned interval, so an assembly gap is never read as repeat loss. Every insertion
or deletion of a haplotype relative to CHM13 is annotated with its repeat content.

Reading the result:
- A repeat's SEQUENCE carries no parent-of-origin label; differences between the
  two haplotypes of a biparental genome are ordinary structural polymorphism
  (contig labels from 04, when present, are reported only for completeness).
- In PWS-DEL/AS-DEL the retained chromosome came from the parent in whom the
  deletion did NOT arise: it samples the population, not the predisposing allele.
- For deletion carriers, rows from the 'primary' assembly inside the CN=1 interval
  (column in_cn1_interval) describe that single retained chromosome; outside it,
  p_ctg mixes both parental chromosomes.

Outputs: results/08_duplicons/haplotype_repeats/
  per_sample/<sample>.{repeat_blocks,sv_repeats}.tsv
  cohort_repeat_blocks.tsv, cohort_sv_repeats.tsv
  repeat_blocks_by_mechanism.tsv, sv_repeats_by_class.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
from duplicon_analysis import cis_analysis, repeat_content  # noqa: E402
from duplicon_analysis.tools import log  # noqa: E402

MIN_SV_BP = 50
MIN_BLOCK_BP = 2_000
HAPLOTYPES = ("hap1", "hap2", "primary")


def gather(folder: Path, samples: list[str], suffix: str) -> pd.DataFrame:
    frames = []
    for s in samples:
        path = folder / f"{s}.{suffix}.tsv"
        if path.is_file() and path.stat().st_size > 1:
            frames.append(pd.read_csv(path, sep="\t"))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["sample"])


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args(argv)
    ref_out = C.REFERENCE_DIR / "window.fa.out"
    if not ref_out.is_file():
        raise SystemExit(f"{ref_out} missing: RepeatMasker must run in 01")
    out = C.HAPLOTYPE_REPEAT_DIR
    per_sample = out / "per_sample"
    per_sample.mkdir(parents=True, exist_ok=True)
    landmarks = C.BREAKPOINT_DIR / "bp_landmarks.bed"
    cohort = cis_analysis().load_cohort(C.METADATA_PATH)
    done = []
    for sample in cohort.samples:
        d = C.ASSEMBLY_DIR / sample
        args = []
        for hap in HAPLOTYPES:
            rm_out, paf = d / f"{sample}.{hap}.15q.fa.out", d / f"{sample}.{hap}.15q.vs_window.paf"
            if rm_out.is_file() and paf.is_file():
                args.append(f"--hap={hap}={rm_out},{paf}")
        if not args:
            log(f"{sample}: no assembled 15q pieces with RepeatMasker output; skipped")
            continue
        args += ["--sample", sample, "--ref-out", str(ref_out), "--sd-bed", str(C.REFERENCE_DIR / "window.sd.bed"),
                 "--min-sv", str(MIN_SV_BP), "--min-block", str(MIN_BLOCK_BP), "--outdir", str(per_sample)]
        if landmarks.is_file():
            args += ["--landmarks-bed", str(landmarks)]
        origin = d / f"{sample}.contig_origin.tsv"
        if origin.is_file():
            args += ["--origin", str(origin)]
        repeat_content.main(args)
        done.append(sample)
    if not done:
        raise SystemExit("no participant had assembled 15q pieces with RepeatMasker output")

    mechanism = {s: cohort.mechanism(s) for s in cohort.samples}
    blocks = gather(per_sample, done, "repeat_blocks")
    svs = gather(per_sample, done, "sv_repeats")
    for table in (blocks, svs):
        table.insert(1, "mechanism", table["sample"].map(mechanism))
    svs["in_cn1_interval"] = False
    if C.STRUCTURAL_EVIDENCE_PATH.is_file() and len(svs):
        ev = pd.read_csv(C.STRUCTURAL_EVIDENCE_PATH, sep="\t").set_index("sample_id")
        for i, row in svs.iterrows():
            if row["sample"] in ev.index and pd.notna(ev.loc[row["sample"], "cn_event_start"]):
                lo, hi = ev.loc[row["sample"], "cn_event_start"], ev.loc[row["sample"], "cn_event_end"]
                svs.at[i, "in_cn1_interval"] = bool(lo <= row["ref_start"] < hi)
    blocks.to_csv(out / "cohort_repeat_blocks.tsv", sep="\t", index=False)
    svs.to_csv(out / "cohort_sv_repeats.tsv", sep="\t", index=False)
    if len(blocks):
        (blocks.groupby(["mechanism", "feature", "class"])
         .agg(haplotypes=("sample", "size"), median_aligned_ref_bp=("aligned_ref_bp", "median"),
              median_delta_bp=("delta_bp", "median"), min_delta_bp=("delta_bp", "min"),
              max_delta_bp=("delta_bp", "max"))
         .reset_index().to_csv(out / "repeat_blocks_by_mechanism.tsv", sep="\t", index=False))
    if len(svs):
        svs["dominant_class"] = svs["repeat_classes"].fillna(".").str.split(":").str[0]
        svs["mostly_repeat"] = svs["repeat_fraction"] >= 0.8
        (svs.groupby(["type", "dominant_class", "mostly_repeat"])
         .agg(events=("sample", "size"), participants=("sample", "nunique"), median_length=("length", "median"))
         .reset_index().to_csv(out / "sv_repeats_by_class.tsv", sep="\t", index=False))
    log(f"{len(done)} participants; outputs in {out}")


if __name__ == "__main__":
    main()
