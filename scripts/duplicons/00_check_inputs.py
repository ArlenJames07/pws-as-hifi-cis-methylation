#!/usr/bin/env python3
"""
00 -- Check every input before the long steps start. Reads only; writes one table.

Per participant in assets/metadata.csv:
  bam          which BAM(s) will be used (<sample>.aligned.bam, a SMRT Link name that
               contains the sample ID, or bams.local.csv), whether it is indexed,
               aligned to the reference in params.local.yml (chr15 length), and
               whether its reads carry MM/ML tags (needed by 04 only)
  assemblies   which of hap1 / hap2 / primary are listed in assemblies.local.csv and
               whether that matches the group (deletion carriers need primary)
  methylation  combined / hap1 / hap2 pb-CpG-tools tracks found (07)
  cnv          HiFiCNV output folder present (analysis 01)
  deletion     PWS-DEL / AS-DEL: status and CN=1 interval from analysis 01, and
               whether the interval lies inside WINDOW

Output: results/08_duplicons/input_check.tsv (and the same table on screen).
Exit status 1 if something would stop or silently distort a later step (a BAM that
cannot be chosen or is aligned to another reference, a missing tool or reference,
a deletion outside WINDOW); --warn-only reports without failing.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
from duplicon_analysis import bams as bam_lookup  # noqa: E402
from duplicon_analysis import cis_analysis  # noqa: E402
from duplicon_analysis.assembly import parse_region  # noqa: E402
from duplicon_analysis.methylation import find_track  # noqa: E402
from duplicon_analysis.tools import log, parse_samples_arg  # noqa: E402

CNV_DIR = C.RESULTS / "05_cnv"


def check_global() -> tuple[list[str], list[str]]:
    blocking, warnings = [], []
    if not C.PARAMS_FILE.is_file():
        blocking.append(f"params file not found: {C.PARAMS_FILE} (set PARAMS_FILE=...)")
    if not C.REFERENCE.is_file():
        blocking.append(f"reference not found: {C.REFERENCE} ('reference:' in {C.PARAMS_FILE.name})")
    elif not Path(f"{C.REFERENCE}.fai").is_file():
        warnings.append(f"{C.REFERENCE}.fai missing: run samtools faidx on the reference "
                        "(needed to compare BAM and reference)")
    for tool, label in ((C.SAMTOOLS, "samtools"), (C.MINIMAP2, "minimap2"), (C.MERYL, "meryl"),
                        (C.REPEATMASKER, "RepeatMasker")):
        if shutil.which(tool) is None and not Path(tool).is_file():
            (blocking if label != "RepeatMasker" else warnings).append(
                f"{label} not found ({tool})" + ("; 03 and 06 need it" if label == "RepeatMasker" else ""))
    if shutil.which(C.PBCPGTOOLS) is None and not Path(C.PBCPGTOOLS).is_file():
        warnings.append(f"pb-CpG-tools not found ({C.PBCPGTOOLS}); only the optional step 04 needs it")
    if not C.ASSEMBLIES_PATH.is_file():
        blocking.append(f"{C.ASSEMBLIES_PATH} not found (03)")
    if not C.STRUCTURAL_EVIDENCE_PATH.is_file():
        blocking.append(f"{C.STRUCTURAL_EVIDENCE_PATH} not found: run "
                        "scripts/analysis/01_build_chr15_evidence_matrix.py first")
    return blocking, warnings


def bam_columns(sample: str, table: dict, blocking: list[str]) -> dict:
    row = {"bam": "", "bam_index": "", "bam_reference": "", "mm_ml_tags": ""}
    try:
        bams = bam_lookup.find_bams(sample, C.ALIGNMENT_DIR, table)
    except bam_lookup.BamLookupError as error:
        row["bam"] = "ERROR"
        blocking.append(str(error))
        return row
    if not bams:
        row["bam"] = "none found"
        return row
    row["bam"] = ";".join(b.name for b in bams)
    row["bam_index"] = ";".join("yes" if bam_lookup.bam_index(b) else "NO" for b in bams)
    refs, tags = [], []
    for bam in bams:
        try:
            ok, message = bam_lookup.reference_check(C.SAMTOOLS, bam, C.REFERENCE, C.CHROM)
        except Exception as error:  # unreadable BAM
            ok, message = False, f"cannot read header ({error})"
        refs.append("ok" if ok else "DIFFERENT")
        if not ok:
            blocking.append(f"{sample}: {bam.name}: {message}")
        if ok and bam_lookup.bam_index(bam):
            present = bam_lookup.has_methylation_tags(C.SAMTOOLS, bam, f"{C.CHROM}:{C.DOMAIN_START}-{C.DOMAIN_END}")
            tags.append({True: "yes", False: "no", None: "?"}[present])
        else:
            tags.append("?")
    row["bam_reference"] = ";".join(refs)
    row["mm_ml_tags"] = ";".join(tags)
    if "NO" in row["bam_index"]:
        blocking.append(f"{sample}: BAM without index; run samtools index on it")
    return row


def assembly_columns(sample: str, mechanism: str, assemblies: dict) -> dict:
    present = sorted(assemblies.get(sample, {}), key=["hap1", "hap2", "primary"].index)
    need = ["primary"] if mechanism in C.MECHANISMS_DELETION else ["hap1", "hap2"]
    status = "ok" if set(need) <= set(present) else f"needs {'+'.join(need)}"
    return {"assemblies": ",".join(present) or "none", "assembly_status": status}


def methylation_columns(sample: str) -> dict:
    found = []
    for kind in ("combined", "hap1", "hap2"):
        try:
            found.append(kind if find_track(C.METHYLATION_DIR, sample, kind) else f"no {kind}")
        except ValueError:
            found.append(f"AMBIGUOUS {kind}")
    return {"methylation_tracks": ",".join(found)}


def deletion_columns(sample: str, mechanism: str, evidence: pd.DataFrame | None, window, blocking) -> dict:
    row = {"cnv_output": "yes" if (CNV_DIR / sample).is_dir() else "no", "deletion": ""}
    if mechanism not in C.MECHANISMS_DELETION or evidence is None:
        return row
    hit = evidence[evidence["sample_id"] == sample]
    if hit.empty:
        row["deletion"] = "not in evidence table"
        return row
    r = hit.iloc[0]
    status = r.get("ic_deletion_status", "")
    if pd.isna(r.get("cn_event_start")):
        row["deletion"] = f"{status}; no CN=1 interval"
        return row
    start, end = int(r["cn_event_start"]), int(r["cn_event_end"])
    inside = window[1] <= start and end <= window[2]
    row["deletion"] = f"{status}; {start:,}-{end:,} ({(end - start) / 1e6:.2f} Mb)" + ("" if inside else " OUTSIDE WINDOW")
    if not inside:
        blocking.append(f"{sample}: deletion {start:,}-{end:,} is not inside WINDOW {C.WINDOW}; "
                        "widen it with DUPLICON_WINDOW and rerun 01 with --force")
    return row


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples")
    ap.add_argument("--warn-only", action="store_true", help="never exit with status 1")
    a = ap.parse_args(argv)
    blocking, warnings = check_global()
    cohort = cis_analysis().load_cohort(C.METADATA_PATH)
    wanted = parse_samples_arg(a.samples)
    table = bam_lookup.read_bam_table(C.BAMS_PATH)
    try:
        import importlib
        assemblies = importlib.import_module("03_place_contigs").read_assemblies(C.ASSEMBLIES_PATH)
    except SystemExit as error:
        blocking.append(str(error))
        assemblies = {}
    evidence = pd.read_csv(C.STRUCTURAL_EVIDENCE_PATH, sep="\t") if C.STRUCTURAL_EVIDENCE_PATH.is_file() else None
    window = parse_region(C.WINDOW)

    rows = []
    for sample in cohort.samples:
        if wanted is not None and sample not in wanted:
            continue
        mechanism = cohort.mechanism(sample)
        row = {"sample": sample, "mechanism": mechanism}
        row.update(bam_columns(sample, table, blocking))
        row.update(assembly_columns(sample, mechanism, assemblies))
        row.update(methylation_columns(sample))
        row.update(deletion_columns(sample, mechanism, evidence, window, blocking))
        rows.append(row)
        if row["bam"] == "none found":
            warnings.append(f"{sample} ({mechanism}): no BAM; 02 and 05 will skip it")
        if row["assembly_status"] != "ok":
            warnings.append(f"{sample} ({mechanism}): assemblies.local.csv {row['assembly_status']}")

    report = pd.DataFrame(rows)
    C.OUT.mkdir(parents=True, exist_ok=True)
    report.to_csv(C.OUT / "input_check.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", None, "display.max_colwidth", 60):
        print(report.to_string(index=False), flush=True)
    print(f"\nBAM folder: {C.ALIGNMENT_DIR}   BAM table: "
          f"{C.BAMS_PATH if C.BAMS_PATH.is_file() else 'none (optional)'}   window: {C.WINDOW}")
    for message in warnings:
        log(f"WARNING {message}")
    for message in blocking:
        log(f"PROBLEM {message}")
    log(f"table: {C.OUT / 'input_check.tsv'}")
    if blocking and not a.warn_only:
        raise SystemExit(f"{len(blocking)} problem(s) above must be fixed before running 01-07 "
                         "(or rerun with --warn-only to continue anyway).")
    log("inputs OK" if not blocking else "problems reported (--warn-only)")


if __name__ == "__main__":
    main()
