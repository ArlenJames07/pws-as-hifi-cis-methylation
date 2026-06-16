#!/usr/bin/env python3
"""
Phase 2 reciprocal cis-architecture analysis for the PWS/AS manuscript.

This script implements the fixed-coordinate Phase 2 prompts:
  - 1 kb per-haplotype methylation profiles across chr15:22,500,000-28,500,000
  - control maternal/paternal reference architecture
  - PWS-DEL retained maternal and AS-DEL retained paternal natural dissection
  - PWS-mUPD maternal/maternal architecture
  - a Figure 2 composite with a gene/IC track and reciprocal-delta boundary panel

The input methylation BEDs in this project are produced by PacBio
aligned_bam_to_cpg_scores. For those files the methylation percentage is column
9 and per-site read coverage is column 6; these are the defaults here because
the Phase 2 prompt explicitly filters on haplotype read coverage.
"""

from __future__ import annotations

import argparse
import ast
import gzip
import math
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter

from paper_vf_q1_pipeline import (
    DEFAULT_CNV_DIR,
    DEFAULT_EXCLUDE,
    DEFAULT_GROUPS,
    DEFAULT_GTF,
    DEFAULT_METADATA_PATH,
    DEFAULT_METHYLATION_DIR,
    DEFAULT_OUTDIR,
    IC_PRIORITY_GENES,
    add_allele_labels,
    assign_allele_labels,
    detect_sample_code,
    discover_methylation_files,
    infer_ic_region_from_gtf,
    load_gtf_genes,
    load_metadata_table,
    mkdir,
    open_maybe_gzip,
    parse_attrs,
    read_methyl_region,
    sample_group,
    weighted_mean,
)


PHASE2_CHROM = "chr15"
PHASE2_START = 17_600_000
PHASE2_END = 28_000_000
PHASE2_WINDOW = 1_000
AS_IC_CENTER = 24_920_000
DEFAULT_PHASE2_SAMPLES = [
    "001P", "002P", "004P", "005P", "006P", "007P",
    "013A", "014A", "016A", "017C", "018C",
]

MATERNAL = "#c7254e"
PATERNAL = "#2266aa"
MUPD_1 = "#d65f2d"
MUPD_2 = "#9b287b"
GREY = "#626262"
LIGHT_GREY = "#d0d0d0"
DEFAULT_ICR_BED = Path("/home/rare/arlen/reference/ICR_t2t.bed")
DEFAULT_SEG_DUP = Path("/home/rare/arlen/reference/dupseg")
DEFAULT_BP_SOURCE_SCRIPT = Path("/home/rare/arlen/scripts/Daniela/pws_chr15_hifi_deletion_analysis.py")
FALLBACK_PWS_BREAKPOINTS = {"BP1": 20_940_000, "BP2": 21_070_000, "BP3": 26_050_000}


def load_pws_breakpoints_from_daniela(script_path: Path = DEFAULT_BP_SOURCE_SCRIPT) -> Dict[str, int]:
    """Read BP1/BP2/BP3 from Daniela's CHM13 deletion-analysis script."""
    try:
        tree = ast.parse(script_path.read_text())
        for node in tree.body:
            if not isinstance(node, ast.Assign):
                continue
            if not any(isinstance(target, ast.Name) and target.id == "PWS_BREAKPOINTS" for target in node.targets):
                continue
            raw = ast.literal_eval(node.value)
            breakpoints = {str(k): int(v) for k, v in raw.items() if str(k) in {"BP1", "BP2", "BP3"}}
            if set(breakpoints) == {"BP1", "BP2", "BP3"}:
                return breakpoints
    except Exception:
        pass
    return dict(FALLBACK_PWS_BREAKPOINTS)


PWS_BREAKPOINTS_T2T = load_pws_breakpoints_from_daniela()
BP_HOTSPOTS_T2T = [(name, PWS_BREAKPOINTS_T2T[name]) for name in ("BP1", "BP2", "BP3")]
BP_CLUSTER_INTERVALS_T2T = [
    {
        "name": "BP1",
        "start": 17_691_439,
        "end": 20_454_275,
        "approx_size": "2.763 Mb",
        "interpretation": "Proximal SD block / BP1 cluster",
    },
    {
        "name": "BP2",
        "start": 20_753_698,
        "end": 21_183_655,
        "approx_size": "430 kb",
        "interpretation": "BP2 SD block",
    },
    {
        "name": "BP3",
        "start": 25_875_912,
        "end": 26_632_507,
        "approx_size": "757 kb",
        "interpretation": "Distal PWS/AS BP3 SD block",
    },
]


def log(msg: str) -> None:
    print(f"[Phase2] {msg}", file=sys.stderr, flush=True)


@dataclass(frozen=True)
class WindowSpec:
    chrom: str
    start: int
    end: int
    size: int

    def frame(self) -> pd.DataFrame:
        starts = np.arange(self.start, self.end, self.size, dtype=int)
        out = pd.DataFrame({"window_start": starts})
        out["window_end"] = np.minimum(out["window_start"] + self.size, self.end).astype(int)
        out["window_mid"] = ((out["window_start"] + out["window_end"]) / 2).astype(int)
        return out


def phase_samples_from_metadata(file_table: pd.DataFrame, requested: Optional[List[str]] = None) -> List[str]:
    if requested:
        present = set(file_table["sample"].unique())
        return [s for s in requested if s in present]
    wanted_groups = {"CONTROL", "PWS_DEL", "PWS_mUPD", "AS_DEL"}
    samples = []
    for sample in sorted(file_table["sample"].unique()):
        if sample_group(sample, DEFAULT_GROUPS) in wanted_groups:
            samples.append(sample)
    return samples


def one_based_column(value, default: Optional[int] = None) -> Optional[int]:
    if value is None:
        return default
    if isinstance(value, str):
        v = value.strip().lower()
        if v in {"", "none", "null", "na"}:
            return None
        if v == "auto":
            return default
        return int(v)
    return int(value)


def read_bed_region_fast(
    path: Path,
    chrom: str,
    start: int,
    end: int,
    bed_meth_col,
    bed_cov_col,
    timeout: int = 240,
) -> pd.DataFrame:
    """Extract a region from plain-text BED methylation files without loading the genome."""
    if path is None or not path.exists() or str(path).endswith(".gz"):
        return pd.DataFrame(columns=["chrom", "start", "end", "meth", "coverage", "mid"])

    meth_col = one_based_column(bed_meth_col, default=9) or 9
    cov_col = one_based_column(bed_cov_col, default=6)
    meth_i = meth_col - 1
    cov_i = cov_col - 1 if cov_col is not None else None
    max_i = max(meth_i, cov_i if cov_i is not None else 0, 2)

    try:
        proc = subprocess.run(
            [
                "awk",
                "-v",
                f"chrom={chrom}",
                "-v",
                f"start={start}",
                "-v",
                f"end={end}",
                '$1==chrom && $2>=start && $2<end {print}',
                str(path),
            ],
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except Exception as exc:
        log(f"WARNING: fast BED extraction failed for {path.name}: {exc}")
        return pd.DataFrame(columns=["chrom", "start", "end", "meth", "coverage", "mid"])

    rows = []
    for line in proc.stdout.splitlines():
        f = line.split()
        if len(f) <= max_i:
            continue
        try:
            s = int(float(f[1]))
            e = int(float(f[2]))
            meth = float(f[meth_i])
            cov = float(f[cov_i]) if cov_i is not None else 1.0
        except Exception:
            continue
        if meth > 1.5:
            meth /= 100.0
        if not np.isfinite(meth):
            continue
        rows.append(
            {
                "chrom": f[0],
                "start": s,
                "end": e,
                "meth": min(1.0, max(0.0, meth)),
                "coverage": cov if np.isfinite(cov) and cov > 0 else 0.0,
                "mid": int((s + e) / 2),
            }
        )
    if not rows:
        return pd.DataFrame(columns=["chrom", "start", "end", "meth", "coverage", "mid"])
    return pd.DataFrame(rows)


def read_phase2_methyl_region(path: Path, chrom: str, start: int, end: int, args) -> pd.DataFrame:
    fmt = args.methyl_format
    if fmt == "auto":
        fmt = "dss" if "dss" in path.name.lower() else "bed"
    if fmt == "bed" and not str(path).endswith(".gz"):
        return read_bed_region_fast(path, chrom, start, end, args.bed_meth_col, args.bed_cov_col)
    return read_methyl_region(
        path,
        chrom,
        start,
        end,
        fmt=fmt,
        bed_meth_col=args.bed_meth_col,
        bed_cov_col=args.bed_cov_col,
        chunksize=args.chunksize,
    )


def read_layer_windows(path: Optional[Path], spec: WindowSpec, args) -> pd.DataFrame:
    base = spec.frame()
    if path is None:
        return base.assign(n_cpg=0, n_reads=0.0, mean_methylation=np.nan)

    m = read_phase2_methyl_region(path, spec.chrom, spec.start, spec.end, args)
    if m.empty:
        return base.assign(n_cpg=0, n_reads=0.0, mean_methylation=np.nan)

    m["window_start"] = ((m["mid"] - spec.start) // spec.size) * spec.size + spec.start
    m["window_end"] = m["window_start"] + spec.size
    rows = []
    for (ws, we), z in m.groupby(["window_start", "window_end"], sort=True):
        rows.append(
            {
                "window_start": int(ws),
                "window_end": int(min(we, spec.end)),
                "n_cpg": int(len(z)),
                "n_reads": float(z["coverage"].sum()),
                "mean_methylation": weighted_mean(z["meth"].values, z["coverage"].values),
            }
        )
    z = pd.DataFrame(rows)
    out = base.merge(z, on=["window_start", "window_end"], how="left")
    out["n_cpg"] = out["n_cpg"].fillna(0).astype(int)
    out["n_reads"] = out["n_reads"].fillna(0.0)
    return out


def sample_layer_paths(file_table: pd.DataFrame, sample: str) -> Dict[str, Optional[Path]]:
    z = file_table[file_table["sample"] == sample].copy()
    out: Dict[str, Optional[Path]] = {"hap1": None, "hap2": None, "combined": None}
    for layer in out:
        hit = z[z["layer"] == layer]
        if not hit.empty:
            out[layer] = Path(hit.iloc[0]["path"])
    return out


def summarize_sample_windows(
    sample: str,
    file_table: pd.DataFrame,
    spec: WindowSpec,
    labels: Dict[str, Dict[str, str]],
    args,
) -> pd.DataFrame:
    paths = sample_layer_paths(file_table, sample)
    h1 = read_layer_windows(paths["hap1"], spec, args).rename(
        columns={
            "n_cpg": "n_CpGs_haplotype1",
            "n_reads": "n_reads_haplotype1",
            "mean_methylation": "mean_meth_haplotype1_raw",
        }
    )
    h2 = read_layer_windows(paths["hap2"], spec, args).rename(
        columns={
            "n_cpg": "n_CpGs_haplotype2",
            "n_reads": "n_reads_haplotype2",
            "mean_methylation": "mean_meth_haplotype2_raw",
        }
    )
    combined = read_layer_windows(paths["combined"], spec, args).rename(
        columns={
            "n_cpg": "n_CpGs_combined",
            "n_reads": "n_reads_combined",
            "mean_methylation": "mean_meth_combined_raw",
        }
    )

    out = h1.merge(h2, on=["window_start", "window_end", "window_mid"], how="outer")
    out = out.merge(combined, on=["window_start", "window_end", "window_mid"], how="outer")
    for c in ["n_CpGs_haplotype1", "n_CpGs_haplotype2", "n_CpGs_combined"]:
        out[c] = out[c].fillna(0).astype(int)
    for c in ["n_reads_haplotype1", "n_reads_haplotype2", "n_reads_combined"]:
        out[c] = out[c].fillna(0.0)

    out["n_CpGs"] = out["n_CpGs_combined"]
    missing_combined = out["n_CpGs"] == 0
    out.loc[missing_combined, "n_CpGs"] = out.loc[
        missing_combined, ["n_CpGs_haplotype1", "n_CpGs_haplotype2"]
    ].max(axis=1)

    out["mean_meth_haplotype1"] = out["mean_meth_haplotype1_raw"]
    out["mean_meth_haplotype2"] = out["mean_meth_haplotype2_raw"]
    out.loc[
        (out["n_reads_haplotype1"] < args.min_reads) | (out["n_CpGs_haplotype1"] < args.min_cpgs),
        "mean_meth_haplotype1",
    ] = np.nan
    out.loc[
        (out["n_reads_haplotype2"] < args.min_reads) | (out["n_CpGs_haplotype2"] < args.min_cpgs),
        "mean_meth_haplotype2",
    ] = np.nan
    out["mean_meth_combined"] = out["mean_meth_combined_raw"]
    out.loc[
        (out["n_reads_combined"] < args.min_reads) | (out["n_CpGs_combined"] < args.min_cpgs),
        "mean_meth_combined",
    ] = np.nan

    lab = labels.get(sample, {})
    out["sample"] = sample
    out["group"] = sample_group(sample, DEFAULT_GROUPS)
    out["parental_label_haplotype1"] = lab.get("hap1", f"{out['group'].iloc[0]}_hap1_unresolved")
    out["parental_label_haplotype2"] = lab.get("hap2", f"{out['group'].iloc[0]}_hap2_unresolved")

    requested_cols = [
        "window_start",
        "window_end",
        "n_CpGs",
        "mean_meth_haplotype1",
        "mean_meth_haplotype2",
        "parental_label_haplotype1",
        "parental_label_haplotype2",
        "n_reads_haplotype1",
        "n_reads_haplotype2",
    ]
    aux_cols = [
        "sample",
        "group",
        "window_mid",
        "n_CpGs_haplotype1",
        "n_CpGs_haplotype2",
        "n_CpGs_combined",
        "mean_meth_combined",
        "n_reads_combined",
    ]
    return out[aux_cols + requested_cols]


def summarize_ic_for_labels(file_table: pd.DataFrame, chrom: str, ic_start: int, ic_end: int, args) -> pd.DataFrame:
    rows = []
    for _, r in file_table[file_table["layer"].isin(["hap1", "hap2"])].iterrows():
        path = Path(r["path"])
        m = read_phase2_methyl_region(path, chrom, ic_start, ic_end, args)
        rows.append(
            {
                "sample": r["sample"],
                "group": sample_group(r["sample"], DEFAULT_GROUPS),
                "layer": r["layer"],
                "n_cpg": int(len(m)),
                "total_coverage": float(m["coverage"].sum()) if not m.empty else 0.0,
                "mean_methylation": weighted_mean(m["meth"].values, m["coverage"].values) if not m.empty else np.nan,
            }
        )
    return pd.DataFrame(rows)


def labels_from_ic_and_group(ic_summary: pd.DataFrame, min_cpgs: int) -> Dict[str, Dict[str, str]]:
    allele_table = assign_allele_labels(ic_summary, min_cpgs=min_cpgs)
    labels: Dict[str, Dict[str, str]] = {}
    for _, r in allele_table.iterrows():
        labels.setdefault(r["sample"], {})[r["layer"]] = r["allele_label"]

    for sample in ic_summary["sample"].unique():
        group = sample_group(sample, DEFAULT_GROUPS)
        labels.setdefault(sample, {})
        if group == "PWS_mUPD":
            labels[sample]["hap1"] = "upd_maternal_like"
            labels[sample]["hap2"] = "upd_maternal_like"
        elif group == "PWS_DEL":
            labels[sample].setdefault("hap1", "pwsdel_hap1_unresolved")
            labels[sample].setdefault("hap2", "pwsdel_hap2_unresolved")
        elif group == "AS_DEL":
            labels[sample].setdefault("hap1", "asdel_hap1_unresolved")
            labels[sample].setdefault("hap2", "asdel_hap2_unresolved")
    return labels


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def nearby_plot_interval(start: int, end: int, spec: WindowSpec, rank: int, side: str, width: int = 8_000) -> Tuple[int, int]:
    if side == "left":
        x0 = spec.start + rank * 15_000
        return x0, min(x0 + width, spec.end)
    x1 = spec.end - rank * 15_000
    return max(spec.start, x1 - width), x1


def assign_plot_coordinates(df: pd.DataFrame, spec: WindowSpec, width: int = 8_000) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    out["within_phase2_window"] = (out["end"] > spec.start) & (out["start"] < spec.end)
    out["plot_start"] = out["start"].clip(lower=spec.start, upper=spec.end)
    out["plot_end"] = out["end"].clip(lower=spec.start, upper=spec.end)
    left_idx = out[out["end"] <= spec.start].sort_values("start").index
    for rank, idx in enumerate(left_idx):
        out.loc[idx, ["plot_start", "plot_end"]] = nearby_plot_interval(
            int(out.loc[idx, "start"]), int(out.loc[idx, "end"]), spec, rank, "left", width
        )
    right_idx = out[out["start"] >= spec.end].sort_values("start").index
    for rank, idx in enumerate(right_idx):
        out.loc[idx, ["plot_start", "plot_end"]] = nearby_plot_interval(
            int(out.loc[idx, "start"]), int(out.loc[idx, "end"]), spec, rank, "right", width
        )
    return out


def read_icr_annotations(path: Path, spec: WindowSpec, flank: int = 1_500_000) -> pd.DataFrame:
    rows = []
    if path is None or not Path(path).exists():
        return pd.DataFrame()
    with open_text(Path(path)) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] != spec.chrom:
                continue
            try:
                start = int(float(f[1]))
                end = int(float(f[2]))
            except Exception:
                continue
            if end <= spec.start - flank or start >= spec.end + flank:
                continue
            rows.append(
                {
                    "chrom": f[0],
                    "start": start,
                    "end": end,
                    "name": f[3] if len(f) > 3 else "ICR",
                    "parent": f[4] if len(f) > 4 else "",
                    "gene": f[5] if len(f) > 5 else "",
                    "source": str(path),
                }
            )
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows).sort_values(["start", "end"]).reset_index(drop=True)
    return assign_plot_coordinates(out, spec, width=7_000)


def read_segdup_annotations(path: Path, spec: WindowSpec) -> Tuple[pd.DataFrame, pd.DataFrame]:
    raw_rows = []
    if path is None or not Path(path).exists():
        return pd.DataFrame(), pd.DataFrame()
    with open_text(Path(path)) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] != spec.chrom:
                continue
            try:
                start = int(float(f[1]))
                end = int(float(f[2]))
            except Exception:
                continue
            if end <= spec.start or start >= spec.end:
                continue
            raw_rows.append(
                {
                    "chrom": f[0],
                    "start": max(start, spec.start),
                    "end": min(end, spec.end),
                    "name": f[3] if len(f) > 3 else "segdup",
                    "target": f[9] if len(f) > 9 else "",
                    "target_start": f[10] if len(f) > 10 else "",
                    "target_end": f[11] if len(f) > 11 else "",
                    "frac_match": f[23] if len(f) > 23 else "",
                    "aln_len": f[15] if len(f) > 15 else "",
                }
            )
    if not raw_rows:
        return pd.DataFrame(), pd.DataFrame()
    raw = pd.DataFrame(raw_rows).sort_values(["start", "end"]).reset_index(drop=True)

    merged = []
    for _, r in raw.sort_values(["start", "end"]).iterrows():
        s, e = int(r["start"]), int(r["end"])
        if not merged or s > merged[-1]["end"]:
            merged.append({"chrom": spec.chrom, "start": s, "end": e, "n_raw_segments": 1})
        else:
            merged[-1]["end"] = max(merged[-1]["end"], e)
            merged[-1]["n_raw_segments"] += 1
    blocks = pd.DataFrame(merged)
    blocks["size_bp"] = blocks["end"] - blocks["start"]
    blocks = blocks[blocks["size_bp"] >= 1_000].reset_index(drop=True)
    return raw, blocks


def nearest_interval(pos: int, intervals: pd.DataFrame, label_col: str) -> Tuple[str, int, str]:
    if intervals is None or intervals.empty:
        return "", np.nan, ""
    z = intervals.copy()
    starts = pd.to_numeric(z["start"], errors="coerce")
    ends = pd.to_numeric(z["end"], errors="coerce")
    dist = np.where((starts <= pos) & (ends >= pos), 0, np.minimum((starts - pos).abs(), (ends - pos).abs()))
    idx = int(np.nanargmin(dist))
    row = z.iloc[idx]
    label = str(row.get(label_col, ""))
    interval = f"{int(row['start'])}-{int(row['end'])}"
    return label, int(dist[idx]), interval


def annotate_deletion_breakpoints(
    deletion_intervals: pd.DataFrame,
    genes: pd.DataFrame,
    icrs: pd.DataFrame,
    segdup_blocks: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    if deletion_intervals.empty:
        return pd.DataFrame()
    gene_intervals = genes[["gene", "start", "end", "strand"]].copy() if not genes.empty else pd.DataFrame()
    icr_intervals = icrs.rename(columns={"name": "icr_name"}).copy() if icrs is not None and not icrs.empty else pd.DataFrame()
    segdup_intervals = segdup_blocks.copy() if segdup_blocks is not None and not segdup_blocks.empty else pd.DataFrame()
    if not segdup_intervals.empty:
        segdup_intervals["segdup_block"] = segdup_intervals["start"].astype(str) + "-" + segdup_intervals["end"].astype(str)

    for _, r in deletion_intervals.iterrows():
        for kind, pos in [("left_or_start", int(r["start"])), ("right_or_end", int(r["end"]))]:
            gene, gene_dist, gene_interval = nearest_interval(pos, gene_intervals, "gene")
            icr, icr_dist, icr_interval = nearest_interval(pos, icr_intervals, "icr_name")
            segdup, segdup_dist, segdup_interval = nearest_interval(pos, segdup_intervals, "segdup_block")
            rows.append(
                {
                    "sample": r["sample"],
                    "group": r["group"],
                    "chrom": r["chrom"],
                    "breakpoint_position": pos,
                    "breakpoint_side": kind,
                    "copy_number": r.get("copy_number", np.nan),
                    "nearest_gene": gene,
                    "distance_to_nearest_gene_bp": gene_dist,
                    "nearest_gene_interval": gene_interval,
                    "nearest_ICR": icr,
                    "distance_to_nearest_ICR_bp": icr_dist,
                    "nearest_ICR_interval": icr_interval,
                    "nearest_or_overlapping_segdup_block": segdup,
                    "distance_to_segdup_block_bp": segdup_dist,
                    "segdup_block_interval": segdup_interval,
                    "source": r.get("source", ""),
                }
            )
    return pd.DataFrame(rows).sort_values(["group", "sample", "breakpoint_position", "breakpoint_side"]).reset_index(drop=True)


def discover_cnv_files(cnv_dir: Path, samples: Iterable[str]) -> pd.DataFrame:
    suffixes = (".bedgraph", ".bedGraph", ".bed", ".bed.gz", ".tsv", ".txt")
    wanted = set(samples)
    rows = []
    if cnv_dir is None or not cnv_dir.exists():
        return pd.DataFrame(columns=["sample", "path"])
    for p in cnv_dir.rglob("*"):
        if not p.is_file() or not any(str(p).endswith(s) for s in suffixes):
            continue
        sample = detect_sample_code(p)
        if sample in wanted:
            rows.append({"sample": sample, "path": str(p)})
    return pd.DataFrame(rows).drop_duplicates().sort_values(["sample", "path"]).reset_index(drop=True)


def read_cnv_deletion_intervals(
    cnv_table: pd.DataFrame,
    spec: WindowSpec,
    copy_threshold: float,
) -> pd.DataFrame:
    rows = []
    for _, r in cnv_table.iterrows():
        path = Path(r["path"])
        try:
            with open_text(path) as fh:
                for line in fh:
                    if not line.strip() or line.startswith("#"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 4 or f[0] != spec.chrom:
                        continue
                    try:
                        s = int(float(f[1]))
                        e = int(float(f[2]))
                        value = float(f[3])
                    except Exception:
                        continue
                    if e <= spec.start or s >= spec.end or value > copy_threshold:
                        continue
                    rows.append(
                        {
                            "sample": r["sample"],
                            "group": sample_group(r["sample"], DEFAULT_GROUPS),
                            "chrom": spec.chrom,
                            "start": max(s, spec.start),
                            "end": min(e, spec.end),
                            "copy_number": value,
                            "source": str(path),
                        }
                    )
        except Exception as exc:
            log(f"WARNING: failed to read CNV file {path}: {exc}")
    if not rows:
        return pd.DataFrame(columns=["sample", "group", "chrom", "start", "end", "copy_number", "source"])
    out = pd.DataFrame(rows)
    out = out[out["end"] > out["start"]].copy()
    return out.sort_values(["sample", "start", "end"]).reset_index(drop=True)


def windows_in_intervals(windows: pd.DataFrame, intervals: pd.DataFrame, sample: Optional[str] = None) -> pd.Series:
    mask = pd.Series(False, index=windows.index)
    if intervals.empty:
        return mask
    z = intervals
    if sample is not None and "sample" in z.columns:
        z = z[z["sample"] == sample]
    for _, r in z.iterrows():
        mask |= (windows["window_end"] > int(r["start"])) & (windows["window_start"] < int(r["end"]))
    return mask


def infer_deletion_intervals_from_windows(all_windows: pd.DataFrame, spec: WindowSpec, min_reads: int) -> pd.DataFrame:
    rows = []
    for sample, z in all_windows.groupby("sample"):
        group = sample_group(sample, DEFAULT_GROUPS)
        if group not in {"PWS_DEL", "AS_DEL"}:
            continue
        zz = z.sort_values("window_start").copy()
        candidate = (
            (zz["n_reads_combined"] >= min_reads)
            & ((zz["n_reads_haplotype1"] < min_reads) | (zz["n_reads_haplotype2"] < min_reads))
        )
        for start_i, end_i in boolean_runs(candidate.to_numpy()):
            start = int(zz.iloc[start_i]["window_start"])
            end = int(zz.iloc[end_i]["window_end"])
            if end - start >= 50_000:
                rows.append(
                    {
                        "sample": sample,
                        "group": group,
                        "chrom": spec.chrom,
                        "start": start,
                        "end": end,
                        "copy_number": np.nan,
                        "source": "coverage_fallback",
                    }
                )
    return pd.DataFrame(rows)


def boolean_runs(mask: np.ndarray) -> List[Tuple[int, int]]:
    runs: List[Tuple[int, int]] = []
    start: Optional[int] = None
    for i, val in enumerate(mask):
        if bool(val) and start is None:
            start = i
        if start is not None and (not bool(val) or i == len(mask) - 1):
            end = i if bool(val) and i == len(mask) - 1 else i - 1
            runs.append((start, end))
            start = None
    return runs


def update_del_labels_by_coverage(
    labels: Dict[str, Dict[str, str]],
    all_windows: pd.DataFrame,
    deletion_intervals: pd.DataFrame,
    fold: float = 1.5,
) -> Dict[str, Dict[str, str]]:
    labels = {s: dict(v) for s, v in labels.items()}
    for sample in sorted(all_windows["sample"].unique()):
        group = sample_group(sample, DEFAULT_GROUPS)
        if group not in {"PWS_DEL", "AS_DEL"}:
            continue
        z = all_windows[all_windows["sample"] == sample].copy()
        mask = windows_in_intervals(z, deletion_intervals, sample=sample)
        if not mask.any():
            continue
        h1 = float(z.loc[mask, "n_reads_haplotype1"].sum())
        h2 = float(z.loc[mask, "n_reads_haplotype2"].sum())
        labels.setdefault(sample, {})
        if max(h1, h2) < 50:
            continue
        if h1 >= h2 * fold:
            retained = "hap1"
            absent = "hap2"
        elif h2 >= h1 * fold:
            retained = "hap2"
            absent = "hap1"
        else:
            continue
        if group == "PWS_DEL":
            labels[sample][retained] = "pwsdel_retained_maternal_like"
            labels[sample][absent] = "pwsdel_deleted_paternal_like"
        else:
            labels[sample][retained] = "asdel_retained_paternal_like"
            labels[sample][absent] = "asdel_deleted_maternal_like"
    return labels


def apply_labels_to_windows(all_windows: pd.DataFrame, labels: Dict[str, Dict[str, str]]) -> pd.DataFrame:
    out = all_windows.copy()
    out["parental_label_haplotype1"] = out["sample"].map(lambda s: labels.get(s, {}).get("hap1", "hap1_unresolved"))
    out["parental_label_haplotype2"] = out["sample"].map(lambda s: labels.get(s, {}).get("hap2", "hap2_unresolved"))
    return out


def write_per_sample_tables(all_windows: pd.DataFrame, outdir: Path) -> None:
    per_dir = mkdir(outdir / "tables" / "phase2_per_sample_windows")
    requested_cols = [
        "window_start",
        "window_end",
        "n_CpGs",
        "mean_meth_haplotype1",
        "mean_meth_haplotype2",
        "parental_label_haplotype1",
        "parental_label_haplotype2",
        "n_reads_haplotype1",
        "n_reads_haplotype2",
    ]
    for sample, z in all_windows.groupby("sample"):
        z = z.sort_values("window_start")
        z[requested_cols].to_csv(
            per_dir / f"Phase2_{sample}_haplotype_1kb_methylation.tsv",
            sep="\t",
            index=False,
            na_rep="NaN",
        )


def value_for_label(z: pd.DataFrame, label_contains: str) -> pd.Series:
    h1 = z["parental_label_haplotype1"].astype(str).str.contains(label_contains, regex=False)
    h2 = z["parental_label_haplotype2"].astype(str).str.contains(label_contains, regex=False)
    out = pd.Series(np.nan, index=z.index, dtype=float)
    out.loc[h1] = z.loc[h1, "mean_meth_haplotype1"]
    out.loc[h2] = z.loc[h2, "mean_meth_haplotype2"]
    return out


def control_reference(all_windows: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    controls = all_windows[all_windows["group"] == "CONTROL"].copy()
    rows = []
    sample_rows = []
    for sample, z in controls.groupby("sample"):
        zz = z.sort_values("window_start")
        maternal = value_for_label(zz, "maternal")
        paternal = value_for_label(zz, "paternal")
        tmp = zz[["sample", "window_start", "window_end", "window_mid"]].copy()
        tmp["control_maternal"] = maternal
        tmp["control_paternal"] = paternal
        sample_rows.append(tmp)
    sample_ref = pd.concat(sample_rows, ignore_index=True) if sample_rows else pd.DataFrame()
    if sample_ref.empty:
        return pd.DataFrame(), sample_ref
    for (ws, we, wm), z in sample_ref.groupby(["window_start", "window_end", "window_mid"]):
        rows.append(
            {
                "window_start": ws,
                "window_end": we,
                "window_mid": wm,
                "control_maternal_mean": float(np.nanmean(z["control_maternal"])),
                "control_paternal_mean": float(np.nanmean(z["control_paternal"])),
                "n_controls_maternal": int(z["control_maternal"].notna().sum()),
                "n_controls_paternal": int(z["control_paternal"].notna().sum()),
            }
        )
    ref = pd.DataFrame(rows).sort_values("window_start").reset_index(drop=True)
    ref["maternal_paternal_delta"] = ref["control_maternal_mean"] - ref["control_paternal_mean"]
    ref["abs_delta"] = ref["maternal_paternal_delta"].abs()
    return ref, sample_ref


def retained_profiles(all_windows: pd.DataFrame, deletion_intervals: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sample, z in all_windows.groupby("sample"):
        group = sample_group(sample, DEFAULT_GROUPS)
        if group not in {"PWS_DEL", "AS_DEL"}:
            continue
        zz = z.sort_values("window_start").copy()
        in_del = windows_in_intervals(zz, deletion_intervals, sample=sample)
        retained_kind = "maternal" if group == "PWS_DEL" else "paternal"
        outside = value_for_label(zz, retained_kind)
        retained = outside.copy()
        retained.loc[in_del] = zz.loc[in_del, "mean_meth_combined"]
        source = pd.Series("assigned_haplotype", index=zz.index, dtype=object)
        source.loc[in_del] = "combined_hemizygous_deletion"
        absent_cov = np.minimum(zz["n_reads_haplotype1"], zz["n_reads_haplotype2"])
        if group == "PWS_DEL":
            hap2_cov = zz["n_reads_haplotype2"]
        else:
            hap2_cov = absent_cov
        tmp = zz[
            [
                "sample",
                "group",
                "window_start",
                "window_end",
                "window_mid",
                "n_reads_haplotype1",
                "n_reads_haplotype2",
                "n_reads_combined",
            ]
        ].copy()
        tmp["retained_parent"] = retained_kind
        tmp["retained_methylation"] = retained
        tmp["retained_source"] = source
        tmp["within_deletion_footprint"] = in_del.to_numpy()
        tmp["absent_haplotype_coverage"] = absent_cov
        tmp["prompt_haplotype2_coverage"] = hap2_cov
        rows.append(tmp)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def pearson_summary(a: pd.Series, b: pd.Series) -> Dict[str, float]:
    m = pd.DataFrame({"a": a, "b": b}).dropna()
    if len(m) < 3 or m["a"].std() == 0 or m["b"].std() == 0:
        r = np.nan
    else:
        r = float(m["a"].corr(m["b"]))
    return {
        "n_windows": int(len(m)),
        "pearson_r": r,
        "mean_abs_delta": float((m["a"] - m["b"]).abs().mean()) if len(m) else np.nan,
        "mean_delta": float((m["a"] - m["b"]).mean()) if len(m) else np.nan,
    }


def compute_correlations(
    retained: pd.DataFrame,
    all_windows: pd.DataFrame,
    control_ref: pd.DataFrame,
    deletion_intervals: pd.DataFrame,
    deviation_threshold: float,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    corr_rows = []
    for sample, z in retained.groupby("sample"):
        group = sample_group(sample, DEFAULT_GROUPS)
        if group not in {"PWS_DEL", "AS_DEL"}:
            continue
        ref_col = "control_maternal_mean" if group == "PWS_DEL" else "control_paternal_mean"
        zz = z.merge(control_ref[["window_start", ref_col]], on="window_start", how="left")
        zz = zz[zz["within_deletion_footprint"]]
        stats = pearson_summary(zz["retained_methylation"], zz[ref_col])
        intervals = deletion_intervals[deletion_intervals["sample"] == sample]
        corr_rows.append(
            {
                "sample": sample,
                "group": group,
                "comparison": f"{sample}_retained_vs_{ref_col}",
                "reference": ref_col,
                "deletion_start": int(intervals["start"].min()) if not intervals.empty else np.nan,
                "deletion_end": int(intervals["end"].max()) if not intervals.empty else np.nan,
                **stats,
            }
        )

    deviations = []
    upd = all_windows[all_windows["sample"] == "004P"].copy()
    if not upd.empty and not control_ref.empty:
        zz = upd.merge(control_ref[["window_start", "control_maternal_mean"]], on="window_start", how="left")
        for hap, col in [("haplotype1", "mean_meth_haplotype1"), ("haplotype2", "mean_meth_haplotype2")]:
            stats = pearson_summary(zz[col], zz["control_maternal_mean"])
            corr_rows.append(
                {
                    "sample": "004P",
                    "group": "PWS_mUPD",
                    "comparison": f"004P_{hap}_vs_control_maternal",
                    "reference": "control_maternal_mean",
                    "deletion_start": np.nan,
                    "deletion_end": np.nan,
                    **stats,
                }
            )
            diff = zz[col] - zz["control_maternal_mean"]
            hit = zz[diff.abs() >= deviation_threshold].copy()
            if not hit.empty:
                hit["haplotype"] = hap
                hit["methylation"] = hit[col]
                hit["control_maternal_mean"] = hit["control_maternal_mean"]
                hit["delta_from_control_maternal"] = diff.loc[hit.index]
                deviations.append(
                    hit[
                        [
                            "sample",
                            "haplotype",
                            "window_start",
                            "window_end",
                            "methylation",
                            "control_maternal_mean",
                            "delta_from_control_maternal",
                        ]
                    ]
                )

    corr = pd.DataFrame(corr_rows)
    dev = pd.concat(deviations, ignore_index=True) if deviations else pd.DataFrame(
        columns=[
            "sample",
            "haplotype",
            "window_start",
            "window_end",
            "methylation",
            "control_maternal_mean",
            "delta_from_control_maternal",
        ]
    )
    return corr, dev


def reciprocal_delta_and_boundaries(
    retained: pd.DataFrame,
    spec: WindowSpec,
    genes: pd.DataFrame,
    icrs: pd.DataFrame,
    segdup_blocks: pd.DataFrame,
    threshold: float = 0.10,
    smooth_window: int = 31,
    min_bp: int = 10_000,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if retained.empty:
        return pd.DataFrame(), pd.DataFrame()

    base = spec.frame()

    def summarize_group(group: str, mean_col: str, n_col: str) -> pd.DataFrame:
        z = retained[retained["group"] == group]
        if z.empty:
            return base[["window_start", "window_end", "window_mid"]].assign(**{mean_col: np.nan, n_col: 0})
        out = (
            z.groupby(["window_start", "window_end", "window_mid"])["retained_methylation"]
            .agg([("mean_value", "mean"), ("n_value", "count")])
            .reset_index()
            .rename(columns={"mean_value": mean_col, "n_value": n_col})
        )
        return base.merge(out, on=["window_start", "window_end", "window_mid"], how="left")

    pws = summarize_group("PWS_DEL", "pwsdel_retained_maternal_mean", "n_pwsdel_samples")
    asdel = summarize_group("AS_DEL", "asdel_retained_paternal_mean", "n_asdel_samples")
    profile = pws.merge(
        asdel[["window_start", "asdel_retained_paternal_mean", "n_asdel_samples"]],
        on="window_start",
        how="left",
    )
    profile = profile.sort_values("window_start").reset_index(drop=True)
    profile["reciprocal_delta"] = profile["pwsdel_retained_maternal_mean"] - profile["asdel_retained_paternal_mean"]
    profile["reciprocal_delta_smoothed"] = smooth_values(profile["reciprocal_delta"], window=smooth_window)
    profile["abs_reciprocal_delta_smoothed"] = profile["reciprocal_delta_smoothed"].abs()
    profile["boundary_threshold"] = threshold
    profile["delta_state"] = "low_delta"
    profile.loc[profile["reciprocal_delta_smoothed"] >= threshold, "delta_state"] = "maternal_higher"
    profile.loc[profile["reciprocal_delta_smoothed"] <= -threshold, "delta_state"] = "paternal_higher"

    gene_intervals = genes[["gene", "start", "end", "strand"]].copy() if genes is not None and not genes.empty else pd.DataFrame()
    icr_intervals = icrs.rename(columns={"name": "icr_name"}).copy() if icrs is not None and not icrs.empty else pd.DataFrame()
    segdup_intervals = segdup_blocks.copy() if segdup_blocks is not None and not segdup_blocks.empty else pd.DataFrame()
    if not segdup_intervals.empty:
        segdup_intervals["segdup_block"] = segdup_intervals["start"].astype(str) + "-" + segdup_intervals["end"].astype(str)

    rows = []
    for state in ["maternal_higher", "paternal_higher"]:
        mask = (profile["delta_state"] == state).fillna(False).to_numpy()
        for start_i, end_i in boolean_runs(mask):
            start = int(profile.iloc[start_i]["window_start"])
            end = int(profile.iloc[end_i]["window_end"])
            if end - start < min_bp:
                continue
            z = profile.iloc[start_i : end_i + 1]
            row = {
                "chrom": spec.chrom,
                "start": start,
                "end": end,
                "size_bp": end - start,
                "delta_state": state,
                "interpretation": "PWS_DEL_maternal_higher_than_AS_DEL_paternal"
                if state == "maternal_higher"
                else "AS_DEL_paternal_higher_than_PWS_DEL_maternal",
                "n_windows": int(end_i - start_i + 1),
                "mean_delta_smoothed": float(z["reciprocal_delta_smoothed"].mean()),
                "max_abs_delta_smoothed": float(z["abs_reciprocal_delta_smoothed"].max()),
                "threshold_abs_delta": threshold,
                "smoothing_windows": smooth_window,
            }
            for prefix, pos in [("left_boundary", start), ("right_boundary", end)]:
                gene, gene_dist, gene_interval = nearest_interval(pos, gene_intervals, "gene")
                icr, icr_dist, icr_interval = nearest_interval(pos, icr_intervals, "icr_name")
                segdup, segdup_dist, segdup_interval = nearest_interval(pos, segdup_intervals, "segdup_block")
                row.update(
                    {
                        f"{prefix}_position": pos,
                        f"{prefix}_nearest_gene": gene,
                        f"{prefix}_distance_to_nearest_gene_bp": gene_dist,
                        f"{prefix}_nearest_gene_interval": gene_interval,
                        f"{prefix}_nearest_ICR": icr,
                        f"{prefix}_distance_to_nearest_ICR_bp": icr_dist,
                        f"{prefix}_nearest_ICR_interval": icr_interval,
                        f"{prefix}_nearest_or_overlapping_segdup_block": segdup,
                        f"{prefix}_distance_to_segdup_block_bp": segdup_dist,
                        f"{prefix}_segdup_block_interval": segdup_interval,
                    }
                )
            rows.append(row)

    boundaries = pd.DataFrame(rows)
    if not boundaries.empty:
        boundaries = boundaries.sort_values(["start", "end", "delta_state"]).reset_index(drop=True)
    return profile, boundaries


def classify_control_regions(
    control_ref: pd.DataFrame,
    divergent_threshold: float,
    convergent_threshold: float,
    min_bp: int,
) -> pd.DataFrame:
    if control_ref.empty:
        return pd.DataFrame()
    rows = []
    z = control_ref.sort_values("window_start").reset_index(drop=True)
    for label, mask in [
        ("imprinted_divergent", (z["abs_delta"] >= divergent_threshold).fillna(False).to_numpy()),
        ("biallelic_convergent", (z["abs_delta"] <= convergent_threshold).fillna(False).to_numpy()),
    ]:
        for start_i, end_i in boolean_runs(mask):
            start = int(z.iloc[start_i]["window_start"])
            end = int(z.iloc[end_i]["window_end"])
            if end - start < min_bp:
                continue
            rows.append(
                {
                    "region_class": label,
                    "chrom": PHASE2_CHROM,
                    "start": start,
                    "end": end,
                    "size_bp": end - start,
                    "mean_abs_delta": float(z.iloc[start_i : end_i + 1]["abs_delta"].mean()),
                    "n_windows": int(end_i - start_i + 1),
                }
            )
    return pd.DataFrame(rows).sort_values(["region_class", "start"]).reset_index(drop=True) if rows else pd.DataFrame()


def feature_from_genes(genes: pd.DataFrame, name: str, patterns: List[str], label: Optional[str] = None) -> Optional[dict]:
    hits = []
    for pat in patterns:
        if pat.endswith("*"):
            prefix = pat[:-1].upper()
            z = genes[genes["gene"].str.upper().str.startswith(prefix)]
        else:
            z = genes[genes["gene"].str.upper() == pat.upper()]
        if not z.empty:
            hits.append(z)
    if not hits:
        return None
    h = pd.concat(hits, ignore_index=True)
    return {
        "feature": label or name,
        "chrom": h["chrom"].iloc[0],
        "start": int(h["start"].min()),
        "end": int(h["end"].max()),
        "strand": h["strand"].iloc[0] if h["strand"].nunique() == 1 else ".",
        "source": "T2T_GTF",
    }


def build_phase2_gene_features(genes: pd.DataFrame, pws_ic_start: int, pws_ic_end: int, spec: WindowSpec) -> pd.DataFrame:
    specs = [
        ("MKRN3", ["MKRN3"], "MKRN3"),
        ("MAGEL2", ["MAGEL2"], "MAGEL2"),
        ("NDN", ["NDN"], "NDN"),
        ("NPAP1", ["NPAP1"], "NPAP1"),
        ("SNURF-SNRPN", ["SNURF", "SNRPN"], "SNURF-SNRPN"),
        ("SNORD116 cluster", ["SNORD116*"], "SNORD116"),
        ("SNORD115 cluster", ["SNORD115*"], "SNORD115"),
        ("UBE3A-ATS", ["SNHG14", "UBE3A-AS1"], "UBE3A-ATS"),
        ("UBE3A", ["UBE3A"], "UBE3A"),
    ]
    rows = [
        {
            "feature": "PWS-IC",
            "chrom": PHASE2_CHROM,
            "start": pws_ic_start,
            "end": pws_ic_end,
            "strand": ".",
            "source": "SNRPN_TSS_proxy",
        },
        {
            "feature": "AS-IC",
            "chrom": PHASE2_CHROM,
            "start": AS_IC_CENTER - 10_000,
            "end": AS_IC_CENTER + 10_000,
            "strand": ".",
            "source": "prompt_coordinate_approximation",
        },
    ]
    for name, patterns, label in specs:
        f = feature_from_genes(genes, name, patterns, label)
        if f is not None:
            rows.append(f)
    out = pd.DataFrame(rows)
    out = out.drop_duplicates("feature", keep="first")
    out = assign_plot_coordinates(out, spec, width=9_000)
    left = out["end"] <= spec.start
    right = out["start"] >= spec.end
    out = out[out["within_phase2_window"] | left | right].copy()
    return out.sort_values(["start", "end"]).reset_index(drop=True)


def add_region_spans(ax, regions: pd.DataFrame, y_text: float = 1.02) -> None:
    if regions.empty:
        return
    div = regions[regions["region_class"] == "imprinted_divergent"].copy()
    if div.empty:
        return
    div = div.sort_values("size_bp", ascending=False).head(3)
    for i, (_, r) in enumerate(div.iterrows()):
        ax.axvspan(r["start"], r["end"], color="#f7c6c7", alpha=0.24, lw=0)
        if i == 0:
            x = (r["start"] + r["end"]) / 2
            ax.text(x, y_text, "imprinted divergence", ha="center", va="bottom", fontsize=7, color="#8a1c2b")


def add_feature_bands(ax, features: pd.DataFrame) -> None:
    if features.empty:
        return
    for _, r in features.iterrows():
        name = str(r["feature"])
        start = int(r["plot_start"]) if "plot_start" in r and pd.notna(r["plot_start"]) else int(r["start"])
        end = int(r["plot_end"]) if "plot_end" in r and pd.notna(r["plot_end"]) else int(r["end"])
        if end <= start:
            continue
        color = "#e7b3c2" if "IC" in name else "#b8c7d9"
        alpha = 0.09 if "IC" in name else 0.035
        ax.axvspan(start, end, color=color, alpha=alpha, lw=0, zorder=0)


def format_mb_axis(ax) -> None:
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x / 1e6:.1f}"))


def smooth_values(values: pd.Series, window: int = 31) -> pd.Series:
    return pd.to_numeric(values, errors="coerce").rolling(window, center=True, min_periods=max(3, window // 5)).median()


def plot_smoothed_track(
    ax,
    x: pd.Series,
    y: pd.Series,
    color: str,
    label: Optional[str] = None,
    lw: float = 1.4,
    alpha: float = 0.95,
    window: int = 31,
    ls: str = "-",
) -> None:
    yy = smooth_values(y, window=window)
    ax.plot(x, yy, color=color, lw=lw, alpha=alpha, label=label, ls=ls)


def plot_breakpoint_lines(ax, intervals: pd.DataFrame, group: Optional[str], color: str) -> None:
    if intervals.empty:
        return
    z = intervals if group is None else intervals[intervals["group"] == group]
    if z.empty:
        return
    positions = sorted(set(z["start"].astype(int).tolist() + z["end"].astype(int).tolist()))
    for pos in positions:
        ax.axvline(pos, color=color, lw=0.55, ls="--", alpha=0.35, zorder=1)


def unique_breakpoint_positions(intervals: pd.DataFrame, group: str, spec: WindowSpec) -> List[int]:
    if intervals.empty:
        return []
    z = intervals[intervals["group"] == group]
    if z.empty:
        return []
    positions = sorted(set(z["start"].astype(int).tolist() + z["end"].astype(int).tolist()))
    return [p for p in positions if spec.start < p < spec.end]


def bp_hotspots_table(spec: WindowSpec) -> pd.DataFrame:
    rows = []
    for name, pos in BP_HOTSPOTS_T2T:
        rows.append(
            {
                "name": name,
                "chrom": spec.chrom,
                "position": pos,
                "position_mb": pos / 1e6,
                "in_plotted_range": spec.start <= pos <= spec.end,
                "source_script": str(DEFAULT_BP_SOURCE_SCRIPT),
                "source_variable": "PWS_BREAKPOINTS",
            }
        )
    return pd.DataFrame(rows)


def bp_cluster_table(spec: WindowSpec) -> pd.DataFrame:
    rows = []
    for cluster in BP_CLUSTER_INTERVALS_T2T:
        start = int(cluster["start"])
        end = int(cluster["end"])
        rows.append(
            {
                "name": cluster["name"],
                "chrom": spec.chrom,
                "start": start,
                "end": end,
                "start_mb": start / 1e6,
                "end_mb": end / 1e6,
                "approx_size": cluster["approx_size"],
                "interpretation": cluster["interpretation"],
                "overlaps_plotted_range": start < spec.end and end > spec.start,
            }
        )
    return pd.DataFrame(rows)


def snrpn_marker_position(gene_features: pd.DataFrame, spec: WindowSpec) -> int:
    if gene_features is not None and not gene_features.empty:
        for feature in ["SNURF-SNRPN", "PWS-IC"]:
            z = gene_features[gene_features["feature"].astype(str) == feature]
            if not z.empty:
                return int(pd.to_numeric(z.iloc[0]["start"], errors="coerce"))
    return min(max(22_560_000, spec.start), spec.end)


def add_panel_vertical_annotations(
    ax,
    deletion_intervals: pd.DataFrame,
    snrpn_pos: int,
    spec: WindowSpec,
    include_labels: bool = True,
) -> None:
    for cluster in BP_CLUSTER_INTERVALS_T2T:
        cluster_start = int(cluster["start"])
        cluster_end = int(cluster["end"])
        plot_start = max(cluster_start, spec.start)
        plot_end = min(cluster_end, spec.end)
        if plot_start >= plot_end:
            continue
        ax.axvspan(plot_start, plot_end, color="#8B6914", alpha=0.075, lw=0, zorder=0)
        ax.axvline(plot_start, color="#8B6914", lw=0.55, ls=":", alpha=0.55, zorder=1)
        ax.axvline(plot_end, color="#8B6914", lw=0.55, ls=":", alpha=0.55, zorder=1)
        if include_labels:
            ax.text(
                (plot_start + plot_end) / 2,
                -0.135,
                f"{cluster['name']} SD cluster\n{cluster_start / 1e6:.3f}-{cluster_end / 1e6:.3f} Mb",
                transform=ax.get_xaxis_transform(),
                rotation=0,
                ha="center",
                va="top",
                fontsize=5.25,
                color="#6e520f",
                linespacing=0.9,
                clip_on=False,
                bbox=dict(fc="white", ec="none", alpha=0.78, pad=0.5),
            )

    if spec.start <= snrpn_pos <= spec.end:
        ax.axvline(snrpn_pos, color="#1f1f1f", lw=0.8, ls="-", alpha=0.58, zorder=3)
        if include_labels:
            ax.text(
                snrpn_pos,
                -0.135,
                f"SNRPN\n{snrpn_pos / 1e6:.3f} Mb",
                transform=ax.get_xaxis_transform(),
                rotation=0,
                ha="center",
                va="top",
                fontsize=5.25,
                color="#1f1f1f",
                linespacing=0.9,
                clip_on=False,
                bbox=dict(fc="white", ec="none", alpha=0.78, pad=0.5),
            )

    for group, color in [
        ("PWS_DEL", MATERNAL),
        ("AS_DEL", PATERNAL),
    ]:
        positions = unique_breakpoint_positions(deletion_intervals, group, spec)
        if not positions:
            continue
        for pos in positions:
            ax.axvline(pos, color=color, lw=0.6, ls="--", alpha=0.34, zorder=2)


def legend_outside(ax, include_annotations: bool = True, **kwargs) -> None:
    handles, labels = ax.get_legend_handles_labels()
    if include_annotations:
        extra = [
            Line2D([0], [0], color="#1f1f1f", lw=0.8, ls="-", alpha=0.72, label="SNRPN"),
            Line2D([0], [0], color=MATERNAL, lw=0.8, ls="--", alpha=0.72, label="PWS-DEL breakpoints"),
            Line2D([0], [0], color=PATERNAL, lw=0.8, ls="--", alpha=0.72, label="AS-DEL breakpoints"),
            Rectangle((0, 0), 1, 1, facecolor="#8B6914", alpha=0.14, edgecolor="#8B6914", label="BP1/BP2/BP3 T2T SD clusters"),
        ]
        existing = set(labels)
        for h in extra:
            label = h.get_label()
            if label not in existing:
                handles.append(h)
                labels.append(label)
                existing.add(label)
    if not handles:
        return
    defaults = {
        "loc": "upper left",
        "bbox_to_anchor": (1.12, 1.0),
        "borderaxespad": 0.0,
        "frameon": False,
        "fontsize": 6.4,
        "handlelength": 2.1,
    }
    defaults.update(kwargs)
    ax.legend(handles, labels, **defaults)


def add_correlation_box(ax, text: str, y: float = 0.02) -> None:
    if not text:
        return
    ax.text(
        1.12,
        0.02 if y is None else y,
        text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=6.5,
        color="#222222",
        bbox=dict(fc="white", ec="#cfcfcf", lw=0.35, alpha=0.82, pad=2.0),
        clip_on=False,
    )


def correlation_label(a: pd.Series, b: pd.Series, prefix: str) -> str:
    stats = pearson_summary(a, b)
    r = stats["pearson_r"]
    n = stats["n_windows"]
    if not np.isfinite(r):
        return f"{prefix}\nr=NA, n={n}"
    return f"{prefix}\nr={r:.3f}, n={n}"


def retained_group_correlation_label(
    retained: pd.DataFrame,
    control_ref: pd.DataFrame,
    group: str,
    deletion_intervals: pd.DataFrame,
) -> str:
    if retained.empty or control_ref.empty:
        return ""
    ref_col = "control_maternal_mean" if group == "PWS_DEL" else "control_paternal_mean"
    label = "retained vs control maternal" if group == "PWS_DEL" else "retained vs control paternal"
    values = []
    for _, z in retained[retained["group"] == group].groupby("sample"):
        zz = z.merge(control_ref[["window_start", ref_col]], on="window_start", how="left")
        zz = zz[zz["within_deletion_footprint"]]
        stats = pearson_summary(zz["retained_methylation"], zz[ref_col])
        if np.isfinite(stats["pearson_r"]):
            values.append(stats["pearson_r"])
    if not values:
        return ""
    return f"Pearson {label}\nmean r={np.mean(values):.3f}\nrange {min(values):.3f}-{max(values):.3f}"


def mupd_correlation_label(all_windows: pd.DataFrame, control_ref: pd.DataFrame) -> str:
    if all_windows.empty or control_ref.empty:
        return ""
    upd = all_windows[all_windows["sample"] == "004P"].copy()
    if upd.empty:
        return ""
    zz = upd.merge(control_ref[["window_start", "control_maternal_mean"]], on="window_start", how="left")
    parts = []
    for hap, col in [("hap1", "mean_meth_haplotype1"), ("hap2", "mean_meth_haplotype2")]:
        stats = pearson_summary(zz[col], zz["control_maternal_mean"])
        if np.isfinite(stats["pearson_r"]):
            parts.append(f"{hap} r={stats['pearson_r']:.3f}")
    return "Pearson vs control maternal\n" + "\n".join(parts) if parts else ""


def reciprocal_overlay_correlation_label(retained: pd.DataFrame) -> str:
    if retained.empty:
        return ""
    pws = (
        retained[retained["group"] == "PWS_DEL"]
        .groupby("window_start")["retained_methylation"]
        .mean()
        .reset_index(name="pws_maternal")
    )
    asdel = (
        retained[retained["group"] == "AS_DEL"]
        .groupby("window_start")["retained_methylation"]
        .mean()
        .reset_index(name="asdel_paternal")
    )
    z = pws.merge(asdel, on="window_start", how="inner")
    return correlation_label(z["pws_maternal"], z["asdel_paternal"], "Pearson PWS maternal vs AS paternal")


def plot_gene_track(
    ax,
    features: pd.DataFrame,
    spec: WindowSpec,
    icrs: Optional[pd.DataFrame] = None,
    segdup_blocks: Optional[pd.DataFrame] = None,
    deletion_intervals: Optional[pd.DataFrame] = None,
) -> None:
    ax.set_xlim(spec.start, spec.end)
    ax.set_ylim(0, 1.08)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="#eeeeee", lw=0.4)
    for y in [0.79, 0.55, 0.34, 0.13]:
        ax.hlines(y, spec.start, spec.end, color="#d6d6d6", lw=0.6)
    for label, y in [("Genes", 0.84), ("ICRs", 0.60), ("Segmental duplications", 0.38), ("CNV breakpoints", 0.18)]:
        ax.text(-0.012, y, label, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7.0, fontweight="bold")

    gene_lanes = {
        "UBE3A-ATS": 0.88,
        "SNURF-SNRPN": 0.72,
        "SNORD116": 0.80,
        "SNORD115": 0.72,
        "UBE3A": 0.80,
    }
    default_levels = [0.73, 0.81, 0.89]
    upstream_gene_labels = []
    regulatory_features = []
    for i, (_, r) in enumerate(features.iterrows()):
        if not bool(r.get("within_phase2_window", True)):
            upstream_gene_labels.append(str(r["feature"]))
            continue
        if "IC" in str(r["feature"]):
            regulatory_features.append(r)
            continue
        y = gene_lanes.get(str(r["feature"]), default_levels[i % len(default_levels)])
        x0 = int(r["plot_start"]) if "plot_start" in r and pd.notna(r["plot_start"]) else int(r["start"])
        x1 = int(r["plot_end"]) if "plot_end" in r and pd.notna(r["plot_end"]) else int(r["end"])
        width = max(1, x1 - x0)
        ax.add_patch(Rectangle((x0, y - 0.020), width, 0.040, fc="#b8c7d9", ec="#333333", lw=0.3))
        ax.text(
            (x0 + x1) / 2,
            y + 0.030,
            str(r["feature"]),
            ha="center",
            va="bottom",
            fontsize=5.2,
            rotation=0,
            bbox=dict(fc="white", ec="none", alpha=0.72, pad=0.5),
        )

    upstream_icr_labels = []
    for i, r in enumerate(regulatory_features):
        x0 = int(r["plot_start"]) if "plot_start" in r and pd.notna(r["plot_start"]) else int(r["start"])
        x1 = int(r["plot_end"]) if "plot_end" in r and pd.notna(r["plot_end"]) else int(r["end"])
        y = 0.60 if str(r["feature"]) == "AS-IC" else 0.52
        ax.add_patch(Rectangle((x0, y - 0.026), max(1, x1 - x0), 0.052, fc="#f2b8c5", ec="#7d3a44", lw=0.3))
        ax.text(
            (x0 + x1) / 2,
            y + 0.038,
            str(r["feature"]),
            ha="center",
            va="bottom",
            fontsize=5.4,
            rotation=28,
            bbox=dict(fc="white", ec="none", alpha=0.74, pad=0.5),
        )

    if icrs is not None and not icrs.empty:
        for i, (_, r) in enumerate(icrs.iterrows()):
            if not bool(r.get("within_phase2_window", True)):
                upstream_icr_labels.append(f"{r['name']} {r['gene']}".strip())
                continue
            x0 = int(r["plot_start"])
            x1 = int(r["plot_end"])
            y = 0.49 + (i % 2) * 0.060
            ax.add_patch(Rectangle((x0, y - 0.025), max(1, x1 - x0), 0.05, fc="#e89aaa", ec="#7d3a44", lw=0.3))
            label = re.sub(r"[*^#]", "", str(r["name"]))
            ax.text(
                (x0 + x1) / 2,
                y + 0.035,
                label,
                ha="center",
                va="bottom",
                fontsize=5.0,
                rotation=28,
                bbox=dict(fc="white", ec="none", alpha=0.72, pad=0.5),
            )

    upstream_bits = []
    if upstream_gene_labels:
        upstream_bits.append("upstream genes: " + ", ".join(upstream_gene_labels))
    if upstream_icr_labels:
        upstream_bits.append("upstream ICRs: " + ", ".join(upstream_icr_labels))
    if upstream_bits:
        ax.text(
            spec.start + 0.004 * (spec.end - spec.start),
            0.985,
            "Outside fixed x-range: " + "; ".join(upstream_bits),
            ha="left",
            va="top",
            fontsize=5.6,
            color="#333333",
            bbox=dict(fc="white", ec="none", alpha=0.80, pad=0.6),
        )

    if segdup_blocks is not None and not segdup_blocks.empty:
        for _, r in segdup_blocks.iterrows():
            ax.add_patch(
                Rectangle(
                    (int(r["start"]), 0.275),
                    max(1, int(r["end"]) - int(r["start"])),
                    0.05,
                    fc="#8f8f8f",
                    ec="none",
                    alpha=0.45,
                )
            )

    if deletion_intervals is not None and not deletion_intervals.empty:
        colors = {"PWS_DEL": MATERNAL, "AS_DEL": PATERNAL}
        for group, z in deletion_intervals.groupby("group"):
            y = 0.09 if group == "PWS_DEL" else 0.15
            for _, r in z.iterrows():
                for pos in [int(r["start"]), int(r["end"])]:
                    ax.vlines(pos, y - 0.035, y + 0.035, color=colors.get(group, GREY), lw=0.65, alpha=0.75)
        ax.text(spec.end, 0.09, "PWS-DEL", ha="right", va="center", fontsize=5.6, color=MATERNAL)
        ax.text(spec.end, 0.15, "AS-DEL", ha="right", va="center", fontsize=5.6, color=PATERNAL)

    ax.text(
        spec.start,
        1.06,
        "T2T annotation tracks: genes, ICRs, segmental duplications, deletion breakpoints",
        ha="left",
        va="top",
        fontsize=8.2,
        fontweight="bold",
    )


def plot_deletion_shading(ax, intervals: pd.DataFrame, group: str) -> None:
    z = intervals[intervals["group"] == group] if not intervals.empty else intervals
    if z.empty:
        return
    for _, r in z.iterrows():
        ax.axvspan(int(r["start"]), int(r["end"]), color="#bdbdbd", alpha=0.10, lw=0)


def plot_phase2_figure(
    outdir: Path,
    spec: WindowSpec,
    all_windows: pd.DataFrame,
    control_ref: pd.DataFrame,
    control_samples: pd.DataFrame,
    retained: pd.DataFrame,
    deletion_intervals: pd.DataFrame,
    regions: pd.DataFrame,
    gene_features: pd.DataFrame,
    reciprocal_delta: Optional[pd.DataFrame] = None,
    reciprocal_boundaries: Optional[pd.DataFrame] = None,
    icrs: Optional[pd.DataFrame] = None,
    segdup_blocks: Optional[pd.DataFrame] = None,
) -> None:
    fig = plt.figure(figsize=(19, 15.8))
    gs = GridSpec(
        5,
        1,
        height_ratios=[1.12, 1.0, 1.0, 0.92, 1.0],
        hspace=0.56,
        figure=fig,
    )
    axes = [fig.add_subplot(gs[i, 0]) for i in range(5)]
    snrpn_pos = snrpn_marker_position(gene_features, spec)

    for ax in axes:
        ax.set_xlim(spec.start, spec.end)
        ax.set_ylim(-0.05, 1.08)
        ax.set_ylabel("Meth.")
        ax.grid(axis="y", color="#ececec", lw=0.5)
        ax.tick_params(axis="x", labelbottom=False)
        ax.spines[["top", "right"]].set_visible(False)
        add_panel_vertical_annotations(ax, deletion_intervals, snrpn_pos, spec)
        format_mb_axis(ax)

    # Panel A: controls, separate lines per haplotype and sample.
    ax = axes[0]
    add_region_spans(ax, regions)
    if not control_samples.empty:
        for sample, z in control_samples.groupby("sample"):
            z = z.sort_values("window_start")
            plot_smoothed_track(ax, z["window_mid"], z["control_maternal"], color=MATERNAL, lw=0.75, alpha=0.28)
            plot_smoothed_track(ax, z["window_mid"], z["control_paternal"], color=PATERNAL, lw=0.75, alpha=0.28)
    if not control_ref.empty:
        plot_smoothed_track(ax, control_ref["window_mid"], control_ref["control_maternal_mean"], color=MATERNAL, lw=2.0, alpha=0.98, label="maternal mean")
        plot_smoothed_track(ax, control_ref["window_mid"], control_ref["control_paternal_mean"], color=PATERNAL, lw=2.0, alpha=0.98, label="paternal mean")
        add_correlation_box(
            ax,
            correlation_label(
                control_ref["control_maternal_mean"],
                control_ref["control_paternal_mean"],
                "Pearson maternal vs paternal",
            ),
        )
    ax.set_title("A. Controls: reciprocal maternal and paternal cis-methylation architecture", loc="left", fontsize=9, fontweight="bold", pad=7)
    legend_outside(ax, fontsize=6.4)

    # Panel B: PWS-DEL retained maternal signal plus coverage track.
    ax = axes[1]
    plot_deletion_shading(ax, deletion_intervals, "PWS_DEL")
    pws = retained[retained["group"] == "PWS_DEL"]
    for sample, z in pws.groupby("sample"):
        z = z.sort_values("window_start")
        plot_smoothed_track(ax, z["window_mid"], z["retained_methylation"], color=MATERNAL, lw=0.65, alpha=0.18)
    if not control_ref.empty:
        plot_smoothed_track(ax, control_ref["window_mid"], control_ref["control_maternal_mean"], color=MATERNAL, lw=1.2, ls=":", alpha=0.85, label="control maternal mean")
    if not pws.empty:
        mean = pws.groupby("window_mid")["retained_methylation"].mean().reset_index()
        plot_smoothed_track(ax, mean["window_mid"], mean["retained_methylation"], color=MATERNAL, lw=2.0, alpha=0.98, label="PWS-DEL retained mean")
    cov_ax = ax.twinx()
    if not pws.empty:
        cov = pws.groupby("window_mid")["prompt_haplotype2_coverage"].mean().reset_index()
        cov_ax.plot(cov["window_mid"], cov["prompt_haplotype2_coverage"], color=GREY, lw=0.65, alpha=0.60)
    cov_ax.set_ylabel("hap2 reads", fontsize=8, color=GREY)
    cov_ax.tick_params(axis="y", labelsize=7, colors=GREY)
    cov_ax.spines["top"].set_visible(False)
    add_correlation_box(ax, retained_group_correlation_label(retained, control_ref, "PWS_DEL", deletion_intervals))
    ax.set_title("B. PWS-DEL: retained maternal-only architecture across the deletion footprint", loc="left", fontsize=9, fontweight="bold", pad=7)
    legend_outside(ax, fontsize=6.4)

    # Panel C: AS-DEL retained paternal signal plus absent-haplotype coverage.
    ax = axes[2]
    plot_deletion_shading(ax, deletion_intervals, "AS_DEL")
    asdel = retained[retained["group"] == "AS_DEL"]
    for sample, z in asdel.groupby("sample"):
        z = z.sort_values("window_start")
        plot_smoothed_track(ax, z["window_mid"], z["retained_methylation"], color=PATERNAL, lw=0.65, alpha=0.18)
    if not control_ref.empty:
        plot_smoothed_track(ax, control_ref["window_mid"], control_ref["control_paternal_mean"], color=PATERNAL, lw=1.2, ls=":", alpha=0.85, label="control paternal mean")
    if not asdel.empty:
        mean = asdel.groupby("window_mid")["retained_methylation"].mean().reset_index()
        plot_smoothed_track(ax, mean["window_mid"], mean["retained_methylation"], color=PATERNAL, lw=2.0, alpha=0.98, label="AS-DEL retained mean")
    cov_ax = ax.twinx()
    if not asdel.empty:
        cov = asdel.groupby("window_mid")["absent_haplotype_coverage"].mean().reset_index()
        cov_ax.plot(cov["window_mid"], cov["absent_haplotype_coverage"], color=GREY, lw=0.65, alpha=0.60)
    cov_ax.set_ylabel("absent-hap reads", fontsize=8, color=GREY)
    cov_ax.tick_params(axis="y", labelsize=7, colors=GREY)
    cov_ax.spines["top"].set_visible(False)
    add_correlation_box(ax, retained_group_correlation_label(retained, control_ref, "AS_DEL", deletion_intervals))
    ax.set_title("C. AS-DEL: retained paternal-only architecture across the deletion footprint", loc="left", fontsize=9, fontweight="bold", pad=7)
    legend_outside(ax, fontsize=6.4)

    # Panel D: PWS-mUPD both haplotypes.
    ax = axes[3]
    upd = all_windows[all_windows["sample"] == "004P"].sort_values("window_start")
    if not upd.empty:
        plot_smoothed_track(ax, upd["window_mid"], upd["mean_meth_haplotype1"], color=MUPD_1, lw=1.35, label="004P haplotype 1")
        plot_smoothed_track(ax, upd["window_mid"], upd["mean_meth_haplotype2"], color=MUPD_2, lw=1.35, label="004P haplotype 2")
    if not control_ref.empty:
        plot_smoothed_track(ax, control_ref["window_mid"], control_ref["control_maternal_mean"], color=MATERNAL, lw=1.1, ls=":", label="control maternal mean")
    add_correlation_box(ax, mupd_correlation_label(all_windows, control_ref))
    ax.set_title("D. PWS-mUPD: both haplotypes test the duplicated maternal architecture", loc="left", fontsize=9, fontweight="bold", pad=7)
    legend_outside(ax, fontsize=6.4)

    # Panel E: reciprocal natural dissection overlay.
    ax = axes[4]
    plot_deletion_shading(ax, deletion_intervals, "PWS_DEL")
    plot_deletion_shading(ax, deletion_intervals, "AS_DEL")
    overlay_rows = []
    for group, color, label in [
        ("PWS_DEL", MATERNAL, "PWS-DEL retained maternal"),
        ("AS_DEL", PATERNAL, "AS-DEL retained paternal"),
    ]:
        z = retained[retained["group"] == group]
        if z.empty:
            continue
        mean = z.groupby("window_mid")["retained_methylation"].mean().reset_index()
        plot_smoothed_track(ax, mean["window_mid"], mean["retained_methylation"], color=color, lw=2.0, label=label)
    ax.set_title("E. Reciprocal natural dissection: maternal-only and paternal-only signals overlaid", loc="left", fontsize=9, fontweight="bold", pad=7)
    add_correlation_box(ax, reciprocal_overlay_correlation_label(retained))
    legend_outside(ax, fontsize=6.4)
    ax.tick_params(axis="x", labelbottom=True)
    ax.set_xlabel(f"{spec.chrom} coordinate (Mb)", labelpad=34)

    fig.suptitle("Figure 2. Reciprocal cis-architecture across the chr15 imprinted domain", y=0.986, fontsize=13, fontweight="bold")
    fig.text(
        0.07,
        0.955,
        "T2T breakpoint SD clusters shaded in all panels: BP1=17.691-20.454 Mb, BP2=20.754-21.184 Mb, BP3=25.876-26.633 Mb.",
        ha="left",
        va="top",
        fontsize=7.2,
        color="#333333",
    )
    fig.subplots_adjust(top=0.925, bottom=0.125, left=0.07, right=0.70)
    fig.savefig(outdir / "figures" / "Figure2_reciprocal_cis_architecture.png", dpi=350, bbox_inches="tight")
    fig.savefig(outdir / "figures" / "Figure2_reciprocal_cis_architecture.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase 2 reciprocal cis-architecture analysis")
    parser.add_argument("--methylation-dir", "--meth-dir", dest="meth_dir", default=DEFAULT_METHYLATION_DIR, type=Path)
    parser.add_argument("--metadata", default=DEFAULT_METADATA_PATH, type=Path)
    parser.add_argument("--gtf", default=DEFAULT_GTF, type=Path)
    parser.add_argument("--cnv-dir", default=DEFAULT_CNV_DIR, type=Path)
    parser.add_argument("--icr-bed", default=DEFAULT_ICR_BED, type=Path, help="T2T imprinting-control-region BED annotations.")
    parser.add_argument("--segdup-bed", default=DEFAULT_SEG_DUP, type=Path, help="T2T segmental-duplication BED-like annotation table.")
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR, type=Path)
    parser.add_argument("--chrom", default=PHASE2_CHROM)
    parser.add_argument("--start", default=PHASE2_START, type=int)
    parser.add_argument("--end", default=PHASE2_END, type=int)
    parser.add_argument("--window-size", default=PHASE2_WINDOW, type=int)
    parser.add_argument("--samples", default=",".join(DEFAULT_PHASE2_SAMPLES),
                        help="Comma-separated sample list for Phase 2; defaults to PWS-DEL, PWS-mUPD, AS-DEL, and controls only.")
    parser.add_argument("--exclude", default=",".join(sorted(DEFAULT_EXCLUDE)))
    parser.add_argument("--methyl-format", default="auto", choices=["auto", "bed", "dss"])
    parser.add_argument("--bed-meth-col", default="9", help="1-based methylation column; pb-CpG BED default is 9.")
    parser.add_argument("--bed-cov-col", default="6", help="1-based read coverage column; pb-CpG BED default is 6.")
    parser.add_argument("--chunksize", default=1_000_000, type=int)
    parser.add_argument("--min-reads", default=10, type=int)
    parser.add_argument("--min-cpgs", default=3, type=int)
    parser.add_argument("--deletion-copy-threshold", default=1.5, type=float)
    parser.add_argument("--divergent-delta", default=0.20, type=float)
    parser.add_argument("--convergent-delta", default=0.10, type=float)
    parser.add_argument("--min-region-bp", default=25_000, type=int)
    parser.add_argument("--deviation-threshold", default=0.25, type=float)
    parser.add_argument("--reciprocal-delta-threshold", default=0.10, type=float)
    parser.add_argument("--reciprocal-smooth-windows", default=31, type=int)
    parser.add_argument("--reciprocal-min-boundary-bp", default=10_000, type=int)
    args = parser.parse_args()

    outdir = mkdir(args.outdir)
    tbldir = mkdir(outdir / "tables")
    figdir = mkdir(outdir / "figures")
    spec = WindowSpec(args.chrom, args.start, args.end, args.window_size)

    exclude = {x.strip() for x in args.exclude.split(",") if x.strip()}

    log("Loading metadata and methylation file inventory")
    load_metadata_table(args.metadata, outdir=outdir)
    file_table = discover_methylation_files(args.meth_dir, exclude)
    file_table["group"] = file_table["sample"].apply(lambda s: sample_group(s, DEFAULT_GROUPS))
    requested_samples = [x.strip() for x in args.samples.split(",") if x.strip()]
    samples = phase_samples_from_metadata(file_table, requested=requested_samples)
    file_table = file_table[file_table["sample"].isin(samples)].copy()
    file_table.to_csv(tbldir / "Phase2_input_methylation_files.tsv", sep="\t", index=False)

    log("Loading T2T gene annotations and IC proxy")
    genes = load_gtf_genes(args.gtf, args.chrom)
    ic_chrom, ic_start, ic_end, ic_gene, _ = infer_ic_region_from_gtf(genes, IC_PRIORITY_GENES, flank=5_000)
    if ic_chrom != args.chrom:
        raise ValueError(f"IC inferred on {ic_chrom}, but Phase 2 is running on {args.chrom}")

    gene_features = build_phase2_gene_features(genes, ic_start, ic_end, spec)
    gene_features.to_csv(tbldir / "Phase2_gene_track_features.tsv", sep="\t", index=False)
    icr_annotations = read_icr_annotations(args.icr_bed, spec)
    if not icr_annotations.empty:
        icr_annotations.to_csv(tbldir / "Phase2_T2T_ICR_annotations.tsv", sep="\t", index=False)
    segdup_raw, segdup_blocks = read_segdup_annotations(args.segdup_bed, spec)
    if not segdup_raw.empty:
        segdup_raw.to_csv(tbldir / "Phase2_T2T_segmental_duplications_raw.tsv.gz", sep="\t", index=False, compression="gzip")
    if not segdup_blocks.empty:
        segdup_blocks.to_csv(tbldir / "Phase2_T2T_segmental_duplication_blocks.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "feature": "phase2_fixed_domain",
                "chrom": args.chrom,
                "start": args.start,
                "end": args.end,
                "window_size": args.window_size,
            },
            {
                "feature": f"PWS_IC_proxy_{ic_gene}",
                "chrom": ic_chrom,
                "start": ic_start,
                "end": ic_end,
                "window_size": np.nan,
            },
            {
                "feature": "AS_IC_prompt_approx",
                "chrom": args.chrom,
                "start": AS_IC_CENTER - 10_000,
                "end": AS_IC_CENTER + 10_000,
                "window_size": np.nan,
            },
        ]
    ).to_csv(tbldir / "Phase2_fixed_coordinates.tsv", sep="\t", index=False)
    bp_hotspots_table(spec).to_csv(tbldir / "Phase2_T2T_BP_hotspots.tsv", sep="\t", index=False)
    bp_cluster_table(spec).to_csv(tbldir / "Phase2_T2T_BP_cluster_annotations.tsv", sep="\t", index=False)

    log("Assigning parental labels from IC methylation and group expectations")
    ic_summary = summarize_ic_for_labels(file_table, args.chrom, ic_start, ic_end, args)
    ic_summary.to_csv(tbldir / "Phase2_IC_haplotype_methylation_for_labeling.tsv", sep="\t", index=False)
    labels = labels_from_ic_and_group(ic_summary, min_cpgs=args.min_cpgs)

    log("Computing 1 kb per-haplotype windows")
    window_tables = []
    for i, sample in enumerate(samples, start=1):
        log(f"Windowing sample {i}/{len(samples)}: {sample}")
        window_tables.append(summarize_sample_windows(sample, file_table, spec, labels, args))
    all_windows = pd.concat(window_tables, ignore_index=True) if window_tables else pd.DataFrame()

    log("Loading CNV deletion footprints")
    cnv_table = discover_cnv_files(args.cnv_dir, samples)
    cnv_table.to_csv(tbldir / "Phase2_input_cnv_files.tsv", sep="\t", index=False)
    deletion_intervals = read_cnv_deletion_intervals(cnv_table, spec, args.deletion_copy_threshold)
    deletion_intervals = deletion_intervals[deletion_intervals["group"].isin(["PWS_DEL", "AS_DEL"])].copy()
    if deletion_intervals.empty:
        log("No CNV deletion intervals found in fixed domain; using haplotype coverage fallback.")
        deletion_intervals = infer_deletion_intervals_from_windows(all_windows, spec, min_reads=args.min_reads)
    deletion_intervals.to_csv(tbldir / "Phase2_deletion_footprints.tsv", sep="\t", index=False)
    breakpoint_annotations = annotate_deletion_breakpoints(deletion_intervals, genes, icr_annotations, segdup_blocks)
    if not breakpoint_annotations.empty:
        breakpoint_annotations.to_csv(tbldir / "Phase2_deletion_breakpoints_T2T_annotated.tsv", sep="\t", index=False)

    labels = update_del_labels_by_coverage(labels, all_windows, deletion_intervals)
    all_windows = apply_labels_to_windows(all_windows, labels)
    write_per_sample_tables(all_windows, outdir)
    all_windows.to_csv(
        tbldir / "Phase2_all_samples_haplotype_1kb_methylation.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
        na_rep="NaN",
    )
    pd.DataFrame(
        [
            {"sample": sample, "haplotype1_label": lab.get("hap1"), "haplotype2_label": lab.get("hap2")}
            for sample, lab in sorted(labels.items())
            if sample in samples
        ]
    ).to_csv(tbldir / "Phase2_haplotype_parental_labels.tsv", sep="\t", index=False)

    log("Building control reference and retained-haplotype profiles")
    ref, control_samples = control_reference(all_windows)
    ref.to_csv(tbldir / "Phase2_control_reference_architecture.tsv", sep="\t", index=False, na_rep="NaN")
    control_samples.to_csv(tbldir / "Phase2_control_sample_architecture.tsv", sep="\t", index=False, na_rep="NaN")
    regions = classify_control_regions(ref, args.divergent_delta, args.convergent_delta, args.min_region_bp)
    regions.to_csv(tbldir / "Phase2_control_divergence_convergence_regions.tsv", sep="\t", index=False)

    retained = retained_profiles(all_windows, deletion_intervals)
    retained.to_csv(tbldir / "Phase2_retained_haplotype_profiles.tsv", sep="\t", index=False, na_rep="NaN")
    corr, dev = compute_correlations(retained, all_windows, ref, deletion_intervals, args.deviation_threshold)
    corr.to_csv(tbldir / "Phase2_retained_and_mUPD_correlations.tsv", sep="\t", index=False, na_rep="NaN")
    dev.to_csv(tbldir / "Phase2_004P_candidate_postzygotic_deviation_windows.tsv", sep="\t", index=False, na_rep="NaN")

    reciprocal_delta, reciprocal_boundaries = reciprocal_delta_and_boundaries(
        retained,
        spec,
        genes,
        icr_annotations,
        segdup_blocks,
        threshold=args.reciprocal_delta_threshold,
        smooth_window=args.reciprocal_smooth_windows,
        min_bp=args.reciprocal_min_boundary_bp,
    )
    reciprocal_delta.to_csv(
        tbldir / "Phase2_reciprocal_delta_profile.tsv",
        sep="\t",
        index=False,
        na_rep="NaN",
    )
    reciprocal_boundaries.to_csv(
        tbldir / "Phase2_reciprocal_delta_boundary_candidates.tsv",
        sep="\t",
        index=False,
        na_rep="NaN",
    )

    log("Rendering Figure 2 composite")
    plot_phase2_figure(
        outdir,
        spec,
        all_windows,
        ref,
        control_samples,
        retained,
        deletion_intervals,
        regions,
        gene_features,
        reciprocal_delta=reciprocal_delta,
        reciprocal_boundaries=reciprocal_boundaries,
        icrs=icr_annotations,
        segdup_blocks=segdup_blocks,
    )

    pd.DataFrame(
        [
            {
                "phase": "Phase2_reciprocal_cis_architecture",
                "domain": f"{args.chrom}:{args.start}-{args.end}",
                "window_size": args.window_size,
                "n_samples": len(samples),
                "n_windows_per_sample": int((args.end - args.start) / args.window_size),
                "min_reads_per_haplotype": args.min_reads,
                "min_cpgs_per_haplotype": args.min_cpgs,
                "figure": str(outdir / "figures" / "Figure2_reciprocal_cis_architecture.png"),
                "all_windows_table": str(tbldir / "Phase2_all_samples_haplotype_1kb_methylation.tsv.gz"),
                "reciprocal_delta_profile_table": str(tbldir / "Phase2_reciprocal_delta_profile.tsv"),
                "reciprocal_boundary_candidate_table": str(tbldir / "Phase2_reciprocal_delta_boundary_candidates.tsv"),
            }
        ]
    ).to_csv(outdir / "phase2_run_summary.tsv", sep="\t", index=False)
    log("Done.")


if __name__ == "__main__":
    main()
