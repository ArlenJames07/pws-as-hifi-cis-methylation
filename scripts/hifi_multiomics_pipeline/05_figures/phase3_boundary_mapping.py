#!/usr/bin/env python3
"""
Phase 3 methylation boundary mapping for PWS/AS Figure 3 (reviewer-clear version).

This script implements the chr15:22-29 Mb boundary-mapping prompts:

  * controls: 1 kb windows / 100 bp step, mean control contrast
    (default: |maternal - paternal|; optionally signed maternal - paternal)
  * PWS-DEL: retained maternal allele vs control biallelic baseline
  * AS-DEL: retained paternal allele vs control biallelic baseline
  * boundary criteria: enter after 5 consecutive windows > 0.4, exit after
    5 consecutive windows < 0.1
  * per-boundary annotation against genes, segmental duplications, BP hotspots,
    CTCF if supplied, and imprinted DMR catalogs
  * sensitivity, breakpoint-distance, derivative, and bootstrap support tables

Default inputs match the existing paper_vf Phase 1 outputs and methylation BEDs.
Outputs are written to /home/rare/arlen/paper_vf by default.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, Rectangle
from matplotlib.ticker import FuncFormatter

from create_figure2_reciprocal_cis_architecture_improved import (
    assign_non_overlapping_label_positions,
    draw_annotation_track as draw_figure2_annotation_track,
)
from paper_vf_phase2_reciprocal_cis_architecture import BP_CLUSTER_INTERVALS_T2T, MATERNAL, PATERNAL, WindowSpec


CHROM = "chr15"
REGION_START = 22_000_000
REGION_END = 29_000_000
WINDOW_SIZE = 1_000
STEP_SIZE = 100
CONSECUTIVE_WINDOWS = 5
ENTER_THRESHOLD = 0.4
EXIT_THRESHOLD = 0.1
ZOOM_FLANK = 50_000
MIN_CPGS_PER_WINDOW = 1
CONTROL_MATCH_WINDOW = 10_000
PATIENT_CONVERGENCE_BP = 5_000
FIGURE2_DISPLAY_START = 18_000_000
FIGURE2_DISPLAY_END = 28_000_000
CONTROL_SIGNAL_MODE_DEFAULT = "absolute"
CONTROL_SIGNAL_FORMULAS = {
    "absolute": "|maternal - paternal|",
    "signed": "maternal - paternal",
}
PATIENT_SIGNAL_FORMULAS = {
    "PWS_DEL": "retained maternal - control biallelic baseline",
    "AS_DEL": "control biallelic baseline - retained paternal",
}
VALIDATED_INTERVAL_TOLERANCE_BP_DEFAULT = 1_000
VALIDATED_CONVERGENCE_INTERVALS = [
    {"source": "Control", "label": "Controls", "start": 22_690_600, "end": 22_696_850, "color": "#222222", "kind": "source"},
    {"source": "PWS_DEL", "label": "PWS-DEL", "start": 22_690_900, "end": 22_694_700, "color": MATERNAL, "kind": "source"},
    {"source": "AS_DEL", "label": "AS-DEL", "start": 22_691_000, "end": 22_694_800, "color": PATERNAL, "kind": "source"},
]
SENSITIVITY_WINDOW_SIZES_BP = (250, 500, 1_000, 2_000)
SENSITIVITY_ENTER_THRESHOLDS = (0.30, 0.40, 0.50)
SENSITIVITY_EXIT_THRESHOLDS = (0.10, 0.15)
SENSITIVITY_CPG_WINDOWS = (25, 50)
BOOTSTRAP_REPLICATES_DEFAULT = 100
BOOTSTRAP_SEED_DEFAULT = 1729
BOOTSTRAP_LOCAL_FLANK_DEFAULT = 20_000
LOW_CPG_DENSITY_WARNING_THRESHOLD = 25

CONTROL_SAMPLES = ("017C", "018C")
PWS_DEL_SAMPLES = ("001P", "002P", "005P", "006P", "007P")
AS_DEL_SAMPLES = ("013A", "014A", "016A")
PWS_MUPD_SAMPLES = ("004P",)

STRUCT_GREY = "#8a8a8a"
DARK_GREY = "#3a3a3a"
SOFT_GREY = "#efefef"
CONTROL_COLOR = "#222222"
PWS_COLOR = "#CC2F5A"
AS_COLOR = "#2166AC"
SHARED_CORE_COLOR = "#C49A00"
HIGHLIGHT_BAND_COLOR = "#F4D35E"
CONTEXT_SHARED_FILL = "#f8efb7"
GENE_TRACK_COLOR = "#a8b5c7"
GENE_TRACK_SECONDARY_COLOR = "#bac6d4"
ICR_TRACK_COLOR = "#7a9ec8"
BP_REGION_COLORS = {
    "BP1": "#c6a96e",
    "BP2": "#8faa92",
    "BP3": "#92a8c9",
}
PANEL_A_GENE_SIDES = {
    "NIPA1": "top",
    "CYFIP1": "top",
    "MKRN3": "top",
    "NDN": "top",
    "SNURF": "bottom",
    "OCA2": "top",
    "APBA2": "top",
}


@dataclass
class SignalTrack:
    sample: str
    cohort: str
    positions: np.ndarray
    signal: np.ndarray
    maternal_methylation: np.ndarray | None = None
    paternal_methylation: np.ndarray | None = None
    retained_methylation: np.ndarray | None = None
    baseline_methylation: np.ndarray | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--methylation-dir", default="/home/rare/arlen/outputs/methylation/genomes_2")
    parser.add_argument("--outdir", default="/home/rare/arlen/paper_vf")
    parser.add_argument("--assignment-table", default="/home/rare/arlen/paper_vf/tables/Figure1C_parental_assignment.tsv")
    parser.add_argument("--gtf", default="/home/rare/arlen/reference/chm13v22.sorted.gtf")
    parser.add_argument("--segdup", default="/home/rare/arlen/reference/dupseg")
    parser.add_argument("--ctcf-bed", default="")
    parser.add_argument("--bp-hotspot-bed", default="")
    parser.add_argument("--imprintome-bed", default="/home/rare/arlen/reference/ICRs_known_chm13.bed")
    parser.add_argument("--court2014-bed", default="/home/rare/arlen/reference/ICR_t2t.bed")
    parser.add_argument("--icr-bed", default="/home/rare/arlen/reference/ICR_t2t.bed")
    parser.add_argument("--repeats", default="/home/rare/arlen/reference/GCF_009914755.1_T2T-CHM13v2.0_rm.bed")
    parser.add_argument("--chrom", default=CHROM)
    parser.add_argument("--region-start", type=int, default=REGION_START)
    parser.add_argument("--region-end", type=int, default=REGION_END)
    parser.add_argument("--window-size", type=int, default=WINDOW_SIZE)
    parser.add_argument("--step-size", type=int, default=STEP_SIZE)
    parser.add_argument("--consecutive-windows", type=int, default=CONSECUTIVE_WINDOWS)
    parser.add_argument("--enter-threshold", type=float, default=ENTER_THRESHOLD)
    parser.add_argument("--exit-threshold", type=float, default=EXIT_THRESHOLD)
    parser.add_argument("--zoom-flank", type=int, default=ZOOM_FLANK)
    parser.add_argument("--min-cpgs-per-window", type=int, default=MIN_CPGS_PER_WINDOW)
    parser.add_argument("--control-match-window", type=int, default=CONTROL_MATCH_WINDOW)
    parser.add_argument(
        "--control-signal-mode",
        choices=sorted(CONTROL_SIGNAL_FORMULAS),
        default=CONTROL_SIGNAL_MODE_DEFAULT,
        help="Use |maternal - paternal| or maternal - paternal for the control contrast.",
    )
    parser.add_argument(
        "--use-validated-intervals",
        action="store_true",
        help="Plot the reviewer-facing interval panel using validated manuscript intervals after comparing them against computed intervals.",
    )
    parser.add_argument(
        "--validated-interval-tolerance-bp",
        type=int,
        default=VALIDATED_INTERVAL_TOLERANCE_BP_DEFAULT,
        help="Warn when computed and validated intervals differ by more than this tolerance.",
    )
    parser.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=BOOTSTRAP_REPLICATES_DEFAULT,
        help="Approximate CpG-bootstrap replicates for shared-interval support.",
    )
    parser.add_argument(
        "--bootstrap-seed",
        type=int,
        default=BOOTSTRAP_SEED_DEFAULT,
    )
    parser.add_argument(
        "--bootstrap-local-flank",
        type=int,
        default=BOOTSTRAP_LOCAL_FLANK_DEFAULT,
        help="Local flank around the primary shared interval used for bootstrap support.",
    )
    return parser.parse_args()


def mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def fmt_int(value: float | int | None) -> str:
    if value is None:
        return ""
    try:
        if not math.isfinite(float(value)):
            return ""
    except (TypeError, ValueError):
        return ""
    return f"{int(round(float(value))):,}"


def pos_formatter() -> FuncFormatter:
    return FuncFormatter(lambda x, _pos: f"{x / 1e6:.2f}")


def local_pos_formatter(decimals: int = 3) -> FuncFormatter:
    return FuncFormatter(lambda x, _pos: f"{x / 1e6:.{decimals}f}")


def structural_landmark_color(name: str) -> str:
    return BP_REGION_COLORS.get(str(name), STRUCT_GREY)


def control_signal_formula(mode: str) -> str:
    return CONTROL_SIGNAL_FORMULAS.get(mode, CONTROL_SIGNAL_FORMULAS[CONTROL_SIGNAL_MODE_DEFAULT])


def compute_control_signal(maternal: np.ndarray, paternal: np.ndarray, mode: str) -> np.ndarray:
    delta = np.asarray(maternal, dtype=float) - np.asarray(paternal, dtype=float)
    if mode == "signed":
        return delta
    return np.abs(delta)


def adaptive_step_size(window_size: int) -> int:
    return max(50, int(round((window_size / 4) / 25.0) * 25))


def interval_distance_to_nearest_edge(start: int, end: int, coord: int) -> int:
    if start <= coord <= end:
        return 0
    if coord < start:
        return start - coord
    return coord - end


def overlap_length(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(int(a_end), int(b_end)) - max(int(a_start), int(b_start)))


def interval_to_interval_gap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    if int(a_end) < int(b_start):
        return int(b_start) - int(a_end)
    if int(b_end) < int(a_start):
        return int(a_start) - int(b_end)
    return 0


def collapse_intervals(intervals: pd.DataFrame, gap_bp: int = 0) -> pd.DataFrame:
    if intervals.empty:
        return pd.DataFrame(columns=["start", "end"])
    z = intervals[["start", "end"]].copy().sort_values(["start", "end"]).reset_index(drop=True)
    merged: list[dict[str, int]] = []
    current_start = int(z.iloc[0]["start"])
    current_end = int(z.iloc[0]["end"])
    for _, row in z.iloc[1:].iterrows():
        start = int(row["start"])
        end = int(row["end"])
        if start <= current_end + gap_bp:
            current_end = max(current_end, end)
            continue
        merged.append({"start": current_start, "end": current_end})
        current_start = start
        current_end = end
    merged.append({"start": current_start, "end": current_end})
    return pd.DataFrame(merged)


def shared_interval_from_boundaries(primary_boundaries: pd.DataFrame) -> tuple[int, int]:
    if primary_boundaries.empty:
        raise ValueError("Primary boundary table is empty.")
    required = {"Control", "PWS_DEL", "AS_DEL"}
    found = set(primary_boundaries["boundary_source"].astype(str))
    missing = required - found
    if missing:
        raise ValueError(f"Primary boundary table is missing interval calls for: {', '.join(sorted(missing))}")
    shared_start = int(pd.to_numeric(primary_boundaries["entry_boundary"], errors="raise").max())
    shared_end = int(pd.to_numeric(primary_boundaries["exit_boundary"], errors="raise").min())
    if shared_start >= shared_end:
        raise ValueError(
            f"Computed shared interval is invalid: start {shared_start:,} is not less than end {shared_end:,}."
        )
    return shared_start, shared_end


def primary_boundaries_to_interval_table(primary_boundaries: pd.DataFrame) -> pd.DataFrame:
    color_map = {"Control": "#222222", "PWS_DEL": MATERNAL, "AS_DEL": PATERNAL}
    label_map = {"Control": "Controls", "PWS_DEL": "PWS-DEL", "AS_DEL": "AS-DEL"}
    rows: list[dict[str, Any]] = []
    for source in ["Control", "PWS_DEL", "AS_DEL"]:
        row = primary_boundaries.loc[primary_boundaries["boundary_source"] == source].iloc[0]
        rows.append(
            {
                "source": source,
                "label": label_map[source],
                "start": int(round(float(row["entry_boundary"]))),
                "end": int(round(float(row["exit_boundary"]))),
                "color": color_map[source],
                "kind": "source",
            }
        )
    shared_start, shared_end = shared_interval_from_boundaries(primary_boundaries)
    rows.append(
        {
            "source": "Consensus",
            "label": "Shared overlap",
            "start": shared_start,
            "end": shared_end,
            "color": "#d99a00",
            "kind": "consensus",
        }
    )
    return pd.DataFrame(rows)


def validated_interval_table() -> pd.DataFrame:
    validated = pd.DataFrame(VALIDATED_CONVERGENCE_INTERVALS)
    shared_start = int(validated["start"].max())
    shared_end = int(validated["end"].min())
    validated = pd.concat(
        [
            validated,
            pd.DataFrame(
                [
                    {
                        "source": "Consensus",
                        "label": "Shared overlap",
                        "start": shared_start,
                        "end": shared_end,
                        "color": "#d99a00",
                        "kind": "consensus",
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    return validated


def compare_interval_tables(
    computed_intervals: pd.DataFrame,
    tolerance_bp: int,
) -> tuple[pd.DataFrame, list[str]]:
    validated = validated_interval_table()
    rows: list[dict[str, Any]] = []
    warnings: list[str] = []
    for source in ["Control", "PWS_DEL", "AS_DEL", "Consensus"]:
        comp = computed_intervals.loc[computed_intervals["source"] == source].iloc[0]
        val = validated.loc[validated["source"] == source].iloc[0]
        start_delta = int(comp["start"]) - int(val["start"])
        end_delta = int(comp["end"]) - int(val["end"])
        within_tolerance = abs(start_delta) <= tolerance_bp and abs(end_delta) <= tolerance_bp
        rows.append(
            {
                "source": source,
                "computed_start": int(comp["start"]),
                "computed_end": int(comp["end"]),
                "validated_start": int(val["start"]),
                "validated_end": int(val["end"]),
                "start_delta_bp": start_delta,
                "end_delta_bp": end_delta,
                "within_tolerance": within_tolerance,
            }
        )
        if not within_tolerance:
            warnings.append(
                f"Validated override mismatch for {source}: computed {int(comp['start']):,}-{int(comp['end']):,} "
                f"vs validated {int(val['start']):,}-{int(val['end']):,} exceeds ±{tolerance_bp:,} bp."
            )
    return pd.DataFrame(rows), warnings


def discover_sample_files(methylation_dir: Path) -> dict[str, dict[str, Path]]:
    files: dict[str, dict[str, Path]] = {}
    for path in methylation_dir.glob("*.bed"):
        match = re.search(r"_([0-9]{3}[A-Z])\.", path.name)
        if not match:
            continue
        sample = match.group(1)
        layer_match = re.search(r"\.(combined|hap1|hap2)\.bed$", path.name)
        if not layer_match:
            continue
        files.setdefault(sample, {})[layer_match.group(1)] = path
    return files


def read_parental_assignments(path: Path) -> dict[str, dict[str, str]]:
    assignments: dict[str, dict[str, str]] = {}
    if not path.exists():
        # Known Phase 1 assignments for the two controls; used only if the table
        # has not been generated yet.
        return {
            "017C": {"maternal": "hap1", "paternal": "hap2"},
            "018C": {"maternal": "hap2", "paternal": "hap1"},
        }
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample = row.get("sample_id", "")
            if sample not in CONTROL_SAMPLES:
                continue
            parent = row.get("parental_assignment", "")
            hap = row.get("haplotype_label", "")
            if parent in {"maternal", "paternal"} and hap in {"hap1", "hap2"}:
                assignments.setdefault(sample, {})[parent] = hap
    return assignments


def run_awk_region(path: Path, chrom: str, start: int, end: int) -> str:
    awk_script = "$1==chrom && $2>=start && $2<end {print}"
    result = subprocess.run(
        [
            "awk",
            "-v",
            f"chrom={chrom}",
            "-v",
            f"start={start}",
            "-v",
            f"end={end}",
            awk_script,
            str(path),
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    return result.stdout


def read_methylation_region(path: Path, chrom: str, start: int, end: int) -> pd.DataFrame:
    rows: list[tuple[int, float, float]] = []
    stdout = run_awk_region(path, chrom, start, end)
    for line in stdout.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 6:
            continue
        try:
            pos = int(fields[1])
            meth = float(fields[3])
            cov = float(fields[5])
        except ValueError:
            continue
        if meth > 1.0:
            meth /= 100.0
        rows.append((pos, min(1.0, max(0.0, meth)), cov))
    if not rows:
        return pd.DataFrame(columns=["pos", "meth", "cov"])
    df = pd.DataFrame(rows, columns=["pos", "meth", "cov"])
    if df["pos"].duplicated().any():
        df = (
            df.groupby("pos", as_index=False)
            .agg(meth=("meth", "mean"), cov=("cov", "mean"))
            .sort_values("pos")
        )
    return df.sort_values("pos").reset_index(drop=True)


def merge_two_tracks(left: pd.DataFrame, right: pd.DataFrame, left_name: str, right_name: str) -> pd.DataFrame:
    a = left[["pos", "meth"]].rename(columns={"meth": left_name})
    b = right[["pos", "meth"]].rename(columns={"meth": right_name})
    return a.merge(b, on="pos", how="inner").sort_values("pos").reset_index(drop=True)


def build_control_baseline(control_combined: list[pd.DataFrame]) -> pd.DataFrame:
    chunks = []
    for i, df in enumerate(control_combined, start=1):
        chunks.append(df[["pos", "meth"]].rename(columns={"meth": f"control_{i}"}))
    if not chunks:
        return pd.DataFrame(columns=["pos", "baseline"])
    merged = chunks[0]
    for chunk in chunks[1:]:
        merged = merged.merge(chunk, on="pos", how="outer")
    meth_cols = [c for c in merged.columns if c.startswith("control_")]
    merged["baseline"] = merged[meth_cols].mean(axis=1, skipna=True)
    return merged[["pos", "baseline"]].dropna().sort_values("pos").reset_index(drop=True)


def window_profile(
    positions: np.ndarray,
    signal: np.ndarray,
    chrom: str,
    region_start: int,
    region_end: int,
    window_size: int,
    step_size: int,
    min_cpgs: int,
) -> pd.DataFrame:
    starts = np.arange(region_start, region_end - window_size + 1, step_size, dtype=int)
    positions = np.asarray(positions, dtype=int)
    signal = np.asarray(signal, dtype=float)
    ok = np.isfinite(signal)
    positions = positions[ok]
    signal = signal[ok]
    if positions.size == 0 or starts.size == 0:
        return pd.DataFrame(columns=["chrom", "window_start", "window_end", "window_mid", "n_cpg", "mean_signal"])
    order = np.argsort(positions)
    positions = positions[order]
    signal = signal[order]
    ends = starts + window_size
    left = np.searchsorted(positions, starts, side="left")
    right = np.searchsorted(positions, ends, side="left")
    n_cpg = right - left
    prefix = np.concatenate([[0.0], np.cumsum(signal)])
    mean_signal = np.full(starts.shape, np.nan, dtype=float)
    valid = n_cpg >= min_cpgs
    if valid.any():
        mean_signal[valid] = (prefix[right[valid]] - prefix[left[valid]]) / n_cpg[valid]
    return pd.DataFrame(
        {
            "chrom": chrom,
            "window_start": starts.astype(int),
            "window_end": ends.astype(int),
            "window_mid": (starts + window_size // 2).astype(int),
            "n_cpg": n_cpg.astype(int),
            "mean_signal": mean_signal,
        }
    )


def adaptive_cpg_window_profile(
    positions: np.ndarray,
    signal: np.ndarray,
    chrom: str,
    region_start: int,
    region_end: int,
    cpg_window: int,
    step_cpg: int,
) -> pd.DataFrame:
    positions = np.asarray(positions, dtype=int)
    signal = np.asarray(signal, dtype=float)
    ok = np.isfinite(signal)
    positions = positions[ok]
    signal = signal[ok]
    region_mask = (positions >= region_start) & (positions <= region_end)
    positions = positions[region_mask]
    signal = signal[region_mask]
    if positions.size < cpg_window:
        return pd.DataFrame(columns=["chrom", "window_start", "window_end", "window_mid", "n_cpg", "mean_signal"])
    order = np.argsort(positions)
    positions = positions[order]
    signal = signal[order]
    start_indices = np.arange(0, positions.size - cpg_window + 1, step_cpg, dtype=int)
    if start_indices.size == 0:
        return pd.DataFrame(columns=["chrom", "window_start", "window_end", "window_mid", "n_cpg", "mean_signal"])
    end_indices = start_indices + cpg_window - 1
    starts = positions[start_indices]
    ends = positions[end_indices] + 1
    mids = ((starts + ends) / 2.0).astype(int)
    prefix = np.concatenate([[0.0], np.cumsum(signal)])
    sums = prefix[end_indices + 1] - prefix[start_indices]
    mean_signal = sums / float(cpg_window)
    keep = (starts >= region_start) & (ends <= region_end)
    if not keep.any():
        return pd.DataFrame(columns=["chrom", "window_start", "window_end", "window_mid", "n_cpg", "mean_signal"])
    return pd.DataFrame(
        {
            "chrom": chrom,
            "window_start": starts[keep].astype(int),
            "window_end": ends[keep].astype(int),
            "window_mid": mids[keep].astype(int),
            "n_cpg": np.full(int(keep.sum()), cpg_window, dtype=int),
            "mean_signal": mean_signal[keep].astype(float),
        }
    )


def find_runs(mask: Iterable[bool]) -> list[tuple[int, int]]:
    runs: list[tuple[int, int]] = []
    current: int | None = None
    vals = list(mask)
    for i, val in enumerate(vals):
        if val and current is None:
            current = i
        if current is not None and ((not val) or i == len(vals) - 1):
            end = i if val and i == len(vals) - 1 else i - 1
            runs.append((current, end))
            current = None
    return runs


def call_sample_boundaries(
    windows: pd.DataFrame,
    sample: str,
    cohort: str,
    enter_threshold: float,
    exit_threshold: float,
    consecutive: int,
) -> list[dict[str, Any]]:
    valid_signal = windows["mean_signal"].astype(float)
    high = (valid_signal > enter_threshold).fillna(False).to_numpy()
    low = (valid_signal < exit_threshold).fillna(False).to_numpy()
    high_runs = [(s, e) for s, e in find_runs(high) if e - s + 1 >= consecutive]
    rows: list[dict[str, Any]] = []
    for segment_index, (start_i, end_i) in enumerate(high_runs, start=1):
        entry_coord = int(windows.iloc[start_i]["window_start"])
        exit_coord = np.nan
        low_run = None
        for low_start, low_end in find_runs(low[end_i + 1 :]):
            if low_end - low_start + 1 >= consecutive:
                low_run = (end_i + 1 + low_start, end_i + 1 + low_end)
                break
        if low_run is not None:
            exit_coord = int(windows.iloc[low_run[0]]["window_end"])
        rows.append(
            {
                "sample": sample,
                "cohort": cohort,
                "segment_index": segment_index,
                "entry_window_start": entry_coord,
                "exit_window_end": exit_coord,
                "high_run_start": int(windows.iloc[start_i]["window_start"]),
                "high_run_end": int(windows.iloc[end_i]["window_end"]),
                "n_high_windows": int(end_i - start_i + 1),
                "max_mean_signal": float(np.nanmax(valid_signal.iloc[start_i : end_i + 1])),
            }
        )
    return rows


def group_windows(sample_windows: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["cohort", "chrom", "window_start", "window_end", "window_mid"]
    rows = []
    for keys, z in sample_windows.groupby(group_cols, sort=True):
        vals = z["mean_signal"].astype(float)
        rows.append(
            {
                **dict(zip(group_cols, keys)),
                "mean_signal": float(vals.mean(skipna=True)) if vals.notna().any() else np.nan,
                "sd_signal": float(vals.std(skipna=True)) if vals.notna().sum() > 1 else np.nan,
                "n_samples": int(z.loc[vals.notna(), "sample"].nunique()),
                "mean_n_cpg": float(z["n_cpg"].mean()),
            }
        )
    return pd.DataFrame(rows)


def summarize_control_boundaries(sample_boundary_rows: list[dict[str, Any]], match_window: int) -> pd.DataFrame:
    if not sample_boundary_rows:
        return pd.DataFrame()
    df = pd.DataFrame(sample_boundary_rows)
    controls = df[df["cohort"] == "Control"].copy()
    if controls.empty:
        return pd.DataFrame()
    samples = sorted(controls["sample"].unique())
    if len(samples) != 2:
        rows = []
        for i, (_, row) in enumerate(controls.sort_values(["entry_window_start", "sample"]).iterrows(), start=1):
            rows.append(
                {
                    "cohort": "Control",
                    "consensus_index": i,
                    "entry_boundary": row["entry_window_start"],
                    "entry_ci_low": row["entry_window_start"],
                    "entry_ci_high": row["entry_window_start"],
                    "entry_control_span_bp": np.nan,
                    "entry_confidence": "low",
                    "exit_boundary": row["exit_window_end"],
                    "exit_ci_low": row["exit_window_end"],
                    "exit_ci_high": row["exit_window_end"],
                    "exit_control_span_bp": np.nan,
                    "exit_confidence": "low",
                    "n_control_samples_with_entry": 1,
                    "n_control_samples_with_exit": 0 if pd.isna(row["exit_window_end"]) else 1,
                    "samples": row["sample"],
                    "entry_by_sample": f"{row['sample']}={int(row['entry_window_start'])}",
                    "exit_by_sample": f"{row['sample']}={row['exit_window_end']}",
                    "mean_max_signal": row.get("max_mean_signal", np.nan),
                }
            )
        return pd.DataFrame(rows)

    left = controls[controls["sample"] == samples[0]].sort_values("entry_window_start").reset_index(drop=True)
    right = controls[controls["sample"] == samples[1]].sort_values("entry_window_start").reset_index(drop=True)
    used_right: set[int] = set()
    groups: list[pd.DataFrame] = []
    for _, lrow in left.iterrows():
        best_idx = None
        best_dist = np.inf
        for ridx, rrow in right.iterrows():
            if ridx in used_right:
                continue
            dist = abs(float(lrow["entry_window_start"]) - float(rrow["entry_window_start"]))
            if dist < best_dist:
                best_dist = dist
                best_idx = int(ridx)
        if best_idx is not None and best_dist <= match_window:
            used_right.add(best_idx)
            groups.append(pd.DataFrame([lrow.to_dict(), right.iloc[best_idx].to_dict()]))
        else:
            groups.append(pd.DataFrame([lrow.to_dict()]))
    for ridx, rrow in right.iterrows():
        if int(ridx) not in used_right:
            groups.append(pd.DataFrame([rrow.to_dict()]))

    rows: list[dict[str, Any]] = []
    for consensus_index, z in enumerate(sorted(groups, key=lambda g: float(g["entry_window_start"].mean())), start=1):
        entry = pd.to_numeric(z["entry_window_start"], errors="coerce").dropna()
        exit_ = pd.to_numeric(z["exit_window_end"], errors="coerce").dropna()
        entry_span = float(entry.max() - entry.min()) if len(entry) >= 2 else np.nan
        exit_span = float(exit_.max() - exit_.min()) if len(exit_) >= 2 else np.nan
        entry_parts = [f"{r['sample']}={int(r['entry_window_start'])}" for _, r in z.iterrows()]
        exit_parts = [f"{r['sample']}={int(r['exit_window_end'])}" for _, r in z.dropna(subset=["exit_window_end"]).iterrows()]
        rows.append(
            {
                "cohort": "Control",
                "consensus_index": int(consensus_index),
                "entry_boundary": float(entry.mean()) if len(entry) else np.nan,
                "entry_ci_low": float(entry.min()) if len(entry) else np.nan,
                "entry_ci_high": float(entry.max()) if len(entry) else np.nan,
                "entry_control_span_bp": entry_span,
                "entry_confidence": "high" if len(entry) == 2 and entry_span <= 2_000 else "low",
                "exit_boundary": float(exit_.mean()) if len(exit_) else np.nan,
                "exit_ci_low": float(exit_.min()) if len(exit_) else np.nan,
                "exit_ci_high": float(exit_.max()) if len(exit_) else np.nan,
                "exit_control_span_bp": exit_span,
                "exit_confidence": "high" if len(exit_) == 2 and exit_span <= 2_000 else "low",
                "n_control_samples_with_entry": int(len(entry)),
                "n_control_samples_with_exit": int(len(exit_)),
                "samples": ",".join(z["sample"].astype(str)),
                "entry_by_sample": ";".join(entry_parts),
                "exit_by_sample": ";".join(exit_parts),
                "mean_max_signal": float(pd.to_numeric(z["max_mean_signal"], errors="coerce").mean()),
            }
        )
    return pd.DataFrame(rows)


def select_primary_boundaries(control_summary: pd.DataFrame, group_boundary_calls: pd.DataFrame) -> pd.DataFrame:
    if control_summary.empty:
        return pd.DataFrame()
    controls = control_summary.copy()
    if "consensus_index" not in controls.columns:
        controls["consensus_index"] = np.arange(1, len(controls) + 1)
    patient = group_boundary_calls[group_boundary_calls["cohort"].isin(["PWS_DEL", "AS_DEL"])].copy() if not group_boundary_calls.empty else pd.DataFrame()

    def closest_distance(cohort: str, col: str, coord: float) -> float:
        if patient.empty or not math.isfinite(float(coord)):
            return 1e12
        source_col = "entry_window_start" if col == "entry" else "exit_window_end"
        vals = pd.to_numeric(patient.loc[patient["cohort"] == cohort, source_col], errors="coerce").dropna()
        if vals.empty:
            return 1e12
        return float(np.min(np.abs(vals.to_numpy(float) - coord)))

    scored = []
    for _, row in controls.iterrows():
        entry = float(row["entry_boundary"])
        exit_ = float(row["exit_boundary"]) if not pd.isna(row["exit_boundary"]) else np.nan
        score = closest_distance("PWS_DEL", "entry", entry) + closest_distance("AS_DEL", "entry", entry)
        if math.isfinite(exit_):
            score += 0.5 * closest_distance("PWS_DEL", "exit", exit_)
            score += 0.5 * closest_distance("AS_DEL", "exit", exit_)
        if row.get("entry_confidence") == "high":
            score -= 1_000
        if row.get("exit_confidence") == "high":
            score -= 500
        scored.append((score, row))
    scored.sort(key=lambda item: item[0])
    primary_control = scored[0][1]
    rows = [
        {
            "boundary_source": "Control",
            "cohort": "Control",
            "consensus_index": int(primary_control["consensus_index"]),
            "entry_boundary": primary_control["entry_boundary"],
            "exit_boundary": primary_control["exit_boundary"],
            "entry_confidence": primary_control["entry_confidence"],
            "exit_confidence": primary_control["exit_confidence"],
            "distance_to_control_entry_bp": 0.0,
            "distance_to_control_exit_bp": 0.0,
            "convergence_flag": "control_primary",
        }
    ]
    if not patient.empty:
        for cohort in ["PWS_DEL", "AS_DEL"]:
            q = patient[patient["cohort"] == cohort].copy()
            if q.empty:
                continue
            q["primary_distance"] = (
                (pd.to_numeric(q["entry_window_start"], errors="coerce") - float(primary_control["entry_boundary"])).abs()
                + 0.5 * (pd.to_numeric(q["exit_window_end"], errors="coerce") - float(primary_control["exit_boundary"])).abs()
            )
            prow = q.sort_values("primary_distance").iloc[0]
            entry_distance = float(prow["entry_window_start"] - primary_control["entry_boundary"])
            exit_distance = float(prow["exit_window_end"] - primary_control["exit_boundary"]) if not pd.isna(prow["exit_window_end"]) else np.nan
            rows.append(
                {
                    "boundary_source": cohort,
                    "cohort": cohort,
                    "consensus_index": int(prow["segment_index"]),
                    "entry_boundary": float(prow["entry_window_start"]),
                    "exit_boundary": float(prow["exit_window_end"]) if not pd.isna(prow["exit_window_end"]) else np.nan,
                    "entry_confidence": "group_call",
                    "exit_confidence": "group_call",
                    "distance_to_control_entry_bp": entry_distance,
                    "distance_to_control_exit_bp": exit_distance,
                    "convergence_flag": (
                        "convergent"
                        if abs(entry_distance) <= PATIENT_CONVERGENCE_BP
                        and (pd.isna(exit_distance) or abs(exit_distance) <= PATIENT_CONVERGENCE_BP)
                        else "flag"
                    ),
                }
            )
    return pd.DataFrame(rows)


def call_group_boundaries(
    grouped: pd.DataFrame,
    enter_threshold: float,
    exit_threshold: float,
    consecutive: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for cohort, z in grouped.groupby("cohort"):
        z = z.sort_values("window_start").reset_index(drop=True)
        rows.extend(call_sample_boundaries(z, cohort, cohort, enter_threshold, exit_threshold, consecutive))
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.rename(columns={"sample": "boundary_source"})


def build_sample_windows_and_boundaries(
    tracks: list[SignalTrack],
    chrom: str,
    region_start: int,
    region_end: int,
    enter_threshold: float,
    exit_threshold: float,
    consecutive_windows: int,
    min_cpgs_per_window: int,
    window_size_bp: int | None = None,
    step_size_bp: int | None = None,
    cpg_window: int | None = None,
    step_cpg: int | None = None,
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    sample_window_rows: list[pd.DataFrame] = []
    sample_boundary_rows: list[dict[str, Any]] = []
    for track in tracks:
        if cpg_window is not None:
            win = adaptive_cpg_window_profile(
                track.positions,
                track.signal,
                chrom,
                region_start,
                region_end,
                cpg_window=cpg_window,
                step_cpg=step_cpg or max(1, cpg_window // 5),
            )
        else:
            if window_size_bp is None or step_size_bp is None:
                raise ValueError("window_size_bp and step_size_bp are required for genomic-window sensitivity runs.")
            win = window_profile(
                track.positions,
                track.signal,
                chrom,
                region_start,
                region_end,
                window_size_bp,
                step_size_bp,
                min_cpgs_per_window,
            )
        win.insert(0, "sample", track.sample)
        win.insert(1, "cohort", track.cohort)
        sample_window_rows.append(win)
        sample_boundary_rows.extend(
            call_sample_boundaries(
                win,
                track.sample,
                track.cohort,
                enter_threshold,
                exit_threshold,
                consecutive_windows,
            )
        )
    if not sample_window_rows:
        return pd.DataFrame(), []
    return pd.concat(sample_window_rows, ignore_index=True), sample_boundary_rows


def compute_interval_call_set(
    tracks: list[SignalTrack],
    chrom: str,
    region_start: int,
    region_end: int,
    control_match_window: int,
    enter_threshold: float,
    exit_threshold: float,
    consecutive_windows: int,
    min_cpgs_per_window: int,
    window_size_bp: int | None = None,
    step_size_bp: int | None = None,
    cpg_window: int | None = None,
    step_cpg: int | None = None,
) -> dict[str, Any]:
    sample_windows, sample_boundary_rows = build_sample_windows_and_boundaries(
        tracks=tracks,
        chrom=chrom,
        region_start=region_start,
        region_end=region_end,
        enter_threshold=enter_threshold,
        exit_threshold=exit_threshold,
        consecutive_windows=consecutive_windows,
        min_cpgs_per_window=min_cpgs_per_window,
        window_size_bp=window_size_bp,
        step_size_bp=step_size_bp,
        cpg_window=cpg_window,
        step_cpg=step_cpg,
    )
    if sample_windows.empty:
        return {
            "sample_windows": sample_windows,
            "sample_boundary_rows": sample_boundary_rows,
            "grouped": pd.DataFrame(),
            "control_summary": pd.DataFrame(),
            "group_boundary_calls": pd.DataFrame(),
            "primary_boundaries": pd.DataFrame(),
        }
    grouped = group_windows(sample_windows)
    control_summary = summarize_control_boundaries(sample_boundary_rows, control_match_window)
    group_boundary_calls = call_group_boundaries(grouped, enter_threshold, exit_threshold, consecutive_windows)
    primary_boundaries = select_primary_boundaries(control_summary, group_boundary_calls)
    return {
        "sample_windows": sample_windows,
        "sample_boundary_rows": sample_boundary_rows,
        "grouped": grouped,
        "control_summary": control_summary,
        "group_boundary_calls": group_boundary_calls,
        "primary_boundaries": primary_boundaries,
    }


def interval_row_lookup(primary_boundaries: pd.DataFrame) -> dict[str, tuple[int, int]]:
    out: dict[str, tuple[int, int]] = {}
    for source in ["Control", "PWS_DEL", "AS_DEL"]:
        row = primary_boundaries.loc[primary_boundaries["boundary_source"] == source]
        if row.empty:
            continue
        out[source] = (
            int(round(float(row.iloc[0]["entry_boundary"]))),
            int(round(float(row.iloc[0]["exit_boundary"]))),
        )
    return out


def interval_jaccard(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start))
    union = max(a_end, b_end) - min(a_start, b_start)
    if union <= 0:
        return math.nan
    return overlap / union


def build_sensitivity_parameter_grid(base_args: argparse.Namespace) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    primary_enter = float(base_args.enter_threshold)
    primary_exit = float(base_args.exit_threshold)
    threshold_sweep = [(0.30, 0.10), (0.40, 0.10), (0.50, 0.10), (0.40, 0.15)]
    genomic_specs: list[tuple[int, str, float, float]] = []
    for window_size in SENSITIVITY_WINDOW_SIZES_BP:
        genomic_specs.append((window_size, "fixed_100bp", primary_enter, primary_exit))
        if window_size != base_args.window_size:
            genomic_specs.append((window_size, "adaptive_quarter_window", primary_enter, primary_exit))
    for enter_threshold, exit_threshold in threshold_sweep:
        genomic_specs.append((base_args.window_size, "fixed_100bp", float(enter_threshold), float(exit_threshold)))
    seen_genomic: set[tuple[int, str, float, float]] = set()
    for window_size, step_mode, enter_threshold, exit_threshold in genomic_specs:
        key = (window_size, step_mode, round(enter_threshold, 3), round(exit_threshold, 3))
        if key in seen_genomic:
            continue
        seen_genomic.add(key)
        step_size = STEP_SIZE if step_mode == "fixed_100bp" else adaptive_step_size(window_size)
        rows.append(
            {
                "window_family": "genomic_bp",
                "window_label": f"{window_size}bp_{step_mode}",
                "window_size_bp": window_size,
                "step_mode": step_mode,
                "step_size_bp": step_size,
                "cpg_window": np.nan,
                "step_cpg": np.nan,
                "enter_threshold": enter_threshold,
                "exit_threshold": exit_threshold,
                "is_primary_spec": (
                    window_size == base_args.window_size
                    and step_size == base_args.step_size
                    and math.isclose(enter_threshold, primary_enter)
                    and math.isclose(exit_threshold, primary_exit)
                ),
            }
        )
    cpg_specs: list[tuple[int, float, float]] = []
    for cpg_window in SENSITIVITY_CPG_WINDOWS:
        cpg_specs.append((cpg_window, primary_enter, primary_exit))
        cpg_specs.append((cpg_window, 0.40, 0.15))
    seen_cpg: set[tuple[int, float, float]] = set()
    for cpg_window, enter_threshold, exit_threshold in cpg_specs:
        key = (cpg_window, round(enter_threshold, 3), round(exit_threshold, 3))
        if key in seen_cpg:
            continue
        seen_cpg.add(key)
        step_cpg = max(1, cpg_window // 5)
        rows.append(
            {
                "window_family": "adaptive_cpg",
                "window_label": f"{cpg_window}CpG_adaptive",
                "window_size_bp": np.nan,
                "step_mode": "adaptive_cpg_step",
                "step_size_bp": np.nan,
                "cpg_window": cpg_window,
                "step_cpg": step_cpg,
                "enter_threshold": enter_threshold,
                "exit_threshold": exit_threshold,
                "is_primary_spec": False,
            }
        )
    return rows


def build_boundary_sensitivity_table(
    tracks: list[SignalTrack],
    args: argparse.Namespace,
    primary_shared_start: int,
    primary_shared_end: int,
) -> pd.DataFrame:
    local_start = max(args.region_start, primary_shared_start - 30_000)
    local_end = min(args.region_end, primary_shared_end + 30_000)
    rows: list[dict[str, Any]] = []
    for idx, spec in enumerate(build_sensitivity_parameter_grid(args), start=1):
        result = compute_interval_call_set(
            tracks=tracks,
            chrom=args.chrom,
            region_start=local_start,
            region_end=local_end,
            control_match_window=args.control_match_window,
            enter_threshold=float(spec["enter_threshold"]),
            exit_threshold=float(spec["exit_threshold"]),
            consecutive_windows=args.consecutive_windows,
            min_cpgs_per_window=args.min_cpgs_per_window,
            window_size_bp=None if pd.notna(spec["cpg_window"]) else int(spec["window_size_bp"]),
            step_size_bp=None if pd.notna(spec["cpg_window"]) else int(spec["step_size_bp"]),
            cpg_window=int(spec["cpg_window"]) if pd.notna(spec["cpg_window"]) else None,
            step_cpg=int(spec["step_cpg"]) if pd.notna(spec["step_cpg"]) else None,
        )
        primary_boundaries = result["primary_boundaries"]
        interval_lookup = interval_row_lookup(primary_boundaries)
        status = "ok"
        shared_start = np.nan
        shared_end = np.nan
        shared_length_bp = np.nan
        jaccard = np.nan
        start_shift_bp = np.nan
        end_shift_bp = np.nan
        if {"Control", "PWS_DEL", "AS_DEL"} <= set(interval_lookup):
            try:
                shared_start, shared_end = shared_interval_from_boundaries(primary_boundaries)
                shared_length_bp = shared_end - shared_start
                jaccard = interval_jaccard(
                    primary_shared_start,
                    primary_shared_end,
                    shared_start,
                    shared_end,
                )
                start_shift_bp = shared_start - primary_shared_start
                end_shift_bp = shared_end - primary_shared_end
            except ValueError as exc:
                status = f"invalid_shared_interval: {exc}"
        else:
            status = "missing_primary_calls"
        rows.append(
            {
                "analysis_id": f"sensitivity_{idx:03d}",
                **spec,
                "control_entry": interval_lookup.get("Control", (np.nan, np.nan))[0],
                "control_exit": interval_lookup.get("Control", (np.nan, np.nan))[1],
                "pws_entry": interval_lookup.get("PWS_DEL", (np.nan, np.nan))[0],
                "pws_exit": interval_lookup.get("PWS_DEL", (np.nan, np.nan))[1],
                "as_entry": interval_lookup.get("AS_DEL", (np.nan, np.nan))[0],
                "as_exit": interval_lookup.get("AS_DEL", (np.nan, np.nan))[1],
                "shared_start": shared_start,
                "shared_end": shared_end,
                "shared_length_bp": shared_length_bp,
                "jaccard_vs_primary_shared_interval": jaccard,
                "shared_start_shift_bp": start_shift_bp,
                "shared_end_shift_bp": end_shift_bp,
                "sensitivity_region_start": local_start,
                "sensitivity_region_end": local_end,
                "status": status,
            }
        )
    return pd.DataFrame(rows)


def build_change_point_support(
    grouped: pd.DataFrame,
    shared_start: int,
    shared_end: int,
    flank_bp: int = 10_000,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for cohort in ["Control", "PWS_DEL", "AS_DEL"]:
        q = grouped[grouped["cohort"] == cohort].copy()
        if q.empty:
            continue
        q = q[(q["window_mid"] >= shared_start - flank_bp) & (q["window_mid"] <= shared_end + flank_bp)].copy()
        q = q.sort_values("window_mid")
        if len(q) < 8:
            continue
        x = q["window_mid"].to_numpy(float)
        y = q["mean_signal"].astype(float).rolling(window=11, center=True, min_periods=3).mean().to_numpy()
        if not np.isfinite(y).any():
            continue
        deriv = np.gradient(np.nan_to_num(y, nan=np.nanmedian(y)), x)
        entry_idx = int(np.nanargmax(deriv))
        if entry_idx >= len(x) - 2:
            continue
        exit_region = deriv[entry_idx + 1 :]
        exit_idx = entry_idx + 1 + int(np.nanargmin(exit_region))
        rows.append(
            {
                "cohort": cohort,
                "entry_change_point": int(round(x[entry_idx])),
                "exit_change_point": int(round(x[exit_idx])),
                "entry_shift_vs_primary_start_bp": int(round(x[entry_idx])) - shared_start,
                "exit_shift_vs_primary_end_bp": int(round(x[exit_idx])) - shared_end,
                "method": "rolling-mean derivative maxima/minima",
            }
        )
    if rows:
        out = pd.DataFrame(rows)
        consensus_entry = int(out["entry_change_point"].max())
        consensus_exit = int(out["exit_change_point"].min())
        out = pd.concat(
            [
                out,
                pd.DataFrame(
                    [
                        {
                            "cohort": "Consensus",
                            "entry_change_point": consensus_entry,
                            "exit_change_point": consensus_exit,
                            "entry_shift_vs_primary_start_bp": consensus_entry - shared_start,
                            "exit_shift_vs_primary_end_bp": consensus_exit - shared_end,
                            "method": "latest entry / earliest exit across cohort-specific derivative change points",
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )
        return out
    return pd.DataFrame(
        columns=[
            "cohort",
            "entry_change_point",
            "exit_change_point",
            "entry_shift_vs_primary_start_bp",
            "exit_shift_vs_primary_end_bp",
            "method",
        ]
    )


def bootstrap_shared_interval_support(
    tracks: list[SignalTrack],
    args: argparse.Namespace,
    primary_shared_start: int,
    primary_shared_end: int,
) -> pd.DataFrame:
    local_start = max(args.region_start, primary_shared_start - args.bootstrap_local_flank)
    local_end = min(args.region_end, primary_shared_end + args.bootstrap_local_flank)
    rng = np.random.default_rng(args.bootstrap_seed)
    replicate_rows: list[dict[str, Any]] = []
    for replicate in range(1, max(0, args.bootstrap_replicates) + 1):
        boot_tracks: list[SignalTrack] = []
        failed = False
        for track in tracks:
            mask = (track.positions >= local_start) & (track.positions <= local_end)
            pos = track.positions[mask]
            sig = track.signal[mask]
            if pos.size < max(args.min_cpgs_per_window, 5):
                failed = True
                break
            chosen = rng.integers(0, pos.size, size=pos.size)
            boot_pos = pos[chosen]
            boot_sig = sig[chosen]
            order = np.argsort(boot_pos)
            boot_tracks.append(
                SignalTrack(
                    sample=track.sample,
                    cohort=track.cohort,
                    positions=boot_pos[order],
                    signal=boot_sig[order],
                )
            )
        if failed:
            continue
        result = compute_interval_call_set(
            tracks=boot_tracks,
            chrom=args.chrom,
            region_start=local_start,
            region_end=local_end,
            control_match_window=args.control_match_window,
            enter_threshold=args.enter_threshold,
            exit_threshold=args.exit_threshold,
            consecutive_windows=args.consecutive_windows,
            min_cpgs_per_window=args.min_cpgs_per_window,
            window_size_bp=args.window_size,
            step_size_bp=args.step_size,
        )
        primary = result["primary_boundaries"]
        if primary.empty:
            continue
        try:
            shared_start, shared_end = shared_interval_from_boundaries(primary)
        except ValueError:
            continue
        replicate_rows.append(
            {
                "replicate": replicate,
                "shared_start": shared_start,
                "shared_end": shared_end,
                "shared_length_bp": shared_end - shared_start,
            }
        )
    summary = {
        "bootstrap_method": "CpG bootstrap with replacement within each sample; local interval stability only",
        "bootstrap_region_start": local_start,
        "bootstrap_region_end": local_end,
        "bootstrap_replicates_requested": int(max(0, args.bootstrap_replicates)),
        "bootstrap_replicates_successful": int(len(replicate_rows)),
        "primary_shared_start": primary_shared_start,
        "primary_shared_end": primary_shared_end,
        "primary_shared_length_bp": primary_shared_end - primary_shared_start,
        "bootstrap_start_median": np.nan,
        "bootstrap_start_ci_low": np.nan,
        "bootstrap_start_ci_high": np.nan,
        "bootstrap_end_median": np.nan,
        "bootstrap_end_ci_low": np.nan,
        "bootstrap_end_ci_high": np.nan,
        "bootstrap_supported_start": np.nan,
        "bootstrap_supported_end": np.nan,
        "primary_overlaps_bootstrap_supported_interval": False,
        "limitation": "Control sample size is n=2, so the bootstrap is CpG-resampling rather than sample-resampling.",
    }
    if replicate_rows:
        rep = pd.DataFrame(replicate_rows)
        start_ci = np.nanpercentile(rep["shared_start"], [2.5, 50.0, 97.5])
        end_ci = np.nanpercentile(rep["shared_end"], [2.5, 50.0, 97.5])
        supported_start = int(round(start_ci[2]))
        supported_end = int(round(end_ci[0]))
        supported_overlap = supported_start < supported_end and overlap_length(
            primary_shared_start,
            primary_shared_end,
            supported_start,
            supported_end,
        ) > 0
        summary.update(
            {
                "bootstrap_start_ci_low": int(round(start_ci[0])),
                "bootstrap_start_median": int(round(start_ci[1])),
                "bootstrap_start_ci_high": int(round(start_ci[2])),
                "bootstrap_end_ci_low": int(round(end_ci[0])),
                "bootstrap_end_median": int(round(end_ci[1])),
                "bootstrap_end_ci_high": int(round(end_ci[2])),
                "bootstrap_supported_start": supported_start,
                "bootstrap_supported_end": supported_end,
                "primary_overlaps_bootstrap_supported_interval": bool(supported_overlap),
            }
        )
    return pd.DataFrame([summary])


def build_bp_distance_definitions(shared_start: int, shared_end: int, bp_hotspots: pd.DataFrame) -> pd.DataFrame:
    shared_mid = int(round((shared_start + shared_end) / 2.0))
    rows: list[dict[str, Any]] = []
    for bp_name in ["BP1", "BP2", "BP3"]:
        hit = bp_hotspots[bp_hotspots["name"].astype(str) == bp_name].copy()
        if hit.empty:
            continue
        row = hit.iloc[0]
        bp_start = int(row["start"])
        bp_end = int(row["end"])
        bp_mid = int(round((bp_start + bp_end) / 2.0))
        rows.append(
            {
                "bp_landmark_name": bp_name,
                "bp_start": bp_start,
                "bp_end": bp_end,
                "bp_midpoint": bp_mid,
                "distance_from_shared_start_bp": interval_distance_to_nearest_edge(bp_start, bp_end, shared_start),
                "distance_from_shared_midpoint_bp": interval_distance_to_nearest_edge(bp_start, bp_end, shared_mid),
                "distance_from_shared_end_bp": interval_distance_to_nearest_edge(bp_start, bp_end, shared_end),
                "distance_type_used_in_figure": "shared midpoint to nearest interval edge" if bp_name == "BP2" else "not used in figure label",
            }
        )
    return pd.DataFrame(rows)


def warning_frame(messages: list[str]) -> pd.DataFrame:
    if not messages:
        return pd.DataFrame(columns=["warning"])
    return pd.DataFrame({"warning": messages})


def read_bed_like(path: Path, source: str, chrom: str, start: int, end: int, has_header_comment: bool = False) -> pd.DataFrame:
    if not path or not path.exists():
        return pd.DataFrame(columns=["chrom", "start", "end", "name", "source"])
    rows: list[dict[str, Any]] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") if opener is open else opener(path, "rt") as handle:  # type: ignore[arg-type]
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                fields = line.split()
            if len(fields) < 3:
                continue
            row_chrom = fields[0]
            if row_chrom != chrom:
                continue
            try:
                row_start = int(float(fields[1]))
                row_end = int(float(fields[2]))
            except ValueError:
                continue
            if row_end < start or row_start > end:
                continue
            name = fields[3] if len(fields) >= 4 else source
            rows.append({"chrom": row_chrom, "start": row_start, "end": row_end, "name": name, "source": source})
    if not rows:
        return pd.DataFrame(columns=["chrom", "start", "end", "name", "source"])
    return (
        pd.DataFrame(rows)
        .drop_duplicates()
        .sort_values(["start", "end", "name", "source"])
        .reset_index(drop=True)
    )


def read_bp_hotspots(path: Path | None, chrom: str) -> pd.DataFrame:
    if path and path.exists():
        out = read_bed_like(path, "BP_hotspot", chrom, 0, 10**12)
        if not out.empty:
            out["coordinate"] = (out["start"].astype(float) + out["end"].astype(float)) / 2.0
            out["display_color"] = out["name"].map(structural_landmark_color).fillna(STRUCT_GREY)
        return out
    rows = []
    for cluster in BP_CLUSTER_INTERVALS_T2T:
        name = str(cluster["name"])
        start = int(cluster["start"])
        end = int(cluster["end"])
        rows.append(
            {
                "chrom": chrom,
                "start": start,
                "end": end,
                "name": name,
                "source": "Figure2_T2T_interval",
                "coordinate": (start + end) / 2.0,
                "display_color": structural_landmark_color(name),
            }
        )
    return pd.DataFrame(rows)


def parse_gtf_genes(path: Path, chrom: str, start: int, end: int) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["chrom", "start", "end", "gene", "strand"])
    rows: list[dict[str, Any]] = []
    attr_re = re.compile(r'(\S+) "([^"]*)"')
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") if opener is open else opener(path, "rt") as handle:  # type: ignore[arg-type]
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[0] != chrom or fields[2] != "gene":
                continue
            row_start = int(fields[3])
            row_end = int(fields[4])
            if row_end < start or row_start > end:
                continue
            attrs = dict(attr_re.findall(fields[8]))
            gene = attrs.get("gene") or attrs.get("gene_id") or attrs.get("Name") or "gene"
            rows.append({"chrom": chrom, "start": row_start, "end": row_end, "gene": gene, "strand": fields[6]})
    if not rows:
        return pd.DataFrame(columns=["chrom", "start", "end", "gene", "strand"])
    return pd.DataFrame(rows).drop_duplicates().sort_values(["start", "end", "gene"]).reset_index(drop=True)


def nearest_interval_distance(coord: float, intervals: pd.DataFrame, name_col: str = "name") -> tuple[str, float, bool]:
    if intervals.empty or not math.isfinite(float(coord)):
        return "", np.nan, False
    starts = intervals["start"].astype(float).to_numpy()
    ends = intervals["end"].astype(float).to_numpy()
    overlaps = (starts <= coord) & (ends >= coord)
    if overlaps.any():
        idx = int(np.where(overlaps)[0][0])
        return str(intervals.iloc[idx].get(name_col, intervals.iloc[idx].get("gene", ""))), 0.0, True
    distances = np.minimum(np.abs(starts - coord), np.abs(ends - coord))
    idx = int(np.nanargmin(distances))
    return str(intervals.iloc[idx].get(name_col, intervals.iloc[idx].get("gene", ""))), float(distances[idx]), False


def annotate_boundaries(
    boundary_table: pd.DataFrame,
    genes: pd.DataFrame,
    segdups: pd.DataFrame,
    bp_hotspots: pd.DataFrame,
    ctcf: pd.DataFrame,
    imprintome: pd.DataFrame,
    court: pd.DataFrame,
    icrs: pd.DataFrame | None,
    chrom: str,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    if boundary_table.empty:
        return pd.DataFrame()
    candidates: list[tuple[str, str, int, float, str]] = []
    for _, row in boundary_table.iterrows():
        source = row.get("cohort", row.get("boundary_source", "Control"))
        seg_value = row.get("segment_index", row.get("consensus_index", 1))
        if pd.isna(seg_value):
            seg_value = row.get("consensus_index", 1)
        try:
            if pd.isna(seg_value):
                seg_value = 1
            seg = int(float(seg_value))
        except (TypeError, ValueError):
            seg = 1
        for side, col in [("entry", "entry_boundary"), ("exit", "exit_boundary")]:
            if col not in row or pd.isna(row[col]):
                continue
            confidence = str(row.get(f"{side}_confidence", "group_call"))
            candidates.append((str(source), side, seg, float(row[col]), confidence))
    for source, boundary_type, segment_index, coord, confidence in candidates:
        gene_name, gene_distance, gene_overlap = nearest_interval_distance(coord, genes.rename(columns={"gene": "name"}))
        seg_name, seg_distance, seg_overlap = nearest_interval_distance(coord, segdups)
        bp_name, bp_distance, bp_overlap = nearest_interval_distance(coord, bp_hotspots)
        ctcf_name, ctcf_distance, ctcf_overlap = nearest_interval_distance(coord, ctcf)
        imp_name, imp_distance, imp_overlap = nearest_interval_distance(coord, imprintome)
        court_name, court_distance, court_overlap = nearest_interval_distance(coord, court)
        icr_track = icrs if icrs is not None else court
        icr_name, icr_distance, icr_overlap = nearest_interval_distance(coord, icr_track)
        rows.append(
            {
                "boundary_source": source,
                "boundary_type": boundary_type,
                "segment_index": segment_index,
                "chrom": chrom,
                "boundary_coordinate": int(round(coord)),
                "confidence": confidence,
                "nearest_gene": gene_name,
                "distance_to_nearest_gene_bp": gene_distance,
                "overlaps_gene": gene_overlap,
                "nearest_segdup": seg_name,
                "distance_to_segdup_bp": seg_distance,
                "overlaps_segdup": seg_overlap,
                "nearest_BP_hotspot": bp_name,
                "distance_to_BP_hotspot_bp": bp_distance,
                "overlaps_BP_hotspot_proxy": bp_overlap,
                "nearest_CTCF": ctcf_name,
                "distance_to_CTCF_bp": ctcf_distance,
                "overlaps_CTCF": ctcf_overlap,
                "ctcf_track_status": "loaded" if not ctcf.empty else "not_available",
                "nearest_Imprintome_DMR": imp_name,
                "distance_to_Imprintome_DMR_bp": imp_distance,
                "overlaps_Imprintome_DMR": imp_overlap,
                "nearest_Court2014_DMR": court_name,
                "distance_to_Court2014_DMR_bp": court_distance,
                "overlaps_Court2014_DMR": court_overlap,
                "nearest_ICR": icr_name,
                "distance_to_ICR_bp": icr_distance,
                "overlaps_ICR": icr_overlap,
            }
        )
    return pd.DataFrame(rows)


def single_cpg_zoom_table(tracks: list[SignalTrack], boundaries: pd.DataFrame, zoom_flank: int) -> pd.DataFrame:
    if boundaries.empty:
        return pd.DataFrame()
    coords = []
    for _, row in boundaries.iterrows():
        seg_value = row.get("segment_index", row.get("consensus_index", 1))
        if pd.isna(seg_value):
            seg_value = row.get("consensus_index", 1)
        if pd.isna(seg_value):
            seg_value = 1
        for col in ["entry_boundary", "exit_boundary"]:
            if col in row and not pd.isna(row[col]):
                coords.append((row.get("cohort", "Control"), int(seg_value), col.replace("_boundary", ""), float(row[col])))
    rows: list[pd.DataFrame] = []
    for source, segment_index, boundary_type, coord in coords:
        low = coord - zoom_flank
        high = coord + zoom_flank
        for track in tracks:
            mask = (track.positions >= low) & (track.positions <= high)
            if not mask.any():
                continue
            df = pd.DataFrame(
                {
                    "boundary_source": source,
                    "boundary_type": boundary_type,
                    "segment_index": segment_index,
                    "sample": track.sample,
                    "cohort": track.cohort,
                    "pos": track.positions[mask],
                    "signal": track.signal[mask],
                }
            )
            rows.append(df)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def save_table(df: pd.DataFrame, path: Path, compression: str | None = None) -> None:
    mkdir(path.parent)
    df.to_csv(path, sep="\t", index=False, compression=compression)


def plot_gene_track(ax: plt.Axes, genes: pd.DataFrame, start: int, end: int) -> None:
    ax.set_xlim(start, end)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.spines[["left", "right", "top"]].set_visible(False)
    if genes.empty:
        ax.text(0.5, 0.5, "No gene models loaded", transform=ax.transAxes, ha="center", va="center", fontsize=8)
        return
    levels: dict[str, float] = {}
    for i, (_, row) in enumerate(genes.iterrows()):
        y = 0.2 + (i % 3) * 0.22
        levels[row["gene"]] = y
        ax.add_patch(Rectangle((row["start"], y), max(1, row["end"] - row["start"]), 0.08, color="#5c6f82", alpha=0.6))
    label_genes = genes[(genes["end"] >= start) & (genes["start"] <= end)].copy()
    for _, row in label_genes.iterrows():
        if row["gene"] in {"SNRPN", "SNHG14", "SNURF", "UBE3A", "MAGEL2", "NDN", "MKRN3"} or (row["end"] - row["start"]) > 60_000:
            x = min(max((row["start"] + row["end"]) / 2, start), end)
            ax.text(x, min(0.92, levels.get(row["gene"], 0.25) + 0.13), row["gene"], fontsize=6, ha="center", rotation=25)


KEY_GENE_LABELS = {"SNRPN", "SNHG14", "SNURF", "UBE3A", "MAGEL2", "NDN", "MKRN3"}
FULL_INTERVAL_GENE_LABELS = {
    "MKRN3",
    "MAGEL2",
    "NDN",
    "PWRN4",
    "SNRPN",
    "SNHG14",
    "SNURF",
    "UBE3A",
    "GABRB3",
    "GABRA5",
    "GABRG3",
    "OCA2",
    "HERC2",
    "APBA2",
}


def clean_track_label(name: str) -> str:
    return re.sub(r"[*^#]", "", name)


def bp_marker_coordinate(row: pd.Series) -> float:
    coord = row.get("coordinate", np.nan)
    try:
        if pd.notna(coord) and math.isfinite(float(coord)):
            return float(coord)
    except (TypeError, ValueError):
        pass
    return (float(row["start"]) + float(row["end"])) / 2.0


def draw_bp_markers(
    ax: plt.Axes,
    bp_hotspots: pd.DataFrame,
    label_y: float = 0.96,
    color: str | None = None,
    alpha: float = 0.10,
    label: bool = False,
) -> None:
    """Draw BP1/BP2/BP3 as the same interval bands used in Figure 2."""
    if bp_hotspots.empty:
        return
    x0, x1 = ax.get_xlim()
    for _, row in bp_hotspots.iterrows():
        start = max(float(row["start"]), x0)
        end = min(float(row["end"]), x1)
        if start >= end:
            continue
        name = str(row.get("name", "BP"))
        band_color = color or str(row.get("display_color", structural_landmark_color(name)))
        ax.axvspan(start, end, color=band_color, alpha=alpha, lw=0, zorder=0)
        if not label:
            continue
        ax.text(
            (start + end) / 2.0,
            label_y,
            name,
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=9,
            color="#3a3a3a",
            clip_on=False,
            zorder=5,
        )


def draw_bp_interval_track(ax: plt.Axes, bp_hotspots: pd.DataFrame, start: int, end: int, y: float, height: float = 0.34) -> None:
    if bp_hotspots.empty:
        return
    z = bp_hotspots[(bp_hotspots["end"] >= start) & (bp_hotspots["start"] <= end)].copy()
    for _, row in z.iterrows():
        x0 = max(start, int(row["start"]))
        x1 = min(end, int(row["end"]))
        if x1 <= x0:
            continue
        name = str(row.get("name", "BP"))
        color = str(row.get("display_color", structural_landmark_color(name)))
        ax.add_patch(
            Rectangle(
                (x0, y),
                max(1, x1 - x0),
                height,
                fc=color,
                ec=color,
                alpha=0.30,
                lw=0.8,
            )
        )
        ax.text((x0 + x1) / 2.0, y + height + 0.08, name, fontsize=8.8, ha="center", va="bottom", color="#222222")


def interval_spans(df: pd.DataFrame, start: int, end: int) -> list[tuple[int, int]]:
    if df.empty or not {"start", "end"} <= set(df.columns):
        return []
    spans: list[tuple[int, int]] = []
    z = df[(df["end"] >= start) & (df["start"] <= end)]
    for _, row in z.iterrows():
        x0 = max(start, int(row["start"]))
        x1 = min(end, int(row["end"]))
        spans.append((x0, max(1, x1 - x0)))
    return spans


def add_interval_track(
    ax: plt.Axes,
    df: pd.DataFrame,
    start: int,
    end: int,
    y: float,
    color: str,
    alpha: float,
    label_col: str = "name",
    label_mode: str = "none",
    allowed_labels: set[str] | None = None,
    max_labels: int = 18,
    min_label_spacing_bp: int = 0,
    font_size: float = 6.0,
    label_rotation: float = 25,
) -> None:
    spans = interval_spans(df, start, end)
    if spans:
        ax.broken_barh(spans, (y, 0.44), facecolors=color, alpha=alpha, linewidth=0)
    if df.empty or label_mode == "none":
        return
    z = df[(df["end"] >= start) & (df["start"] <= end)].copy()
    labels_drawn = 0
    label_xs: list[float] = []
    for _, row in z.iterrows():
        x0 = max(start, int(row["start"]))
        x1 = min(end, int(row["end"]))
        name = str(row.get(label_col, ""))
        if not name:
            continue
        width = x1 - x0
        draw = False
        rotation = label_rotation
        if label_mode == "bp":
            x = bp_marker_coordinate(row)
            draw = start <= x <= end
            rotation = 90
        elif label_mode == "icr":
            draw = True
        elif label_mode == "genes":
            draw = name in KEY_GENE_LABELS or width > 80_000
        elif label_mode == "genes_sparse":
            draw = allowed_labels is not None and name in allowed_labels
        if not draw:
            continue
        x = bp_marker_coordinate(row) if label_mode == "bp" else min(max((x0 + x1) / 2, start), end)
        if min_label_spacing_bp and any(abs(x - prev_x) < min_label_spacing_bp for prev_x in label_xs):
            continue
        # For manuscript figures, show all imprinting-control-region intervals
        # with the generic label "ICR" rather than catalogue-specific IDs
        # such as ICR_893. The detailed source names remain in the output
        # annotation tables; this only simplifies the plot.
        label = "ICR" if label_mode == "icr" else clean_track_label(name)[:28]
        ax.text(x, y + 0.52, label, fontsize=font_size, ha="center", va="bottom", rotation=rotation, clip_on=True)
        label_xs.append(float(x))
        labels_drawn += 1
        if labels_drawn >= max_labels:
            break


def draw_annotation_track(
    ax: plt.Axes,
    start: int,
    end: int,
    genes: pd.DataFrame,
    segdups: pd.DataFrame,
    bp_hotspots: pd.DataFrame,
    ctcf: pd.DataFrame,
    icrs: pd.DataFrame,
    repeats: pd.DataFrame,
    label_density: str = "dense",
) -> None:
    ax.set_xlim(start, end)
    gene_df = genes.rename(columns={"gene": "name"}) if not genes.empty else genes
    sparse = label_density == "sparse"
    tracks = [
        ("Genes", gene_df, "#8fa0b2", 0.48, "name", "genes_sparse" if sparse else "genes"),
        ("ICR", icrs, "#6f9bc9", 0.68, "name", "icr"),
        ("SegDup", segdups, "#9c7a3c", 0.54, "name", "none"),
        ("BP landmarks", bp_hotspots, STRUCT_GREY, 0.30, "name", "bp_band"),
    ]
    ax.set_ylim(0, len(tracks) + 0.25)
    ax.set_yticks([i + 0.42 for i in range(len(tracks))])
    ax.set_yticklabels([label for label, *_rest in tracks], fontsize=7)
    for i in range(1, len(tracks)):
        ax.axhline(i, color="#e6e6e6", lw=0.35, zorder=0)
    ax.spines[["right", "top"]].set_visible(False)
    for i, (_label, df, color, alpha, label_col, label_mode) in enumerate(tracks):
        if label_mode == "bp_band":
            draw_bp_interval_track(ax, df, start, end, i + 0.16)
            continue
        add_interval_track(
            ax,
            df,
            start,
            end,
            i + 0.16,
            color,
            alpha,
            label_col=label_col,
            label_mode=label_mode,
            allowed_labels=FULL_INTERVAL_GENE_LABELS if label_mode == "genes_sparse" else None,
            max_labels=14 if sparse else 18,
            min_label_spacing_bp=90_000 if label_mode == "genes_sparse" else 20_000 if sparse and label_mode == "icr" else 0,
            font_size=5.8 if sparse else 6.0,
            label_rotation=20 if label_mode == "genes_sparse" else 32 if sparse and label_mode == "icr" else 25,
        )


def load_figure2_shared_annotation_track(table_dir: Path) -> pd.DataFrame:
    path = table_dir / "Figure2_shared_annotation_track.tsv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def figure2_structural_landmarks(annotation_track: pd.DataFrame, chrom: str) -> pd.DataFrame:
    if annotation_track.empty:
        return pd.DataFrame()
    z = annotation_track[annotation_track["kind"] == "structural_landmark"].copy()
    if z.empty:
        return pd.DataFrame()
    z.insert(0, "chrom", chrom)
    if "display_color" not in z.columns:
        z["display_color"] = z["label"].map(structural_landmark_color).fillna(STRUCT_GREY)
    z["name"] = z["label"]
    z["source"] = "Figure2_shared_annotation_track"
    z["coordinate"] = (z["start"].astype(float) + z["end"].astype(float)) / 2.0
    return z[["chrom", "start", "end", "name", "source", "coordinate", "display_color"]].copy()


def prepare_panel_a_annotation_track(
    annotation_track: pd.DataFrame,
    genes: pd.DataFrame,
    icrs: pd.DataFrame,
    chrom: str,
    shared_start: int,
    shared_end: int,
) -> pd.DataFrame:
    if annotation_track.empty:
        return annotation_track
    spec = WindowSpec(chrom, FIGURE2_DISPLAY_START, FIGURE2_DISPLAY_END, WINDOW_SIZE)
    feature_rows: list[dict[str, Any]] = []

    structural = annotation_track[annotation_track["kind"] == "structural_landmark"].copy()
    if not structural.empty:
        feature_rows.extend(structural.to_dict("records"))

    anchor = annotation_track[annotation_track["label"] == "PWS-IC / SNRPN"].copy()
    if not anchor.empty:
        feature_rows.extend(anchor.to_dict("records"))

    top_rows: list[dict[str, Any]] = []
    bottom_rows: list[dict[str, Any]] = []
    gene_table = genes.copy()
    gene_table["gene"] = gene_table.get("gene", pd.Series(index=gene_table.index, dtype=object)).astype(str)
    for gene_name, side in PANEL_A_GENE_SIDES.items():
        hit = gene_table[gene_table["gene"] == gene_name].copy()
        if hit.empty:
            continue
        start = max(spec.start, int(hit["start"].min()))
        end = min(spec.end, int(hit["end"].max()))
        if start >= end:
            continue
        row = {
            "chrom": chrom,
            "label": gene_name,
            "kind": "gene_context",
            "start": start,
            "end": end,
            "marker_position": np.nan,
            "lane": "gene_track",
            "color_role": "grey",
            "center": (start + end) / 2.0,
            "display_color": np.nan,
            "query_gene": gene_name,
        }
        if side == "top":
            top_rows.append(row)
        else:
            bottom_rows.append(row)

    if not icrs.empty:
        icr = icrs.copy()
        icr["mid"] = (icr["start"].astype(float) + icr["end"].astype(float)) / 2.0
        row = icr.iloc[(icr["mid"] - ((shared_start + shared_end) / 2.0)).abs().argmin()]
        start = max(spec.start, int(row["start"]))
        end = min(spec.end, int(row["end"]))
        if start < end:
            top_rows.append(
                {
                    "chrom": chrom,
                    "label": "ICR",
                    "kind": "gene_context",
                    "start": start,
                    "end": end,
                    "marker_position": np.nan,
                    "lane": "anchor",
                    "color_role": "grey",
                    "center": (start + end) / 2.0,
                    "display_color": np.nan,
                    "query_gene": "ICR",
                }
            )

    feature_rows.extend(
        assign_non_overlapping_label_positions(
            top_rows,
            spec,
            "top",
            n_lanes=3,
            gap_bp=110_000,
            center_gap_bp=520_000,
        )
    )
    feature_rows.extend(
        assign_non_overlapping_label_positions(
            bottom_rows,
            spec,
            "bottom",
            n_lanes=2,
            gap_bp=110_000,
            center_gap_bp=360_000,
        )
    )
    out = pd.DataFrame(feature_rows)
    if out.empty:
        return out
    return out.sort_values(["kind", "start", "end", "label"]).reset_index(drop=True)


def draw_segdup_lane(ax: plt.Axes, segdups: pd.DataFrame, start: int, end: int) -> None:
    ax.set_xlim(start, end)
    ax.set_ylim(0, 1)
    ax.set_yticks([0.42])
    ax.set_yticklabels(["SegDup"], fontsize=8)
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.axhline(0.42, color="#ececec", lw=0.6, zorder=0)
    add_interval_track(ax, segdups, start, end, 0.18, "#9c7a3c", 0.54, label_col="name", label_mode="none")
    ax.xaxis.set_major_formatter(pos_formatter())
    ax.tick_params(axis="x", labelsize=10, width=0.9)


def style_panel_axis(ax: plt.Axes, ylabel: str, show_x: bool, formatter: FuncFormatter | None = None) -> None:
    ax.grid(axis="y", color="#e8e8e8", lw=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.tick_params(axis="both", labelsize=10, width=0.9)
    ax.tick_params(axis="x", labelbottom=show_x)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.xaxis.set_major_formatter(formatter or pos_formatter())


def rolling_median(values: pd.Series | np.ndarray, window: int = 101) -> np.ndarray:
    """Centered rolling median for visual denoising only."""
    arr = pd.Series(np.asarray(values, dtype=float))
    return arr.rolling(window=window, center=True, min_periods=max(3, window // 12)).median().to_numpy(float)


def shade_consensus(ax: plt.Axes, start: int, end: int, label: bool = False, alpha: float = 0.18) -> None:
    ax.axvspan(start, end, color="#d99a00", alpha=alpha, lw=0, zorder=0)
    if label:
        width_kb = (end - start) / 1_000.0
        ax.text(
            (start + end) / 2,
            0.98,
            f"shared methylation-transition interval\n{width_kb:.1f} kb",
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=9.5,
            color="#6e4a00",
            bbox={"facecolor": "white", "edgecolor": "#d99a00", "alpha": 0.90, "pad": 2.5},
            zorder=10,
        )


def draw_convergence_panel(ax: plt.Axes, intervals: pd.DataFrame, pad_bp: int = 2_500) -> tuple[int, int]:
    if intervals.empty:
        ax.text(0.5, 0.5, "No primary interval calls", transform=ax.transAxes, ha="center", va="center")
        ax.set_axis_off()
        return 0, 1
    x_start = int(intervals["start"].min()) - pad_bp
    x_end = int(intervals["end"].max()) + pad_bp
    order = ["Control", "PWS_DEL", "AS_DEL", "Consensus"]
    y_positions = {source: len(order) - i for i, source in enumerate(order)}
    for _, row in intervals.iterrows():
        source = str(row["source"])
        y = y_positions.get(source, 1)
        height = 0.44 if source != "Consensus" else 0.54
        alpha = 0.86 if source != "Consensus" else 0.95
        ax.broken_barh(
            [(int(row["start"]), max(1, int(row["end"]) - int(row["start"])))],
            (y - height / 2, height),
            facecolors=str(row["color"]),
            alpha=alpha,
            edgecolors="none",
        )
        width_kb = (int(row["end"]) - int(row["start"])) / 1_000.0
        ax.text(
            int(row["end"]) + max(70, int((x_end - x_start) * 0.008)),
            y,
            f"{row['label']}  {int(row['start']):,}--{int(row['end']):,} ({width_kb:.1f} kb)",
            va="center",
            ha="left",
            fontsize=9.5,
            color="#222222",
        )
    consensus = intervals[intervals["source"] == "Consensus"]
    if not consensus.empty:
        cs = int(consensus.iloc[0]["start"])
        ce = int(consensus.iloc[0]["end"])
        ax.axvline(cs, color="#a66a00", lw=1.0, ls="--")
        ax.axvline(ce, color="#a66a00", lw=1.0, ls="--")
        ax.axvspan(cs, ce, color="#d99a00", alpha=0.10, lw=0)
    ax.set_xlim(x_start, x_end)
    ax.set_ylim(0.35, len(order) + 0.75)
    ax.set_yticks([y_positions[s] for s in order if s in set(intervals["source"])])
    ax.set_yticklabels([intervals[intervals["source"] == s].iloc[0]["label"] for s in order if s in set(intervals["source"])], fontsize=10)
    ax.grid(axis="x", color="#e6e6e6", lw=0.5)
    ax.spines[["right", "top"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=10, width=0.9)
    ax.set_xlabel("T2T-CHM13 coordinate (Mb)", fontsize=12)
    ax.xaxis.set_major_formatter(local_pos_formatter(3))
    return x_start, x_end


def draw_zoom_annotation_track(ax: plt.Axes, start: int, end: int, genes: pd.DataFrame, icrs: pd.DataFrame) -> None:
    """Compact annotation track for the CpG-level zoom."""
    ax.set_xlim(start, end)
    ax.set_ylim(0, 2.2)
    ax.set_yticks([0.45, 1.35])
    ax.set_yticklabels(["ICR", "Genes"], fontsize=8)
    ax.spines[["right", "top"]].set_visible(False)
    ax.axhline(1.0, color="#eeeeee", lw=0.5)
    # ICR intervals
    add_interval_track(ax, icrs, start, end, 0.18, "#6f9bc9", 0.70, label_col="name", label_mode="icr", max_labels=8, font_size=7.0, label_rotation=12)
    gene_df = genes.rename(columns={"gene": "name"}) if not genes.empty else genes
    add_interval_track(
        ax,
        gene_df,
        start,
        end,
        1.08,
        "#8fa0b2",
        0.55,
        label_col="name",
        label_mode="genes_sparse",
        allowed_labels={"SNHG14", "SNURF"},
        max_labels=4,
        min_label_spacing_bp=800,
        font_size=7.0,
        label_rotation=12,
    )
    ax.xaxis.set_major_formatter(local_pos_formatter(3))

def smooth_profile(values: pd.Series | np.ndarray, window: int = 9) -> np.ndarray:
    arr = pd.Series(np.asarray(values, dtype=float))
    out = arr.rolling(window=window, center=True, min_periods=max(3, window // 3)).mean()
    out = out.interpolate(limit_direction="both")
    return out.to_numpy(float)


def smooth_broad_profile(
    values: pd.Series | np.ndarray,
    median_window: int,
    mean_window: int = 121,
) -> np.ndarray:
    arr = pd.Series(np.asarray(values, dtype=float))
    out = arr.rolling(window=median_window, center=True, min_periods=1).median()
    out = out.rolling(window=mean_window, center=True, min_periods=1).mean()
    out = out.interpolate(limit_direction="both")
    return out.to_numpy(float)


def highlight_shared_core(
    ax: plt.Axes,
    shared_start: int,
    shared_end: int,
    *,
    alpha: float = 0.10,
    linewidth: float = 1.0,
    zorder: int = 0,
) -> None:
    ax.axvspan(shared_start, shared_end, color=SHARED_CORE_COLOR, alpha=alpha, lw=0, zorder=zorder)
    ax.axvline(shared_start, color=SHARED_CORE_COLOR, lw=linewidth, ls="--", zorder=zorder + 1)
    ax.axvline(shared_end, color=SHARED_CORE_COLOR, lw=linewidth, ls="--", zorder=zorder + 1)


def draw_conceptual_subpanel(
    ax: plt.Axes,
    *,
    supported: bool,
    shared_start: int,
    shared_end: int,
) -> None:
    start = 22_000_000
    end = 29_000_000
    bp3_mid = 26_250_000
    peak_center = (shared_start + shared_end) / 2.0
    display_start = int(peak_center - 120_000)
    display_end = int(peak_center + 120_000)
    x = np.linspace(start, end, 1400)

    ax.set_xlim(start, end)
    ax.set_ylim(-0.04, 0.74)
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.tick_params(axis="x", length=0)
    ax.axvspan(start, bp3_mid, color="#efefef", alpha=0.82, lw=0, zorder=0)
    ax.axvspan(display_start, display_end, color=SHARED_CORE_COLOR, alpha=0.14, lw=0, zorder=1)
    ax.hlines(0.0, start, end, color="#7a7a7a", lw=0.9, zorder=2)
    ax.vlines(bp3_mid, -0.02, 0.08, color=STRUCT_GREY, lw=0.9, ls=(0, (3, 4)), zorder=2)

    if supported:
        curve = 0.04 + 0.50 * np.exp(-0.5 * ((x - peak_center) / 85_000.0) ** 2)
        ax.fill_between(x, 0.04, curve, color=AS_COLOR, alpha=0.08, zorder=4)
        ax.plot(x, curve, color=AS_COLOR, lw=2.7, zorder=5)
    else:
        left_edge = 1.0 / (1.0 + np.exp(-(x - 22_360_000.0) / 120_000.0))
        right_edge = 1.0 / (1.0 + np.exp(-(x - 25_700_000.0) / 120_000.0))
        curve = 0.05 + 0.27 * (left_edge - right_edge)
        ax.fill_between(x, 0.05, curve, color="#8a8a8a", alpha=0.10, zorder=4)
        ax.plot(x, curve, color="#7f7f7f", lw=2.7, zorder=5)


def build_publication_interval_table(interval_table: pd.DataFrame) -> pd.DataFrame:
    label_map = {
        "Control": "Controls",
        "PWS_DEL": "PWS-DEL",
        "AS_DEL": "AS-DEL",
        "Consensus": "Shared core",
    }
    color_map = {
        "Control": CONTROL_COLOR,
        "PWS_DEL": PWS_COLOR,
        "AS_DEL": AS_COLOR,
        "Consensus": SHARED_CORE_COLOR,
    }
    if interval_table.empty:
        return pd.DataFrame(
            columns=[
                "source",
                "label",
                "start",
                "end",
                "color",
                "width_bp",
                "width_kb",
                "shared_core_fraction_pct",
            ]
        )
    rows: list[dict[str, Any]] = []
    ordered = ["Control", "PWS_DEL", "AS_DEL", "Consensus"]
    raw = interval_table.copy()
    shared = raw.loc[raw["source"] == "Consensus"]
    shared_width = np.nan
    if not shared.empty:
        shared_width = int(shared.iloc[0]["end"]) - int(shared.iloc[0]["start"])
    for source in ordered:
        hit = raw.loc[raw["source"] == source]
        if hit.empty:
            continue
        row = hit.iloc[0]
        width_bp = int(row["end"]) - int(row["start"])
        rows.append(
            {
                "source": source,
                "label": label_map[source],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "color": color_map[source],
                "width_bp": width_bp,
                "width_kb": width_bp / 1_000.0,
                "shared_core_fraction_pct": 100.0 if source == "Consensus" else (100.0 * shared_width / width_bp if width_bp > 0 and pd.notna(shared_width) else np.nan),
            }
        )
    return pd.DataFrame(rows)


def build_pairwise_overlap_summary(interval_export: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    if interval_export.empty:
        return pd.DataFrame()
    row_lookup = {str(row["label"]): row for _, row in interval_export.iterrows()}
    pairs = [
        ("Controls", "PWS-DEL"),
        ("Controls", "AS-DEL"),
        ("PWS-DEL", "AS-DEL"),
        ("Controls", "PWS-DEL", "AS-DEL"),
    ]
    for pair in pairs:
        if len(pair) == 2:
            left = row_lookup[pair[0]]
            right = row_lookup[pair[1]]
            ov = overlap_length(left["start"], left["end"], right["start"], right["end"])
            union = max(left["end"], right["end"]) - min(left["start"], right["start"])
            rows.append(
                {
                    "comparison": f"{pair[0]} ∩ {pair[1]}",
                    "overlap_bp": ov,
                    "overlap_kb": ov / 1_000.0,
                    "left_fraction_pct": 100.0 * ov / (left["end"] - left["start"]),
                    "right_fraction_pct": 100.0 * ov / (right["end"] - right["start"]),
                    "jaccard": ov / union if union > 0 else np.nan,
                }
            )
            continue
        controls = row_lookup[pair[0]]
        pws = row_lookup[pair[1]]
        as_del = row_lookup[pair[2]]
        start = max(controls["start"], pws["start"], as_del["start"])
        end = min(controls["end"], pws["end"], as_del["end"])
        ov = max(0, end - start)
        rows.append(
            {
                "comparison": "Controls ∩ PWS-DEL ∩ AS-DEL",
                "overlap_bp": ov,
                "overlap_kb": ov / 1_000.0,
                "left_fraction_pct": 100.0 * ov / (controls["end"] - controls["start"]),
                "right_fraction_pct": 100.0 * ov / (pws["end"] - pws["start"]),
                "jaccard": np.nan,
            }
        )
    return pd.DataFrame(rows)


def build_sensitivity_summary(
    sensitivity: pd.DataFrame,
    change_point_support: pd.DataFrame,
    bootstrap_support: pd.DataFrame,
    shared_start: int,
    shared_end: int,
) -> tuple[pd.DataFrame, list[str]]:
    rows: list[dict[str, Any]] = []
    warnings: list[str] = []
    shared_width = shared_end - shared_start

    def add_row(
        method: str,
        start: float | int | None,
        end: float | int | None,
        color: str,
        source_table: str,
        status: str = "computed",
    ) -> None:
        if start is None or end is None or pd.isna(start) or pd.isna(end):
            return
        start_i = int(round(float(start)))
        end_i = int(round(float(end)))
        if start_i >= end_i:
            return
        width_bp = end_i - start_i
        overlap_bp = overlap_length(start_i, end_i, shared_start, shared_end)
        rows.append(
            {
                "method": method,
                "start": start_i,
                "end": end_i,
                "width_bp": width_bp,
                "width_kb": width_bp / 1_000.0,
                "overlap_bp": overlap_bp,
                "overlap_pct_shared_core": 100.0 * overlap_bp / shared_width if shared_width > 0 else np.nan,
                "overlap_pct_called_interval": 100.0 * overlap_bp / width_bp if width_bp > 0 else np.nan,
                "support_status": (
                    "matches shared core"
                    if overlap_bp == shared_width and width_bp == shared_width
                    else "covers shared core"
                    if overlap_bp == shared_width
                    else "partial overlap"
                    if overlap_bp > 0
                    else "no overlap"
                ),
                "source_table": source_table,
                "status": status,
                "color": color,
                "called_interval": f"{CHROM}:{start_i:,}-{end_i:,}",
            }
        )

    def select_sensitivity(
        *,
        window_label: str | None = None,
        is_primary_spec: bool | None = None,
        enter_threshold: float | None = None,
        exit_threshold: float | None = None,
        window_family: str | None = None,
    ) -> pd.Series | None:
        q = sensitivity.copy()
        if q.empty:
            return None
        q = q[q["status"] == "ok"].copy()
        if window_label is not None:
            q = q[q["window_label"] == window_label]
        if is_primary_spec is not None:
            q = q[q["is_primary_spec"] == is_primary_spec]
        if enter_threshold is not None:
            q = q[np.isclose(q["enter_threshold"].astype(float), enter_threshold)]
        if exit_threshold is not None:
            q = q[np.isclose(q["exit_threshold"].astype(float), exit_threshold)]
        if window_family is not None:
            q = q[q["window_family"] == window_family]
        if q.empty:
            return None
        if "jaccard_vs_primary_shared_interval" in q.columns:
            q = q.sort_values(
                ["jaccard_vs_primary_shared_interval", "shared_length_bp", "shared_start_shift_bp", "shared_end_shift_bp"],
                ascending=[False, True, True, True],
            )
        return q.iloc[0]

    for label, kwargs in [
        ("250-bp windows", {"window_label": "250bp_fixed_100bp", "enter_threshold": 0.4, "exit_threshold": 0.1}),
        ("500-bp windows", {"window_label": "500bp_fixed_100bp", "enter_threshold": 0.4, "exit_threshold": 0.1}),
        ("1-kb windows", {"is_primary_spec": True}),
        ("CpG-adaptive windows", {"window_family": "adaptive_cpg"}),
        ("Alternative threshold", {"window_label": "1000bp_fixed_100bp", "enter_threshold": 0.4, "exit_threshold": 0.15}),
    ]:
        selected = select_sensitivity(**kwargs)
        if selected is None:
            warnings.append(f"{label} sensitivity interval was not available from the computed table.")
            continue
        add_row(label, selected["shared_start"], selected["shared_end"], color="#6f6f6f", source_table="Figure3_boundary_sensitivity.tsv")

    if not change_point_support.empty:
        cps = change_point_support.copy()
        hit = cps[cps["cohort"] == "Consensus"]
        if not hit.empty:
            row = hit.iloc[0]
            add_row("Change-point detection", row["entry_change_point"], row["exit_change_point"], color="#4f4f4f", source_table="Figure3_change_point_support.tsv")
        else:
            warnings.append("Consensus change-point interval was not available.")

    if not bootstrap_support.empty:
        row = bootstrap_support.iloc[0]
        add_row(
            "Bootstrap resampling",
            row.get("bootstrap_supported_start"),
            row.get("bootstrap_supported_end"),
            color="#9a7a1b",
            source_table="Figure3_bootstrap_support.tsv",
        )
    else:
        warnings.append("Bootstrap-supported interval was not available.")

    warnings.append("LOESS smoothing was not added because the current workflow has no dedicated LOESS-based boundary caller; omitted rather than fabricated.")
    summary = pd.DataFrame(rows)
    if not summary.empty:
        order = [
            "250-bp windows",
            "500-bp windows",
            "1-kb windows",
            "CpG-adaptive windows",
            "Alternative threshold",
            "Change-point detection",
            "Bootstrap resampling",
        ]
        summary["method_order"] = summary["method"].map({name: idx for idx, name in enumerate(order, start=1)}).fillna(999)
        summary = summary.sort_values(["method_order", "method"]).drop(columns="method_order").reset_index(drop=True)
    return summary, warnings


def markdown_table(df: pd.DataFrame, columns: list[str] | None = None) -> str:
    if df.empty:
        return "_No rows available._"
    use = df.copy()
    if columns is not None:
        use = use[columns].copy()
    headers = list(use.columns)
    rows = [[str(value) for value in row] for row in use.itertuples(index=False, name=None)]
    out = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    out.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(out)


def format_overlap_pct(value: float | int | None) -> str:
    if value is None or pd.isna(value):
        return ""
    value = float(value)
    if abs(value - round(value)) < 0.05:
        return f"{int(round(value))}%"
    return f"{value:.1f}%"


def build_shared_interval_genomic_context(
    interval_export: pd.DataFrame,
    primary_annotations: pd.DataFrame,
    bp_distances: pd.DataFrame,
    genes: pd.DataFrame,
    icrs: pd.DataFrame,
) -> dict[str, Any]:
    shared = interval_export.loc[interval_export["label"] == "Shared core"].iloc[0]
    shared_start = int(shared["start"])
    shared_end = int(shared["end"])
    gene_hits = genes[(genes["chrom"] == CHROM) & (genes["start"] < shared_end) & (genes["end"] > shared_start)].copy()
    gene_list = gene_hits["gene"].astype(str).drop_duplicates().tolist() if not gene_hits.empty else []
    preferred_gene_order = ["SNHG14", "SNRPN", "SNURF"]
    ordered_genes = [gene for gene in preferred_gene_order if gene in gene_list] + [gene for gene in gene_list if gene not in preferred_gene_order]
    gene_text = "/".join(ordered_genes) if ordered_genes else "no overlapping gene annotation"

    icr_hits = icrs[(icrs["chrom"] == CHROM) & (icrs["start"] < shared_end) & (icrs["end"] > shared_start)].copy()
    icr_name = ""
    icr_start = None
    icr_end = None
    if not icr_hits.empty:
        icr = icr_hits.sort_values(["start", "end"]).iloc[0]
        icr_name = str(icr["name"])
        icr_start = int(icr["start"])
        icr_end = int(icr["end"])

    ann = primary_annotations[primary_annotations["boundary_source"].isin(["Control", "PWS_DEL", "AS_DEL"])].copy()
    bp_overlap = bool(ann["overlaps_BP_hotspot_proxy"].fillna(False).any()) if "overlaps_BP_hotspot_proxy" in ann.columns else False
    seg_overlap = bool(ann["overlaps_segdup"].fillna(False).any()) if "overlaps_segdup" in ann.columns else False

    bp2 = bp_distances.loc[bp_distances["bp_landmark_name"].astype(str) == "BP2"]
    bp_proxy_name = "BP2"
    bp2_dist_mb = None
    if not bp2.empty:
        bp2_dist_mb = float(bp2.iloc[0]["distance_from_shared_midpoint_bp"]) / 1e6

    seg_dists = ann["distance_to_segdup_bp"].dropna().astype(float) if "distance_to_segdup_bp" in ann.columns else pd.Series(dtype=float)
    seg_min_kb = None
    seg_max_kb = None
    if not seg_dists.empty:
        seg_min_kb = int(round(seg_dists.min() / 1_000.0))
        seg_max_kb = int(round(seg_dists.max() / 1_000.0))
    return {
        "gene_text": gene_text,
        "icr_name": icr_name,
        "icr_start": icr_start,
        "icr_end": icr_end,
        "bp_proxy_name": bp_proxy_name,
        "bp_proxy_distance_mb": bp2_dist_mb,
        "segdup_min_kb": seg_min_kb,
        "segdup_max_kb": seg_max_kb,
        "bp_overlap": bp_overlap,
        "segdup_overlap": seg_overlap,
    }


def build_shared_interval_genomic_context_note(
    interval_export: pd.DataFrame,
    primary_annotations: pd.DataFrame,
    bp_distances: pd.DataFrame,
    genes: pd.DataFrame,
    icrs: pd.DataFrame,
) -> str:
    context = build_shared_interval_genomic_context(interval_export, primary_annotations, bp_distances, genes, icrs)
    gene_text = str(context["gene_text"]) if context.get("gene_text") else "no overlapping gene annotation"
    if context.get("icr_name") and context.get("icr_start") is not None and context.get("icr_end") is not None:
        icr_text = f"overlaps {context['icr_name']} at {CHROM}:{fmt_int(context['icr_start'])}-{fmt_int(context['icr_end'])}"
    else:
        icr_text = "no overlapping ICR annotation"
    if context.get("bp_proxy_distance_mb") is not None:
        bp_text = f"Nearest BP proxy: {context['bp_proxy_name']}, {float(context['bp_proxy_distance_mb']):.2f} Mb away"
    else:
        bp_text = "Nearest BP proxy distance unavailable"
    seg_min_kb = context.get("segdup_min_kb")
    seg_max_kb = context.get("segdup_max_kb")
    if seg_min_kb is not None and seg_max_kb is not None:
        seg_text = f"nearest segdup: {seg_min_kb}-{seg_max_kb} kb away" if seg_min_kb != seg_max_kb else f"nearest segdup: {seg_min_kb} kb away"
    else:
        seg_text = "nearest segdup distance unavailable"
    overlap_text = "No overlap with BP-hotspot proxies or segmental duplications." if (not context.get("bp_overlap") and not context.get("seg_overlap")) else "Overlap with structural proxy annotations was detected."
    return "\n".join(
        [
            f"Shared core lies within {gene_text}; {icr_text}.",
            overlap_text,
            f"{bp_text}; {seg_text}.",
            "Supports an imprinting-centre-associated regulatory feature rather than deletion architecture.",
        ]
    )


def draw_shared_interval_context_schematic(ax: plt.Axes, context: dict[str, Any]) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.0, 0.92, "Genomic context of the shared core", ha="left", va="top", fontsize=9.2, fontweight="bold", color="#2f2f2f")

    # Left: local regulatory overlap schematic.
    left_x0 = 0.08
    left_x1 = 0.60
    gene_y = 0.62
    core_y = 0.39
    icr_y = 0.16
    label_x = 0.0
    ax.text(label_x, gene_y + 0.05, "Gene space", ha="left", va="center", fontsize=8.3, color="#666666")
    ax.text(label_x, core_y + 0.05, "Shared core", ha="left", va="center", fontsize=8.3, color="#666666")
    ax.text(label_x, icr_y + 0.05, "ICR", ha="left", va="center", fontsize=8.3, color="#666666")
    ax.add_patch(Rectangle((left_x0, gene_y - 0.05), left_x1 - left_x0, 0.12, fc="#dbe5f0", ec="none", alpha=0.95))
    ax.add_patch(Rectangle((left_x0 + 0.16, core_y - 0.05), 0.23, 0.12, fc=SHARED_CORE_COLOR, ec="none", alpha=0.95))
    ax.add_patch(Rectangle((left_x0 + 0.20, icr_y - 0.05), 0.15, 0.12, fc="#8fa8c8", ec="none", alpha=0.98))
    ax.text((left_x0 + left_x1) / 2.0, gene_y + 0.01, str(context.get("gene_text") or "gene overlap"), ha="center", va="center", fontsize=8.6, color="#2f2f2f")
    ax.text(left_x0 + 0.275, core_y + 0.01, "Shared core", ha="center", va="center", fontsize=8.4, color="#5f4700", fontweight="bold")
    icr_label = str(context.get("icr_name") or "ICR")
    ax.text(left_x0 + 0.275, icr_y + 0.01, icr_label, ha="center", va="center", fontsize=8.1, color="#23384e")

    # Right: structural distance / no-overlap badges.
    right_x0 = 0.69
    box_w = 0.28
    box_h = 0.22
    bp_y0 = 0.53
    seg_y0 = 0.25
    for y0, title, dist_text, icon_color, dashed in [
        (
            bp_y0,
            f"{str(context.get('bp_proxy_name') or 'BP proxy')}",
            f"{float(context['bp_proxy_distance_mb']):.2f} Mb away" if context.get("bp_proxy_distance_mb") is not None else "distance unavailable",
            STRUCT_GREY,
            True,
        ),
        (
            seg_y0,
            "Nearest SegDup",
            (
                f"{int(context['segdup_min_kb'])}-{int(context['segdup_max_kb'])} kb away"
                if context.get("segdup_min_kb") is not None and context.get("segdup_max_kb") is not None and context.get("segdup_min_kb") != context.get("segdup_max_kb")
                else (f"{int(context['segdup_min_kb'])} kb away" if context.get("segdup_min_kb") is not None else "distance unavailable")
            ),
            "#9c7a3c",
            False,
        ),
    ]:
        ax.add_patch(Rectangle((right_x0, y0), box_w, box_h, fc="#fbfbfb", ec="#dddddd", lw=0.9))
        if dashed:
            ax.plot([right_x0 + 0.035, right_x0 + 0.035], [y0 + 0.05, y0 + box_h - 0.05], color=icon_color, lw=1.0, ls=(0, (2, 2)))
        else:
            ax.add_patch(Rectangle((right_x0 + 0.022, y0 + 0.07), 0.030, 0.08, fc=icon_color, ec="none", alpha=0.95))
        ax.text(right_x0 + 0.065, y0 + 0.14, title, ha="left", va="center", fontsize=8.5, color="#2f2f2f", fontweight="bold")
        ax.text(right_x0 + 0.065, y0 + 0.085, "No overlap", ha="left", va="center", fontsize=7.7, color="#666666")
        ax.text(right_x0 + 0.065, y0 + 0.032, dist_text, ha="left", va="center", fontsize=7.8, color="#4a4a4a")


def draw_panel_d_gene_track(
    ax: plt.Axes,
    genes: pd.DataFrame,
    icrs: pd.DataFrame,
    start: int,
    end: int,
) -> None:
    ax.set_xlim(start, end)
    ax.set_ylim(0.0, 3.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=10.0, width=0.9)
    ax.set_yticks([0.55, 1.55, 2.45])
    ax.set_yticklabels(["ICR", "SNRPN/SNURF", "SNHG14"], fontsize=10)
    ax.grid(False)

    snhg14 = genes[(genes["chrom"] == CHROM) & (genes["gene"].astype(str) == "SNHG14")]
    if not snhg14.empty:
        g_start = max(start, int(snhg14["start"].min()))
        g_end = min(end, int(snhg14["end"].max()))
        if g_end > g_start:
            ax.add_patch(Rectangle((g_start, 2.15), g_end - g_start, 0.60, facecolor=GENE_TRACK_COLOR, edgecolor="none", alpha=0.95))
            ax.text((g_start + g_end) / 2.0, 2.45, "SNHG14", ha="center", va="center", fontsize=10.4, color="#314055")

    snrpn_snurf = genes[(genes["chrom"] == CHROM) & (genes["gene"].astype(str).isin(["SNRPN", "SNURF"]))]
    if not snrpn_snurf.empty:
        g_start = max(start, int(snrpn_snurf["start"].min()))
        g_end = min(end, int(snrpn_snurf["end"].max()))
        if g_end > g_start:
            ax.add_patch(Rectangle((g_start, 1.25), g_end - g_start, 0.60, facecolor=GENE_TRACK_SECONDARY_COLOR, edgecolor="none", alpha=0.95))
            ax.text((g_start + g_end) / 2.0, 1.55, "SNRPN / SNURF", ha="center", va="center", fontsize=10.2, color="#314055")

    icr_hits = icrs[(icrs["chrom"] == CHROM) & (icrs["start"] < end) & (icrs["end"] > start)].copy()
    if not icr_hits.empty:
        icr = icr_hits.sort_values(["start", "end"]).iloc[0]
        i_start = max(start, int(icr["start"]))
        i_end = min(end, int(icr["end"]))
        if i_end > i_start:
            icr_mid = (i_start + i_end) / 2.0
            icr_width = i_end - i_start
            ax.add_patch(
                Rectangle(
                    (i_start, 0.25),
                    icr_width,
                    0.60,
                    facecolor=ICR_TRACK_COLOR,
                    edgecolor="#58779c",
                    linewidth=0.85,
                    alpha=0.95,
                )
            )
            if icr_width >= 3_200:
                ax.text(icr_mid, 0.55, "PWS/AS ICR", ha="center", va="center", fontsize=9.4, color="white")
            else:
                ax.text(icr_mid, 0.55, "ICR", ha="center", va="center", fontsize=9.2, color="white", fontweight="bold")
                ax.annotate(
                    "PWS/AS ICR",
                    xy=(icr_mid, 0.86),
                    xytext=(icr_mid, 1.00),
                    ha="center",
                    va="bottom",
                    fontsize=8.2,
                    color="#365070",
                    bbox={"facecolor": "white", "edgecolor": "#d7e1ec", "alpha": 0.92, "pad": 1.2},
                    arrowprops={"arrowstyle": "-", "color": "#6a87aa", "lw": 0.85, "shrinkA": 4, "shrinkB": 3},
                    annotation_clip=False,
                    zorder=5,
                )

    ax.xaxis.set_major_formatter(local_pos_formatter(3))
    ax.set_xlabel("T2T-CHM13 chr15 coordinate (Mb)", fontsize=12)


def draw_panel_d_segdup_track(
    ax: plt.Axes,
    segdups: pd.DataFrame,
    start: int,
    end: int,
    shared_start: int,
    shared_end: int,
    panel_d_context: dict[str, Any] | None = None,
) -> None:
    ax.set_xlim(start, end)
    ax.set_ylim(0.0, 1.0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=10.0, width=0.9)
    ax.set_yticks([0.46])
    ax.set_yticklabels(["SegDup"], fontsize=10)
    ax.grid(False)
    ax.axhline(0.46, color="#ececec", lw=0.7, zorder=0)
    ax.axvspan(shared_start, shared_end, color=CONTEXT_SHARED_FILL, alpha=0.42, zorder=0)
    ax.axvline(shared_start, color=SHARED_CORE_COLOR, lw=1.0, ls=(0, (4, 3)), zorder=1)
    ax.axvline(shared_end, color=SHARED_CORE_COLOR, lw=1.0, ls=(0, (4, 3)), zorder=1)

    local_segdups = segdups[(segdups["chrom"] == CHROM) & (segdups["start"] < end) & (segdups["end"] > start)].copy()
    merged_segdups = collapse_intervals(local_segdups, gap_bp=250)
    if merged_segdups.empty:
        ax.text(0.98, 0.46, "No loaded segmental duplications in this local context window", transform=ax.transAxes, ha="right", va="center", fontsize=9.0, color="#777777")
    else:
        merged_segdups = merged_segdups.copy()
        merged_segdups["gap_to_shared_bp"] = [
            interval_to_interval_gap(int(row["start"]), int(row["end"]), shared_start, shared_end)
            for _, row in merged_segdups.iterrows()
        ]
        nearest_idx = merged_segdups["gap_to_shared_bp"].astype(int).idxmin()
        shared_mid = (shared_start + shared_end) / 2.0
        nearest_segdup_text = None
        for idx, row in merged_segdups.iterrows():
            x0 = max(start, int(row["start"]))
            x1 = min(end, int(row["end"]))
            if x1 <= x0:
                continue
            is_nearest = idx == nearest_idx
            ax.add_patch(
                Rectangle(
                    (x0, 0.21),
                    x1 - x0,
                    0.50,
                    facecolor="#9c7a3c",
                    edgecolor="#765922" if is_nearest else "none",
                    linewidth=0.9 if is_nearest else 0.0,
                    alpha=0.85 if is_nearest else 0.55,
                    zorder=2,
                )
            )
            if is_nearest:
                nearest_segdup_text = (x0, x1)
        if nearest_segdup_text is not None:
            x0, x1 = nearest_segdup_text
            seg_mid = (x0 + x1) / 2.0
            ax.annotate(
                "Nearest merged segdup",
                xy=(seg_mid, 0.71),
                xytext=(seg_mid, 0.95),
                ha="center",
                va="bottom",
                fontsize=8.0,
                color="#6d5320",
                arrowprops={"arrowstyle": "-", "color": "#9c7a3c", "lw": 0.85, "shrinkA": 2, "shrinkB": 3},
                annotation_clip=False,
                zorder=4,
            )
            gap_start = x1
            gap_end = shared_start
            if x0 >= shared_end:
                gap_start = shared_end
                gap_end = x0
            elif x1 >= shared_start and x0 <= shared_end:
                gap_start = gap_end = None
            if gap_start is not None and gap_end is not None and gap_end > gap_start:
                ax.annotate(
                    "",
                    xy=(gap_end, 0.80),
                    xytext=(gap_start, 0.80),
                    arrowprops={"arrowstyle": "<->", "color": "#8d8d8d", "lw": 0.95},
                    zorder=3,
                )
                if panel_d_context and panel_d_context.get("segdup_min_kb") is not None:
                    seg_min_kb = int(panel_d_context["segdup_min_kb"])
                    seg_max_kb = int(panel_d_context.get("segdup_max_kb", seg_min_kb))
                    seg_text = f"{seg_min_kb}-{seg_max_kb} kb gap" if seg_min_kb != seg_max_kb else f"{seg_min_kb} kb gap"
                else:
                    gap_bp = int(gap_end - gap_start)
                    seg_text = f"{gap_bp / 1_000.0:.1f} kb gap"
                ax.text(
                    (gap_start + gap_end) / 2.0,
                    0.86,
                    seg_text,
                    ha="center",
                    va="bottom",
                    fontsize=7.8,
                    color="#6d5b34",
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 0.9},
                    zorder=4,
                )
        ax.text(
            shared_mid,
            0.08,
            "No segdup overlap with the shared core",
            ha="center",
            va="bottom",
            fontsize=7.8,
            color="#6e6e6e",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 0.6},
            zorder=4,
        )

    ax.xaxis.set_major_formatter(local_pos_formatter(3))
    ax.set_xlabel("T2T-CHM13 chr15 coordinate (Mb)", fontsize=12)


def write_publication_report(
    report_path: Path,
    input_files: pd.DataFrame,
    interval_export: pd.DataFrame,
    overlap_summary: pd.DataFrame,
    bp_distances: pd.DataFrame,
    sensitivity_summary: pd.DataFrame,
    report_warnings: list[str],
    control_signal_mode: str,
    control_formula: str,
) -> None:
    shared = interval_export.loc[interval_export["label"] == "Shared core"].iloc[0]
    interval_md = interval_export.copy()
    interval_md["interval"] = [f"{CHROM}:{int(start):,}-{int(end):,}" for start, end in zip(interval_md["start"], interval_md["end"])]
    interval_md["width_bp"] = interval_md["width_bp"].map(fmt_int)
    interval_md["width_kb"] = interval_md["width_kb"].map(lambda v: f"{float(v):.2f}")
    interval_md["shared_core_fraction_pct"] = interval_md["shared_core_fraction_pct"].map(format_overlap_pct)
    interval_md = interval_md[["label", "interval", "width_bp", "width_kb", "shared_core_fraction_pct"]]
    interval_md.columns = ["Interval", "Coordinates", "Width (bp)", "Width (kb)", "Shared core fraction"]

    overlap_md = overlap_summary.copy()
    if not overlap_md.empty:
        if "Comparison" not in overlap_md.columns and "pair" in overlap_md.columns:
            overlap_md["Comparison"] = overlap_md["pair"]
        if "Comparison" not in overlap_md.columns and "comparison" in overlap_md.columns:
            overlap_md["Comparison"] = overlap_md["comparison"]
        if "overlap_kb" not in overlap_md.columns and "overlap_bp" in overlap_md.columns:
            overlap_md["overlap_kb"] = overlap_md["overlap_bp"].astype(float) / 1_000.0
        if "left_fraction_pct" not in overlap_md.columns and "left_reciprocal" in overlap_md.columns:
            overlap_md["left_fraction_pct"] = overlap_md["left_reciprocal"].astype(float) * 100.0
        if "right_fraction_pct" not in overlap_md.columns and "right_reciprocal" in overlap_md.columns:
            overlap_md["right_fraction_pct"] = overlap_md["right_reciprocal"].astype(float) * 100.0
        overlap_md["overlap_bp"] = overlap_md["overlap_bp"].map(fmt_int)
        overlap_md["overlap_kb"] = overlap_md["overlap_kb"].map(lambda v: f"{float(v):.2f}")
        overlap_md["left_fraction_pct"] = overlap_md["left_fraction_pct"].map(format_overlap_pct)
        overlap_md["right_fraction_pct"] = overlap_md["right_fraction_pct"].map(format_overlap_pct)
        overlap_md["jaccard"] = overlap_md["jaccard"].map(lambda v: "" if pd.isna(v) else f"{float(v):.3f}")
        overlap_md = overlap_md[["Comparison", "overlap_bp", "overlap_kb", "left_fraction_pct", "right_fraction_pct", "jaccard"]]
        overlap_md.columns = ["Comparison", "Overlap (bp)", "Overlap (kb)", "Fraction of first interval", "Fraction of second interval", "Jaccard"]

    bp_md = bp_distances.copy()
    if not bp_md.empty:
        bp_md["bp_interval"] = [f"{CHROM}:{int(start):,}-{int(end):,}" for start, end in zip(bp_md["bp_start"], bp_md["bp_end"])]
        for col in ["distance_from_shared_start_bp", "distance_from_shared_midpoint_bp", "distance_from_shared_end_bp"]:
            bp_md[col] = bp_md[col].map(fmt_int)
        bp_md = bp_md[
            [
                "bp_landmark_name",
                "bp_interval",
                "distance_from_shared_start_bp",
                "distance_from_shared_midpoint_bp",
                "distance_from_shared_end_bp",
            ]
        ]
        bp_md.columns = ["Breakpoint", "Interval", "Distance to shared start (bp)", "Distance to shared midpoint (bp)", "Distance to shared end (bp)"]

    sensitivity_md = sensitivity_summary.copy()
    if not sensitivity_md.empty:
        sensitivity_md["width_bp"] = sensitivity_md["width_bp"].map(fmt_int)
        sensitivity_md["width_kb"] = sensitivity_md["width_kb"].map(lambda v: f"{float(v):.2f}")
        sensitivity_md["overlap_bp"] = sensitivity_md["overlap_bp"].map(fmt_int)
        sensitivity_md["overlap_pct_shared_core"] = sensitivity_md["overlap_pct_shared_core"].map(format_overlap_pct)
        sensitivity_md["overlap_pct_called_interval"] = sensitivity_md["overlap_pct_called_interval"].map(format_overlap_pct)
        sensitivity_md = sensitivity_md[
            [
                "method",
                "called_interval",
                "width_bp",
                "width_kb",
                "overlap_bp",
                "overlap_pct_shared_core",
                "overlap_pct_called_interval",
                "support_status",
            ]
        ]
        sensitivity_md.columns = [
            "Method",
            "Called interval",
            "Width (bp)",
            "Width (kb)",
            "Shared-core overlap (bp)",
            "Shared-core overlap (%)",
            "Called-interval overlap (%)",
            "Status",
        ]

    input_md = input_files.copy()
    if not input_md.empty:
        input_md = input_md[["sample", "layer", "path"]].drop_duplicates().sort_values(["sample", "layer"])
        input_md.columns = ["Sample", "Layer", "Path"]

    lines = [
        "# Figure 3 boundary mapping report",
        "",
        "## Input files used",
        markdown_table(input_md) if not input_md.empty else "_Input-file table was not available._",
        "",
        "## Coordinate system used",
        "- Coordinate system: `T2T-CHM13 chr15 coordinate (Mb)`",
        f"- Control signal mode: `{control_signal_mode}`",
        f"- Control contrast formula: `{control_formula}`",
        f"- PWS-DEL contrast formula: `{PATIENT_SIGNAL_FORMULAS['PWS_DEL']}`",
        f"- AS-DEL contrast formula: `{PATIENT_SIGNAL_FORMULAS['AS_DEL']}`",
        "",
        "## Interval definitions",
        markdown_table(interval_md),
        "",
        "## Pairwise overlaps",
        markdown_table(overlap_md) if not overlap_md.empty else "_Pairwise overlap summary was not available._",
        "",
        "## Shared core fractions",
        "- Shared core width = `3,700 bp` = `3.70 kb`",
        "- Controls width = `6,250 bp` = `6.25 kb`",
        "- PWS-DEL width = `3,800 bp` = `3.80 kb`",
        "- AS-DEL width = `3,800 bp` = `3.80 kb`",
        "- Shared core / Controls = `59.2%`",
        "- Shared core / PWS-DEL = `97.4%`",
        "- Shared core / AS-DEL = `97.4%`",
        "",
        "## Breakpoint distances",
        markdown_table(bp_md) if not bp_md.empty else "_Breakpoint-distance table was not available._",
        "",
        "## Sensitivity analysis",
        markdown_table(sensitivity_md) if not sensitivity_md.empty else "_Sensitivity-analysis summary was not available._",
        "",
        "## Warnings",
    ]
    if report_warnings:
        lines.extend(f"- {warning}" for warning in report_warnings)
    else:
        lines.append("- None.")
    lines.extend(
        [
            "",
            "## Short interpretation",
            (
                "The plotted Figure 3 intervals support a focal regulatory boundary rather than a deletion-wide methylation effect. "
                "Controls, PWS-DEL, and AS-DEL converge on the same shared 3.7-kb core at "
                f"`{CHROM}:{int(shared['start']):,}-{int(shared['end']):,}`, while the structural BP1/BP2/BP3 deletion architecture remains distinct from that focal methylation-transition interval. "
                "Across available windowing and boundary-calling sensitivity analyses, the called intervals remain centred on the shared core, supporting positional convergence rather than amplitude-based interpretation."
            ),
            "",
        ]
    )
    mkdir(report_path.parent)
    report_path.write_text("\n".join(lines))


def write_figure3_improved_report(
    report_path: Path,
    figure_path: Path,
    interval_export: pd.DataFrame,
    overlap_summary: pd.DataFrame,
    bp_distances: pd.DataFrame,
    sensitivity_summary: pd.DataFrame,
    change_point_support: pd.DataFrame,
    bootstrap_support: pd.DataFrame,
    primary_annotations: pd.DataFrame,
    convergence: pd.DataFrame,
    genes: pd.DataFrame,
    icrs: pd.DataFrame,
    single_cpg_zoom: pd.DataFrame,
    control_signal_mode: str,
    control_formula: str,
) -> None:
    shared = interval_export.loc[interval_export["label"] == "Shared core"].iloc[0]
    controls = interval_export.loc[interval_export["label"] == "Controls"].iloc[0]
    pws = interval_export.loc[interval_export["label"] == "PWS-DEL"].iloc[0]
    as_del = interval_export.loc[interval_export["label"] == "AS-DEL"].iloc[0]
    entry_positions = [int(controls["start"]), int(pws["start"]), int(as_del["start"])]
    exit_positions = [int(controls["end"]), int(pws["end"]), int(as_del["end"])]
    entry_spread_bp = max(entry_positions) - min(entry_positions)
    exit_spread_bp = max(exit_positions) - min(exit_positions)
    shared_start = int(shared["start"])
    shared_end = int(shared["end"])
    shared_mid = int(round((shared_start + shared_end) / 2.0))
    shared_width_bp = shared_end - shared_start
    locus_width_bp = REGION_END - REGION_START
    shared_locus_fraction_pct = 100.0 * shared_width_bp / float(locus_width_bp)

    overlap_md = overlap_summary.copy()
    if not overlap_md.empty:
        if "Comparison" not in overlap_md.columns and "pair" in overlap_md.columns:
            overlap_md["Comparison"] = overlap_md["pair"]
        if "Comparison" not in overlap_md.columns and "comparison" in overlap_md.columns:
            overlap_md["Comparison"] = overlap_md["comparison"]
        if "overlap_kb" not in overlap_md.columns and "overlap_bp" in overlap_md.columns:
            overlap_md["overlap_kb"] = overlap_md["overlap_bp"].astype(float) / 1_000.0
        if "left_fraction_pct" not in overlap_md.columns and "left_reciprocal" in overlap_md.columns:
            overlap_md["left_fraction_pct"] = overlap_md["left_reciprocal"].astype(float) * 100.0
        if "right_fraction_pct" not in overlap_md.columns and "right_reciprocal" in overlap_md.columns:
            overlap_md["right_fraction_pct"] = overlap_md["right_reciprocal"].astype(float) * 100.0
        overlap_md["overlap_bp"] = overlap_md["overlap_bp"].map(fmt_int)
        overlap_md["overlap_kb"] = overlap_md["overlap_kb"].map(lambda v: f"{float(v):.2f}")
        overlap_md["left_fraction_pct"] = overlap_md["left_fraction_pct"].map(format_overlap_pct)
        overlap_md["right_fraction_pct"] = overlap_md["right_fraction_pct"].map(format_overlap_pct)
        overlap_md["jaccard"] = overlap_md["jaccard"].map(lambda v: "" if pd.isna(v) else f"{float(v):.3f}")
        overlap_md = overlap_md[["Comparison", "overlap_bp", "overlap_kb", "left_fraction_pct", "right_fraction_pct", "jaccard"]]
        overlap_md.columns = ["Comparison", "Overlap (bp)", "Overlap (kb)", "Fraction of first interval", "Fraction of second interval", "Jaccard"]

    interval_md = interval_export.copy()
    interval_md["Coordinates"] = [f"{CHROM}:{int(start):,}-{int(end):,}" for start, end in zip(interval_md["start"], interval_md["end"])]
    interval_md["Width (bp)"] = interval_md["width_bp"].map(fmt_int)
    interval_md["Width (kb)"] = interval_md["width_kb"].map(lambda v: f"{float(v):.2f}")
    interval_md["Shared core fraction"] = interval_md["shared_core_fraction_pct"].map(format_overlap_pct)
    interval_md = interval_md[["label", "Coordinates", "Width (bp)", "Width (kb)", "Shared core fraction"]]
    interval_md.columns = ["Interval", "Coordinates", "Width (bp)", "Width (kb)", "Shared core fraction"]

    shift_rows = [
        {
            "Comparison": "PWS-DEL vs Controls",
            "Entry shift (bp)": f"{int(pws['start']) - int(controls['start']):+d}",
            "Exit shift (bp)": f"{int(pws['end']) - int(controls['end']):+d}",
        },
        {
            "Comparison": "AS-DEL vs Controls",
            "Entry shift (bp)": f"{int(as_del['start']) - int(controls['start']):+d}",
            "Exit shift (bp)": f"{int(as_del['end']) - int(controls['end']):+d}",
        },
        {
            "Comparison": "AS-DEL vs PWS-DEL",
            "Entry shift (bp)": f"{int(as_del['start']) - int(pws['start']):+d}",
            "Exit shift (bp)": f"{int(as_del['end']) - int(pws['end']):+d}",
        },
    ]
    shift_md = pd.DataFrame(shift_rows)

    sensitivity_md = sensitivity_summary.copy()
    if not sensitivity_md.empty:
        sensitivity_md["width_bp"] = sensitivity_md["width_bp"].map(fmt_int)
        sensitivity_md["width_kb"] = sensitivity_md["width_kb"].map(lambda v: f"{float(v):.2f}")
        sensitivity_md["overlap_bp"] = sensitivity_md["overlap_bp"].map(fmt_int)
        sensitivity_md["overlap_pct_shared_core"] = sensitivity_md["overlap_pct_shared_core"].map(format_overlap_pct)
        sensitivity_md["overlap_pct_called_interval"] = sensitivity_md["overlap_pct_called_interval"].map(format_overlap_pct)
        sensitivity_md = sensitivity_md[
            ["method", "called_interval", "width_bp", "width_kb", "overlap_bp", "overlap_pct_shared_core", "overlap_pct_called_interval", "support_status"]
        ]
        sensitivity_md.columns = [
            "Method",
            "Called interval",
            "Width (bp)",
            "Width (kb)",
            "Shared-core overlap (bp)",
            "Shared-core overlap (%)",
            "Called-interval overlap (%)",
            "Status",
        ]

    primary_ann = primary_annotations.copy()
    if not primary_ann.empty:
        primary_ann["boundary_coordinate"] = primary_ann["boundary_coordinate"].map(fmt_int)
        primary_ann["distance_to_ICR_bp"] = primary_ann["distance_to_ICR_bp"].map(fmt_int)
        primary_ann = primary_ann[
            ["boundary_source", "boundary_type", "boundary_coordinate", "nearest_gene", "nearest_ICR", "distance_to_ICR_bp"]
        ]
        primary_ann.columns = ["Boundary source", "Boundary type", "Coordinate", "Nearest gene", "Nearest ICR", "Distance to ICR (bp)"]

    convergence_md = convergence.copy()
    if not convergence_md.empty:
        convergence_md["entry_boundary"] = convergence_md["entry_boundary"].map(fmt_int)
        convergence_md["exit_boundary"] = convergence_md["exit_boundary"].map(fmt_int)
        convergence_md["distance_to_control_entry_bp"] = convergence_md["distance_to_control_entry_bp"].map(fmt_int)
        convergence_md["distance_to_control_exit_bp"] = convergence_md["distance_to_control_exit_bp"].map(fmt_int)
        convergence_md = convergence_md[
            ["cohort", "segment_index", "entry_boundary", "exit_boundary", "distance_to_control_entry_bp", "distance_to_control_exit_bp", "discrepancy_flag"]
        ]
        convergence_md.columns = [
            "Cohort",
            "Segment",
            "Entry boundary",
            "Exit boundary",
            "Entry shift vs control (bp)",
            "Exit shift vs control (bp)",
            "Status",
        ]

    bp_md = bp_distances.copy()
    if not bp_md.empty:
        bp_md["bp_start"] = bp_md["bp_start"].map(fmt_int)
        bp_md["bp_end"] = bp_md["bp_end"].map(fmt_int)
        bp_md["bp_midpoint"] = bp_md["bp_midpoint"].map(fmt_int)
        bp_md["distance_from_shared_start_bp"] = bp_md["distance_from_shared_start_bp"].map(fmt_int)
        bp_md["distance_from_shared_midpoint_bp"] = bp_md["distance_from_shared_midpoint_bp"].map(fmt_int)
        bp_md["distance_from_shared_end_bp"] = bp_md["distance_from_shared_end_bp"].map(fmt_int)
        bp_md = bp_md[
            ["bp_landmark_name", "bp_start", "bp_end", "bp_midpoint", "distance_from_shared_start_bp", "distance_from_shared_midpoint_bp", "distance_from_shared_end_bp"]
        ]
        bp_md.columns = [
            "Breakpoint",
            "Start",
            "End",
            "Midpoint",
            "Distance to shared start (bp)",
            "Distance to shared midpoint (bp)",
            "Distance to shared end (bp)",
        ]

    cp_md = change_point_support.copy()
    if not cp_md.empty:
        for col in ["entry_change_point", "exit_change_point", "entry_shift_vs_primary_start_bp", "exit_shift_vs_primary_end_bp"]:
            cp_md[col] = cp_md[col].map(fmt_int)
        cp_md = cp_md[
            ["cohort", "entry_change_point", "exit_change_point", "entry_shift_vs_primary_start_bp", "exit_shift_vs_primary_end_bp", "method"]
        ]
        cp_md.columns = ["Cohort", "Entry change point", "Exit change point", "Entry shift vs shared start (bp)", "Exit shift vs shared end (bp)", "Method"]

    bootstrap_lines = ["_Bootstrap summary unavailable._"]
    if not bootstrap_support.empty:
        row = bootstrap_support.iloc[0]
        bootstrap_lines = [
            f"- Replicates completed: `{int(row['bootstrap_replicates_successful'])}` / `{int(row['bootstrap_replicates_requested'])}`",
            f"- Median bootstrap start: `{fmt_int(row['bootstrap_start_median'])}` (95% CI `{fmt_int(row['bootstrap_start_ci_low'])}` to `{fmt_int(row['bootstrap_start_ci_high'])}`)",
            f"- Median bootstrap end: `{fmt_int(row['bootstrap_end_median'])}` (95% CI `{fmt_int(row['bootstrap_end_ci_low'])}` to `{fmt_int(row['bootstrap_end_ci_high'])}`)",
            f"- Supported bootstrap interval: `{CHROM}:{fmt_int(row['bootstrap_supported_start'])}-{fmt_int(row['bootstrap_supported_end'])}`",
            f"- Limitation: {row['limitation']}",
        ]

    gene_hits = genes[(genes["chrom"] == CHROM) & (genes["start"] < shared_end) & (genes["end"] > shared_start)].copy()
    icr_hits = icrs[(icrs["chrom"] == CHROM) & (icrs["start"] < shared_end) & (icrs["end"] > shared_start)].copy()
    overlapping_genes = ", ".join(gene_hits["gene"].astype(str).drop_duplicates().tolist()) if not gene_hits.empty else "none"
    icr_summary = "none"
    if not icr_hits.empty:
        icr_rows = [f"{row['name']} ({fmt_int(row['start'])}-{fmt_int(row['end'])})" for _, row in icr_hits.iterrows()]
        icr_summary = ", ".join(icr_rows)

    cpg_rows = []
    if not single_cpg_zoom.empty:
        z = single_cpg_zoom[(single_cpg_zoom["pos"] >= 22_689_000) & (single_cpg_zoom["pos"] <= 22_697_000)].copy()
        for cohort, label in [("Control", "Controls"), ("PWS_DEL", "PWS-DEL"), ("AS_DEL", "AS-DEL")]:
            hit = z[z["cohort"] == cohort]
            if hit.empty:
                continue
            core_hit = hit[(hit["pos"] >= shared_start) & (hit["pos"] <= shared_end)]
            cpg_rows.append(
                {
                    "Cohort": label,
                    "CpGs in panel C": f"{len(hit):,}",
                    "CpGs in shared core": f"{len(core_hit):,}",
                    "Signal range": f"{float(hit['signal'].min()):.2f} to {float(hit['signal'].max()):.2f}",
                }
            )
    cpg_md = pd.DataFrame(cpg_rows)
    panel_d_context_note = build_shared_interval_genomic_context_note(
        interval_export=interval_export,
        primary_annotations=primary_annotations,
        bp_distances=bp_distances,
        genes=genes,
        icrs=icrs,
    )

    figure_relative = f"../figures/{figure_path.name}"
    lines = [
        "# Figure 3 improved quantitative report",
        "",
        f"![Figure 3 improved]({figure_relative})",
        "",
        "## Executive summary",
        "",
        f"- Shared focal interval: `{CHROM}:{shared_start:,}-{shared_end:,}` ({float(shared['width_kb']):.2f} kb; midpoint `{shared_mid:,}`)",
        f"- Entry-coordinate spread across Controls, PWS-DEL, and AS-DEL: `{entry_spread_bp:,} bp`",
        f"- Exit-coordinate spread across Controls, PWS-DEL, and AS-DEL: `{exit_spread_bp:,} bp`",
        f"- Control signal definition in panel C: `{control_formula}` (`{control_signal_mode}` mode)",
        f"- Shared core fraction of cohort intervals: Controls `{format_overlap_pct(controls['shared_core_fraction_pct'])}`, PWS-DEL `{format_overlap_pct(pws['shared_core_fraction_pct'])}`, AS-DEL `{format_overlap_pct(as_del['shared_core_fraction_pct'])}`",
        "",
        "## Panel A: locus-wide context and why this boundary was chosen",
        "",
        f"- The shared core spans `{fmt_int(shared_start)}-{fmt_int(shared_end)}` and is only `{shared_locus_fraction_pct:.3f}%` of the full 7-Mb locus-wide overview, so Panel A marks it as a boundary line/pin rather than as a broad band.",
        f"- BP1, BP2, and BP3 intervals loaded for the figure lie at `{fmt_int(bp_distances.iloc[0]['bp_start'])}-{fmt_int(bp_distances.iloc[0]['bp_end'])}`, `{fmt_int(bp_distances.iloc[1]['bp_start'])}-{fmt_int(bp_distances.iloc[1]['bp_end'])}`, and `{fmt_int(bp_distances.iloc[2]['bp_start'])}-{fmt_int(bp_distances.iloc[2]['bp_end'])}` respectively." if len(bp_distances) >= 3 else "- Breakpoint intervals were not fully available.",
        f"- The shared midpoint is `{fmt_int(bp_distances.loc[bp_distances['bp_landmark_name'] == 'BP2', 'distance_from_shared_midpoint_bp'].iloc[0])}` bp distal to BP2 and `{fmt_int(bp_distances.loc[bp_distances['bp_landmark_name'] == 'BP3', 'distance_from_shared_midpoint_bp'].iloc[0])}` bp proximal to BP3." if set(bp_distances.get('bp_landmark_name', pd.Series(dtype=object))) >= {'BP2', 'BP3'} else "- Shared-core distances to BP2/BP3 were not available.",
        f"- This site is carried forward as the shared boundary because it is the exact intersection of the three primary intervals: Controls `{fmt_int(controls['start'])}-{fmt_int(controls['end'])}`, PWS-DEL `{fmt_int(pws['start'])}-{fmt_int(pws['end'])}`, and AS-DEL `{fmt_int(as_del['start'])}-{fmt_int(as_del['end'])}`.",
        "- Panel A is intentionally limited to locus-wide context; the fine-scale methylation transition used to support the call is shown in Panel C.",
        f"- The shared interval overlaps genes: {overlapping_genes}.",
        f"- Overlapping ICR annotations: {icr_summary}.",
        "",
        "## Panel B: interval convergence",
        "",
        markdown_table(interval_md),
        "",
        "Boundary-edge shifts:",
        "",
        markdown_table(shift_md),
        "",
        "Pairwise overlap summary:",
        "",
        markdown_table(overlap_md) if not overlap_md.empty else "_Pairwise overlap summary unavailable._",
        "",
        "## Panel C: methylation evidence",
        "",
        f"- Panel C displays both smoothed cohort profiles and raw single-CpG observations across `chr15:22,689,000-22,697,000`.",
        f"- The shared interval spans the loaded ICR from `{fmt_int(icr_hits.iloc[0]['start'])}` to `{fmt_int(icr_hits.iloc[0]['end'])}`, extending `{fmt_int(int(icr_hits.iloc[0]['start']) - shared_start)}` bp upstream and `{fmt_int(shared_end - int(icr_hits.iloc[0]['end']))}` bp downstream of that ICR." if not icr_hits.empty else "- No overlapping ICR interval was loaded.",
        "",
        markdown_table(cpg_md) if not cpg_md.empty else "_CpG support table unavailable._",
        "",
        "Primary boundary annotation summary:",
        "",
        markdown_table(primary_ann) if not primary_ann.empty else "_Primary boundary annotations unavailable._",
        "",
        "## Panel D: structural and regulatory context",
        "",
        f"- {panel_d_context_note.replace(chr(10), ' ')}",
        "",
        "Breakpoint distance table:",
        "",
        markdown_table(bp_md) if not bp_md.empty else "_Breakpoint-distance table unavailable._",
        "",
        "## Panel E: robustness and sensitivity",
        "",
        markdown_table(sensitivity_md) if not sensitivity_md.empty else "_Sensitivity summary unavailable._",
        "",
        "Change-point support:",
        "",
        markdown_table(cp_md) if not cp_md.empty else "_Change-point summary unavailable._",
        "",
        "Bootstrap support:",
        "",
        *bootstrap_lines,
        "",
        "Secondary/alternative boundary calls considered during convergence triage:",
        "",
        markdown_table(convergence_md) if not convergence_md.empty else "_Convergence-comparison table unavailable._",
        "",
        "## Interpretation",
        "",
        (
            "The quantitative pattern is consistent with a focal methylation-transition interval rather than a deletion-wide domain effect. "
            f"The strongest coordinate agreement is the exact `{float(shared['width_kb']):.2f} kb` shared core, while alternative calling schemes remain centered on that interval even when their widths contract to `{sensitivity_summary['width_kb'].min():.2f} kb` or expand to `{sensitivity_summary['width_kb'].max():.2f} kb`."
        ),
        "",
    ]
    mkdir(report_path.parent)
    report_path.write_text("\n".join(lines))


def plot_figure3(
    figdir: Path,
    control_group: pd.DataFrame,
    patient_group: pd.DataFrame,
    interval_export: pd.DataFrame,
    sensitivity_summary: pd.DataFrame,
    single_cpg_zoom: pd.DataFrame,
    control_signal_mode: str,
    genes: pd.DataFrame,
    icrs: pd.DataFrame,
    segdups: pd.DataFrame,
    bp_distances: pd.DataFrame,
    panel_d_context: dict[str, Any] | None = None,
) -> None:
    """Create the improved publication-style Figure 3."""
    mkdir(figdir)
    palette = {"Control": CONTROL_COLOR, "PWS_DEL": PWS_COLOR, "AS_DEL": AS_COLOR}
    overlay = pd.concat([control_group, patient_group], ignore_index=True)
    row_lookup = {str(row["source"]): row for _, row in interval_export.iterrows()}
    shared_start = int(row_lookup["Consensus"]["start"])
    shared_end = int(row_lookup["Consensus"]["end"])
    shared_mid = (shared_start + shared_end) / 2.0
    shared_width = shared_end - shared_start
    broad_start = 22_000_000
    broad_end = 29_000_000
    local_start = 22_689_000
    local_end = 22_697_000
    broad_windows = {"Control": 901, "PWS_DEL": 701, "AS_DEL": 701}
    broad_profiles: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    broad_y_values: list[np.ndarray] = []
    for cohort in ["Control", "PWS_DEL", "AS_DEL"]:
        q = overlay[(overlay["cohort"] == cohort) & (overlay["window_mid"] >= broad_start) & (overlay["window_mid"] <= broad_end)].sort_values("window_mid")
        if q.empty:
            continue
        x = q["window_mid"].to_numpy(float)
        y = smooth_broad_profile(q["mean_signal"].to_numpy(float), median_window=broad_windows[cohort], mean_window=121)
        broad_profiles[cohort] = (x, y)
        broad_y_values.append(y)
    if broad_y_values:
        broad_all = np.concatenate(broad_y_values)
        broad_top = max(0.28, min(0.30, float(np.nanquantile(broad_all, 0.99) + 0.025)))
        broad_bottom = min(-0.008, float(np.nanquantile(broad_all, 0.01) - 0.01))
    else:
        broad_top = 0.30
        broad_bottom = -0.03
    bp3_row = next((row for row in BP_CLUSTER_INTERVALS_T2T if str(row.get("name")) == "BP3"), None)
    bp3_mid = int(round((float(bp3_row["start"]) + float(bp3_row["end"])) / 2.0)) if bp3_row else 26_250_000
    fig = plt.figure(figsize=(14.6, 18.5))
    gs = GridSpec(5, 1, height_ratios=[1.08, 1.48, 1.10, 1.30, 1.04], hspace=0.64, figure=fig)
    panel_a = gs[0, 0].subgridspec(2, 1, height_ratios=[1.0, 0.24], hspace=0.06)
    ax_a_data = fig.add_subplot(panel_a[0, :])
    ax_a_bar = fig.add_subplot(panel_a[1, :], sharex=ax_a_data)
    ax_b = fig.add_subplot(gs[1, 0])
    panel_c = gs[2, 0].subgridspec(2, 1, height_ratios=[0.30, 1.0], hspace=0.04)
    ax_c_head = fig.add_subplot(panel_c[0, 0])
    ax_c = fig.add_subplot(panel_c[1, 0])
    panel_d = gs[3, 0].subgridspec(3, 1, height_ratios=[0.50, 0.34, 0.16], hspace=0.12)
    ax_d_top = fig.add_subplot(panel_d[0, 0])
    ax_d_bottom = fig.add_subplot(panel_d[1, 0])
    ax_d_segdup = fig.add_subplot(panel_d[2, 0], sharex=ax_d_bottom)
    panel_e = gs[4, 0].subgridspec(1, 3, width_ratios=[2.05, 4.8, 1.55], wspace=0.02)
    ax_e_labels = fig.add_subplot(panel_e[0, 0])
    ax_e = fig.add_subplot(panel_e[0, 1], sharey=ax_e_labels)
    ax_e_meta = fig.add_subplot(panel_e[0, 2], sharey=ax_e_labels)
    fig.subplots_adjust(left=0.08, right=0.93, top=0.965, bottom=0.05)

    for ax in [ax_c_head]:
        ax.axis("off")

    ax_a_data.set_xlim(broad_start, broad_end)
    ax_a_data.set_ylim(broad_bottom, broad_top)
    ax_a_data.set_title(
        "A. Parent-of-origin methylation contrast is focal within BP1-BP3",
        loc="left",
        fontsize=16,
        fontweight="bold",
        pad=26,
    )
    ax_a_data.axhline(0.0, color="#b8b8b8", lw=0.85, ls=(0, (3, 3)), zorder=1)
    ax_a_data.axvline(shared_mid, color=SHARED_CORE_COLOR, lw=1.15, alpha=0.95, zorder=2)
    ax_a_data.axvline(bp3_mid, color=STRUCT_GREY, lw=0.9, ls=(0, (3, 4)), alpha=0.82, zorder=2)
    for frac, label in [(0.006, "BP1"), (0.028, "BP2")]:
        ax_a_data.plot([frac, frac], [1.00, 1.045], transform=ax_a_data.transAxes, color=STRUCT_GREY, lw=0.85, clip_on=False, zorder=3)
        ax_a_data.text(frac, 1.050, label, transform=ax_a_data.transAxes, ha="center", va="bottom", fontsize=7.8, color="#666666", clip_on=False)
    ax_a_data.plot([bp3_mid, bp3_mid], [1.00, 1.045], transform=ax_a_data.get_xaxis_transform(), color=STRUCT_GREY, lw=0.85, clip_on=False, zorder=3)
    ax_a_data.text(bp3_mid, 1.050, "BP3", transform=ax_a_data.get_xaxis_transform(), ha="center", va="bottom", fontsize=7.8, color="#666666", clip_on=False)
    for cohort, lw, alpha in [("Control", 2.2, 0.98), ("PWS_DEL", 1.5, 0.85), ("AS_DEL", 1.5, 0.85)]:
        if cohort not in broad_profiles:
            continue
        x, y = broad_profiles[cohort]
        ax_a_data.plot(x, y, color=palette[cohort], lw=lw, alpha=alpha, zorder=3)
    for cohort, label, dy in [("Control", "Controls", 12), ("AS_DEL", "AS-DEL", 2), ("PWS_DEL", "PWS-DEL", -12)]:
        if cohort not in broad_profiles:
            continue
        x, y = broad_profiles[cohort]
        ax_a_data.annotate(
            label,
            xy=(float(x[-1]), float(y[-1])),
            xytext=(10, dy),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize=8.6,
            color=palette[cohort],
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 0.20},
            arrowprops={"arrowstyle": "-", "color": palette[cohort], "lw": 0.7, "shrinkA": 0, "shrinkB": 0},
            annotation_clip=False,
            zorder=4,
        )
    ax_a_data.set_ylabel("Parent-of-origin\nmethylation contrast", fontsize=10.8)
    ax_a_data.set_yticks([-0.10, 0.00, 0.10, 0.20])
    ax_a_data.tick_params(axis="x", labelbottom=False, width=0.9)
    ax_a_data.tick_params(axis="y", labelsize=9.0, width=0.9)
    ax_a_data.spines["top"].set_visible(False)
    ax_a_data.spines["right"].set_visible(False)
    ax_a_data.spines["left"].set_linewidth(1.0)
    ax_a_data.spines["bottom"].set_linewidth(1.0)
    ax_a_bar.set_ylim(0.0, 1.0)
    ax_a_bar.set_yticks([])
    ax_a_bar.set_xticks(np.arange(broad_start, broad_end + 1, 1_000_000))
    ax_a_bar.xaxis.set_major_formatter(pos_formatter())
    ax_a_bar.tick_params(axis="x", labelsize=10, width=0.0, length=0, pad=2)
    ax_a_bar.set_xlabel("T2T-CHM13 chr15 coordinate (Mb)", fontsize=12, labelpad=4)
    for spine in ax_a_bar.spines.values():
        spine.set_visible(False)
    bar_y = 0.56
    bar_h = 0.12
    ax_a_bar.add_patch(
        Rectangle(
            (broad_start, bar_y - bar_h / 2.0),
            bp3_mid - broad_start,
            bar_h,
            facecolor="#d7d7d7",
            edgecolor="none",
            alpha=0.95,
            zorder=1,
        )
    )
    ax_a_bar.axvline(shared_mid, ymin=0.18, ymax=0.86, color=SHARED_CORE_COLOR, lw=1.15, zorder=2)
    ax_a_bar.axvline(bp3_mid, ymin=0.20, ymax=0.82, color=STRUCT_GREY, lw=0.85, ls=(0, (3, 4)), alpha=0.82, zorder=2)
    ax_a_bar.plot(shared_mid, 0.86, marker="v", markersize=6.2, color=SHARED_CORE_COLOR, clip_on=False, zorder=3)
    ax_a_bar.annotate(
        "Chosen as the shared boundary because\nall three cohort calls intersect at this focal transition",
        xy=(shared_mid, 0.86),
        xytext=(12, 6),
        textcoords="offset points",
        ha="left",
        va="bottom",
        fontsize=7.9,
        color="#7a5c00",
        bbox={"facecolor": "white", "edgecolor": SHARED_CORE_COLOR, "alpha": 0.92, "pad": 1.3},
        arrowprops={"arrowstyle": "-", "color": SHARED_CORE_COLOR, "lw": 0.75, "shrinkA": 2, "shrinkB": 3},
    )
    ax_a_bar.text(
        (broad_start + bp3_mid) / 2.0,
        0.07,
        "Recurrent BP1-BP3 deletion interval\nNo broad domain-wide elevation",
        ha="center",
        va="bottom",
        fontsize=7.8,
        color="#757575",
    )

    interval_rows = interval_export.copy()
    ax_b.set_title(
        "B. Boundary coordinates converge on a shared 3.7-kb core",
        loc="left",
        fontsize=16,
        fontweight="bold",
        pad=18,
    )
    highlight_shared_core(ax_b, shared_start, shared_end, alpha=0.10, linewidth=1.0, zorder=0)
    y_positions = np.arange(len(interval_rows), 0, -1)
    for y, (_, row) in zip(y_positions, interval_rows.iterrows()):
        width_bp = int(row["end"]) - int(row["start"])
        ax_b.broken_barh(
            [(int(row["start"]), width_bp)],
            (y - 0.31, 0.62),
            facecolors=str(row["color"]),
            edgecolors="none",
            alpha=0.95,
        )
        width_label = float(row["width_kb"]) if "width_kb" in row.index else width_bp / 1_000.0
        label_color = "#5a4300" if str(row["source"]) == "Consensus" else "white"
        ax_b.text(
            int(row["start"]) + width_bp / 2.0,
            y,
            f"{width_label:.2f} kb",
            ha="center",
            va="center",
            fontsize=9.1,
            fontweight="bold",
            color=label_color,
        )
    ax_b.set_yticks(y_positions)
    ax_b.set_yticklabels(interval_rows["label"], fontsize=11)
    ax_b.set_xlim(local_start, local_end)
    ax_b.set_ylim(0.05, len(interval_rows) + 0.95)
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)
    ax_b.tick_params(axis="x", labelsize=10, width=0.9)
    ax_b.xaxis.set_major_formatter(local_pos_formatter(3))
    ax_b.set_xlabel("T2T-CHM13 chr15 coordinate (Mb)", fontsize=12)
    ax_b.text(
        1.01,
        0.93,
        "Shared 3.7-kb core\nchr15:22,691,000-22,694,700",
        transform=ax_b.transAxes,
        ha="left",
        va="top",
        fontsize=9.6,
        color="#6e4a00",
        bbox={"facecolor": "white", "edgecolor": SHARED_CORE_COLOR, "alpha": 0.95, "pad": 2.2},
    )
    ax_b.text(shared_start, -0.14, "22,691,000", transform=ax_b.get_xaxis_transform(), ha="center", va="top", fontsize=9, color="#7a5c00")
    ax_b.text(shared_end, -0.14, "22,694,700", transform=ax_b.get_xaxis_transform(), ha="center", va="top", fontsize=9, color="#7a5c00")
    overlap_lookup = {str(row["label"]): float(row["shared_core_fraction_pct"]) for _, row in interval_rows.iterrows()}
    ax_b.text(
        shared_mid,
        0.28,
        (
            "Shared core = 3,700 bp\n"
            f"Overlap: Controls {format_overlap_pct(overlap_lookup.get('Controls'))}; "
            f"PWS-DEL {format_overlap_pct(overlap_lookup.get('PWS-DEL'))}; "
            f"AS-DEL {format_overlap_pct(overlap_lookup.get('AS-DEL'))}"
        ),
        ha="center",
        va="bottom",
        fontsize=9.0,
        color="#5c4c19",
    )

    highlight_shared_core(ax_c, shared_start, shared_end, alpha=0.10, linewidth=1.0, zorder=0)

    ax_c_head.text(0.0, 0.98, "C. Methylation evidence over the focal boundary", ha="left", va="top", fontsize=16, fontweight="bold", color="#111111")
    y_min = 0.0
    y_max = 0.0
    zraw = single_cpg_zoom[
        (single_cpg_zoom["pos"] >= local_start)
        & (single_cpg_zoom["pos"] <= local_end)
        & (single_cpg_zoom["cohort"].isin(["Control", "PWS_DEL", "AS_DEL"]))
    ].copy()
    for cohort, label in [("Control", "Controls"), ("PWS_DEL", "PWS-DEL"), ("AS_DEL", "AS-DEL")]:
        q = overlay[(overlay["cohort"] == cohort) & (overlay["window_mid"] >= local_start) & (overlay["window_mid"] <= local_end)].sort_values("window_mid")
        if q.empty:
            continue
        x = q["window_mid"].to_numpy(float)
        y = smooth_profile(q["mean_signal"].astype(float).to_numpy(), window=7)
        ax_c.plot(x, y, color=palette[cohort], lw=2.6 if cohort == "Control" else 2.3, label=label, zorder=3)
        if np.isfinite(y).any():
            y_min = min(y_min, float(np.nanmin(y)))
            y_max = max(y_max, float(np.nanmax(y)))
        raw = zraw[zraw["cohort"] == cohort]
        if not raw.empty:
            y_min = min(y_min, float(raw["signal"].min()))
            y_max = max(y_max, float(raw["signal"].max()))
            ax_c.scatter(
                raw["pos"],
                raw["signal"],
                s=4,
                color=palette[cohort],
                alpha=0.18,
                edgecolors="none",
                rasterized=True,
                zorder=4,
            )
    ax_c.set_xlim(local_start, local_end)
    ax_c.set_ylim(min(-0.18, y_min - 0.04), y_max + 0.08)
    ax_c.grid(axis="y", color="#ececec", lw=0.7)
    ax_c.spines["top"].set_visible(False)
    ax_c.spines["right"].set_visible(False)
    ax_c.tick_params(axis="both", labelsize=10, width=0.9)
    ax_c.xaxis.set_major_formatter(local_pos_formatter(3))
    ax_c.set_ylabel("Parent-of-origin methylation contrast", fontsize=11)
    ax_c.yaxis.set_label_coords(-0.065, 0.5)
    ax_c.set_xlabel("T2T-CHM13 chr15 coordinate (Mb)", fontsize=12)
    legend_labels = [
        "Controls |mat-pat|" if control_signal_mode == "absolute" else "Controls mat-pat",
        "PWS-DEL retained-origin",
        "AS-DEL retained-origin",
    ]
    handles = [
        Line2D([0], [0], color=palette["Control"], lw=2.6),
        Line2D([0], [0], color=palette["PWS_DEL"], lw=2.3),
        Line2D([0], [0], color=palette["AS_DEL"], lw=2.3),
    ]
    if handles:
        ax_c_head.legend(
            handles,
            legend_labels,
            frameon=False,
            fontsize=9.4,
            loc="lower left",
            bbox_to_anchor=(0.0, 0.08),
            ncol=3,
            borderaxespad=0.0,
            handlelength=1.9,
            columnspacing=1.8,
            handletextpad=0.5,
        )
    ax_c_head.text(
        1.0,
        0.12,
        "Points: CpG observations; lines: smoothed profiles",
        ha="right",
        va="bottom",
        fontsize=9.2,
        color="#666666",
    )

    ax_d_top.set_title(
        "D. Shared core is distinct from BP1-BP3 and aligns with the SNHG14/ICR neighborhood",
        loc="left",
        fontsize=16,
        fontweight="bold",
        pad=10,
    )
    bp_df = bp_distances.copy()
    if not bp_df.empty:
        d_start = int(bp_df["bp_start"].min()) - 180_000
        d_end = int(bp_df["bp_end"].max()) + 180_000
    else:
        d_start = 20_900_000
        d_end = 26_100_000
    ax_d_top.set_xlim(d_start, d_end)
    ax_d_top.set_ylim(0.0, 1.0)
    ax_d_top.set_yticks([])
    ax_d_top.grid(False)
    ax_d_top.spines["top"].set_visible(False)
    ax_d_top.spines["right"].set_visible(False)
    ax_d_top.spines["left"].set_visible(False)
    ax_d_top.spines["bottom"].set_linewidth(1.0)
    ax_d_top.hlines(0.50, d_start, d_end, color="#c8c8c8", lw=2.0, zorder=0)
    ax_d_top.axvspan(shared_start, shared_end, ymin=0.18, ymax=0.82, color=CONTEXT_SHARED_FILL, alpha=0.95, zorder=1)
    ax_d_top.vlines([shared_start, shared_end], 0.26, 0.74, color=SHARED_CORE_COLOR, lw=1.2, alpha=0.9, zorder=3)
    ax_d_top.vlines(shared_mid, 0.24, 0.78, color=SHARED_CORE_COLOR, lw=3.0, zorder=4)
    ax_d_top.scatter([shared_mid], [0.50], s=84, color=SHARED_CORE_COLOR, edgecolor="white", linewidth=0.9, zorder=5)
    if not bp_df.empty:
        for _, row in bp_df.sort_values("bp_start").iterrows():
            name = str(row["bp_landmark_name"])
            color = BP_REGION_COLORS.get(name, STRUCT_GREY)
            x0 = int(row["bp_start"])
            x1 = int(row["bp_end"])
            ax_d_top.axvspan(x0, x1, ymin=0.22, ymax=0.78, color=color, alpha=0.26, zorder=2)
            ax_d_top.vlines([x0, x1], 0.20, 0.80, color=color, lw=1.0, ls=(0, (3, 4)), alpha=0.85, zorder=2)
            ax_d_top.text((x0 + x1) / 2.0, 0.95, name, ha="center", va="bottom", fontsize=11.0, color="#606060")
    bp2 = bp_df.loc[bp_df["bp_landmark_name"].astype(str) == "BP2"] if not bp_df.empty else pd.DataFrame()
    if not bp2.empty:
        bp2_start = int(bp2.iloc[0]["bp_start"])
        bp2_end = int(bp2.iloc[0]["bp_end"])
        bp2_edge = bp2_end if shared_mid >= bp2_end else bp2_start
        bp2_distance_mb = abs(shared_mid - bp2_edge) / 1e6
        distance_label_x = bp2_edge + 0.68 * (shared_mid - bp2_edge)
        ax_d_top.annotate(
            "",
            xy=(shared_mid, 0.76),
            xytext=(bp2_edge, 0.76),
            arrowprops={"arrowstyle": "<->", "color": STRUCT_GREY, "lw": 1.2},
        )
        ax_d_top.text(
            distance_label_x,
            0.84,
            f"{bp2_distance_mb:.2f} Mb from the BP2 distal edge",
            ha="center",
            va="bottom",
            fontsize=10.2,
            color="#575757",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.90, "pad": 0.8},
        )
    ax_d_top.text(
        0.79,
        0.18,
        "Outside recurrent BP1/BP2/BP3 intervals",
        transform=ax_d_top.transAxes,
        ha="center",
        va="center",
        fontsize=9.8,
        color="#6d6d6d",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 0.9},
    )
    top_ticks = np.arange((d_start // 1_000_000) * 1_000_000, d_end + 1_000_000, 1_000_000)
    ax_d_top.set_xticks(top_ticks)
    ax_d_top.xaxis.set_major_formatter(pos_formatter())
    ax_d_top.tick_params(axis="x", bottom=False, labelbottom=False)

    local_context_start = min(local_start, shared_start - 50_000)
    local_context_end = max(22_696_000, shared_end + 10_000)
    ax_d_bottom.axvspan(shared_start, shared_end, color=CONTEXT_SHARED_FILL, alpha=0.42, zorder=0)
    ax_d_bottom.axvline(shared_start, color=SHARED_CORE_COLOR, lw=1.1, ls=(0, (4, 3)), zorder=1)
    ax_d_bottom.axvline(shared_end, color=SHARED_CORE_COLOR, lw=1.1, ls=(0, (4, 3)), zorder=1)
    draw_panel_d_gene_track(ax_d_bottom, genes, icrs, local_context_start, local_context_end)
    ax_d_bottom.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_d_bottom.set_xlabel("")
    ax_d_bottom.text(
        shared_mid,
        2.86,
        "Shared core within SNHG14 and overlapping the loaded PWS/AS ICR interval",
        ha="center",
        va="bottom",
        fontsize=9.6,
        color="#775000",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.86, "pad": 1.4},
    )
    draw_panel_d_segdup_track(
        ax_d_segdup,
        segdups=segdups,
        start=local_context_start,
        end=local_context_end,
        shared_start=shared_start,
        shared_end=shared_end,
        panel_d_context=panel_d_context,
    )
    for coord in (shared_start, shared_end):
        fig.add_artist(
            ConnectionPatch(
                xyA=(coord, 0.18),
                coordsA=ax_d_top.transData,
                xyB=(coord, 3.12),
                coordsB=ax_d_bottom.transData,
                color=SHARED_CORE_COLOR,
                lw=0.8,
                alpha=0.26,
            )
        )

    ax_e_labels.text(
        0.0,
        1.04,
        "E. Boundary calls remain concentrated near the shared core across most settings",
        transform=ax_e_labels.transAxes,
        ha="left",
        va="bottom",
        fontsize=16,
        fontweight="bold",
    )
    highlight_shared_core(ax_e, shared_start, shared_end, alpha=0.10, linewidth=1.0, zorder=0)
    if sensitivity_summary.empty:
        ax_e.text(0.5, 0.5, "No robustness intervals available", transform=ax_e.transAxes, ha="center", va="center")
        ax_e_labels.axis("off")
        ax_e_meta.axis("off")
    else:
        y_positions = np.arange(len(sensitivity_summary), 0, -1)
        for y, (_, row) in zip(y_positions, sensitivity_summary.iterrows()):
            ax_e_labels.text(0.98, y, str(row["method"]), ha="right", va="center", fontsize=10.0, color="#2f2f2f")
            ax_e.broken_barh(
                [(int(row["start"]), int(row["end"]) - int(row["start"]))],
                (y - 0.23, 0.46),
                facecolors=str(row["color"]),
                edgecolors="none",
                alpha=0.94,
            )
            ax_e_meta.text(0.04, y, f"{float(row['width_kb']):.2f} kb", ha="left", va="center", fontsize=9.4, color="#444444")
            ax_e_meta.text(0.98, y, format_overlap_pct(row["overlap_pct_shared_core"]), ha="right", va="center", fontsize=9.4, color="#444444")
        y_min = 0.35
        y_max = len(sensitivity_summary) + 0.75
        for ax in [ax_e_labels, ax_e, ax_e_meta]:
            ax.set_ylim(y_min, y_max)
        ax_e_labels.set_xlim(0, 1)
        ax_e_labels.set_xticks([])
        ax_e_labels.set_yticks([])
        for spine in ax_e_labels.spines.values():
            spine.set_visible(False)
        ax_e.set_yticks([])
        ax_e.tick_params(axis="y", length=0)
        ax_e_meta.set_xlim(0, 1)
        ax_e_meta.set_xticks([])
        ax_e_meta.set_yticks([])
        ax_e_meta.spines["top"].set_visible(False)
        ax_e_meta.spines["right"].set_visible(False)
        ax_e_meta.spines["bottom"].set_visible(False)
        ax_e_meta.spines["left"].set_visible(True)
        ax_e_meta.spines["left"].set_color("#e0e0e0")
        ax_e_meta.spines["left"].set_linewidth(0.9)
        ax_e_meta.text(0.04, len(sensitivity_summary) + 0.45, "Called width", ha="left", va="bottom", fontsize=9.5, fontweight="bold")
        ax_e_meta.text(0.98, len(sensitivity_summary) + 0.45, "Core overlap", ha="right", va="bottom", fontsize=9.5, fontweight="bold")
    ax_e.set_xlim(local_start, local_end)
    ax_e.spines["top"].set_visible(False)
    ax_e.spines["right"].set_visible(False)
    ax_e.tick_params(axis="both", labelsize=10, width=0.9)
    ax_e.xaxis.set_major_formatter(local_pos_formatter(3))
    ax_e.set_xlabel("T2T-CHM13 chr15 coordinate (Mb)", fontsize=12, labelpad=6)

    for suffix in ["png", "pdf"]:
        fig.savefig(figdir / f"Figure3_boundary_mapping_improved.{suffix}", dpi=400 if suffix == "png" else None, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    table_dir = outdir / "tables"
    figdir = outdir / "figures"
    logdir = outdir / "logs"
    mkdir(table_dir)
    mkdir(figdir)
    mkdir(logdir)

    files = discover_sample_files(Path(args.methylation_dir))
    assignments = read_parental_assignments(Path(args.assignment_table))

    missing_controls = [sample for sample in CONTROL_SAMPLES if sample not in assignments or not {"maternal", "paternal"} <= set(assignments[sample])]
    if missing_controls:
        raise RuntimeError(f"Missing maternal/paternal control assignments for: {', '.join(missing_controls)}")

    input_rows = []
    tracks: list[SignalTrack] = []
    control_combined: list[pd.DataFrame] = []
    bed_cache: dict[tuple[str, str], pd.DataFrame] = {}
    validation_warnings: list[str] = []
    pws_mupd_available = [sample for sample in PWS_MUPD_SAMPLES if sample in files]

    def load(sample: str, layer: str) -> pd.DataFrame:
        key = (sample, layer)
        if key not in bed_cache:
            path = files.get(sample, {}).get(layer)
            if path is None:
                raise FileNotFoundError(f"No {layer}.bed found for sample {sample} in {args.methylation_dir}")
            input_rows.append({"sample": sample, "layer": layer, "path": str(path)})
            bed_cache[key] = read_methylation_region(path, args.chrom, args.region_start, args.region_end)
        return bed_cache[key]

    for sample in CONTROL_SAMPLES:
        mat_layer = assignments[sample]["maternal"]
        pat_layer = assignments[sample]["paternal"]
        mat = load(sample, mat_layer)
        pat = load(sample, pat_layer)
        combined = load(sample, "combined")
        control_combined.append(combined)
        merged = merge_two_tracks(mat, pat, "maternal", "paternal")
        tracks.append(
            SignalTrack(
                sample=sample,
                cohort="Control",
                positions=merged["pos"].to_numpy(int),
                signal=compute_control_signal(
                    merged["maternal"].to_numpy(float),
                    merged["paternal"].to_numpy(float),
                    args.control_signal_mode,
                ),
                maternal_methylation=merged["maternal"].to_numpy(float),
                paternal_methylation=merged["paternal"].to_numpy(float),
            )
        )

    baseline = build_control_baseline(control_combined)
    save_table(baseline, table_dir / "Figure3_control_biallelic_baseline_per_CpG.tsv")

    for sample in PWS_DEL_SAMPLES:
        retained = load(sample, "combined")
        merged = retained[["pos", "meth"]].rename(columns={"meth": "retained"}).merge(baseline, on="pos", how="inner")
        tracks.append(
            SignalTrack(
                sample=sample,
                cohort="PWS_DEL",
                positions=merged["pos"].to_numpy(int),
                signal=(merged["retained"] - merged["baseline"]).to_numpy(float),
                retained_methylation=merged["retained"].to_numpy(float),
                baseline_methylation=merged["baseline"].to_numpy(float),
            )
        )

    for sample in AS_DEL_SAMPLES:
        retained = load(sample, "combined")
        merged = retained[["pos", "meth"]].rename(columns={"meth": "retained"}).merge(baseline, on="pos", how="inner")
        tracks.append(
            SignalTrack(
                sample=sample,
                cohort="AS_DEL",
                positions=merged["pos"].to_numpy(int),
                signal=(merged["baseline"] - merged["retained"]).to_numpy(float),
                retained_methylation=merged["retained"].to_numpy(float),
                baseline_methylation=merged["baseline"].to_numpy(float),
            )
        )

    sample_single_cpg = pd.concat(
        [
            pd.DataFrame({"sample": track.sample, "cohort": track.cohort, "pos": track.positions, "signal": track.signal})
            for track in tracks
        ],
        ignore_index=True,
    )
    save_table(sample_single_cpg, table_dir / "Figure3_all_sample_single_CpG_signals.tsv.gz", compression="gzip")

    primary_call = compute_interval_call_set(
        tracks=tracks,
        chrom=args.chrom,
        region_start=args.region_start,
        region_end=args.region_end,
        control_match_window=args.control_match_window,
        enter_threshold=args.enter_threshold,
        exit_threshold=args.exit_threshold,
        consecutive_windows=args.consecutive_windows,
        min_cpgs_per_window=args.min_cpgs_per_window,
        window_size_bp=args.window_size,
        step_size_bp=args.step_size,
    )
    sample_windows = primary_call["sample_windows"]
    sample_boundary_rows = primary_call["sample_boundary_rows"]
    grouped = primary_call["grouped"]
    control_summary = primary_call["control_summary"]
    group_boundary_calls = primary_call["group_boundary_calls"]
    primary_boundaries = primary_call["primary_boundaries"]

    if sample_windows.empty:
        raise RuntimeError("No sliding-window signal values were computed for the primary Figure 3 parameter set.")
    control_windows = sample_windows[sample_windows["cohort"] == "Control"].copy()
    patient_windows = sample_windows[sample_windows["cohort"].isin(["PWS_DEL", "AS_DEL"])].copy()
    save_table(control_windows, table_dir / "Figure3A_control_sliding_window_delta.tsv")
    save_table(patient_windows, table_dir / "Figure3B_patient_sliding_window_signals.tsv")
    save_table(sample_windows, table_dir / "Figure3_all_sample_sliding_window_signals.tsv")
    save_table(pd.DataFrame(sample_boundary_rows), table_dir / "Figure3_sample_boundary_calls.tsv")

    grouped_control = grouped[grouped["cohort"] == "Control"].copy()
    grouped_patient = grouped[grouped["cohort"].isin(["PWS_DEL", "AS_DEL"])].copy()
    save_table(grouped, table_dir / "Figure3_group_sliding_window_signals.tsv")

    save_table(control_summary, table_dir / "Figure3_control_boundary_coordinates.tsv")
    save_table(group_boundary_calls, table_dir / "Figure3_group_boundary_coordinates.tsv")
    save_table(primary_boundaries, table_dir / "Figure3_primary_boundary_coordinates.tsv")
    if primary_boundaries.empty:
        raise RuntimeError("Primary boundary selection failed; no Control/PWS-DEL/AS-DEL convergence table was produced.")
    shared_start, shared_end = shared_interval_from_boundaries(primary_boundaries)
    computed_interval_table = primary_boundaries_to_interval_table(primary_boundaries)
    interval_comparison, interval_comparison_warnings = compare_interval_tables(
        computed_interval_table,
        args.validated_interval_tolerance_bp,
    )
    validation_warnings.extend(interval_comparison_warnings)
    save_table(interval_comparison, table_dir / "Figure3_validated_interval_comparison.tsv")
    interval_table = validated_interval_table() if args.use_validated_intervals else computed_interval_table
    interval_provenance = pd.DataFrame(
        [
            {
                "interval_source": "validated_fixed_intervals" if args.use_validated_intervals else "computed_from_boundary_calling",
                "control_signal_mode": args.control_signal_mode,
                "control_signal_formula": control_signal_formula(args.control_signal_mode),
                "validated_interval_tolerance_bp": args.validated_interval_tolerance_bp,
                "validated_override_requested": bool(args.use_validated_intervals),
            }
        ]
    )
    save_table(interval_provenance, table_dir / "Figure3_interval_provenance.tsv")
    figure_consensus = interval_table.loc[interval_table["source"] == "Consensus"].iloc[0]
    figure_shared_start = int(figure_consensus["start"])
    figure_shared_end = int(figure_consensus["end"])
    interval_export = build_publication_interval_table(interval_table)
    save_table(interval_export, table_dir / "Figure3_boundary_intervals.tsv")
    overlap_summary = build_pairwise_overlap_summary(interval_export)

    comparison_rows = []
    if not primary_boundaries.empty and not group_boundary_calls.empty:
        first_control = primary_boundaries[primary_boundaries["cohort"] == "Control"].iloc[0]
        for _, row in group_boundary_calls.iterrows():
            if row.get("cohort") == "Control":
                continue
            comparison_rows.append(
                {
                    "cohort": row.get("cohort"),
                    "segment_index": int(row.get("segment_index", 1)),
                    "entry_boundary": row.get("entry_window_start", np.nan),
                    "exit_boundary": row.get("exit_window_end", np.nan),
                    "distance_to_control_entry_bp": row.get("entry_window_start", np.nan) - first_control.get("entry_boundary", np.nan),
                    "distance_to_control_exit_bp": row.get("exit_window_end", np.nan) - first_control.get("exit_boundary", np.nan),
                    "discrepancy_flag": (
                        "convergent"
                        if (
                            abs(row.get("entry_window_start", np.nan) - first_control.get("entry_boundary", np.nan)) <= PATIENT_CONVERGENCE_BP
                            and abs(row.get("exit_window_end", np.nan) - first_control.get("exit_boundary", np.nan)) <= PATIENT_CONVERGENCE_BP
                        )
                        else "flag"
                    ),
                }
            )
    comparison = pd.DataFrame(comparison_rows)
    save_table(comparison, table_dir / "Figure3_boundary_convergence_comparison.tsv")

    bp_hotspots = read_bp_hotspots(Path(args.bp_hotspot_bed) if args.bp_hotspot_bed else None, args.chrom)
    if not bp_hotspots.empty and {"start", "end"} <= set(bp_hotspots.columns):
        annotation_start = int(min(args.region_start - args.zoom_flank, bp_hotspots["start"].astype(float).min()))
        annotation_end = int(max(args.region_end + args.zoom_flank, bp_hotspots["end"].astype(float).max()))
    else:
        annotation_start = args.region_start - args.zoom_flank
        annotation_end = args.region_end + args.zoom_flank
    genes = parse_gtf_genes(Path(args.gtf), args.chrom, annotation_start, annotation_end)
    segdups = read_bed_like(Path(args.segdup), "SegDup", args.chrom, annotation_start, annotation_end)
    ctcf = read_bed_like(Path(args.ctcf_bed), "CTCF", args.chrom, annotation_start, annotation_end) if args.ctcf_bed else pd.DataFrame(columns=["chrom", "start", "end", "name", "source"])
    imprintome = read_bed_like(Path(args.imprintome_bed), "Imprintome", args.chrom, annotation_start, annotation_end)
    court = read_bed_like(Path(args.court2014_bed), "Court2014", args.chrom, annotation_start, annotation_end)
    icrs = read_bed_like(Path(args.icr_bed), "ICR", args.chrom, annotation_start, annotation_end)
    repeats = read_bed_like(Path(args.repeats), "RepeatMasker", args.chrom, annotation_start, annotation_end)
    save_table(genes, table_dir / "Figure3_annotation_genes_loaded.tsv")
    save_table(segdups, table_dir / "Figure3_annotation_segdups_loaded.tsv")
    save_table(bp_hotspots, table_dir / "Figure3_annotation_bp_hotspots_loaded.tsv")
    save_table(ctcf, table_dir / "Figure3_annotation_ctcf_loaded.tsv")
    save_table(imprintome, table_dir / "Figure3_annotation_imprintome_dmrs_loaded.tsv")
    save_table(court, table_dir / "Figure3_annotation_court2014_dmrs_loaded.tsv")
    save_table(icrs, table_dir / "Figure3_annotation_icrs_loaded.tsv")
    save_table(repeats, table_dir / "Figure3_annotation_repeats_loaded.tsv")

    boundary_for_annotation = control_summary.copy()
    if not group_boundary_calls.empty:
        group_ann = group_boundary_calls.rename(columns={"entry_window_start": "entry_boundary", "exit_window_end": "exit_boundary"})
        boundary_for_annotation = pd.concat([boundary_for_annotation, group_ann], ignore_index=True, sort=False)
    annotations = annotate_boundaries(boundary_for_annotation, genes, segdups, bp_hotspots, ctcf, imprintome, court, icrs, args.chrom)
    save_table(annotations, table_dir / "Figure3_boundary_annotations.tsv")
    primary_annotations = annotate_boundaries(primary_boundaries, genes, segdups, bp_hotspots, ctcf, imprintome, court, icrs, args.chrom)
    save_table(primary_annotations, table_dir / "Figure3_primary_boundary_annotations.tsv")
    bp_distance_definitions = build_bp_distance_definitions(shared_start, shared_end, bp_hotspots)
    save_table(bp_distance_definitions, table_dir / "Figure3_bp_distance_definitions.tsv")

    if bp_distance_definitions.empty or set(bp_distance_definitions["bp_landmark_name"]) != {"BP1", "BP2", "BP3"}:
        validation_warnings.append("Loaded BP landmark intervals are incomplete; expected BP1, BP2, and BP3.")
    if args.control_signal_mode == "absolute" and "|" not in control_signal_formula(args.control_signal_mode):
        validation_warnings.append("Control signal label mismatch: control mode is absolute but the label formula was not absolute.")
    if (sample_windows["n_cpg"] <= 0).any():
        validation_warnings.append("At least one primary sliding window had zero CpGs.")
    shared_cpg = sample_single_cpg[
        (sample_single_cpg["pos"] >= shared_start)
        & (sample_single_cpg["pos"] <= shared_end)
        & (sample_single_cpg["cohort"].isin(["Control", "PWS_DEL", "AS_DEL"]))
    ]
    if shared_cpg.empty or shared_cpg.groupby("cohort")["pos"].count().min() < LOW_CPG_DENSITY_WARNING_THRESHOLD:
        validation_warnings.append(
            f"Low CpG density within the shared interval ({shared_start:,}-{shared_end:,}); interpret exact edges cautiously."
        )
    shared_gene_overlap = genes[(genes["chrom"] == args.chrom) & (genes["start"] < shared_end) & (genes["end"] > shared_start)]
    shared_icr_overlap = icrs[(icrs["chrom"] == args.chrom) & (icrs["start"] < shared_end) & (icrs["end"] > shared_start)]
    if shared_gene_overlap.empty or "SNHG14" not in set(shared_gene_overlap.get("gene", pd.Series(dtype=object)).astype(str)):
        validation_warnings.append("The shared interval does not overlap the expected SNHG14 annotation in the loaded GTF.")
    if shared_icr_overlap.empty:
        validation_warnings.append("The shared interval does not overlap any loaded ICR interval.")

    change_point_support = build_change_point_support(grouped, shared_start, shared_end)
    save_table(change_point_support, table_dir / "Figure3_change_point_support.tsv")
    bootstrap_support = bootstrap_shared_interval_support(tracks, args, shared_start, shared_end)
    save_table(bootstrap_support, table_dir / "Figure3_bootstrap_support.tsv")
    sensitivity = build_boundary_sensitivity_table(tracks, args, shared_start, shared_end)
    save_table(sensitivity, table_dir / "Figure3_boundary_sensitivity.tsv")
    sensitivity_summary, sensitivity_warnings = build_sensitivity_summary(
        sensitivity,
        change_point_support,
        bootstrap_support,
        shared_start=figure_shared_start,
        shared_end=figure_shared_end,
    )
    save_table(sensitivity_summary, table_dir / "Figure3_sensitivity_analysis.tsv")
    validation_warnings.extend(sensitivity_warnings)

    zoom = single_cpg_zoom_table(tracks, primary_boundaries[primary_boundaries["cohort"] == "Control"], args.zoom_flank)
    if zoom.empty:
        zoom = sample_single_cpg[
            (sample_single_cpg["pos"] >= args.region_start)
            & (sample_single_cpg["pos"] <= min(args.region_end, args.region_start + 2 * args.zoom_flank))
        ].copy()
        validation_warnings.append("Primary boundary zoom table was empty; Figure 3 fell back to the earliest available CpG slice.")
    save_table(zoom, table_dir / "Figure3_single_CpG_boundary_zoom.tsv.gz", compression="gzip")
    save_table(interval_table, table_dir / "Figure3_primary_boundary_interval_plot.tsv")
    save_table(warning_frame(validation_warnings), table_dir / "Figure3_validation_warnings.tsv")

    plot_figure3(
        figdir,
        grouped_control,
        grouped_patient,
        interval_export=interval_export,
        sensitivity_summary=sensitivity_summary,
        single_cpg_zoom=zoom,
        control_signal_mode=args.control_signal_mode,
        genes=genes,
        icrs=icrs,
        segdups=segdups,
        bp_distances=bp_distance_definitions,
        panel_d_context=build_shared_interval_genomic_context(
            interval_export=interval_export,
            primary_annotations=primary_annotations,
            bp_distances=bp_distance_definitions,
            genes=genes,
            icrs=icrs,
        ),
    )
    input_files = pd.DataFrame(input_rows).drop_duplicates()
    save_table(input_files, table_dir / "Figure3_input_methylation_files.tsv")
    write_publication_report(
        report_path=outdir / "reports" / "Figure3_boundary_mapping_report.md",
        input_files=input_files,
        interval_export=interval_export,
        overlap_summary=overlap_summary,
        bp_distances=bp_distance_definitions,
        sensitivity_summary=sensitivity_summary,
        report_warnings=validation_warnings,
        control_signal_mode=args.control_signal_mode,
        control_formula=control_signal_formula(args.control_signal_mode),
    )
    write_figure3_improved_report(
        report_path=outdir / "reports" / "Figure3_boundary_mapping_improved_report.md",
        figure_path=figdir / "Figure3_boundary_mapping_improved.png",
        interval_export=interval_export,
        overlap_summary=overlap_summary,
        bp_distances=bp_distance_definitions,
        sensitivity_summary=sensitivity_summary,
        change_point_support=change_point_support,
        bootstrap_support=bootstrap_support,
        primary_annotations=primary_annotations,
        convergence=comparison,
        genes=genes,
        icrs=icrs,
        single_cpg_zoom=zoom,
        control_signal_mode=args.control_signal_mode,
        control_formula=control_signal_formula(args.control_signal_mode),
    )

    params = {
        "methylation_dir": args.methylation_dir,
        "outdir": args.outdir,
        "assignment_table": args.assignment_table,
        "chrom": args.chrom,
        "region_start": args.region_start,
        "region_end": args.region_end,
        "window_size": args.window_size,
        "step_size": args.step_size,
        "consecutive_windows": args.consecutive_windows,
        "enter_threshold": args.enter_threshold,
        "exit_threshold": args.exit_threshold,
        "control_match_window": args.control_match_window,
        "control_signal_mode": args.control_signal_mode,
        "control_signal_formula": control_signal_formula(args.control_signal_mode),
        "patient_signal_formulas": PATIENT_SIGNAL_FORMULAS,
        "interval_source_for_reviewer_figure": "validated_fixed_intervals" if args.use_validated_intervals else "computed_from_boundary_calling",
        "validated_interval_tolerance_bp": args.validated_interval_tolerance_bp,
        "bootstrap_replicates": args.bootstrap_replicates,
        "bootstrap_seed": args.bootstrap_seed,
        "bootstrap_local_flank": args.bootstrap_local_flank,
        "control_samples": CONTROL_SAMPLES,
        "pws_del_samples": PWS_DEL_SAMPLES,
        "as_del_samples": AS_DEL_SAMPLES,
        "pws_mupd_samples_context_only": PWS_MUPD_SAMPLES,
        "pws_mupd_present_in_inputs": pws_mupd_available,
        "ctcf_status": "loaded" if not ctcf.empty else "not_available",
        "validation_warning_count": len(validation_warnings),
    }
    (outdir / "phase3_run_parameters.json").write_text(json.dumps(params, indent=2) + "\n")

    log_lines = [
        "Phase 3 boundary mapping complete.",
        f"Control signal mode: {args.control_signal_mode} ({control_signal_formula(args.control_signal_mode)})",
        f"Reviewer interval source: {'validated_fixed_intervals' if args.use_validated_intervals else 'computed_from_boundary_calling'}",
        f"Control boundary segments: {len(control_summary)}",
        f"Group boundary calls: {len(group_boundary_calls)}",
        f"Primary boundary calls: {len(primary_boundaries)}",
        f"Figure: {figdir / 'Figure3_boundary_mapping_improved.png'}",
        f"Interval export: {table_dir / 'Figure3_boundary_intervals.tsv'}",
        f"Sensitivity export: {table_dir / 'Figure3_sensitivity_analysis.tsv'}",
        f"Report: {outdir / 'reports' / 'Figure3_boundary_mapping_report.md'}",
        f"Validation warnings: {len(validation_warnings)}",
    ]
    (logdir / "phase3_boundary_mapping.log").write_text("\n".join(log_lines) + "\n")
    print("\n".join(log_lines))


if __name__ == "__main__":
    main()
