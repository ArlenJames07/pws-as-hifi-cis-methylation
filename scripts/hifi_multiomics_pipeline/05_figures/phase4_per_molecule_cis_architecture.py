#!/usr/bin/env python3
"""
Render the restructured Figure 4 from direct BAM-derived per-molecule methylation.

Panels:
A. Control paternal single-molecule SNORD116 heatmap
B. Control maternal single-molecule SNORD116 heatmap
C. Control boundary-crossing molecule heatmap centered on the SNHG14 boundary
D. Per-molecule entropy across regions
E. Per-molecule mean methylation contrast at the SNORD116 display window
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from scipy import stats


DEFAULT_BAM_DIR = Path("/mnt/diskrare/arlenb/08/hiphase_results/bamfiles")
DEFAULT_GTF = Path("/home/rare/arlen/reference/chm13v22.sorted.gtf")
DEFAULT_ASSIGNMENT = Path("/home/rare/arlen/paper_vf/tables/Figure1C_parental_assignment.tsv")
DEFAULT_OUTDIR = Path("/home/rare/arlen/paper_vf")

CONTROL_SAMPLES = ["017C", "018C"]
GROUP_ORDER = ["Control", "PWS-DEL", "PWS-mUPD", "AS-DEL"]
PANEL_E_ORDER = [
    "Control paternal",
    "AS-DEL retained paternal",
    "Control maternal",
    "PWS-DEL retained maternal",
    "PWS-mUPD hap1",
    "PWS-mUPD hap2",
]
PANEL_D_GROUP_ORDER = ["Control", "AS-DEL", "PWS-DEL", "PWS-mUPD"]
PANEL_D_REGION_ORDER = ["PWS-IC", "SNORD116 display window", "Downstream control"]
BOUNDARY_PARENT_ORDER = ["paternal", "maternal"]

HEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "methylation_probability",
    ["#2f6fb0", "#f5f4ef", "#d24a43"],
)
HEATMAP_CMAP.set_bad("#d7d7d7")

PANEL_E_PALETTE = {
    "Control paternal": "#2f6fb0",
    "AS-DEL retained paternal": "#5f93c9",
    "Control maternal": "#bf3d3d",
    "PWS-DEL retained maternal": "#d9796b",
    "PWS-mUPD hap1": "#e59a92",
    "PWS-mUPD hap2": "#cf6c67",
}
GROUP_PALETTE = {
    "Control": "#6e6e6e",
    "AS-DEL": "#2f6fb0",
    "PWS-DEL": "#bf3d3d",
    "PWS-mUPD": "#d9796b",
}
BOUNDARY_PARENT_PALETTE = {
    "paternal": "#2f6fb0",
    "maternal": "#bf3d3d",
}


@dataclass(frozen=True)
class RegionSpec:
    name: str
    chrom: str
    start: int
    end: int
    min_cpgs: int

    @property
    def span_bp(self) -> int:
        return self.end - self.start

    @property
    def label(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


SNORD116_DISPLAY = RegionSpec("SNORD116 display window", "chr15", 22_808_000, 22_845_000, 20)
PWS_IC = RegionSpec("PWS-IC", "chr15", 22_691_258, 22_693_494, 5)
DOWNSTREAM_CONTROL = RegionSpec("Downstream control", "chr15", 23_175_683, 23_275_682, 10)
BOUNDARY_CENTER = 22_693_725
BOUNDARY_REGION = RegionSpec("SNHG14 boundary window", "chr15", BOUNDARY_CENTER - 15_000, BOUNDARY_CENTER + 15_000, 20)
CONTROL_ENTRY_BOUNDARY = 22_690_600
CONTROL_EXIT_BOUNDARY = 22_696_850
ICR_START = 22_691_258
ICR_END = 22_693_494

PANEL_OUTPUTS = {
    "panel_a": "Figure4_panelA_control_paternal_SNORD116",
    "panel_b": "Figure4_panelB_control_maternal_SNORD116",
    "panel_c": "Figure4_panelC_boundary_crossing",
}


def log(message: str) -> None:
    print(f"[Phase4 Figure4] {message}", file=sys.stderr, flush=True)


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def detect_sample(path: Path) -> str | None:
    match = re.search(r"(?<!\d)(\d{3}[A-Z])(?![A-Za-z0-9])", path.name)
    return match.group(1) if match else None


def fnum(value: float | int | str | None, digits: int = 3) -> str:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(x):
        return "NA"
    return f"{x:.{digits}f}"


def fmt_pvalue(value: float | int | str | None) -> str:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(x):
        return "NA"
    if x < 1e-4:
        return f"{x:.2e}"
    return f"{x:.4f}"


def shannon_entropy(binary_states: np.ndarray) -> float:
    if binary_states.size == 0:
        return math.nan
    p = float(binary_states.mean())
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return float(-(p * math.log2(p) + (1.0 - p) * math.log2(1.0 - p)))


def discover_bam_paths(bam_dir: Path, sample_ids: list[str], outdir: Path) -> dict[str, Path]:
    table_path = outdir / "tables" / "Figure4_input_bam_files.tsv"
    out: dict[str, Path] = {}
    if table_path.exists():
        table = pd.read_csv(table_path, sep="\t")
        if {"sample_id", "bam_path"}.issubset(table.columns):
            for row in table.itertuples(index=False):
                sample_id = str(row.sample_id)
                bam_path = Path(str(row.bam_path))
                if sample_id in sample_ids and bam_path.exists():
                    out[sample_id] = bam_path
    if len(out) == len(sample_ids):
        return out
    for path in sorted(bam_dir.glob("*.bam")):
        sample = detect_sample(path)
        if sample is None or sample not in sample_ids:
            continue
        if sample not in out or ("v2" not in path.name.lower() and "v2" in out[sample].name.lower()):
            out[sample] = path
    return out


def load_assignment_table(path: Path) -> tuple[pd.DataFrame, dict[str, str], dict[str, dict[int, str]]]:
    table = pd.read_csv(path, sep="\t")
    required = {"sample_id", "molecular_mechanism", "haplotype_label", "parental_assignment"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"Assignment table missing columns: {sorted(missing)}")

    sample_group_map: dict[str, str] = {}
    control_hp_map: dict[str, dict[int, str]] = {}
    for row in table.itertuples(index=False):
        sample_id = str(row.sample_id)
        mechanism = str(row.molecular_mechanism)
        sample_group_map[sample_id] = mechanism
        hap_label = str(row.haplotype_label)
        assignment = str(row.parental_assignment)
        if mechanism == "Control" and hap_label in {"hap1", "hap2"} and assignment in {"maternal", "paternal"}:
            hp = 1 if hap_label == "hap1" else 2
            control_hp_map.setdefault(sample_id, {})[hp] = assignment
    return table, sample_group_map, control_hp_map


def build_bam_input_table(sample_group_map: dict[str, str], bam_paths: dict[str, Path]) -> pd.DataFrame:
    rows = []
    for sample_id in sorted(bam_paths):
        bam_path = bam_paths[sample_id]
        rows.append(
            {
                "sample_id": sample_id,
                "sample_group": sample_group_map.get(sample_id, "NA"),
                "bam_path": str(bam_path),
                "bai_exists": bam_path.with_suffix(bam_path.suffix + ".bai").exists() or Path(str(bam_path) + ".bai").exists(),
            }
        )
    return pd.DataFrame(rows)


def resolve_analysis_label(
    sample_id: str,
    sample_group: str,
    hp_tag: int | None,
    control_hp_map: dict[str, dict[int, str]],
) -> tuple[str, str] | tuple[None, None]:
    if sample_group == "Control":
        if hp_tag is None:
            return None, None
        parent = control_hp_map.get(sample_id, {}).get(hp_tag)
        if parent == "maternal":
            return "Control maternal", "maternal"
        if parent == "paternal":
            return "Control paternal", "paternal"
        return None, None
    if sample_group == "AS-DEL":
        return "AS-DEL retained paternal", "paternal"
    if sample_group == "PWS-DEL":
        return "PWS-DEL retained maternal", "maternal"
    if sample_group == "PWS-mUPD":
        if hp_tag == 1:
            return "PWS-mUPD hap1", "hap1"
        if hp_tag == 2:
            return "PWS-mUPD hap2", "hap2"
    return None, None


def pick_modified_base_key(read: pysam.AlignedSegment) -> tuple[str, int, str] | None:
    if not read.modified_bases:
        return None
    for key in read.modified_bases:
        if len(key) >= 3 and key[0] == "C" and key[2] == "m":
            return key
    return None


def extract_region_data(
    bam_path: Path,
    sample_id: str,
    sample_group: str,
    region: RegionSpec,
    control_hp_map: dict[str, dict[int, str]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    call_rows: list[dict[str, object]] = []
    read_rows: list[dict[str, object]] = []
    bam = pysam.AlignmentFile(bam_path, "rb")
    for read in bam.fetch(region.chrom, region.start, region.end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < 20:
            continue
        try:
            hp_tag = int(read.get_tag("HP"))
        except Exception:
            hp_tag = None
        analysis_label, parental_label = resolve_analysis_label(sample_id, sample_group, hp_tag, control_hp_map)
        if analysis_label is None:
            continue
        key = pick_modified_base_key(read)
        if key is None:
            continue
        q2r = {
            query_pos: ref_pos
            for query_pos, ref_pos in read.get_aligned_pairs(matches_only=False)
            if query_pos is not None and ref_pos is not None
        }
        positions: list[int] = []
        probs: list[float] = []
        for query_pos, ml in read.modified_bases[key]:
            ref_pos = q2r.get(query_pos)
            if ref_pos is None:
                continue
            ref_1based = ref_pos + 1
            if region.start <= ref_1based <= region.end:
                prob = float(ml) / 255.0
                positions.append(ref_1based)
                probs.append(prob)
        if len(positions) < region.min_cpgs:
            continue
        order = np.argsort(np.array(positions))
        positions = [positions[i] for i in order]
        probs_arr = np.array([probs[i] for i in order], dtype=float)
        binary_arr = (probs_arr >= 0.5).astype(float)
        read_key = f"{sample_id}|{read.query_name}"
        read_rows.append(
            {
                "sample_id": sample_id,
                "sample_group": sample_group,
                "analysis_label": analysis_label,
                "parental_label": parental_label,
                "region": region.name,
                "chrom": region.chrom,
                "region_start": region.start,
                "region_end": region.end,
                "read_id": read.query_name,
                "read_key": read_key,
                "hp_tag": hp_tag if hp_tag is not None else "NA",
                "alignment_start": int(read.reference_start + 1),
                "alignment_end": int(read.reference_end),
                "n_cpgs": int(len(positions)),
                "mean_methylation_probability": float(probs_arr.mean()),
                "mean_methylation_binary": float(binary_arr.mean()),
                "methylation_entropy": shannon_entropy(binary_arr),
            }
        )
        for pos, prob, binary in zip(positions, probs_arr.tolist(), binary_arr.tolist()):
            call_rows.append(
                {
                    "sample_id": sample_id,
                    "sample_group": sample_group,
                    "analysis_label": analysis_label,
                    "parental_label": parental_label,
                    "region": region.name,
                    "chrom": region.chrom,
                    "region_start": region.start,
                    "region_end": region.end,
                    "read_id": read.query_name,
                    "read_key": read_key,
                    "hp_tag": hp_tag if hp_tag is not None else "NA",
                    "cpg_position": int(pos),
                    "methylation_probability": float(prob),
                    "methylation_binary": int(binary),
                }
            )
    bam.close()
    return read_rows, call_rows


def add_boundary_metrics(read_df: pd.DataFrame, call_df: pd.DataFrame) -> pd.DataFrame:
    if read_df.empty:
        return read_df
    out = read_df.copy()
    extras = []
    grouped = call_df.sort_values(["read_key", "cpg_position"]).groupby("read_key")
    for read_key, group in grouped:
        positions = group["cpg_position"].to_numpy(dtype=int)
        probs = group["methylation_probability"].to_numpy(dtype=float)
        binary = group["methylation_binary"].to_numpy(dtype=float)
        left_mask = positions < BOUNDARY_CENTER
        right_mask = positions >= BOUNDARY_CENTER
        n_left = int(left_mask.sum())
        n_right = int(right_mask.sum())
        left_mean = float(probs[left_mask].mean()) if n_left else math.nan
        right_mean = float(probs[right_mask].mean()) if n_right else math.nan
        crossing = n_left >= 5 and n_right >= 5
        transition_position = math.nan
        transition_offset = math.nan
        transition_delta = math.nan
        if crossing:
            best_score = -1.0
            best_position = None
            best_delta = None
            for idx in range(5, len(positions) - 5):
                left_state = binary[:idx]
                right_state = binary[idx:]
                delta = float(right_state.mean() - left_state.mean())
                score = abs(delta)
                if score > best_score:
                    best_score = score
                    best_position = int(positions[idx])
                    best_delta = delta
            if best_position is not None:
                transition_position = best_position
                transition_offset = best_position - BOUNDARY_CENTER
                transition_delta = best_delta
        extras.append(
            {
                "read_key": read_key,
                "n_left_cpgs": n_left,
                "n_right_cpgs": n_right,
                "left_mean_probability": left_mean,
                "right_mean_probability": right_mean,
                "delta_right_minus_left": right_mean - left_mean if crossing else math.nan,
                "boundary_crossing": crossing,
                "transition_position": transition_position,
                "transition_offset_bp": transition_offset,
                "transition_delta_binary": transition_delta,
            }
        )
    extra_df = pd.DataFrame(extras)
    return out.merge(extra_df, on="read_key", how="left")


def pivot_heatmap(
    call_df: pd.DataFrame,
    read_df: pd.DataFrame,
    value_col: str,
    order_cols: list[str],
    sort_cols: list[str],
    ascending: list[bool],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if read_df.empty or call_df.empty:
        return pd.DataFrame(), pd.DataFrame()
    meta = read_df.copy().sort_values(sort_cols, ascending=ascending).reset_index(drop=True)
    matrix = call_df.pivot_table(index="read_key", columns="cpg_position", values=value_col, aggfunc="mean")
    matrix = matrix.reindex(columns=sorted(matrix.columns.tolist()))
    meta = meta.drop_duplicates("read_key").set_index("read_key")
    matrix = matrix.reindex(meta.index.tolist())
    meta = meta.reset_index()[order_cols]
    meta = meta.set_index("read_key")
    return matrix, meta


def nearest_column_index(columns: np.ndarray, position: float) -> int:
    return int(np.argmin(np.abs(columns - position)))


def tick_labels_from_positions(
    positions: np.ndarray,
    n_ticks: int = 6,
    relative_to: int | None = None,
) -> tuple[list[int], list[str]]:
    if positions.size == 0:
        return [], []
    if positions.size <= n_ticks:
        idx = list(range(positions.size))
    else:
        idx = np.linspace(0, positions.size - 1, n_ticks).round().astype(int).tolist()
    labels = []
    for i in idx:
        pos = int(positions[i])
        if relative_to is None:
            labels.append(f"{pos / 1_000_000:.3f}")
        else:
            labels.append(f"{(pos - relative_to) / 1000:.0f}")
    return idx, labels


def fixed_relative_ticks(
    positions: np.ndarray,
    center: int,
    tick_kb: list[int],
) -> tuple[list[int], list[str]]:
    if positions.size == 0:
        return [], []
    idx = []
    labels = []
    for kb in tick_kb:
        target = center + kb * 1000
        idx.append(nearest_column_index(positions, target))
        labels.append(str(kb))
    unique = []
    out_labels = []
    seen = set()
    for i, label in zip(idx, labels):
        if i in seen:
            continue
        seen.add(i)
        unique.append(i)
        out_labels.append(label)
    return unique, out_labels


def write_gzip_tsv(path: Path, df: pd.DataFrame) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        df.to_csv(handle, sep="\t", index=False)


def save_heatmap_tables(
    tbldir: Path,
    basename: str,
    matrix: pd.DataFrame,
    meta: pd.DataFrame,
) -> None:
    if matrix.empty or meta.empty:
        return
    out_matrix = matrix.copy()
    out_matrix.index.name = "read_key"
    write_gzip_tsv(tbldir / f"{basename}_matrix.tsv.gz", out_matrix.reset_index())
    meta.reset_index().to_csv(tbldir / f"{basename}_rows.tsv", sep="\t", index=False)


def load_gtf_features(gtf_path: Path, regions: list[tuple[str, RegionSpec]]) -> pd.DataFrame:
    feature_rows: list[dict[str, object]] = []
    wanted_names = {
        "snord116": lambda gene: gene.startswith("SNORD116"),
        "boundary": lambda gene: gene in {"SNHG14", "SNURF"},
    }
    with open(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene" or fields[0] != "chr15":
                continue
            start = int(fields[3])
            end = int(fields[4])
            attrs = {}
            for part in fields[8].split(";"):
                part = part.strip()
                if not part or " " not in part:
                    continue
                key, value = part.split(" ", 1)
                attrs[key] = value.strip().strip('"')
            gene_name = attrs.get("gene_name") or attrs.get("gene_id") or "NA"
            for panel_name, region in regions:
                if end < region.start or start > region.end:
                    continue
                if wanted_names[panel_name](gene_name):
                    feature_rows.append(
                        {
                            "panel": panel_name,
                            "feature_name": gene_name,
                            "chrom": fields[0],
                            "start": start,
                            "end": end,
                        }
                    )
    # Add synthetic regulatory bar for the boundary panel.
    feature_rows.extend(
        [
            {
                "panel": "boundary",
                "feature_name": "Boundary span",
                "chrom": BOUNDARY_REGION.chrom,
                "start": CONTROL_ENTRY_BOUNDARY,
                "end": CONTROL_EXIT_BOUNDARY,
            },
        ]
    )
    return pd.DataFrame(feature_rows)


def smooth_values(values: np.ndarray, window: int = 9) -> np.ndarray:
    if values.size == 0:
        return values
    return pd.Series(values).rolling(window=window, center=True, min_periods=1).mean().to_numpy(dtype=float)


def compress_matrix_for_display(matrix: pd.DataFrame, max_columns: int) -> pd.DataFrame:
    if matrix.empty or matrix.shape[1] <= max_columns:
        return matrix.copy()
    n_cols = matrix.shape[1]
    bin_ids = np.floor(np.arange(n_cols) * max_columns / n_cols).astype(int)
    grouped = matrix.T.groupby(bin_ids).mean().T
    source_positions = matrix.columns.to_numpy(dtype=float)
    display_positions = []
    for bin_id in sorted(np.unique(bin_ids)):
        pos = float(np.mean(source_positions[bin_ids == bin_id]))
        display_positions.append(pos)
    grouped.columns = display_positions
    return grouped


def compute_column_profile(matrix: pd.DataFrame) -> np.ndarray:
    if matrix.empty:
        return np.array([], dtype=float)
    values = matrix.to_numpy(dtype=float)
    counts = np.isfinite(values).sum(axis=0)
    sums = np.nansum(values, axis=0)
    return np.divide(sums, counts, out=np.full(counts.shape, np.nan, dtype=float), where=counts > 0)


def render_mean_track(
    ax: plt.Axes,
    matrix: pd.DataFrame,
    line_color: str,
    ylabel: str = "Mean methylation",
    smooth_window: int = 11,
    line_width: float = 1.7,
    fill_alpha: float = 0.16,
) -> np.ndarray:
    if matrix.empty:
        ax.set_xticks([])
        ax.set_yticks([])
        return np.array([], dtype=float)
    profile = compute_column_profile(matrix)
    profile = smooth_values(profile, window=smooth_window)
    xs = np.arange(matrix.shape[1], dtype=float)
    ax.fill_between(xs, profile, color=line_color, alpha=fill_alpha)
    ax.plot(xs, profile, color=line_color, linewidth=line_width)
    ax.set_ylim(0.0, 1.0)
    ax.set_yticks([0.0, 1.0])
    ax.set_yticklabels(["0", "1"], fontsize=7)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#ececec", linewidth=0.6)
    ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    return profile


def render_feature_track(
    ax: plt.Axes,
    columns: np.ndarray,
    feature_df: pd.DataFrame,
    region: RegionSpec,
    relative_to: int | None = None,
    fixed_relative_tick_kb: list[int] | None = None,
) -> None:
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.grid(False)
    if columns.size == 0:
        ax.set_xticks([])
        return
    panel = feature_df["panel"].iloc[0] if not feature_df.empty and "panel" in feature_df.columns else "unknown"
    y0 = 0.28 if panel == "snord116" else 0.20
    height = 0.24 if panel == "snord116" else 0.42
    for idx, row in feature_df.reset_index(drop=True).iterrows():
        start_idx = nearest_column_index(columns, max(row["start"], region.start))
        end_idx = nearest_column_index(columns, min(row["end"], region.end))
        width = max(end_idx - start_idx + 1, 1)
        color = "#7f7f7f"
        alpha = 0.9
        if row["feature_name"] == "Boundary span":
            color = "#4C78A8"
            alpha = 0.28
        ax.add_patch(Rectangle((start_idx - 0.5, y0), width, height, facecolor=color, edgecolor="none", alpha=alpha))
        if panel == "boundary" and row["feature_name"] in {"SNHG14", "SNURF", "Boundary span"}:
            ax.text(start_idx, 0.80, row["feature_name"], ha="left", va="bottom", fontsize=8.5)
    if panel == "snord116":
        ax.text(0.00, 0.78, "SNORD116 copies", transform=ax.transAxes, ha="left", va="bottom", fontsize=7.5, color="#666666")
    if relative_to is not None and fixed_relative_tick_kb is not None:
        ticks, labels = fixed_relative_ticks(columns, relative_to, fixed_relative_tick_kb)
    else:
        ticks, labels = tick_labels_from_positions(columns, relative_to=relative_to)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=8)


def render_snord116_panel(
    fig: plt.Figure,
    spec,
    matrix: pd.DataFrame,
    meta: pd.DataFrame,
    region: RegionSpec,
    title: str,
    line_color: str,
    feature_df: pd.DataFrame,
    show_ylabel: bool,
) -> matplotlib.image.AxesImage:
    sub = spec.subgridspec(3, 1, height_ratios=[0.32, 1.0, 0.14], hspace=0.04)
    ax_top = fig.add_subplot(sub[0])
    ax_main = fig.add_subplot(sub[1], sharex=ax_top)
    ax_feat = fig.add_subplot(sub[2], sharex=ax_top)

    render_mean_track(ax_top, matrix, line_color, ylabel="Mean meth.", smooth_window=17, line_width=1.95, fill_alpha=0.18)
    ax_top.set_title(title, loc="left", fontsize=13.5, fontweight="bold", color=line_color, pad=3.5)

    image = ax_main.imshow(
        matrix.to_numpy(dtype=float),
        aspect="auto",
        interpolation="nearest",
        cmap=HEATMAP_CMAP,
        vmin=0.0,
        vmax=1.0,
    )
    positions = matrix.columns.to_numpy(dtype=float)
    ticks, labels = tick_labels_from_positions(positions)
    ax_main.set_xticks(ticks)
    ax_main.set_xticklabels(labels, fontsize=9)
    ax_main.set_yticks([])
    ax_main.set_ylabel("Rows = molecules\n(sorted by mean methylation)" if show_ylabel else "", fontsize=10.5)
    ax_main.grid(False)
    ax_main.text(
        0.012,
        0.04,
        "rows = molecules; columns = CpGs",
        transform=ax_main.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.4,
        color="#5a5a5a",
        bbox={"facecolor": "white", "alpha": 0.74, "edgecolor": "none", "pad": 1.2},
    )
    ax_main.text(
        0.995,
        0.98,
        (
            f"n = {len(meta)}\n"
            f"median mean meth. = "
            f"{fnum(pd.to_numeric(meta['mean_methylation_probability'], errors='coerce').median())}"
        ),
        transform=ax_main.transAxes,
        ha="right",
        va="top",
        fontsize=8.8,
        bbox={"facecolor": "white", "alpha": 0.86, "edgecolor": "none", "pad": 2.0},
    )

    render_feature_track(ax_feat, positions, feature_df, region)
    ax_feat.set_xlabel("Shared CpG positions across SNORD116 (chr15 Mb, T2T-CHM13)", fontsize=10.5)
    return image


def render_boundary_panel(
    fig: plt.Figure,
    spec,
    matrix: pd.DataFrame,
    meta: pd.DataFrame,
    region: RegionSpec,
    feature_df: pd.DataFrame,
    title: str = "C. Boundary-centered molecules at the SNHG14 boundary",
) -> None:
    positions = matrix.columns.to_numpy(dtype=float)
    boundary_idx = nearest_column_index(positions, BOUNDARY_CENTER)
    span_start_idx = nearest_column_index(positions, CONTROL_ENTRY_BOUNDARY)
    span_end_idx = nearest_column_index(positions, CONTROL_EXIT_BOUNDARY)
    pat_meta = meta[meta["parental_label"] == "paternal"].copy()
    mat_meta = meta[meta["parental_label"] == "maternal"].copy()
    pat_matrix = matrix.reindex(pat_meta.index.tolist())
    mat_matrix = matrix.reindex(mat_meta.index.tolist())

    sub = spec.subgridspec(5, 1, height_ratios=[0.34, 0.82, 0.06, 0.82, 0.19], hspace=0.03)
    ax_top = fig.add_subplot(sub[0])
    ax_pat = fig.add_subplot(sub[1], sharex=ax_top)
    ax_gap = fig.add_subplot(sub[2], sharex=ax_top)
    ax_mat = fig.add_subplot(sub[3], sharex=ax_top)
    ax_feat = fig.add_subplot(sub[4], sharex=ax_top)
    ax_gap.axis("off")

    if not pat_matrix.empty:
        pat_profile = smooth_values(compute_column_profile(pat_matrix), window=19)
        ax_top.plot(np.arange(pat_matrix.shape[1]), pat_profile, color=BOUNDARY_PARENT_PALETTE["paternal"], linewidth=2.0, label="Paternal-like")
    if not mat_matrix.empty:
        mat_profile = smooth_values(compute_column_profile(mat_matrix), window=19)
        ax_top.plot(np.arange(mat_matrix.shape[1]), mat_profile, color=BOUNDARY_PARENT_PALETTE["maternal"], linewidth=2.0, label="Maternal-like")
    for ax in [ax_top, ax_pat, ax_mat, ax_feat]:
        ax.axvspan(span_start_idx - 0.5, span_end_idx + 0.5, color="#8fa5bb", alpha=0.12, zorder=-5)
    ax_top.axvline(boundary_idx, color="#303030", linewidth=1.45, linestyle="--", alpha=0.95)
    ax_top.set_ylim(0.0, 1.0)
    ax_top.set_yticks([0.0, 1.0])
    ax_top.set_yticklabels(["0", "1"], fontsize=7)
    ax_top.set_ylabel("Mean meth.", fontsize=8)
    ax_top.set_title(title, loc="left", fontsize=13.5, fontweight="bold")
    ax_top.legend(frameon=False, fontsize=9, loc="upper right", bbox_to_anchor=(1.0, 1.26), ncol=2, handlelength=2.6)
    ax_top.grid(axis="y", color="#ececec", linewidth=0.6)
    ax_top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_top.text(
        boundary_idx,
        1.06,
        "SNHG14 transition center",
        transform=ax_top.get_xaxis_transform(),
        ha="center",
        va="bottom",
        fontsize=9.4,
        fontweight="bold",
        color="#303030",
    )

    def _draw_boundary_heatmap(ax: plt.Axes, sub_matrix: pd.DataFrame, sub_meta: pd.DataFrame, label: str, color: str) -> None:
        ax.imshow(
            sub_matrix.to_numpy(dtype=float),
            aspect="auto",
            interpolation="nearest",
            cmap=HEATMAP_CMAP,
            vmin=0.0,
            vmax=1.0,
            zorder=1,
        )
        ax.axvline(boundary_idx, color="#303030", linewidth=1.45, linestyle="--", alpha=0.95, zorder=3)
        ax.set_yticks([])
        ax.grid(False)
        ax.set_ylabel(label, fontsize=10.5, color=color, fontweight="bold")
        ax.set_xticks([])
        ax.add_patch(Rectangle((-0.055, 0), 0.024, 1.0, transform=ax.transAxes, facecolor=color, edgecolor="none", clip_on=False))
        if not sub_meta.empty:
            trans = sub_meta["transition_position"].dropna()
            if not trans.empty:
                xs = [nearest_column_index(positions, float(pos)) for pos in trans.to_numpy(dtype=float)]
                ys = np.arange(len(xs))
                ax.scatter(xs, ys, s=12, c="#111111", alpha=0.72, linewidths=0, zorder=4)
            ax.text(
                0.995,
                0.96,
                f"n = {len(sub_meta)}",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=9.4,
                bbox={"facecolor": "white", "alpha": 0.90, "edgecolor": "none", "pad": 1.8},
            )

    _draw_boundary_heatmap(ax_pat, pat_matrix, pat_meta, "Paternal-like\nmolecules", BOUNDARY_PARENT_PALETTE["paternal"])
    _draw_boundary_heatmap(ax_mat, mat_matrix, mat_meta, "Maternal-like\nmolecules", BOUNDARY_PARENT_PALETTE["maternal"])

    render_feature_track(
        ax_feat,
        positions,
        feature_df,
        region,
        relative_to=BOUNDARY_CENTER,
        fixed_relative_tick_kb=[-15, -10, -5, 0, 5, 10, 15],
    )
    ax_feat.axvline(boundary_idx, color="#303030", linewidth=1.45, linestyle="--", alpha=0.95)
    ax_feat.set_xlabel("Position relative to SNHG14 transition center (kb)", fontsize=10.5)


def save_expanded_boundary_panel(
    figdir: Path,
    matrix: pd.DataFrame,
    meta: pd.DataFrame,
    region: RegionSpec,
    feature_df: pd.DataFrame,
) -> None:
    fig = plt.figure(figsize=(15.5, 8.7))
    outer = fig.add_gridspec(1, 1, left=0.07, right=0.98, top=0.93, bottom=0.14)
    render_boundary_panel(
        fig,
        outer[0, 0],
        matrix,
        meta,
        region,
        feature_df,
        title="Expanded boundary-centered molecules at the SNHG14 boundary",
    )
    cax = fig.add_axes([0.24, 0.055, 0.52, 0.020])
    scalar = plt.cm.ScalarMappable(cmap=HEATMAP_CMAP, norm=plt.Normalize(0.0, 1.0))
    scalar.set_array([])
    cbar = fig.colorbar(scalar, cax=cax, orientation="horizontal")
    cbar.set_label("CpG methylation value (0 = unmethylated, 1 = methylated)", fontsize=10)
    cbar.ax.tick_params(labelsize=9)
    base = figdir / "Figure4_panelC_boundary_transition_expanded"
    fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def build_entropy_stats(entropy_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for region in PANEL_D_REGION_ORDER:
        z = entropy_df[entropy_df["region"] == region].copy()
        groups = []
        medians = {}
        counts = {}
        for group in PANEL_D_GROUP_ORDER:
            vals = pd.to_numeric(
                z.loc[z["sample_group"] == group, "methylation_entropy"],
                errors="coerce",
            ).dropna()
            groups.append(vals)
            medians[group] = float(vals.median()) if not vals.empty else math.nan
            counts[group] = int(vals.size)
        valid_groups = [vals for vals in groups if len(vals) > 0]
        if len(valid_groups) >= 2:
            stat, pvalue = stats.kruskal(*valid_groups)
        else:
            stat, pvalue = math.nan, math.nan
        rows.append(
            {
                "region": region,
                "test": "Kruskal-Wallis across sample groups",
                "statistic_H": stat,
                "p_value": pvalue,
                "n_Control": counts["Control"],
                "n_PWS-DEL": counts["PWS-DEL"],
                "n_PWS-mUPD": counts["PWS-mUPD"],
                "n_AS-DEL": counts["AS-DEL"],
                "median_entropy_Control": medians["Control"],
                "median_entropy_PWS-DEL": medians["PWS-DEL"],
                "median_entropy_PWS-mUPD": medians["PWS-mUPD"],
                "median_entropy_AS-DEL": medians["AS-DEL"],
            }
        )
    return pd.DataFrame(rows)


def build_panel_e_stats(snord116_df: pd.DataFrame) -> tuple[pd.DataFrame, float, float]:
    rows = []
    valid_groups = []
    for label in PANEL_E_ORDER:
        vals = pd.to_numeric(
            snord116_df.loc[snord116_df["analysis_label"] == label, "mean_methylation_probability"],
            errors="coerce",
        ).dropna()
        valid_groups.append(vals)
        sample_count = snord116_df.loc[snord116_df["analysis_label"] == label, "sample_id"].nunique()
        rows.append(
            {
                "analysis_label": label,
                "n_reads": int(vals.size),
                "n_samples": int(sample_count),
                "median_mean_methylation": float(vals.median()) if not vals.empty else math.nan,
                "mean_mean_methylation": float(vals.mean()) if not vals.empty else math.nan,
                "q1_mean_methylation": float(vals.quantile(0.25)) if not vals.empty else math.nan,
                "q3_mean_methylation": float(vals.quantile(0.75)) if not vals.empty else math.nan,
            }
        )
    valid_groups = [vals for vals in valid_groups if len(vals) > 0]
    if len(valid_groups) >= 2:
        stat, pvalue = stats.kruskal(*valid_groups)
    else:
        stat, pvalue = math.nan, math.nan
    return pd.DataFrame(rows), stat, pvalue


def build_boundary_stats(boundary_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for parent in BOUNDARY_PARENT_ORDER:
        z = boundary_df[boundary_df["parental_label"] == parent].copy()
        rows.append(
            {
                "parental_label": parent,
                "n_reads": int(len(z)),
                "median_left_mean_probability": float(pd.to_numeric(z["left_mean_probability"], errors="coerce").median()),
                "median_right_mean_probability": float(pd.to_numeric(z["right_mean_probability"], errors="coerce").median()),
                "median_delta_right_minus_left": float(pd.to_numeric(z["delta_right_minus_left"], errors="coerce").median()),
                "median_transition_offset_bp": float(pd.to_numeric(z["transition_offset_bp"], errors="coerce").median()),
            }
        )
    return pd.DataFrame(rows)


def overlay_grouped_sample_medians(
    ax: plt.Axes,
    sample_medians: pd.DataFrame,
    x_order: list[str],
    hue_order: list[str],
    palette: dict[str, str],
    x_col: str,
    hue_col: str,
    y_col: str,
    total_width: float = 0.72,
) -> None:
    if sample_medians.empty:
        return
    group_width = total_width / max(len(hue_order), 1)
    for x_idx, x_value in enumerate(x_order):
        for hue_idx, hue_value in enumerate(hue_order):
            sub = sample_medians[(sample_medians[x_col] == x_value) & (sample_medians[hue_col] == hue_value)].copy()
            if sub.empty:
                continue
            center = x_idx - total_width / 2 + group_width / 2 + hue_idx * group_width
            if len(sub) == 1:
                xs = [center]
            else:
                xs = np.linspace(center - 0.045, center + 0.045, len(sub))
            ax.scatter(
                xs,
                sub[y_col].to_numpy(dtype=float),
                s=44,
                facecolor="white",
                edgecolor=palette[hue_value],
                linewidth=1.2,
                zorder=5,
            )


def overlay_category_sample_medians(
    ax: plt.Axes,
    sample_medians: pd.DataFrame,
    order: list[str],
    palette: dict[str, str],
    x_col: str,
    y_col: str,
) -> None:
    if sample_medians.empty:
        return
    for idx, label in enumerate(order):
        sub = sample_medians[sample_medians[x_col] == label].copy()
        if sub.empty:
            continue
        if len(sub) == 1:
            xs = [idx]
        else:
            xs = np.linspace(idx - 0.08, idx + 0.08, len(sub))
        ax.scatter(
            xs,
            sub[y_col].to_numpy(dtype=float),
            s=46,
            facecolor="white",
            edgecolor=palette[label],
            linewidth=1.2,
            zorder=5,
        )


def render_figure(
    outdir: Path,
    gtf_path: Path,
    all_read_df: pd.DataFrame,
    all_call_df: pd.DataFrame,
    feature_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    figdir = mkdir(outdir / "figures")
    tbldir = mkdir(outdir / "tables")

    snord116_reads = all_read_df[all_read_df["region"] == SNORD116_DISPLAY.name].copy()
    snord116_calls = all_call_df[all_call_df["region"] == SNORD116_DISPLAY.name].copy()
    entropy_reads = all_read_df[all_read_df["region"].isin(PANEL_D_REGION_ORDER)].copy()
    boundary_reads = all_read_df[all_read_df["region"] == BOUNDARY_REGION.name].copy()
    boundary_calls = all_call_df[all_call_df["region"] == BOUNDARY_REGION.name].copy()

    boundary_reads = add_boundary_metrics(boundary_reads, boundary_calls)
    boundary_reads = boundary_reads[boundary_reads["boundary_crossing"] == True].copy()
    boundary_calls = boundary_calls[boundary_calls["read_key"].isin(boundary_reads["read_key"])].copy()
    boundary_reads["parental_label"] = pd.Categorical(
        boundary_reads["parental_label"],
        categories=BOUNDARY_PARENT_ORDER,
        ordered=True,
    )
    boundary_reads["parent_order"] = boundary_reads["parental_label"].map(
        {label: idx for idx, label in enumerate(BOUNDARY_PARENT_ORDER)}
    )

    panel_a_reads = snord116_reads[snord116_reads["analysis_label"] == "Control paternal"].copy()
    panel_a_calls = snord116_calls[snord116_calls["analysis_label"] == "Control paternal"].copy()
    matrix_a, meta_a = pivot_heatmap(
        panel_a_calls,
        panel_a_reads,
        value_col="methylation_probability",
        order_cols=[
            "read_key",
            "sample_id",
            "sample_group",
            "analysis_label",
            "parental_label",
            "alignment_start",
            "alignment_end",
            "n_cpgs",
            "mean_methylation_probability",
            "mean_methylation_binary",
            "methylation_entropy",
        ],
        sort_cols=["mean_methylation_probability", "sample_id", "alignment_start"],
        ascending=[True, True, True],
    )

    panel_b_reads = snord116_reads[snord116_reads["analysis_label"] == "Control maternal"].copy()
    panel_b_calls = snord116_calls[snord116_calls["analysis_label"] == "Control maternal"].copy()
    matrix_b, meta_b = pivot_heatmap(
        panel_b_calls,
        panel_b_reads,
        value_col="methylation_probability",
        order_cols=[
            "read_key",
            "sample_id",
            "sample_group",
            "analysis_label",
            "parental_label",
            "alignment_start",
            "alignment_end",
            "n_cpgs",
            "mean_methylation_probability",
            "mean_methylation_binary",
            "methylation_entropy",
        ],
        sort_cols=["mean_methylation_probability", "sample_id", "alignment_start"],
        ascending=[True, True, True],
    )
    union_snord116_cols = sorted(set(matrix_a.columns.tolist()).union(matrix_b.columns.tolist()))
    if union_snord116_cols:
        matrix_a = matrix_a.reindex(columns=union_snord116_cols)
        matrix_b = matrix_b.reindex(columns=union_snord116_cols)
    save_heatmap_tables(tbldir, PANEL_OUTPUTS["panel_a"], matrix_a, meta_a)
    save_heatmap_tables(tbldir, PANEL_OUTPUTS["panel_b"], matrix_b, meta_b)

    matrix_c, meta_c = pivot_heatmap(
        boundary_calls,
        boundary_reads,
        value_col="methylation_binary",
        order_cols=[
            "read_key",
            "sample_id",
            "sample_group",
            "analysis_label",
            "parental_label",
            "alignment_start",
            "alignment_end",
            "n_cpgs",
            "mean_methylation_probability",
            "mean_methylation_binary",
            "methylation_entropy",
            "n_left_cpgs",
            "n_right_cpgs",
            "left_mean_probability",
            "right_mean_probability",
            "delta_right_minus_left",
            "transition_position",
            "transition_offset_bp",
            "transition_delta_binary",
        ],
        sort_cols=["parent_order", "mean_methylation_probability", "transition_offset_bp", "sample_id", "alignment_start"],
        ascending=[True, True, True, True, True],
    )
    if not matrix_c.empty:
        matrix_c = matrix_c.reindex(columns=sorted(matrix_c.columns.tolist()))
    save_heatmap_tables(tbldir, PANEL_OUTPUTS["panel_c"], matrix_c, meta_c)
    save_expanded_boundary_panel(figdir, matrix_c, meta_c, BOUNDARY_REGION, feature_df[feature_df["panel"] == "boundary"].copy())

    matrix_a_display = compress_matrix_for_display(matrix_a, max_columns=320)
    matrix_b_display = compress_matrix_for_display(matrix_b, max_columns=320)
    matrix_c_display = compress_matrix_for_display(matrix_c, max_columns=300)

    entropy_df = entropy_reads.copy()
    entropy_df["sample_group"] = pd.Categorical(entropy_df["sample_group"], PANEL_D_GROUP_ORDER, ordered=True)
    entropy_df["region"] = pd.Categorical(entropy_df["region"], PANEL_D_REGION_ORDER, ordered=True)
    entropy_stats = build_entropy_stats(entropy_df)
    entropy_sample_medians = (
        entropy_df.groupby(["region", "sample_group", "sample_id"], observed=True, as_index=False)["methylation_entropy"]
        .median()
        .rename(columns={"methylation_entropy": "sample_median_entropy"})
    )
    entropy_df.to_csv(tbldir / "Figure4_panelD_entropy_plot_input.tsv", sep="\t", index=False)
    entropy_stats.to_csv(tbldir / "Figure4_panelD_entropy_statistics.tsv", sep="\t", index=False)

    panel_e_df = snord116_reads[snord116_reads["analysis_label"].isin(PANEL_E_ORDER)].copy()
    panel_e_df["analysis_label"] = pd.Categorical(panel_e_df["analysis_label"], PANEL_E_ORDER, ordered=True)
    panel_e_sample_medians = (
        panel_e_df.groupby(["analysis_label", "sample_id"], observed=True, as_index=False)["mean_methylation_probability"]
        .median()
        .rename(columns={"mean_methylation_probability": "sample_median_mean_methylation"})
    )
    panel_e_df.to_csv(tbldir / "Figure4_panelE_SNORD116_per_read_summary.tsv", sep="\t", index=False)
    panel_e_stats, panel_e_statistic, panel_e_pvalue = build_panel_e_stats(panel_e_df)
    panel_e_stats["overall_kruskal_H"] = panel_e_statistic
    panel_e_stats["overall_kruskal_p"] = panel_e_pvalue
    panel_e_stats.to_csv(tbldir / "Figure4_panelE_SNORD116_group_statistics.tsv", sep="\t", index=False)

    boundary_stats = build_boundary_stats(boundary_reads)
    boundary_reads.to_csv(tbldir / "Figure4_panelC_boundary_crossing_summary.tsv", sep="\t", index=False)
    boundary_stats.to_csv(tbldir / "Figure4_panelC_boundary_transition_statistics.tsv", sep="\t", index=False)

    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 12.5,
            "axes.labelsize": 10.5,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "axes.facecolor": "#fcfcfc",
            "grid.color": "#e8e8e8",
            "grid.linewidth": 0.7,
        }
    )
    fig = plt.figure(figsize=(17.8, 15.8))
    fig.subplots_adjust(left=0.06, right=0.985, top=0.93, bottom=0.105)
    outer = fig.add_gridspec(3, 2, height_ratios=[1.00, 1.24, 0.95], hspace=0.45, wspace=0.24)

    feature_a = feature_df[feature_df["panel"] == "snord116"].copy()
    feature_c = feature_df[feature_df["panel"] == "boundary"].copy()

    image = render_snord116_panel(
        fig,
        outer[0, 0],
        matrix_a_display,
        meta_a,
        SNORD116_DISPLAY,
        "A. Control paternal molecules at SNORD116",
        BOUNDARY_PARENT_PALETTE["paternal"],
        feature_a,
        show_ylabel=True,
    )
    render_snord116_panel(
        fig,
        outer[0, 1],
        matrix_b_display,
        meta_b,
        SNORD116_DISPLAY,
        "B. Control maternal molecules at SNORD116",
        BOUNDARY_PARENT_PALETTE["maternal"],
        feature_a,
        show_ylabel=False,
    )

    render_boundary_panel(fig, outer[1, :], matrix_c_display, meta_c, BOUNDARY_REGION, feature_c)

    ax_d = fig.add_subplot(outer[2, 0])
    sns.boxplot(
        data=entropy_df,
        x="region",
        y="methylation_entropy",
        hue="sample_group",
        order=PANEL_D_REGION_ORDER,
        hue_order=PANEL_D_GROUP_ORDER,
        palette=GROUP_PALETTE,
        width=0.68,
        linewidth=0.9,
        showfliers=False,
        ax=ax_d,
    )
    overlay_grouped_sample_medians(
        ax_d,
        entropy_sample_medians,
        x_order=PANEL_D_REGION_ORDER,
        hue_order=PANEL_D_GROUP_ORDER,
        palette=GROUP_PALETTE,
        x_col="region",
        hue_col="sample_group",
        y_col="sample_median_entropy",
        total_width=0.68,
    )
    ax_d.set_title("D. Binary methylation entropy across regions", loc="left", fontsize=13.5, fontweight="bold", pad=8.0)
    ax_d.set_xlabel("")
    ax_d.set_ylabel("Binary methylation entropy", fontsize=10.5)
    ax_d.set_ylim(0.0, 1.08)
    ax_d.set_xticks(range(len(PANEL_D_REGION_ORDER)))
    ax_d.set_xticklabels(["PWS-IC", "SNORD116\nwindow", "Downstream\ncontrol"], fontsize=10)
    legend_d = ax_d.legend(title="", frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.34, 1.20), columnspacing=1.3, handletextpad=0.5)
    for text in legend_d.get_texts():
        text.set_fontsize(9)
    ax_d.grid(axis="y", color="#e5e5e5", linewidth=0.7)
    ax_d.set_axisbelow(True)
    entropy_lookup = entropy_stats.set_index("region")["p_value"].to_dict()
    ax_d.text(
        0.985,
        0.93,
        "\n".join(
            [
                f"PWS-IC p = {fmt_pvalue(entropy_lookup.get('PWS-IC'))}",
                f"SNORD116 p = {fmt_pvalue(entropy_lookup.get('SNORD116 display window'))}",
                f"Downstream p = {fmt_pvalue(entropy_lookup.get('Downstream control'))}",
            ]
        ),
        transform=ax_d.transAxes,
        ha="right",
        va="top",
        fontsize=8.5,
        bbox={"facecolor": "white", "alpha": 0.86, "edgecolor": "none", "pad": 1.8},
    )
    sns.despine(ax=ax_d)

    ax_e = fig.add_subplot(outer[2, 1])
    sns.violinplot(
        data=panel_e_df,
        x="analysis_label",
        y="mean_methylation_probability",
        hue="analysis_label",
        order=PANEL_E_ORDER,
        hue_order=PANEL_E_ORDER,
        palette=PANEL_E_PALETTE,
        dodge=False,
        cut=0,
        inner=None,
        linewidth=0.8,
        ax=ax_e,
    )
    if ax_e.legend_ is not None:
        ax_e.legend_.remove()
    sns.boxplot(
        data=panel_e_df,
        x="analysis_label",
        y="mean_methylation_probability",
        order=PANEL_E_ORDER,
        width=0.22,
        showcaps=True,
        boxprops={"facecolor": "white", "edgecolor": "#333333", "linewidth": 0.8, "zorder": 3},
        whiskerprops={"color": "#333333", "linewidth": 0.8},
        medianprops={"color": "#111111", "linewidth": 1.0},
        showfliers=False,
        ax=ax_e,
    )
    overlay_category_sample_medians(
        ax_e,
        panel_e_sample_medians,
        order=PANEL_E_ORDER,
        palette=PANEL_E_PALETTE,
        x_col="analysis_label",
        y_col="sample_median_mean_methylation",
    )
    ax_e.axvspan(-0.5, 1.5, color=BOUNDARY_PARENT_PALETTE["paternal"], alpha=0.05, zorder=0)
    ax_e.axvspan(1.5, len(PANEL_E_ORDER) - 0.5, color=BOUNDARY_PARENT_PALETTE["maternal"], alpha=0.04, zorder=0)
    ax_e.axvline(1.5, color="#bcbcbc", linewidth=1.2)
    ax_e.set_title(
        "E. SNORD116 methylation separates paternal-like and maternal-like configurations",
        loc="left",
        fontsize=13.5,
        fontweight="bold",
        pad=8.0,
    )
    ax_e.set_xlabel("")
    ax_e.set_ylabel("Mean methylation probability per read", fontsize=10.5)
    ax_e.set_ylim(0.0, 1.08)
    ax_e.set_xticks(range(len(PANEL_E_ORDER)))
    ax_e.set_xticklabels(
        ["Ctrl\npat", "AS-DEL\nret pat", "Ctrl\nmat", "PWS-DEL\nret mat", "mUPD\nhap1", "mUPD\nhap2"],
        fontsize=9,
    )
    ax_e.grid(axis="y", color="#e5e5e5", linewidth=0.7)
    ax_e.set_axisbelow(True)
    ax_e.text(
        0.985,
        0.98,
        f"Kruskal-Wallis p = {fmt_pvalue(panel_e_pvalue)}",
        transform=ax_e.transAxes,
        ha="right",
        va="top",
        fontsize=8.8,
        bbox={"facecolor": "white", "alpha": 0.86, "edgecolor": "none", "pad": 1.8},
    )
    ax_e.text(
        0.5,
        1.045,
        "Paternal-like",
        transform=ax_e.transData,
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="bold",
        color=BOUNDARY_PARENT_PALETTE["paternal"],
    )
    ax_e.text(
        3.5,
        1.045,
        "Maternal-like",
        transform=ax_e.transData,
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="bold",
        color=BOUNDARY_PARENT_PALETTE["maternal"],
    )
    bracket_y = 0.935
    ax_e.plot([1, 1, 3, 3], [bracket_y - 0.02, bracket_y, bracket_y, bracket_y - 0.02], color="#444444", linewidth=1.0, clip_on=False)
    ax_e.text(
        2.0,
        bracket_y + 0.015,
        "Reciprocal retained haplotypes",
        ha="center",
        va="bottom",
        fontsize=8.5,
        color="#444444",
    )
    med_lookup = panel_e_stats.set_index("analysis_label")["median_mean_methylation"].to_dict()
    for idx, label in enumerate(PANEL_E_ORDER):
        median_value = med_lookup.get(label)
        if median_value is None or not np.isfinite(median_value):
            continue
        ax_e.text(
            idx,
            min(float(median_value) + 0.035, 0.965),
            f"{float(median_value):.2f}",
            ha="center",
            va="bottom",
            fontsize=8.6,
            bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "none", "pad": 1.0},
        )
    sns.despine(ax=ax_e)

    fig.suptitle(
        "Figure 4. Single-molecule methylation architecture at SNORD116 and the SNHG14 boundary",
        fontsize=17,
        fontweight="bold",
        x=0.06,
        ha="left",
        y=0.985,
    )

    cax = fig.add_axes([0.26, 0.045, 0.48, 0.015])
    cbar = fig.colorbar(image, cax=cax, orientation="horizontal")
    cbar.set_label("CpG methylation value (0 = unmethylated, 1 = methylated)", fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    for base in [figdir / "Figure4", figdir / "Figure4_per_molecule_cis_architecture"]:
        fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight")
        fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return entropy_stats, panel_e_stats, boundary_stats


def write_run_summary(
    outdir: Path,
    all_read_df: pd.DataFrame,
    entropy_stats: pd.DataFrame,
    panel_e_stats: pd.DataFrame,
    boundary_stats: pd.DataFrame,
) -> None:
    summary_rows = [
        ["samples", str(all_read_df["sample_id"].nunique())],
        ["reads_direct_extracted", str(len(all_read_df))],
        ["snord116_display_window", SNORD116_DISPLAY.label],
        ["snord116_display_window_bp", str(SNORD116_DISPLAY.span_bp)],
        ["boundary_center", str(BOUNDARY_CENTER)],
        ["boundary_window", BOUNDARY_REGION.label],
        ["pws_ic_region", PWS_IC.label],
        ["downstream_control_region", DOWNSTREAM_CONTROL.label],
        ["panelA_reads_control_paternal", str(int(panel_e_stats.loc[panel_e_stats["analysis_label"] == "Control paternal", "n_reads"].sum()))],
        ["panelB_reads_control_maternal", str(int(panel_e_stats.loc[panel_e_stats["analysis_label"] == "Control maternal", "n_reads"].sum()))],
        ["panelC_boundary_crossing_reads", str(int(boundary_stats["n_reads"].sum()))],
        ["panelE_overall_categories", str(len(PANEL_E_ORDER))],
    ]
    run_summary = pd.DataFrame(summary_rows, columns=["metric", "value"])
    run_summary.to_csv(outdir / "phase4_run_summary.tsv", sep="\t", index=False)
    logdir = mkdir(outdir / "logs")
    run_summary.to_csv(logdir / "phase4_run_summary.tsv", sep="\t", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam-dir", type=Path, default=DEFAULT_BAM_DIR)
    parser.add_argument("--gtf", type=Path, default=DEFAULT_GTF)
    parser.add_argument("--parental-assignment", type=Path, default=DEFAULT_ASSIGNMENT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--fasta", type=Path, default=None, help="Accepted for compatibility with older wrappers.")
    parser.add_argument("--chrom", default="chr15", help="Accepted for compatibility with older wrappers.")
    parser.add_argument("--min-overlap-bp", type=int, default=5000, help="Accepted for compatibility with older wrappers.")
    parser.add_argument("--min-cpgs", type=int, default=5, help="Accepted for compatibility with older wrappers.")
    parser.add_argument("--min-ml-probability", type=float, default=0.70, help="Accepted for compatibility with older wrappers.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = mkdir(args.outdir)
    tbldir = mkdir(outdir / "tables")

    assignment_table, sample_group_map, control_hp_map = load_assignment_table(args.parental_assignment)
    sample_ids = sorted(sample_group_map)
    bam_paths = discover_bam_paths(args.bam_dir, sample_ids, outdir)
    missing = sorted(set(sample_ids).difference(bam_paths))
    if missing:
        raise FileNotFoundError(f"Missing BAMs for samples: {missing}")

    bam_input_df = build_bam_input_table(sample_group_map, bam_paths)
    bam_input_df.to_csv(tbldir / "Figure4_input_bam_files.tsv", sep="\t", index=False)

    region_specs = [SNORD116_DISPLAY, PWS_IC, DOWNSTREAM_CONTROL, BOUNDARY_REGION]
    read_rows: list[dict[str, object]] = []
    call_rows: list[dict[str, object]] = []
    for region in region_specs:
        log(f"Extracting region {region.name}: {region.label}")
        for sample_id in sorted(bam_paths):
            sample_group = sample_group_map[sample_id]
            if region == BOUNDARY_REGION and sample_group != "Control":
                continue
            rr, cr = extract_region_data(
                bam_path=bam_paths[sample_id],
                sample_id=sample_id,
                sample_group=sample_group,
                region=region,
                control_hp_map=control_hp_map,
            )
            read_rows.extend(rr)
            call_rows.extend(cr)

    all_read_df = pd.DataFrame(read_rows)
    all_call_df = pd.DataFrame(call_rows)
    if all_read_df.empty or all_call_df.empty:
        raise RuntimeError("No qualifying methylation calls were extracted from the BAMs.")

    write_gzip_tsv(tbldir / "Figure4_direct_region_read_summary.tsv.gz", all_read_df)
    write_gzip_tsv(tbldir / "Figure4_direct_region_cpg_calls.tsv.gz", all_call_df)

    feature_df = load_gtf_features(
        args.gtf,
        regions=[
            ("snord116", SNORD116_DISPLAY),
            ("boundary", BOUNDARY_REGION),
        ],
    )
    feature_df.to_csv(tbldir / "Figure4_gene_features.tsv", sep="\t", index=False)

    entropy_stats, panel_e_stats, boundary_stats = render_figure(
        outdir=outdir,
        gtf_path=args.gtf,
        all_read_df=all_read_df,
        all_call_df=all_call_df,
        feature_df=feature_df,
    )
    write_run_summary(outdir, all_read_df, entropy_stats, panel_e_stats, boundary_stats)
    log("Done.")


if __name__ == "__main__":
    main()
