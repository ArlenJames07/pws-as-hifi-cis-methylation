#!/usr/bin/env python3
"""
Create an improved Figure 2 for the Q1 genomics paper from existing Phase 2 tables.

This second-version figure keeps the original Phase 2 computations intact and
re-renders the interpretation as a cleaner reciprocal parental-architecture
figure with:
  - a schematic logic panel,
  - a shared T2T annotation track,
  - control, PWS-DEL, AS-DEL, and PWS-mUPD architecture panels,
  - a delta-architecture panel,
  - panel-specific plotting tables,
  - and a Markdown report.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import re
import sys
import textwrap
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.ticker import FuncFormatter

from paper_vf_phase2_reciprocal_cis_architecture import (
    BP_CLUSTER_INTERVALS_T2T,
    GREY,
    LIGHT_GREY,
    MATERNAL,
    MUPD_1,
    MUPD_2,
    PATERNAL,
    WindowSpec,
    pearson_summary,
    smooth_values,
    snrpn_marker_position,
)


DEFAULT_OUTDIR = Path("/home/rare/arlen/paper_vf")
DEFAULT_GTF = Path("/home/rare/arlen/reference/chm13v22.sorted.gtf")
DEFAULT_ICR_BED = Path("/home/rare/arlen/reference/ICR_t2t.bed")
DISPLAY_START = 18_000_000
DISPLAY_END = 28_000_000
DISPLAY_WINDOW = 1_000
STRUCT_GREY = "#8a8a8a"
DARK_GREY = "#3a3a3a"
SOFT_GREY = "#efefef"
BP_REGION_COLORS = {
    "BP1": "#c6a96e",
    "BP2": "#8faa92",
    "BP3": "#92a8c9",
}

ANNOTATION_GENE_SPECS = [
    ("NIPA1", "NIPA1", "top"),
    ("NIPA2", "NIPA2", "bottom"),
    ("CYFIP1", "CYFIP1", "top"),
    ("TUBGCP5", "TUBGCP5", "bottom"),
    ("MKRN3", "MKRN3", "top"),
    ("MAGEL2", "MAGEL2", "bottom"),
    ("NDN", "NDN", "top"),
    ("NPAP1", "NPAP1", "bottom"),
    ("SNORD116", "SNORD116", "top"),
    ("SNORD115", "SNORD115", "bottom"),
    ("UBE3A-ATS", "SNHG14", "top"),
    ("UBE3A", "UBE3A", "bottom"),
    ("ATP10A", "ATP10A", "top"),
    ("GABRB3", "GABRB3", "bottom"),
    ("GABRA5", "GABRA5", "top"),
    ("GABRG3", "GABRG3", "bottom"),
    ("OCA2", "OCA2", "top"),
    ("HERC2", "HERC2", "bottom"),
    ("APBA2", "APBA2", "top"),
    ("TJP1", "TJP1", "bottom"),
]


@dataclass(frozen=True)
class PanelStats:
    title: str
    stats_box: str
    badge: str


def log(message: str) -> None:
    print(f"[Figure2-improved] {message}", file=sys.stderr, flush=True)


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


def parse_gtf_attrs(attrs: str) -> Dict[str, str]:
    return {m.group(1): m.group(2) for m in re.finditer(r'(\S+) "([^"]*)";', attrs)}


def mb_formatter(value: float, _pos: float) -> str:
    return f"{value / 1e6:.1f}"


def overlap_filter(df: pd.DataFrame, start_col: str, end_col: str, spec: WindowSpec) -> pd.DataFrame:
    mask = (pd.to_numeric(df[start_col], errors="coerce") < spec.end) & (
        pd.to_numeric(df[end_col], errors="coerce") > spec.start
    )
    return df.loc[mask].copy()


def clipped_interval(start: float, end: float, spec: WindowSpec) -> Tuple[int, int]:
    return max(int(start), spec.start), min(int(end), spec.end)


def add_smoothed_columns(df: pd.DataFrame, columns: List[str], window: int) -> pd.DataFrame:
    out = df.copy()
    for col in columns:
        out[f"{col}_smoothed"] = smooth_values(out[col], window=window)
    return out


def structural_landmark_color(name: str) -> str:
    return BP_REGION_COLORS.get(name, STRUCT_GREY)


def logic_schematic_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "state": "Control",
                "allele_1": "maternal",
                "allele_2": "paternal",
                "interpretation": "Biparental reference",
            },
            {
                "state": "PWS-DEL",
                "allele_1": "retained_maternal",
                "allele_2": "paternal_deleted",
                "interpretation": "Retained maternal, paternal deleted",
            },
            {
                "state": "AS-DEL",
                "allele_1": "maternal_deleted",
                "allele_2": "retained_paternal",
                "interpretation": "Retained paternal, maternal deleted",
            },
            {
                "state": "PWS-mUPD",
                "allele_1": "maternal_haplotype_1",
                "allele_2": "maternal_haplotype_2",
                "interpretation": "Copy-neutral duplicated maternal state",
            },
        ]
    )


def load_annotation_genes(gtf_path: Path, spec: WindowSpec) -> pd.DataFrame:
    rows = []
    with open(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[0] != spec.chrom or fields[2] != "gene":
                continue
            start = int(fields[3])
            end = int(fields[4])
            if end < spec.start or start > spec.end:
                continue
            attrs = parse_gtf_attrs(fields[8])
            gene = attrs.get("gene") or attrs.get("gene_name") or attrs.get("gene_id", "")
            biotype = attrs.get("gene_biotype", "")
            description = attrs.get("description", "")
            rows.append(
                {
                    "gene": gene,
                    "start": start,
                    "end": end,
                    "strand": fields[6],
                    "gene_biotype": biotype,
                    "description": description,
                }
            )
    genes = pd.DataFrame(rows)
    if genes.empty:
        return genes
    genes = genes[~genes["gene"].str.startswith("LOC")].copy()
    genes = genes[~genes["gene_biotype"].str.contains("pseudogene", case=False, na=False)].copy()
    return genes.sort_values(["start", "end", "gene"]).reset_index(drop=True)


def load_pws_ic_anchor(icr_path: Path, spec: WindowSpec) -> pd.DataFrame:
    rows = []
    with open(icr_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6 or fields[0] != spec.chrom:
                continue
            start = int(fields[1])
            end = int(fields[2])
            if end < spec.start or start > spec.end:
                continue
            rows.append(
                {
                    "chrom": fields[0],
                    "start": start,
                    "end": end,
                    "name": fields[3],
                    "parent": fields[4],
                    "gene": fields[5],
                }
            )
    if not rows:
        return pd.DataFrame()
    icrs = pd.DataFrame(rows)
    icrs = icrs[icrs["gene"].astype(str).str.contains("SNRPN", na=False)].copy()
    return icrs.sort_values(["start", "end"]).reset_index(drop=True)


def select_gene_interval(genes: pd.DataFrame, query: str) -> pd.DataFrame:
    if genes.empty:
        return pd.DataFrame()
    if query in {"SNORD116", "SNORD115"}:
        return genes[genes["gene"].astype(str).str.startswith(f"{query}-")].copy()
    return genes[genes["gene"].astype(str) == query].copy()


def estimate_label_width_bp(label: str) -> int:
    return max(220_000, 24_000 * len(label) + 110_000)


def assign_non_overlapping_label_positions(
    feature_rows: List[Dict[str, object]],
    spec: WindowSpec,
    side: str,
    n_lanes: int = 4,
    gap_bp: int = 80_000,
    center_gap_bp: int = 450_000,
) -> List[Dict[str, object]]:
    if not feature_rows:
        return []
    lane_right = [spec.start - gap_bp] * n_lanes
    lane_last_center = [spec.start - center_gap_bp * 2] * n_lanes
    placed: List[Dict[str, object]] = []
    for row in sorted(feature_rows, key=lambda x: x["center"]):
        width = estimate_label_width_bp(str(row["label"]))
        chosen_lane = None
        chosen_left = None
        for lane in range(n_lanes):
            center_spacing_ok = (float(row["center"]) - lane_last_center[lane]) >= center_gap_bp
            left = max(int(row["center"] - width / 2), int(lane_right[lane] + gap_bp))
            if left + width <= spec.end and center_spacing_ok:
                chosen_lane = lane
                chosen_left = left
                break
        if chosen_lane is None:
            chosen_lane = int(np.argmin(lane_right))
            chosen_left = min(max(int(row["center"] - width / 2), spec.start), spec.end - width)
        lane_right[chosen_lane] = chosen_left + width
        lane_last_center[chosen_lane] = float(row["center"])
        out = dict(row)
        out["label_side"] = side
        out["lane_index"] = chosen_lane
        out["label_x"] = chosen_left + width / 2
        out["label_width_bp"] = width
        placed.append(out)
    return placed


def build_annotation_track_table(
    gene_features: pd.DataFrame,
    gtf_path: Path,
    icr_path: Path,
    spec: WindowSpec,
) -> pd.DataFrame:
    genes = load_annotation_genes(gtf_path, spec)
    pws_ic = load_pws_ic_anchor(icr_path, spec)
    feature_rows = []
    for cluster in BP_CLUSTER_INTERVALS_T2T:
        start, end = clipped_interval(cluster["start"], cluster["end"], spec)
        if start >= end:
            continue
        feature_rows.append(
            {
                "label": cluster["name"],
                "kind": "structural_landmark",
                "start": start,
                "end": end,
                "marker_position": np.nan,
                "lane": "structural",
                "color_role": "grey",
                "center": (start + end) / 2,
                "display_color": structural_landmark_color(cluster["name"]),
            }
        )

    snrpn_pos = snrpn_marker_position(gene_features, spec)
    if not pws_ic.empty:
        icr_row = pws_ic.iloc[0]
        start, end = clipped_interval(icr_row["start"], icr_row["end"], spec)
        if start < end:
            feature_rows.append(
                {
                    "label": "PWS-IC / SNRPN",
                    "kind": "regulatory_anchor",
                    "start": start,
                    "end": end,
                    "marker_position": snrpn_pos,
                    "lane": "anchor",
                    "color_role": "grey",
                    "center": (start + end) / 2,
                }
            )

    top_rows: List[Dict[str, object]] = []
    bottom_rows: List[Dict[str, object]] = []
    for label, query, side in ANNOTATION_GENE_SPECS:
        hit = select_gene_interval(genes, query)
        if hit.empty:
            continue
        start, end = clipped_interval(hit["start"].min(), hit["end"].max(), spec)
        if start >= end:
            continue
        row = {
            "label": label,
            "kind": "gene_context",
            "start": start,
            "end": end,
            "marker_position": np.nan,
            "lane": "gene_track",
            "color_role": "grey",
            "center": (start + end) / 2,
            "query_gene": query,
        }
        if side == "top":
            top_rows.append(row)
        else:
            bottom_rows.append(row)

    feature_rows.extend(assign_non_overlapping_label_positions(top_rows, spec, "top"))
    feature_rows.extend(assign_non_overlapping_label_positions(bottom_rows, spec, "bottom"))
    out = pd.DataFrame(feature_rows)
    if not out.empty:
        out = out.sort_values(["kind", "start", "end", "label"]).reset_index(drop=True)
    return out


def build_panel_tables(
    display_spec: WindowSpec,
    control_ref: pd.DataFrame,
    retained: pd.DataFrame,
    all_windows: pd.DataFrame,
    reciprocal_delta: pd.DataFrame,
    smooth_window: int,
) -> Dict[str, pd.DataFrame]:
    tables: Dict[str, pd.DataFrame] = {}

    panel_b = overlap_filter(control_ref, "window_start", "window_end", display_spec)
    panel_b = add_smoothed_columns(
        panel_b,
        ["control_maternal_mean", "control_paternal_mean"],
        smooth_window,
    )
    tables["Figure2B_control_reference_plot.tsv"] = panel_b

    pws_mean = (
        retained[retained["group"] == "PWS_DEL"]
        .groupby(["window_start", "window_end", "window_mid"], as_index=False)
        .agg(
            pwsdel_retained_maternal_mean=("retained_methylation", "mean"),
            n_pwsdel_samples=("sample", "nunique"),
        )
    )
    panel_c = overlap_filter(
        control_ref.merge(pws_mean, on=["window_start", "window_end", "window_mid"], how="left"),
        "window_start",
        "window_end",
        display_spec,
    )
    panel_c = add_smoothed_columns(
        panel_c,
        ["control_maternal_mean", "pwsdel_retained_maternal_mean"],
        smooth_window,
    )
    tables["Figure2C_pwsdel_retained_maternal_plot.tsv"] = panel_c

    as_mean = (
        retained[retained["group"] == "AS_DEL"]
        .groupby(["window_start", "window_end", "window_mid"], as_index=False)
        .agg(
            asdel_retained_paternal_mean=("retained_methylation", "mean"),
            n_asdel_samples=("sample", "nunique"),
        )
    )
    panel_d = overlap_filter(
        control_ref.merge(as_mean, on=["window_start", "window_end", "window_mid"], how="left"),
        "window_start",
        "window_end",
        display_spec,
    )
    panel_d = add_smoothed_columns(
        panel_d,
        ["control_paternal_mean", "asdel_retained_paternal_mean"],
        smooth_window,
    )
    tables["Figure2D_asdel_retained_paternal_plot.tsv"] = panel_d

    panel_e = overlap_filter(
        all_windows[all_windows["sample"] == "004P"][
            [
                "window_start",
                "window_end",
                "window_mid",
                "mean_meth_haplotype1",
                "mean_meth_haplotype2",
            ]
        ].merge(
            control_ref[["window_start", "control_maternal_mean"]],
            on="window_start",
            how="left",
        ),
        "window_start",
        "window_end",
        display_spec,
    )
    panel_e = add_smoothed_columns(
        panel_e,
        ["mean_meth_haplotype1", "mean_meth_haplotype2", "control_maternal_mean"],
        smooth_window,
    )
    tables["Figure2E_pwsmupd_maternal_plot.tsv"] = panel_e

    panel_f = overlap_filter(reciprocal_delta, "window_start", "window_end", display_spec)
    panel_f = add_smoothed_columns(
        panel_f,
        ["pwsdel_retained_maternal_mean", "asdel_retained_paternal_mean", "reciprocal_delta"],
        smooth_window,
    )
    panel_f["panel_label"] = "shared methylation scaffold"
    panel_f.loc[panel_f["delta_state"] == "paternal_higher", "panel_label"] = "localized_parental_divergence"
    tables["Figure2F_reciprocal_delta_architecture_plot.tsv"] = panel_f
    return tables


def panel_stats(
    control_ref: pd.DataFrame,
    corr: pd.DataFrame,
    reciprocal_delta: pd.DataFrame,
) -> Dict[str, PanelStats]:
    control_stats = pearson_summary(
        control_ref["control_maternal_mean"],
        control_ref["control_paternal_mean"],
    )
    pws = corr[corr["group"] == "PWS_DEL"].copy()
    asdel = corr[corr["group"] == "AS_DEL"].copy()
    upd = corr[corr["group"] == "PWS_mUPD"].copy().sort_values("comparison")
    reciprocal_nonnull = reciprocal_delta.dropna(
        subset=["pwsdel_retained_maternal_mean", "asdel_retained_paternal_mean"]
    )
    reciprocal_stats = pearson_summary(
        reciprocal_nonnull["pwsdel_retained_maternal_mean"],
        reciprocal_nonnull["asdel_retained_paternal_mean"],
    )

    return {
        "A": PanelStats(
            title="A. Controls",
            stats_box=f"r = {control_stats['pearson_r']:.3f}\nn = {control_stats['n_windows']:,} informative windows",
            badge="Biparental reference",
        ),
        "B": PanelStats(
            title="B. PWS-DEL",
            stats_box=(
                f"mean r = {pws['pearson_r'].mean():.3f}\n"
                f"range = {pws['pearson_r'].min():.3f}-{pws['pearson_r'].max():.3f}"
            ),
            badge="Maternal-only architecture exposed by paternal deletion",
        ),
        "C": PanelStats(
            title="C. AS-DEL",
            stats_box=(
                f"mean r = {asdel['pearson_r'].mean():.3f}\n"
                f"range = {asdel['pearson_r'].min():.3f}-{asdel['pearson_r'].max():.3f}"
            ),
            badge="Paternal-only architecture exposed by maternal deletion",
        ),
        "D": PanelStats(
            title="D. PWS-mUPD",
            stats_box=(
                f"hap1 r = {upd.iloc[0]['pearson_r']:.3f}\n"
                f"hap2 r = {upd.iloc[1]['pearson_r']:.3f}"
            ),
            badge="Copy-neutral duplicated maternal architecture",
        ),
        "E": PanelStats(
            title="E. Delta architecture",
            stats_box=(
                f"r = {reciprocal_stats['pearson_r']:.3f}\n"
                f"n = {reciprocal_stats['n_windows']:,} shared windows"
            ),
            badge="Shared methylation scaffold",
        ),
    }


def add_shared_axis_style(ax, spec: WindowSpec, ylabel: str, show_x: bool) -> None:
    ax.set_xlim(spec.start, spec.end)
    ax.grid(axis="y", color="#e8e8e8", lw=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.tick_params(axis="both", labelsize=10, width=0.9)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.xaxis.set_major_formatter(FuncFormatter(mb_formatter))
    ax.tick_params(axis="x", labelbottom=show_x)


def add_stats_box(ax, text: str) -> None:
    ax.text(
        1.015,
        0.98,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        color="#222222",
        clip_on=False,
        bbox=dict(fc="white", ec="#d0d0d0", lw=0.7, alpha=0.92, boxstyle="round,pad=0.35"),
    )


def add_badge(ax, text: str, edgecolor: str) -> None:
    ax.text(
        1.015,
        0.64,
        textwrap.fill(text, width=34),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        color="#1f1f1f",
        clip_on=False,
        bbox=dict(fc="white", ec=edgecolor, lw=0.9, alpha=0.92, boxstyle="round,pad=0.35"),
    )


def add_major_landmarks(ax, spec: WindowSpec, snrpn_pos: int) -> None:
    for cluster in BP_CLUSTER_INTERVALS_T2T:
        start, end = clipped_interval(cluster["start"], cluster["end"], spec)
        if start >= end:
            continue
        ax.axvspan(start, end, color=structural_landmark_color(cluster["name"]), alpha=0.10, lw=0, zorder=0)
    if spec.start <= snrpn_pos <= spec.end:
        ax.axvline(snrpn_pos, color=DARK_GREY, lw=1.0, alpha=0.65, zorder=2)


def add_logic_block(ax, x: float, title: str, top_label: str, bottom_label: str, top_style: str, bottom_style: str) -> None:
    width = 0.18
    box_h = 0.18
    ax.text(x + width / 2, 0.93, title, ha="center", va="top", fontsize=13, fontweight="bold")
    for y0, label, style in [
        (0.60, top_label, top_style),
        (0.35, bottom_label, bottom_style),
    ]:
        fc = "white"
        ec = STRUCT_GREY
        hatch = None
        text_color = "#222222"
        if style in {"maternal", "retained_maternal", "maternal_haplotype_1", "maternal_haplotype_2"}:
            fc = MATERNAL if style in {"maternal", "retained_maternal"} else MUPD_1 if style == "maternal_haplotype_1" else MUPD_2
            ec = fc
            text_color = "white"
        elif style in {"paternal", "retained_paternal"}:
            fc = PATERNAL
            ec = fc
            text_color = "white"
        elif "deleted" in style:
            fc = "#f4f4f4"
            ec = STRUCT_GREY
            hatch = "////"
        rect = FancyBboxPatch(
            (x, y0),
            width,
            box_h,
            boxstyle="round,pad=0.02,rounding_size=0.025",
            fc=fc,
            ec=ec,
            lw=1.0,
            hatch=hatch,
        )
        ax.add_patch(rect)
        ax.text(x + width / 2, y0 + box_h / 2, label, ha="center", va="center", fontsize=11, color=text_color)


def draw_logic_panel(ax) -> None:
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(
        0.0,
        1.02,
        "A. Reciprocal study logic",
        ha="left",
        va="bottom",
        fontsize=16,
        fontweight="bold",
    )
    blocks = [
        ("Control", "maternal", "paternal", "maternal", "paternal"),
        ("PWS-DEL", "retained maternal", "paternal deleted", "retained_maternal", "paternal_deleted"),
        ("AS-DEL", "maternal deleted", "retained paternal", "maternal_deleted", "retained_paternal"),
        ("PWS-mUPD", "maternal haplotype 1", "maternal haplotype 2", "maternal_haplotype_1", "maternal_haplotype_2"),
    ]
    xs = [0.03, 0.285, 0.54, 0.795]
    for x, (title, top_label, bottom_label, top_style, bottom_style) in zip(xs, blocks):
        add_logic_block(ax, x, title, top_label, bottom_label, top_style, bottom_style)
    ax.text(
        0.5,
        0.07,
        "Reciprocal deletions expose parental cis-methylation architecture, while PWS-mUPD confirms duplicated maternal identity.",
        ha="center",
        va="center",
        fontsize=12,
        color="#333333",
    )


def draw_annotation_track(ax, annotation_track: pd.DataFrame, spec: WindowSpec) -> None:
    ax.set_xlim(spec.start, spec.end)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.tick_params(axis="x", labelbottom=False)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.xaxis.set_major_formatter(FuncFormatter(mb_formatter))
    ax.hlines(0.73, spec.start, spec.end, color="#d8d8d8", lw=0.8)
    ax.hlines(0.58, spec.start, spec.end, color="#d8d8d8", lw=0.8)
    ax.hlines(0.43, spec.start, spec.end, color="#d8d8d8", lw=0.8)
    ax.text(
        0.0,
        1.06,
        "Shared T2T annotation track",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11.5,
        fontweight="bold",
    )

    lane_y = {
        "structural": 0.73,
        "anchor": 0.58,
        "gene_track": 0.43,
    }
    top_label_y = {0: 0.92, 1: 1.02, 2: 1.12, 3: 1.22}
    bottom_label_y = {0: 0.28, 1: 0.18, 2: 0.08, 3: -0.02}
    xtrans = ax.get_xaxis_transform()
    for _, row in annotation_track.iterrows():
        y = lane_y[row["lane"]]
        start, end = int(row["start"]), int(row["end"])
        if row["kind"] == "structural_landmark":
            rect_fc = row.get("display_color", STRUCT_GREY)
            rect_ec = row.get("display_color", STRUCT_GREY)
            alpha = 0.30
            height = 0.10
        elif row["kind"] == "regulatory_anchor":
            rect_fc = DARK_GREY
            rect_ec = DARK_GREY
            alpha = 0.25
            height = 0.10
        else:
            rect_fc = "#bdbdbd"
            rect_ec = "#7f7f7f"
            alpha = 0.45
            height = 0.09
        ax.add_patch(
            Rectangle(
                (start, y - height / 2),
                max(1, end - start),
                height,
                fc=rect_fc,
                ec=rect_ec,
                alpha=alpha,
                lw=0.8,
            )
        )
        center = float(row.get("center", (start + end) / 2))
        if row["kind"] in {"structural_landmark", "regulatory_anchor"}:
            label_y = 0.84 if row["kind"] == "structural_landmark" else 0.69
            ax.plot([center, center], [y + height / 2, label_y - 0.03], color="#757575", lw=0.7, alpha=0.8, transform=xtrans)
            ax.text(
                center,
                label_y,
                row["label"],
                transform=xtrans,
                ha="center",
                va="bottom",
                fontsize=10,
                color="#222222",
                clip_on=False,
            )
        else:
            lane_index = int(row.get("lane_index", 0))
            if row.get("label_side") == "top":
                label_y = top_label_y.get(lane_index, 1.20)
                line_end = label_y - 0.035
                va = "bottom"
            else:
                label_y = bottom_label_y.get(lane_index, 0.04)
                line_end = label_y + 0.035
                va = "top"
            label_x = float(row.get("label_x", center))
            ax.plot([center, center], [y + height / 2, 0.50], color="#969696", lw=0.55, alpha=0.75, transform=xtrans, clip_on=False)
            ax.plot([center, label_x], [0.50, line_end], color="#969696", lw=0.55, alpha=0.75, transform=xtrans, clip_on=False)
            ax.text(
                label_x,
                label_y,
                row["label"],
                transform=xtrans,
                ha="center",
                va=va,
                fontsize=8.6,
                color="#222222",
                clip_on=False,
            )
        if pd.notna(row["marker_position"]):
            marker = int(row["marker_position"])
            if spec.start <= marker <= spec.end:
                ax.axvline(marker, color=DARK_GREY, lw=1.1, alpha=0.75)

    ax.text(spec.start, 0.765, "BP landmarks", ha="left", va="bottom", fontsize=9.5, color="#5a5a5a")
    ax.text(spec.start, 0.615, "PWS imprinting center", ha="left", va="bottom", fontsize=9.5, color="#5a5a5a")
    ax.text(spec.start, 0.465, "Broader gene context", ha="left", va="bottom", fontsize=9.5, color="#5a5a5a")


def plot_series(ax, x: pd.Series, y: pd.Series, color: str, label: str, lw: float = 2.0, ls: str = "-") -> None:
    ax.plot(x, y, color=color, lw=lw, ls=ls, label=label)


def boundary_highlights(boundaries: pd.DataFrame, gene_features: pd.DataFrame) -> List[Tuple[int, int, str]]:
    if boundaries.empty or gene_features.empty:
        return []
    features = gene_features.set_index("feature")
    regions: List[Tuple[int, int, str]] = []
    for label, feature_names in [
        ("localized parental divergence", ["SNURF-SNRPN", "UBE3A-ATS"]),
        ("localized parental divergence", ["SNORD116", "UBE3A-ATS"]),
    ]:
        feature_intervals = []
        for name in feature_names:
            if name in features.index:
                row = features.loc[name]
                feature_intervals.append((int(row["start"]), int(row["end"])))
        if not feature_intervals:
            continue
        fstart = min(x[0] for x in feature_intervals)
        fend = max(x[1] for x in feature_intervals)
        overlap = boundaries[(boundaries["start"] < fend) & (boundaries["end"] > fstart)].copy()
        if overlap.empty:
            continue
        best = overlap.sort_values("max_abs_delta_smoothed", ascending=False).iloc[0]
        regions.append((int(best["start"]), int(best["end"]), label))
    unique = []
    seen = set()
    for start, end, label in regions:
        key = (start, end, label)
        if key not in seen:
            unique.append((start, end, label))
            seen.add(key)
    return unique


def render_figure(
    outdir: Path,
    display_spec: WindowSpec,
    annotation_track: pd.DataFrame,
    panel_tables: Dict[str, pd.DataFrame],
    boundaries: pd.DataFrame,
    gene_features: pd.DataFrame,
    stats: Dict[str, PanelStats],
) -> None:
    fig = plt.figure(figsize=(19.5, 15.9))
    gs = GridSpec(
        6,
        1,
        height_ratios=[1.10, 1.14, 1.14, 1.14, 1.14, 1.22],
        hspace=0.34,
        figure=fig,
    )

    anno_ax = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0], sharex=anno_ax)
    ax_c = fig.add_subplot(gs[2, 0], sharex=anno_ax)
    ax_d = fig.add_subplot(gs[3, 0], sharex=anno_ax)
    ax_e = fig.add_subplot(gs[4, 0], sharex=anno_ax)
    ax_f = fig.add_subplot(gs[5, 0], sharex=anno_ax)

    draw_annotation_track(anno_ax, annotation_track, display_spec)
    snrpn_markers = annotation_track.loc[annotation_track["kind"] == "regulatory_anchor", "marker_position"].dropna()
    if snrpn_markers.empty:
        raise ValueError("Annotation track is missing the PWS-IC/SNRPN anchor row.")
    snrpn_pos = int(snrpn_markers.iloc[0])

    # Panel A
    table_b = panel_tables["Figure2B_control_reference_plot.tsv"]
    add_shared_axis_style(ax_b, display_spec, "Mean methylation", show_x=False)
    ax_b.set_ylim(-0.02, 1.04)
    add_major_landmarks(ax_b, display_spec, snrpn_pos)
    plot_series(ax_b, table_b["window_mid"], table_b["control_maternal_mean_smoothed"], MATERNAL, "control maternal")
    plot_series(ax_b, table_b["window_mid"], table_b["control_paternal_mean_smoothed"], PATERNAL, "control paternal")
    ax_b.set_title(stats["A"].title, loc="left", fontsize=15, fontweight="bold", pad=8)
    add_badge(ax_b, stats["A"].badge, MATERNAL)
    add_stats_box(ax_b, stats["A"].stats_box)
    ax_b.legend(loc="lower left", bbox_to_anchor=(1.015, 0.03), fontsize=10, frameon=False, borderaxespad=0.0)

    # Panel B
    table_c = panel_tables["Figure2C_pwsdel_retained_maternal_plot.tsv"]
    add_shared_axis_style(ax_c, display_spec, "Mean methylation", show_x=False)
    ax_c.set_ylim(-0.02, 1.04)
    add_major_landmarks(ax_c, display_spec, snrpn_pos)
    plot_series(ax_c, table_c["window_mid"], table_c["control_maternal_mean_smoothed"], MATERNAL, "control maternal reference", lw=1.6, ls="--")
    plot_series(ax_c, table_c["window_mid"], table_c["pwsdel_retained_maternal_mean_smoothed"], MATERNAL, "PWS-DEL retained maternal", lw=2.2)
    ax_c.set_title(stats["B"].title, loc="left", fontsize=15, fontweight="bold", pad=8)
    add_badge(ax_c, stats["B"].badge, MATERNAL)
    add_stats_box(ax_c, stats["B"].stats_box)
    ax_c.legend(loc="lower left", bbox_to_anchor=(1.015, 0.03), fontsize=10, frameon=False, borderaxespad=0.0)

    # Panel C
    table_d = panel_tables["Figure2D_asdel_retained_paternal_plot.tsv"]
    add_shared_axis_style(ax_d, display_spec, "Mean methylation", show_x=False)
    ax_d.set_ylim(-0.02, 1.04)
    add_major_landmarks(ax_d, display_spec, snrpn_pos)
    plot_series(ax_d, table_d["window_mid"], table_d["control_paternal_mean_smoothed"], PATERNAL, "control paternal reference", lw=1.6, ls="--")
    plot_series(ax_d, table_d["window_mid"], table_d["asdel_retained_paternal_mean_smoothed"], PATERNAL, "AS-DEL retained paternal", lw=2.2)
    ax_d.set_title(stats["C"].title, loc="left", fontsize=15, fontweight="bold", pad=8)
    add_badge(ax_d, stats["C"].badge, PATERNAL)
    add_stats_box(ax_d, stats["C"].stats_box)
    ax_d.legend(loc="lower left", bbox_to_anchor=(1.015, 0.03), fontsize=10, frameon=False, borderaxespad=0.0)

    # Panel D
    table_e = panel_tables["Figure2E_pwsmupd_maternal_plot.tsv"]
    add_shared_axis_style(ax_e, display_spec, "Mean methylation", show_x=False)
    ax_e.set_ylim(-0.02, 1.04)
    add_major_landmarks(ax_e, display_spec, snrpn_pos)
    plot_series(ax_e, table_e["window_mid"], table_e["control_maternal_mean_smoothed"], MATERNAL, "control maternal reference", lw=1.5, ls="--")
    plot_series(ax_e, table_e["window_mid"], table_e["mean_meth_haplotype1_smoothed"], MUPD_1, "004P haplotype 1", lw=2.0)
    plot_series(ax_e, table_e["window_mid"], table_e["mean_meth_haplotype2_smoothed"], MUPD_2, "004P haplotype 2", lw=2.0)
    ax_e.set_title(stats["D"].title, loc="left", fontsize=15, fontweight="bold", pad=8)
    add_badge(ax_e, stats["D"].badge, MUPD_2)
    add_stats_box(ax_e, stats["D"].stats_box)
    ax_e.legend(loc="lower left", bbox_to_anchor=(1.015, 0.03), fontsize=10, frameon=False, borderaxespad=0.0)

    # Panel E
    table_f = panel_tables["Figure2F_reciprocal_delta_architecture_plot.tsv"]
    add_shared_axis_style(ax_f, display_spec, "Delta methylation\n(PWS-DEL maternal - AS-DEL paternal)", show_x=True)
    ax_f.set_ylim(-0.33, 0.18)
    add_major_landmarks(ax_f, display_spec, snrpn_pos)
    ax_f.axhline(0.0, color=DARK_GREY, lw=1.0, alpha=0.8)
    delta = table_f["reciprocal_delta_smoothed"]
    ax_f.fill_between(table_f["window_mid"], 0, delta, where=delta >= 0, color=MATERNAL, alpha=0.16, lw=0)
    ax_f.fill_between(table_f["window_mid"], 0, delta, where=delta < 0, color=PATERNAL, alpha=0.18, lw=0)
    plot_series(ax_f, table_f["window_mid"], delta, DARK_GREY, "PWS-DEL maternal minus AS-DEL paternal", lw=2.1)

    for start, end, _label in boundary_highlights(boundaries, gene_features):
        ax_f.axvspan(start, end, color="#ffd9a8", alpha=0.22, lw=0)

    ax_f.annotate(
        "shared methylation scaffold",
        xy=(24_250_000, 0.01),
        xytext=(24_250_000, 0.12),
        ha="center",
        va="bottom",
        fontsize=10,
        color="#333333",
        arrowprops=dict(arrowstyle="-", color="#666666", lw=0.8),
    )
    ax_f.set_title(stats["E"].title, loc="left", fontsize=15, fontweight="bold", pad=8)
    add_badge(ax_f, stats["E"].badge, PATERNAL)
    add_stats_box(ax_f, stats["E"].stats_box)
    ax_f.legend(loc="lower left", bbox_to_anchor=(1.015, 0.03), fontsize=10, frameon=False, borderaxespad=0.0)
    ax_f.set_xlabel("chr15 coordinate, T2T-CHM13v2.0 (Mb)", fontsize=14)

    fig.subplots_adjust(top=0.985, bottom=0.06, left=0.08, right=0.80)

    figdir = mkdir(outdir / "figures")
    fig.savefig(figdir / "Figure2_reciprocal_cis_architecture_improved.png", dpi=400, bbox_inches="tight")
    fig.savefig(figdir / "Figure2_reciprocal_cis_architecture_improved.pdf", bbox_inches="tight")
    fig.savefig(figdir / "Figure2_reciprocal_cis_architecture_improved.svg", bbox_inches="tight")
    plt.close(fig)


def report_text(
    outdir: Path,
    report_path: Path,
    input_files: List[str],
    output_files: List[str],
    control_ref: pd.DataFrame,
    corr: pd.DataFrame,
    reciprocal_delta: pd.DataFrame,
    boundaries: pd.DataFrame,
    smooth_window: int,
) -> str:
    timestamp = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %z")
    control_stats = pearson_summary(
        control_ref["control_maternal_mean"],
        control_ref["control_paternal_mean"],
    )
    pws = corr[corr["group"] == "PWS_DEL"].copy()
    asdel = corr[corr["group"] == "AS_DEL"].copy()
    upd = corr[corr["group"] == "PWS_mUPD"].copy().sort_values("comparison")
    reciprocal_nonnull = reciprocal_delta.dropna(
        subset=["pwsdel_retained_maternal_mean", "asdel_retained_paternal_mean"]
    )
    reciprocal_stats = pearson_summary(
        reciprocal_nonnull["pwsdel_retained_maternal_mean"],
        reciprocal_nonnull["asdel_retained_paternal_mean"],
    )

    rel = lambda path: path.relative_to(report_path.parent.parent).as_posix()
    lines = [
        "# Figure 2 Report: Reciprocal cis-methylation architecture across the 15q11-q13 imprinted domain",
        "",
        f"Generated: {timestamp}",
        "",
        "## 1. Purpose of Figure 2",
        "Figure 2 tests whether reciprocal deletion genomes reconstruct the parental cis-methylation architecture across the T2T chr15 imprinted domain: controls provide the biparental reference, PWS-DEL exposes the retained maternal-only profile, AS-DEL exposes the retained paternal-only profile, and PWS-mUPD confirms duplicated maternal identity.",
        "",
        "## 2. Input files",
    ]
    lines.extend([f"- `{path}`" for path in input_files])
    lines.extend(
        [
            "",
            "## 3. Coordinate system",
            "- Reference: T2T-CHM13v2.0",
            "- Analysis domain: `chr15:17,600,000-28,000,000`",
            "- Display domain in the improved figure: `chr15:18,000,000-28,000,000`",
            "- Major landmarks retained in the main figure: BP1, BP2, BP3, and PWS-IC/SNRPN",
            "",
            "## 4. Window size and smoothing parameters",
            f"- Base window size: 1 kb",
            f"- Smoothed display tracks: centered rolling median, {smooth_window} windows",
            f"- Approximate smoothing span: {smooth_window} kb",
            "",
            "## 5. Number of informative windows per panel",
            f"- Shared annotation track: structural and gene context only, no methylation windows",
            f"- Panel A (Controls): {control_stats['n_windows']:,} informative windows with both control maternal and control paternal means",
            f"- Panel B (PWS-DEL): {int(pws['n_windows'].sum()):,} within-footprint sample-windows across 5 deletion genomes",
            f"- Panel C (AS-DEL): {int(asdel['n_windows'].sum()):,} within-footprint sample-windows across 3 deletion genomes",
            f"- Panel D (PWS-mUPD): {int(upd.iloc[0]['n_windows']):,} windows for haplotype 1 and {int(upd.iloc[1]['n_windows']):,} windows for haplotype 2",
            f"- Panel E (Delta architecture): {reciprocal_stats['n_windows']:,} shared windows with both PWS-DEL retained maternal and AS-DEL retained paternal means",
            "",
            "## 6. Correlation statistics",
            f"- Controls: `r = {control_stats['pearson_r']:.3f}` across {control_stats['n_windows']:,} informative windows",
            f"- PWS-DEL retained maternal versus control maternal: mean `r = {pws['pearson_r'].mean():.3f}`, range `{pws['pearson_r'].min():.3f}-{pws['pearson_r'].max():.3f}`",
            f"- AS-DEL retained paternal versus control paternal: mean `r = {asdel['pearson_r'].mean():.3f}`, range `{asdel['pearson_r'].min():.3f}-{asdel['pearson_r'].max():.3f}`",
            f"- PWS-mUPD haplotype 1 versus control maternal: `r = {upd.iloc[0]['pearson_r']:.3f}`",
            f"- PWS-mUPD haplotype 2 versus control maternal: `r = {upd.iloc[1]['pearson_r']:.3f}`",
            f"- Reciprocal overlay (PWS-DEL retained maternal versus AS-DEL retained paternal): `r = {reciprocal_stats['pearson_r']:.3f}` across {reciprocal_stats['n_windows']:,} shared windows",
            "",
            "## 7. Interpretation of each panel",
            "- Shared annotation track: provides T2T structural and gene context for BP1, BP2, BP3, the PWS imprinting center, and broader 15q11-q13 genes.",
            "- Panel A: controls provide the biparental reference and show the reciprocal maternal and paternal methylation architecture used throughout the rest of the figure.",
            "- Panel B: PWS-DEL retained methylation closely follows the control maternal reference, supporting exposure of a maternal-only architecture by paternal deletion.",
            "- Panel C: AS-DEL retained methylation closely follows the control paternal reference, supporting exposure of a paternal-only architecture by maternal deletion.",
            "- Panel D: both 004P haplotypes track the control maternal reference, consistent with copy-neutral duplicated maternal architecture rather than biparental inheritance.",
            "- Panel E: the delta track is near zero across broad intervals, indicating a shared methylation scaffold, but shows localized parental divergence centered on the SNRPN/SNHG14 and SNORD116 interval.",
            "",
            "## 8. BP1/BP2/BP3 are structural landmarks, not methylation-boundary calls",
            "BP1, BP2, and BP3 are shown only as T2T structural breakpoint-cluster landmarks. They are not interpreted here as methylation boundaries. The localized methylation divergences highlighted in Panel E are data-driven parental differences within the imprinted gene domain, not calls anchored to BP1/BP2/BP3.",
            "",
            "## 9. Main conclusion",
            "Reciprocal PWS and AS deletion genomes reconstruct the parental cis-methylation architecture across the chr15 imprinted domain. PWS-DEL reveals the retained maternal scaffold, AS-DEL reveals the retained paternal scaffold, PWS-mUPD independently confirms duplicated maternal identity, and controls supply the biparental reference needed to interpret localized parental divergence near SNRPN/SNHG14 and SNORD116.",
            "",
            "## 10. Output file list",
        ]
    )
    lines.extend([f"- `{path}`" for path in output_files])
    return "\n".join(lines) + "\n"


def write_panel_tables(outdir: Path, annotation_track: pd.DataFrame, panel_tables: Dict[str, pd.DataFrame]) -> List[str]:
    tbldir = mkdir(outdir / "tables")
    outputs = []
    annotation_name = "Figure2_shared_annotation_track.tsv"
    annotation_track.to_csv(tbldir / annotation_name, sep="\t", index=False, na_rep="NaN")
    outputs.append(f"tables/{annotation_name}")
    for name, table in panel_tables.items():
        table.to_csv(tbldir / name, sep="\t", index=False, na_rep="NaN")
        outputs.append(f"tables/{name}")
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description="Create the improved Figure 2 from existing Phase 2 outputs")
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR, type=Path)
    parser.add_argument("--display-start", default=DISPLAY_START, type=int)
    parser.add_argument("--display-end", default=DISPLAY_END, type=int)
    parser.add_argument("--smooth-window", default=31, type=int)
    parser.add_argument("--gtf", default=DEFAULT_GTF, type=Path)
    parser.add_argument("--icr-bed", default=DEFAULT_ICR_BED, type=Path)
    args = parser.parse_args()

    outdir = args.outdir
    tbldir = outdir / "tables"
    reports_dir = mkdir(outdir / "reports")

    log("Loading existing Phase 2 tables")
    run_summary = read_tsv(outdir / "phase2_run_summary.tsv")
    window_size = int(run_summary.iloc[0]["window_size"])
    analysis_spec = WindowSpec(
        chrom=str(run_summary.iloc[0]["domain"]).split(":")[0],
        start=int(str(run_summary.iloc[0]["domain"]).split(":")[1].split("-")[0]),
        end=int(str(run_summary.iloc[0]["domain"]).split(":")[1].split("-")[1]),
        size=window_size,
    )
    display_spec = WindowSpec(analysis_spec.chrom, args.display_start, args.display_end, DISPLAY_WINDOW)

    control_ref = read_tsv(tbldir / "Phase2_control_reference_architecture.tsv")
    retained = read_tsv(tbldir / "Phase2_retained_haplotype_profiles.tsv")
    all_windows = read_tsv(tbldir / "Phase2_all_samples_haplotype_1kb_methylation.tsv.gz")
    corr = read_tsv(tbldir / "Phase2_retained_and_mUPD_correlations.tsv")
    reciprocal_delta = read_tsv(tbldir / "Phase2_reciprocal_delta_profile.tsv")
    boundaries = read_tsv(tbldir / "Phase2_reciprocal_delta_boundary_candidates.tsv")
    gene_features = read_tsv(tbldir / "Phase2_gene_track_features.tsv")

    annotation_track = build_annotation_track_table(gene_features, args.gtf, args.icr_bed, display_spec)
    panel_tables = build_panel_tables(
        display_spec,
        control_ref,
        retained,
        all_windows,
        reciprocal_delta,
        args.smooth_window,
    )
    stats = panel_stats(control_ref, corr, reciprocal_delta)

    log("Writing panel plotting tables")
    output_table_paths = write_panel_tables(outdir, annotation_track, panel_tables)

    log("Rendering improved Figure 2")
    render_figure(
        outdir,
        display_spec,
        annotation_track,
        panel_tables,
        boundaries,
        gene_features,
        stats,
    )

    report_path = reports_dir / "Figure2_reciprocal_cis_architecture_report.md"
    input_files = [
        "phase2_run_summary.tsv",
        "tables/Phase2_control_reference_architecture.tsv",
        "tables/Phase2_retained_haplotype_profiles.tsv",
        "tables/Phase2_all_samples_haplotype_1kb_methylation.tsv.gz",
        "tables/Phase2_retained_and_mUPD_correlations.tsv",
        "tables/Phase2_reciprocal_delta_profile.tsv",
        "tables/Phase2_reciprocal_delta_boundary_candidates.tsv",
        "tables/Phase2_gene_track_features.tsv",
        str(args.gtf),
        str(args.icr_bed),
        "scripts/paper_vf/create_figure2_reciprocal_cis_architecture_improved.py",
    ]
    output_files = [
        "figures/Figure2_reciprocal_cis_architecture_improved.png",
        "figures/Figure2_reciprocal_cis_architecture_improved.pdf",
        "figures/Figure2_reciprocal_cis_architecture_improved.svg",
        "reports/Figure2_reciprocal_cis_architecture_report.md",
    ] + output_table_paths
    report = report_text(
        outdir,
        report_path,
        input_files,
        output_files,
        control_ref,
        corr,
        reciprocal_delta,
        boundaries,
        args.smooth_window,
    )
    report_path.write_text(report)
    log("Done")


if __name__ == "__main__":
    main()
