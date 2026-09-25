#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle


# ============================== CONFIGURATION ==============================
PROJECT_ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = PROJECT_ROOT / "results" / "analysis"
CIS_DIR = ANALYSIS_DIR / "03_cis_architecture"
INPUTS = {
    "windows": CIS_DIR / "parent_associated_windows.tsv.gz",
    "window_participants": CIS_DIR / "parent_window_participant_values.tsv.gz",
    "focal": CIS_DIR / "focal_intervals.tsv",
    "regions": CIS_DIR / "regional_parent_contrasts.tsv",
    "regions_by_estimator": CIS_DIR / "regional_parent_contrasts_by_estimator.tsv",
    "participants": CIS_DIR / "regional_participant_values.tsv.gz",
    "asm": CIS_DIR / "regional_phase_invariant_asm.tsv",
    "missingness": CIS_DIR / "participant_missingness.tsv",
    "core": ANALYSIS_DIR / "01_evidence_matrix" / "common_reciprocal_cn1_core.tsv",
    "annotation": ANALYSIS_DIR / "annotation_genes.tsv",
}
OUTPUT_DIR = PROJECT_ROOT / "results" / "07_figures" / "figure_2"
MAIN_STEM = "Figure2"
SUPP_PROFILE_STEM = "Supplementary_Figure_S2A_retained_copy_profiles"
SUPP_DEPTH_STEM = "Supplementary_Figure_S2B_depth_robustness"
FORMATS = ("png", "pdf", "svg")
PNG_DPI = 600
PAD_INCHES = 0.04

TITLE = "Reciprocal deletions reveal focal parent-associated methylation divergence across 15q11–q13"
PANEL_C_TITLE = "Phase-invariant allelic methylation on intact chromosome 15"
FIGURE_WIDTH_IN = 7.0
MAX_EXPORT_WIDTH_IN = 7.09
FIGURE_HEIGHT_IN = 7.7
BASE_FONT = 7.0
INSET_RATIO = 2.0
A_YLIM_FACTOR = 1.6
FOCAL_LABEL_MERGE_MB = 0.02
JITTER = 0.07
JITTER_SEED = 11
NONEVALUABLE_MIN_BP = 5_000

MATERNAL = "#C8472B"
PATERNAL = "#2F6DB0"
TRACE = "#6A2C91"
INCONCLUSIVE = "#8A8A8A"
RAW = "#CFCFCF"
INK = "#222222"
MUTED = "#6B6B6B"
CONTROL = "#A86F0C"
DIGEORGE = "#0E8A74"
ESTIMATOR_STYLE = {
    "full_depth": ("Full depth", "#222222", "o"),
    "common_depth_downsampled": ("Common-depth downsampled", "#6A2C91", "s"),
    "capped_effective_coverage": ("Capped effective coverage", "#A86F0C", "D"),
    "shared_cpg": ("Shared CpGs", "#0E8A74", "^"),
}
REGION_SHORT = {
    "MAGEL2/NDN": "MAGEL2/NDN",
    "PWS/AS imprinting centre": "IC",
    "SNRPN/SNHG14": "SNRPN",
    "SNORD116": "SNORD116",
    "UBE3A": "UBE3A",
    "GABRB3/GABA receptor cluster": "GABA-R cluster",
    "OCA2 downstream control": "OCA2",
}
IC_REGION = "PWS/AS imprinting centre"
# ===========================================================================

REQUIRED_COLUMNS = {
    "windows": {
        "start", "end", "mid", "in_common_cn1", "evaluable", "segment_id", "delta_beta",
        "delta_rolling_median", "rolling_ci_low", "rolling_ci_high", "pws_n", "as_n",
        "mean_beta_pws_maternal_retained", "mean_beta_as_paternal_retained",
        "diploid_combined_mean_beta_descriptive", "missing_fraction",
        "delta_beta_common_depth_downsampled", "delta_beta_shared_cpg",
    },
    "window_participants": {"window_id", "start", "sample_id", "estimator", "beta", "mechanism"},
    "focal": {"start", "end", "direction", "reproducible"},
    "regions": {
        "display_order", "region_id", "pws_n", "as_n", "delta_beta", "ci_low", "ci_high",
        "status", "coverage_sensitivity", "shared_cpg_sensitivity", "robust",
    },
    "regions_by_estimator": {"region_id", "estimator", "delta_beta", "ci_low", "ci_high", "ci_width", "pws_n", "as_n"},
    "participants": {"analysis", "region_id", "sample_id", "cohort", "estimator", "value"},
    "asm": {"region_id", "cohort", "estimator", "n_participants", "mean_absolute_asm", "ci_low", "ci_high"},
    "missingness": {"sample_id", "mechanism", "estimator", "missing_fraction", "median_retained_copy_depth"},
    "core": {"start", "end"},
    "annotation": {"label", "start", "end"},
}


def load_inputs() -> dict[str, pd.DataFrame]:
    tables: dict[str, pd.DataFrame] = {}
    for name, path in INPUTS.items():
        if not path.is_file():
            raise FileNotFoundError(f"Missing Figure 2 input ({name}): {path}. Run scripts/analysis/run_analysis.py first.")
        table = pd.read_csv(path, sep="\t", low_memory=False)
        missing = REQUIRED_COLUMNS[name] - set(table.columns)
        if missing:
            raise ValueError(f"{path.name} lacks columns: {sorted(missing)}")
        tables[name] = table
    for name in ("windows", "regions", "asm", "regions_by_estimator"):
        if tables[name].empty:
            raise ValueError(f"{INPUTS[name].name} is empty")
    return tables


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": BASE_FONT,
            "axes.titlesize": BASE_FONT + 1,
            "axes.labelsize": BASE_FONT,
            "xtick.labelsize": BASE_FONT - 0.5,
            "ytick.labelsize": BASE_FONT - 0.5,
            "legend.fontsize": BASE_FONT - 0.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "axes.edgecolor": MUTED,
            "axes.labelcolor": INK,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "text.color": INK,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
        }
    )


def save(fig: plt.Figure, stem: str) -> list[Path]:
    fig.canvas.draw()
    extent = fig.get_tightbbox(fig.canvas.get_renderer())
    width = extent.width + 2 * PAD_INCHES
    if width > MAX_EXPORT_WIDTH_IN + 1e-6:
        raise RuntimeError(f"{stem} exports at {width:.2f} in, wider than {MAX_EXPORT_WIDTH_IN} in")
    outdir = OUTPUT_DIR / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    for fmt in FORMATS:
        path = outdir / f"{stem}.{fmt}"
        fig.savefig(path, dpi=PNG_DPI if fmt == "png" else None, bbox_inches="tight", pad_inches=PAD_INCHES)
        paths.append(path)
    plt.close(fig)
    return paths


def panel_label(ax: plt.Axes, letter: str, title: str, x: float = -0.07) -> None:
    ax.annotate(letter, (x, 1.04), xycoords="axes fraction", fontsize=BASE_FONT + 3, fontweight="bold",
                va="bottom", ha="left")
    ax.annotate(title, (x, 1.04), xycoords="axes fraction", xytext=(12, 0), textcoords="offset points",
                fontsize=BASE_FONT + 1, va="bottom", ha="left")


def mb(value: float) -> float:
    return value / 1e6


def segment_runs(frame: pd.DataFrame, value_col: str, segment_col: str):
    valid = frame[frame[value_col].notna() & frame[segment_col].ge(0)]
    for _, run in valid.groupby(segment_col, sort=True):
        yield run


def nonevaluable_spans(windows: pd.DataFrame) -> list[tuple[int, int]]:
    spans: list[tuple[int, int]] = []
    current: list[int] | None = None
    for row in windows.itertuples():
        if not row.evaluable:
            if current is None:
                current = [int(row.start), int(row.end)]
            else:
                current[1] = int(row.end)
        elif current is not None:
            spans.append(tuple(current))
            current = None
    if current is not None:
        spans.append(tuple(current))
    return [span for span in spans if span[1] - span[0] >= NONEVALUABLE_MIN_BP]


def direction_color(status: str) -> str:
    return {"maternal_retained_higher": MATERNAL, "paternal_retained_higher": PATERNAL}.get(status, INCONCLUSIVE)


def draw_panel_a(ax: plt.Axes, track: plt.Axes, tables: dict[str, pd.DataFrame], common: tuple[int, int]) -> None:
    windows = tables["windows"]
    windows = windows[windows["in_common_cn1"].astype(bool)].sort_values("start").reset_index(drop=True)
    focal = tables["focal"]
    focal = focal[focal["reproducible"].astype(bool)] if not focal.empty else focal
    smooth_extent = np.nanmax(
        np.abs(windows[["delta_rolling_median", "rolling_ci_low", "rolling_ci_high"]].to_numpy(float))
    )
    limit = float(np.clip(max(A_YLIM_FACTOR * smooth_extent, 0.3), 0.3, 1.0))
    limit = np.ceil(limit * 10) / 10
    groups: list[list] = []
    for record in focal.sort_values("start").itertuples():
        color = direction_color(record.direction)
        ax.axvspan(mb(record.start), mb(record.end), color=color, alpha=0.13, lw=0, zorder=0)
        if groups and mb(record.start) - groups[-1][1] < FOCAL_LABEL_MERGE_MB and groups[-1][3] == color:
            groups[-1][1] = mb(record.end)
            groups[-1][2].append(record.focal_id)
        else:
            groups.append([mb(record.start), mb(record.end), [record.focal_id], color])
    for start, end, labels, color in groups:
        ax.text(end + 0.02, limit * 0.93, "\u2013".join([labels[0], labels[-1]]) if len(labels) > 1 else labels[0],
                ha="left", va="top", fontsize=BASE_FONT - 0.5, color=color, fontweight="bold")
    for start, end in nonevaluable_spans(windows):
        ax.add_patch(Rectangle((mb(start), -limit), mb(end) - mb(start), limit * 0.035, color="#BDBDBD", lw=0, zorder=1))
    ax.axhline(0, color=MUTED, lw=0.6, zorder=1)
    for run in segment_runs(windows, "delta_beta", "segment_id"):
        x = mb(run["mid"].to_numpy())
        y = np.clip(run["delta_beta"].to_numpy(), -limit, limit)
        if len(run) == 1:
            ax.plot(x, y, ".", color=RAW, ms=1.5, zorder=2)
        else:
            ax.plot(x, y, color=RAW, lw=0.45, zorder=2, solid_joinstyle="round")
    clipped = windows[windows["delta_beta"].abs() > limit]
    for sign, marker in ((1, "^"), (-1, "v")):
        subset = clipped[np.sign(clipped["delta_beta"]) == sign]
        if not subset.empty:
            ax.plot(mb(subset["mid"]), np.full(len(subset), sign * limit * 0.975), marker=marker, ls="none",
                    ms=2.2, color=MUTED, zorder=3)
    for run in segment_runs(windows, "delta_rolling_median", "segment_id"):
        x = mb(run["mid"].to_numpy())
        ribbon = run["rolling_ci_low"].notna() & run["rolling_ci_high"].notna()
        if ribbon.any():
            ax.fill_between(x, run["rolling_ci_low"], run["rolling_ci_high"], where=ribbon.to_numpy(),
                            color=TRACE, alpha=0.2, lw=0, zorder=3)
        ax.plot(x, run["delta_rolling_median"], color=TRACE, lw=1.1, zorder=4, solid_capstyle="round")
    annotation = tables["annotation"]
    ic = annotation[annotation["label"].eq("IC")]
    if not ic.empty:
        ic_mid = mb((int(ic["start"].iloc[0]) + int(ic["end"].iloc[0])) / 2)
        ax.axvline(ic_mid, color=INK, lw=0.6, ls=(0, (2, 2)), zorder=1)
        ax.text(ic_mid, -limit * 0.97, " IC", ha="left", va="bottom", fontsize=BASE_FONT - 0.5, color=INK)
    ax.set_xlim(mb(common[0]), mb(common[1]))
    ax.set_ylim(-limit, limit)
    ax.set_ylabel("Δβ (PWS − AS)")
    ax.tick_params(labelbottom=False)
    ax.text(1.005, 0.97, "maternal-\nretained\nhigher", transform=ax.transAxes, color=MATERNAL, fontsize=BASE_FONT - 1,
            va="top", ha="left")
    ax.text(1.005, 0.03, "paternal-\nretained\nhigher", transform=ax.transAxes, color=PATERNAL, fontsize=BASE_FONT - 1,
            va="bottom", ha="left")
    legend = [
        Line2D([], [], color=RAW, lw=1.2, label="1-kb \u0394\u03b2"),
        Line2D([], [], color=TRACE, lw=1.2, label="21-kb median"),
        Patch(color=TRACE, alpha=0.2, label="participant bootstrap 95% CI"),
        Patch(color="#BDBDBD", label="not evaluable"),
        Line2D([], [], marker="^", ls="none", color=MUTED, ms=3, label="beyond axis"),
    ]
    ax.legend(handles=legend, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=5, fontsize=BASE_FONT - 1,
              handlelength=1.4, columnspacing=1.0, borderaxespad=0.2)
    draw_track(track, tables, common)


def repel(centres: np.ndarray, widths: np.ndarray, low: float, high: float) -> np.ndarray:
    positions = centres.astype(float).copy()
    for _ in range(500):
        moved = False
        for index in range(1, len(positions)):
            need = (widths[index] + widths[index - 1]) / 2
            gap = positions[index] - positions[index - 1]
            if gap < need:
                shift = (need - gap) / 2
                positions[index - 1] -= shift
                positions[index] += shift
                moved = True
        positions = np.clip(positions, low + widths / 2, high - widths / 2)
        if not moved:
            break
    return positions


def draw_track(ax: plt.Axes, tables: dict[str, pd.DataFrame], common: tuple[int, int]) -> None:
    annotation = tables["annotation"]
    regions = tables["regions"].sort_values("display_order")
    low, high = mb(common[0]), mb(common[1])
    items: list[tuple[float, str, str]] = []
    covered = []
    for record in regions.itertuples():
        if not np.isfinite(record.analysed_start):
            continue
        start, end = mb(record.analysed_start), mb(record.analysed_end)
        covered.append((record.analysed_start, record.analysed_end))
        ax.add_patch(Rectangle((start, 0.66), max(end - start, 0.003), 0.22, color="#4A4A4A", lw=0))
        items.append(((start + end) / 2, REGION_SHORT.get(record.region_id, record.region_id), INK))
    for record in annotation.itertuples():
        if record.label == "IC" or record.end < common[0] or record.start > common[1]:
            continue
        if any(record.start < c_end and record.end > c_start for c_start, c_end in covered):
            continue
        start, end = mb(max(record.start, common[0])), mb(min(record.end, common[1]))
        ax.add_patch(Rectangle((start, 0.70), end - start, 0.14, color="#C9C9C9", lw=0))
        items.append(((start + end) / 2, record.label, MUTED))
    items.sort()
    bbox = ax.get_position()
    axis_inches = bbox.width * ax.figure.get_figwidth()
    char_width = (BASE_FONT - 1) * 0.6 / 72 * (high - low) / axis_inches
    centres = np.array([item[0] for item in items])
    widths = np.array([len(item[1]) * char_width + 3 * char_width for item in items])
    positions = repel(centres, widths, low, high)
    for (centre, label, color), position in zip(items, positions):
        ax.text(position, 0.08, label, ha="center", va="bottom", fontsize=BASE_FONT - 1, color=color)
        if abs(position - centre) > char_width:
            ax.plot([centre, centre, position], [0.64, 0.5, 0.36], color="#9A9A9A", lw=0.4)
    ic = annotation[annotation["label"].eq("IC")]
    if not ic.empty:
        ic_mid = mb((int(ic["start"].iloc[0]) + int(ic["end"].iloc[0])) / 2)
        ax.plot([ic_mid], [0.96], marker="v", color=INK, ms=3.2, clip_on=False)
    ax.set_xlim(low, high)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("T2T-CHM13 chr15 coordinate (Mb)")
    ax.text(-0.01, 0.77, "regions\n/ genes", transform=ax.transAxes, ha="right", va="center",
            fontsize=BASE_FONT - 1, color=MUTED)


def split_axes(fig: plt.Figure, spec, use_inset: bool, width_ratio: float = 0.24):
    if not use_inset:
        return fig.add_subplot(spec), None
    grid = GridSpecFromSubplotSpec(1, 2, subplot_spec=spec, width_ratios=[1 - width_ratio, width_ratio], wspace=0.08)
    main = fig.add_subplot(grid[0])
    inset = fig.add_subplot(grid[1], sharey=main)
    inset.tick_params(labelleft=False, left=False)
    inset.spines["left"].set_visible(False)
    inset.set_facecolor("#F6F6F6")
    return main, inset


def needs_inset(values: pd.Series, others: pd.Series) -> bool:
    other_extent = np.nanmax(np.abs(others.to_numpy(float))) if len(others) else np.nan
    extent = np.nanmax(np.abs(values.to_numpy(float))) if len(values) else np.nan
    return bool(np.isfinite(extent) and np.isfinite(other_extent) and extent > INSET_RATIO * other_extent)


def padded_limits(values: np.ndarray, include_zero: bool = True, pad: float = 0.12) -> tuple[float, float]:
    values = values[np.isfinite(values)]
    if include_zero:
        values = np.append(values, 0.0)
    low, high = values.min(), values.max()
    span = max(high - low, 0.05)
    return low - pad * span, high + pad * span


def draw_panel_b(fig: plt.Figure, spec, tables: dict[str, pd.DataFrame]) -> tuple[plt.Axes, list[str]]:
    regions = tables["regions"].sort_values("display_order").reset_index(drop=True)
    order = regions["region_id"].tolist()
    ypos = {region: index for index, region in enumerate(order)}
    ic_rows = regions[regions["region_id"].eq(IC_REGION)]
    others = regions[~regions["region_id"].eq(IC_REGION)]
    use_inset = needs_inset(ic_rows[["ci_low", "ci_high", "delta_beta"]].stack(),
                            others[["ci_low", "ci_high", "delta_beta"]].stack())
    main, inset = split_axes(fig, spec, use_inset)
    for axis in filter(None, (main, inset)):
        axis.axvline(0, color=MUTED, lw=0.6)
    for record in regions.itertuples():
        axis = inset if (use_inset and record.region_id == IC_REGION) else main
        y = ypos[record.region_id]
        if record.status == "not_evaluable":
            main.text(0, y, "not evaluable", ha="center", va="center", fontsize=BASE_FONT - 1, color=MUTED)
            continue
        color = direction_color(record.status)
        robust_inputs = record.coverage_sensitivity == "pass" and record.shared_cpg_sensitivity == "pass"
        axis.plot([record.ci_low, record.ci_high], [y, y], color=color, lw=1.6, solid_capstyle="round", zorder=2)
        axis.plot(record.delta_beta, y, marker="o", ms=5.5, mec=color, mew=1.2,
                  mfc=color if robust_inputs else "white", zorder=3)
        if use_inset and record.region_id == IC_REGION:
            main.text(1.0, y, "→", transform=main.get_yaxis_transform(),
                      ha="right", va="center", fontsize=BASE_FONT, color=MUTED)
    main_values = others[["ci_low", "ci_high"]].to_numpy(float).ravel() if use_inset else regions[["ci_low", "ci_high"]].to_numpy(float).ravel()
    main.set_xlim(*padded_limits(main_values))
    if use_inset:
        inset_values = ic_rows[["ci_low", "ci_high"]].to_numpy(float).ravel()
        low, high = padded_limits(inset_values, include_zero=False, pad=0.6)
        inset.set_xlim(low, high)
        inset.xaxis.set_major_locator(plt.MaxNLocator(2))
        inset.set_title("IC scale", fontsize=BASE_FONT - 1, color=MUTED, pad=2)
    main.set_yticks(range(len(order)))
    main.set_yticklabels(order)
    main.set_ylim(len(order) - 0.5, -0.7)
    main.set_xlabel("Δβ (PWS maternal-retained − AS paternal-retained)")
    main.xaxis.set_major_locator(plt.MaxNLocator(4))
    right = inset if use_inset else main
    right.text(1.04, -0.75, "PWS/AS\nn", transform=right.get_yaxis_transform(), ha="left", va="bottom",
               fontsize=BASE_FONT - 1, color=MUTED)
    for record in regions.itertuples():
        flag = "" if record.status == "not_evaluable" or (
            record.coverage_sensitivity == "pass" and record.shared_cpg_sensitivity == "pass") else "†"
        right.text(1.04, ypos[record.region_id], f"{record.pws_n}/{record.as_n}{flag}",
                   transform=right.get_yaxis_transform(), ha="left", va="center", fontsize=BASE_FONT - 1, color=INK)
    legend = [
        Line2D([], [], marker="o", color=MATERNAL, ms=4.5, lw=1.4, label="maternal-retained higher"),
        Line2D([], [], marker="o", color=PATERNAL, ms=4.5, lw=1.4, label="paternal-retained higher"),
        Line2D([], [], marker="o", color=INCONCLUSIVE, ms=4.5, lw=1.4, label="inconclusive"),
        Line2D([], [], marker="o", color=MUTED, mfc="white", ms=4.5, lw=0, label="† fails depth or shared-CpG check"),
    ]
    main.legend(handles=legend, loc="upper left", bbox_to_anchor=(-0.02, -0.17), ncol=1, fontsize=BASE_FONT - 1,
                handlelength=1.3)
    return main, order


def draw_panel_c(fig: plt.Figure, spec, tables: dict[str, pd.DataFrame], order: list[str]) -> plt.Axes:
    asm = tables["asm"][tables["asm"]["estimator"].eq("full_depth")]
    participants = tables["participants"]
    participants = participants[
        participants["analysis"].eq("phase_invariant_asm") & participants["estimator"].eq("full_depth")
    ]
    ypos = {region: index for index, region in enumerate(order)}
    ic_values = participants.loc[participants["region_id"].eq(IC_REGION), "value"]
    other_values = participants.loc[~participants["region_id"].eq(IC_REGION), "value"]
    use_inset = needs_inset(ic_values, other_values)
    main, inset = split_axes(fig, spec, use_inset)
    rng = np.random.default_rng(JITTER_SEED)
    cohorts = (("Control", CONTROL, "o", -0.17, "Unaffected control"), ("DiGeorge", DIGEORGE, "^", 0.17, "DiGeorge"))
    totals = {}
    for cohort, color, marker, offset, _ in cohorts:
        subset = participants[participants["cohort"].eq(cohort)]
        totals[cohort] = subset["sample_id"].nunique()
        for region in order:
            axis = inset if (use_inset and region == IC_REGION) else main
            y = ypos[region] + offset
            values = subset.loc[subset["region_id"].eq(region), "value"].dropna().to_numpy()
            jitter = rng.uniform(-JITTER, JITTER, len(values))
            axis.scatter(values, y + jitter, s=11, marker=marker, facecolor="white", edgecolor=color, lw=0.8, zorder=3)
            summary = asm[asm["cohort"].eq(cohort) & asm["region_id"].eq(region)]
            if summary.empty or not np.isfinite(summary["mean_absolute_asm"].iloc[0]):
                continue
            record = summary.iloc[0]
            if np.isfinite(record["ci_low"]):
                axis.plot([record["ci_low"], record["ci_high"]], [y, y], color=color, lw=1.4, zorder=2,
                          solid_capstyle="round")
            axis.plot([record["mean_absolute_asm"]] * 2, [y - 0.13, y + 0.13], color=color, lw=2.0, zorder=4)
    for axis in filter(None, (main, inset)):
        axis.set_xlim(left=0)
    other = participants[~participants["region_id"].eq(IC_REGION)] if use_inset else participants
    upper = np.nanmax(other["value"].to_numpy(float))
    main.set_xlim(0, upper * 1.15)
    if use_inset:
        values = ic_values.dropna().to_numpy()
        span = max(values.max() - values.min(), 0.05)
        inset.set_xlim(values.min() - 0.6 * span, min(1.0, values.max() + 0.6 * span))
        inset.xaxis.set_major_locator(plt.MaxNLocator(2))
        inset.set_title("IC scale", fontsize=BASE_FONT - 1, color=MUTED, pad=2)
        main.text(1.0, ypos[IC_REGION], "→", transform=main.get_yaxis_transform(), ha="right", va="center",
                  fontsize=BASE_FONT, color=MUTED)
    main.set_yticks(range(len(order)))
    main.set_yticklabels([])
    main.set_ylim(len(order) - 0.5, -0.7)
    main.set_xlabel("regional mean |β$_{H1}$ − β$_{H2}$| per participant")
    main.xaxis.set_major_locator(plt.MaxNLocator(4))
    right = inset if use_inset else main
    right.text(1.04, -0.75, "Ctrl/DG\nn", transform=right.get_yaxis_transform(), ha="left", va="bottom",
               fontsize=BASE_FONT - 1, color=MUTED)
    for region in order:
        counts = [
            int(asm.loc[asm["cohort"].eq(cohort) & asm["region_id"].eq(region), "n_participants"].sum())
            for cohort, *_ in cohorts
        ]
        right.text(1.04, ypos[region], f"{counts[0]}/{counts[1]}", transform=right.get_yaxis_transform(),
                   ha="left", va="center", fontsize=BASE_FONT - 1, color=INK)
    legend = [
        Line2D([], [], marker="o", ls="none", mfc="white", mec=CONTROL, ms=4, label=f"Unaffected control (n={totals['Control']})"),
        Line2D([], [], marker="^", ls="none", mfc="white", mec=DIGEORGE, ms=4, label=f"DiGeorge (n={totals['DiGeorge']})"),
        Line2D([], [], color=MUTED, lw=2.0, marker="|", ms=0, label="group mean (equal weight)"),
        Line2D([], [], color=DIGEORGE, lw=1.4, label="95% participant bootstrap,\ndescriptive (only n≥3)"),
    ]
    main.legend(handles=legend, loc="upper left", bbox_to_anchor=(-0.02, -0.17), ncol=1, fontsize=BASE_FONT - 1,
                handlelength=1.3)
    return main


def render_main(tables: dict[str, pd.DataFrame], common: tuple[int, int]) -> list[Path]:
    fig = plt.figure(figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN))
    outer = GridSpec(4, 1, figure=fig, height_ratios=[2.0, 0.5, 1.05, 2.75], hspace=0.0,
                     left=0.235, right=0.87, top=0.9, bottom=0.16)
    ax_a = fig.add_subplot(outer[0])
    track = fig.add_subplot(outer[1], sharex=ax_a)
    lower = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[3], width_ratios=[1.0, 0.8], wspace=0.32)
    draw_panel_a(ax_a, track, tables, common)
    windows = tables["windows"]
    n_pws = int(windows["pws_n"].max())
    n_as = int(windows["as_n"].max())
    ax_b, order = draw_panel_b(fig, lower[0, 0], tables)
    ax_c = draw_panel_c(fig, lower[0, 1], tables, order)
    left = 0.012
    a_top = ax_a.get_position().y1
    b_top = ax_b.get_position().y1
    fig.text(left, a_top + 0.03, "a", fontsize=BASE_FONT + 3, fontweight="bold", va="bottom")
    fig.text(left + 0.028, a_top + 0.03,
             f"Direct parent-associated contrast in the common reciprocal CN=1 interval (PWS n={n_pws}, AS n={n_as})",
             fontsize=BASE_FONT + 1, va="bottom")
    fig.text(left, b_top + 0.035, "b", fontsize=BASE_FONT + 3, fontweight="bold", va="bottom")
    fig.text(left + 0.028, b_top + 0.035, "Prespecified regional parent-associated effects",
             fontsize=BASE_FONT + 1, va="bottom")
    c_left = ax_c.get_position().x0 - 0.035
    fig.text(c_left, b_top + 0.035, "c", fontsize=BASE_FONT + 3, fontweight="bold", va="bottom")
    fig.text(c_left + 0.028, b_top + 0.035, PANEL_C_TITLE.replace(" on intact", "\non intact"),
             fontsize=BASE_FONT + 1, va="bottom", linespacing=1.1)
    fig.suptitle(TITLE, x=left, y=0.985, ha="left", fontsize=BASE_FONT + 2, fontweight="bold")
    return save(fig, MAIN_STEM)


def heatmap(ax: plt.Axes, tables: dict[str, pd.DataFrame], common: tuple[int, int]) -> None:
    data = tables["window_participants"]
    data = data[data["estimator"].eq("full_depth")]
    wide = data.pivot(index="sample_id", columns="start", values="beta")
    mechanism = data.drop_duplicates("sample_id").set_index("sample_id")["mechanism"]
    rows = sorted(wide.index, key=lambda sample: (mechanism[sample] != "PWS-DEL", sample))
    wide = wide.loc[rows]
    cmap = plt.get_cmap("Purples").copy()
    cmap.set_bad("white")
    image = ax.imshow(np.ma.masked_invalid(wide.to_numpy(float)), aspect="auto", interpolation="nearest",
                      cmap=cmap, vmin=0, vmax=1,
                      extent=(mb(wide.columns.min()), mb(wide.columns.max() + 1000), len(rows) - 0.5, -0.5))
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{sample} ({'PWS mat.' if mechanism[sample] == 'PWS-DEL' else 'AS pat.'})" for sample in rows])
    boundary = sum(mechanism[sample] == "PWS-DEL" for sample in rows) - 0.5
    ax.axhline(boundary, color=INK, lw=0.8)
    bar = plt.colorbar(image, cax=ax.inset_axes([1.01, 0.0, 0.012, 1.0]))
    bar.set_label("retained-copy β", fontsize=BASE_FONT - 1)
    ax.set_xlim(mb(common[0]), mb(common[1]))


def profile_lines(ax: plt.Axes, windows: pd.DataFrame, column: str, color: str, label: str) -> None:
    first = True
    for run in segment_runs(windows.assign(seg=windows["segment_id"]), column, "seg"):
        ax.plot(mb(run["mid"]), run[column], color=color, lw=0.35, alpha=0.85, label=label if first else None)
        first = False


def diploid_lines(ax: plt.Axes, windows: pd.DataFrame, column: str) -> None:
    valid = windows[column].notna().to_numpy()
    starts, ends = windows["start"].to_numpy(), windows["end"].to_numpy()
    first = True
    run_id = np.cumsum(np.r_[True, (starts[1:] != ends[:-1]) | ~valid[1:] | ~valid[:-1]])
    for _, run in windows[valid].groupby(run_id[valid]):
        ax.plot(mb(run["mid"]), run[column], color=MUTED, lw=0.35,
                label="Control + DiGeorge combined-track mean β (descriptive)" if first else None)
        first = False


def render_profiles(tables: dict[str, pd.DataFrame], common: tuple[int, int]) -> list[Path]:
    windows = tables["windows"]
    windows = windows[windows["in_common_cn1"].astype(bool)].sort_values("start").reset_index(drop=True)
    fig, axes = plt.subplots(4, 1, figsize=(FIGURE_WIDTH_IN - 0.1, 7.6), sharex=True,
                             gridspec_kw={"height_ratios": [1.6, 1.1, 1.0, 0.8], "hspace": 0.28})
    heatmap(axes[0], tables, common)
    panel_label(axes[0], "a", "Participant-by-window retained-copy β (full depth; white = not evaluable)", x=-0.12)
    profile_lines(axes[1], windows, "mean_beta_pws_maternal_retained", MATERNAL, "PWS maternal-retained mean β")
    profile_lines(axes[1], windows, "mean_beta_as_paternal_retained", PATERNAL, "AS paternal-retained mean β")
    axes[1].set_ylim(0, 1)
    axes[1].set_ylabel("mean β")
    axes[1].legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, fontsize=BASE_FONT - 1, borderaxespad=0.1)
    panel_label(axes[1], "b", "Raw retained-copy group means", x=-0.12)
    diploid_lines(axes[2], windows, "diploid_combined_mean_beta_descriptive")
    axes[2].set_ylim(0, 1)
    axes[2].set_ylabel("mean β")
    axes[2].legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), fontsize=BASE_FONT - 1, borderaxespad=0.1)
    panel_label(axes[2], "c", "Pooled diploid combined-track mean (descriptive only)", x=-0.12)
    axes[3].step(mb(windows["mid"]), windows["pws_n"], where="mid", color=MATERNAL, lw=0.5, label="PWS evaluable n")
    axes[3].step(mb(windows["mid"]), windows["as_n"], where="mid", color=PATERNAL, lw=0.5, label="AS evaluable n")
    axes[3].set_ylabel("participants")
    axes[3].set_ylim(-0.3, max(windows["pws_n"].max(), windows["as_n"].max()) + 0.8)
    axes[3].legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, fontsize=BASE_FONT - 1, borderaxespad=0.1)
    axes[3].set_xlabel("T2T-CHM13 chr15 coordinate (Mb)")
    panel_label(axes[3], "d", "Per-window participant support", x=-0.12)
    fig.suptitle("Supplementary Figure S2A. Retained-copy methylation in the common reciprocal CN=1 interval",
                 x=0.02, y=0.96, ha="left", fontsize=BASE_FONT + 1.5, fontweight="bold")
    return save(fig, SUPP_PROFILE_STEM)


def render_depth(tables: dict[str, pd.DataFrame]) -> list[Path]:
    by_estimator = tables["regions_by_estimator"]
    regions = tables["regions"].sort_values("display_order")
    order = regions["region_id"].tolist()
    ypos = {region: index for index, region in enumerate(order)}
    estimators = list(ESTIMATOR_STYLE)
    offsets = np.linspace(-0.27, 0.27, len(estimators))
    fig = plt.figure(figsize=(FIGURE_WIDTH_IN - 0.5, 8.2))
    grid = GridSpec(3, 3, figure=fig, height_ratios=[1.25, 1.0, 1.0], hspace=0.62, wspace=0.5,
                    left=0.25, right=0.97, top=0.9, bottom=0.07)
    ax = fig.add_subplot(grid[0, :2])
    ax.axvline(0, color=MUTED, lw=0.6)
    for offset, estimator in zip(offsets, estimators):
        label, color, marker = ESTIMATOR_STYLE[estimator]
        subset = by_estimator[by_estimator["estimator"].eq(estimator)]
        for record in subset.itertuples():
            y = ypos[record.region_id] + offset
            ax.plot([record.ci_low, record.ci_high], [y, y], color=color, lw=1.0)
            ax.plot(record.delta_beta, y, marker=marker, color=color, ms=3.5, ls="none",
                    label=label if record.region_id == order[0] else None)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order)
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.set_xlabel("Δβ (PWS − AS) with participant-bootstrap 95% CI")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=BASE_FONT - 1)
    panel_label(ax, "a", "Regional effect agreement across depth handling", x=-0.45)
    ax = fig.add_subplot(grid[1, 0])
    full = by_estimator[by_estimator["estimator"].eq("full_depth")].set_index("region_id")
    ax.axvline(1, color=MUTED, lw=0.6)
    for offset, estimator in zip(offsets[1:], estimators[1:]):
        label, color, marker = ESTIMATOR_STYLE[estimator]
        subset = by_estimator[by_estimator["estimator"].eq(estimator)].set_index("region_id")
        ratio = subset["ci_width"] / full["ci_width"]
        ax.plot(ratio.values, [ypos[region] + offset for region in ratio.index], marker=marker, color=color, ls="none", ms=3.5)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order)
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.set_xlabel("CI width / full-depth CI width")
    panel_label(ax, "b", "Change in uncertainty", x=-0.9)
    ax = fig.add_subplot(grid[1, 1])
    for index, estimator in enumerate(estimators):
        subset = by_estimator[by_estimator["estimator"].eq(estimator)].set_index("region_id")
        for region in order:
            record = subset.loc[region]
            complete = record["pws_n"] >= full.loc[region, "pws_n"] and record["as_n"] >= full.loc[region, "as_n"]
            ax.text(index, ypos[region], f"{int(record['pws_n'])}/{int(record['as_n'])}", ha="center", va="center",
                    fontsize=BASE_FONT - 1, color=INK if complete else MATERNAL,
                    fontweight="normal" if complete else "bold")
    ax.set_xlim(-0.5, len(estimators) - 0.5)
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.set_xticks(range(len(estimators)))
    ax.set_xticklabels(["full", "down-\nsampled", "capped", "shared\nCpG"])
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    panel_label(ax, "c", "Participants retained (PWS/AS)", x=-0.05)
    ax = fig.add_subplot(grid[1, 2])
    asm = tables["asm"]
    for offset, estimator in zip(offsets, estimators):
        label, color, marker = ESTIMATOR_STYLE[estimator]
        subset = asm[asm["estimator"].eq(estimator) & asm["cohort"].eq("DiGeorge")]
        ax.plot(subset["mean_absolute_asm"], [ypos[region] + offset for region in subset["region_id"]],
                marker=marker, color=color, ls="none", ms=3.5)
    ax.set_xscale("log")
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("DiGeorge mean |β$_{H1}$ − β$_{H2}$| (log)")
    panel_label(ax, "d", "Phase-invariant ASM by estimator", x=-0.05)
    ax = fig.add_subplot(grid[2, :2])
    missing = tables["missingness"]
    samples = (
        missing.drop_duplicates("sample_id")
        .sort_values(["mechanism", "median_retained_copy_depth"], ascending=[False, True])["sample_id"].tolist()
    )
    for offset, estimator in zip(offsets, estimators):
        label, color, marker = ESTIMATOR_STYLE[estimator]
        subset = missing[missing["estimator"].eq(estimator)].set_index("sample_id").loc[samples]
        ax.plot(np.arange(len(samples)) + offset * 0.8, subset["missing_fraction"], marker=marker, color=color,
                ls="none", ms=3.5)
    depth = missing.drop_duplicates("sample_id").set_index("sample_id").loc[samples]
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels([f"{sample}\n{'PWS' if depth.loc[sample, 'mechanism'] == 'PWS-DEL' else 'AS'}\n{depth.loc[sample, 'median_retained_copy_depth']:.0f}×"
                        for sample in samples], fontsize=BASE_FONT - 1)
    ax.set_ylabel("fraction of CN=1 windows\nnot evaluable (label: depth)")
    ax.set_ylim(0, min(1.0, 1.12 * missing["missing_fraction"].max()))
    panel_label(ax, "e", "Missingness by participant and estimator", x=-0.14)
    windows = tables["windows"]
    windows = windows[windows["evaluable"].astype(bool)]
    ax = fig.add_subplot(grid[2, 2])
    pairs = windows[["delta_beta", "delta_beta_common_depth_downsampled"]].dropna()
    ax.hexbin(pairs["delta_beta"], pairs["delta_beta_common_depth_downsampled"], gridsize=40, cmap="Purples",
              mincnt=1, bins="log", linewidths=0)
    limit = np.nanpercentile(np.abs(pairs.to_numpy()), 99.8)
    ax.plot([-limit, limit], [-limit, limit], color=MUTED, lw=0.6, ls="--")
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    r = pairs.corr().iloc[0, 1]
    ax.text(0.04, 0.96, f"r = {r:.2f}\n{len(pairs):,} windows", transform=ax.transAxes, va="top", fontsize=BASE_FONT - 1)
    ax.set_xlabel("full-depth 1-kb Δβ")
    ax.set_ylabel("downsampled 1-kb Δβ")
    panel_label(ax, "f", "Window concordance", x=-0.3)
    fig.suptitle("Supplementary Figure S2B. Sequencing-depth and CpG-composition robustness", x=0.02, y=0.975,
                 ha="left", fontsize=BASE_FONT + 1.5, fontweight="bold")
    fig.text(0.02, 0.945, "Depth handling evaluates robustness only; it does not create biological evidence. "
             "Intervals resample participants, not windows.", fontsize=BASE_FONT - 0.5, color=MUTED)
    return save(fig, SUPP_DEPTH_STEM)


def write_render_log(paths: list[Path], tables: dict[str, pd.DataFrame], common: tuple[int, int]) -> Path:
    regions = tables["regions"].sort_values("display_order")
    rows = [{"item": "common_reciprocal_cn1_interval", "value": f"chr15:{common[0]}-{common[1]}"}]
    rows += [{"item": "output", "value": str(path)} for path in paths]
    for record in regions.itertuples():
        rows.append(
            {
                "item": f"panel_b:{record.region_id}",
                "value": f"{record.delta_beta:+.3f} [{record.ci_low:+.3f},{record.ci_high:+.3f}] "
                f"PWS/AS {record.pws_n}/{record.as_n} {record.status} coverage={record.coverage_sensitivity} "
                f"shared_cpg={record.shared_cpg_sensitivity}",
            }
        )
    path = OUTPUT_DIR / "reports" / "figure2_render_log.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def main() -> None:
    setup_style()
    tables = load_inputs()
    common = (int(tables["core"]["start"].iloc[0]), int(tables["core"]["end"].iloc[0]))
    windows = tables["windows"]
    if windows.loc[~windows["in_common_cn1"].astype(bool), "delta_beta"].notna().any():
        raise ValueError("Window table contains parental contrasts outside the common CN=1 interval")
    paths = render_main(tables, common)
    paths += render_profiles(tables, common)
    paths += render_depth(tables)
    paths.append(write_render_log(paths, tables, common))
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
