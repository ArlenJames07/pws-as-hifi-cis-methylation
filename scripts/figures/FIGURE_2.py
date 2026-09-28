#!/usr/bin/env python3
"""
Figure 2 -- Reciprocal parental methylation architecture across the chromosome 15 imprinted domain.

Built from the new analysis outputs (nothing is recomputed from raw data):
  a  retained-copy methylation per participant and 1-kb window in the common reciprocal CN=1
     interval (analysis 03), with the PWS-mUPD genome and the biparental mean below
     (analysis 01 window matrix)
  b  parental contrast at two scales on one axis: 1-kb CpG windows, PWS-DEL maternal-retained
     minus AS-DEL paternal-retained, 21-kb median with participant bootstrap (analysis 03), and
     repeat elements in 250-kb bins, maternal minus paternal with the envelope of all label
     permutations (scripts/duplicons 07); the largest methylation domain (duplicons 08) is shaded
  c  the domain tested in each genome: (inside - outside) minus the biparental median, for
     CpG windows (analysis 01) and repeat elements (duplicons 07); PWS-mUPD took no part in
     finding the domain
  d  prespecified regional contrasts with bootstrap intervals and sensitivity flags (analysis 03)
  e  phase-invariant allelic methylation |H1 - H2| on intact chromosomes (analysis 03)
Supplementary S2A (retained-copy profiles) and S2B (depth and CpG-composition robustness).

Inputs : results/analysis/{01_evidence_matrix,03_cis_architecture}/..., analysis/annotation_genes.tsv,
         results/08_duplicons/{repeat_methylation,followup}/... (panels b and c; skipped if absent)
Outputs: results/07_figures/figure_2/{figures,tables,reports}/

Usage: python3 scripts/figures/FIGURE_2.py [--results DIR] [--outdir DIR] [--domain START-END]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figlib import cohort as cohort_lib  # noqa: E402
from figlib import paths as paths_lib  # noqa: E402
from figlib import stats  # noqa: E402
from figlib import style  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

# set in main() from the command line
INPUTS: dict[str, Path] = {}
OUTPUT_DIR = Path(".")
ANALYSIS_REPORT_PATH = Path(".")
DOMAINS = pd.DataFrame()
REPEAT_BINS = pd.DataFrame()      # duplicons 07: 250-kb bins, raw maternal - paternal and envelope
SHIFTS = pd.DataFrame()           # per genome: level (CpG / repeat), group, inside, outside, shift
EXTRA_ROWS = pd.DataFrame()       # analysis 01: PWS-mUPD and biparental-mean combined beta per window
REPEAT_COLOR = "#8C510A"

MAIN_STEM = "Figure2"
SUPP_PROFILE_STEM = "Supplementary_Figure_S2A_retained_copy_profiles"
SUPP_DEPTH_STEM = "Supplementary_Figure_S2B_depth_robustness"
FORMATS = style.FORMATS
PNG_DPI = style.DPI
PAD_INCHES = 0.04

TITLE = "Reciprocal parental methylation architecture across the chromosome 15 imprinted domain"
PANEL_C_TITLE = "Phase-invariant allelic methylation on intact chromosome 15"
FIGURE_WIDTH_IN = 7.0
MAX_EXPORT_WIDTH_IN = 7.2
FIGURE_HEIGHT_IN = 7.7
BASE_FONT = style.BASE_FONT
INSET_RATIO = 2.0
A_YLIM_FACTOR = 1.6
FOCAL_LABEL_MERGE_MB = 0.02
JITTER = 0.07
JITTER_SEED = 11
NONEVALUABLE_MIN_BP = 5_000

MATERNAL = style.MATERNAL
PATERNAL = style.PATERNAL
TRACE = style.TRACE
INCONCLUSIVE = "#8A8A8A"
RAW = "#CFCFCF"
INK = style.INK
MUTED = style.MUTED
CONTROL = cohort_lib.COLOR["Control"]
DIGEORGE = cohort_lib.COLOR["DiGeorge"]
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
BASELINE_REGION = "OCA2 downstream control"
REPORT_NAME = "Figure2_report.md"
ASM_EXCESS_THRESHOLD = 0.02
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
    style.setup()


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
    for record in DOMAINS.itertuples():
        lo, hi = max(record.start, common[0]), min(record.end, common[1])
        if hi > lo:
            ax.axvspan(mb(lo), mb(hi), color=style.DOMAIN, alpha=0.8, lw=0, zorder=0)
            direction = str(getattr(record, "direction", "")).replace("_", "-")
            ax.text(mb(lo) + 0.02, limit * 0.93, f"{direction} repeat domain" if direction not in ("", "nan")
                    else "repeat methylation domain", ha="left", va="top", fontsize=BASE_FONT - 1.5, color=MUTED)
    if len(REPEAT_BINS):
        rb = REPEAT_BINS[(REPEAT_BINS["bin_end"] > common[0]) & (REPEAT_BINS["bin_start"] < common[1])]
        x = mb((rb["bin_start"] + rb["bin_end"]) / 2)
        if rb["env_low"].notna().any():
            ax.fill_between(x, rb["env_low"], rb["env_high"], step="mid", color=REPEAT_COLOR, alpha=0.12, lw=0, zorder=2)
        ax.step(x, rb["raw_delta"], where="mid", color=REPEAT_COLOR, lw=1.0, zorder=5)
        sig = rb["p"] <= 0.05
        ax.scatter(x[sig.to_numpy()], rb.loc[sig, "raw_delta"], s=14, marker="s", color=REPEAT_COLOR, zorder=6)
        ax.scatter(x[~sig.to_numpy()], rb.loc[~sig, "raw_delta"], s=14, marker="s", facecolor="white",
                   edgecolor=REPEAT_COLOR, lw=0.8, zorder=6)
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
    ax.set_ylabel("maternal − paternal\nmethylation")
    ax.tick_params(labelbottom=False)
    ax.text(1.005, 0.97, "maternal-\nretained\nhigher", transform=ax.transAxes, color=MATERNAL, fontsize=BASE_FONT - 1,
            va="top", ha="left")
    ax.text(1.005, 0.03, "paternal-\nretained\nhigher", transform=ax.transAxes, color=PATERNAL, fontsize=BASE_FONT - 1,
            va="bottom", ha="left")
    legend = [
        Line2D([], [], color=RAW, lw=1.2, label="CpG 1-kb \u0394\u03b2"),
        Line2D([], [], color=TRACE, lw=1.2, label="21-kb median (bootstrap CI)"),
        Patch(color="#BDBDBD", label="not evaluable"),
    ]
    if len(REPEAT_BINS):
        legend += [Line2D([], [], color=REPEAT_COLOR, marker="s", ms=3.5, lw=1.0,
                          label="repeat elements, 250 kb (filled p \u2264 0.05)"),
                   Patch(color=REPEAT_COLOR, alpha=0.12, label="label-permutation envelope")]
    ax.legend(handles=legend, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=3, fontsize=BASE_FONT - 1.5,
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
    ax.set_xlim(low, high)
    centres = np.array([item[0] for item in items])
    texts = [ax.text(c, 0.08, label, ha="center", va="bottom", fontsize=BASE_FONT - 1, color=color)
             for c, label, color in items]
    # measured label widths in data units; shrink the font when the labels cannot all fit
    renderer = ax.figure.canvas.get_renderer()
    data_per_px = (high - low) / max(1.0, ax.get_window_extent(renderer).width)
    for _ in range(4):
        widths = np.array([t.get_window_extent(renderer).width * data_per_px for t in texts]) * 1.12
        if widths.sum() <= (high - low) * 0.98 or texts[0].get_fontsize() <= BASE_FONT - 2.5:
            break
        for t in texts:
            t.set_fontsize(t.get_fontsize() - 0.5)
    positions = repel(centres, widths, low, high)
    for t, centre, position, width in zip(texts, centres, positions, widths):
        t.set_x(position)
        if abs(position - centre) > width * 0.05:
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
    o, v = np.abs(others.to_numpy(float)), np.abs(values.to_numpy(float))
    other_extent = np.max(o[np.isfinite(o)]) if np.isfinite(o).any() else np.nan
    extent = np.max(v[np.isfinite(v)]) if np.isfinite(v).any() else np.nan
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
    cohorts = (("Control", CONTROL, cohort_lib.MARKER["Control"], -0.17, "Control"),
               ("DiGeorge", DIGEORGE, cohort_lib.MARKER["DiGeorge"], 0.17, "22q11.2DS"))
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
        Line2D([], [], marker=cohort_lib.MARKER["Control"], ls="none", mfc="white", mec=CONTROL, ms=4,
               label=f"Control (n={totals['Control']})"),
        Line2D([], [], marker=cohort_lib.MARKER["DiGeorge"], ls="none", mfc="white", mec=DIGEORGE, ms=4,
               label=f"22q11.2DS (n={totals['DiGeorge']})"),
        Line2D([], [], color=MUTED, lw=2.0, marker="|", ms=0, label="group mean (equal weight)"),
        Line2D([], [], color=DIGEORGE, lw=1.4, label="95% participant bootstrap,\ndescriptive (only n≥3)"),
    ]
    main.legend(handles=legend, loc="upper left", bbox_to_anchor=(-0.02, -0.17), ncol=1, fontsize=BASE_FONT - 1,
                handlelength=1.3)
    return main


def draw_validation(fig, spec) -> list:
    """Panel c: per-genome shift of the domain at two levels."""
    levels = [lv for lv in ("CpG windows (analysis 01)", "repeat elements (duplicons 07)")
              if len(SHIFTS) and (SHIFTS["level"] == lv).any()]
    if not levels:
        ax = fig.add_subplot(spec)
        ax.text(0.5, 0.5, "domain validation not available\n(run scripts/duplicons 07-08)", ha="center",
                va="center", transform=ax.transAxes, color=MUTED)
        ax.set_axis_off()
        return [ax]
    sub = GridSpecFromSubplotSpec(1, len(levels), subplot_spec=spec, wspace=0.35)
    order = ["Biparental", "PWS-mUPD", "PWS-DEL", "AS-DEL"]
    colors = {"Biparental": style.BIPARENTAL, "PWS-mUPD": cohort_lib.COLOR["PWS-mUPD"],
              "PWS-DEL": cohort_lib.COLOR["PWS-DEL"], "AS-DEL": cohort_lib.COLOR["AS-DEL"]}
    markers = {"Biparental": "o", "PWS-mUPD": cohort_lib.MARKER["PWS-mUPD"], "PWS-DEL": cohort_lib.MARKER["PWS-DEL"],
               "AS-DEL": cohort_lib.MARKER["AS-DEL"]}
    rng = np.random.default_rng(5)
    axes = []
    for k, lv in enumerate(levels):
        ax = fig.add_subplot(sub[k], sharey=axes[0] if axes else None)
        d = SHIFTS[SHIFTS["level"] == lv]
        for i, g in enumerate(order):
            v = d.loc[d["group"] == g, "shift"].to_numpy()
            if not len(v):
                continue
            ax.scatter(i + rng.uniform(-0.12, 0.12, len(v)), v, s=16, marker=markers[g], color=colors[g],
                       edgecolor="white", lw=0.4, zorder=3)
            if len(v) > 1:
                ax.hlines(np.median(v), i - 0.28, i + 0.28, color=INK, lw=1.1, zorder=4)
        ax.axhline(0, color=MUTED, lw=0.6)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(["Bipar.", "mUPD", "PWS", "AS"], fontsize=BASE_FONT - 1)
        ax.set_xlim(-0.6, len(order) - 0.4)
        st = SHIFTS.attrs.get(lv, {})
        ax.set_title(f"{lv}\nPWS vs AS p = {style.fmt_p(st.get('p'))}; mUPD < "
                     f"{st.get('mupd_below', 0)}/{st.get('n_biparental', 0)} biparental", fontsize=BASE_FONT - 1.5)
        if k == 0:
            ax.set_ylabel("domain shift vs\nbiparental genomes")
        axes.append(ax)
    return axes


def render_main(tables: dict[str, pd.DataFrame], common: tuple[int, int], labels: dict) -> list[Path]:
    fig = plt.figure(figsize=(FIGURE_WIDTH_IN, 10.5))
    outer = GridSpec(8, 1, figure=fig, height_ratios=[1.3, 0.62, 1.85, 0.5, 1.0, 1.3, 0.85, 2.3], hspace=0.0,
                     left=0.2, right=0.9, top=0.94, bottom=0.075)
    ax_h = fig.add_subplot(outer[0])
    heatmap(ax_h, tables, common, extra=True, labels=labels)
    ax_h.tick_params(labelbottom=False)
    ax_a = fig.add_subplot(outer[2], sharex=ax_h)
    track = fig.add_subplot(outer[3], sharex=ax_h)
    draw_panel_a(ax_a, track, tables, common)
    mid = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[5], width_ratios=[1.0, 0.45], wspace=0.15)
    val_axes = draw_validation(fig, mid[0, 0])
    leg = fig.add_subplot(mid[0, 1])
    leg.set_axis_off()
    handles = [Line2D([], [], marker="o", ls="", ms=5, color=style.BIPARENTAL, label="biparental (Control + 22q11.2DS)")]
    handles += style.group_legend_handles(["PWS-mUPD", "PWS-DEL", "AS-DEL"])
    leg.legend(handles=handles, loc="center left", fontsize=BASE_FONT - 1, title="genome", title_fontsize=BASE_FONT - 1)
    lower = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[7], width_ratios=[1.0, 0.8], wspace=0.32)
    ax_b, order = draw_panel_b(fig, lower[0, 0], tables)
    ax_c = draw_panel_c(fig, lower[0, 1], tables, order)
    windows = tables["windows"]
    n_pws, n_as = int(windows["pws_n"].max()), int(windows["as_n"].max())
    left = 0.012
    for ax, letter, title, dy in (
            (ax_h, "a", f"Retained-copy methylation, common CN=1 interval (PWS n={n_pws}, AS n={n_as}); "
                        f"PWS-mUPD and biparental mean below", 0.008),
            (ax_a, "b", "Parental contrast of CpGs (1 kb) and repeat elements (250 kb)", 0.042),
            (val_axes[0], "c", "The domain in each genome: (inside − outside) relative to biparental genomes", 0.035),
            (ax_b, "d", "Prespecified regional parent-associated effects", 0.032)):
        top = ax.get_position().y1
        fig.text(left, top + dy, letter, fontsize=BASE_FONT + 3, fontweight="bold", va="bottom")
        fig.text(left + 0.028, top + dy, title, fontsize=BASE_FONT + 0.5, va="bottom")
    c_left = ax_c.get_position().x0 - 0.035
    b_top = ax_b.get_position().y1
    fig.text(c_left, b_top + 0.032, "e", fontsize=BASE_FONT + 3, fontweight="bold", va="bottom")
    fig.text(c_left + 0.028, b_top + 0.032, PANEL_C_TITLE.replace(" on intact", "\non intact"),
             fontsize=BASE_FONT + 0.5, va="bottom", linespacing=1.1)
    fig.suptitle(TITLE, x=left, y=0.99, ha="left", fontsize=BASE_FONT + 1.5, fontweight="bold")
    return save(fig, MAIN_STEM)


def heatmap(ax: plt.Axes, tables: dict[str, pd.DataFrame], common: tuple[int, int], extra: bool = False,
            labels: dict | None = None) -> None:
    data = tables["window_participants"]
    data = data[data["estimator"].eq("full_depth")]
    wide = data.pivot(index="sample_id", columns="start", values="beta")
    mechanism = data.drop_duplicates("sample_id").set_index("sample_id")["mechanism"]
    rows = sorted(wide.index, key=lambda sample: (mechanism[sample] != "PWS-DEL", sample))
    wide = wide.loc[rows]
    names = [f"{(labels or {}).get(r, r)} ({'mat.' if mechanism[r] == 'PWS-DEL' else 'pat.'})" for r in rows]
    if extra and len(EXTRA_ROWS):
        ex = EXTRA_ROWS.pivot(index="row", columns="start", values="beta").reindex(columns=wide.columns)
        wide = pd.concat([wide, ex])
        names += list(ex.index)
    cmap = plt.get_cmap("Purples").copy()
    cmap.set_bad("white")
    image = ax.imshow(np.ma.masked_invalid(wide.to_numpy(float)), aspect="auto", interpolation="nearest",
                      cmap=cmap, vmin=0, vmax=1,
                      extent=(mb(wide.columns.min()), mb(wide.columns.max() + 1000), len(wide) - 0.5, -0.5))
    ax.set_yticks(range(len(wide)))
    ax.set_yticklabels(names, fontsize=BASE_FONT - 1.5)
    n_pws = sum(mechanism[sample] == "PWS-DEL" for sample in rows)
    ax.axhline(n_pws - 0.5, color=INK, lw=0.8)
    if len(wide) > len(rows):
        ax.axhline(len(rows) - 0.5, color=INK, lw=1.4)
    for tick, name in zip(ax.get_yticklabels(), names):
        tick.set_color(MATERNAL if "(mat." in name or "mUPD" in name else PATERNAL if "(pat." in name else MUTED)
    bar = plt.colorbar(image, cax=ax.inset_axes([1.01, 0.0, 0.012, 1.0]))
    bar.set_label("β", fontsize=BASE_FONT - 1)
    bar.ax.tick_params(labelsize=BASE_FONT - 2)
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


def render_profiles(tables: dict[str, pd.DataFrame], common: tuple[int, int], labels: dict) -> list[Path]:
    windows = tables["windows"]
    windows = windows[windows["in_common_cn1"].astype(bool)].sort_values("start").reset_index(drop=True)
    fig, axes = plt.subplots(4, 1, figsize=(FIGURE_WIDTH_IN - 0.1, 7.6), sharex=True,
                             gridspec_kw={"height_ratios": [1.6, 1.1, 1.0, 0.8], "hspace": 0.28})
    heatmap(axes[0], tables, common, labels=labels)
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


def render_depth(tables: dict[str, pd.DataFrame], labels: dict) -> list[Path]:
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
    ax.set_xticklabels([f"{labels.get(sample, sample)}\n{'PWS' if depth.loc[sample, 'mechanism'] == 'PWS-DEL' else 'AS'}\n{depth.loc[sample, 'median_retained_copy_depth']:.0f}×"
                        for sample in samples], fontsize=BASE_FONT - 1)
    ax.set_ylabel("fraction of CN=1 windows\nnot evaluable (label: depth)")
    ax.set_ylim(0, min(1.0, max(0.05, 1.12 * float(np.nan_to_num(missing["missing_fraction"].max())))))
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

STATUS_TEXT = {
    "maternal_retained_higher": "maternal-retained higher",
    "paternal_retained_higher": "paternal-retained higher",
    "inconclusive": "inconclusive",
    "not_evaluable": "not evaluable",
}


def fmt(value: float, digits: int = 3, signed: bool = True) -> str:
    if value is None or not np.isfinite(value):
        return "NA"
    return f"{value:+.{digits}f}" if signed else f"{value:.{digits}f}"


def md_table(header: list[str], rows: list[list[object]]) -> str:
    def cell(value: object) -> str:
        return str(value).replace("|", "\\|")
    lines = ["| " + " | ".join(cell(item) for item in header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(cell(item) for item in row) + " |" for row in rows]
    return "\n".join(lines)


def overlapping_labels(start: int, end: int, tables: dict[str, pd.DataFrame]) -> str:
    labels = [
        record.label for record in tables["annotation"].itertuples()
        if record.start < end and record.end > start
    ]
    labels += [
        REGION_SHORT.get(record.region_id, record.region_id) for record in tables["regions"].itertuples()
        if np.isfinite(record.analysed_start) and record.analysed_start < end and record.analysed_end > start
    ]
    return ", ".join(dict.fromkeys(labels)) or "no annotated gene or region"


def asm_excess(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    values = tables["participants"]
    values = values[values["analysis"].eq("phase_invariant_asm") & values["estimator"].eq("full_depth")]
    wide = values.pivot_table(index=["sample_id", "cohort"], columns="region_id", values="value")
    if BASELINE_REGION not in wide:
        return pd.DataFrame()
    excess = wide.sub(wide[BASELINE_REGION], axis=0).drop(columns=BASELINE_REGION)
    long = excess.reset_index().melt(id_vars=["sample_id", "cohort"], var_name="region_id", value_name="excess")
    return (
        long.dropna(subset=["excess"])
        .groupby(["region_id", "cohort"])["excess"]
        .agg(["mean", "min", "max", "count"])
        .reset_index()
    )


def write_markdown_report(tables: dict[str, pd.DataFrame], common: tuple[int, int], figure_paths: list[Path]) -> Path:
    windows = tables["windows"]
    windows = windows[windows["in_common_cn1"].astype(bool)]
    regions = tables["regions"].sort_values("display_order")
    order = regions["region_id"].tolist()
    focal = tables["focal"]
    asm = tables["asm"][tables["asm"]["estimator"].eq("full_depth")]
    n_pws, n_as = int(windows["pws_n"].max()), int(windows["as_n"].max())
    cohort_n = (
        tables["participants"][tables["participants"]["analysis"].eq("phase_invariant_asm")]
        .groupby("cohort")["sample_id"].nunique().to_dict()
    )
    out: list[str] = []
    add = out.append
    add(f"# Figure 2 report\n\n**{TITLE}**\n")
    add("Generated automatically by `scripts/figures/FIGURE_2.py` from the analysis tables in "
        "`results/analysis/03_cis_architecture/`. Every number below is read from those tables.\n")
    png = next((path for path in figure_paths if path.name == f"{MAIN_STEM}.png"), None)
    if png is not None:
        add(f"![Figure 2](../figures/{png.name})\n")

    add("## 1. Design\n")
    add(md_table(
        ["Item", "Value"],
        [
            ["Common reciprocal CN=1 interval", f"chr15:{common[0]:,}–{common[1]:,} ({(common[1] - common[0]) / 1e6:.2f} Mb)"],
            ["PWS deletion (maternal copy retained)", f"n = {n_pws}"],
            ["AS deletion (paternal copy retained)", f"n = {n_as}"],
            ["Unaffected controls", f"n = {cohort_n.get('Control', 0)}"],
            ["DiGeorge (independent diploid disease controls)", f"n = {cohort_n.get('DiGeorge', 0)}"],
            ["1-kb windows in the interval", f"{len(windows):,}"],
            ["Evaluable windows", f"{int(windows['evaluable'].sum()):,} ({windows['evaluable'].mean():.1%})"],
        ],
    ))
    add("\nΔβ = equally weighted mean β of PWS (maternal-retained) − AS (paternal-retained). "
        "Positive: maternal-retained copy more methylated; negative: paternal-retained copy more methylated. "
        "Parental direction comes only from the deletions, and only inside the common CN=1 interval. "
        "Control and DiGeorge haplotypes (H1/H2) are unoriented, so panel c has no direction.\n")

    add("## 2. Panel b: CpG-level parental contrast and focal intervals\n")
    smooth = windows["delta_rolling_median"].dropna()
    add(f"- Smoothed trace available for {len(smooth):,} windows; median |Δβ| = {smooth.abs().median():.3f}; "
        f"{(smooth.abs() >= 0.10).mean():.1%} of smoothed windows reach |Δβ| ≥ 0.10.")
    reproducible = focal[focal["reproducible"].astype(bool)] if not focal.empty else focal
    add(f"- Candidate intervals: {len(focal)}; reproducible focal intervals: {len(reproducible)}.\n")
    if len(reproducible):
        add(md_table(
            ["ID", "Coordinates", "Length", "Direction", "Interval Δβ", "Leave-one-out min |Δβ|", "Overlaps"],
            [
                [
                    record.focal_id,
                    f"chr15:{int(record.start):,}–{int(record.end):,}",
                    f"{(record.end - record.start) / 1e3:.0f} kb",
                    STATUS_TEXT.get(record.direction, record.direction),
                    fmt(getattr(record, "interval_delta", record.mean_smoothed_delta)),
                    fmt(getattr(record, "loo_min_abs_delta", np.nan), signed=False),
                    overlapping_labels(int(record.start), int(record.end), tables),
                ]
                for record in reproducible.itertuples()
            ],
        ))
        add("")
    failed = focal[~focal["reproducible"].astype(bool)] if not focal.empty else focal
    if len(failed):
        add("Candidates that did not meet the prespecified criteria:\n")
        add(md_table(
            ["Candidate", "Coordinates", "Length", "Failed criteria"],
            [
                [record.candidate_id, f"chr15:{int(record.start):,}–{int(record.end):,}",
                 f"{(record.end - record.start) / 1e3:.0f} kb", record.failed_criteria]
                for record in failed.itertuples()
            ],
        ))
        add("")

    add("## 3. Panels b and c: repeat elements and the domain in each genome\n")
    if len(REPEAT_BINS):
        rb = REPEAT_BINS
        add(f"Repeat-element bins (scripts/duplicons 07, offset {REPEAT_BINS.attrs.get('offset', 0):+.3f} added back): "
            f"{int((rb['p'] <= 0.05).sum())} of {len(rb)} bins with label-permutation p <= 0.05.\n")
        add(md_table(["bin", "maternal − paternal", "envelope", "p"],
                     [[f"{int(r.bin_start):,}–{int(r.bin_end):,}", fmt(r.raw_delta), f"[{fmt(r.env_low)}, {fmt(r.env_high)}]",
                       fmt(r.p, signed=False)] for r in rb.itertuples()]))
        add("")
    if len(SHIFTS):
        add("Domain shift per genome = (inside − outside) − median of the biparental genomes; PWS-mUPD did not "
            "take part in finding the domain.\n")
        add("```\n" + SHIFTS.attrs.get("note", "") + "\n```\n")
        add(md_table(["level", "genome", "group", "inside", "outside", "shift"],
                     [[r.level.split(" (")[0], r.sample_id, r.group, fmt(r.inside, signed=False), fmt(r.outside, signed=False),
                       fmt(r.shift)] for r in SHIFTS.itertuples()]))
        add("")
    else:
        add("No scripts/duplicons 07-08 outputs: panels b (repeat bins) and c are not drawn.\n")

    add("## 4. Panel d: prespecified regional effects\n")
    rows = []
    for record in regions.itertuples():
        welch = f"[{fmt(getattr(record, 'welch_ci_low', np.nan))}, {fmt(getattr(record, 'welch_ci_high', np.nan))}]"
        rows.append([
            record.region_id, f"{record.pws_n}/{record.as_n}", fmt(record.delta_beta),
            f"[{fmt(record.ci_low)}, {fmt(record.ci_high)}]", welch,
            STATUS_TEXT.get(record.status, record.status), record.coverage_sensitivity,
            record.shared_cpg_sensitivity, "yes" if bool(record.robust) else "no",
        ])
    add(md_table(
        ["Region", "PWS/AS n", "Δβ", "Bootstrap 95% CI", "Welch 95% CI", "Call", "Depth check", "Shared-CpG check", "Robust"],
        rows,
    ))
    add("")
    for record in regions.itertuples():
        call = STATUS_TEXT.get(record.status, record.status)
        if record.status in {"maternal_retained_higher", "paternal_retained_higher"}:
            robust = "and passes both sensitivity checks" if bool(record.robust) else "but fails a sensitivity check"
            add(f"- **{record.region_id}**: {call} (Δβ {fmt(record.delta_beta)}), {robust}.")
        elif record.status == "inconclusive":
            reason = (
                "the CI excludes 0 but |Δβ| is below the effect-size floor"
                if np.isfinite(record.ci_low) and (record.ci_low > 0 or record.ci_high < 0)
                else "the CI includes 0"
            )
            add(f"- **{record.region_id}**: inconclusive (Δβ {fmt(record.delta_beta)}); {reason}.")
        else:
            add(f"- **{record.region_id}**: not evaluable in the common CN=1 interval.")
    add("")

    add("## 5. Panel e: phase-invariant allelic methylation on intact chromosome 15\n")
    excess = asm_excess(tables)
    rows = []
    for region in order:
        cells = [region]
        for cohort in ("Control", "DiGeorge"):
            record = asm[asm["region_id"].eq(region) & asm["cohort"].eq(cohort)]
            if record.empty:
                cells += ["NA", "NA"]
                continue
            record = record.iloc[0]
            interval = (
                f" [{fmt(record['ci_low'], signed=False)}, {fmt(record['ci_high'], signed=False)}]"
                if np.isfinite(record["ci_low"]) else ""
            )
            cells.append(f"{fmt(record['mean_absolute_asm'], signed=False)}{interval} (n={int(record['n_participants'])})")
            match = excess[excess["region_id"].eq(region) & excess["cohort"].eq(cohort)] if not excess.empty else excess
            cells.append(fmt(match["mean"].iloc[0]) if len(match) else ("baseline" if region == BASELINE_REGION else "NA"))
        rows.append(cells)
    add(md_table(
        ["Region", "Control mean |H1−H2|", "Control excess over OCA2", "DiGeorge mean |H1−H2|", "DiGeorge excess over OCA2"],
        rows,
    ))
    add(f"\n|H1−H2| is above zero even without allelic methylation, so each participant's value is compared with "
        f"their own {BASELINE_REGION} value (excess). Intervals are descriptive participant bootstraps, "
        f"shown only for n ≥ 3; controls are descriptive only.\n")

    add("## 6. Agreement between reciprocal deletions and intact chromosomes\n")
    rows = []
    for record in regions.itertuples():
        if record.region_id == BASELINE_REGION:
            continue
        match = excess[excess["region_id"].eq(record.region_id)] if not excess.empty else excess
        elevated = bool(len(match)) and bool((match["mean"] >= ASM_EXCESS_THRESHOLD).all())
        directional = record.status in {"maternal_retained_higher", "paternal_retained_higher"}
        if directional and elevated:
            verdict = "concordant: parental divergence and allelic methylation"
        elif not directional and not elevated:
            verdict = "concordant: no divergence detected"
        elif elevated:
            verdict = "allelic methylation without net parental Δβ (possible mixed-sign divergence)"
        else:
            verdict = "parental Δβ without elevated allelic methylation"
        rows.append([record.region_id, STATUS_TEXT.get(record.status, record.status),
                     "elevated" if elevated else "at baseline", verdict])
    add(md_table(["Region", "Panel b call", f"Panel c (excess ≥ {ASM_EXCESS_THRESHOLD} in both cohorts)", "Reading"], rows))
    add("")

    add("## 7. Sequencing-depth and CpG-composition robustness\n")
    by_estimator = tables["regions_by_estimator"]
    rows = []
    for region in order:
        cells = [region]
        for estimator in ESTIMATOR_STYLE:
            record = by_estimator[by_estimator["region_id"].eq(region) & by_estimator["estimator"].eq(estimator)]
            cells.append(fmt(record["delta_beta"].iloc[0]) if len(record) else "NA")
        rows.append(cells)
    add(md_table(["Region"] + [label for label, _, _ in ESTIMATOR_STYLE.values()], rows))
    add("")
    evaluable = tables["windows"][tables["windows"]["evaluable"].astype(bool)]
    for column, label in (("delta_beta_common_depth_downsampled", "downsampled"), ("delta_beta_shared_cpg", "shared-CpG")):
        pairs = evaluable[["delta_beta", column]].dropna()
        if len(pairs) > 2:
            add(f"- Window-level agreement, full depth vs {label}: r = {pairs.corr().iloc[0, 1]:.2f} ({len(pairs):,} windows).")
    missing = tables["missingness"]
    full = missing[missing["estimator"].eq("full_depth")]
    if len(full):
        add(f"- Full-depth missingness per participant: {full['missing_fraction'].min():.1%}–{full['missing_fraction'].max():.1%} "
            f"of CN=1 windows (retained-copy depth {full['median_retained_copy_depth'].min():.0f}–"
            f"{full['median_retained_copy_depth'].max():.0f}×).")
    add("- Depth handling tests robustness only; it does not add biological evidence.\n")

    if ANALYSIS_REPORT_PATH.is_file():
        checks = pd.read_csv(ANALYSIS_REPORT_PATH, sep="\t")
        checks = checks[checks["status"].isin(["pass", "fail"])]
        add("## 8. Analysis checks\n")
        add(md_table(["Section", "Check", "Result", "Status"],
                     [[r.section, r.item,
                       Path(str(r.value)).name if r.section == "input_validation" and "/" in str(r.value) else str(r.value)[:140],
                       r.status] for r in checks.itertuples()]))
        add("")

    add("## 9. Limitations\n")
    add(f"- Small groups (PWS n = {n_pws}, AS n = {n_as}, controls n = {cohort_n.get('Control', 0)}): percentile "
        "bootstrap intervals are descriptive and can be narrow; compare them with the Welch intervals.")
    add("- The PWS vs AS contrast also compares two diagnoses and two sets of individuals. Check the flanking "
        "diploid sequence (expected Δβ ≈ 0) and PWS-mUPD agreement before calling a pattern strictly parent-of-origin.")
    add("- No parental DNA: direction is assigned only from the deletions and only inside the common CN=1 interval.")
    add("- DiGeorge participants are an independent diploid disease-control cohort, not parental references.")
    add("- Focal intervals need a full 21-kb run of evaluable windows; shorter peaks (e.g. promoter DMRs) are "
        "visible in the 1-kb trace but cannot be called here.")
    path = OUTPUT_DIR / "reports" / REPORT_NAME
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(out) + "\n", encoding="utf-8")
    return path


def load_domains(paths, arg: str | None = None) -> pd.DataFrame:
    """The paternal/maternal repeat-methylation domain (scripts/duplicons 08), largest first."""
    if arg:
        a, b = arg.replace(",", "").split("-")
        return pd.DataFrame([{"start": int(a), "end": int(b)}])
    path = paths.duplicons / "followup" / "methylation_domains.tsv"
    if path.is_file() and path.stat().st_size:
        d = pd.read_csv(path, sep="\t")
        if {"start", "end"} <= set(d.columns) and len(d):
            return d.assign(size=d["end"] - d["start"]).sort_values("size", ascending=False)
    return pd.DataFrame(columns=["start", "end"])


def load_repeat_bins(paths) -> pd.DataFrame:
    d = paths.duplicons / "repeat_methylation"
    reg_path, val_path = d / "direct_delta_by_region.tsv", d / "validation.tsv"
    if not reg_path.is_file():
        return pd.DataFrame()
    reg = pd.read_csv(reg_path, sep="\t")
    offset = 0.0
    if val_path.is_file():
        val = pd.read_csv(val_path, sep="\t")
        hit = val[val["quantity"].astype(str).str.startswith("parental offset removed")]
        if len(hit):
            offset = float(hit["value"].iloc[0])
    out = pd.DataFrame({"bin_start": reg["bin_start"], "bin_end": reg["bin_end"],
                        "raw_delta": reg["observed_median_delta"] + offset})
    out["env_low"] = reg["permutation_min"] + offset if "permutation_min" in reg else np.nan
    out["env_high"] = reg["permutation_max"] + offset if "permutation_max" in reg else np.nan
    out["p"] = reg["permutation_p_two_sided"] if "permutation_p_two_sided" in reg else np.nan
    out.attrs["offset"] = offset
    return out.dropna(subset=["raw_delta"])


GROUP_OF = {"Control": "Biparental", "DiGeorge": "Biparental", "PWS-mUPD": "PWS-mUPD", "PWS-DEL": "PWS-DEL",
            "AS-DEL": "AS-DEL"}


def _shifts(per: pd.DataFrame, level: str) -> pd.DataFrame:
    """per: sample_id, group, inside, outside -> + shift relative to the biparental median."""
    per = per.dropna(subset=["inside", "outside"]).copy()
    per["in_minus_out"] = per["inside"] - per["outside"]
    ref = per.loc[per["group"] == "Biparental", "in_minus_out"].median()
    per["shift"] = per["in_minus_out"] - ref
    per["level"] = level
    return per


def load_matrix(paths) -> pd.DataFrame:
    path = paths.analysis / "01_evidence_matrix" / "chr15_window_evidence_matrix.tsv.gz"
    if not path.is_file():
        return pd.DataFrame()
    cols = ["sample_id", "mechanism", "track_kind", "start", "end", "n_cpg", "beta_site_mean"]
    m = pd.read_csv(path, sep="\t", usecols=cols, low_memory=False)
    return m[(m["track_kind"] == "combined") & (m["n_cpg"] >= 3)]


def cpg_domain_shifts(matrix: pd.DataFrame, core, domain) -> pd.DataFrame:
    if matrix.empty or domain is None:
        return pd.DataFrame()
    ic0, ic1 = 22_691_258 - 5_000, 22_693_494 + 5_000
    m = matrix[(matrix["start"] >= core[0]) & (matrix["end"] <= core[1])
               & ~((matrix["end"] > ic0) & (matrix["start"] < ic1))].copy()
    inside = (m["start"] >= domain[0]) & (m["end"] <= domain[1])
    per = (m.assign(region=np.where(inside, "inside", "outside"))
           .groupby(["sample_id", "mechanism", "region"])["beta_site_mean"].mean().unstack("region").reset_index())
    per["group"] = per["mechanism"].map(GROUP_OF)
    return _shifts(per.dropna(subset=["group"]), "CpG windows (analysis 01)")


def repeat_domain_shifts(paths, domain) -> pd.DataFrame:
    d = paths.duplicons / "repeat_methylation"
    ebp_path, summ_path = d / "element_by_participant.tsv.gz", d / "element_summary.tsv"
    if domain is None or not ebp_path.is_file() or not summ_path.is_file():
        return pd.DataFrame()
    summ = pd.read_csv(summ_path, sep="\t", usecols=["element_id", "start", "end", "testable", "in_common_core"])
    ok = summ["testable"].astype(str).isin(["True", "true", "1"]) & summ["in_common_core"].astype(str).isin(["True", "true", "1"])
    summ = summ[ok]
    region = np.where((summ["start"] >= domain[0]) & (summ["end"] <= domain[1]), "inside",
                      np.where((summ["end"] <= domain[0]) | (summ["start"] >= domain[1]), "outside", "edge"))
    reg = pd.Series(region, index=summ["element_id"])
    e = pd.read_csv(ebp_path, sep="\t", low_memory=False)
    keep = (e["mechanism"].isin(["PWS-DEL", "AS-DEL"]) & (e["measure"] == "retained_combined")) | \
           (~e["mechanism"].isin(["PWS-DEL", "AS-DEL"]) & (e["measure"] == "scaffold_combined"))
    e = e[keep & e["element_id"].isin(reg.index)].dropna(subset=["value"])
    e["region"] = e["element_id"].map(reg)
    e = e[e["region"] != "edge"]
    per = e.groupby(["sample_id", "mechanism", "region"])["value"].median().unstack("region").reset_index()
    per["group"] = per["mechanism"].map(GROUP_OF)
    return _shifts(per.dropna(subset=["group"]), "repeat elements (duplicons 07)")


def extra_rows(matrix: pd.DataFrame, cohort, core) -> pd.DataFrame:
    if matrix.empty:
        return pd.DataFrame()
    m = matrix[(matrix["start"] >= core[0]) & (matrix["end"] <= core[1])]
    rows = []
    for s in cohort.of("PWS-mUPD"):
        d = m[m["sample_id"] == s]
        rows.append(pd.DataFrame({"row": f"{cohort.label(s)} (2 × mat.)", "start": d["start"], "beta": d["beta_site_mean"]}))
    bi = m[m["mechanism"].isin(["Control", "DiGeorge"])]
    if len(bi):
        n = bi["sample_id"].nunique()
        g = bi.groupby("start")["beta_site_mean"].mean().reset_index()
        rows.append(pd.DataFrame({"row": f"biparental mean (n = {n})", "start": g["start"], "beta": g["beta_site_mean"]}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def shift_summary(shifts: pd.DataFrame, domain) -> str:
    lines = [f"Domain chr15:{int(domain[0]):,}–{int(domain[1]):,}" if domain is not None else "No domain"]
    for lv, d in shifts.groupby("level", sort=False):
        a = d.loc[d["group"] == "PWS-DEL", "shift"]
        b = d.loc[d["group"] == "AS-DEL", "shift"]
        u = d.loc[d["group"] == "PWS-mUPD", "shift"]
        bip = d.loc[d["group"] == "Biparental", "shift"]
        diff, p, n = stats.label_permutation(a, b)
        shifts.attrs[lv] = {"p": p, "diff": diff, "labellings": n}
        below = int((bip > u.median()).sum()) if len(u) else 0
        shifts.attrs[lv].update({"mupd_below": below, "n_biparental": len(bip)})
        lines.append(f"{lv.split(' (')[0]}: PWS-DEL {np.median(a):+.3f}, AS-DEL {np.median(b):+.3f}, "
                     f"PWS-mUPD {u.median():+.3f} (below {below} of {len(bip)} biparental); "
                     f"PWS vs AS p = {style.fmt_p(p)} ({n} labellings)")
    return "\n".join(lines)


def main(argv=None) -> None:
    global INPUTS, OUTPUT_DIR, ANALYSIS_REPORT_PATH, DOMAINS, REPEAT_BINS, SHIFTS, EXTRA_ROWS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    paths_lib.add_common_arguments(ap)
    ap.add_argument("--domain", help="START-END of the domain tested in panel c (default: scripts/duplicons 08)")
    a = ap.parse_args(argv)
    paths = paths_lib.resolve(a)
    cis = paths.analysis / "03_cis_architecture"
    INPUTS = {
        "windows": cis / "parent_associated_windows.tsv.gz",
        "window_participants": cis / "parent_window_participant_values.tsv.gz",
        "focal": cis / "focal_intervals.tsv",
        "regions": cis / "regional_parent_contrasts.tsv",
        "regions_by_estimator": cis / "regional_parent_contrasts_by_estimator.tsv",
        "participants": cis / "regional_participant_values.tsv.gz",
        "asm": cis / "regional_phase_invariant_asm.tsv",
        "missingness": cis / "participant_missingness.tsv",
        "core": paths.analysis / "01_evidence_matrix" / "common_reciprocal_cn1_core.tsv",
        "annotation": paths.analysis / "annotation_genes.tsv",
    }
    ANALYSIS_REPORT_PATH = cis / "figure2_analysis_report.tsv"
    OUTPUT_DIR = paths.figure_dir(2, a.outdir)
    setup_style()
    tables = load_inputs()
    common = (int(tables["core"]["start"].iloc[0]), int(tables["core"]["end"].iloc[0]))
    windows = tables["windows"]
    if windows.loc[~windows["in_common_cn1"].astype(bool), "delta_beta"].notna().any():
        raise ValueError("Window table contains parental contrasts outside the common CN=1 interval")
    cohort = cohort_lib.load(paths.metadata)
    labels = {s: cohort.label(s) for s in cohort.samples}
    DOMAINS = load_domains(paths, a.domain)
    domain = (int(DOMAINS["start"].iloc[0]), int(DOMAINS["end"].iloc[0])) if len(DOMAINS) else None
    REPEAT_BINS = load_repeat_bins(paths)
    matrix = load_matrix(paths)
    EXTRA_ROWS = extra_rows(matrix, cohort, common)
    parts = [cpg_domain_shifts(matrix, common, domain), repeat_domain_shifts(paths, domain)]
    parts = [x for x in parts if len(x)]
    SHIFTS = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    if len(SHIFTS):
        SHIFTS.attrs["note"] = shift_summary(SHIFTS, domain)
    style.write_table(windows[windows["in_common_cn1"].astype(bool)], OUTPUT_DIR, "Figure2b_cpg_window_contrast.tsv")
    style.write_table(REPEAT_BINS, OUTPUT_DIR, "Figure2b_repeat_bins.tsv")
    style.write_table(SHIFTS, OUTPUT_DIR, "Figure2c_domain_shifts.tsv")
    style.write_table(tables["regions"], OUTPUT_DIR, "Figure2d_regional_contrasts.tsv")
    style.write_table(tables["asm"], OUTPUT_DIR, "Figure2e_phase_invariant_asm.tsv")
    paths_out = render_main(tables, common, labels)
    paths_out += render_profiles(tables, common, labels)
    paths_out += render_depth(tables, labels)
    paths_out.append(write_render_log(paths_out, tables, common))
    paths_out.append(write_markdown_report(tables, common, paths_out))
    for path in paths_out:
        print(path)


if __name__ == "__main__":
    main()
