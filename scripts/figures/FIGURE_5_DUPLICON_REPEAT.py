#!/usr/bin/env python3
"""Render Figure 5 from existing duplicon and repeat-methylation summary tables.

Run the analysis/01, duplicons/07, and duplicons/08 stages first. No BAM, CpG
BED, or assembly is read here. The old FIGURE_5.py remains independently runnable.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
FOLLOWUP = ROOT / "results/08_duplicons/followup"
CORE_PATH = ROOT / "results/analysis/01_evidence_matrix/common_reciprocal_cn1_core.tsv"
OUTPUT = ROOT / "results/07_figures/figure_5"
PATHS = {
    "breakpoints": FOLLOWUP / "breakpoint_status.tsv",
    "regions": FOLLOWUP / "methylation_regions.tsv",
    "core": CORE_PATH,
}
REQUIRED = {
    "breakpoints": {"sample", "mechanism", "status", "nahr_consistent",
                    "junction_p_at_crossover", "assembled_fusion", "split_reads",
                    "split_reads_panel_max"},
    "regions": {"bin_start", "bin_end", "raw_median_delta", "permutation_p_two_sided",
                "pws_del_minus_biparental", "as_del_minus_biparental",
                "mupd_minus_biparental"},
    "core": {"start", "end"},
}

MATERNAL = "#C8472B"
PATERNAL = "#2F6DB0"
UPD = "#704B94"
INK = "#252525"
FAINT = "#E5E5E5"
NO_DATA = "#9A9A9A"
STATUS_SHORT = {
    "confirmed NAHR": "NAHR supported",
    "confirmed by assembly": "Assembly; mechanism open",
    "junction resolved at bp level": "Junction resolved",
    "compatible with NAHR": "NAHR compatible",
    "edges near one SD pair, not homologous": "Edges not homologous",
    "edge outside WINDOW": "Edge outside window",
    "unresolved": "Unresolved",
}


def load() -> dict[str, pd.DataFrame]:
    tables = {}
    for name, path in PATHS.items():
        if not path.is_file():
            raise FileNotFoundError(f"Missing {path}; run analysis/01, duplicons/07 and duplicons/08 first")
        table = pd.read_csv(path, sep="\t")
        missing = REQUIRED[name] - set(table.columns)
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        tables[name] = table
    if len(tables["core"]) != 1 or tables["breakpoints"].empty or tables["regions"].empty:
        raise ValueError("Expected one common CN=1 core, deletion carrier rows and repeat bins")
    if tables["breakpoints"]["sample"].duplicated().any():
        raise ValueError("Breakpoint table has duplicate participants")
    if tables["regions"]["bin_start"].duplicated().any():
        raise ValueError("Repeat table has duplicate bins")
    allowed = {"PWS-DEL", "AS-DEL"}
    if not set(tables["breakpoints"]["mechanism"]).issubset(allowed):
        raise ValueError("Breakpoint table must contain only PWS-DEL and AS-DEL carriers")
    return tables


def truth(value: object) -> bool | None:
    if pd.isna(value):
        return None
    if value is True or str(value).lower() in {"true", "1"}:
        return True
    if value is False or str(value).lower() in {"false", "0"}:
        return False
    return None


def evidence(row: pd.Series) -> list[bool | None]:
    p = pd.to_numeric(row["junction_p_at_crossover"], errors="coerce")
    homology = pd.to_numeric(row.get("fusion_homolog_distance_bp", np.nan), errors="coerce")
    split = pd.to_numeric(row["split_reads"], errors="coerce")
    panel = pd.to_numeric(row["split_reads_panel_max"], errors="coerce")
    fused = truth(row["assembled_fusion"])
    return [
        truth(row["nahr_consistent"]),
        bool(p < 1e-3) if pd.notna(p) else None,
        bool(homology <= 20_000) if pd.notna(homology) else (False if fused is False else None),
        bool(split >= 2 and panel == 0) if pd.notna(split) and pd.notna(panel) else None,
    ]


def plot_breakpoints(ax: plt.Axes, rows: pd.DataFrame) -> None:
    rows = rows.assign(group=rows.mechanism.map({"PWS-DEL": 0, "AS-DEL": 1}))
    rows = rows.sort_values(["group", "sample"]).reset_index(drop=True)
    ax.set_xlim(-0.6, 6.2)
    ax.set_ylim(len(rows) - 0.25, -1.5)
    headings = ["SD pair\ngeometry", "Junction\nexcess", "Homologous\ncontig", "Split-read\njunction"]
    ax.set_xticks(range(4), headings)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(rows)), rows["sample"])
    ax.set_ylabel("Deletion carrier")
    ax.tick_params(length=0, pad=4)
    for i, row in rows.iterrows():
        color = MATERNAL if row.mechanism == "PWS-DEL" else PATERNAL
        ax.axhspan(i - 0.5, i + 0.5, color="#FAFAFA" if i % 2 else "white", zorder=0)
        for j, value in enumerate(evidence(row)):
            if value is True:
                ax.scatter(j, i, s=51, marker="o", c=color, edgecolors="none", zorder=3)
            elif value is False:
                ax.scatter(j, i, s=43, marker="o", facecolors="white", edgecolors="#888888", lw=0.8, zorder=3)
            else:
                ax.text(j, i, "–", ha="center", va="center", color=NO_DATA, fontsize=9)
        status = STATUS_SHORT.get(str(row.status), str(row.status))
        ax.text(3.6, i, status, va="center", ha="left", fontsize=7, color=INK)
    ax.text(3.6, -0.95, "Combined verdict", ha="left", va="center", weight="bold", fontsize=7)
    n_pws = sum(rows.mechanism == "PWS-DEL")
    if 0 < n_pws < len(rows):
        ax.axhline(n_pws - 0.5, color=INK, lw=0.75)
    ax.set_xlabel("Filled = supporting evidence; open = assessed, absent; dash = unavailable", labelpad=9)
    ax.xaxis.set_label_position("bottom")
    for spine in ax.spines.values():
        spine.set_visible(False)


def numeric(data: pd.DataFrame, column: str) -> np.ndarray:
    return pd.to_numeric(data[column], errors="coerce").to_numpy(dtype=float)


def shared_coordinate(ax: plt.Axes, core: tuple[float, float]) -> None:
    ax.axvspan(*core, color="#F4F1EA", zorder=-2)
    for edge in core:
        ax.axvline(edge, lw=0.6, color="#9F9A91", ls=(0, (2, 2)), zorder=-1)
    ax.axhline(0, color=INK, lw=0.65, zorder=0)
    ax.grid(axis="y", color=FAINT, lw=0.55, zorder=-1)
    ax.spines[["top", "right"]].set_visible(False)


def plot_repeat_contrast(ax: plt.Axes, regions: pd.DataFrame, core: tuple[float, float]) -> None:
    x = (numeric(regions, "bin_start") + numeric(regions, "bin_end")) / 2e6
    width = (numeric(regions, "bin_end") - numeric(regions, "bin_start")) / 1e6
    delta = numeric(regions, "raw_median_delta")
    p = numeric(regions, "permutation_p_two_sided")
    shared_coordinate(ax, core)
    for xpos, w, val, pv in zip(x, width, delta, p):
        if not np.isfinite(val) or xpos - w / 2 < core[0] or xpos + w / 2 > core[1]:
            continue
        color = MATERNAL if val >= 0 else PATERNAL
        ax.bar(xpos, val, width=w * 0.85, color=color, alpha=0.75, edgecolor="none", zorder=2)
        # Nominal bin-level permutation support; avoid stars that imply corrected significance.
        if np.isfinite(pv) and pv <= 0.05 and abs(val) >= 0.02:
            ax.plot(xpos, val, marker="o", markersize=3.1, color=INK, linestyle="none", zorder=3)
    ax.set_ylabel("Repeat median Δβ\nPWS maternal − AS paternal")
    ax.set_title("b  Repeat methylation: direct contrast in the common CN=1 interval", loc="left", weight="bold")
    ax.text(0.99, 0.98, "Dot: nominal permutation p ≤ 0.05 and |Δβ| ≥ 0.02",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.5, color="#555555")


def plot_comparator(ax: plt.Axes, regions: pd.DataFrame, core: tuple[float, float]) -> None:
    x = (numeric(regions, "bin_start") + numeric(regions, "bin_end")) / 2e6
    shared_coordinate(ax, core)
    series = [
        ("pws_del_minus_biparental", MATERNAL, "PWS-DEL retained maternal", "o"),
        ("as_del_minus_biparental", PATERNAL, "AS-DEL retained paternal", "s"),
        ("mupd_minus_biparental", UPD, "PWS-mUPD (one participant)", "D"),
    ]
    for column, color, label, symbol in series:
        y = numeric(regions, column)
        ok = (np.isfinite(x) & np.isfinite(y) &
              (numeric(regions, "bin_start") / 1e6 >= core[0]) &
              (numeric(regions, "bin_end") / 1e6 <= core[1]))
        # Points, without joining missing bins, make each 250-kb summary explicit.
        ax.scatter(x[ok], y[ok], s=13 if symbol != "D" else 17, marker=symbol,
                   color=color, label=label, zorder=2)
    ax.set_xlabel("T2T-CHM13v2.0 chr15 position (Mb); shading = common reciprocal CN=1 core")
    ax.set_ylabel("Repeat β minus biparental\ncombined median")
    ax.set_title("c  Descriptive check against Control + DiGeorge combined tracks", loc="left", weight="bold", pad=27)
    ax.legend(frameon=False, fontsize=6.5, ncol=3, loc="lower left",
              bbox_to_anchor=(0, 1.005), borderaxespad=0, handletextpad=0.3, columnspacing=0.8)


def run() -> list[Path]:
    tables = load()
    rows = tables["breakpoints"]
    regions = tables["regions"].sort_values("bin_start").reset_index(drop=True)
    core_row = tables["core"].iloc[0]
    core = (float(core_row.start) / 1e6, float(core_row.end) / 1e6)
    if core[0] >= core[1]:
        raise ValueError("Invalid common CN=1 interval")
    positions = numeric(regions, "bin_start")
    if not np.all(np.isfinite(positions)):
        raise ValueError("Non-numeric repeat bin coordinates")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7,
                         "axes.titlesize": 8, "axes.labelsize": 7,
                         "xtick.labelsize": 7, "ytick.labelsize": 7,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(7.05, 8.1), layout="constrained")
    grid = fig.add_gridspec(3, 1, height_ratios=[2.05, 1.0, 1.15], hspace=0.13)
    ax_a, ax_b, ax_c = (fig.add_subplot(grid[i, 0]) for i in range(3))
    plot_breakpoints(ax_a, rows)
    ax_a.set_title("a  Breakpoint mechanism: evidence by deletion carrier", loc="left", weight="bold", pad=9)
    plot_repeat_contrast(ax_b, regions, core)
    plot_comparator(ax_c, regions, core)
    limits = (core[0] - 0.15, core[1] + 0.15)
    ax_b.set_xlim(*limits)
    ax_c.set_xlim(*limits)
    ax_b.tick_params(labelbottom=False)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    paths = [OUTPUT / f"Figure5_duplicon_repeat.{ext}" for ext in ("png", "pdf", "svg")]
    for path in paths:
        fig.savefig(path, dpi=600 if path.suffix == ".png" else None, facecolor="white")
    plt.close(fig)
    counts = "\n".join(f"- {status}: {n}" for status, n in rows["status"].value_counts().items())
    (OUTPUT / "Figure5_duplicon_repeat_report.md").write_text(
        "# Figure 5 data audit\n\n"
        f"Deletion carriers with status rows: {len(rows)}. "
        f"Repeat bins: {len(regions)}. Common CN=1: chr15:{int(core_row.start):,}–{int(core_row.end):,}.\n\n"
        f"## Breakpoint verdicts\n\n{counts}\n\n"
        "Panel a shows the current 08_followup verdict; support types are distinct. "
        "Panel b shows the raw PWS-DEL minus AS-DEL repeat-element bin median; "
        "dots are nominal permutation support and are not multiplicity corrected. "
        "Panel c is a descriptive comparison of group medians, including one mUPD participant. "
        "Panels b and c must not be read as independently replicated element-level tests.\n"
    )
    return paths


if __name__ == "__main__":
    for output in run():
        print(output)
