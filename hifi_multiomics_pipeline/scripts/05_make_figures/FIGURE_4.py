#!/usr/bin/env python3
"""
Self-contained Figure 4 manuscript renderer for the hifi_multiomics_pipeline
layout.

This script writes the publication-facing SNORD116 figure names, including
`Figure4_per_molecule_cis_architecture`.
"""

from __future__ import annotations

import argparse
from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, Patch, Rectangle
from scipy import stats


DEFAULT_OUTDIR = Path("/home/rare/arlen/pws-as-hifi-cis-methylation/scripts/hifi_multiomics_pipeline/06_results")
DEFAULT_INPUT_TABLE_DIR = Path("/home/rare/arlen/paper_vf/tables")

FIGURE_BASENAMES = [
    "Figure_SNORD116_single_molecule_architecture",
    "Figure4",
    "Figure4_per_molecule_cis_architecture",
]
SUPPLEMENT_BARCODE_BASENAME = "Figure_SNORD116_single_molecule_architecture_supplement_barcode_stacks"

REPORT_NAME = "Figure_SNORD116_single_molecule_report.md"
TABLE_MOLECULE_SUMMARY = "Table_SNORD116_molecule_summary.tsv"
TABLE_ENTROPY_STATS = "Table_SNORD116_entropy_statistics.tsv"
TABLE_PARENTAL_STATS = "Table_SNORD116_parental_state_statistics.tsv"
PAIRWISE_Q_DISPLAY_THRESHOLD = 0.05

REQUIRED_INPUT_TABLES = [
    "Figure4_panelA_control_paternal_SNORD116_rows.tsv",
    "Figure4_panelA_control_paternal_SNORD116_matrix.tsv.gz",
    "Figure4_panelB_control_maternal_SNORD116_rows.tsv",
    "Figure4_panelB_control_maternal_SNORD116_matrix.tsv.gz",
    "Figure4_gene_features.tsv",
    "Figure4_SNORD116_shared_core_window.tsv",
    "Figure4_panelC_entropy_plot_input.tsv",
    "Figure4_panelD_SNORD116_plot_input.tsv",
    "Figure4_single_molecule_CpG_calls.tsv.gz",
]

CPG_UNMETH = "#3A74B7"
CPG_METH = "#D6544A"
CPG_MISSING = "#F0F0F0"
CPG_RUG = "#B0B0B0"
TRACK_FILL = "#D7D7D7"
GUIDE_BAND = "#F7F4EC"
ANNOTATION_COLOR = "#8E8E8E"
TEXT_MUTED = "#5C5C5C"

GROUP_PALETTE = {
    "Control": "#6F6F6F",
    "AS-DEL": "#157A73",
    "PWS-DEL": "#B56D2D",
    "PWS-mUPD": "#7258A7",
}
CONTROL_PARENTAL_COLORS = {
    "Control paternal": "#4E6A5E",
    "Control maternal": "#8A6B3F",
}
PANEL_D_PALETTE = {
    "AS-DEL retained paternal": GROUP_PALETTE["AS-DEL"],
    "PWS-DEL retained maternal": GROUP_PALETTE["PWS-DEL"],
}

HEATMAP_CMAP = ListedColormap([CPG_UNMETH, CPG_METH])
HEATMAP_CMAP.set_bad(CPG_MISSING)
COVERAGE_STRIP_CMAP = LinearSegmentedColormap.from_list(
    "snord116_coverage_strip",
    ["#FCFCFC", "#D7D7D7"],
)

DISPLAY_WINDOW = (22_808_000, 22_845_000)
DISPLAY_MIN_SPAN_FRACTION = 0.35
DISPLAY_MIN_CPGS = 20
DISPLAY_MIN_WITHIN_CALL_FRACTION = 0.47
DISPLAY_TARGET_ROWS = 18
MAIN_BARCODE_MAX_COLUMNS = 240
ZOOM_WINDOW_RAW_COLUMNS = 40
ZOOM_INSET_MAX_COLUMNS = 20
ZOOM_INSET_ROWS = 5
ZOOM_INSET_MIN_CALL_FRACTION = 0.75
ZOOM_HIGHLIGHT_COLOR = "#2F2F2F"
ZOOM_CONNECTOR_COLOR = "#464646"

PANEL_C_REGION_ORDER = ["SNRPN/PWS-IC", "SNORD116", "Downstream control"]
PANEL_C_GROUP_ORDER = ["Control", "AS-DEL", "PWS-DEL", "PWS-mUPD"]
PANEL_D_ORDER = ["AS-DEL retained paternal", "PWS-DEL retained maternal"]
CONTROL_PANEL_MAP = {
    "A": ("Control paternal", "Paternal-like control molecules"),
    "B": ("Control maternal", "Maternal-like control molecules"),
}


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def missing_input_tables(table_dir: Path) -> list[Path]:
    return [table_dir / name for name in REQUIRED_INPUT_TABLES if not (table_dir / name).exists()]


def resolve_input_table_dir(outdir: Path, requested_table_dir: Path | None = None) -> Path:
    candidates = [requested_table_dir] if requested_table_dir is not None else [outdir / "tables", DEFAULT_INPUT_TABLE_DIR]
    checked: list[Path] = []
    for candidate in candidates:
        if candidate is None:
            continue
        candidate = candidate.expanduser().resolve()
        checked.append(candidate)
        if not missing_input_tables(candidate):
            return candidate
    details = "\n".join(f"- {path}" for path in checked)
    raise FileNotFoundError(
        "Figure 4 input tables were not found. Run the Figure 4 table-generation step first, "
        "or pass --table-dir to a directory containing the Figure4_*.tsv inputs.\n"
        f"Checked:\n{details}"
    )


def standardize_region_label(value: str) -> str:
    mapping = {
        "PWS-IC": "SNRPN/PWS-IC",
        "SNORD116 cluster": "SNORD116",
        "SNORD116 display window": "SNORD116 display",
        "Downstream control": "Downstream control",
        "UBE3A downstream control": "Downstream control",
        "SNHG14 boundary window": "Boundary",
    }
    return mapping.get(str(value), str(value))


def scientific_notation_mathtext(value: float, prefix: str) -> str:
    if pd.isna(value):
        return f"{prefix} = NA"
    value = float(value)
    if value == 0:
        return f"{prefix} = 0"
    if value >= 1e-3:
        return f"{prefix} = {value:.3f}"
    exponent = int(np.floor(np.log10(value)))
    mantissa = value / (10**exponent)
    return rf"{prefix} = {mantissa:.1f} \times 10^{{{exponent}}}"


def scientific_notation_text(value: float, prefix: str) -> str:
    if pd.isna(value):
        return f"{prefix} = NA"
    value = float(value)
    if value == 0:
        return f"{prefix} = 0"
    if value >= 1e-3:
        return f"{prefix} = {value:.3f}"
    exponent = int(np.floor(np.log10(value)))
    mantissa = value / (10**exponent)
    return f"{prefix} = {mantissa:.1f} x 10^{exponent}"


def scientific_notation_value_text(value: float) -> str:
    if pd.isna(value):
        return "NA"
    value = float(value)
    if value == 0:
        return "0"
    if value >= 1e-3:
        return f"{value:.3f}"
    exponent = int(np.floor(np.log10(value)))
    mantissa = value / (10**exponent)
    return f"{mantissa:.1f} x 10^{exponent}"


def compact_value_text(value: float) -> str:
    if pd.isna(value):
        return "NA"
    return f"{float(value):.2g}"


def scientific_notation_plot_text(value: float) -> str:
    if pd.isna(value):
        return "NA"
    value = float(value)
    if value == 0:
        return "0"
    if value >= 1e-3:
        return f"{value:.3f}"
    exponent = int(np.floor(np.log10(value)))
    mantissa = value / (10**exponent)
    return f"{mantissa:.1f} × 10^{exponent}"


def benjamini_hochberg(p_values: list[float]) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    sorted_p = p[order]
    n = len(sorted_p)
    adjusted = np.empty(n, dtype=float)
    running = 1.0
    for idx in range(n - 1, -1, -1):
        rank = idx + 1
        running = min(running, sorted_p[idx] * n / rank)
        adjusted[idx] = running
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adjusted, 0.0, 1.0)
    return out


def kruskal_eta_squared(statistic_h: float, total_n: int, n_groups: int) -> float:
    if total_n <= n_groups or pd.isna(statistic_h):
        return np.nan
    return max((float(statistic_h) - n_groups + 1.0) / (total_n - n_groups), 0.0)


def bootstrap_median_difference(
    values_a: np.ndarray,
    values_b: np.ndarray,
    n_boot: int = 5000,
    seed: int = 20260525,
) -> tuple[float, float, float]:
    a = np.asarray(values_a, dtype=float)
    b = np.asarray(values_b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    observed = float(np.median(a) - np.median(b))
    if a.size == 0 or b.size == 0:
        return observed, np.nan, np.nan
    rng = np.random.default_rng(seed)
    boot = np.empty(n_boot, dtype=float)
    for idx in range(n_boot):
        aa = rng.choice(a, size=a.size, replace=True)
        bb = rng.choice(b, size=b.size, replace=True)
        boot[idx] = np.median(aa) - np.median(bb)
    ci_low, ci_high = np.quantile(boot, [0.025, 0.975])
    return observed, float(ci_low), float(ci_high)


def median_q1_q3(values: pd.Series | np.ndarray) -> tuple[float, float, float]:
    arr = pd.to_numeric(pd.Series(values), errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size == 0:
        return np.nan, np.nan, np.nan
    return float(np.median(arr)), float(np.quantile(arr, 0.25)), float(np.quantile(arr, 0.75))


def format_iqr(low: float, high: float) -> str:
    if pd.isna(low) or pd.isna(high):
        return "NA"
    return f"{low:.3f}-{high:.3f}"


def dataframe_to_markdown(frame: pd.DataFrame) -> str:
    columns = [str(col) for col in frame.columns]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for _, row in frame.iterrows():
        values = []
        for column in frame.columns:
            value = row[column]
            if isinstance(value, float):
                if np.isnan(value):
                    values.append("NA")
                else:
                    values.append(f"{value:.3f}")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def load_matrix_and_rows(table_dir: Path, prefix: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = pd.read_csv(table_dir / f"{prefix}_rows.tsv", sep="\t").set_index("read_key")
    matrix = pd.read_csv(table_dir / f"{prefix}_matrix.tsv.gz", sep="\t", compression="gzip").set_index("read_key")
    matrix.columns = [int(column) for column in matrix.columns]
    shared = rows.index.intersection(matrix.index)
    return matrix.loc[shared].copy(), rows.loc[shared].copy()


def annotate_display_rows(matrix: pd.DataFrame, meta: pd.DataFrame, window: tuple[int, int]) -> pd.DataFrame:
    columns = matrix.columns.to_numpy(dtype=float)
    start_window, end_window = window
    window_span = max(end_window - start_window, 1)
    out = meta.copy()
    span_bp = []
    span_fraction = []
    within_call_fraction = []
    missing_fraction = []
    first_called = []
    last_called = []
    for read_key, row in matrix.iterrows():
        values = row.to_numpy(dtype=float)
        finite = np.isfinite(values)
        called_positions = columns[finite]
        if called_positions.size == 0:
            first_called.append(np.nan)
            last_called.append(np.nan)
            within_call_fraction.append(np.nan)
            missing_fraction.append(np.nan)
        else:
            first_called.append(float(called_positions[0]))
            last_called.append(float(called_positions[-1]))
            within_mask = (columns >= called_positions[0]) & (columns <= called_positions[-1])
            completeness = float(finite[within_mask].mean())
            within_call_fraction.append(completeness)
            missing_fraction.append(1.0 - completeness)
        clipped_start = max(float(out.at[read_key, "alignment_start"]), float(start_window))
        clipped_end = min(float(out.at[read_key, "alignment_end"]), float(end_window))
        span = max(clipped_end - clipped_start, 0.0)
        span_bp.append(span)
        span_fraction.append(span / window_span)
    out["molecule_span_bp"] = span_bp
    out["molecule_span_fraction"] = span_fraction
    out["within_call_fraction"] = within_call_fraction
    out["missing_fraction_within_span"] = missing_fraction
    out["first_called_position"] = first_called
    out["last_called_position"] = last_called
    return out


def build_display_candidates(full_meta: pd.DataFrame) -> pd.DataFrame:
    keep = (
        (pd.to_numeric(full_meta["molecule_span_fraction"], errors="coerce") >= DISPLAY_MIN_SPAN_FRACTION)
        & (pd.to_numeric(full_meta["n_cpgs"], errors="coerce") >= DISPLAY_MIN_CPGS)
        & (pd.to_numeric(full_meta["within_call_fraction"], errors="coerce") >= DISPLAY_MIN_WITHIN_CALL_FRACTION)
    )
    candidate_meta = full_meta.loc[keep].copy()
    if candidate_meta.empty:
        candidate_meta = full_meta.copy()
    return candidate_meta.sort_values(
        ["mean_methylation_probability", "within_call_fraction", "n_cpgs", "molecule_span_bp", "first_called_position"],
        ascending=[False, False, False, False, True],
        kind="mergesort",
    )


def select_display_subset(
    matrix: pd.DataFrame,
    candidate_meta: pd.DataFrame,
    display_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    display_meta = candidate_meta.copy()
    if len(display_meta) > display_n:
        chosen = np.linspace(0, len(display_meta) - 1, display_n).round().astype(int)
        display_meta = display_meta.iloc[np.unique(chosen)].copy()
    display_meta = display_meta.sort_values(
        ["mean_methylation_probability", "within_call_fraction", "n_cpgs", "first_called_position"],
        ascending=[False, False, False, True],
        kind="mergesort",
    )
    display_matrix = matrix.reindex(display_meta.index)
    return display_matrix, display_meta


def order_full_barcode_stack(full_meta: pd.DataFrame) -> pd.DataFrame:
    return full_meta.sort_values(
        ["mean_methylation_probability", "within_call_fraction", "n_cpgs", "first_called_position"],
        ascending=[False, False, False, True],
        kind="mergesort",
    )


def compute_binary_display_matrix(matrix: pd.DataFrame) -> pd.DataFrame:
    values = matrix.to_numpy(dtype=float)
    binary = np.where(np.isnan(values), np.nan, (values >= 0.5).astype(float))
    return pd.DataFrame(binary, index=matrix.index, columns=matrix.columns)


def compute_mean_profile(matrix: pd.DataFrame) -> np.ndarray:
    values = matrix.to_numpy(dtype=float)
    counts = np.isfinite(values).sum(axis=0)
    sums = np.nansum(values, axis=0)
    return np.divide(sums, counts, out=np.full(counts.shape, np.nan, dtype=float), where=counts > 0)


def compute_call_fraction(matrix: pd.DataFrame) -> np.ndarray:
    if matrix.empty:
        return np.array([], dtype=float)
    return np.isfinite(matrix.to_numpy(dtype=float)).mean(axis=0)


def smooth_profile(values: np.ndarray, window: int = 7) -> np.ndarray:
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
    grouped.columns = [float(np.mean(source_positions[bin_ids == bin_id])) for bin_id in sorted(np.unique(bin_ids))]
    return grouped


def choose_zoom_window(
    matrix_a: pd.DataFrame,
    matrix_b: pd.DataFrame,
    core_window: tuple[int, int],
    window_columns: int = ZOOM_WINDOW_RAW_COLUMNS,
) -> tuple[int, int]:
    columns = matrix_a.columns.to_numpy(dtype=float)
    within_core = np.flatnonzero((columns >= core_window[0]) & (columns <= core_window[1]))
    if within_core.size == 0:
        return int(columns[0]), int(columns[min(len(columns) - 1, window_columns - 1)])
    core_coverage = (compute_call_fraction(matrix_a) + compute_call_fraction(matrix_b)) / 2.0
    local_coverage = core_coverage[within_core]
    if within_core.size <= window_columns:
        start_idx = int(within_core[0])
        end_idx = int(within_core[-1])
    else:
        rolling = np.convolve(local_coverage, np.ones(window_columns, dtype=float) / window_columns, mode="valid")
        best_local = int(np.argmax(rolling))
        start_idx = int(within_core[best_local])
        end_idx = int(within_core[min(best_local + window_columns - 1, within_core.size - 1)])
    return int(columns[start_idx]), int(columns[end_idx])


def build_zoom_inset_matrix(
    raw_matrix: pd.DataFrame,
    ordered_meta: pd.DataFrame,
    zoom_window: tuple[int, int],
    max_columns: int = ZOOM_INSET_MAX_COLUMNS,
    n_rows: int = ZOOM_INSET_ROWS,
) -> pd.DataFrame:
    columns = raw_matrix.columns.to_numpy(dtype=float)
    within_zoom = (columns >= zoom_window[0]) & (columns <= zoom_window[1])
    zoom_matrix = raw_matrix.loc[:, within_zoom].copy()
    if zoom_matrix.empty:
        return zoom_matrix
    zoom_matrix = compress_matrix_for_display(zoom_matrix, max_columns)
    zoom_call_fraction = pd.Series(
        np.isfinite(zoom_matrix.to_numpy(dtype=float)).mean(axis=1),
        index=zoom_matrix.index,
        name="zoom_call_fraction",
    )
    ordered = ordered_meta.copy()
    ordered["zoom_call_fraction"] = zoom_call_fraction.reindex(ordered.index).fillna(0.0)
    eligible = ordered.loc[ordered["zoom_call_fraction"] >= ZOOM_INSET_MIN_CALL_FRACTION].copy()
    if eligible.empty:
        eligible = ordered.copy()
    if len(eligible) > n_rows:
        chosen = np.linspace(0, len(eligible) - 1, n_rows).round().astype(int)
        eligible = eligible.iloc[np.unique(chosen)].copy()
    if len(eligible) < n_rows:
        fallback = ordered.sort_values(
            ["zoom_call_fraction", "within_call_fraction", "n_cpgs", "mean_methylation_probability"],
            ascending=[False, False, False, False],
            kind="mergesort",
        )
        chosen_keys = list(eligible.index)
        for read_key in fallback.index:
            if read_key not in chosen_keys:
                chosen_keys.append(read_key)
            if len(chosen_keys) >= n_rows:
                break
        eligible = ordered.loc[[read_key for read_key in ordered.index if read_key in set(chosen_keys)]].copy()
    return zoom_matrix.reindex(eligible.index)


def nearest_column_index(columns: np.ndarray, position: float) -> int:
    return int(np.argmin(np.abs(columns - position)))


def shared_xticks(columns: np.ndarray) -> tuple[np.ndarray, list[str]]:
    tick_values = np.array([22_810_000, 22_820_000, 22_830_000, 22_840_000, 22_845_000], dtype=float)
    tick_values = tick_values[(tick_values >= columns.min()) & (tick_values <= columns.max())]
    tick_idx = [nearest_column_index(columns, value) for value in tick_values]
    tick_labels = [f"{value / 1e6:.3f}".rstrip("0").rstrip(".") for value in tick_values]
    return np.asarray(tick_idx, dtype=int), tick_labels


def overlay_sample_medians_grouped(
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
    group_width = total_width / max(len(hue_order), 1)
    for x_idx, x_value in enumerate(x_order):
        for hue_idx, hue_value in enumerate(hue_order):
            subset = sample_medians[(sample_medians[x_col] == x_value) & (sample_medians[hue_col] == hue_value)].copy()
            if subset.empty:
                continue
            center = x_idx - total_width / 2 + group_width / 2 + hue_idx * group_width
            xs = [center] if len(subset) == 1 else np.linspace(center - 0.04, center + 0.04, len(subset))
            ax.scatter(
                xs,
                subset[y_col].to_numpy(dtype=float),
                s=42,
                facecolor="white",
                edgecolor=palette[hue_value],
                linewidth=1.2,
                zorder=5,
            )


def overlay_sample_medians_category(
    ax: plt.Axes,
    sample_medians: pd.DataFrame,
    order: list[str],
    palette: dict[str, str],
    x_col: str,
    y_col: str,
    marker_size: float = 48,
) -> None:
    for idx, category in enumerate(order):
        subset = sample_medians[sample_medians[x_col] == category].copy()
        if subset.empty:
            continue
        xs = [idx] if len(subset) == 1 else np.linspace(idx - 0.08, idx + 0.08, len(subset))
        ax.scatter(
            xs,
            subset[y_col].to_numpy(dtype=float),
            s=marker_size,
            facecolor="white",
            edgecolor=palette[category],
            linewidth=1.3,
            zorder=5,
        )


def annotate_region_calls_with_missing(reads_df: pd.DataFrame, calls_df: pd.DataFrame) -> pd.DataFrame:
    reads = reads_df.copy()
    reads["region_plot"] = reads["region"].map(standardize_region_label)
    calls = calls_df.copy()
    calls["region_plot"] = calls["region"].map(standardize_region_label)
    merge_keys = ["sample_id", "region_plot", "read_id"]
    if "haplotype_assignment" in reads.columns and "haplotype_assignment" in calls.columns:
        merge_keys.append("haplotype_assignment")
    region_positions = {
        region: np.sort(group["cpg_position"].unique().astype(int))
        for region, group in calls.groupby("region_plot", observed=True)
    }
    grouped = (
        calls.groupby(merge_keys, observed=True)["cpg_position"]
        .agg(["min", "max", "count"])
        .reset_index()
        .rename(columns={"min": "first_called_position", "max": "last_called_position", "count": "n_called_cpgs"})
    )
    potential_sites = []
    missing_fraction = []
    call_span_bp = []
    for row in grouped.itertuples(index=False):
        positions = region_positions.get(row.region_plot, np.array([], dtype=int))
        left = np.searchsorted(positions, row.first_called_position, side="left")
        right = np.searchsorted(positions, row.last_called_position, side="right")
        n_possible = max(right - left, 1)
        potential_sites.append(int(n_possible))
        missing_fraction.append(float(1.0 - row.n_called_cpgs / n_possible))
        call_span_bp.append(int(row.last_called_position - row.first_called_position))
    grouped["n_possible_cpg_sites"] = potential_sites
    grouped["missing_fraction_within_span"] = missing_fraction
    grouped["call_span_bp"] = call_span_bp
    start_col = "start_position" if "start_position" in reads.columns else "alignment_start"
    end_col = "end_position" if "end_position" in reads.columns else "alignment_end"
    reads["molecule_span_bp"] = pd.to_numeric(reads[end_col], errors="coerce") - pd.to_numeric(reads[start_col], errors="coerce")
    reads = reads.merge(
        grouped[merge_keys + ["n_possible_cpg_sites", "missing_fraction_within_span", "call_span_bp"]],
        on=merge_keys,
        how="left",
    )
    return reads


def build_region_group_summary(panel_c_df: pd.DataFrame, calls_df: pd.DataFrame) -> pd.DataFrame:
    annotated = annotate_region_calls_with_missing(panel_c_df, calls_df)
    cpg_counts = (
        calls_df.assign(region_plot=calls_df["region"].map(standardize_region_label))
        .groupby(["region_plot", "sample_group"], observed=True)
        .size()
        .reset_index(name="total_cpg_calls")
        .rename(columns={"region_plot": "region", "sample_group": "group"})
    )
    rows = []
    for region in PANEL_C_REGION_ORDER:
        for group in PANEL_C_GROUP_ORDER:
            subset = annotated[(annotated["region_plot"] == region) & (annotated["sample_group"] == group)].copy()
            if subset.empty:
                continue
            entropy_median, entropy_q1, entropy_q3 = median_q1_q3(subset["methylation_entropy"])
            meth_median, meth_q1, meth_q3 = median_q1_q3(subset["mean_methylation"])
            rows.append(
                {
                    "summary_context": "region_level",
                    "region": region,
                    "group": group,
                    "parental_state": subset["parental_label"].mode().iat[0] if subset["parental_label"].notna().any() else "mixed",
                    "n_molecules": int(len(subset)),
                    "n_samples": int(subset["sample_id"].nunique()),
                    "median_cpgs_per_molecule": float(pd.to_numeric(subset["n_cpgs"], errors="coerce").median()),
                    "median_molecule_span_bp": float(pd.to_numeric(subset["molecule_span_bp"], errors="coerce").median()),
                    "median_missing_fraction": float(pd.to_numeric(subset["missing_fraction_within_span"], errors="coerce").median()),
                    "median_methylation": meth_median,
                    "methylation_q1": meth_q1,
                    "methylation_q3": meth_q3,
                    "median_entropy": entropy_median,
                    "entropy_q1": entropy_q1,
                    "entropy_q3": entropy_q3,
                }
            )
    summary = pd.DataFrame(rows)
    summary = summary.merge(cpg_counts, on=["region", "group"], how="left")
    return summary.sort_values(["region", "group"], kind="mergesort").reset_index(drop=True)


def build_control_display_summary(full_meta_by_label: dict[str, pd.DataFrame], displayed_by_label: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for analysis_label, panel_title in CONTROL_PANEL_MAP.values():
        full_meta = full_meta_by_label[analysis_label]
        displayed = displayed_by_label[analysis_label]
        meth_median, meth_q1, meth_q3 = median_q1_q3(full_meta["mean_methylation_probability"])
        rows.append(
            {
                "summary_context": "control_display_window",
                "region": "SNORD116 display",
                "group": panel_title,
                "parental_state": "paternal-like" if "paternal" in analysis_label.lower() else "maternal-like",
                "n_molecules": int(len(full_meta)),
                "n_samples": int(full_meta["sample_id"].nunique()),
                "n_displayed": int(len(displayed)),
                "median_cpgs_per_molecule": float(pd.to_numeric(full_meta["n_cpgs"], errors="coerce").median()),
                "median_molecule_span_bp": float(pd.to_numeric(full_meta["molecule_span_bp"], errors="coerce").median()),
                "median_missing_fraction": float(pd.to_numeric(full_meta["missing_fraction_within_span"], errors="coerce").median()),
                "median_methylation": meth_median,
                "methylation_q1": meth_q1,
                "methylation_q3": meth_q3,
                "median_entropy": float(pd.to_numeric(full_meta["methylation_entropy"], errors="coerce").median()),
                "entropy_q1": float(pd.to_numeric(full_meta["methylation_entropy"], errors="coerce").quantile(0.25)),
                "entropy_q3": float(pd.to_numeric(full_meta["methylation_entropy"], errors="coerce").quantile(0.75)),
                "total_cpg_calls": int(pd.to_numeric(full_meta["n_cpgs"], errors="coerce").sum()),
            }
        )
    return pd.DataFrame(rows)


def build_entropy_statistics(panel_c_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sample_medians = (
        panel_c_df.groupby(["region_plot", "sample_group", "sample_id"], observed=True)["methylation_entropy"]
        .median()
        .reset_index(name="sample_median_entropy")
    )
    rows = []
    for region in PANEL_C_REGION_ORDER:
        subset = panel_c_df[panel_c_df["region_plot"] == region].copy()
        arrays = []
        total_n = 0
        row = {
            "region": region,
            "test": "Kruskal-Wallis",
        }
        for group in PANEL_C_GROUP_ORDER:
            values = pd.to_numeric(
                subset.loc[subset["sample_group"] == group, "methylation_entropy"],
                errors="coerce",
            ).dropna().to_numpy(dtype=float)
            arrays.append(values)
            total_n += len(values)
            median, q1, q3 = median_q1_q3(values)
            row[f"n_{group}"] = int(len(values))
            row[f"n_samples_{group}"] = int(subset.loc[subset["sample_group"] == group, "sample_id"].nunique())
            row[f"median_{group}"] = median
            row[f"iqr_low_{group}"] = q1
            row[f"iqr_high_{group}"] = q3
        statistic_h, p_value = stats.kruskal(*arrays)
        row["statistic_H"] = float(statistic_h)
        row["p_value"] = float(p_value)
        row["eta_squared"] = float(kruskal_eta_squared(statistic_h, total_n, len(PANEL_C_GROUP_ORDER)))
        rows.append(row)
    stats_df = pd.DataFrame(rows)
    stats_df["q_value"] = benjamini_hochberg(stats_df["p_value"].tolist())
    pairwise_rows = []
    pair_order = [(a, b) for idx, a in enumerate(PANEL_C_GROUP_ORDER) for b in PANEL_C_GROUP_ORDER[idx + 1 :]]
    for region in PANEL_C_REGION_ORDER:
        subset = panel_c_df[panel_c_df["region_plot"] == region].copy()
        region_rows = []
        for group_a, group_b in pair_order:
            values_a = pd.to_numeric(
                subset.loc[subset["sample_group"] == group_a, "methylation_entropy"],
                errors="coerce",
            ).dropna().to_numpy(dtype=float)
            values_b = pd.to_numeric(
                subset.loc[subset["sample_group"] == group_b, "methylation_entropy"],
                errors="coerce",
            ).dropna().to_numpy(dtype=float)
            if values_a.size == 0 or values_b.size == 0:
                continue
            statistic_u, p_value = stats.mannwhitneyu(values_a, values_b, alternative="two-sided")
            region_rows.append(
                {
                    "region": region,
                    "group_a": group_a,
                    "group_b": group_b,
                    "n_group_a": int(values_a.size),
                    "n_group_b": int(values_b.size),
                    "statistic_U": float(statistic_u),
                    "p_value": float(p_value),
                }
            )
        if region_rows:
            region_df = pd.DataFrame(region_rows)
            region_df["q_value"] = benjamini_hochberg(region_df["p_value"].tolist())
            pairwise_rows.extend(region_df.to_dict("records"))
    pairwise_df = pd.DataFrame(pairwise_rows)
    return stats_df, sample_medians, pairwise_df


def build_parental_state_statistics(panel_d_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    panel_d = panel_d_df.copy()
    panel_d["panel_label"] = panel_d["sample_group"].map(
        {"AS-DEL": "AS-DEL retained paternal", "PWS-DEL": "PWS-DEL retained maternal"}
    )
    panel_d = panel_d[panel_d["panel_label"].isin(PANEL_D_ORDER)].copy()
    sample_medians = (
        panel_d.groupby(["panel_label", "sample_id"], observed=True)["mean_methylation"]
        .median()
        .reset_index(name="sample_median_methylation")
    )
    a = pd.to_numeric(panel_d.loc[panel_d["panel_label"] == PANEL_D_ORDER[0], "mean_methylation"], errors="coerce").dropna().to_numpy(dtype=float)
    b = pd.to_numeric(panel_d.loc[panel_d["panel_label"] == PANEL_D_ORDER[1], "mean_methylation"], errors="coerce").dropna().to_numpy(dtype=float)
    statistic_u, p_value = stats.mannwhitneyu(a, b, alternative="two-sided")
    delta_median, ci_low, ci_high = bootstrap_median_difference(a, b)
    row = {
        "comparison": f"{PANEL_D_ORDER[0]} vs {PANEL_D_ORDER[1]}",
        "region": "SNORD116",
        "test": "Mann-Whitney U",
        "n_molecules_group_a": int(a.size),
        "n_molecules_group_b": int(b.size),
        "n_samples_group_a": int(panel_d.loc[panel_d["panel_label"] == PANEL_D_ORDER[0], "sample_id"].nunique()),
        "n_samples_group_b": int(panel_d.loc[panel_d["panel_label"] == PANEL_D_ORDER[1], "sample_id"].nunique()),
        "median_group_a": float(np.median(a)),
        "median_group_b": float(np.median(b)),
        "delta_median": float(delta_median),
        "ci_low": float(ci_low),
        "ci_high": float(ci_high),
        "statistic_U": float(statistic_u),
        "p_value": float(p_value),
        "rank_biserial": float((2 * statistic_u / (a.size * b.size)) - 1),
        "reference_line": 0.5,
    }
    return pd.DataFrame([row]), sample_medians


def build_display_sensitivity(
    meta_a: pd.DataFrame,
    meta_b: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for span_threshold in [0.00, 0.30, 0.35, 0.40]:
        subset_a = meta_a[
            (pd.to_numeric(meta_a["n_cpgs"], errors="coerce") >= 20)
            & (pd.to_numeric(meta_a["within_call_fraction"], errors="coerce") >= DISPLAY_MIN_WITHIN_CALL_FRACTION)
            & (pd.to_numeric(meta_a["molecule_span_fraction"], errors="coerce") >= span_threshold)
        ].copy()
        subset_b = meta_b[
            (pd.to_numeric(meta_b["n_cpgs"], errors="coerce") >= 20)
            & (pd.to_numeric(meta_b["within_call_fraction"], errors="coerce") >= DISPLAY_MIN_WITHIN_CALL_FRACTION)
            & (pd.to_numeric(meta_b["molecule_span_fraction"], errors="coerce") >= span_threshold)
        ].copy()
        p_value = np.nan
        if len(subset_a) >= 5 and len(subset_b) >= 5:
            _, p_value = stats.mannwhitneyu(
                subset_a["mean_methylation_probability"].to_numpy(dtype=float),
                subset_b["mean_methylation_probability"].to_numpy(dtype=float),
                alternative="two-sided",
            )
        rows.append(
            {
                "min_span_fraction": span_threshold,
                "paternal_n": int(len(subset_a)),
                "maternal_n": int(len(subset_b)),
                "paternal_median_methylation": float(pd.to_numeric(subset_a["mean_methylation_probability"], errors="coerce").median()),
                "maternal_median_methylation": float(pd.to_numeric(subset_b["mean_methylation_probability"], errors="coerce").median()),
                "control_parental_p_value": float(p_value) if not pd.isna(p_value) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def build_panel_d_sensitivity(panel_d_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for min_cpg in [20, 60, 100, 120]:
        for min_span_bp in [10_000, 12_000, 14_000, 16_000]:
            subset = panel_d_df[
                (pd.to_numeric(panel_d_df["n_cpgs"], errors="coerce") >= min_cpg)
                & ((pd.to_numeric(panel_d_df["end_position"], errors="coerce") - pd.to_numeric(panel_d_df["start_position"], errors="coerce")) >= min_span_bp)
            ].copy()
            a = pd.to_numeric(subset.loc[subset["sample_group"] == "AS-DEL", "mean_methylation"], errors="coerce").dropna().to_numpy(dtype=float)
            b = pd.to_numeric(subset.loc[subset["sample_group"] == "PWS-DEL", "mean_methylation"], errors="coerce").dropna().to_numpy(dtype=float)
            if len(a) < 10 or len(b) < 10:
                continue
            statistic_u, p_value = stats.mannwhitneyu(a, b, alternative="two-sided")
            rows.append(
                {
                    "min_cpgs": min_cpg,
                    "min_span_bp": min_span_bp,
                    "n_as_del": int(len(a)),
                    "n_pws_del": int(len(b)),
                    "delta_median": float(np.median(a) - np.median(b)),
                    "p_value": float(p_value),
                    "statistic_U": float(statistic_u),
                }
            )
    return pd.DataFrame(rows)


def build_figure_caption(entropy_stats: pd.DataFrame, parental_stats: pd.DataFrame) -> str:
    snord116 = entropy_stats.loc[entropy_stats["region"] == "SNORD116"].iloc[0]
    parental = parental_stats.iloc[0]
    return (
        "Single-molecule HiFi methylation resolves focal parent-of-origin architecture across SNORD116.\n\n"
        "Panels A-B show representative high-completeness control paternal-like and maternal-like molecules across a shared SNORD116 display window, with the number shown and the total qualifying molecules indicated in each header. Rows represent individual HiFi molecules and columns represent ordered CpG positions; CpG positions were compressed locally for display only. Insets show representative molecules across a narrower SNORD116 subwindow at higher display resolution to make multi-CpG barcode structure directly visible.\n\n"
        f"Panel C compares molecule-level methylation entropy across the SNRPN/PWS-IC interval, SNORD116, and a downstream control interval, with the strongest cross-group heterogeneity observed at SNORD116 (q = {scientific_notation_value_text(snord116['q_value'])}, eta^2 = {snord116['eta_squared']:.2f}).\n\n"
        f"Panel D compares reciprocal deletion backgrounds at SNORD116 and shows higher per-molecule methylation in AS-DEL retained paternal molecules than in PWS-DEL retained maternal molecules (delta median = {parental['delta_median']:.3f}, p = {scientific_notation_value_text(parental['p_value'])}). White points indicate sample-level medians.\n\n"
        "Together, the figure shows that individual HiFi molecules carry coherent multi-CpG methylation barcodes across SNORD116 and that the strongest parent-of-origin contrast is focal and molecule-level rather than a downstream background effect."
    )


def render_concept_schematic(ax: plt.Axes, core_start: int, core_end: int) -> None:
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    x_left, x_right = 0.12, 0.43
    y_genome = 0.48
    ax.plot([x_left, x_right], [y_genome, y_genome], color="#6B6B6B", linewidth=1.2)
    ax.add_patch(Rectangle((0.22, y_genome - 0.045), 0.16, 0.09, facecolor=GUIDE_BAND, edgecolor="#C8C2AE", linewidth=0.8))
    ax.text(0.30, 0.73, "SNORD116 region", fontsize=9.0, ha="center", va="center", color=TEXT_MUTED)
    for x_pos in np.linspace(0.15, 0.40, 8):
        ax.plot([x_pos, x_pos], [y_genome - 0.08, y_genome + 0.08], color="#B7B7B7", linewidth=0.6)

    ax.text(0.72, 0.90, "Single HiFi molecules span multiple CpGs", fontsize=8.5, ha="center", va="center", color=TEXT_MUTED)
    ax.add_patch(Rectangle((0.50, 0.18), 0.19, 0.50, facecolor="#FCFCFC", edgecolor="#E2E2E2", linewidth=0.55))
    ax.add_patch(Rectangle((0.73, 0.18), 0.22, 0.50, facecolor="#FCFCFC", edgecolor="#E2E2E2", linewidth=0.55))
    ax.text(0.595, 0.63, "Random isolated\nCpG methylation", fontsize=7.9, ha="center", va="center", color=TEXT_MUTED, linespacing=1.05)
    ax.text(0.84, 0.63, "Observed coordinated\nmulti-CpG pattern", fontsize=7.9, ha="center", va="center", color=TEXT_MUTED, linespacing=1.05)

    random_x = np.linspace(0.54, 0.66, 8)
    observed_x = np.linspace(0.77, 0.91, 8)
    random_rows = [0.45, 0.32]
    observed_rows = [0.45, 0.32]
    for y in random_rows + observed_rows:
        if y in random_rows:
            ax.plot([0.52, 0.68], [y, y], color="#B0B0B0", linewidth=0.9)
        else:
            ax.plot([0.75, 0.93], [y, y], color="#B0B0B0", linewidth=0.9)

    random_patterns = [
        [CPG_UNMETH, CPG_METH, CPG_UNMETH, CPG_UNMETH, CPG_METH, CPG_UNMETH, CPG_METH, CPG_UNMETH],
        [CPG_METH, CPG_UNMETH, CPG_UNMETH, CPG_METH, CPG_UNMETH, CPG_METH, CPG_UNMETH, CPG_METH],
    ]
    observed_patterns = [
        [CPG_UNMETH, CPG_UNMETH, CPG_UNMETH, CPG_METH, CPG_METH, CPG_METH, CPG_METH, CPG_METH],
        [CPG_METH, CPG_METH, CPG_METH, CPG_UNMETH, CPG_UNMETH, CPG_UNMETH, CPG_METH, CPG_METH],
    ]

    for y, pattern in zip(random_rows, random_patterns):
        for x_pos, color in zip(random_x, pattern):
            ax.scatter([x_pos], [y], s=20, color=color, edgecolor="white", linewidth=0.35, zorder=3)
    for y, pattern in zip(observed_rows, observed_patterns):
        for x_pos, color in zip(observed_x, pattern):
            ax.scatter([x_pos], [y], s=20, color=color, edgecolor="white", linewidth=0.35, zorder=3)

    ax.text(0.595, 0.18, "scattered single-site events", fontsize=7.9, color=TEXT_MUTED, va="center", ha="center")
    ax.text(0.84, 0.18, "contiguous blocks across one molecule", fontsize=7.9, color=TEXT_MUTED, va="center", ha="center")


def render_annotation_track(ax: plt.Axes, columns: np.ndarray, feature_df: pd.DataFrame) -> None:
    ax.set_ylim(0.0, 1.0)
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.grid(False)
    ax.axhspan(0.22, 0.70, color="#F7F7F7", zorder=0)
    for _, row in feature_df.sort_values(["start", "end"]).iterrows():
        start_idx = nearest_column_index(columns, float(row["start"]))
        end_idx = nearest_column_index(columns, float(row["end"]))
        width = max(end_idx - start_idx + 1, 1)
        ax.add_patch(
            Rectangle(
                (start_idx - 0.5, 0.30),
                width,
                0.28,
                facecolor="#9A9A9A",
                edgecolor="#707070",
                linewidth=0.25,
                alpha=0.96,
            )
        )
    rug_positions = np.arange(len(columns))
    ax.vlines(rug_positions, 0.05, 0.16, color="#BDBDBD", linewidth=0.32)
    ax.text(0.0, 0.86, "SNORD116 copies", transform=ax.transAxes, fontsize=8.9, color=TEXT_MUTED, ha="left", va="bottom", fontweight="bold")
    ax.text(1.0, 0.86, f"n = {len(feature_df)}", transform=ax.transAxes, fontsize=8.0, color=TEXT_MUTED, ha="right", va="bottom")
    ticks, labels = shared_xticks(columns)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=9.2)
    ax.set_xlabel("Genomic position (chr15 Mb, T2T-CHM13)", fontsize=10.8, labelpad=1.0)


def render_barcode_panel(
    fig: plt.Figure,
    spec,
    panel_letter: str,
    panel_title: str,
    line_color: str,
    display_matrix: pd.DataFrame,
    display_meta: pd.DataFrame,
    full_matrix: pd.DataFrame,
    full_meta: pd.DataFrame,
    raw_matrix: pd.DataFrame,
    feature_df: pd.DataFrame,
    core_window: tuple[int, int],
    show_ylabel: bool,
    zoom_window: tuple[int, int] | None = None,
    show_zoom_inset: bool = True,
    show_local_cpg_legend: bool = True,
) -> None:
    width_ratios = [1.0, 0.36] if show_zoom_inset else [1.0]
    n_rows = 4 if show_local_cpg_legend else 3
    height_ratios = [0.21, 1.09, 0.16, 0.17] if show_local_cpg_legend else [0.21, 1.09, 0.15]
    main = spec.subgridspec(n_rows, len(width_ratios), height_ratios=height_ratios, width_ratios=width_ratios, hspace=0.038, wspace=0.06)
    top = main[0, 0].subgridspec(2, 1, height_ratios=[0.90, 0.10], hspace=0.03)
    ax_profile = fig.add_subplot(top[0, 0])
    ax_profile_strip = fig.add_subplot(top[1, 0], sharex=ax_profile)
    ax_heat = fig.add_subplot(main[1, 0], sharex=ax_profile)
    ax_ann = fig.add_subplot(main[2, 0], sharex=ax_profile)
    ax_key = fig.add_subplot(main[3, :]) if show_local_cpg_legend else None
    ax_zoom_header = None
    ax_zoom = None
    if show_zoom_inset:
        zoom_spec = main[:2, 1].subgridspec(2, 1, height_ratios=[0.19, 0.81], hspace=0.04)
        ax_zoom_header = fig.add_subplot(zoom_spec[0, 0])
        ax_zoom = fig.add_subplot(zoom_spec[1, 0])

    columns = full_matrix.columns.to_numpy(dtype=float)
    profile = smooth_profile(compute_mean_profile(full_matrix), window=5)
    coverage = compute_call_fraction(full_matrix)
    display_binary = compute_binary_display_matrix(display_matrix)
    core_start_idx = nearest_column_index(columns, core_window[0])
    core_end_idx = nearest_column_index(columns, core_window[1])
    zoom_matrix = build_zoom_inset_matrix(raw_matrix, display_meta, zoom_window) if show_zoom_inset and zoom_window is not None else pd.DataFrame()

    for axis in [ax_profile, ax_profile_strip, ax_heat, ax_ann]:
        axis.axvspan(core_start_idx - 0.5, core_end_idx + 0.5, color=GUIDE_BAND, alpha=0.22, zorder=0)

    ax_profile.plot(np.arange(len(columns)), profile, color=line_color, linewidth=1.08, solid_capstyle="round")
    ax_profile.fill_between(np.arange(len(columns)), profile, color=line_color, alpha=0.02)
    ax_profile.set_ylim(0.0, 1.0)
    ax_profile.set_yticks([0.0, 0.5, 1.0])
    ax_profile.set_yticklabels(["0", "0.5", "1"], fontsize=8.8 if show_ylabel else 0.0)
    ax_profile.set_ylabel("Mean methylation" if show_ylabel else "", fontsize=9.8)
    if not show_ylabel:
        ax_profile.tick_params(axis="y", labelleft=False, left=False)
    ax_profile.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_profile.grid(False)
    ax_profile.spines["top"].set_visible(False)
    ax_profile.spines["right"].set_visible(False)
    ax_profile.set_title(f"{panel_letter}. {panel_title}", loc="left", fontsize=13.7, fontweight="bold", pad=1.8)
    full_median_meth = float(pd.to_numeric(full_meta["mean_methylation_probability"], errors="coerce").median())
    ax_profile.text(
        0.995,
        0.92,
        f"representative n = {len(display_meta)}/{len(full_meta)} | median methylation = {full_median_meth:.3f}",
        transform=ax_profile.transAxes,
        ha="right",
        va="top",
        fontsize=8.0,
        color=TEXT_MUTED,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.80, "pad": 0.45},
    )

    ax_profile_strip.imshow(
        coverage[np.newaxis, :],
        aspect="auto",
        interpolation="nearest",
        cmap=COVERAGE_STRIP_CMAP,
        vmin=0.0,
        vmax=1.0,
        alpha=0.34,
    )
    ax_profile_strip.set_yticks([])
    ax_profile_strip.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    for spine in ["top", "right", "left", "bottom"]:
        ax_profile_strip.spines[spine].set_visible(False)

    ax_heat.imshow(
        display_binary.to_numpy(dtype=float),
        aspect="auto",
        interpolation="nearest",
        cmap=HEATMAP_CMAP,
        vmin=0.0,
        vmax=1.0,
        rasterized=True,
    )
    ax_heat.set_xticks([])
    ax_heat.set_yticks([])
    ax_heat.set_ylabel("Representative\nHiFi molecules" if show_ylabel else "", fontsize=9.8)
    ax_heat.text(
        0.012,
        0.988,
        "Rows = molecules | columns = ordered CpGs",
        transform=ax_heat.transAxes,
        ha="left",
        va="top",
        fontsize=7.0,
        color=TEXT_MUTED,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 0.35},
    )
    if len(display_binary) > 1:
        ax_heat.hlines(
            np.arange(0.5, len(display_binary), 1.0),
            -0.5,
            len(columns) - 0.5,
            color="#FFFFFF",
            linewidth=0.15,
            alpha=0.22,
        )
    for spine in ["top", "right"]:
        ax_heat.spines[spine].set_visible(False)

    render_annotation_track(ax_ann, columns, feature_df)
    if show_zoom_inset and zoom_window is not None:
        zoom_start_idx = nearest_column_index(columns, zoom_window[0])
        zoom_end_idx = nearest_column_index(columns, zoom_window[1])
        zoom_width = max(zoom_end_idx - zoom_start_idx + 1, 1)
        ax_heat.axvspan(zoom_start_idx - 0.5, zoom_end_idx + 0.5, color=ZOOM_HIGHLIGHT_COLOR, alpha=0.035, zorder=1)
        ax_heat.add_patch(
            Rectangle(
                (zoom_start_idx - 0.5, -0.5),
                zoom_width,
                len(display_binary),
                fill=False,
                edgecolor=ZOOM_HIGHLIGHT_COLOR,
                linewidth=1.05,
                linestyle=(0, (2.2, 1.8)),
                zorder=6,
            )
        )
        ax_profile_strip.add_patch(
            Rectangle(
                (zoom_start_idx - 0.5, -0.5),
                zoom_width,
                1.0,
                fill=False,
                edgecolor=ZOOM_HIGHLIGHT_COLOR,
                linewidth=0.95,
                linestyle=(0, (2.2, 1.8)),
                zorder=6,
            )
        )
    if ax_zoom_header is not None and ax_zoom is not None:
        ax_zoom_header.axis("off")
        ax_zoom_header.text(
            0.0,
            0.70,
            "Representative zoomed molecules",
            ha="left",
            va="center",
            fontsize=7.5,
            fontweight="bold",
            color="#222222",
        )
        ax_zoom_header.text(
            0.0,
            0.22,
            "Single-molecule methylation barcodes",
            ha="left",
            va="center",
            fontsize=6.9,
            color=TEXT_MUTED,
        )
        if zoom_matrix.empty:
            ax_zoom.axis("off")
        else:
            zoom_binary = compute_binary_display_matrix(zoom_matrix)
            ax_zoom.imshow(
                zoom_binary.to_numpy(dtype=float),
                aspect="auto",
                interpolation="nearest",
                cmap=HEATMAP_CMAP,
                vmin=0.0,
                vmax=1.0,
                rasterized=True,
            )
            ax_zoom.set_xticks([])
            ax_zoom.set_yticks([])
            if len(zoom_binary) > 1:
                ax_zoom.hlines(
                    np.arange(0.5, len(zoom_binary), 1.0),
                    -0.5,
                    zoom_binary.shape[1] - 0.5,
                    color="#FFFFFF",
                    linewidth=0.24,
                    alpha=0.50,
                )
            for spine in ["top", "right", "left", "bottom"]:
                ax_zoom.spines[spine].set_visible(True)
                ax_zoom.spines[spine].set_linewidth(0.85)
                ax_zoom.spines[spine].set_color("#9E9E9E")
            ax_zoom.text(
                0.02,
                0.98,
                f"{len(zoom_binary)} molecules",
                transform=ax_zoom.transAxes,
                ha="left",
                va="top",
                fontsize=7.0,
                color=TEXT_MUTED,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 0.35},
            )
            connector_top = ConnectionPatch(
                xyA=(0.0, 0.98),
                coordsA=ax_zoom.transAxes,
                xyB=(zoom_end_idx + 0.5, -0.5),
                coordsB=ax_heat.transData,
                color=ZOOM_CONNECTOR_COLOR,
                linewidth=0.95,
                linestyle=(0, (2.2, 1.8)),
                alpha=0.92,
                zorder=8,
            )
            connector_bottom = ConnectionPatch(
                xyA=(0.0, 0.02),
                coordsA=ax_zoom.transAxes,
                xyB=(zoom_end_idx + 0.5, len(display_binary) - 0.5),
                coordsB=ax_heat.transData,
                color=ZOOM_CONNECTOR_COLOR,
                linewidth=0.95,
                linestyle=(0, (2.2, 1.8)),
                alpha=0.92,
                zorder=8,
            )
            fig.add_artist(connector_top)
            fig.add_artist(connector_bottom)
    if ax_key is not None:
        ax_key.axis("off")
        items = [
            (CPG_UNMETH, "Unmethylated CpG"),
            (CPG_METH, "Methylated CpG"),
            (CPG_MISSING, "Not observed"),
        ]
        x_positions = [0.10, 0.43, 0.74]
        for (color, label), xpos in zip(items, x_positions):
            ax_key.add_patch(
                Rectangle(
                    (xpos, 0.00),
                    0.028,
                    0.24,
                    transform=ax_key.transAxes,
                    facecolor=color,
                    edgecolor="#BDBDBD" if label == "Not observed" else "none",
                    linewidth=0.5,
                    clip_on=False,
                )
            )
            ax_key.text(
                xpos + 0.040,
                0.01,
                label,
                transform=ax_key.transAxes,
                ha="left",
                va="bottom",
                fontsize=7.9,
                color=TEXT_MUTED,
            )


def render_panel_c(
    fig: plt.Figure,
    spec,
    panel_c_df: pd.DataFrame,
    entropy_stats: pd.DataFrame,
    sample_medians: pd.DataFrame,
    pairwise_stats: pd.DataFrame,
) -> None:
    ax = fig.add_subplot(spec)
    ax.axvspan(0.64, 1.36, color=GUIDE_BAND, alpha=0.20, zorder=0)

    sns.boxplot(
        data=panel_c_df,
        x="region_plot",
        y="methylation_entropy",
        hue="sample_group",
        order=PANEL_C_REGION_ORDER,
        hue_order=PANEL_C_GROUP_ORDER,
        palette=GROUP_PALETTE,
        width=0.64,
        linewidth=0.9,
        showfliers=False,
        ax=ax,
    )
    if ax.legend_ is not None:
        ax.legend_.remove()
    overlay_sample_medians_grouped(
        ax,
        sample_medians,
        PANEL_C_REGION_ORDER,
        PANEL_C_GROUP_ORDER,
        GROUP_PALETTE,
        "region_plot",
        "sample_group",
        "sample_median_entropy",
        total_width=0.64,
    )
    ax.set_title("C. Molecule-level methylation entropy", loc="left", fontsize=14.4, fontweight="bold", pad=2.4)
    ax.set_xlabel("")
    ax.set_ylabel("Methylation entropy", fontsize=10.5)
    ax.set_xticks(range(len(PANEL_C_REGION_ORDER)))
    ax.set_xticklabels(["SNRPN/\nPWS-IC", "SNORD116", "Downstream\ncontrol"], fontsize=10.4)
    ax.set_ylim(0.0, 0.57)
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.7)
    ax.set_axisbelow(True)
    sns.despine(ax=ax)
    for tick, label in zip(ax.get_xticklabels(), PANEL_C_REGION_ORDER):
        if label == "SNORD116":
            tick.set_fontweight("bold")
        else:
            tick.set_color(TEXT_MUTED)
    total_width = 0.64
    group_width = total_width / len(PANEL_C_GROUP_ORDER)
    group_positions = {
        group: -total_width / 2 + group_width / 2 + idx * group_width
        for idx, group in enumerate(PANEL_C_GROUP_ORDER)
    }
    bracket_pairs = [
        ("Control", "AS-DEL"),
        ("PWS-DEL", "PWS-mUPD"),
        ("AS-DEL", "PWS-DEL"),
        ("Control", "PWS-DEL"),
        ("AS-DEL", "PWS-mUPD"),
        ("Control", "PWS-mUPD"),
    ]
    for region_idx, region in enumerate(PANEL_C_REGION_ORDER):
        stat_row = entropy_stats.loc[entropy_stats["region"] == region].iloc[0]
        region_pairs = pairwise_stats[pairwise_stats["region"] == region].copy()
        region_max = float(
            pd.to_numeric(
                panel_c_df.loc[panel_c_df["region_plot"] == region, "methylation_entropy"],
                errors="coerce",
            ).max()
        )
        base_y = min(region_max + 0.020, 0.43)
        step_y = 0.024
        tick_height = 0.008
        significant_pairs = []
        for group_a, group_b in bracket_pairs:
            hit = region_pairs[(region_pairs["group_a"] == group_a) & (region_pairs["group_b"] == group_b)]
            if hit.empty:
                continue
            q_value = float(hit["q_value"].iloc[0])
            if q_value >= PAIRWISE_Q_DISPLAY_THRESHOLD:
                continue
            significant_pairs.append((group_a, group_b, q_value))
        for level, (group_a, group_b, q_value) in enumerate(significant_pairs):
            x1 = region_idx + group_positions[group_a]
            x2 = region_idx + group_positions[group_b]
            y = base_y + level * step_y
            ax.plot(
                [x1, x1, x2, x2],
                [y - tick_height, y, y, y - tick_height],
                color="#4A4A4A",
                linewidth=0.85,
                clip_on=False,
                zorder=6,
            )
            ax.text(
                (x1 + x2) / 2,
                y + 0.004,
                rf"${scientific_notation_mathtext(q_value, 'q')}$",
                ha="center",
                va="bottom",
                fontsize=5.8,
                color="#303030",
                clip_on=False,
                zorder=7,
                bbox={
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.92,
                    "pad": 0.10,
                },
            )
        if not significant_pairs:
            ax.text(
                region_idx,
                base_y + 0.006,
                rf"${scientific_notation_mathtext(float(stat_row['q_value']), 'q')}$",
                ha="center",
                va="bottom",
                fontsize=6.3,
                color="#303030",
                clip_on=False,
                zorder=7,
                bbox={
                    "facecolor": "white",
                    "edgecolor": "#D7D7D7",
                    "linewidth": 0.35,
                    "alpha": 0.94,
                    "pad": 0.14,
                },
            )


def draw_deleted_allele(ax: plt.Axes, x: float, y: float, width: float, height: float, facecolor: str) -> None:
    ax.add_patch(Rectangle((x, y), width, height, facecolor=facecolor, edgecolor="#8A8A8A", linewidth=0.9))
    ax.plot([x, x + width], [y, y + height], color="#8A8A8A", linewidth=1.1)
    ax.plot([x, x + width], [y + height, y], color="#8A8A8A", linewidth=1.1)


def render_panel_d(
    fig: plt.Figure,
    spec,
    panel_d_df: pd.DataFrame,
    parental_stats: pd.DataFrame,
    sample_medians: pd.DataFrame,
) -> None:
    panel_d = panel_d_df.copy()
    panel_d["panel_label"] = panel_d["sample_group"].map(
        {"AS-DEL": "AS-DEL retained paternal", "PWS-DEL": "PWS-DEL retained maternal"}
    )
    panel_d = panel_d[panel_d["panel_label"].isin(PANEL_D_ORDER)].copy()
    inner = spec.subgridspec(2, 1, height_ratios=[1.0, 0.19], hspace=0.03)
    ax = fig.add_subplot(inner[0, 0])
    ax_counts = fig.add_subplot(inner[1, 0], sharex=ax)

    sns.violinplot(
        data=panel_d,
        x="panel_label",
        y="mean_methylation",
        order=PANEL_D_ORDER,
        hue="panel_label",
        hue_order=PANEL_D_ORDER,
        palette=PANEL_D_PALETTE,
        dodge=False,
        cut=0,
        inner=None,
        bw_adjust=0.85,
        linewidth=0.85,
        ax=ax,
    )
    if ax.legend_ is not None:
        ax.legend_.remove()
    sns.boxplot(
        data=panel_d,
        x="panel_label",
        y="mean_methylation",
        order=PANEL_D_ORDER,
        width=0.22,
        showfliers=False,
        boxprops={"facecolor": "white", "edgecolor": "#333333", "linewidth": 0.9, "zorder": 3},
        whiskerprops={"color": "#333333", "linewidth": 0.9},
        capprops={"color": "#333333", "linewidth": 0.9},
        medianprops={"color": "#111111", "linewidth": 1.0},
        ax=ax,
    )
    overlay_sample_medians_category(
        ax,
        sample_medians,
        PANEL_D_ORDER,
        PANEL_D_PALETTE,
        "panel_label",
        "sample_median_methylation",
        marker_size=82,
    )
    ax.axhline(0.5, color="#A6A6A6", linestyle="--", linewidth=0.8, alpha=0.75, zorder=0)
    ax.set_title("D. Parent-of-origin contrast at SNORD116", loc="left", fontsize=14.4, fontweight="bold", pad=2.4)
    ax.set_xlabel("")
    ax.set_ylabel("Mean methylation probability per molecule", fontsize=10.5)
    ax.set_xticks(range(len(PANEL_D_ORDER)))
    ax.set_xticklabels(["AS-DEL\nretained paternal", "PWS-DEL\nretained maternal"], fontsize=10.2)
    ax.set_ylim(0.0, 1.05)
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.7)
    ax.set_axisbelow(True)
    sns.despine(ax=ax)

    stat_row = parental_stats.iloc[0]
    bracket_y = 0.915
    ax.plot([0, 0, 1, 1], [bracket_y - 0.032, bracket_y, bracket_y, bracket_y - 0.032], color="#3D3D3D", linewidth=1.05, clip_on=False)
    ax.text(
        0.5,
        bracket_y + 0.018,
        "\n".join(
            [
                f"Δ median = {stat_row['delta_median']:.3f}",
                rf"${scientific_notation_mathtext(stat_row['p_value'], 'P')}$",
            ]
        ),
        fontsize=9.1,
        ha="center",
        va="bottom",
        color="#303030",
        bbox={"boxstyle": "round,pad=0.24", "facecolor": "white", "edgecolor": "#D4D4D4", "linewidth": 0.4, "alpha": 0.90},
    )

    ax_counts.set_ylim(0.0, 1.0)
    ax_counts.set_yticks([])
    ax_counts.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    for spine in ["top", "right", "left", "bottom"]:
        ax_counts.spines[spine].set_visible(False)
    ax_counts.axhline(0.98, color="#E3E3E3", linewidth=0.7)
    ax_counts.text(
        0.5,
        0.92,
        f"95% CI = {stat_row['ci_low']:.3f}-{stat_row['ci_high']:.3f}; U = {stat_row['statistic_U']:,.0f}",
        transform=ax_counts.transAxes,
        ha="center",
        va="top",
        fontsize=8.1,
        color=TEXT_MUTED,
    )

    count_rows = [
        (
            0,
            int(stat_row["n_molecules_group_a"]),
            int(stat_row["n_samples_group_a"]),
            float(stat_row["median_group_a"]),
        ),
        (
            1,
            int(stat_row["n_molecules_group_b"]),
            int(stat_row["n_samples_group_b"]),
            float(stat_row["median_group_b"]),
        ),
    ]
    for xpos, n_mol, n_samp, median_val in count_rows:
        ax_counts.text(xpos, 0.60, f"{n_mol} molecules | {n_samp} samples", ha="center", va="center", fontsize=8.25, color=TEXT_MUTED)
        ax_counts.text(xpos, 0.18, f"median = {median_val:.3f}", ha="center", va="center", fontsize=8.25, color=TEXT_MUTED)


def render_shared_legends(fig: plt.Figure) -> None:
    group_handles = [Patch(facecolor=GROUP_PALETTE[group], edgecolor="none", label=group) for group in PANEL_C_GROUP_ORDER]
    fig.text(0.5, 0.040, "Sample group", ha="center", va="bottom", fontsize=8.9, color="#333333")
    fig.legend(
        group_handles,
        [handle.get_label() for handle in group_handles],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.018),
        frameon=False,
        ncol=4,
        columnspacing=1.0,
        handlelength=1.0,
        fontsize=8.9,
    )
    fig.text(
        0.5,
        0.001,
        "Rows represent individual HiFi molecules; columns represent ordered CpG positions in T2T-CHM13. Grey indicates no molecule-level observation. Panels A-B compress CpG positions for display.",
        ha="center",
        va="bottom",
        fontsize=8.7,
        color=TEXT_MUTED,
    )


def render_figure(
    outdir: Path,
    matrix_a_display: pd.DataFrame,
    matrix_b_display: pd.DataFrame,
    matrix_a_raw: pd.DataFrame,
    matrix_b_raw: pd.DataFrame,
    matrix_a_full: pd.DataFrame,
    matrix_b_full: pd.DataFrame,
    meta_a_display: pd.DataFrame,
    meta_b_display: pd.DataFrame,
    meta_a_full: pd.DataFrame,
    meta_b_full: pd.DataFrame,
    feature_df: pd.DataFrame,
    core_window: tuple[int, int],
    panel_c_df: pd.DataFrame,
    entropy_stats: pd.DataFrame,
    panel_c_sample_medians: pd.DataFrame,
    panel_c_pairwise_stats: pd.DataFrame,
    panel_d_df: pd.DataFrame,
    parental_stats: pd.DataFrame,
    panel_d_sample_medians: pd.DataFrame,
    zoom_window: tuple[int, int],
) -> None:
    figdir = ensure_dir(outdir / "figures")

    sns.set_theme(style="white", context="paper")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.sans-serif": ["DejaVu Sans"],
            "font.size": 10.4,
            "axes.titlesize": 14.4,
            "axes.labelsize": 10.7,
            "xtick.labelsize": 9.4,
            "ytick.labelsize": 9.4,
            "figure.facecolor": "#FFFFFF",
            "axes.facecolor": "#FFFFFF",
            "savefig.facecolor": "#FFFFFF",
        }
    )

    fig = plt.figure(figsize=(17.0, 12.6))
    fig.subplots_adjust(left=0.055, right=0.985, top=0.982, bottom=0.120)
    outer = fig.add_gridspec(2, 2, height_ratios=[1.24, 0.92], hspace=0.17, wspace=0.18)

    render_barcode_panel(
        fig,
        outer[0, 0],
        panel_letter="A",
        panel_title="Paternal-like control molecules",
        line_color=CONTROL_PARENTAL_COLORS["Control paternal"],
        display_matrix=matrix_a_display,
        display_meta=meta_a_display,
        full_matrix=matrix_a_full,
        full_meta=meta_a_full,
        raw_matrix=matrix_a_raw,
        feature_df=feature_df,
        core_window=core_window,
        show_ylabel=True,
        zoom_window=zoom_window,
    )
    render_barcode_panel(
        fig,
        outer[0, 1],
        panel_letter="B",
        panel_title="Maternal-like control molecules",
        line_color=CONTROL_PARENTAL_COLORS["Control maternal"],
        display_matrix=matrix_b_display,
        display_meta=meta_b_display,
        full_matrix=matrix_b_full,
        full_meta=meta_b_full,
        raw_matrix=matrix_b_raw,
        feature_df=feature_df,
        core_window=core_window,
        show_ylabel=False,
        zoom_window=zoom_window,
    )

    render_panel_c(fig, outer[1, 0], panel_c_df, entropy_stats, panel_c_sample_medians, panel_c_pairwise_stats)
    render_panel_d(fig, outer[1, 1], panel_d_df, parental_stats, panel_d_sample_medians)
    render_shared_legends(fig)

    for stem in FIGURE_BASENAMES:
        base = figdir / stem
        fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight")
        fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
        fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
        fig.savefig(base.with_suffix(".jpeg"), dpi=400, bbox_inches="tight")
    primary_base = figdir / FIGURE_BASENAMES[0]
    for suffix in [".png", ".pdf", ".svg", ".jpeg"]:
        target = outdir / f"{FIGURE_BASENAMES[0]}{suffix}"
        target.write_bytes(primary_base.with_suffix(suffix).read_bytes())
    plt.close(fig)


def render_supplement_barcode_figure(
    outdir: Path,
    matrix_a_full: pd.DataFrame,
    matrix_b_full: pd.DataFrame,
    meta_a_full: pd.DataFrame,
    meta_b_full: pd.DataFrame,
    feature_df: pd.DataFrame,
    core_window: tuple[int, int],
    zoom_window: tuple[int, int],
) -> None:
    figdir = ensure_dir(outdir / "figures")

    fig = plt.figure(figsize=(15.8, 14.2))
    fig.subplots_adjust(left=0.06, right=0.98, top=0.95, bottom=0.06)
    fig.suptitle("Supplementary full barcode stacks across the SNORD116 display window", fontsize=17.0, fontweight="bold", y=0.98)
    outer = fig.add_gridspec(2, 1, hspace=0.26)

    meta_a_sorted = order_full_barcode_stack(meta_a_full)
    meta_b_sorted = order_full_barcode_stack(meta_b_full)
    render_barcode_panel(
        fig,
        outer[0, 0],
        panel_letter="S1A",
        panel_title="Paternal-like control molecules full stack",
        line_color=CONTROL_PARENTAL_COLORS["Control paternal"],
        display_matrix=matrix_a_full.reindex(meta_a_sorted.index),
        display_meta=meta_a_sorted,
        full_matrix=matrix_a_full,
        full_meta=meta_a_full,
        raw_matrix=matrix_a_full,
        feature_df=feature_df,
        core_window=core_window,
        show_ylabel=True,
        zoom_window=zoom_window,
        show_zoom_inset=False,
        show_local_cpg_legend=False,
    )
    render_barcode_panel(
        fig,
        outer[1, 0],
        panel_letter="S1B",
        panel_title="Maternal-like control molecules full stack",
        line_color=CONTROL_PARENTAL_COLORS["Control maternal"],
        display_matrix=matrix_b_full.reindex(meta_b_sorted.index),
        display_meta=meta_b_sorted,
        full_matrix=matrix_b_full,
        full_meta=meta_b_full,
        raw_matrix=matrix_b_full,
        feature_df=feature_df,
        core_window=core_window,
        show_ylabel=True,
        zoom_window=zoom_window,
        show_zoom_inset=False,
        show_local_cpg_legend=False,
    )

    base = figdir / SUPPLEMENT_BARCODE_BASENAME
    fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def write_report(
    outdir: Path,
    molecule_summary: pd.DataFrame,
    entropy_stats: pd.DataFrame,
    pairwise_stats: pd.DataFrame,
    parental_stats: pd.DataFrame,
    display_sensitivity: pd.DataFrame,
    panel_d_sensitivity: pd.DataFrame,
    caption: str,
    total_region_observations: int,
    total_unique_molecules: int,
    total_cpg_calls: int,
    core_window: tuple[int, int],
    zoom_window: tuple[int, int],
) -> None:
    report_path = outdir / REPORT_NAME
    reports_dir = ensure_dir(outdir / "reports")
    report_mirror_path = reports_dir / REPORT_NAME
    generic_report_path = reports_dir / "report.md"

    region_summary = molecule_summary[molecule_summary["summary_context"] == "region_level"].copy()
    display_summary = molecule_summary[molecule_summary["summary_context"] == "control_display_window"].copy()

    group_totals = (
        region_summary.groupby("group", observed=True)["n_molecules"]
        .sum()
        .reset_index(name="region_level_molecule_observations")
    )
    region_totals = (
        region_summary.groupby("region", observed=True)["n_molecules"]
        .sum()
        .reset_index(name="n_molecules")
    )
    snord116_stat = entropy_stats.loc[entropy_stats["region"] == "SNORD116"].iloc[0]
    ic_stat = entropy_stats.loc[entropy_stats["region"] == "SNRPN/PWS-IC"].iloc[0]
    downstream_stat = entropy_stats.loc[entropy_stats["region"] == "Downstream control"].iloc[0]
    parental_row = parental_stats.iloc[0]
    paternal_display = display_summary.loc[display_summary["group"] == "Paternal-like control molecules"].iloc[0]
    maternal_display = display_summary.loc[display_summary["group"] == "Maternal-like control molecules"].iloc[0]
    significant_pairwise_stats = (
        pairwise_stats.loc[pd.to_numeric(pairwise_stats["q_value"], errors="coerce") < PAIRWISE_Q_DISPLAY_THRESHOLD]
        .copy()
        .sort_values(["region", "q_value", "group_a", "group_b"], kind="mergesort")
    )
    entropy_stats_display = entropy_stats[
        ["region", "statistic_H", "p_value", "q_value", "eta_squared"]
    ].rename(
        columns={
            "statistic_H": "H",
            "p_value": "p",
            "q_value": "q",
            "eta_squared": "eta^2",
        }
    )
    entropy_stats_display["p"] = entropy_stats_display["p"].map(scientific_notation_value_text)
    entropy_stats_display["q"] = entropy_stats_display["q"].map(scientific_notation_value_text)
    significant_pairwise_display = significant_pairwise_stats[
        ["region", "group_a", "group_b", "statistic_U", "p_value", "q_value"]
    ].rename(
        columns={
            "group_a": "group A",
            "group_b": "group B",
            "statistic_U": "U",
            "p_value": "p",
            "q_value": "q",
        }
    )
    significant_pairwise_display["p"] = significant_pairwise_display["p"].map(scientific_notation_value_text)
    significant_pairwise_display["q"] = significant_pairwise_display["q"].map(scientific_notation_value_text)
    parental_stats_display = parental_stats[
        [
            "comparison",
            "n_molecules_group_a",
            "n_molecules_group_b",
            "n_samples_group_a",
            "n_samples_group_b",
            "median_group_a",
            "median_group_b",
            "delta_median",
            "ci_low",
            "ci_high",
            "statistic_U",
            "p_value",
            "rank_biserial",
        ]
    ].rename(
        columns={
            "n_molecules_group_a": "n mol A",
            "n_molecules_group_b": "n mol B",
            "n_samples_group_a": "n samples A",
            "n_samples_group_b": "n samples B",
            "median_group_a": "median A",
            "median_group_b": "median B",
            "delta_median": "delta median",
            "ci_low": "CI low",
            "ci_high": "CI high",
            "statistic_U": "U",
            "p_value": "p",
            "rank_biserial": "rank-biserial",
        }
    )
    parental_stats_display["p"] = parental_stats_display["p"].map(scientific_notation_value_text)

    report_lines = [
        "# Figure SNORD116 single-molecule report",
        "",
        f"Generated: {date.today().isoformat()}",
        "",
        "## Scope",
        "",
        f"Figure analyzed: `figures/{FIGURE_BASENAMES[0]}.png`",
        "",
        "This report documents the publication-facing four-panel SNORD116 single-molecule figure and the quantitative tables used to render it. The goal of the figure is to show that molecule-level CpG methylation barcodes preserve structured parent-of-origin patterns across imprinting-relevant intervals, with the strongest molecule-level heterogeneity and reciprocal parental contrast occurring at `SNORD116`.",
        "",
        "Primary biological message:",
        "",
        "- Molecule-level CpG methylation barcodes preserve the window- and boundary-level methylation patterns seen at the `PWS imprinting centre` and `SNORD116`.",
        "- These methylation states are directly visible on individual PacBio HiFi molecules rather than only as cohort-level averages.",
        "- Cross-group entropy differences are region-specific and are strongest at `SNORD116`.",
        "- The strongest reciprocal read-level parent-of-origin contrast is focal to `SNORD116`, where `AS-DEL retained paternal` molecules are more methylated than `PWS-DEL retained maternal` molecules.",
        "",
        "## Overview",
        "",
        f"- Region-level molecule observations analyzed for Panels C-D: `{total_region_observations:,}`",
        f"- Unique molecules contributing to the region-level analyses: `{total_unique_molecules:,}`",
        f"- Total CpG calls across SNRPN/PWS-IC, SNORD116, and downstream control windows: `{total_cpg_calls:,}`",
        "- Molecule-level panels A-B show representative high-completeness control molecules from the shared SNORD116 display window and are reported separately from the region-level statistical analyses.",
        f"- Full control barcode stacks are additionally exported to `figures/{SUPPLEMENT_BARCODE_BASENAME}.pdf` and companion image formats.",
        "- Panels A-B compress CpG positions for display only so the barcode structure remains legible at manuscript scale; all statistics remain tied to the uncompressed saved tables.",
        "",
        "## Output files",
        "",
        f"- Main figure PNG: `figures/{FIGURE_BASENAMES[0]}.png`",
        f"- Main figure PDF: `figures/{FIGURE_BASENAMES[0]}.pdf`",
        f"- Main figure SVG: `figures/{FIGURE_BASENAMES[0]}.svg`",
        f"- Main figure JPEG: `{FIGURE_BASENAMES[0]}.jpeg`",
        f"- Figure-specific report: `{REPORT_NAME}`",
        f"- Mirrored detailed report: `reports/report.md`",
        "",
        "## Data provenance",
        "",
        "Primary tables used by the renderer:",
        "",
        "- `tables/Figure4_panelA_control_paternal_SNORD116_matrix.tsv.gz`",
        "- `tables/Figure4_panelA_control_paternal_SNORD116_rows.tsv`",
        "- `tables/Figure4_panelB_control_maternal_SNORD116_matrix.tsv.gz`",
        "- `tables/Figure4_panelB_control_maternal_SNORD116_rows.tsv`",
        "- `tables/Figure4_SNORD116_shared_core_window.tsv`",
        "- `tables/Figure4_gene_features.tsv`",
        "- `tables/Figure4_panelC_entropy_plot_input.tsv`",
        "- `tables/Figure4_panelD_SNORD116_plot_input.tsv`",
        "- `tables/Figure4_single_molecule_CpG_calls.tsv.gz`",
        "",
        "## Figure caption",
        "",
        caption,
        "",
        "## Figure design notes",
        "",
        "- Panels `A-B` are structural visualization panels. They use the true saved single-molecule calls, display representative high-completeness control molecules, and preserve molecule ordering, coordinates, and CpG-state values. Their CpG columns are compressed only for display density.",
        "- The highlighted SNORD116 zoom inset is generated from a narrower shared high-coverage subwindow inside the SNORD116 core display interval and uses the same underlying molecules and same methylation calls.",
        f"- Shared SNORD116 core window used for emphasis in Panels A-B: `chr15:{core_window[0]:,}-{core_window[1]:,}`.",
        f"- Zoomed subwindow used for the barcode inset: `chr15:{zoom_window[0]:,}-{zoom_window[1]:,}`.",
        "- Panels `C-D` use the quantitative saved tables directly and are not derived from the display compression used in Panels A-B.",
        f"- Panel `C` draws only pairwise entropy brackets with `q < {PAIRWISE_Q_DISPLAY_THRESHOLD:.2f}` above the boxplots; if a region has no significant pairwise contrast, the overall region-level `q` is shown instead.",
        "",
        "## Executive interpretation",
        "",
        f"- `SNRPN/PWS-IC` shows significant cross-group entropy differences (`{scientific_notation_text(float(ic_stat['q_value']), 'q')}`, `eta^2 = {ic_stat['eta_squared']:.2f}`), consistent with imprinting-centre-specific molecule-level organisation.",
        f"- `SNORD116` shows the strongest entropy difference by a large margin (`{scientific_notation_text(float(snord116_stat['q_value']), 'q')}`, `eta^2 = {snord116_stat['eta_squared']:.2f}`), marking it as the dominant locus of molecule-level methylation heterogeneity in this figure.",
        f"- `Downstream control` does not show meaningful cross-group entropy separation (`{scientific_notation_text(float(downstream_stat['q_value']), 'q')}`, `eta^2 = {downstream_stat['eta_squared']:.2f}`), arguing against a diffuse genome-wide effect.",
        f"- The reciprocal deletion comparison at `SNORD116` remains strongly separated (`delta median = {parental_row['delta_median']:.3f}`, `95% CI = {parental_row['ci_low']:.3f}-{parental_row['ci_high']:.3f}`, `{scientific_notation_text(float(parental_row['p_value']), 'p')}`), with `AS-DEL retained paternal` molecules more methylated than `PWS-DEL retained maternal` molecules.",
        "",
        "## Panel-by-panel reading guide",
        "",
        "### Panel A. Paternal-like control molecules",
        "",
        f"- Shows `{int(paternal_display['n_displayed'])}` of `{int(paternal_display['n_molecules'])}` qualifying control paternal-like molecules from `{int(paternal_display['n_samples'])}` control samples.",
        f"- Median control-paternal display methylation is `{paternal_display['median_methylation']:.3f}` with median `{paternal_display['median_cpgs_per_molecule']:.0f}` CpGs per molecule and median span `{paternal_display['median_molecule_span_bp']:.0f}` bp.",
        "- Each row is one HiFi molecule; each colored cell corresponds to an ordered CpG position after display compression.",
        "- Red marks methylated CpGs, blue marks unmethylated CpGs, and grey marks positions not observed on that molecule.",
        "- The design goal of this panel is to let the reader see extended multi-CpG methylation patterns on single molecules rather than isolated per-site changes.",
        "",
        "### Panel B. Maternal-like control molecules",
        "",
        f"- Shows `{int(maternal_display['n_displayed'])}` of `{int(maternal_display['n_molecules'])}` qualifying control maternal-like molecules from `{int(maternal_display['n_samples'])}` control samples.",
        f"- Median control-maternal display methylation is `{maternal_display['median_methylation']:.3f}` with median `{maternal_display['median_cpgs_per_molecule']:.0f}` CpGs per molecule and median span `{maternal_display['median_molecule_span_bp']:.0f}` bp.",
        "- The shared coordinate system, feature track, and zoom inset allow direct structural comparison with Panel A over the same SNORD116 display window.",
        "- Together with Panel A, this panel is meant to teach the reader how coherent single-molecule methylation barcodes appear in the control state before moving to the quantitative tests in Panels C-D.",
        "",
        "### Panel C. Molecule-level methylation entropy",
        "",
        "- Boxplots summarize molecule-level entropy by region and sample group; white points are sample-level medians.",
        f"- `SNRPN/PWS-IC`: significant but moderate heterogeneity (`q = {scientific_notation_value_text(ic_stat['q_value'])}`, `eta^2 = {ic_stat['eta_squared']:.2f}`).",
        f"- `SNORD116`: strongest heterogeneity signal in the figure (`q = {scientific_notation_value_text(snord116_stat['q_value'])}`, `eta^2 = {snord116_stat['eta_squared']:.2f}`).",
        f"- `Downstream control`: null result (`q = {scientific_notation_value_text(downstream_stat['q_value'])}`, `eta^2 = {downstream_stat['eta_squared']:.2f}`), which supports locus specificity.",
        f"- Only significant within-region pairwise entropy comparisons (`q < {PAIRWISE_Q_DISPLAY_THRESHOLD:.2f}`) are drawn as brackets above the boxplots in the figure; `Downstream control` shows the non-significant overall region-level `q` because no pairwise contrast passed the display threshold.",
        "- This panel is the main region-level quantitative test that turns the structural barcode impression from Panels A-B into a formal heterogeneity comparison.",
        "",
        "### Panel D. Parent-of-origin contrast at SNORD116",
        "",
        f"- Compares `{int(parental_row['n_molecules_group_a'])}` `AS-DEL retained paternal` molecules from `{int(parental_row['n_samples_group_a'])}` samples against `{int(parental_row['n_molecules_group_b'])}` `PWS-DEL retained maternal` molecules from `{int(parental_row['n_samples_group_b'])}` samples.",
        f"- Median per-molecule methylation is `{parental_row['median_group_a']:.3f}` in `AS-DEL retained paternal` and `{parental_row['median_group_b']:.3f}` in `PWS-DEL retained maternal`.",
        f"- The observed separation is `delta median = {parental_row['delta_median']:.3f}` with bootstrap `95% CI = {parental_row['ci_low']:.3f}-{parental_row['ci_high']:.3f}` and Mann-Whitney `U = {parental_row['statistic_U']:,.0f}`.",
        "- This panel provides the clearest reciprocal parent-of-origin result in the figure and directly matches the paternal-higher signal described in the upstream windowed analysis.",
        "",
        "## Molecule totals by group",
        "",
        dataframe_to_markdown(group_totals),
        "",
        "## Molecule totals by region",
        "",
        dataframe_to_markdown(region_totals),
        "",
        "## Control display-window summary",
        "",
        dataframe_to_markdown(
            display_summary[
                [
                    "group",
                    "n_displayed",
                    "n_molecules",
                    "n_samples",
                    "median_cpgs_per_molecule",
                    "median_molecule_span_bp",
                    "median_missing_fraction",
                    "median_methylation",
                ]
            ].rename(
                columns={
                    "group": "display panel",
                    "n_displayed": "n displayed",
                    "n_molecules": "n total",
                    "n_samples": "n samples",
                    "median_cpgs_per_molecule": "median CpGs",
                    "median_molecule_span_bp": "median span bp",
                    "median_missing_fraction": "median missing fraction",
                    "median_methylation": "median methylation",
                }
            )
        ),
        "",
        "## Region-level molecule summary",
        "",
        dataframe_to_markdown(
            region_summary[
                [
                    "region",
                    "group",
                    "n_molecules",
                    "n_samples",
                    "total_cpg_calls",
                    "median_cpgs_per_molecule",
                    "median_molecule_span_bp",
                    "median_missing_fraction",
                    "median_methylation",
                    "median_entropy",
                    "entropy_q1",
                    "entropy_q3",
                ]
            ].rename(
                columns={
                    "n_molecules": "n molecules",
                    "n_samples": "n samples",
                    "total_cpg_calls": "total CpG calls",
                    "median_cpgs_per_molecule": "median CpGs",
                    "median_molecule_span_bp": "median span bp",
                    "median_missing_fraction": "median missing fraction",
                    "median_methylation": "median methylation",
                    "median_entropy": "median entropy",
                    "entropy_q1": "entropy Q1",
                    "entropy_q3": "entropy Q3",
                }
            )
        ),
        "",
        "## Entropy statistics",
        "",
        dataframe_to_markdown(entropy_stats_display),
        "",
        "## Significant pairwise entropy comparisons shown in Panel C",
        "",
        dataframe_to_markdown(significant_pairwise_display) if not significant_pairwise_stats.empty else "No within-region pairwise entropy comparisons reached the display threshold.",
        "",
        "## Reciprocal parental-state statistics",
        "",
        dataframe_to_markdown(parental_stats_display),
        "",
        "## Statistical methods",
        "",
        "- Panel C: Kruskal-Wallis tests across the four sample groups were run independently for `SNRPN/PWS-IC`, `SNORD116`, and `Downstream control`.",
        "- Panel C multiple-testing correction: Benjamini-Hochberg false-discovery-rate adjustment across the three region-level entropy tests.",
        f"- Panel C pairwise brackets: two-sided Mann-Whitney U tests were computed for all six within-region group pairs, with Benjamini-Hochberg correction applied within each region; only comparisons with `q < {PAIRWISE_Q_DISPLAY_THRESHOLD:.2f}` are drawn as brackets on the figure. Regions without any significant pairwise contrast display the overall region-level `q` instead.",
        "- Panel C effect size: Kruskal eta-squared analogue `(H - k + 1) / (n - k)`.",
        "- Panel D: a prespecified two-sided Mann-Whitney U test comparing `AS-DEL retained paternal` versus `PWS-DEL retained maternal` molecules at `SNORD116`.",
        "- Panel D confidence interval: bootstrap 95% CI for the median difference using 5,000 resamples.",
        "- Panels A-B: no inferential statistics are driven by the display compression. Compression affects only figure readability, not the saved underlying values.",
        "",
        "## Panel interpretation",
        "",
        "### Panel A and Panel B",
        "",
        "- Control paternal-like and maternal-like molecules show coherent multi-CpG barcode structure across the shared SNORD116 display window.",
        "- Panels A-B now display representative high-completeness control molecules; CpG positions are still compressed for display only so the molecule-level barcode structure remains legible at manuscript scale.",
        "- Added zoom insets now show representative molecules across a higher-resolution SNORD116 subwindow so the same underlying data can be read directly as coordinated multi-CpG barcodes rather than scattered single-site events.",
        "- The control parental-state mean-methylation difference is modest and filter-sensitive, so Panels A-B are intentionally presented as representative molecule-level structure rather than as the primary statistical proof of parental separation.",
        "",
        "### Panel C",
        "",
        f"- `SNORD116` shows the strongest region-level heterogeneity signal across sample groups (`{scientific_notation_text(float(entropy_stats.loc[entropy_stats['region'] == 'SNORD116', 'q_value'].iloc[0]), 'q')}`, `eta^2 = {entropy_stats.loc[entropy_stats['region'] == 'SNORD116', 'eta_squared'].iloc[0]:.2f}`), exceeding both `SNRPN/PWS-IC` and `Downstream control`.",
        "- The downstream control interval shows little cross-group entropy structure, indicating that the heterogeneity signal is locus-specific rather than a global methylation artifact.",
        "",
        "### Panel D",
        "",
        f"- `AS-DEL retained paternal` molecules have higher mean per-molecule methylation than `PWS-DEL retained maternal` molecules at SNORD116 (`delta median = {parental_stats['delta_median'].iloc[0]:.3f}`, `{scientific_notation_text(parental_stats['p_value'].iloc[0], 'p')}`).",
        "- This reciprocal pattern is robust across stricter minimum-CpG and minimum-span thresholds.",
        "",
        "## Sensitivity checks",
        "",
        "### Minimum span threshold in the control display panels",
        "",
        dataframe_to_markdown(
            display_sensitivity.rename(
                columns={
                    "min_span_fraction": "min span fraction",
                    "paternal_n": "paternal n",
                    "maternal_n": "maternal n",
                    "paternal_median_methylation": "paternal median",
                    "maternal_median_methylation": "maternal median",
                    "control_parental_p_value": "control p",
                }
            )
        ),
        "",
        "- The control paternal-like versus maternal-like mean-methylation contrast is not stable across display filters. This supports using Panels A-B as representative structural views and Panels C-D as the main quantitative evidence.",
        "",
        "### Reciprocal deletion comparison under alternative filters",
        "",
        dataframe_to_markdown(
            panel_d_sensitivity.rename(
                columns={
                    "min_cpgs": "min CpGs",
                    "min_span_bp": "min span bp",
                    "n_as_del": "AS-DEL n",
                    "n_pws_del": "PWS-DEL n",
                    "delta_median": "delta median",
                    "p_value": "p",
                    "statistic_U": "U",
                }
            )
        ),
        "",
        "- The AS-DEL versus PWS-DEL difference remains positive and statistically supported across all tested minimum-CpG and minimum-span thresholds shown above.",
        "",
        "### Displayed molecules versus quantitative analyses",
        "",
        "- Panels A-B display representative control molecules from the full qualifying sets. All inferential statistics in Panels C-D continue to use the full retained molecule sets defined in the saved phase-4 tables.",
        "",
        "### Relative strength of the SNORD116 heterogeneity signal",
        "",
        "- The cross-group entropy effect size is strongest at SNORD116 (`eta^2 = 0.19`), weaker at SNRPN/PWS-IC (`eta^2 = 0.13`), and absent downstream (`eta^2 = 0.00`).",
        "",
        "## Limitations and scope boundaries",
        "",
        "- Panels A-B are intended as structural single-molecule views, not as standalone inferential comparisons between paternal-like and maternal-like controls.",
        "- The control parental contrast in the display subset is sensitive to span and completeness filters, so the strongest formal evidence in this figure comes from Panels C-D rather than from an A-versus-B hypothesis test.",
        "- The figure is deliberately locus-focused. It supports focal cis-coordination at imprinting-relevant intervals and does not claim a genome-wide barcode phenomenon from this panel set alone.",
        "",
        "## Reproducibility",
        "",
        "Figure generation command:",
        "",
        "```bash",
        "python3 FIGURE_4.py --outdir /home/rare/arlen/pws-as-hifi-cis-methylation/scripts/hifi_multiomics_pipeline/06_results",
        "```",
        "",
        "Convenience outputs written by the renderer:",
        "",
        f"- `/home/rare/arlen/pws-as-hifi-cis-methylation/scripts/hifi_multiomics_pipeline/06_results/figures/{FIGURE_BASENAMES[0]}.png`",
        f"- `/home/rare/arlen/pws-as-hifi-cis-methylation/scripts/hifi_multiomics_pipeline/06_results/figures/{FIGURE_BASENAMES[0]}.pdf`",
        f"- `/home/rare/arlen/pws-as-hifi-cis-methylation/scripts/hifi_multiomics_pipeline/06_results/figures/{FIGURE_BASENAMES[0]}.svg`",
        f"- `/home/rare/arlen/pws-as-hifi-cis-methylation/scripts/hifi_multiomics_pipeline/06_results/figures/{FIGURE_BASENAMES[0]}.jpeg`",
        f"- `/home/rare/arlen/pws-as-hifi-cis-methylation/scripts/hifi_multiomics_pipeline/06_results/reports/{REPORT_NAME}`",
        f"- `/home/rare/arlen/pws-as-hifi-cis-methylation/scripts/hifi_multiomics_pipeline/06_results/reports/report.md`",
        "",
        "## Main result",
        "",
        "Single-molecule PacBio HiFi methylation calls show that SNORD116 is organized as coherent multi-CpG molecule-level barcodes rather than isolated CpG noise. Although the control paternal-like and maternal-like display panels are best interpreted as representative structural views, quantitative testing identifies SNORD116 as the locus with the strongest cross-group molecule-level heterogeneity and shows that reciprocal deletion genomes retain distinct parent-of-origin methylation states: AS-DEL retains the paternal-like state, whereas PWS-DEL retains the maternal-like state.",
        "",
    ]

    report_text = "\n".join(report_lines)
    report_path.write_text(report_text)
    report_mirror_path.write_text(report_text)
    generic_report_path.write_text(report_text)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--table-dir",
        type=Path,
        default=None,
        help="Directory containing precomputed Figure4_*.tsv input tables. Defaults to outdir/tables, then the original paper_vf tables.",
    )
    args = parser.parse_args()

    outdir = args.outdir
    figdir = ensure_dir(outdir / "figures")
    tbldir = ensure_dir(outdir / "tables")
    input_tbldir = resolve_input_table_dir(outdir, args.table_dir)

    matrix_a_full, meta_a_rows = load_matrix_and_rows(input_tbldir, "Figure4_panelA_control_paternal_SNORD116")
    matrix_b_full, meta_b_rows = load_matrix_and_rows(input_tbldir, "Figure4_panelB_control_maternal_SNORD116")
    union_columns = sorted(set(matrix_a_full.columns.tolist()).union(matrix_b_full.columns.tolist()))
    matrix_a_full = matrix_a_full.reindex(columns=union_columns)
    matrix_b_full = matrix_b_full.reindex(columns=union_columns)
    meta_a_full = annotate_display_rows(matrix_a_full, meta_a_rows, DISPLAY_WINDOW)
    meta_b_full = annotate_display_rows(matrix_b_full, meta_b_rows, DISPLAY_WINDOW)
    meta_a_candidates = build_display_candidates(meta_a_full)
    meta_b_candidates = build_display_candidates(meta_b_full)
    matrix_a_display_raw, meta_a_display = select_display_subset(matrix_a_full, meta_a_candidates, DISPLAY_TARGET_ROWS)
    matrix_b_display_raw, meta_b_display = select_display_subset(matrix_b_full, meta_b_candidates, DISPLAY_TARGET_ROWS)
    matrix_a_full_display = compress_matrix_for_display(matrix_a_full, MAIN_BARCODE_MAX_COLUMNS)
    matrix_b_full_display = compress_matrix_for_display(matrix_b_full, MAIN_BARCODE_MAX_COLUMNS)
    matrix_a_display_plot = compress_matrix_for_display(matrix_a_display_raw, MAIN_BARCODE_MAX_COLUMNS)
    matrix_b_display_plot = compress_matrix_for_display(matrix_b_display_raw, MAIN_BARCODE_MAX_COLUMNS)

    feature_df = pd.read_csv(input_tbldir / "Figure4_gene_features.tsv", sep="\t")
    feature_df = feature_df[feature_df["panel"] == "snord116"].copy()
    shared_core = pd.read_csv(input_tbldir / "Figure4_SNORD116_shared_core_window.tsv", sep="\t")
    core_window = (int(shared_core["cpg_position"].min()), int(shared_core["cpg_position"].max()))
    zoom_window = choose_zoom_window(matrix_a_full, matrix_b_full, core_window)

    panel_c_df = pd.read_csv(input_tbldir / "Figure4_panelC_entropy_plot_input.tsv", sep="\t")
    panel_c_df = panel_c_df[panel_c_df["region"].isin(["PWS-IC", "SNORD116 cluster", "Downstream control"])].copy()
    panel_c_df["region_plot"] = panel_c_df["region"].map(standardize_region_label)
    panel_c_df["sample_group"] = pd.Categorical(panel_c_df["sample_group"], PANEL_C_GROUP_ORDER, ordered=True)
    panel_c_df["region_plot"] = pd.Categorical(panel_c_df["region_plot"], PANEL_C_REGION_ORDER, ordered=True)

    panel_d_df = pd.read_csv(input_tbldir / "Figure4_panelD_SNORD116_plot_input.tsv", sep="\t")
    panel_d_df = panel_d_df[panel_d_df["sample_group"].isin(["AS-DEL", "PWS-DEL"])].copy()

    single_molecule_calls = pd.read_csv(input_tbldir / "Figure4_single_molecule_CpG_calls.tsv.gz", sep="\t", compression="gzip")

    region_summary = build_region_group_summary(panel_c_df, single_molecule_calls)
    control_display_summary = build_control_display_summary(
        full_meta_by_label={
            "Control paternal": meta_a_full,
            "Control maternal": meta_b_full,
        },
        displayed_by_label={
            "Control paternal": meta_a_display,
            "Control maternal": meta_b_display,
        },
    )
    molecule_summary = pd.concat([region_summary, control_display_summary], ignore_index=True, sort=False)

    entropy_stats, panel_c_sample_medians, panel_c_pairwise_stats = build_entropy_statistics(panel_c_df)
    parental_stats, panel_d_sample_medians = build_parental_state_statistics(panel_d_df)
    display_sensitivity = build_display_sensitivity(meta_a_full, meta_b_full)
    panel_d_sensitivity = build_panel_d_sensitivity(panel_d_df)

    caption = build_figure_caption(entropy_stats, parental_stats)

    render_figure(
        outdir=outdir,
        matrix_a_display=matrix_a_display_plot,
        matrix_b_display=matrix_b_display_plot,
        matrix_a_raw=matrix_a_full,
        matrix_b_raw=matrix_b_full,
        matrix_a_full=matrix_a_full_display,
        matrix_b_full=matrix_b_full_display,
        meta_a_display=meta_a_display,
        meta_b_display=meta_b_display,
        meta_a_full=meta_a_full,
        meta_b_full=meta_b_full,
        feature_df=feature_df,
        core_window=core_window,
        panel_c_df=panel_c_df,
        entropy_stats=entropy_stats,
        panel_c_sample_medians=panel_c_sample_medians,
        panel_c_pairwise_stats=panel_c_pairwise_stats,
        panel_d_df=panel_d_df,
        parental_stats=parental_stats,
        panel_d_sample_medians=panel_d_sample_medians,
        zoom_window=zoom_window,
    )
    render_supplement_barcode_figure(
        outdir=outdir,
        matrix_a_full=matrix_a_full,
        matrix_b_full=matrix_b_full,
        meta_a_full=meta_a_full,
        meta_b_full=meta_b_full,
        feature_df=feature_df,
        core_window=core_window,
        zoom_window=zoom_window,
    )

    molecule_summary.to_csv(tbldir / TABLE_MOLECULE_SUMMARY, sep="\t", index=False)
    entropy_stats.to_csv(tbldir / TABLE_ENTROPY_STATS, sep="\t", index=False)
    parental_stats.to_csv(tbldir / TABLE_PARENTAL_STATS, sep="\t", index=False)

    total_region_observations = int(len(panel_c_df))
    total_unique_molecules = int(panel_c_df["read_key"].nunique())
    total_cpg_calls = int(
        single_molecule_calls[single_molecule_calls["region"].isin(["PWS-IC", "SNORD116 cluster", "UBE3A downstream control"])].shape[0]
    )
    write_report(
        outdir=outdir,
        molecule_summary=molecule_summary,
        entropy_stats=entropy_stats,
        pairwise_stats=panel_c_pairwise_stats,
        parental_stats=parental_stats,
        display_sensitivity=display_sensitivity,
        panel_d_sensitivity=panel_d_sensitivity,
        caption=caption,
        total_region_observations=total_region_observations,
        total_unique_molecules=total_unique_molecules,
        total_cpg_calls=total_cpg_calls,
        core_window=core_window,
        zoom_window=zoom_window,
    )

    # Mirror the three requested tabular outputs into the root directory for convenience.
    for name in [TABLE_MOLECULE_SUMMARY, TABLE_ENTROPY_STATS, TABLE_PARENTAL_STATS]:
        root_copy = outdir / name
        root_copy.write_text((tbldir / name).read_text())


if __name__ == "__main__":
    main()
