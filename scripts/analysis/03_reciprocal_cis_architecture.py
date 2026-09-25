#!/usr/bin/env python3
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from cis_analysis import load_cohort


# ============================== CONFIGURATION ==============================
PROJECT_ROOT = Path(__file__).resolve().parents[2]
METADATA_PATH = PROJECT_ROOT / "assets" / "metadata.csv"
EVIDENCE_DIR = PROJECT_ROOT / "results" / "analysis" / "01_evidence_matrix"
EVIDENCE_MATRIX_PATH = EVIDENCE_DIR / "chr15_window_evidence_matrix.tsv.gz"
TRACK_INVENTORY_PATH = EVIDENCE_DIR / "methylation_track_inventory.tsv"
COMMON_CN1_PATH = EVIDENCE_DIR / "common_reciprocal_cn1_core.tsv"
REGIONS_PATH = PROJECT_ROOT / "results" / "analysis" / "prespecified_regions.tsv"
OUTPUT_DIR = PROJECT_ROOT / "results" / "analysis" / "03_cis_architecture"

PWS_GROUP = "PWS-DEL"
AS_GROUP = "AS-DEL"
DIPLOID_GROUPS = ("Control", "DiGeorge")
REQUIRED_REGIONS = (
    "MAGEL2/NDN",
    "PWS/AS imprinting centre",
    "SNRPN/SNHG14",
    "SNORD116",
    "UBE3A",
    "GABRB3/GABA receptor cluster",
    "OCA2 downstream control",
)

ESTIMATORS = {
    "full_depth": ("beta_site_mean", "n_cpg", "effective_obs"),
    "common_depth_downsampled": ("beta_downsampled", "n_cpg_downsampled", "effective_obs_downsampled"),
    "capped_effective_coverage": ("beta_depth_cap_15", "n_cpg", "effective_obs"),
    "shared_cpg": ("beta_shared_cpg", "n_cpg_shared", "effective_obs_shared"),
}
PRIMARY_ESTIMATOR = "full_depth"
COVERAGE_ESTIMATORS = ("common_depth_downsampled", "capped_effective_coverage")
SHARED_ESTIMATORS = ("shared_cpg",)

MIN_CPGS = 3
MIN_EFFECTIVE_OBS = 12.0
ESTIMATOR_WINDOW_MINIMUMS = {"shared_cpg": (1, 4.0)}
MIN_PWS_PARTICIPANTS = 3
MIN_AS_PARTICIPANTS = 2
MIN_REGION_WINDOW_FRACTION = 0.25
MIN_REGION_CPGS = 10
MIN_ASM_PARTICIPANTS_FOR_INTERVAL = 3

WINDOW_BOOTSTRAP_REPLICATES = 1_000
REGION_BOOTSTRAP_REPLICATES = 10_000
BOOTSTRAP_MIN_FINITE_FRACTION = 0.95
RANDOM_SEED = 20260924

ROLLING_WINDOWS = 21
ROLLING_MIN_WINDOWS = 21

FOCAL_MIN_ABS_DELTA = 0.15
FOCAL_MIN_WINDOWS = 10
REGION_MIN_ABS_DELTA = 0.05
SENSITIVITY_RETENTION = 0.5
SENSITIVITY_ABS_TOLERANCE = 0.05
# ===========================================================================

MATRIX_COLUMNS = {
    "sample_id",
    "mechanism",
    "track_kind",
    "window_id",
    "start",
    "end",
    "mid",
    "inside_common_cn1",
    "parental_direction_observed",
    "parental_state",
    "evidence_class",
} | {column for spec in ESTIMATORS.values() for column in spec}
INVENTORY_COLUMNS = {"sample_id", "mechanism", "track_kind", "status", "reference_depth_common_cn1"}
CORE_COLUMNS = {"chrom", "start", "end"}
REGION_COLUMNS = {"region_id", "chrom", "start", "end"}


def write_table(table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if path.suffix == ".gz" else None
    table.to_csv(path, sep="\t", index=False, na_rep="NA", compression=compression)


def read_validated(path: Path, required: set[str], label: str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")
    table = pd.read_csv(path, sep="\t", low_memory=False)
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"{label} lacks columns: {sorted(missing)}")
    if table.empty:
        raise ValueError(f"{label} is empty: {path}")
    return table


def nanmean(values: np.ndarray, axis: int) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(values, axis=axis)


def nanmedian(values: np.ndarray, axis: int) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmedian(values, axis=axis)


def percentile_interval(samples: np.ndarray, axis: int = -1) -> tuple[np.ndarray, np.ndarray]:
    finite = np.isfinite(samples).mean(axis=axis)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        low, high = np.nanquantile(samples, [0.025, 0.975], axis=axis)
    usable = finite >= BOOTSTRAP_MIN_FINITE_FRACTION
    return np.where(usable, low, np.nan), np.where(usable, high, np.nan)


def participant_matrices(
    table: pd.DataFrame,
    windows: pd.DataFrame,
    samples: list[str],
    estimator: str,
) -> dict[str, np.ndarray]:
    value_col, n_col, eff_col = ESTIMATORS[estimator]
    min_cpgs, min_effective = ESTIMATOR_WINDOW_MINIMUMS.get(estimator, (MIN_CPGS, MIN_EFFECTIVE_OBS))
    valid = (
        table[n_col].ge(min_cpgs)
        & table[eff_col].ge(min_effective)
        & table[value_col].notna()
    )
    frame = table.assign(
        value=table[value_col].where(valid),
        support=table[n_col].where(valid),
        effective=table[eff_col].where(valid),
    )
    output: dict[str, np.ndarray] = {}
    for column in ("value", "support", "effective"):
        wide = frame.pivot(index="window_id", columns="sample_id", values=column)
        output[column] = wide.reindex(index=windows["window_id"], columns=samples).to_numpy(float)
    return output


def contiguous_segments(evaluable: np.ndarray, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    segment = np.full(len(evaluable), -1, dtype=int)
    current = -1
    previous_end = None
    for index, ok in enumerate(evaluable):
        if not ok:
            previous_end = None
            continue
        if previous_end is None or starts[index] != previous_end:
            current += 1
        segment[index] = current
        previous_end = ends[index]
    return segment


def gap_safe_rolling_median(values: np.ndarray, segment: np.ndarray) -> np.ndarray:
    matrix = values.reshape(len(values), -1)
    output = np.full(matrix.shape, np.nan)
    half = ROLLING_WINDOWS // 2
    min_half = ROLLING_MIN_WINDOWS // 2
    for segment_id in np.unique(segment[segment >= 0]):
        rows = np.flatnonzero(segment == segment_id)
        if len(rows) < ROLLING_MIN_WINDOWS:
            continue
        block = matrix[rows]
        for position in range(len(rows)):
            reach = min(half, position, len(rows) - 1 - position)
            if reach < min_half:
                continue
            output[rows[position]] = nanmedian(block[position - reach : position + reach + 1], axis=0)
    return output.reshape(values.shape)


def group_window_summary(values: np.ndarray, support: np.ndarray, effective: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "n": np.isfinite(values).sum(axis=1),
        "mean": nanmean(values, axis=1),
        "median_cpgs": nanmedian(support, axis=1),
        "effective": np.nansum(effective, axis=1),
    }


def window_bootstrap(
    pws: np.ndarray,
    as_: np.ndarray,
    rows: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    replicates = WINDOW_BOOTSTRAP_REPLICATES
    pws_index = rng.integers(0, pws.shape[1], size=(replicates, pws.shape[1]))
    as_index = rng.integers(0, as_.shape[1], size=(replicates, as_.shape[1]))
    output = np.full((len(rows), replicates), np.nan, dtype=np.float32)
    for offset in range(0, len(rows), 400):
        chunk = rows[offset : offset + 400]
        boot_pws = nanmean(pws[chunk][:, pws_index], axis=2)
        boot_as = nanmean(as_[chunk][:, as_index], axis=2)
        output[offset : offset + len(chunk)] = boot_pws - boot_as
    return output


def reciprocal_windows(
    table: pd.DataFrame,
    windows: pd.DataFrame,
    pws_samples: list[str],
    as_samples: list[str],
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, dict[str, np.ndarray]]]:
    direct = table[
        table["mechanism"].isin([PWS_GROUP, AS_GROUP])
        & table["track_kind"].eq("combined")
        & table["parental_direction_observed"].astype(bool)
        & table["inside_common_cn1"].astype(bool)
    ]
    samples = pws_samples + as_samples
    out = windows.copy()
    starts = out["start"].to_numpy()
    ends = out["end"].to_numpy()
    in_common = out["in_common_cn1"].to_numpy(bool)
    matrices: dict[str, dict[str, np.ndarray]] = {}
    for estimator in ESTIMATORS:
        data = participant_matrices(direct, windows, samples, estimator)
        matrices[estimator] = data
        pws = {key: value[:, : len(pws_samples)] for key, value in data.items()}
        as_ = {key: value[:, len(pws_samples) :] for key, value in data.items()}
        pws_summary = group_window_summary(pws["value"], pws["support"], pws["effective"])
        as_summary = group_window_summary(as_["value"], as_["support"], as_["effective"])
        evaluable = (
            in_common
            & (pws_summary["n"] >= MIN_PWS_PARTICIPANTS)
            & (as_summary["n"] >= MIN_AS_PARTICIPANTS)
        )
        delta = np.where(evaluable, pws_summary["mean"] - as_summary["mean"], np.nan)
        segment = contiguous_segments(evaluable, starts, ends)
        smoothed = gap_safe_rolling_median(delta, segment)
        suffix = "" if estimator == PRIMARY_ESTIMATOR else f"_{estimator}"
        out[f"pws_n{suffix}"] = pws_summary["n"]
        out[f"as_n{suffix}"] = as_summary["n"]
        out[f"evaluable{suffix}"] = evaluable
        out[f"delta_beta{suffix}"] = delta
        out[f"segment_id{suffix}"] = segment
        out[f"delta_rolling_median{suffix}"] = smoothed
        if estimator != PRIMARY_ESTIMATOR:
            continue
        out["mean_beta_pws_maternal_retained"] = np.where(in_common, pws_summary["mean"], np.nan)
        out["mean_beta_as_paternal_retained"] = np.where(in_common, as_summary["mean"], np.nan)
        out["pws_median_cpgs"] = pws_summary["median_cpgs"]
        out["as_median_cpgs"] = as_summary["median_cpgs"]
        out["genomic_support_cpgs"] = nanmedian(data["support"], axis=1)
        out["effective_observations"] = pws_summary["effective"] + as_summary["effective"]
        out["pws_missing_fraction"] = np.where(in_common, 1 - pws_summary["n"] / len(pws_samples), np.nan)
        out["as_missing_fraction"] = np.where(in_common, 1 - as_summary["n"] / len(as_samples), np.nan)
        out["missing_fraction"] = np.where(
            in_common, 1 - (pws_summary["n"] + as_summary["n"]) / len(samples), np.nan
        )
        rows = np.flatnonzero(evaluable)
        boot = np.full((len(out), WINDOW_BOOTSTRAP_REPLICATES), np.nan, dtype=np.float32)
        boot[rows] = window_bootstrap(pws["value"], as_["value"], rows, rng)
        low, high = percentile_interval(boot[rows], axis=1)
        out["delta_ci_low"] = np.nan
        out["delta_ci_high"] = np.nan
        out.loc[rows, "delta_ci_low"] = low
        out.loc[rows, "delta_ci_high"] = high
        out["ci_estimable"] = out["delta_ci_low"].notna()
        boot_smoothed = gap_safe_rolling_median(boot.astype(float), segment)
        smoothed_rows = np.flatnonzero(np.isfinite(smoothed))
        out["rolling_ci_low"] = np.nan
        out["rolling_ci_high"] = np.nan
        if len(smoothed_rows):
            low, high = percentile_interval(boot_smoothed[smoothed_rows], axis=1)
            out.loc[smoothed_rows, "rolling_ci_low"] = low
            out.loc[smoothed_rows, "rolling_ci_high"] = high
    participant_long = []
    for estimator, data in matrices.items():
        frame = pd.DataFrame(
            {
                "window_id": np.repeat(windows["window_id"].to_numpy(), len(samples)),
                "start": np.repeat(starts, len(samples)),
                "sample_id": np.tile(samples, len(windows)),
                "estimator": estimator,
                "beta": data["value"].ravel(),
                "n_cpg": data["support"].ravel(),
                "effective_obs": data["effective"].ravel(),
            }
        )
        frame = frame[np.repeat(in_common, len(samples))]
        participant_long.append(frame)
    participants = pd.concat(participant_long, ignore_index=True)
    participants["mechanism"] = np.where(participants["sample_id"].isin(pws_samples), PWS_GROUP, AS_GROUP)
    participants["retained_copy"] = np.where(
        participants["mechanism"].eq(PWS_GROUP), "maternal_retained", "paternal_retained"
    )
    return out, participants, matrices


def diploid_combined_track(table: pd.DataFrame, windows: pd.DataFrame, samples: list[str]) -> pd.DataFrame:
    rows = table[table["mechanism"].isin(DIPLOID_GROUPS) & table["track_kind"].eq("combined")]
    data = participant_matrices(rows, windows, samples, PRIMARY_ESTIMATOR)
    return pd.DataFrame(
        {
            "window_id": windows["window_id"],
            "diploid_combined_mean_beta_descriptive": nanmean(data["value"], axis=1),
            "diploid_combined_n": np.isfinite(data["value"]).sum(axis=1),
        }
    )


def bootstrap_mean_difference(
    a: np.ndarray, b: np.ndarray, rng: np.random.Generator
) -> tuple[float, float]:
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    replicates = REGION_BOOTSTRAP_REPLICATES
    boot_a = a[rng.integers(0, len(a), size=(replicates, len(a)))].mean(axis=1)
    boot_b = b[rng.integers(0, len(b), size=(replicates, len(b)))].mean(axis=1)
    low, high = np.quantile(boot_a - boot_b, [0.025, 0.975])
    return float(low), float(high)


def bootstrap_mean(values: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    if len(values) < 2:
        return np.nan, np.nan
    boot = values[rng.integers(0, len(values), size=(REGION_BOOTSTRAP_REPLICATES, len(values)))].mean(axis=1)
    low, high = np.quantile(boot, [0.025, 0.975])
    return float(low), float(high)


def welch_interval(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    from scipy import stats

    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    se = np.sqrt(va + vb)
    if se == 0:
        return float(a.mean() - b.mean()), float(a.mean() - b.mean())
    df = (va + vb) ** 2 / (va**2 / (len(a) - 1) + vb**2 / (len(b) - 1))
    margin = stats.t.ppf(0.975, df) * se
    difference = a.mean() - b.mean()
    return float(difference - margin), float(difference + margin)


def analysed_interval(region: pd.Series, common: tuple[int, int]) -> tuple[int, int] | None:
    start, end = max(int(region["start"]), common[0]), min(int(region["end"]), common[1])
    return (start, end) if end > start else None


def region_participant_values(
    values: np.ndarray,
    support: np.ndarray,
    effective: np.ndarray,
    rows: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    block = values[rows]
    n_windows = np.isfinite(block).sum(axis=0)
    n_cpg = np.nansum(support[rows], axis=0)
    n_effective = np.nansum(effective[rows], axis=0)
    needed = max(1, int(np.ceil(MIN_REGION_WINDOW_FRACTION * len(rows))))
    ok = (n_windows >= needed) & (n_cpg >= MIN_REGION_CPGS)
    mean = np.where(ok, nanmean(block, axis=0), np.nan)
    return mean, n_windows, n_cpg, n_effective


def sensitivity_status(primary: float, value: float, n_pws: int, n_as: int) -> str:
    if n_pws < MIN_PWS_PARTICIPANTS or n_as < MIN_AS_PARTICIPANTS:
        return "fail_participant_loss"
    if not np.isfinite(primary) or not np.isfinite(value):
        return "fail_not_estimable"
    if abs(value - primary) <= SENSITIVITY_ABS_TOLERANCE:
        return "pass"
    if np.sign(value) == np.sign(primary) and abs(value) >= SENSITIVITY_RETENTION * abs(primary):
        return "pass"
    return "fail_effect_changed"


def combine_status(statuses: list[str]) -> str:
    failures = [status for status in statuses if status != "pass"]
    return failures[0] if failures else "pass"


def classify_region(delta: float, low: float, high: float, evaluable: bool) -> str:
    if not evaluable:
        return "not_evaluable"
    if np.isfinite(low) and low > 0 and abs(delta) >= REGION_MIN_ABS_DELTA:
        return "maternal_retained_higher"
    if np.isfinite(high) and high < 0 and abs(delta) >= REGION_MIN_ABS_DELTA:
        return "paternal_retained_higher"
    return "inconclusive"


def regional_contrasts(
    windows: pd.DataFrame,
    regions: pd.DataFrame,
    matrices: dict[str, dict[str, np.ndarray]],
    pws_samples: list[str],
    as_samples: list[str],
    common: tuple[int, int],
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    mids = windows["mid"].to_numpy()
    in_common = windows["in_common_cn1"].to_numpy(bool)
    samples = pws_samples + as_samples
    groups = np.array([PWS_GROUP] * len(pws_samples) + [AS_GROUP] * len(as_samples))
    participant_rows: list[dict[str, object]] = []
    estimator_rows: list[dict[str, object]] = []
    for _, region in regions.iterrows():
        interval = analysed_interval(region, common)
        rows = (
            np.flatnonzero((mids >= interval[0]) & (mids < interval[1]) & in_common)
            if interval
            else np.array([], dtype=int)
        )
        for estimator, data in matrices.items():
            if len(rows):
                mean, n_windows, n_cpg, n_effective = region_participant_values(
                    data["value"], data["support"], data["effective"], rows
                )
            else:
                mean = np.full(len(samples), np.nan)
                n_windows = np.zeros(len(samples), dtype=int)
                n_cpg = np.zeros(len(samples))
                n_effective = np.zeros(len(samples))
            for index, sample in enumerate(samples):
                participant_rows.append(
                    {
                        "analysis": "parent_associated_contrast",
                        "region_id": region["region_id"],
                        "sample_id": sample,
                        "cohort": groups[index],
                        "retained_copy": "maternal_retained" if groups[index] == PWS_GROUP else "paternal_retained",
                        "estimator": estimator,
                        "value": mean[index],
                        "value_type": "retained_copy_mean_beta",
                        "evaluable_windows": int(n_windows[index]),
                        "region_windows": int(len(rows)),
                        "n_cpg": float(n_cpg[index]),
                        "effective_obs": float(n_effective[index]),
                        "analysed_start": interval[0] if interval else np.nan,
                        "analysed_end": interval[1] if interval else np.nan,
                    }
                )
            a = mean[: len(pws_samples)]
            b = mean[len(pws_samples) :]
            a, b = a[np.isfinite(a)], b[np.isfinite(b)]
            evaluable = bool(interval) and len(a) >= MIN_PWS_PARTICIPANTS and len(b) >= MIN_AS_PARTICIPANTS
            delta = float(a.mean() - b.mean()) if evaluable else np.nan
            low, high = bootstrap_mean_difference(a, b, rng) if evaluable else (np.nan, np.nan)
            welch_low, welch_high = welch_interval(a, b) if evaluable else (np.nan, np.nan)
            used = np.isfinite(mean)
            estimator_rows.append(
                {
                    "region_id": region["region_id"],
                    "estimator": estimator,
                    "pws_n": int(len(a)),
                    "as_n": int(len(b)),
                    "pws_mean_beta": float(a.mean()) if len(a) else np.nan,
                    "as_mean_beta": float(b.mean()) if len(b) else np.nan,
                    "delta_beta": delta,
                    "ci_low": low,
                    "ci_high": high,
                    "ci_width": high - low if np.isfinite(low) else np.nan,
                    "welch_ci_low": welch_low,
                    "welch_ci_high": welch_high,
                    "evaluable": evaluable,
                    "region_windows": int(len(rows)),
                    "median_evaluable_windows": float(np.median(n_windows[used])) if used.any() else 0.0,
                    "total_cpg_observations": float(n_cpg[used].sum()),
                    "total_effective_obs": float(n_effective[used].sum()),
                }
            )
    by_estimator = pd.DataFrame(estimator_rows)
    summary_rows: list[dict[str, object]] = []
    for _, region in regions.iterrows():
        subset = by_estimator[by_estimator["region_id"].eq(region["region_id"])].set_index("estimator")
        primary = subset.loc[PRIMARY_ESTIMATOR]
        interval = analysed_interval(region, common)
        statuses = {
            estimator: sensitivity_status(
                float(primary["delta_beta"]),
                float(subset.loc[estimator, "delta_beta"]),
                int(subset.loc[estimator, "pws_n"]),
                int(subset.loc[estimator, "as_n"]),
            )
            for estimator in COVERAGE_ESTIMATORS + SHARED_ESTIMATORS
        }
        status = classify_region(
            float(primary["delta_beta"]),
            float(primary["ci_low"]),
            float(primary["ci_high"]),
            bool(primary["evaluable"]),
        )
        coverage = combine_status([statuses[name] for name in COVERAGE_ESTIMATORS])
        shared = combine_status([statuses[name] for name in SHARED_ESTIMATORS])
        if status == "not_evaluable":
            coverage = shared = "not_evaluable"
        region_length = int(region["end"]) - int(region["start"])
        row: dict[str, object] = {
            "display_order": int(region["display_order"]),
            "region_id": region["region_id"],
            "chrom": region["chrom"],
            "start": int(region["start"]),
            "end": int(region["end"]),
            "analysed_start": interval[0] if interval else np.nan,
            "analysed_end": interval[1] if interval else np.nan,
            "fraction_in_common_cn1": (interval[1] - interval[0]) / region_length if interval else 0.0,
            "pws_n": int(primary["pws_n"]),
            "as_n": int(primary["as_n"]),
            "pws_maternal_retained_mean_beta": primary["pws_mean_beta"],
            "as_paternal_retained_mean_beta": primary["as_mean_beta"],
            "delta_beta": primary["delta_beta"],
            "ci_low": primary["ci_low"],
            "ci_high": primary["ci_high"],
            "ci_method": "participant bootstrap, stratified by deletion group",
            "welch_ci_low": primary["welch_ci_low"],
            "welch_ci_high": primary["welch_ci_high"],
            "region_windows": primary["region_windows"],
            "median_evaluable_windows_per_participant": primary["median_evaluable_windows"],
            "total_cpg_observations": primary["total_cpg_observations"],
            "status": status,
            "coverage_sensitivity": coverage,
            "shared_cpg_sensitivity": shared,
            "robust": status.endswith("_higher") and coverage == "pass" and shared == "pass",
        }
        for estimator in COVERAGE_ESTIMATORS + SHARED_ESTIMATORS:
            row[f"delta_{estimator}"] = subset.loc[estimator, "delta_beta"]
            row[f"ci_low_{estimator}"] = subset.loc[estimator, "ci_low"]
            row[f"ci_high_{estimator}"] = subset.loc[estimator, "ci_high"]
            row[f"pws_n_{estimator}"] = int(subset.loc[estimator, "pws_n"])
            row[f"as_n_{estimator}"] = int(subset.loc[estimator, "as_n"])
            row[f"sensitivity_{estimator}"] = statuses[estimator]
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows).sort_values("display_order")
    return summary, by_estimator, pd.DataFrame(participant_rows)


def focal_intervals(
    out: pd.DataFrame,
    matrices: dict[str, dict[str, np.ndarray]],
    pws_samples: list[str],
    as_samples: list[str],
) -> pd.DataFrame:
    smooth = out["delta_rolling_median"].to_numpy()
    low = out["rolling_ci_low"].to_numpy()
    high = out["rolling_ci_high"].to_numpy()
    segment = out["segment_id"].to_numpy()
    sign = np.where(
        (smooth >= FOCAL_MIN_ABS_DELTA) & (low > 0),
        1,
        np.where((smooth <= -FOCAL_MIN_ABS_DELTA) & (high < 0), -1, 0),
    )
    candidates: list[tuple[int, int, int]] = []
    index = 0
    while index < len(sign):
        if sign[index] == 0:
            index += 1
            continue
        stop = index
        while (
            stop + 1 < len(sign)
            and sign[stop + 1] == sign[index]
            and segment[stop + 1] == segment[index]
        ):
            stop += 1
        candidates.append((index, stop, int(sign[index])))
        index = stop + 1
    rows: list[dict[str, object]] = []
    n_pws = len(pws_samples)
    for number, (first, last, direction) in enumerate(candidates, start=1):
        rows_index = np.arange(first, last + 1)
        length = len(rows_index)
        record: dict[str, object] = {
            "candidate_id": f"C{number:03d}",
            "start": int(out.loc[first, "start"]),
            "end": int(out.loc[last, "end"]),
            "n_windows": length,
            "direction": "maternal_retained_higher" if direction > 0 else "paternal_retained_higher",
            "mean_smoothed_delta": float(np.mean(smooth[rows_index])),
            "min_pws_n": int(out.loc[rows_index, "pws_n"].min()),
            "min_as_n": int(out.loc[rows_index, "as_n"].min()),
        }
        interval_values = {}
        for estimator, data in matrices.items():
            participant, _, _, _ = region_participant_values(
                data["value"], data["support"], data["effective"], rows_index
            )
            interval_values[estimator] = participant
        participant = interval_values[PRIMARY_ESTIMATOR]
        a, b = participant[:n_pws], participant[n_pws:]
        loo = []
        for drop in range(len(participant)):
            if not np.isfinite(participant[drop]):
                continue
            keep = np.isfinite(participant) & (np.arange(len(participant)) != drop)
            ka, kb = keep[:n_pws], keep[n_pws:]
            if ka.sum() < 1 or kb.sum() < 1:
                continue
            loo.append(float(a[ka].mean() - b[kb].mean()))
        record["pws_n"] = int(np.isfinite(a).sum())
        record["as_n"] = int(np.isfinite(b).sum())
        record["interval_delta"] = float(np.nanmean(a) - np.nanmean(b)) if record["pws_n"] and record["as_n"] else np.nan
        record["loo_min_abs_delta"] = float(min(abs(value) for value in loo)) if loo else np.nan
        record["loo_sign_consistent"] = bool(loo) and all(np.sign(value) == direction for value in loo)
        sensitivity_pass = True
        for estimator in COVERAGE_ESTIMATORS + SHARED_ESTIMATORS:
            values = interval_values[estimator]
            va, vb = values[:n_pws], values[n_pws:]
            n_a, n_b = int(np.isfinite(va).sum()), int(np.isfinite(vb).sum())
            mean_delta = float(np.nanmean(va) - np.nanmean(vb)) if n_a and n_b else np.nan
            record[f"delta_{estimator}"] = mean_delta
            record[f"window_retention_{estimator}"] = float(
                np.isfinite(out.loc[rows_index, f"delta_beta_{estimator}"]).mean()
            )
            passed = (
                n_a >= MIN_PWS_PARTICIPANTS
                and n_b >= MIN_AS_PARTICIPANTS
                and np.isfinite(mean_delta)
                and np.sign(mean_delta) == direction
                and abs(mean_delta) >= SENSITIVITY_RETENTION * FOCAL_MIN_ABS_DELTA
            )
            record[f"passes_{estimator}"] = bool(passed)
            sensitivity_pass &= bool(passed)
        reasons = []
        if length < FOCAL_MIN_WINDOWS:
            reasons.append(f"fewer_than_{FOCAL_MIN_WINDOWS}_windows")
        if record["pws_n"] < MIN_PWS_PARTICIPANTS or record["as_n"] < MIN_AS_PARTICIPANTS:
            reasons.append("participant_support")
        if not record["loo_sign_consistent"] or record["loo_min_abs_delta"] < SENSITIVITY_RETENTION * FOCAL_MIN_ABS_DELTA:
            reasons.append("leave_one_participant_out")
        if not sensitivity_pass:
            reasons.append("depth_or_shared_cpg_sensitivity")
        record["reproducible"] = not reasons
        record["failed_criteria"] = ";".join(reasons) if reasons else "none"
        rows.append(record)
    columns = ["candidate_id", "start", "end", "n_windows", "direction", "reproducible", "failed_criteria"]
    table = pd.DataFrame(rows, columns=None if rows else columns)
    if table.empty:
        return pd.DataFrame(columns=columns + ["focal_id"])
    table["focal_id"] = ""
    passing = table.index[table["reproducible"]]
    table.loc[passing, "focal_id"] = [f"F{i + 1}" for i in range(len(passing))]
    return table


def phase_invariant_asm(
    table: pd.DataFrame,
    windows: pd.DataFrame,
    regions: pd.DataFrame,
    cohort_samples: dict[str, list[str]],
    common: tuple[int, int],
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    diploid = table[table["mechanism"].isin(DIPLOID_GROUPS)]
    samples = [sample for group in DIPLOID_GROUPS for sample in cohort_samples[group]]
    cohorts = np.array([group for group in DIPLOID_GROUPS for _ in cohort_samples[group]])
    mids = windows["mid"].to_numpy()
    participant_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for estimator in ESTIMATORS:
        haps = {
            kind: participant_matrices(diploid[diploid["track_kind"].eq(kind)], windows, samples, estimator)
            for kind in ("hap1", "hap2")
        }
        absolute = np.abs(haps["hap1"]["value"] - haps["hap2"]["value"])
        support = np.fmin(haps["hap1"]["support"], haps["hap2"]["support"])
        support = np.where(np.isfinite(absolute), support, np.nan)
        effective = np.where(
            np.isfinite(absolute), np.fmin(haps["hap1"]["effective"], haps["hap2"]["effective"]), np.nan
        )
        for _, region in regions.iterrows():
            interval = analysed_interval(region, common)
            scope = "common_cn1_clipped" if interval else "full_region"
            interval = interval or (int(region["start"]), int(region["end"]))
            rows = np.flatnonzero((mids >= interval[0]) & (mids < interval[1]))
            if len(rows):
                mean, n_windows, n_cpg, n_effective = region_participant_values(absolute, support, effective, rows)
            else:
                mean = np.full(len(samples), np.nan)
                n_windows = np.zeros(len(samples), dtype=int)
                n_cpg = n_effective = np.zeros(len(samples))
            for index, sample in enumerate(samples):
                participant_rows.append(
                    {
                        "analysis": "phase_invariant_asm",
                        "region_id": region["region_id"],
                        "sample_id": sample,
                        "cohort": cohorts[index],
                        "retained_copy": "not_applicable_unoriented_haplotypes",
                        "estimator": estimator,
                        "value": mean[index],
                        "value_type": "mean_abs_beta_H1_minus_H2",
                        "evaluable_windows": int(n_windows[index]),
                        "region_windows": int(len(rows)),
                        "n_cpg": float(n_cpg[index]),
                        "effective_obs": float(n_effective[index]),
                        "analysed_start": interval[0],
                        "analysed_end": interval[1],
                    }
                )
            for cohort in DIPLOID_GROUPS:
                values = mean[(cohorts == cohort) & np.isfinite(mean)]
                show_interval = len(values) >= MIN_ASM_PARTICIPANTS_FOR_INTERVAL
                low, high = bootstrap_mean(values, rng) if show_interval else (np.nan, np.nan)
                summary_rows.append(
                    {
                        "display_order": int(region["display_order"]),
                        "region_id": region["region_id"],
                        "cohort": cohort,
                        "estimator": estimator,
                        "n_participants": int(len(values)),
                        "n_participants_total": int((cohorts == cohort).sum()),
                        "mean_absolute_asm": float(values.mean()) if len(values) else np.nan,
                        "median_absolute_asm": float(np.median(values)) if len(values) else np.nan,
                        "min_absolute_asm": float(values.min()) if len(values) else np.nan,
                        "max_absolute_asm": float(values.max()) if len(values) else np.nan,
                        "ci_low": low,
                        "ci_high": high,
                        "interval_label": (
                            "descriptive participant bootstrap 95% interval"
                            if show_interval
                            else f"no interval (n<{MIN_ASM_PARTICIPANTS_FOR_INTERVAL}); descriptive"
                        ),
                        "analysed_scope": scope,
                        "analysed_start": interval[0],
                        "analysed_end": interval[1],
                        "haplotype_labels": "unoriented H1/H2; no parental assignment",
                    }
                )
    return pd.DataFrame(summary_rows), pd.DataFrame(participant_rows)


def participant_missingness(
    participants: pd.DataFrame, inventory: pd.DataFrame, windows_in_common: int
) -> pd.DataFrame:
    counts = (
        participants.assign(evaluable=participants["beta"].notna())
        .groupby(["sample_id", "mechanism", "estimator"], as_index=False)["evaluable"]
        .sum()
    )
    counts["windows_in_common_cn1"] = windows_in_common
    counts["evaluable_fraction"] = counts["evaluable"] / windows_in_common
    counts["missing_fraction"] = 1 - counts["evaluable_fraction"]
    depth = inventory[inventory["track_kind"].eq("combined")][
        ["sample_id", "reference_depth_common_cn1", "downsample_fraction"]
    ]
    return counts.merge(depth, on="sample_id", how="left").rename(
        columns={"evaluable": "evaluable_windows", "reference_depth_common_cn1": "median_retained_copy_depth"}
    )


def verify_gap_safety(out: pd.DataFrame) -> tuple[bool, str]:
    problems = 0
    checked = 0
    for estimator in ESTIMATORS:
        suffix = "" if estimator == PRIMARY_ESTIMATOR else f"_{estimator}"
        smoothed = out[f"delta_rolling_median{suffix}"].to_numpy()
        evaluable = out[f"evaluable{suffix}"].to_numpy(bool)
        segment = out[f"segment_id{suffix}"].to_numpy()
        problems += int((np.isfinite(smoothed) & ~evaluable).sum())
        half = ROLLING_WINDOWS // 2
        for index in np.flatnonzero(np.isfinite(smoothed)):
            checked += 1
            rows = np.flatnonzero(segment == segment[index])
            position = int(np.searchsorted(rows, index))
            reach = min(half, position, len(rows) - 1 - position)
            used = rows[position - reach : position + reach + 1]
            if 2 * reach + 1 < ROLLING_MIN_WINDOWS or np.any(np.diff(used) != 1) or not evaluable[used].all():
                problems += 1
        for segment_id in np.unique(segment[segment >= 0]):
            rows = np.flatnonzero(segment == segment_id)
            if np.any(out["start"].to_numpy()[rows[1:]] != out["end"].to_numpy()[rows[:-1]]):
                problems += 1
    return problems == 0, f"{checked} smoothed estimates checked; {problems} gap violations"


def run() -> dict[str, Path]:
    cohort = load_cohort(METADATA_PATH)
    table = read_validated(EVIDENCE_MATRIX_PATH, MATRIX_COLUMNS, "Evidence matrix")
    inventory = read_validated(TRACK_INVENTORY_PATH, INVENTORY_COLUMNS, "Track inventory")
    core = read_validated(COMMON_CN1_PATH, CORE_COLUMNS, "Common CN=1 interval")
    regions = read_validated(REGIONS_PATH, REGION_COLUMNS | {"display_order"}, "Prespecified regions")
    missing_regions = set(REQUIRED_REGIONS) - set(regions["region_id"])
    if missing_regions:
        raise ValueError(f"Prespecified regions missing: {sorted(missing_regions)}")
    regions = regions[regions["region_id"].isin(REQUIRED_REGIONS)].sort_values("display_order").reset_index(drop=True)
    if len(core) != 1:
        raise ValueError("Common CN=1 table must contain exactly one interval")
    common = (int(core.loc[0, "start"]), int(core.loc[0, "end"]))
    matrix_groups = table.drop_duplicates("sample_id").set_index("sample_id")["mechanism"]
    for sample, mechanism in matrix_groups.items():
        if cohort.mechanism(sample) != mechanism:
            raise ValueError(f"Mechanism mismatch for {sample}: matrix {mechanism}, metadata {cohort.mechanism(sample)}")
    cohort_samples = {
        group: sorted(sample for sample in cohort.samples_for(group) if sample in matrix_groups.index)
        for group in (PWS_GROUP, AS_GROUP) + DIPLOID_GROUPS
    }
    for group, members in cohort_samples.items():
        if not members:
            raise ValueError(f"No {group} participants in evidence matrix")
    diploid_rows = table[table["mechanism"].isin(DIPLOID_GROUPS)]
    labelled = diploid_rows["parental_state"].notna() | diploid_rows["parental_direction_observed"].astype(bool)
    if labelled.any():
        raise ValueError("Parental labels found on control or DiGeorge rows")
    domain_start, domain_end = int(table["start"].min()), int(table["end"].max())
    outside = regions[(regions["start"] < domain_start) | (regions["end"] > domain_end)]
    if not outside.empty:
        raise ValueError(f"Regions extend beyond the evidence-matrix domain: {outside['region_id'].tolist()}")
    windows = (
        table[["window_id", "chrom", "start", "end", "mid"]]
        .drop_duplicates("window_id")
        .sort_values("start")
        .reset_index(drop=True)
    )
    windows["in_common_cn1"] = windows["start"].ge(common[0]) & windows["end"].le(common[1])
    flagged = set(table.loc[table["inside_common_cn1"].astype(bool), "window_id"])
    if flagged != set(windows.loc[windows["in_common_cn1"], "window_id"]):
        raise ValueError("Evidence-matrix common CN=1 flags disagree with the common interval table")
    rng = np.random.default_rng(RANDOM_SEED)
    pws_samples, as_samples = cohort_samples[PWS_GROUP], cohort_samples[AS_GROUP]
    window_table, participants, matrices = reciprocal_windows(table, windows, pws_samples, as_samples, rng)
    diploid_samples = cohort_samples["Control"] + cohort_samples["DiGeorge"]
    window_table = window_table.merge(diploid_combined_track(table, windows, diploid_samples), on="window_id")
    focal = focal_intervals(window_table, matrices, pws_samples, as_samples)
    window_table["focal_id"] = ""
    for record in focal[focal["reproducible"].astype(bool)].itertuples():
        inside = window_table["start"].ge(record.start) & window_table["end"].le(record.end)
        window_table.loc[inside, "focal_id"] = record.focal_id
    regional, regional_long, regional_participants = regional_contrasts(
        window_table, regions, matrices, pws_samples, as_samples, common, rng
    )
    asm_summary, asm_participants = phase_invariant_asm(table, windows, regions, cohort_samples, common, rng)
    missingness = participant_missingness(participants, inventory, int(windows["in_common_cn1"].sum()))
    outdir = OUTPUT_DIR
    paths = {
        "windows": outdir / "parent_associated_windows.tsv.gz",
        "window_participants": outdir / "parent_window_participant_values.tsv.gz",
        "focal": outdir / "focal_intervals.tsv",
        "regions": outdir / "regional_parent_contrasts.tsv",
        "regions_by_estimator": outdir / "regional_parent_contrasts_by_estimator.tsv",
        "participants": outdir / "regional_participant_values.tsv.gz",
        "asm": outdir / "regional_phase_invariant_asm.tsv",
        "missingness": outdir / "participant_missingness.tsv",
        "report": outdir / "figure2_analysis_report.tsv",
        "settings": outdir / "analysis_settings.json",
    }
    outside_delta = int(window_table.loc[~window_table["in_common_cn1"], "delta_beta"].notna().sum())
    if outside_delta:
        raise RuntimeError("Parental contrast estimated outside the common CN=1 interval")
    gap_ok, gap_detail = verify_gap_safety(window_table)
    if not gap_ok:
        raise RuntimeError(f"Smoothing bridged a non-evaluable interval: {gap_detail}")
    report = build_report(
        table, window_table, regional, asm_summary, focal, inventory, cohort_samples, common, gap_detail
    )
    write_table(window_table, paths["windows"])
    write_table(participants, paths["window_participants"])
    write_table(focal, paths["focal"])
    write_table(regional, paths["regions"])
    write_table(regional_long, paths["regions_by_estimator"])
    write_table(pd.concat([regional_participants, asm_participants], ignore_index=True), paths["participants"])
    write_table(asm_summary, paths["asm"])
    write_table(missingness, paths["missingness"])
    write_table(report, paths["report"])
    settings = {key: value for key, value in globals().items() if key.isupper() and isinstance(value, (int, float, str, tuple, dict))}
    paths["settings"].write_text(json.dumps(settings, indent=2, default=str) + "\n")
    return paths


def build_report(
    table: pd.DataFrame,
    windows: pd.DataFrame,
    regional: pd.DataFrame,
    asm: pd.DataFrame,
    focal: pd.DataFrame,
    inventory: pd.DataFrame,
    cohort_samples: dict[str, list[str]],
    common: tuple[int, int],
    gap_detail: str,
) -> pd.DataFrame:
    rows: list[tuple[str, str, str, str]] = []

    def add(section: str, item: str, value: object, status: str = "info") -> None:
        rows.append((section, item, str(value), status))

    for label, path in (
        ("metadata", METADATA_PATH),
        ("evidence_matrix", EVIDENCE_MATRIX_PATH),
        ("track_inventory", TRACK_INVENTORY_PATH),
        ("common_cn1_interval", COMMON_CN1_PATH),
        ("prespecified_regions", REGIONS_PATH),
    ):
        add("input_validation", label, path, "pass")
    add("input_validation", "required_matrix_columns", ",".join(sorted(MATRIX_COLUMNS)), "pass")
    for group, members in cohort_samples.items():
        add("cohort", f"{group}_participants", f"{len(members)}: {','.join(members)}")
    add("parental_contrast", "common_reciprocal_cn1_interval", f"chr15:{common[0]}-{common[1]}", "pass")
    add("parental_contrast", "windows_in_common_cn1", int(windows["in_common_cn1"].sum()))
    add("parental_contrast", "evaluable_windows", int(windows["evaluable"].sum()))
    add(
        "parental_contrast",
        "estimates_outside_common_cn1",
        int(windows.loc[~windows["in_common_cn1"], "delta_beta"].notna().sum()),
        "pass",
    )
    add(
        "parental_contrast",
        "bootstrap_unit",
        f"participants resampled within PWS-DEL and AS-DEL; windows and CpGs never resampled; "
        f"{WINDOW_BOOTSTRAP_REPLICATES} window-level and {REGION_BOOTSTRAP_REPLICATES} regional replicates",
        "pass",
    )
    add("parental_contrast", "group_weighting", "equal weight per participant; no depth weighting of beta", "pass")
    add("smoothing", "gap_safe_rolling_median", f"symmetric centred median of {ROLLING_MIN_WINDOWS}-{ROLLING_WINDOWS} contiguous evaluable 1-kb windows; never crosses a non-evaluable window; {gap_detail}", "pass")
    diploid = table[table["mechanism"].isin(DIPLOID_GROUPS)]
    add(
        "phase_invariant_asm",
        "parental_labels_on_diploid_haplotypes",
        int(diploid["parental_state"].notna().sum()),
        "pass",
    )
    add("phase_invariant_asm", "statistic", "regional mean of |beta_H1 - beta_H2| per participant; H1/H2 unoriented")
    for set_name, frame in inventory.dropna(subset=["track_set"]).groupby("track_set"):
        add(
            "depth",
            f"{set_name}_common_depth_target",
            f"{frame['downsample_target_depth'].iloc[0]:.1f}x; fractions {frame['downsample_fraction'].min():.2f}-{frame['downsample_fraction'].max():.2f}",
        )
    reproducible = focal[focal["reproducible"].astype(bool)] if not focal.empty else focal
    add("focal_intervals", "candidates", len(focal))
    add("focal_intervals", "reproducible", len(reproducible))
    for record in reproducible.itertuples():
        add(
            "focal_intervals",
            record.focal_id,
            f"chr15:{record.start}-{record.end} ({record.n_windows} kb) {record.direction}; "
            f"mean smoothed delta {record.mean_smoothed_delta:+.3f}",
        )
    for record in regional.itertuples():
        add(
            "regional_contrast",
            record.region_id,
            f"delta {record.delta_beta:+.3f} [{record.ci_low:+.3f}, {record.ci_high:+.3f}]; "
            f"PWS n={record.pws_n}, AS n={record.as_n}; {record.status}; "
            f"coverage {record.coverage_sensitivity}; shared-CpG {record.shared_cpg_sensitivity}",
            "robust" if record.robust else ("inconclusive" if record.status == "inconclusive" else "not_robust"),
        )
    primary_asm = asm[asm["estimator"].eq(PRIMARY_ESTIMATOR)]
    for record in primary_asm.itertuples():
        add(
            "phase_invariant_asm",
            f"{record.region_id} | {record.cohort}",
            f"mean |H1-H2| {record.mean_absolute_asm:.3f} (n={record.n_participants}); {record.interval_label}",
            "descriptive",
        )
    add(
        "limitations",
        "small_groups",
        f"AS-DEL n={len(cohort_samples[AS_GROUP])} and Control n={len(cohort_samples['Control'])}; "
        "percentile bootstrap intervals are descriptive and can be narrow; Welch intervals are exported for comparison",
    )
    add("limitations", "depth_correction", "depth sensitivity analyses test robustness only and do not create biological evidence")
    add("limitations", "digeorge", "DiGeorge participants are an independent diploid disease-control cohort, not parental references")
    return pd.DataFrame(rows, columns=["section", "item", "value", "status"])


def main() -> None:
    paths = run()
    for name, path in paths.items():
        print(f"{name}\t{path}")


if __name__ == "__main__":
    main()
