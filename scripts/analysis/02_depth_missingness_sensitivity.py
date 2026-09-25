#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
EVIDENCE_MATRIX_PATH = (
    PROJECT_ROOT
    / "results"
    / "analysis"
    / "01_evidence_matrix"
    / "chr15_window_evidence_matrix.tsv.gz"
)
OUTPUT_DIR = PROJECT_ROOT / "results" / "analysis" / "02_depth_sensitivity"
MIN_CPG_GRID = (3, 5, 10)
PRIMARY_MIN_CPGS = 3
MIN_PWS_PARTICIPANTS = 3
MIN_AS_PARTICIPANTS = 2


def write_table(table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=False, na_rep="NA")


def safe_correlation(x: pd.Series, y: pd.Series, method: str) -> float:
    valid = x.notna() & y.notna()
    if valid.sum() < 3 or x[valid].nunique() < 2 or y[valid].nunique() < 2:
        return np.nan
    return float(x[valid].corr(y[valid], method=method))


def contrast_windows(
    table: pd.DataFrame,
    estimator: str,
    min_cpgs: int,
    min_pws: int,
    min_as: int,
) -> pd.DataFrame:
    direct = table[
        table["parental_direction_observed"].astype(bool)
        & table["inside_common_cn1"].astype(bool)
        & table["track_kind"].eq("combined")
        & table["n_cpg"].ge(min_cpgs)
    ].copy()
    participant = (
        direct.groupby(
            ["window_id", "start", "end", "mid", "sample_id", "mechanism"], as_index=False
        )[estimator]
        .mean()
        .dropna(subset=[estimator])
    )
    summaries = (
        participant.groupby(["window_id", "start", "end", "mid", "mechanism"])[estimator]
        .agg(["mean", "count"])
        .reset_index()
    )
    means = summaries.pivot(index=["window_id", "start", "end", "mid"], columns="mechanism", values="mean")
    counts = summaries.pivot(index=["window_id", "start", "end", "mid"], columns="mechanism", values="count")
    out = means.join(counts, lsuffix="_mean", rsuffix="_n").reset_index()
    for column in ("PWS-DEL_mean", "AS-DEL_mean", "PWS-DEL_n", "AS-DEL_n"):
        if column not in out:
            out[column] = np.nan
    out["delta"] = out["PWS-DEL_mean"] - out["AS-DEL_mean"]
    out["supported"] = out["PWS-DEL_n"].ge(min_pws) & out["AS-DEL_n"].ge(min_as)
    out["estimator"] = estimator
    out["min_cpgs"] = min_cpgs
    return out


def run() -> dict[str, Path]:
    if not EVIDENCE_MATRIX_PATH.is_file():
        raise FileNotFoundError(f"Evidence matrix not found: {EVIDENCE_MATRIX_PATH}")
    table = pd.read_csv(EVIDENCE_MATRIX_PATH, sep="\t", low_memory=False)
    required = {
        "sample_id",
        "mechanism",
        "track_kind",
        "window_id",
        "n_cpg",
        "median_depth",
        "beta_site_mean",
        "beta_read_weighted",
        "parental_direction_observed",
        "inside_common_cn1",
    }
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Evidence matrix lacks columns: {sorted(missing)}")
    estimators = [
        column
        for column in table.columns
        if column in {"beta_site_mean", "beta_read_weighted"}
        or column.startswith("beta_depth_cap_")
    ]
    sample_qc = (
        table.groupby(["sample_id", "mechanism", "track_kind"], as_index=False)
        .agg(
            n_windows=("window_id", "nunique"),
            windows_with_3_cpg=("n_cpg", lambda x: int((x >= 3).sum())),
            windows_with_5_cpg=("n_cpg", lambda x: int((x >= 5).sum())),
            windows_with_10_cpg=("n_cpg", lambda x: int((x >= 10).sum())),
            median_cpgs=("n_cpg", "median"),
            median_depth=("median_depth", "median"),
            beta_site_mean=("beta_site_mean", "mean"),
            beta_read_weighted=("beta_read_weighted", "mean"),
        )
    )
    for threshold in (3, 5, 10):
        sample_qc[f"recovery_fraction_{threshold}_cpg"] = (
            sample_qc[f"windows_with_{threshold}_cpg"] / sample_qc["n_windows"]
        )
    sample_qc["absolute_estimator_shift"] = (
        sample_qc["beta_site_mean"] - sample_qc["beta_read_weighted"]
    ).abs()
    group_qc = (
        sample_qc.groupby(["mechanism", "track_kind"], as_index=False)
        .agg(
            n_participants=("sample_id", "nunique"),
            median_depth=("median_depth", "median"),
            depth_min=("median_depth", "min"),
            depth_max=("median_depth", "max"),
            median_recovery_3_cpg=("recovery_fraction_3_cpg", "median"),
            median_recovery_5_cpg=("recovery_fraction_5_cpg", "median"),
            median_recovery_10_cpg=("recovery_fraction_10_cpg", "median"),
            median_estimator_shift=("absolute_estimator_shift", "median"),
        )
    )
    relationships: list[dict[str, object]] = []
    for kind, frame in sample_qc.groupby("track_kind"):
        for recovery in (
            "recovery_fraction_3_cpg",
            "recovery_fraction_5_cpg",
            "recovery_fraction_10_cpg",
            "absolute_estimator_shift",
        ):
            relationships.append(
                {
                    "track_kind": kind,
                    "outcome": recovery,
                    "n": int(frame[["median_depth", recovery]].dropna().shape[0]),
                    "pearson_r": safe_correlation(frame["median_depth"], frame[recovery], "pearson"),
                    "spearman_rho": safe_correlation(frame["median_depth"], frame[recovery], "spearman"),
                }
            )
    contrast_tables: list[pd.DataFrame] = []
    for estimator in estimators:
        for threshold in MIN_CPG_GRID:
            contrast_tables.append(
                contrast_windows(
                    table,
                    estimator,
                    threshold,
                    MIN_PWS_PARTICIPANTS,
                    MIN_AS_PARTICIPANTS,
                )
            )
    contrast = pd.concat(contrast_tables, ignore_index=True)
    primary = contrast[
        contrast["estimator"].eq("beta_site_mean")
        & contrast["min_cpgs"].eq(PRIMARY_MIN_CPGS)
        & contrast["supported"]
    ][["window_id", "delta"]].rename(columns={"delta": "primary_delta"})
    concordance_rows: list[dict[str, object]] = []
    for (estimator, threshold), frame in contrast.groupby(["estimator", "min_cpgs"]):
        joined = frame[frame["supported"]].merge(primary, on="window_id", how="inner")
        difference = joined["delta"] - joined["primary_delta"]
        concordance_rows.append(
            {
                "estimator": estimator,
                "min_cpgs": int(threshold),
                "n_common_windows": len(joined),
                "pearson_r": safe_correlation(joined["delta"], joined["primary_delta"], "pearson"),
                "spearman_rho": safe_correlation(joined["delta"], joined["primary_delta"], "spearman"),
                "rmse": float(np.sqrt(np.mean(np.square(difference)))) if len(joined) else np.nan,
                "median_absolute_difference": float(np.median(np.abs(difference))) if len(joined) else np.nan,
                "sign_concordance": (
                    float(np.mean(np.sign(joined["delta"]) == np.sign(joined["primary_delta"])))
                    if len(joined)
                    else np.nan
                ),
            }
        )
    outdir = OUTPUT_DIR
    paths = {
        "sample_qc": outdir / "depth_missingness_by_sample.tsv",
        "group_qc": outdir / "depth_missingness_by_group.tsv",
        "relationships": outdir / "depth_recovery_relationships.tsv",
        "contrast": outdir / "reciprocal_contrast_sensitivity.tsv",
        "concordance": outdir / "estimator_concordance.tsv",
    }
    write_table(sample_qc, paths["sample_qc"])
    write_table(group_qc, paths["group_qc"])
    write_table(pd.DataFrame(relationships), paths["relationships"])
    write_table(contrast, paths["contrast"])
    write_table(pd.DataFrame(concordance_rows), paths["concordance"])
    return paths


def main() -> None:
    paths = run()
    for name, path in paths.items():
        print(f"{name}\t{path}")


if __name__ == "__main__":
    main()
