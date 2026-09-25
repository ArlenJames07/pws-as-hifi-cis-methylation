from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from .methylation import MethylationTrack


def summarize_track_windows(
    track: MethylationTrack,
    windows: pd.DataFrame,
    depth_caps: Sequence[float] = (5.0, 10.0, 15.0, 20.0),
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    positions = track.position
    for window in windows.itertuples(index=False):
        left = int(np.searchsorted(positions, int(window.start), side="left"))
        right = int(np.searchsorted(positions, int(window.end), side="left"))
        beta = track.beta[left:right]
        coverage = track.coverage[left:right]
        row: dict[str, object] = {
            "chrom": window.chrom,
            "window_id": window.window_id,
            "start": int(window.start),
            "end": int(window.end),
            "mid": float(window.mid),
            "n_cpg": int(len(beta)),
            "median_depth": float(np.median(coverage)) if len(beta) else np.nan,
            "mean_depth": float(np.mean(coverage)) if len(beta) else np.nan,
            "total_depth": float(np.sum(coverage)) if len(beta) else 0.0,
            "beta_site_mean": float(np.mean(beta)) if len(beta) else np.nan,
            "beta_read_weighted": (
                float(np.average(beta, weights=coverage)) if len(beta) and coverage.sum() > 0 else np.nan
            ),
        }
        for cap in depth_caps:
            weights = np.minimum(coverage, float(cap))
            label = str(float(cap)).rstrip("0").rstrip(".").replace(".", "p")
            row[f"beta_depth_cap_{label}"] = (
                float(np.average(beta, weights=weights)) if len(beta) and weights.sum() > 0 else np.nan
            )
        rows.append(row)
    return pd.DataFrame(rows)


def bootstrap_group_contrast(
    table: pd.DataFrame,
    value_col: str,
    sample_col: str,
    group_col: str,
    group_a: str,
    group_b: str,
    replicates: int,
    seed: int,
) -> pd.DataFrame:
    required = {"window_id", value_col, sample_col, group_col}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Contrast input lacks columns: {sorted(missing)}")
    rng = np.random.default_rng(seed)
    output: list[dict[str, object]] = []
    for window_id, frame in table.groupby("window_id", sort=False):
        vectors: dict[str, np.ndarray] = {}
        for group in (group_a, group_b):
            values = (
                frame.loc[frame[group_col] == group]
                .groupby(sample_col, sort=False)[value_col]
                .mean()
                .dropna()
                .to_numpy(dtype=float)
            )
            vectors[group] = values
        a = vectors[group_a]
        b = vectors[group_b]
        base = frame.iloc[0]
        row: dict[str, object] = {
            "window_id": window_id,
            "start": int(base["start"]),
            "end": int(base["end"]),
            "mid": float(base["mid"]),
            "group_a": group_a,
            "group_b": group_b,
            "n_a": int(len(a)),
            "n_b": int(len(b)),
            "mean_a": float(np.mean(a)) if len(a) else np.nan,
            "mean_b": float(np.mean(b)) if len(b) else np.nan,
            "delta": float(np.mean(a) - np.mean(b)) if len(a) and len(b) else np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
        }
        if len(a) and len(b) and replicates > 0:
            boot_a = a[rng.integers(0, len(a), size=(replicates, len(a)))].mean(axis=1)
            boot_b = b[rng.integers(0, len(b), size=(replicates, len(b)))].mean(axis=1)
            row["ci_low"], row["ci_high"] = np.quantile(boot_a - boot_b, [0.025, 0.975])
        output.append(row)
    return pd.DataFrame(output)
