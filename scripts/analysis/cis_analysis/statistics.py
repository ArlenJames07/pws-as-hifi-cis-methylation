from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from .methylation import MethylationTrack


def _window_index(positions: np.ndarray, windows: pd.DataFrame) -> np.ndarray:
    starts = windows["start"].to_numpy(dtype=np.int64)
    ends = windows["end"].to_numpy(dtype=np.int64)
    index = np.searchsorted(starts, positions, side="right") - 1
    valid = (index >= 0) & (positions < ends[np.clip(index, 0, len(ends) - 1)])
    return np.where(valid, index, -1)


def _cap_label(cap: float) -> str:
    return str(float(cap)).rstrip("0").rstrip(".").replace(".", "p")


def _site_summary(
    track: MethylationTrack,
    windows: pd.DataFrame,
    effective_cap: float,
    mask: np.ndarray | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    beta, coverage, positions = track.beta, track.coverage, track.position
    if mask is not None:
        beta, coverage, positions = beta[mask], coverage[mask], positions[mask]
    index = _window_index(positions, windows)
    keep = index >= 0
    frame = pd.DataFrame(
        {
            "window": index[keep],
            "beta": beta[keep],
            "coverage": coverage[keep],
            "capped": np.minimum(coverage[keep], float(effective_cap)),
        }
    )
    grouped = frame.groupby("window", sort=True)
    out = pd.DataFrame(index=pd.RangeIndex(len(windows)))
    out["n_cpg"] = grouped.size().reindex(out.index, fill_value=0).astype(int)
    out["beta_site_mean"] = grouped["beta"].mean().reindex(out.index)
    out["effective_obs"] = grouped["capped"].sum().reindex(out.index, fill_value=0.0)
    out["median_depth"] = grouped["coverage"].median().reindex(out.index)
    out["mean_depth"] = grouped["coverage"].mean().reindex(out.index)
    out["total_depth"] = grouped["coverage"].sum().reindex(out.index, fill_value=0.0)
    frame["weighted"] = frame["beta"] * frame["coverage"]
    out["beta_read_weighted"] = (
        frame.groupby("window")["weighted"].sum().reindex(out.index) / out["total_depth"].replace(0, np.nan)
    )
    return out, frame


def summarize_track_windows(
    track: MethylationTrack,
    windows: pd.DataFrame,
    depth_caps: Sequence[float] = (5.0, 10.0, 15.0, 20.0),
    effective_cap: float = 15.0,
    downsampled: MethylationTrack | None = None,
    shared_mask: np.ndarray | None = None,
) -> pd.DataFrame:
    base, frame = _site_summary(track, windows, effective_cap)
    out = windows[["chrom", "window_id", "start", "end", "mid"]].reset_index(drop=True).copy()
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)
    out["mid"] = out["mid"].astype(float)
    for column in (
        "n_cpg",
        "median_depth",
        "mean_depth",
        "total_depth",
        "beta_site_mean",
        "beta_read_weighted",
        "effective_obs",
    ):
        out[column] = base[column].to_numpy()
    for cap in depth_caps:
        weights = np.minimum(frame["coverage"].to_numpy(), float(cap))
        numerator = pd.Series(frame["beta"].to_numpy() * weights).groupby(frame["window"].to_numpy()).sum()
        denominator = pd.Series(weights).groupby(frame["window"].to_numpy()).sum()
        out[f"beta_depth_cap_{_cap_label(cap)}"] = (
            (numerator / denominator.replace(0, np.nan)).reindex(out.index).to_numpy()
        )
    if downsampled is not None:
        thinned, _ = _site_summary(downsampled, windows, effective_cap)
        out["n_cpg_downsampled"] = thinned["n_cpg"].to_numpy()
        out["beta_downsampled"] = thinned["beta_site_mean"].to_numpy()
        out["effective_obs_downsampled"] = thinned["effective_obs"].to_numpy()
    if shared_mask is not None:
        shared, _ = _site_summary(track, windows, effective_cap, shared_mask)
        out["n_cpg_shared"] = shared["n_cpg"].to_numpy()
        out["beta_shared_cpg"] = shared["beta_site_mean"].to_numpy()
        out["effective_obs_shared"] = shared["effective_obs"].to_numpy()
    return out


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
