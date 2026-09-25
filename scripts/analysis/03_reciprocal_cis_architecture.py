#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from cis_analysis import bootstrap_group_contrast


PROJECT_ROOT = Path(__file__).resolve().parents[2]
EVIDENCE_MATRIX_PATH = (
    PROJECT_ROOT
    / "results"
    / "analysis"
    / "01_evidence_matrix"
    / "chr15_window_evidence_matrix.tsv.gz"
)
REGIONS_PATH = PROJECT_ROOT / "results" / "analysis" / "prespecified_regions.tsv"
OUTPUT_DIR = PROJECT_ROOT / "results" / "analysis" / "03_cis_architecture"
MIN_CPGS = 3
MIN_PWS_PARTICIPANTS = 3
MIN_AS_PARTICIPANTS = 2
BOOTSTRAP_REPLICATES = 2_000
RANDOM_SEED = 20260924
FOCAL_EFFECT_THRESHOLD = 0.10
EQUIVALENCE_MARGIN = 0.10


def write_table(table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if path.suffix == ".gz" else None
    table.to_csv(path, sep="\t", index=False, na_rep="NA", compression=compression)


def participant_scaffold(table: pd.DataFrame, min_cpgs: int) -> pd.DataFrame:
    data = table[
        table["mechanism"].isin(["Control", "DiGeorge"])
        & table["track_kind"].eq("combined")
        & table["n_cpg"].ge(min_cpgs)
    ].copy()
    participant = (
        data.groupby(
            ["window_id", "start", "end", "mid", "sample_id", "mechanism"], as_index=False
        )["beta_site_mean"]
        .mean()
        .dropna(subset=["beta_site_mean"])
    )
    return (
        participant.groupby(["window_id", "start", "end", "mid"], as_index=False)
        .agg(
            scaffold_beta=("beta_site_mean", "mean"),
            scaffold_sd=("beta_site_mean", "std"),
            scaffold_n=("sample_id", "nunique"),
            scaffold_control_n=("mechanism", lambda x: int((x == "Control").sum())),
            scaffold_digeorge_n=("mechanism", lambda x: int((x == "DiGeorge").sum())),
        )
    )


def phase_invariant_asm(table: pd.DataFrame, min_cpgs: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    data = table[
        table["mechanism"].isin(["Control", "DiGeorge"])
        & table["track_kind"].isin(["hap1", "hap2"])
        & table["n_cpg"].ge(min_cpgs)
    ].copy()
    wide = data.pivot_table(
        index=["window_id", "start", "end", "mid", "sample_id", "mechanism"],
        columns="track_kind",
        values="beta_site_mean",
        aggfunc="mean",
    ).reset_index()
    if "hap1" not in wide or "hap2" not in wide:
        raise ValueError("Both hap1 and hap2 tracks are required for phase-invariant ASM")
    wide["absolute_asm"] = (wide["hap1"] - wide["hap2"]).abs()
    wide["signed_hp1_minus_hp2"] = wide["hap1"] - wide["hap2"]
    summary = (
        wide.dropna(subset=["absolute_asm"])
        .groupby(["window_id", "start", "end", "mid"], as_index=False)
        .agg(
            absolute_asm_mean=("absolute_asm", "mean"),
            absolute_asm_sd=("absolute_asm", "std"),
            absolute_asm_n=("sample_id", "nunique"),
            signed_hp_difference_mean=("signed_hp1_minus_hp2", "mean"),
        )
    )
    return wide, summary


def bootstrap_scaffold(table: pd.DataFrame, replicates: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    for window_id, frame in table.groupby("window_id", sort=False):
        values = frame.groupby("sample_id")["beta_site_mean"].mean().dropna().to_numpy(float)
        base = frame.iloc[0]
        low = high = np.nan
        if len(values) and replicates > 0:
            samples = values[rng.integers(0, len(values), size=(replicates, len(values)))]
            low, high = np.quantile(samples.mean(axis=1), [0.025, 0.975])
        rows.append(
            {
                "window_id": window_id,
                "scaffold_ci_low": low,
                "scaffold_ci_high": high,
                "start": int(base["start"]),
                "end": int(base["end"]),
                "mid": float(base["mid"]),
            }
        )
    return pd.DataFrame(rows)


def classify_delta(table: pd.DataFrame, focal_threshold: float, equivalence_margin: float) -> pd.Series:
    conditions = [
        table["ci_low"].gt(focal_threshold),
        table["ci_high"].lt(-focal_threshold),
        table["ci_low"].gt(-equivalence_margin) & table["ci_high"].lt(equivalence_margin),
    ]
    labels = ["maternal_higher", "paternal_higher", "shared_equivalent"]
    return pd.Series(np.select(conditions, labels, default="uncertain"), index=table.index)


def region_summary(windows: pd.DataFrame, regions: pd.DataFrame) -> pd.DataFrame:
    required = {"region_id", "chrom", "start", "end"}
    missing = required - set(regions.columns)
    if missing:
        raise ValueError(f"Region table lacks columns: {sorted(missing)}")
    rows: list[dict[str, object]] = []
    for region in regions.itertuples(index=False):
        overlap = windows[(windows["mid"] >= int(region.start)) & (windows["mid"] < int(region.end))]
        weights = (overlap["end"] - overlap["start"]).to_numpy(float)
        valid = overlap["delta_maternal_minus_paternal"].notna().to_numpy()
        delta = (
            float(
                np.average(
                    overlap.loc[valid, "delta_maternal_minus_paternal"],
                    weights=weights[valid],
                )
            )
            if valid.any()
            else np.nan
        )
        rows.append(
            {
                "region_id": region.region_id,
                "chrom": region.chrom,
                "start": int(region.start),
                "end": int(region.end),
                "n_windows": int(len(overlap)),
                "delta_maternal_minus_paternal": delta,
                "median_scaffold_beta": float(overlap["scaffold_beta"].median()) if len(overlap) else np.nan,
                "median_absolute_asm": float(overlap["absolute_asm_mean"].median()) if len(overlap) else np.nan,
                "fraction_maternal_higher": float((overlap["state"] == "maternal_higher").mean()) if len(overlap) else np.nan,
                "fraction_paternal_higher": float((overlap["state"] == "paternal_higher").mean()) if len(overlap) else np.nan,
                "fraction_shared_equivalent": float((overlap["state"] == "shared_equivalent").mean()) if len(overlap) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def run() -> dict[str, Path]:
    if not EVIDENCE_MATRIX_PATH.is_file():
        raise FileNotFoundError(f"Evidence matrix not found: {EVIDENCE_MATRIX_PATH}")
    if not REGIONS_PATH.is_file():
        raise FileNotFoundError(f"Prespecified regions not found: {REGIONS_PATH}")
    table = pd.read_csv(EVIDENCE_MATRIX_PATH, sep="\t", low_memory=False)
    direct = table[
        table["parental_direction_observed"].astype(bool)
        & table["inside_common_cn1"].astype(bool)
        & table["track_kind"].eq("combined")
        & table["n_cpg"].ge(MIN_CPGS)
    ].copy()
    reciprocal = bootstrap_group_contrast(
        direct,
        "beta_site_mean",
        "sample_id",
        "mechanism",
        "PWS-DEL",
        "AS-DEL",
        BOOTSTRAP_REPLICATES,
        RANDOM_SEED,
    )
    reciprocal["supported"] = reciprocal["n_a"].ge(
        MIN_PWS_PARTICIPANTS
    ) & reciprocal["n_b"].ge(MIN_AS_PARTICIPANTS)
    reciprocal.loc[~reciprocal["supported"], ["delta", "ci_low", "ci_high"]] = np.nan
    reciprocal["state"] = classify_delta(
        reciprocal, FOCAL_EFFECT_THRESHOLD, EQUIVALENCE_MARGIN
    )
    reciprocal.loc[~reciprocal["supported"], "state"] = "insufficient_support"
    scaffold = participant_scaffold(table, MIN_CPGS)
    scaffold_input = table[
        table["mechanism"].isin(["Control", "DiGeorge"])
        & table["track_kind"].eq("combined")
        & table["n_cpg"].ge(MIN_CPGS)
    ]
    scaffold_ci = bootstrap_scaffold(
        scaffold_input,
        BOOTSTRAP_REPLICATES,
        RANDOM_SEED + 1,
    )
    scaffold = scaffold.merge(scaffold_ci[["window_id", "scaffold_ci_low", "scaffold_ci_high"]], on="window_id", how="left")
    asm_participant, asm_summary = phase_invariant_asm(table, MIN_CPGS)
    architecture = reciprocal.merge(scaffold, on=["window_id", "start", "end", "mid"], how="outer")
    architecture = architecture.merge(
        asm_summary, on=["window_id", "start", "end", "mid"], how="outer"
    ).sort_values(["start", "end"])
    architecture = architecture.rename(
        columns={
            "mean_a": "maternal_retained_mean",
            "mean_b": "paternal_retained_mean",
            "n_a": "pws_n",
            "n_b": "as_n",
            "delta": "delta_maternal_minus_paternal",
            "ci_low": "delta_ci_low",
            "ci_high": "delta_ci_high",
        }
    )
    outdir = OUTPUT_DIR
    paths = {
        "windows": outdir / "cis_architecture_windows.tsv.gz",
        "asm": outdir / "phase_invariant_asm_by_participant.tsv.gz",
        "regions": outdir / "prespecified_region_effects.tsv",
    }
    write_table(architecture, paths["windows"])
    write_table(asm_participant, paths["asm"])
    regions = pd.read_csv(REGIONS_PATH, sep="\t")
    region_table = region_summary(architecture, regions)
    write_table(region_table, paths["regions"])
    return paths


def main() -> None:
    paths = run()
    for name, path in paths.items():
        print(f"{name}\t{path}")


if __name__ == "__main__":
    main()
