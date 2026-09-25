#!/usr/bin/env python3
from __future__ import annotations

import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from cis_analysis import (
    MethylationTrack,
    classify_evidence,
    downsample_track,
    find_track,
    fixed_windows,
    load_cohort,
    load_deletion_map,
    read_track,
    read_track_metadata,
    summarize_track_windows,
)


# ============================== CONFIGURATION ==============================
PROJECT_ROOT = Path(__file__).resolve().parents[2]
METADATA_PATH = PROJECT_ROOT / "assets" / "metadata.csv"
METHYLATION_DIR = PROJECT_ROOT / "results" / "06_methylation"
DELETION_PROVENANCE_PATH = (
    PROJECT_ROOT
    / "results"
    / "07_figures"
    / "figure_1"
    / "supplementary"
    / "Supplementary_chr15_deletion_classification.tsv"
)
OUTPUT_DIR = PROJECT_ROOT / "results" / "analysis" / "01_evidence_matrix"
CHROM = "chr15"
DOMAIN_START = 20_500_000
DOMAIN_END = 28_000_000
WINDOW_SIZE = 1_000
CN1_BREAKPOINT_BUFFER = 75_000
MIN_CPGS = 3
DEPTH_CAPS = (5.0, 10.0, 15.0, 20.0)
EFFECTIVE_COVERAGE_CAP = 15.0
DOWNSAMPLE_TARGET_DEPTH: float | None = None
DOWNSAMPLE_MIN_COVERAGE = 4
DOWNSAMPLE_SEED = 20260925
TRACK_SETS = {
    "reciprocal_deletion": {
        "mechanisms": ("PWS-DEL", "AS-DEL"),
        "kinds": ("combined",),
        "shared_scope": "all_tracks",
    },
    "diploid_haplotypes": {
        "mechanisms": ("Control", "DiGeorge"),
        "kinds": ("hap1", "hap2"),
        "shared_scope": "within_participant",
    },
    "diploid_combined": {
        "mechanisms": ("Control", "DiGeorge"),
        "kinds": ("combined",),
        "shared_scope": "all_tracks",
    },
}
# ===========================================================================


def write_table(table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if path.suffix == ".gz" else None
    table.to_csv(path, sep="\t", index=False, na_rep="NA", compression=compression)


def track_set_for(mechanism: str, kind: str) -> str | None:
    for name, spec in TRACK_SETS.items():
        if mechanism in spec["mechanisms"] and kind in spec["kinds"]:
            return name
    return None


def reference_depth(track: MethylationTrack, interval: tuple[int, int]) -> float:
    inside = (track.position >= interval[0]) & (track.position < interval[1])
    values = track.coverage[inside] if inside.any() else track.coverage
    return float(np.median(values)) if len(values) else np.nan


def shared_positions(
    tracks: dict[tuple[str, str], MethylationTrack],
    members: list[tuple[str, str]],
    scope: str,
) -> dict[tuple[str, str], np.ndarray]:
    if scope == "all_tracks":
        common = None
        for key in members:
            positions = tracks[key].position
            common = positions if common is None else np.intersect1d(common, positions, assume_unique=True)
        return {key: common if common is not None else np.array([], dtype=np.int64) for key in members}
    if scope == "within_participant":
        by_sample: dict[str, list[tuple[str, str]]] = {}
        for key in members:
            by_sample.setdefault(key[0], []).append(key)
        output: dict[tuple[str, str], np.ndarray] = {}
        for keys in by_sample.values():
            common = None
            for key in keys:
                positions = tracks[key].position
                common = positions if common is None else np.intersect1d(common, positions, assume_unique=True)
            for key in keys:
                output[key] = common if len(keys) > 1 else np.array([], dtype=np.int64)
        return output
    raise ValueError(f"Unknown shared-CpG scope: {scope}")


def build_matrix() -> dict[str, Path]:
    if not METADATA_PATH.is_file():
        raise FileNotFoundError(f"Metadata not found: {METADATA_PATH}")
    if not METHYLATION_DIR.is_dir():
        raise FileNotFoundError(f"Methylation directory not found: {METHYLATION_DIR}")
    if not DELETION_PROVENANCE_PATH.is_file():
        raise FileNotFoundError(
            f"Figure 1 deletion provenance not found: {DELETION_PROVENANCE_PATH}"
        )
    cohort = load_cohort(METADATA_PATH)
    deletion_samples = cohort.samples_for("PWS-DEL") + cohort.samples_for("AS-DEL")
    deletion_map = load_deletion_map(
        DELETION_PROVENANCE_PATH,
        deletion_samples,
        CN1_BREAKPOINT_BUFFER,
    )
    windows = fixed_windows(CHROM, DOMAIN_START, DOMAIN_END, WINDOW_SIZE)
    common_interval = (deletion_map.common_start, deletion_map.common_end)
    if common_interval[0] < DOMAIN_START or common_interval[1] > DOMAIN_END:
        raise ValueError(
            f"Common CN=1 interval {common_interval} extends beyond the analysed domain "
            f"{DOMAIN_START}-{DOMAIN_END}; widen DOMAIN_START/DOMAIN_END"
        )
    tracks: dict[tuple[str, str], MethylationTrack] = {}
    inventory: list[dict[str, object]] = []
    for sample_id in cohort.samples:
        mechanism = cohort.mechanism(sample_id)
        for kind in ("combined", "hap1", "hap2"):
            path = find_track(METHYLATION_DIR, sample_id, kind)
            required = kind == "combined" or mechanism in {"Control", "DiGeorge", "PWS-mUPD"}
            if path is None:
                inventory.append(
                    {
                        "sample_id": sample_id,
                        "mechanism": mechanism,
                        "track_kind": kind,
                        "path": None,
                        "required": required,
                        "status": "missing",
                    }
                )
                if required:
                    raise FileNotFoundError(f"Missing required {kind} methylation track for {sample_id}")
                continue
            track = read_track(path, CHROM, DOMAIN_START, DOMAIN_END)
            track_metadata = read_track_metadata(path)
            tracks[(sample_id, kind)] = track
            inventory.append(
                {
                    "sample_id": sample_id,
                    "mechanism": mechanism,
                    "track_kind": kind,
                    "path": str(path),
                    "required": required,
                    "status": "loaded",
                    "domain_cpgs": len(track.position),
                    "domain_median_depth": (
                        float(np.median(track.coverage)) if len(track.coverage) else np.nan
                    ),
                    "reference_depth_common_cn1": reference_depth(track, common_interval),
                    "track_set": track_set_for(mechanism, kind),
                    "pbcpgtools_version": track_metadata.get("pb_cpg_tools_version"),
                    "pileup_mode": track_metadata.get("pileup_mode"),
                    "modsites_mode": track_metadata.get("modsites_mode"),
                    "min_coverage": track_metadata.get("min_coverage"),
                    "min_mapq": track_metadata.get("min_mapq"),
                }
            )
    inventory_table = pd.DataFrame(inventory)
    for setting in ("pileup_mode", "modsites_mode", "min_coverage", "min_mapq"):
        values = inventory_table[setting].dropna().astype(str).unique()
        if len(values) > 1:
            raise ValueError(f"Mixed pb-CpG-tools {setting} values: {sorted(values)}")
    loaded = inventory_table[inventory_table["status"].eq("loaded")]
    shared: dict[tuple[str, str], np.ndarray] = {}
    fractions: dict[tuple[str, str], float] = {}
    targets: dict[str, float] = {}
    for set_name, spec in TRACK_SETS.items():
        members_table = loaded[loaded["track_set"].eq(set_name)]
        members = list(zip(members_table["sample_id"], members_table["track_kind"]))
        if not members:
            continue
        shared.update(shared_positions(tracks, members, str(spec["shared_scope"])))
        depths = members_table["reference_depth_common_cn1"].astype(float)
        target = float(DOWNSAMPLE_TARGET_DEPTH) if DOWNSAMPLE_TARGET_DEPTH else float(depths.min())
        targets[set_name] = target
        for key, depth in zip(members, depths):
            fractions[key] = min(1.0, target / depth) if depth > 0 else 1.0
    inventory_table["downsample_target_depth"] = inventory_table["track_set"].map(targets)
    inventory_table["downsample_fraction"] = [
        fractions.get((sample, kind), np.nan)
        for sample, kind in zip(inventory_table["sample_id"], inventory_table["track_kind"])
    ]
    inventory_table["shared_cpgs"] = [
        len(shared[(sample, kind)]) if (sample, kind) in shared else np.nan
        for sample, kind in zip(inventory_table["sample_id"], inventory_table["track_kind"])
    ]
    rng = np.random.default_rng(DOWNSAMPLE_SEED)
    rows: list[pd.DataFrame] = []
    for sample_id in cohort.samples:
        mechanism = cohort.mechanism(sample_id)
        deletion_interval = deletion_map.interval(sample_id)
        for kind in ("combined", "hap1", "hap2"):
            key = (sample_id, kind)
            if key not in tracks:
                continue
            track = tracks[key]
            thinned = None
            mask = None
            if key in fractions:
                thinned = downsample_track(track, fractions[key], DOWNSAMPLE_MIN_COVERAGE, rng)
                mask = np.isin(track.position, shared[key], assume_unique=True)
            summary = summarize_track_windows(
                track, windows, DEPTH_CAPS, EFFECTIVE_COVERAGE_CAP, thinned, mask
            )
            for column in (
                "n_cpg_downsampled",
                "beta_downsampled",
                "effective_obs_downsampled",
                "n_cpg_shared",
                "beta_shared_cpg",
                "effective_obs_shared",
            ):
                if column not in summary:
                    summary[column] = np.nan
            annotations = [
                classify_evidence(
                    mechanism,
                    kind,
                    int(window.start),
                    int(window.end),
                    deletion_interval,
                    common_interval,
                )
                for window in summary.itertuples(index=False)
            ]
            annotation_table = pd.DataFrame(annotations)
            summary.insert(0, "sample_id", sample_id)
            summary.insert(1, "mechanism", mechanism)
            summary.insert(2, "track_kind", kind)
            summary["track_set"] = track_set_for(mechanism, kind)
            summary["downsample_fraction"] = fractions.get(key, np.nan)
            summary["source_path"] = str(track.source)
            summary = pd.concat([summary.reset_index(drop=True), annotation_table], axis=1)
            summary["measurement_available"] = summary["n_cpg"].ge(MIN_CPGS)
            summary["primary_measurement"] = (
                summary["eligible_primary"] & summary["measurement_available"]
            )
            rows.append(summary)
    matrix = pd.concat(rows, ignore_index=True)
    duplicated = matrix.duplicated(["sample_id", "track_kind", "window_id"])
    if duplicated.any():
        raise RuntimeError("Evidence matrix contains duplicate sample-track-window rows")
    outdir = OUTPUT_DIR
    outdir.mkdir(parents=True, exist_ok=True)
    paths = {
        "matrix": outdir / "chr15_window_evidence_matrix.tsv.gz",
        "inventory": outdir / "methylation_track_inventory.tsv",
        "deletions": outdir / "validated_cn1_intervals.tsv",
        "core": outdir / "common_reciprocal_cn1_core.tsv",
        "depth": outdir / "sample_track_depth_qc.tsv",
        "manifest": outdir / "analysis_manifest.json",
    }
    write_table(matrix, paths["matrix"])
    write_table(inventory_table, paths["inventory"])
    write_table(deletion_map.table, paths["deletions"])
    write_table(
        pd.DataFrame(
            [
                {
                    "chrom": CHROM,
                    "start": deletion_map.common_start,
                    "end": deletion_map.common_end,
                    "cn1_buffer": deletion_map.buffer,
                    "n_pws": len(cohort.samples_for("PWS-DEL")),
                    "n_as": len(cohort.samples_for("AS-DEL")),
                }
            ]
        ),
        paths["core"],
    )
    depth = (
        matrix.groupby(["sample_id", "mechanism", "track_kind"], as_index=False)
        .agg(
            observed_windows=("measurement_available", "sum"),
            total_windows=("window_id", "nunique"),
            median_cpgs_per_window=("n_cpg", "median"),
            median_depth_per_cpg=("median_depth", "median"),
            mean_beta_site=("beta_site_mean", "mean"),
            mean_beta_read_weighted=("beta_read_weighted", "mean"),
        )
    )
    depth["window_recovery_fraction"] = depth["observed_windows"] / depth["total_windows"]
    depth["estimator_absolute_difference"] = (
        depth["mean_beta_site"] - depth["mean_beta_read_weighted"]
    ).abs()
    write_table(depth, paths["depth"])
    manifest = {
        "analysis": "chr15_evidence_matrix",
        "coordinate_system": "T2T-CHM13v2.0; zero-based half-open",
        "chrom": CHROM,
        "start": DOMAIN_START,
        "end": DOMAIN_END,
        "window_size": WINDOW_SIZE,
        "minimum_cpgs": MIN_CPGS,
        "depth_caps": list(DEPTH_CAPS),
        "effective_coverage_cap": EFFECTIVE_COVERAGE_CAP,
        "downsample_targets": targets,
        "downsample_min_coverage": DOWNSAMPLE_MIN_COVERAGE,
        "downsample_seed": DOWNSAMPLE_SEED,
        "downsample_method": "per-CpG binomial read thinning with hypergeometric methylated-read draw",
        "shared_cpg_scopes": {name: spec["shared_scope"] for name, spec in TRACK_SETS.items()},
        "primary_beta": "unweighted mean of emitted CpG modification scores",
        "parental_direction_rule": "diagnosis-derived only for combined tracks wholly inside validated buffered CN=1 intervals",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "metadata": str(METADATA_PATH.resolve()),
        "methylation_dir": str(METHYLATION_DIR.resolve()),
        "deletion_provenance": str(DELETION_PROVENANCE_PATH.resolve()),
        "outputs": {key: str(path.resolve()) for key, path in paths.items() if key != "manifest"},
        "command": " ".join(sys.argv),
    }
    paths["manifest"].write_text(json.dumps(manifest, indent=2) + "\n")
    return paths


def main() -> None:
    paths = build_matrix()
    for name, path in paths.items():
        print(f"{name}\t{path}")


if __name__ == "__main__":
    main()
