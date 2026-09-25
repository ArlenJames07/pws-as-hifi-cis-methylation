#!/usr/bin/env python3
from __future__ import annotations

import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from cis_analysis import (
    classify_evidence,
    find_track,
    fixed_windows,
    load_cohort,
    load_deletion_map,
    read_track,
    read_track_metadata,
    summarize_track_windows,
)


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
DOMAIN_START = 22_000_000
DOMAIN_END = 28_000_000
WINDOW_SIZE = 1_000
CN1_BREAKPOINT_BUFFER = 75_000
MIN_CPGS = 3
DEPTH_CAPS = (5.0, 10.0, 15.0, 20.0)


def write_table(table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if path.suffix == ".gz" else None
    table.to_csv(path, sep="\t", index=False, na_rep="NA", compression=compression)


def track_kinds(mechanism: str) -> tuple[str, ...]:
    if mechanism in {"Control", "DiGeorge", "PWS-mUPD"}:
        return "combined", "hap1", "hap2"
    return "combined", "hap1", "hap2"


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
    rows: list[pd.DataFrame] = []
    inventory: list[dict[str, object]] = []
    for sample_id in cohort.samples:
        mechanism = cohort.mechanism(sample_id)
        deletion_interval = deletion_map.interval(sample_id)
        for kind in track_kinds(mechanism):
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
                    "pbcpgtools_version": track_metadata.get("pb_cpg_tools_version"),
                    "pileup_mode": track_metadata.get("pileup_mode"),
                    "modsites_mode": track_metadata.get("modsites_mode"),
                    "min_coverage": track_metadata.get("min_coverage"),
                    "min_mapq": track_metadata.get("min_mapq"),
                }
            )
            summary = summarize_track_windows(track, windows, DEPTH_CAPS)
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
            summary["source_path"] = str(path)
            summary = pd.concat([summary.reset_index(drop=True), annotation_table], axis=1)
            summary["measurement_available"] = summary["n_cpg"].ge(MIN_CPGS)
            summary["primary_measurement"] = (
                summary["eligible_primary"] & summary["measurement_available"]
            )
            rows.append(summary)
    matrix = pd.concat(rows, ignore_index=True)
    inventory_table = pd.DataFrame(inventory)
    for setting in ("pileup_mode", "modsites_mode", "min_coverage", "min_mapq"):
        values = inventory_table[setting].dropna().astype(str).unique()
        if len(values) > 1:
            raise ValueError(f"Mixed pb-CpG-tools {setting} values: {sorted(values)}")
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
