#!/usr/bin/env python3
from __future__ import annotations

import gzip
import json
import platform
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

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
CNV_DIR = PROJECT_ROOT / "results" / "05_cnv"
PHASING_DIR = PROJECT_ROOT / "results" / "04_phasing"
OUTPUT_DIR = PROJECT_ROOT / "results" / "analysis" / "01_evidence_matrix"

# Canonical structural evidence produced by analysis. Figures consume this table;
# figures must never be the source of structural/deletion provenance.
STRUCTURAL_EVIDENCE_PATH = OUTPUT_DIR / "chr15_structural_evidence.tsv"
CN_SEGMENTS_PATH = OUTPUT_DIR / "chr15_copy_number_segments.tsv.gz"

CHROM = "chr15"
DOMAIN_START = 22_000_000
DOMAIN_END = 28_000_000
WINDOW_SIZE = 1_000

PWS_IC_START = 22_691_258
PWS_IC_END = 22_693_494

CN1_BREAKPOINT_BUFFER = 75_000
MIN_CPGS = 3
DEPTH_CAPS = (5.0, 10.0, 15.0, 20.0)

# HiFiCNV thresholds used only to establish structural provenance.
CN_DELETION_THRESHOLD = 1.35
CN1_MEAN_MIN = 0.65
CN1_MEAN_MAX = 1.35
CN_MERGE_GAP_BP = 150_000
BP_NEAREST_MAX_DISTANCE = 500_000
BREAKPOINT_LANDMARKS = {
    "BP1": 20_940_000,
    "BP2": 21_070_000,
    "BP3": 26_050_000,
    "BP4": 26_460_000,
    "BP5": 31_840_000,
}


@dataclass(frozen=True)
class CNSegment:
    chrom: str
    start: int
    end: int
    copy_number: float
    source: str


def write_table(table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if path.suffix == ".gz" else None
    table.to_csv(path, sep="\t", index=False, na_rep="NA", compression=compression)


def _open_text_maybe_gzip(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open("rt")


def _choose_file(files: list[Path], sample_id: str) -> Path | None:
    if not files:
        return None
    exact = [
        path
        for path in files
        if re.search(rf"(?:^|[_-]){re.escape(sample_id)}(?:\.|_|-|$)", path.name)
    ]
    candidates = exact or files
    return sorted(candidates, key=lambda p: (len(p.name), p.name))[0]


def _find_sample_file(directory: Path, sample_id: str, suffixes: tuple[str, ...]) -> Path | None:
    if not directory.exists():
        return None
    sample_dir = directory / sample_id
    search_root = sample_dir if sample_dir.is_dir() else directory
    matches: list[Path] = []
    for suffix in suffixes:
        matches.extend(search_root.rglob(f"*{sample_id}*{suffix}"))
    return _choose_file(matches, sample_id)


def find_hificnv_cn_track(directory: Path, sample_id: str) -> Path | None:
    return _find_sample_file(
        directory,
        sample_id,
        (
            ".copynum.bedgraph",
            ".copynum.bedgraph.gz",
            ".cnv.bed",
            ".cnv.bed.gz",
            ".bedgraph",
            ".bedgraph.gz",
        ),
    )


def find_sv_vcf(directory: Path, sample_id: str) -> Path | None:
    return _find_sample_file(
        directory,
        sample_id,
        (
            ".sv.phased.vcf.gz",
            ".sv.pass.vcf.gz",
            ".sv.vcf.gz",
            ".sv.phased.vcf",
            ".sv.pass.vcf",
            ".sv.vcf",
        ),
    )


def _extract_cn_from_fields(fields: list[str], path: Path) -> float | None:
    tail = fields[3:]
    tail_text = " ".join(tail)

    match = re.search(
        r"(?:^|[;\s])(?:CN|copy[_ -]?number|copynum)\s*[=:]\s*([0-9]+(?:\.[0-9]+)?)",
        tail_text,
        flags=re.IGNORECASE,
    )
    if match:
        return float(match.group(1))

    lower_name = path.name.lower()
    if ("copynum" in lower_name or "bedgraph" in lower_name) and len(fields) >= 4:
        try:
            value = float(fields[3])
        except ValueError:
            value = float("nan")
        if np.isfinite(value) and 0 <= value <= 10:
            return value

    # Legacy HiFiCNV/CNV BED fallback. Restrict the heuristic to dosage files.
    if "hificnv" in lower_name or "cnv" in lower_name:
        for token in tail:
            try:
                value = float(token.strip(",;"))
            except ValueError:
                continue
            if 0 <= value <= 6 and abs(value - round(value)) <= 0.05:
                return float(value)
    return None


def read_hificnv_cn_segments(path: Path | None) -> list[CNSegment]:
    if path is None or not path.exists():
        return []

    segments: list[CNSegment] = []
    try:
        with _open_text_maybe_gzip(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip().split()
                if len(fields) < 4 or fields[0] != CHROM:
                    continue
                try:
                    start = int(float(fields[1]))
                    end = int(float(fields[2]))
                except ValueError:
                    continue
                if end <= start:
                    continue
                cn = _extract_cn_from_fields(fields, path)
                if cn is None or not np.isfinite(cn):
                    continue
                segments.append(
                    CNSegment(
                        chrom=CHROM,
                        start=start,
                        end=end,
                        copy_number=float(cn),
                        source=str(path),
                    )
                )
    except (OSError, EOFError):
        return []

    return sorted(segments, key=lambda s: (s.start, s.end))


def _merge_deleted_cn_segments(segments: list[CNSegment]) -> list[dict[str, Any]]:
    loss = [seg for seg in segments if seg.copy_number <= CN_DELETION_THRESHOLD]
    if not loss:
        return []

    events: list[dict[str, Any]] = []
    current: dict[str, Any] = {
        "start": loss[0].start,
        "end": loss[0].end,
        "segments": [loss[0]],
    }
    for seg in loss[1:]:
        if seg.start <= int(current["end"]) + CN_MERGE_GAP_BP:
            current["end"] = max(int(current["end"]), seg.end)
            current["segments"].append(seg)
        else:
            events.append(current)
            current = {"start": seg.start, "end": seg.end, "segments": [seg]}
    events.append(current)

    for event in events:
        segs: list[CNSegment] = event["segments"]
        lengths = np.asarray([seg.end - seg.start for seg in segs], dtype=float)
        cns = np.asarray([seg.copy_number for seg in segs], dtype=float)
        event["mean_cn"] = (
            float(np.average(cns, weights=lengths))
            if lengths.sum() > 0
            else float(np.mean(cns))
        )
        event["min_cn"] = float(np.min(cns))
        event["n_segments"] = len(segs)
    return events


def _parse_info_field(info: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for token in info.split(";"):
        if not token:
            continue
        if "=" in token:
            key, value = token.split("=", 1)
            parsed[key] = value
        else:
            parsed[token] = "True"
    return parsed


def sv_vcf_deletion_overlaps_ic(path: Path | None) -> tuple[bool | None, str]:
    if path is None or not path.exists():
        return None, "SV VCF unavailable"
    try:
        with _open_text_maybe_gzip(path) as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) < 8 or fields[0] != CHROM:
                    continue
                try:
                    pos = int(fields[1])
                except ValueError:
                    continue
                info = _parse_info_field(fields[7])
                if info.get("SVTYPE", "").upper() != "DEL":
                    continue
                try:
                    if "END" in info:
                        end = int(info["END"])
                    else:
                        svlen = abs(int(info.get("SVLEN", "0").split(",")[0]))
                        end = pos + svlen
                except ValueError:
                    continue
                left, right = min(pos, end), max(pos, end)
                if left <= PWS_IC_START and right >= PWS_IC_END:
                    return True, f"SV DEL {CHROM}:{left}-{right}"
        return False, "No SV DEL spanning complete IC"
    except (OSError, EOFError) as exc:
        return None, f"SV VCF parse error: {exc}"


def cnv_copy_number_evidence(
    path: Path | None,
) -> tuple[bool | None, str, dict[str, Any] | None, list[CNSegment]]:
    if path is None or not path.exists():
        return None, "HiFiCNV copy-number track unavailable", None, []

    segments = read_hificnv_cn_segments(path)
    if not segments:
        return None, "No parseable chr15 copy-number records", None, []

    events = _merge_deleted_cn_segments(segments)
    spanning = [
        event
        for event in events
        if int(event["start"]) <= PWS_IC_START and int(event["end"]) >= PWS_IC_END
    ]
    if spanning:
        event = max(spanning, key=lambda e: int(e["end"]) - int(e["start"]))
        mean_cn = float(event["mean_cn"])
        cn1_ok = CN1_MEAN_MIN <= mean_cn <= CN1_MEAN_MAX
        note = (
            f"HiFiCNV CN-loss {CHROM}:{event['start']}-{event['end']} "
            f"(mean CN={mean_cn:.3f}; required CN1 range "
            f"{CN1_MEAN_MIN:.2f}-{CN1_MEAN_MAX:.2f})"
        )
        return cn1_ok, note, event, segments

    overlapping = [
        seg for seg in segments if seg.start < PWS_IC_END and seg.end > PWS_IC_START
    ]
    if overlapping:
        lengths = np.asarray(
            [
                max(0, min(seg.end, PWS_IC_END) - max(seg.start, PWS_IC_START))
                for seg in overlapping
            ],
            dtype=float,
        )
        cns = np.asarray([seg.copy_number for seg in overlapping], dtype=float)
        mean_cn = (
            float(np.average(cns, weights=lengths))
            if lengths.sum() > 0
            else float(np.mean(cns))
        )
        return False, f"IC mean CN={mean_cn:.3f}; no CN1 event spans complete IC", None, segments

    return None, "CN track has chr15 records but no bins overlapping the IC", None, segments


def _nearest_breakpoint(position: int, candidates: tuple[str, ...]) -> tuple[str | None, int | None]:
    distances = {
        name: abs(position - BREAKPOINT_LANDMARKS[name])
        for name in candidates
    }
    if not distances:
        return None, None
    name = min(distances, key=distances.get)
    distance = distances[name]
    if distance > BP_NEAREST_MAX_DISTANCE:
        return None, distance
    return name, distance


def classify_deletion_from_cn_event(event: dict[str, Any] | None) -> dict[str, Any]:
    if event is None:
        return {
            "deletion_type": "no chr15 CN1 deletion",
            "left_landmark": "",
            "right_landmark": "",
            "left_distance_bp": np.nan,
            "right_distance_bp": np.nan,
        }

    left_name, left_dist = _nearest_breakpoint(int(event["start"]), ("BP1", "BP2"))
    right_name, right_dist = _nearest_breakpoint(int(event["end"]), ("BP3", "BP4", "BP5"))
    if left_name == "BP1" and right_name == "BP3":
        deletion_type = "Type I-like (BP1-BP3)"
    elif left_name == "BP2" and right_name == "BP3":
        deletion_type = "Type II-like (BP2-BP3)"
    elif left_name and right_name:
        deletion_type = f"atypical/extended {left_name}-{right_name}-like"
    else:
        deletion_type = "atypical/extended deletion"

    return {
        "deletion_type": deletion_type,
        "left_landmark": left_name or "unassigned",
        "right_landmark": right_name or "unassigned",
        "left_distance_bp": left_dist if left_dist is not None else np.nan,
        "right_distance_bp": right_dist if right_dist is not None else np.nan,
    }


def _display_labels(cohort) -> dict[str, str]:
    prefixes = {
        "PWS-DEL": "PW",
        "AS-DEL": "AS",
        "PWS-mUPD": "UPD",
        "DiGeorge": "DC",
        "Disease control": "DC",
        "Control": "CTRL",
    }
    counts: dict[str, int] = {}
    labels: dict[str, str] = {}
    for sample_id in cohort.samples:
        mechanism = cohort.mechanism(sample_id)
        counts[mechanism] = counts.get(mechanism, 0) + 1
        labels[sample_id] = f"{prefixes.get(mechanism, mechanism)}-{counts[mechanism]}"
    return labels


def build_structural_evidence(cohort) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Derive structural provenance from primary CNV/SV outputs.

    Crucially, parent-of-origin interpretation in downstream cis analysis is
    permitted only when HiFiCNV independently confirms a CN1 interval spanning
    the IC. pbsv is orthogonal breakpoint support, not a substitute for CN1.
    """
    labels = _display_labels(cohort)
    rows: list[dict[str, Any]] = []
    segment_rows: list[dict[str, Any]] = []

    for sample_id in cohort.samples:
        mechanism = cohort.mechanism(sample_id)
        expected_deletion = mechanism in {"PWS-DEL", "AS-DEL"}
        cn_track = find_hificnv_cn_track(CNV_DIR, sample_id)
        sv_vcf = find_sv_vcf(PHASING_DIR, sample_id)

        cn_ok, cn_note, cn_event, segments = cnv_copy_number_evidence(cn_track)
        sv_ok, sv_note = sv_vcf_deletion_overlaps_ic(sv_vcf)

        cn1_eligible = bool(cn_ok is True and cn_event is not None)
        if expected_deletion:
            if cn1_eligible:
                status = "confirmed"
            elif cn_ok is None:
                status = "unavailable"
            else:
                status = "not_confirmed"
        else:
            status = "not_expected"

        if cn1_eligible and sv_ok is True:
            evidence_basis = "concordant HiFiCNV CN1 + pbsv DEL"
        elif cn1_eligible:
            evidence_basis = "HiFiCNV CN1 dosage support; pbsv DEL not required"
        elif sv_ok is True:
            evidence_basis = "pbsv DEL support only; CN1 not independently confirmed"
        else:
            evidence_basis = "no confirmed chr15 CN1 deletion"

        deletion_class = classify_deletion_from_cn_event(cn_event if cn1_eligible else None)
        row = {
            "sample_id": sample_id,
            "display_label": labels[sample_id],
            "clinical_diagnosis": "",
            "molecular_mechanism": mechanism,
            "cn_track": str(cn_track) if cn_track else "",
            "sv_vcf": str(sv_vcf) if sv_vcf else "",
            "expected_ic_deletion": expected_deletion,
            "ic_deletion_status": status,
            "evidence_basis": evidence_basis,
            "sv_support": sv_ok,
            "sv_note": sv_note,
            "cnv_support": cn_ok,
            "cnv_note": cn_note,
            "cn_event_start": int(cn_event["start"]) if cn1_eligible else np.nan,
            "cn_event_end": int(cn_event["end"]) if cn1_eligible else np.nan,
            "cn_event_size_mb": (
                (int(cn_event["end"]) - int(cn_event["start"])) / 1e6
                if cn1_eligible
                else np.nan
            ),
            "cn_event_mean_cn": float(cn_event["mean_cn"]) if cn1_eligible else np.nan,
            "cn1_eligible": cn1_eligible,
            **deletion_class,
            # Compatibility aliases for existing Figure 1/report code.
            "pbsv_support": sv_ok,
            "pbsv_note": sv_note,
            "hificnv_support": cn_ok,
            "hificnv_note": cn_note,
        }
        rows.append(row)

        for seg in segments:
            segment_rows.append(
                {
                    "sample_id": sample_id,
                    "display_label": labels[sample_id],
                    "molecular_mechanism": mechanism,
                    "chrom": seg.chrom,
                    "start": seg.start,
                    "end": seg.end,
                    "copy_number": seg.copy_number,
                    "source_file": seg.source,
                }
            )

    structural = pd.DataFrame(rows)
    cn_segments = pd.DataFrame(segment_rows)

    # All deletion participants must have independently validated CN1 before
    # their retained combined methylation can be interpreted directionally.
    deletion_mask = structural["molecular_mechanism"].isin(["PWS-DEL", "AS-DEL"])
    invalid = structural.loc[deletion_mask & ~structural["cn1_eligible"].astype(bool)]
    if not invalid.empty:
        details = "; ".join(
            f"{row.sample_id} ({row.ic_deletion_status}: {row.cnv_note})"
            for row in invalid.itertuples(index=False)
        )
        raise RuntimeError(
            "Directional retained-allele analysis requires independently confirmed "
            f"HiFiCNV CN1 intervals for every PWS/AS deletion sample. Failed: {details}"
        )

    write_table(structural, STRUCTURAL_EVIDENCE_PATH)
    write_table(cn_segments, CN_SEGMENTS_PATH)
    return structural, cn_segments


def track_kinds(mechanism: str) -> tuple[str, ...]:
    # Keep all three tracks in the evidence inventory. Deletion haplotype tracks
    # remain technical-only evidence and are never parentally interpreted.
    return "combined", "hap1", "hap2"


def build_matrix() -> dict[str, Path]:
    if not METADATA_PATH.is_file():
        raise FileNotFoundError(f"Metadata not found: {METADATA_PATH}")
    if not METHYLATION_DIR.is_dir():
        raise FileNotFoundError(f"Methylation directory not found: {METHYLATION_DIR}")
    if not CNV_DIR.is_dir():
        raise FileNotFoundError(f"CNV directory not found: {CNV_DIR}")

    cohort = load_cohort(METADATA_PATH)

    # Analysis, not Figure 1, is now the authoritative producer of structural
    # provenance. Figure 1 reads STRUCTURAL_EVIDENCE_PATH generated here.
    build_structural_evidence(cohort)

    deletion_samples = cohort.samples_for("PWS-DEL") + cohort.samples_for("AS-DEL")
    deletion_map = load_deletion_map(
        STRUCTURAL_EVIDENCE_PATH,
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
                    raise FileNotFoundError(
                        f"Missing required {kind} methylation track for {sample_id}"
                    )
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
            summary = pd.concat(
                [summary.reset_index(drop=True), annotation_table], axis=1
            )
            summary["measurement_available"] = summary["n_cpg"].ge(MIN_CPGS)
            summary["primary_measurement"] = (
                summary["eligible_primary"] & summary["measurement_available"]
            )
            rows.append(summary)

    if not rows:
        raise RuntimeError("No methylation evidence rows were generated")

    matrix = pd.concat(rows, ignore_index=True)
    inventory_table = pd.DataFrame(inventory)

    for setting in ("pileup_mode", "modsites_mode", "min_coverage", "min_mapq"):
        if setting not in inventory_table.columns:
            continue
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
        "structural": STRUCTURAL_EVIDENCE_PATH,
        "cn_segments": CN_SEGMENTS_PATH,
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
        "structural_provenance_rule": (
            "HiFiCNV CN1 spanning the IC is required for directional retained-allele "
            "interpretation; pbsv is orthogonal breakpoint support"
        ),
        "parental_direction_rule": (
            "diagnosis-derived only for combined tracks wholly inside independently "
            "validated buffered CN1 intervals"
        ),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "metadata": str(METADATA_PATH.resolve()),
        "methylation_dir": str(METHYLATION_DIR.resolve()),
        "cnv_dir": str(CNV_DIR.resolve()),
        "phasing_dir": str(PHASING_DIR.resolve()),
        "structural_evidence": str(STRUCTURAL_EVIDENCE_PATH.resolve()),
        "outputs": {
            key: str(path.resolve())
            for key, path in paths.items()
            if key != "manifest"
        },
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
