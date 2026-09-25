from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class DeletionMap:
    table: pd.DataFrame
    common_start: int
    common_end: int
    buffer: int

    def interval(self, sample_id: str) -> tuple[int, int] | None:
        rows = self.table[self.table["sample_id"] == sample_id]
        if rows.empty:
            return None
        row = rows.iloc[0]
        return int(row["analysis_start"]), int(row["analysis_end"])


def load_deletion_map(path: Path, deletion_samples: tuple[str, ...], buffer: int) -> DeletionMap:
    if not path.is_file():
        raise FileNotFoundError(path)
    table = pd.read_csv(path, sep="\t")
    required = {
        "sample_id",
        "ic_deletion_status",
        "cn_event_start",
        "cn_event_end",
        "cn_event_mean_cn",
        "evidence_basis",
    }
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Deletion provenance lacks columns: {sorted(missing)}")
    table = table[table["sample_id"].isin(deletion_samples)].copy()
    if set(table["sample_id"]) != set(deletion_samples) or table["sample_id"].duplicated().any():
        raise ValueError("Deletion provenance must contain exactly one row per PWS/AS deletion sample")
    for column in ("cn_event_start", "cn_event_end", "cn_event_mean_cn"):
        table[column] = pd.to_numeric(table[column], errors="raise")
    invalid = (
        table["ic_deletion_status"].ne("confirmed")
        | table["cn_event_start"].ge(table["cn_event_end"])
        | table["cn_event_mean_cn"].lt(0.65)
        | table["cn_event_mean_cn"].gt(1.35)
        | table["evidence_basis"].isna()
    )
    if invalid.any():
        bad = table.loc[invalid, "sample_id"].tolist()
        raise ValueError(f"Invalid or unconfirmed CN=1 provenance: {bad}")
    table["analysis_start"] = table["cn_event_start"].astype(int) + int(buffer)
    table["analysis_end"] = table["cn_event_end"].astype(int) - int(buffer)
    if table["analysis_start"].ge(table["analysis_end"]).any():
        raise ValueError("Deletion buffer removes one or more CN=1 intervals")
    common_start = int(table["analysis_start"].max())
    common_end = int(table["analysis_end"].min())
    if common_end <= common_start:
        raise ValueError("No common reciprocal CN=1 interval remains")
    return DeletionMap(
        table=table.sort_values("sample_id").reset_index(drop=True),
        common_start=common_start,
        common_end=common_end,
        buffer=int(buffer),
    )


def classify_evidence(
    mechanism: str,
    track_kind: str,
    window_start: int,
    window_end: int,
    deletion_interval: tuple[int, int] | None,
    common_interval: tuple[int, int],
) -> dict[str, object]:
    in_sample_cn1 = False
    if deletion_interval is not None:
        in_sample_cn1 = window_start >= deletion_interval[0] and window_end <= deletion_interval[1]
    in_common_cn1 = window_start >= common_interval[0] and window_end <= common_interval[1]
    evidence_class = "not_primary"
    parental_state: str | None = None
    primary = False
    if mechanism == "PWS-DEL" and track_kind == "combined" and in_sample_cn1:
        evidence_class = "parental_observed_cn1"
        parental_state = "maternal_retained"
        primary = True
    elif mechanism == "AS-DEL" and track_kind == "combined" and in_sample_cn1:
        evidence_class = "parental_observed_cn1"
        parental_state = "paternal_retained"
        primary = True
    elif mechanism in {"Control", "DiGeorge"} and track_kind in {"hap1", "hap2"}:
        evidence_class = "diploid_unoriented"
        primary = True
    elif mechanism in {"Control", "DiGeorge"} and track_kind == "combined":
        evidence_class = "diploid_background"
        primary = True
    elif mechanism == "PWS-mUPD" and track_kind in {"hap1", "hap2"}:
        evidence_class = "maternal_origin_homologue_unoriented"
        parental_state = "maternal_origin"
        primary = True
    elif mechanism in {"PWS-DEL", "AS-DEL"} and track_kind in {"hap1", "hap2"}:
        evidence_class = "deletion_haplotype_technical_only"
    elif mechanism in {"PWS-DEL", "AS-DEL"} and track_kind == "combined":
        evidence_class = "outside_cn1_combined_unoriented"
    elif mechanism == "PWS-mUPD" and track_kind == "combined":
        evidence_class = "maternal_upd_background"
        parental_state = "maternal_origin"
        primary = True
    return {
        "inside_sample_cn1": bool(in_sample_cn1),
        "inside_common_cn1": bool(in_common_cn1),
        "evidence_class": evidence_class,
        "parental_state": parental_state,
        "eligible_primary": bool(primary),
        "parental_direction_observed": evidence_class == "parental_observed_cn1",
    }
