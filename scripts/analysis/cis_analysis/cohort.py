from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


MECHANISM_ALIASES = {
    "pws-del": "PWS-DEL",
    "pws_del": "PWS-DEL",
    "pws del": "PWS-DEL",
    "as-del": "AS-DEL",
    "as_del": "AS-DEL",
    "as del": "AS-DEL",
    "pws-mupd": "PWS-mUPD",
    "pws_mupd": "PWS-mUPD",
    "pws mupd": "PWS-mUPD",
    "digeorge": "DiGeorge",
    "disease control": "DiGeorge",
    "control": "Control",
    "unaffected control": "Control",
}


def normalize_mechanism(value: object) -> str:
    raw = str(value).strip()
    key = raw.lower()
    if key not in MECHANISM_ALIASES:
        raise ValueError(f"Unsupported molecular mechanism: {raw}")
    return MECHANISM_ALIASES[key]


@dataclass(frozen=True)
class Cohort:
    table: pd.DataFrame

    @property
    def samples(self) -> tuple[str, ...]:
        return tuple(self.table["sample_id"].astype(str))

    def mechanism(self, sample_id: str) -> str:
        rows = self.table.loc[self.table["sample_id"] == sample_id, "mechanism"]
        if len(rows) != 1:
            raise KeyError(sample_id)
        return str(rows.iloc[0])

    def samples_for(self, mechanism: str) -> tuple[str, ...]:
        return tuple(self.table.loc[self.table["mechanism"] == mechanism, "sample_id"])


def load_cohort(path: Path) -> Cohort:
    table = pd.read_csv(path)
    sample_col = "sample_id" if "sample_id" in table.columns else "sample"
    mechanism_col = next(
        (name for name in ("molecular_mechanism", "mechanism", "group") if name in table.columns),
        None,
    )
    if mechanism_col is None:
        raise ValueError("Metadata needs molecular_mechanism, mechanism, or group")
    out = table.copy()
    out["sample_id"] = out[sample_col].astype(str).str.strip()
    out["mechanism"] = out[mechanism_col].map(normalize_mechanism)
    if out["sample_id"].eq("").any() or out["sample_id"].duplicated().any():
        raise ValueError("Metadata sample IDs must be non-empty and unique")
    required = {"PWS-DEL", "AS-DEL", "Control", "DiGeorge"}
    missing = required - set(out["mechanism"])
    if missing:
        raise ValueError(f"Metadata lacks required cohort groups: {sorted(missing)}")
    keep = ["sample_id", "mechanism"] + [
        column for column in table.columns if column not in {sample_col, mechanism_col}
    ]
    return Cohort(out[keep].reset_index(drop=True))
