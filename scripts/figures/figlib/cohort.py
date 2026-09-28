"""Cohort groups with the Figure 1 palette, markers and display labels, so every figure
names and colours participants the same way."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .paths import import_analysis_package

GROUPS = ("PWS-DEL", "AS-DEL", "PWS-mUPD", "DiGeorge", "Control")          # Figure 1 order
DISPLAY = {"PWS-DEL": "PWS-DEL", "AS-DEL": "AS-DEL", "PWS-mUPD": "PWS-mUPD",
           "DiGeorge": "22q11.2DS", "Control": "Control"}
COLOR = {"PWS-DEL": "#D55E00", "AS-DEL": "#0072B2", "PWS-mUPD": "#CC79A7",
         "DiGeorge": "#E69F00", "Control": "#4D4D4D"}
MARKER = {"PWS-DEL": "o", "AS-DEL": "s", "PWS-mUPD": "D", "DiGeorge": "^", "Control": "P"}
PREFIX = {"PWS-DEL": "PW", "AS-DEL": "AS", "PWS-mUPD": "UPD", "DiGeorge": "DC", "Control": "CTRL"}
RETAINED = {"PWS-DEL": "maternal", "AS-DEL": "paternal"}   # the copy left by the deletion
BIPARENTAL = ("Control", "DiGeorge")
DELETIONS = ("PWS-DEL", "AS-DEL")


@dataclass(frozen=True)
class Cohort:
    table: pd.DataFrame        # sample_id, mechanism, label

    @property
    def samples(self) -> tuple[str, ...]:
        return tuple(self.table["sample_id"])

    def group(self, sample: str) -> str:
        return str(self.table.set_index("sample_id").at[sample, "mechanism"])

    def label(self, sample: str) -> str:
        return str(self.table.set_index("sample_id").at[sample, "label"])

    def of(self, *groups: str) -> tuple[str, ...]:
        return tuple(self.table.loc[self.table["mechanism"].isin(groups), "sample_id"])

    def groups_present(self) -> list[str]:
        return [g for g in GROUPS if (self.table["mechanism"] == g).any()]


def load(metadata: Path) -> Cohort:
    """Metadata through scripts/analysis (same normalisation of group names); display
    labels PW-1, AS-1, UPD-1, DC-1, CTRL-1 in metadata order, as in analysis 01."""
    table = import_analysis_package().load_cohort(Path(metadata)).table[["sample_id", "mechanism"]].copy()
    counts: dict[str, int] = {}
    labels = []
    for mech in table["mechanism"]:
        counts[mech] = counts.get(mech, 0) + 1
        labels.append(f"{PREFIX.get(mech, mech)}-{counts[mech]}")
    table["label"] = labels
    table["order"] = table["mechanism"].map({g: i for i, g in enumerate(GROUPS)})
    return Cohort(table.sort_values(["order", "sample_id"]).drop(columns="order").reset_index(drop=True))
