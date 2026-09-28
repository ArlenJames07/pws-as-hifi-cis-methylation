"""
Repeat annotation readers shared by analysis scripts 05 and 06.

Accepted formats (auto-detected):
  * RepeatMasker .out (query names may be samtools-faidx windows 'chr15:S-E',
    which are converted to chromosome coordinates)
  * UCSC rmsk.txt (hs1 / CHM13)
  * BED: chrom start end name [class [family]]; a numeric column 5 is taken as a
    score, and the class is then read from column 7 when present
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

CLASSES = ("SINE", "LINE", "LTR", "DNA", "Simple_repeat", "Low_complexity", "Satellite", "other")
COLUMNS = ["start", "end", "name", "class", "family", "strand", "div"]


def class_group(value: str) -> str:
    head = str(value).split("/")[0].replace("?", "")
    return head if head in CLASSES else "other"


def _numeric(text: str) -> bool:
    return re.match(r"^-?[\d.]+(e-?\d+)?$", str(text)) is not None


def read_repeatmasker_out(path: Path | str) -> pd.DataFrame:
    """RepeatMasker .out -> DataFrame(seq, start (0-based), end, name, class, family, strand, div)."""
    rows = []
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 11 or not fields[0].isdigit():
                continue
            class_family = fields[10].split("/")
            rows.append((fields[4], int(fields[5]) - 1, int(fields[6]), fields[9],
                         class_group(fields[10]),
                         class_family[1] if len(class_family) > 1 else class_family[0],
                         "+" if fields[8] == "+" else "-", float(fields[1])))
    return pd.DataFrame(rows, columns=["seq"] + COLUMNS)


def read_repeat_annotation(path: Path | str, chrom: str, start: int, end: int) -> pd.DataFrame:
    """Elements overlapping chrom:[start, end) in chromosome coordinates."""
    path = Path(path)
    with open(path) as handle:
        head = handle.readline()
    if "SW" in head and "perc" in head:
        table = read_repeatmasker_out(path)
        converted = []
        for seq, a, b in zip(table["seq"], table["start"], table["end"]):
            match = re.match(r"^(.+):(\d+)-(\d+)$", seq)
            if match:
                offset = int(match.group(2)) - 1
                converted.append((match.group(1), a + offset, b + offset))
            else:
                converted.append((seq, a, b))
        table[["chrom", "start", "end"]] = pd.DataFrame(converted, index=table.index)
    else:
        rows = []
        with open(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) >= 17 and f[5].startswith("chr") and _numeric(f[6]):     # UCSC rmsk.txt
                    rows.append((f[5], int(f[6]), int(f[7]), f[10], class_group(f[11]), f[12],
                                 f[9], int(f[2]) / 10.0))
                    continue
                if len(f) < 3:
                    continue
                name = f[3] if len(f) > 3 else "."
                if len(f) > 4 and not _numeric(f[4]):
                    cls, fam = f[4], f[5] if len(f) > 5 else f[4]
                elif len(f) > 6 and not _numeric(f[6]):
                    cls, fam = f[6], f[7] if len(f) > 7 else f[6]
                else:
                    cls, fam = "unclassified", "unclassified"
                strand = f[5] if len(f) > 5 and f[5] in "+-" else "."
                rows.append((f[0], int(f[1]), int(f[2]), name, class_group(cls) if cls != "unclassified" else cls,
                             fam, strand, float("nan")))
        table = pd.DataFrame(rows, columns=["chrom"] + COLUMNS)
    table = table[(table["chrom"] == chrom) & (table["end"] > start) & (table["start"] < end)]
    table = table.sort_values(["start", "end"]).reset_index(drop=True)
    return table[["chrom"] + COLUMNS]
