#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path


# ============================== CONFIGURATION ==============================
PROJECT_ROOT = Path(__file__).resolve().parents[2]
GTF_CANDIDATES = (
    Path("/home/rare/arlen/reference/chm13v22.sorted.gtf"),
    PROJECT_ROOT / "reference" / "genes.gtf",
)
OUTPUT_PATH = PROJECT_ROOT / "results" / "analysis" / "prespecified_regions.tsv"
ANNOTATION_FILENAME = "annotation_genes.tsv"
CHROM = "chr15"
REGION_PADDING = 5_000
IC_START = 22_691_258
IC_END = 22_693_494
IC_EXCLUSION_FLANK = 2_000
IC_REGION_ID = "PWS/AS imprinting centre"

REGION_DEFINITIONS = (
    ("MAGEL2/NDN", ("MAGEL2", "NDN"), "imprinted_gene_block"),
    ("SNRPN/SNHG14", ("SNURF", "SNRPN"), "SNRPN_body_downstream_of_IC"),
    ("SNORD116", ("SNORD116",), "snoRNA_cluster"),
    ("UBE3A", ("UBE3A",), "imprinted_gene"),
    ("GABRB3/GABA receptor cluster", ("GABRB3", "GABRA5", "GABRG3"), "regional_context"),
    ("OCA2 downstream control", ("OCA2",), "non_imprinted_control"),
)
REGION_ORDER = (
    "MAGEL2/NDN",
    IC_REGION_ID,
    "SNRPN/SNHG14",
    "SNORD116",
    "UBE3A",
    "GABRB3/GABA receptor cluster",
    "OCA2 downstream control",
)
ANNOTATION_GENES = (
    ("MKRN3", ("MKRN3",)),
    ("MAGEL2", ("MAGEL2",)),
    ("NDN", ("NDN",)),
    ("SNRPN", ("SNRPN",)),
    ("SNORD116", ("SNORD116",)),
    ("SNORD115", ("SNORD115",)),
    ("UBE3A", ("UBE3A",)),
    ("ATP10A", ("ATP10A",)),
    ("GABRB3", ("GABRB3",)),
    ("GABRA5", ("GABRA5",)),
    ("GABRG3", ("GABRG3",)),
    ("OCA2", ("OCA2",)),
    ("HERC2", ("HERC2",)),
)
# ===========================================================================

GTF_PATH = next((path for path in GTF_CANDIDATES if path.is_file()), GTF_CANDIDATES[0])


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def attribute_value(attributes: str, key: str) -> str | None:
    match = re.search(rf'(?:^|;)\s*{re.escape(key)}\s+["\']?([^;"\']+)', attributes)
    return match.group(1).strip() if match else None


def load_genes() -> dict[str, tuple[int, int]]:
    genes: dict[str, tuple[int, int]] = {}
    with open_text(GTF_PATH) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[0] != CHROM or fields[2] not in {"gene", "transcript"}:
                continue
            name = attribute_value(fields[8], "gene_name") or attribute_value(fields[8], "gene_id")
            if not name:
                continue
            start = int(fields[3]) - 1
            end = int(fields[4])
            if name in genes:
                previous_start, previous_end = genes[name]
                genes[name] = min(previous_start, start), max(previous_end, end)
            else:
                genes[name] = start, end
    if not genes:
        raise ValueError(f"No {CHROM} genes found in {GTF_PATH}")
    return genes


def matching_genes(
    genes: dict[str, tuple[int, int]],
    queries: tuple[str, ...],
) -> list[tuple[str, int, int]]:
    selected: list[tuple[str, int, int]] = []
    for name, (start, end) in genes.items():
        keep = any(
            name.startswith(query) if query in {"SNORD116", "SNORD115"} else name == query
            for query in queries
        )
        if keep:
            selected.append((name, start, end))
    return sorted(selected)


def carve_ic(row: dict[str, object]) -> dict[str, object]:
    low, high = IC_START - IC_EXCLUSION_FLANK, IC_END + IC_EXCLUSION_FLANK
    start, end = int(row["start"]), int(row["end"])
    if end <= low or start >= high:
        return row
    left, right = (start, low), (high, end)
    keep = max((left, right), key=lambda interval: interval[1] - interval[0])
    if keep[1] <= keep[0]:
        raise ValueError(f"{row['region_id']} lies entirely within the imprinting centre")
    return {**row, "start": keep[0], "end": keep[1]}


def build_regions() -> list[dict[str, object]]:
    if not GTF_PATH.is_file():
        raise FileNotFoundError(f"T2T-CHM13 GTF not found: {GTF_PATH}")
    genes = load_genes()
    rows: list[dict[str, object]] = [
        {
            "region_id": IC_REGION_ID,
            "chrom": CHROM,
            "start": IC_START,
            "end": IC_END,
            "region_class": "imprinting_control_region",
            "member_genes": "IC",
            "source": "T2T-CHM13 coordinate specified by study",
        }
    ]
    for region_id, queries, region_class in REGION_DEFINITIONS:
        selected = matching_genes(genes, queries)
        if not selected:
            raise ValueError(f"No GTF feature found for prespecified region {region_id}")
        row = {
            "region_id": region_id,
            "chrom": CHROM,
            "start": max(0, min(start for _, start, _ in selected) - REGION_PADDING),
            "end": max(end for _, _, end in selected) + REGION_PADDING,
            "region_class": region_class,
            "member_genes": ",".join(name for name, _, _ in selected),
            "source": str(GTF_PATH),
        }
        rows.append(carve_ic(row))
    rows.sort(key=lambda row: (int(row["start"]), int(row["end"])))
    identifiers = [str(row["region_id"]) for row in rows]
    invalid = len(identifiers) != len(set(identifiers)) or any(
        int(row["start"]) >= int(row["end"]) for row in rows
    )
    if invalid:
        raise RuntimeError("Invalid prespecified region catalog")
    for previous, current in zip(rows, rows[1:]):
        if int(current["start"]) < int(previous["end"]):
            raise ValueError(
                f"Prespecified regions overlap: {previous['region_id']} and {current['region_id']}"
            )
    for row in rows:
        row["display_order"] = REGION_ORDER.index(str(row["region_id"])) + 1
    return rows


def build_annotation() -> list[dict[str, object]]:
    genes = load_genes()
    rows: list[dict[str, object]] = []
    for label, queries in ANNOTATION_GENES:
        selected = matching_genes(genes, queries)
        if not selected:
            continue
        rows.append(
            {
                "label": label,
                "chrom": CHROM,
                "start": min(start for _, start, _ in selected),
                "end": max(end for _, _, end in selected),
                "n_features": len(selected),
            }
        )
    rows.append(
        {"label": "IC", "chrom": CHROM, "start": IC_START, "end": IC_END, "n_features": 1}
    )
    return sorted(rows, key=lambda row: int(row["start"]))


def main() -> None:
    rows = build_regions()
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = (
        "region_id",
        "chrom",
        "start",
        "end",
        "region_class",
        "member_genes",
        "source",
        "display_order",
    )
    with OUTPUT_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    annotation_path = OUTPUT_PATH.with_name(ANNOTATION_FILENAME)
    with annotation_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("label", "chrom", "start", "end", "n_features"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(build_annotation())
    print(OUTPUT_PATH)
    print(annotation_path)


if __name__ == "__main__":
    main()
