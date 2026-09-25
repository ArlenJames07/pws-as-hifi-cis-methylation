#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
GTF_PATH = Path("/home/rare/arlen/reference/chm13v22.sorted.gtf")
OUTPUT_PATH = PROJECT_ROOT / "results" / "analysis" / "prespecified_regions.tsv"
CHROM = "chr15"
REGION_PADDING = 5_000
IC_START = 22_691_258
IC_END = 22_693_494

REGION_DEFINITIONS = (
    ("MAGEL2/NDN", ("MAGEL2", "NDN"), "imprinted_gene_block"),
    ("SNRPN/SNHG14", ("SNRPN", "SNHG14"), "imprinting_domain"),
    ("SNORD116", ("SNORD116",), "snoRNA_cluster"),
    ("UBE3A", ("UBE3A",), "imprinted_gene"),
    ("GABRB3/GABA receptor", ("GABRB3", "GABRA5", "GABRG3"), "regional_context"),
    ("OCA2", ("OCA2",), "non_imprinted_control"),
)


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
            name.startswith(query) if query == "SNORD116" else name == query
            for query in queries
        )
        if keep:
            selected.append((name, start, end))
    return sorted(selected)


def build_regions() -> list[dict[str, object]]:
    if not GTF_PATH.is_file():
        raise FileNotFoundError(f"T2T-CHM13 GTF not found: {GTF_PATH}")
    genes = load_genes()
    rows: list[dict[str, object]] = [
        {
            "region_id": "PWS/AS imprinting centre",
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
        rows.append(
            {
                "region_id": region_id,
                "chrom": CHROM,
                "start": max(0, min(start for _, start, _ in selected) - REGION_PADDING),
                "end": max(end for _, _, end in selected) + REGION_PADDING,
                "region_class": region_class,
                "member_genes": ",".join(name for name, _, _ in selected),
                "source": str(GTF_PATH),
            }
        )
    rows.sort(key=lambda row: (int(row["start"]), int(row["end"])))
    identifiers = [str(row["region_id"]) for row in rows]
    invalid = len(identifiers) != len(set(identifiers)) or any(
        int(row["start"]) >= int(row["end"]) for row in rows
    )
    if invalid:
        raise RuntimeError("Invalid prespecified region catalog")
    return rows


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
    )
    with OUTPUT_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(OUTPUT_PATH)


if __name__ == "__main__":
    main()
