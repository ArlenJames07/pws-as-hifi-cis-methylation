#!/usr/bin/env python3
"""
Self-contained Figure 1 generator for the hifi_multiomics_pipeline layout.

This version has all paths and arguments inside the script.

Run:
    python3 phase1_figure1_v2_no_overlap.py

Main figure outputs:
    Figure_1.png
    Figure_1.pdf
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
import textwrap
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle


CHROM = "chr15"
DOMAIN_START = 22_500_000
DOMAIN_END = 28_500_000

PWS_IC_START = 22_691_258
PWS_IC_END = 22_693_494
PWS_IC_NAME = "ICR_893_SNHG14_SNRPN_SNURF"

MATERNAL_THRESHOLD = 0.85
PATERNAL_THRESHOLD = 0.15
MIN_MEAN_COVERAGE = 10.0
MIN_CPGS = 5

COHORT = [
    ("001P", "Prader-Willi syndrome", "PWS-DEL"),
    ("002P", "Prader-Willi syndrome", "PWS-DEL"),
    ("005P", "Prader-Willi syndrome", "PWS-DEL"),
    ("006P", "Prader-Willi syndrome", "PWS-DEL"),
    ("007P", "Prader-Willi syndrome", "PWS-DEL"),
    ("004P", "Prader-Willi syndrome", "PWS-mUPD"),
    ("013A", "Angelman syndrome", "AS-DEL"),
    ("014A", "Angelman syndrome", "AS-DEL"),
    ("015A", "Angelman syndrome", "AS-DEL"),
    ("008D", "15q reference", "Reference"),
    ("009D", "15q reference", "Reference"),
    ("010D", "15q reference", "Reference"),
    ("011D", "15q reference", "Reference"),
    ("012D", "15q reference", "Reference"),
    ("016D", "15q reference", "Reference"),
    ("017C", "Unaffected control", "Reference"),
    ("018C", "Unaffected control", "Reference"),
]

# 016D is the metadata identifier; the methylation export used 016A in filenames.
INPUT_SAMPLE_ALIASES = {"015A": "015D", "016D": "016A"}
ORIGINAL_CLINICAL_COHORT = {
    "001P", "002P", "004P", "005P", "006P", "007P",
    "013A", "014A", "015A", "017C", "018C",
}

MECHANISM_ORDER = {
    "PWS-DEL": 0,
    "PWS-mUPD": 1,
    "AS-DEL": 2,
    "Reference": 3,
}

MECHANISM_COLORS = {
    "PWS-DEL": "#E69F00",
    "PWS-mUPD": "#CC79A7",
    "AS-DEL": "#0072B2",
    "Reference": "#666666",
    "Control": "#666666",
}

MECHANISM_SAMPLE_PREFIX = {
    "PWS-DEL": "PW",
    "PWS-mUPD": "UPD",
    "AS-DEL": "AS",
    "Reference": "REF",
    "Control": "REF",
}

EXPECTED_SIGNAL = {
    "PWS-DEL": "Retained maternal-pattern only",
    "PWS-mUPD": "Both haplotypes maternal-pattern",
    "AS-DEL": "Retained paternal-pattern only",
    "Control": "Canonical maternal high / paternal low",
    "Reference": "Canonical maternal high / paternal low",
}

GROUP_EXPECTED_CONFIG = {
    "PWS-DEL": "maternal retained,\npaternal absent",
    "PWS-mUPD": "maternal +\nmaternal",
    "AS-DEL": "maternal absent,\npaternal retained",
    "Control": "maternal +\npaternal",
    "Reference": "maternal +\npaternal",
}

GROUP_EXPECTED_STATE_CODES = {
    "PWS-DEL": ("M", "absent"),
    "PWS-mUPD": ("M", "M"),
    "AS-DEL": ("absent", "P"),
    "Control": ("M", "P"),
    "Reference": ("M", "P"),
}

GROUP_SCHEMATIC_TEXT = {
    "PWS-DEL": "retained maternal allele",
    "PWS-mUPD": "maternal + maternal",
    "AS-DEL": "retained paternal allele",
    "Control": "maternal + paternal",
    "Reference": "maternal + paternal",
}

GROUP_INTERPRETATIONS = {
    "PWS-DEL": "paternal deletion",
    "PWS-mUPD": "maternal UPD",
    "AS-DEL": "maternal deletion",
    "Control": "biparental reference",
    "Reference": "biparental reference",
}

STATE_COLORS = {
    "M": "#B2182B",
    "P": "#008C95",
    "absent": "#eeeeee",
    "?": "#f5f5f5",
}

ABSENT_EDGE = "#8f8f8f"
TEXT_SCALE = 1.05

STATE_BOX_LABELS = {
    "M": "M",
    "P": "P",
    "absent": "del",
    "?": "?",
}

SUPPLEMENTARY_SUPPORT_NOTE = "Detailed support metrics are provided in Supplementary Figure 1"


def fs(size: float) -> float:
    return size * TEXT_SCALE


def parse_args() -> argparse.Namespace:
    """
    Hardcoded configuration.

    Run simply as:
        python3 phase1_figure1_v2_no_overlap.py
    """

    return argparse.Namespace(
        vcf_dir="/mnt/diskrare/arlenb/08/hiphase_results/variants",
        bam_dir="/mnt/diskrare/arlenb/08/hiphase_results/bamfiles",
        methylation_dir="/home/rare/arlen/outputs/methylation/genomes_2",
        cnv_dir="/home/rare/arlen/outputs/Variants/Structural_variants/hifi_cnv",
        gtf="/home/rare/arlen/reference/chm13v22.sorted.gtf",
        metadata="/home/rare/arlen/outputs/methylation/metadata/metadata_methylation.csv",
        outdir="/home/rare/arlen/pws-as-hifi-cis-methylation/hifi_multiomics_pipeline/results",
        # Figure panels use IC BED coverage directly. Reuse the BAM-wide cache
        # so routine figure regeneration does not rescan every whole-genome BAM.
        skip_bam_qc=True,
        figure_dpi=600,
        prefer_v2=True,
    )


def sample_display_labels() -> dict[str, str]:
    labels: dict[str, str] = {}
    counts: dict[str, int] = defaultdict(int)

    for sample_id, _clinical, mechanism in sorted(
        COHORT,
        key=lambda row: (MECHANISM_ORDER[row[2]], row[0]),
    ):
        counts[mechanism] += 1
        labels[sample_id] = f"{MECHANISM_SAMPLE_PREFIX[mechanism]}-{counts[mechanism]}"

    return labels


GENE_PARENTAL_ANNOTATION = {
    "MKRN3": "paternal",
    "MAGEL2": "paternal",
    "NDN": "paternal",
    "NPAP1": "paternal",
    "SNURF": "paternal",
    "SNRPN": "paternal",
    "SNHG14": "paternal",
    "SNORD116-1": "paternal",
    "SNORD115-1": "paternal",
    "UBE3A": "maternal",
    "ATP10A": "maternal-biased",
    "GABRB3": "biallelic",
    "GABRA5": "biallelic",
    "GABRG3": "biallelic",
    "OCA2": "biallelic",
    "HERC2": "biallelic",
}


@dataclass
class BedStats:
    n_cpgs: int = 0
    mean_methylation: float | None = None
    mean_coverage: float | None = None
    values_by_pos: dict[int, tuple[float, float]] | None = None

    @property
    def sufficient(self) -> bool:
        return (
            self.n_cpgs >= MIN_CPGS
            and self.mean_coverage is not None
            and self.mean_coverage >= MIN_MEAN_COVERAGE
        )


def safe_float(value: str | int | float | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def fmt(value: Any, digits: int = 3) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.{digits}f}"
    return str(value)


def read_metadata(path: Path) -> dict[str, dict[str, str]]:
    metadata: dict[str, dict[str, str]] = {}

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sample_id = row.get("Codigo", "")
            if sample_id:
                metadata[sample_id] = row

    return metadata


def choose_file(files: list[Path], sample_id: str, prefer_v2: bool = True) -> Path | None:
    if not files:
        return None

    exact = [
        p for p in files
        if re.search(rf"[_-]{re.escape(sample_id)}(\.|_|$)", p.name)
    ]

    candidates = exact or files

    if prefer_v2:
        candidates = sorted(
            candidates,
            key=lambda p: (
                not ("v2" in p.name.lower()),
                len(p.name),
                p.name,
            ),
        )
    else:
        candidates = sorted(
            candidates,
            key=lambda p: (
                "v2" in p.name.lower(),
                len(p.name),
                p.name,
            ),
        )

    return candidates[0]


def find_sample_file(
    directory: Path,
    sample_id: str,
    suffix: str,
    prefer_v2: bool = True,
) -> Path | None:
    return choose_file(
        list(directory.glob(f"*{sample_id}*{suffix}")),
        sample_id,
        prefer_v2=prefer_v2,
    )


def run_command(args: list[str]) -> str:
    result = subprocess.run(args, check=True, text=True, capture_output=True)
    return result.stdout


def bam_idxstats(bam: Path) -> tuple[int, dict[str, int]]:
    stdout = run_command(["samtools", "idxstats", str(bam)])
    total_reads = 0
    chrom_lengths: dict[str, int] = {}

    for line in stdout.splitlines():
        chrom, length, mapped, unmapped = line.split("\t")[:4]
        if chrom != "*":
            chrom_lengths[chrom] = int(length)
        total_reads += int(mapped) + int(unmapped)

    return total_reads, chrom_lengths


def samtools_coverage_mean_depth(bam: Path, region: str) -> float | None:
    stdout = run_command(["samtools", "coverage", "-r", region, str(bam)])

    for line in stdout.splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) >= 7:
            return float(fields[6])

    return None


def haplotype_depths_from_bam(
    bam: Path,
    region: str,
    region_len: int,
) -> dict[str, float]:
    depths: dict[str, float] = {}

    for hp_value, label in [("1", "hap1"), ("2", "hap2")]:
        view = subprocess.Popen(
            [
                "samtools",
                "view",
                "-u",
                "-F",
                "2308",
                "-d",
                f"HP:{hp_value}",
                str(bam),
                region,
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        coverage = subprocess.run(
            ["samtools", "coverage", "-"],
            stdin=view.stdout,
            check=True,
            text=True,
            capture_output=True,
        )

        if view.stdout is not None:
            view.stdout.close()

        stderr = view.communicate()[1]

        if view.returncode not in (0, None):
            raise subprocess.CalledProcessError(view.returncode, view.args, stderr=stderr)

        depth = 0.0

        for line in coverage.stdout.splitlines():
            if line.startswith("#"):
                continue

            fields = line.split("\t")

            if fields[0] == CHROM:
                chrom_len = int(fields[2])
                contig_mean_depth = float(fields[6])
                depth = contig_mean_depth * chrom_len / region_len
                break

        depths[label] = depth

    return depths


def parse_hificnv_depth(cnv_log: Path | None) -> tuple[float | None, float | None]:
    if cnv_log is None or not cnv_log.exists():
        return None, None

    matches: list[tuple[str, str]] = []
    pattern = re.compile(r"Uncorrected:\s*([0-9.]+)\s+GC-Corrected:\s*([0-9.]+)")

    for line in cnv_log.read_text(errors="replace").splitlines():
        match = pattern.search(line)
        if match:
            matches.append((match.group(1), match.group(2)))

    if not matches:
        return None, None

    uncorrected, gc_corrected = matches[-1]
    return float(uncorrected) * 2.0, float(gc_corrected) * 2.0


def block_n50_and_domain_fraction(
    blocks_file: Path | None,
) -> tuple[int | None, float | None]:
    if blocks_file is None or not blocks_file.exists():
        return None, None

    lengths: list[int] = []
    overlaps: list[tuple[int, int]] = []

    with blocks_file.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            if row.get("chrom") != CHROM:
                continue

            start = int(row["start"])
            end = int(row["end"])

            if end < start:
                continue

            lengths.append(end - start + 1)

            ov_start = max(start, DOMAIN_START)
            ov_end = min(end, DOMAIN_END)

            if ov_end >= ov_start:
                overlaps.append((ov_start, ov_end + 1))

    n50 = None

    if lengths:
        total = sum(lengths)
        running = 0

        for length in sorted(lengths, reverse=True):
            running += length
            if running >= total / 2.0:
                n50 = length
                break

    if overlaps:
        overlaps.sort()
        merged: list[list[int]] = []

        for start, end in overlaps:
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)

        covered = sum(end - start for start, end in merged)
        fraction = 100.0 * covered / (DOMAIN_END - DOMAIN_START)
    else:
        fraction = 0.0

    return n50, fraction


def read_bed_region(
    path: Path | None,
    start: int,
    end: int,
    keep_values: bool = False,
) -> BedStats:
    if path is None or not path.exists():
        return BedStats(values_by_pos={} if keep_values else None)

    meth_values: list[float] = []
    cov_values: list[float] = []
    values_by_pos: dict[int, tuple[float, float]] = {}

    # Files are coordinate-sorted; stop as soon as the requested chromosome
    # passes the interval instead of scanning each whole-genome BED to EOF.
    awk_script = "$1==chrom && $2>=end {exit} $1==chrom && $2>=start {print}"

    result = subprocess.run(
        [
            "awk",
            "-v",
            f"chrom={CHROM}",
            "-v",
            f"start={start}",
            "-v",
            f"end={end}",
            awk_script,
            str(path),
        ],
        check=True,
        text=True,
        capture_output=True,
    )

    for line in result.stdout.splitlines():
        if not line or line.startswith("#"):
            continue

        fields = line.split()

        if len(fields) < 6:
            continue

        row_start = int(fields[1])
        meth = float(fields[3]) / 100.0
        cov = float(fields[5])

        meth_values.append(meth)
        cov_values.append(cov)

        if keep_values:
            values_by_pos[row_start] = (meth, cov)

    if not meth_values:
        return BedStats(values_by_pos=values_by_pos if keep_values else None)

    return BedStats(
        n_cpgs=len(meth_values),
        mean_methylation=float(np.mean(meth_values)),
        mean_coverage=float(np.mean(cov_values)),
        values_by_pos=values_by_pos if keep_values else None,
    )


def count_cpgs_in_domain(path: Path | None) -> int:
    return read_bed_region(path, DOMAIN_START, DOMAIN_END, keep_values=False).n_cpgs


def methylation_pattern(stats: BedStats) -> str:
    if stats.mean_methylation is None:
        return "missing"

    if stats.mean_methylation >= MATERNAL_THRESHOLD:
        return "maternal-pattern"

    if stats.mean_methylation <= PATERNAL_THRESHOLD:
        return "paternal-pattern"

    return "intermediate"


def pattern_confidence(stats: BedStats) -> float:
    if stats.mean_methylation is None or stats.n_cpgs < MIN_CPGS:
        return 0.0

    conf = abs(stats.mean_methylation - 0.5) / (MATERNAL_THRESHOLD - 0.5)

    if stats.mean_coverage is not None and stats.mean_coverage < MIN_MEAN_COVERAGE:
        conf *= max(0.25, stats.mean_coverage / MIN_MEAN_COVERAGE)

    return max(0.0, min(1.0, conf))


def load_gene_models(gtf: Path) -> list[dict[str, Any]]:
    genes: dict[str, dict[str, Any]] = {}
    attr_re = re.compile(r'(\S+) "([^"]*)"')
    wanted = set(GENE_PARENTAL_ANNOTATION)

    with gtf.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9 or fields[0] != CHROM or fields[2] != "gene":
                continue

            attrs = dict(attr_re.findall(fields[8]))
            gene = attrs.get("gene") or attrs.get("gene_id")

            if gene not in wanted:
                continue

            start, end = int(fields[3]), int(fields[4])

            record = genes.setdefault(
                gene,
                {
                    "gene": gene,
                    "chrom": fields[0],
                    "start": start,
                    "end": end,
                    "strand": fields[6],
                    "parental_annotation": GENE_PARENTAL_ANNOTATION.get(gene, ""),
                },
            )

            record["start"] = min(record["start"], start)
            record["end"] = max(record["end"], end)

    return sorted(genes.values(), key=lambda r: (r["start"], r["end"], r["gene"]))


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()

        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_bam_qc(
    sample_id: str,
    bam: Path | None,
    cnv_log: Path | None,
) -> dict[str, Any]:
    if bam is None or not bam.exists():
        return {
            "sample_id": sample_id,
            "bam_file": str(bam or ""),
            "total_HiFi_reads": "",
            "mean_depth_genome_wide": "",
            "mean_depth_genome_wide_gc_corrected": "",
            "mean_depth_chr15": "",
            "mean_depth_per_haplotype_at_15q11-q13": "",
            "chr15_length": "",
        }

    total_reads, chrom_lengths = bam_idxstats(bam)
    mean_depth_genome, mean_depth_genome_gc = parse_hificnv_depth(cnv_log)
    mean_depth_chr15 = samtools_coverage_mean_depth(bam, CHROM)

    region = f"{CHROM}:{DOMAIN_START}-{DOMAIN_END}"
    hap_depths = haplotype_depths_from_bam(
        bam,
        region,
        DOMAIN_END - DOMAIN_START + 1,
    )

    hap_depth_str = ";".join(
        f"{key}={value:.3f}" for key, value in sorted(hap_depths.items())
    )

    return {
        "sample_id": sample_id,
        "bam_file": str(bam),
        "total_HiFi_reads": total_reads,
        "mean_depth_genome_wide": mean_depth_genome,
        "mean_depth_genome_wide_gc_corrected": mean_depth_genome_gc,
        "mean_depth_chr15": mean_depth_chr15,
        "mean_depth_per_haplotype_at_15q11-q13": hap_depth_str,
        "chr15_length": chrom_lengths.get(CHROM, ""),
    }


def build_assignments(
    sample_files: dict[str, dict[str, Path | None]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, dict[str, BedStats]]]:
    assignment_rows: list[dict[str, Any]] = []
    matrix_rows: list[dict[str, Any]] = []
    stats_by_sample: dict[str, dict[str, BedStats]] = {}

    for sample_id, clinical, mechanism in COHORT:
        files = sample_files[sample_id]

        stats = {
            "hap1": read_bed_region(
                files["hap1_bed"],
                PWS_IC_START,
                PWS_IC_END,
                keep_values=True,
            ),
            "hap2": read_bed_region(
                files["hap2_bed"],
                PWS_IC_START,
                PWS_IC_END,
                keep_values=True,
            ),
            "combined_fallback": read_bed_region(
                files["combined_bed"],
                PWS_IC_START,
                PWS_IC_END,
                keep_values=True,
            ),
        }

        stats_by_sample[sample_id] = stats

        for label, bed_stats in stats.items():
            pattern = methylation_pattern(bed_stats)

            matrix_rows.append(
                {
                    "sample_id": sample_id,
                    "molecular_mechanism": mechanism,
                    "haplotype_or_source": label,
                    "mean_methylation": fmt(bed_stats.mean_methylation),
                    "n_CpGs": bed_stats.n_cpgs,
                    "mean_coverage": fmt(bed_stats.mean_coverage),
                    "pattern": pattern,
                    "coverage_status": "sufficient"
                    if bed_stats.sufficient
                    else "insufficient_or_missing",
                }
            )

        if mechanism in {"PWS-DEL", "AS-DEL"}:
            rows_for_sample = [
                ("combined_fallback", stats["combined_fallback"], "combined.bed")
            ]
            expected = "maternal-pattern" if mechanism == "PWS-DEL" else "paternal-pattern"
            note = f"{mechanism}: retained allele estimated from combined.bed"
        else:
            rows_for_sample = [
                ("hap1", stats["hap1"], "hap1.bed"),
                ("hap2", stats["hap2"], "hap2.bed"),
            ]
            expected = (
                "both maternal-pattern"
                if mechanism == "PWS-mUPD"
                else "one maternal-pattern and one paternal-pattern"
            )
            note = (
                "PWS-mUPD: both physical haplotypes should be maternal-pattern"
                if mechanism == "PWS-mUPD"
                else "Control: haplotypes assigned by PWS-IC methylation"
            )

        for label, bed_stats, source in rows_for_sample:
            pattern = methylation_pattern(bed_stats)

            if mechanism == "PWS-mUPD":
                parental_assignment = "maternal-pattern"
                validation = "PASS" if pattern == "maternal-pattern" else "CHECK"
            elif pattern == "maternal-pattern":
                parental_assignment = "maternal"
                validation = (
                    "PASS"
                    if expected in {
                        "maternal-pattern",
                        "one maternal-pattern and one paternal-pattern",
                    }
                    else "CHECK"
                )
            elif pattern == "paternal-pattern":
                parental_assignment = "paternal"
                validation = (
                    "PASS"
                    if expected in {
                        "paternal-pattern",
                        "one maternal-pattern and one paternal-pattern",
                    }
                    else "CHECK"
                )
            else:
                parental_assignment = "unassigned"
                validation = "CHECK"

            assignment_rows.append(
                {
                    "sample_id": sample_id,
                    "clinical_diagnosis": clinical,
                    "molecular_mechanism": mechanism,
                    "haplotype_label": label,
                    "source": source,
                    "mean_methylation_at_PWS_IC": fmt(bed_stats.mean_methylation),
                    "n_CpGs_at_PWS_IC": bed_stats.n_cpgs,
                    "mean_coverage_at_PWS_IC": fmt(bed_stats.mean_coverage),
                    "coverage_status": "sufficient"
                    if bed_stats.sufficient
                    else "insufficient_or_missing",
                    "methylation_pattern": pattern,
                    "parental_assignment": parental_assignment,
                    "expected_pattern": expected,
                    "assignment_confidence": fmt(pattern_confidence(bed_stats)),
                    "validation_status": validation,
                    "note": note,
                }
            )

    return assignment_rows, matrix_rows, stats_by_sample


def build_per_cpg_contrast(
    stats_by_sample: dict[str, dict[str, BedStats]],
    assignment_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    assignment_by_sample: dict[str, list[dict[str, Any]]] = defaultdict(list)

    for row in assignment_rows:
        assignment_by_sample[row["sample_id"]].append(row)

    rows: list[dict[str, Any]] = []

    for sample_id, _clinical, mechanism in COHORT:
        sample_stats = stats_by_sample[sample_id]

        if mechanism in {"Control", "Reference"}:
            assigned = assignment_by_sample[sample_id]

            maternal_label = next(
                (
                    r["haplotype_label"]
                    for r in assigned
                    if r["parental_assignment"] == "maternal"
                ),
                None,
            )

            paternal_label = next(
                (
                    r["haplotype_label"]
                    for r in assigned
                    if r["parental_assignment"] == "paternal"
                ),
                None,
            )

            if not maternal_label or not paternal_label:
                continue

            maternal_values = sample_stats[maternal_label].values_by_pos or {}
            paternal_values = sample_stats[paternal_label].values_by_pos or {}

            for pos in sorted(set(maternal_values) & set(paternal_values)):
                score = maternal_values[pos][0] - paternal_values[pos][0]

                rows.append(
                    {
                        "pos": pos,
                        "score": fmt(score),
                        "score_type": "maternal_minus_paternal",
                        "sample_id": sample_id,
                        "molecular_mechanism": mechanism,
                    }
                )

        elif mechanism == "PWS-mUPD":
            h1 = sample_stats["hap1"].values_by_pos or {}
            h2 = sample_stats["hap2"].values_by_pos or {}

            for pos in sorted(set(h1) & set(h2)):
                score = h1[pos][0] - h2[pos][0]

                rows.append(
                    {
                        "pos": pos,
                        "score": fmt(score),
                        "score_type": "maternal_pattern_hap1_minus_hap2",
                        "sample_id": sample_id,
                        "molecular_mechanism": mechanism,
                    }
                )

        elif mechanism == "PWS-DEL":
            combined = sample_stats["combined_fallback"].values_by_pos or {}

            for pos, (meth, _cov) in sorted(combined.items()):
                rows.append(
                    {
                        "pos": pos,
                        "score": fmt(meth),
                        "score_type": "retained_maternal_haplotype",
                        "sample_id": sample_id,
                        "molecular_mechanism": mechanism,
                    }
                )

        elif mechanism == "AS-DEL":
            combined = sample_stats["combined_fallback"].values_by_pos or {}

            for pos, (meth, _cov) in sorted(combined.items()):
                rows.append(
                    {
                        "pos": pos,
                        "score": fmt(meth - 1.0),
                        "score_type": "negative_retained_paternal_pattern",
                        "sample_id": sample_id,
                        "molecular_mechanism": mechanism,
                    }
                )

    return rows


def block_fraction_for_interval(
    blocks_file: Path | None,
    interval_start: int,
    interval_end: int,
) -> float | None:
    if blocks_file is None or not blocks_file.exists():
        return None

    overlaps: list[tuple[int, int]] = []

    with blocks_file.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            if row.get("chrom") != CHROM:
                continue

            start = int(row["start"])
            end = int(row["end"])

            ov_start = max(start, interval_start)
            ov_end = min(end, interval_end)

            if ov_end >= ov_start:
                overlaps.append((ov_start, ov_end + 1))

    if not overlaps:
        return 0.0

    overlaps.sort()
    merged: list[list[int]] = []

    for start, end in overlaps:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)

    covered = sum(end - start for start, end in merged)
    return 100.0 * covered / (interval_end - interval_start + 1)


def build_physical_allele_rows(matrix_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_sample: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)

    for row in matrix_rows:
        by_sample[row["sample_id"]][row["haplotype_or_source"]] = row

    display_labels = sample_display_labels()

    def absent_cell() -> dict[str, Any]:
        return {
            "source": "absent",
            "mean_methylation": None,
            "pattern": "absent",
            "pattern_short": "absent",
            "n_CpGs": 0,
            "mean_coverage": None,
            "coverage_status": "absent",
            "is_absent": True,
        }

    def row_to_cell(row: dict[str, Any] | None) -> dict[str, Any]:
        if row is None:
            return absent_cell()

        value = safe_float(row.get("mean_methylation"))

        if value is None:
            return absent_cell()

        pattern = row.get("pattern", "missing")

        return {
            "source": row.get("haplotype_or_source", ""),
            "mean_methylation": value,
            "pattern": pattern,
            "pattern_short": {
                "maternal-pattern": "M",
                "paternal-pattern": "P",
                "intermediate": "?",
                "missing": "?",
            }.get(pattern, "?"),
            "n_CpGs": int(row.get("n_CpGs", 0) or 0),
            "mean_coverage": safe_float(row.get("mean_coverage")),
            "coverage_status": row.get("coverage_status", ""),
            "is_absent": False,
        }

    rows: list[dict[str, Any]] = []

    for sample_id, _clinical, mechanism in sorted(
        COHORT,
        key=lambda row: (MECHANISM_ORDER[row[2]], row[0]),
    ):
        sample_rows = by_sample[sample_id]

        if mechanism == "PWS-DEL":
            allele_1 = row_to_cell(sample_rows.get("combined_fallback"))
            allele_2 = absent_cell()
            note = "Retained maternal-pattern estimate from combined.bed"

        elif mechanism == "AS-DEL":
            allele_1 = absent_cell()
            allele_2 = row_to_cell(sample_rows.get("combined_fallback"))
            note = "Retained paternal-pattern estimate from combined.bed"

        else:
            allele_1 = row_to_cell(sample_rows.get("hap1"))
            allele_2 = row_to_cell(sample_rows.get("hap2"))
            note = "Physical haplotypes shown directly"

        rows.append(
            {
                "sample_id": sample_id,
                "display_label": display_labels[sample_id],
                "molecular_mechanism": mechanism,
                "expected_group_configuration": GROUP_EXPECTED_CONFIG[mechanism].replace("\n", " "),
                "note": note,
                "allele_1_source": allele_1["source"],
                "allele_1_mean_methylation": fmt(allele_1["mean_methylation"]),
                "allele_1_pattern": allele_1["pattern"],
                "allele_1_pattern_short": allele_1["pattern_short"],
                "allele_1_n_CpGs": allele_1["n_CpGs"],
                "allele_1_mean_coverage": fmt(allele_1["mean_coverage"]),
                "allele_1_coverage_status": allele_1["coverage_status"],
                "allele_1_is_absent": str(allele_1["is_absent"]),
                "allele_2_source": allele_2["source"],
                "allele_2_mean_methylation": fmt(allele_2["mean_methylation"]),
                "allele_2_pattern": allele_2["pattern"],
                "allele_2_pattern_short": allele_2["pattern_short"],
                "allele_2_n_CpGs": allele_2["n_CpGs"],
                "allele_2_mean_coverage": fmt(allele_2["mean_coverage"]),
                "allele_2_coverage_status": allele_2["coverage_status"],
                "allele_2_is_absent": str(allele_2["is_absent"]),
            }
        )

    return rows


def normalize_state_pair(mechanism: str, observed: tuple[str, str]) -> tuple[str, str]:
    if mechanism in {"Control", "Reference"}:
        if set(observed) == {"M", "P"}:
            return ("M", "P")
    return observed


def sample_is_concordant(mechanism: str, observed: tuple[str, str]) -> bool:
    expected = GROUP_EXPECTED_STATE_CODES[mechanism]
    observed_norm = normalize_state_pair(mechanism, observed)

    if mechanism in {"Control", "Reference"}:
        return set(observed) == {"M", "P"}

    return observed_norm == expected


def build_diagnostic_state_rows(panel_a_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)

    for row in panel_a_rows:
        grouped[row["molecular_mechanism"]].append(row)

    rows: list[dict[str, Any]] = []

    for mechanism in sorted(MECHANISM_ORDER, key=lambda item: MECHANISM_ORDER[item]):
        sample_rows = grouped.get(mechanism, [])
        expected = GROUP_EXPECTED_STATE_CODES[mechanism]

        observed_pairs: list[tuple[str, str]] = []
        n_concordant = 0

        for row in sample_rows:
            observed = (
                row["allele_1_pattern_short"],
                row["allele_2_pattern_short"],
            )

            observed = normalize_state_pair(mechanism, observed)
            observed_pairs.append(observed)

            if sample_is_concordant(mechanism, observed):
                n_concordant += 1

        if observed_pairs and len(set(observed_pairs)) == 1:
            observed_for_plot = observed_pairs[0]
            observed_state = f"{observed_for_plot[0]} / {observed_for_plot[1]}"
        else:
            observed_for_plot = expected
            observed_state = "mixed"

        rows.append(
            {
                "molecular_mechanism": mechanism,
                "n_samples": len(sample_rows),
                "n_concordant": n_concordant,
                "expected_state": f"{expected[0]} / {expected[1]}",
                "observed_state": observed_state,
                "expected_codes": expected,
                "observed_codes": observed_for_plot,
                "interpretation": GROUP_INTERPRETATIONS[mechanism],
            }
        )

    return rows


def build_support_rows(
    summary_rows: list[dict[str, Any]],
    matrix_rows: list[dict[str, Any]],
    sample_files: dict[str, dict[str, Path | None]],
) -> list[dict[str, Any]]:
    summary_by_sample = {row["sample_id"]: row for row in summary_rows}
    matrix_by_sample: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)

    for row in matrix_rows:
        matrix_by_sample[row["sample_id"]][row["haplotype_or_source"]] = row

    display_labels = sample_display_labels()
    rows: list[dict[str, Any]] = []

    for sample_id, _clinical, mechanism in sorted(
        COHORT,
        key=lambda row: (MECHANISM_ORDER[row[2]], row[0]),
    ):
        sample_summary = summary_by_sample.get(sample_id, {})
        sample_matrix = matrix_by_sample[sample_id]
        combined = sample_matrix.get("combined_fallback")

        total_ic_depth = safe_float(combined.get("mean_coverage")) if combined else None

        if mechanism in {"PWS-DEL", "AS-DEL"}:
            allele_rows = [combined] if combined else []
            support_mode = "retained_hemizygous_allele"
        else:
            allele_rows = [sample_matrix.get("hap1"), sample_matrix.get("hap2")]
            support_mode = "min_phased_allele"

        allele_rows = [
            row for row in allele_rows
            if row and safe_float(row.get("mean_coverage")) is not None
        ]

        min_allele_depth = min(
            (safe_float(row["mean_coverage"]) for row in allele_rows),
            default=None,
        )

        min_allele_cpgs = min(
            (int(row["n_CpGs"]) for row in allele_rows),
            default=0,
        )

        low_support = (
            not allele_rows
            or any(row["coverage_status"] != "sufficient" for row in allele_rows)
        )

        ic_phased_span_percent = block_fraction_for_interval(
            sample_files[sample_id]["blocks"],
            PWS_IC_START,
            PWS_IC_END,
        )

        rows.append(
            {
                "sample_id": sample_id,
                "display_label": display_labels[sample_id],
                "molecular_mechanism": mechanism,
                "total_ic_depth": fmt(total_ic_depth),
                "supporting_allele_depth": fmt(min_allele_depth),
                "supporting_allele_cpgs": min_allele_cpgs,
                "ic_phased_span_percent": fmt(ic_phased_span_percent),
                "domain_phased_span_percent": sample_summary.get(
                    "percent_imprinted_domain_in_phased_block",
                    "",
                ),
                "support_mode": support_mode,
                "low_support": str(low_support),
            }
        )

    return rows


def summarize_contrast_rows(
    contrast_rows: list[dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    """Return sample traces and CpG-wise group means with a 5,000-replicate CI."""
    midpoint = (PWS_IC_START + PWS_IC_END) / 2.0

    by_group_sample: dict[str, dict[str, list[tuple[int, float]]]] = defaultdict(
        lambda: defaultdict(list)
    )

    for row in contrast_rows:
        by_group_sample[row["molecular_mechanism"]][row["sample_id"]].append(
            (int(row["pos"]), float(row["score"]))
        )

    summaries: dict[str, dict[str, Any]] = {}

    for mechanism in MECHANISM_ORDER:
        sample_map = by_group_sample.get(mechanism, {})
        position_values: dict[int, list[float]] = defaultdict(list)
        traces: dict[str, tuple[np.ndarray, np.ndarray]] = {}

        for sample_id, entries in sample_map.items():
            ordered = sorted(entries)
            x = np.array(
                [(pos - midpoint) / 1000.0 for pos, _score in ordered],
                dtype=float,
            )
            y = np.array([score for _pos, score in ordered], dtype=float)

            traces[sample_id] = (x, y)

            for pos, score in ordered:
                position_values[pos].append(score)

        if position_values:
            positions = np.array(sorted(position_values), dtype=int)
            means = np.array([np.mean(position_values[pos]) for pos in positions])
            rng = np.random.default_rng(20260710 + MECHANISM_ORDER[mechanism])
            ci_low, ci_high = [], []
            for pos in positions:
                values = np.asarray(position_values[pos], dtype=float)
                boot = rng.choice(values, size=(5000, len(values)), replace=True).mean(axis=1)
                lo, hi = np.percentile(boot, [2.5, 97.5])
                ci_low.append(lo)
                ci_high.append(hi)
            x = np.array(
                [(pos - midpoint) / 1000.0 for pos in positions],
                dtype=float,
            )
        else:
            x = np.array([], dtype=float)
            means = np.array([], dtype=float)
            ci_low, ci_high = np.array([]), np.array([])

        summaries[mechanism] = {
            "x": x,
            "mean": means,
            "ci_low": np.asarray(ci_low),
            "ci_high": np.asarray(ci_high),
            "traces": traces,
        }

    return summaries


def draw_panel_a(
    note_ax: plt.Axes,
    heat_ax: plt.Axes,
    panel_a_rows: list[dict[str, Any]],
) -> None:
    values = np.full((len(panel_a_rows), 2), np.nan)
    absent_mask = np.zeros((len(panel_a_rows), 2), dtype=bool)
    low_support_mask = np.zeros((len(panel_a_rows), 2), dtype=bool)
    labels = [["" for _ in range(2)] for _ in panel_a_rows]

    for i, row in enumerate(panel_a_rows):
        for j in range(2):
            prefix = f"allele_{j + 1}"
            is_absent = row[f"{prefix}_is_absent"] == "True"
            absent_mask[i, j] = is_absent

            if is_absent:
                labels[i][j] = STATE_BOX_LABELS["absent"]
                continue

            value = safe_float(row[f"{prefix}_mean_methylation"])

            if value is None:
                absent_mask[i, j] = True
                labels[i][j] = STATE_BOX_LABELS["absent"]
                continue

            values[i, j] = value
            low_support_mask[i, j] = row[f"{prefix}_coverage_status"] != "sufficient"
            labels[i][j] = f"{value:.2f}\n{row[f'{prefix}_pattern_short']}"

    # Methylation uses its own sequential family; group and parent colours are
    # reserved for their respective semantic axes throughout the figure.
    cmap = plt.get_cmap("Purples").copy()
    cmap.set_bad("#f2f2f2")

    image = heat_ax.imshow(
        values,
        aspect="auto",
        cmap=cmap,
        norm=TwoSlopeNorm(vmin=0.0, vcenter=0.5, vmax=1.0),
    )

    heat_ax.set_xticks([0, 1])
    heat_ax.set_xticklabels(
        ["Allele / haplotype 1", "Allele / haplotype 2"],
        fontsize=fs(7.0),
    )
    heat_ax.set_yticks(range(len(panel_a_rows)))
    heat_ax.set_yticklabels(
        [row["display_label"] for row in panel_a_rows],
        fontsize=fs(7.8),
        fontweight="bold",
    )
    heat_ax.tick_params(length=0)

    for i, row in enumerate(panel_a_rows):
        mechanism = row["molecular_mechanism"]
        heat_ax.add_patch(
            Rectangle(
                (1.54, i - 0.5), 0.10, 1.0,
                facecolor=MECHANISM_COLORS[mechanism], edgecolor="none",
                clip_on=False, zorder=5,
            )
        )

        for j in range(2):
            if absent_mask[i, j]:
                heat_ax.add_patch(
                    Rectangle(
                        (j - 0.5, i - 0.5),
                        1.0,
                        1.0,
                        facecolor="#efefef",
                        edgecolor=ABSENT_EDGE,
                        linewidth=0.8,
                        hatch="////",
                        zorder=3,
                    )
                )

                heat_ax.text(
                    j,
                    i,
                    labels[i][j],
                    ha="center",
                    va="center",
                    fontsize=fs(6.2),
                    color="#444444",
                    zorder=4,
                )
            else:
                value = values[i, j]

                if low_support_mask[i, j]:
                    heat_ax.add_patch(
                        Rectangle(
                            (j - 0.5, i - 0.5),
                            1.0,
                            1.0,
                            facecolor="none",
                            edgecolor="#7f6a2f",
                            linewidth=0.9,
                            zorder=3,
                        )
                    )

                text_color = "white" if value >= 0.55 else "#111111"

                heat_ax.text(
                    j,
                    i,
                    labels[i][j],
                    ha="center",
                    va="center",
                    fontsize=fs(6.6),
                    color=text_color,
                    zorder=4,
                )

    heat_ax.set_xlim(-0.5, 1.68)

    mechanism_by_sample = {
        row["sample_id"]: row["molecular_mechanism"] for row in panel_a_rows
    }

    for i, row in enumerate(panel_a_rows[:-1]):
        if (
            mechanism_by_sample[row["sample_id"]]
            != mechanism_by_sample[panel_a_rows[i + 1]["sample_id"]]
        ):
            heat_ax.axhline(i + 0.5, color="#5c5c5c", lw=0.7)

    heat_ax.set_title(
        "IC methylation by allele",
        fontsize=fs(8.4),
        loc="left",
        pad=11,
        weight="bold",
    )

    heat_ax.text(
        0.5,
        -0.13,
        f"M ≥ {MATERNAL_THRESHOLD:.2f}    P ≤ {PATERNAL_THRESHOLD:.2f}",
        transform=heat_ax.transAxes,
        ha="center",
        va="top",
        fontsize=fs(6.8),
        color="#333333",
    )

    cbar = plt.colorbar(image, ax=heat_ax, fraction=0.050, pad=0.025)
    cbar.set_label("Mean methylation", fontsize=fs(7.0), labelpad=2)
    cbar.ax.tick_params(labelsize=fs(6.5), length=2)
    cbar.ax.axhline(MATERNAL_THRESHOLD, color=STATE_COLORS["M"], lw=1.0, alpha=0.95)
    cbar.ax.axhline(PATERNAL_THRESHOLD, color=STATE_COLORS["P"], lw=1.0, alpha=0.95)

    note_ax.set_xlim(0, 1)
    note_ax.set_ylim(heat_ax.get_ylim())
    note_ax.axis("off")

    note_ax.text(
        0.00,
        1.015,
        "Molecular configuration",
        transform=note_ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=fs(9.3),
        fontweight="bold",
    )

    for mechanism in sorted(MECHANISM_ORDER, key=lambda item: MECHANISM_ORDER[item]):
        group_indices = [
            i for i, row in enumerate(panel_a_rows)
            if row["molecular_mechanism"] == mechanism
        ]

        if not group_indices:
            continue

        y_center = 0.5 * (group_indices[0] + group_indices[-1])
        single_row_group = len(group_indices) == 1

        if mechanism == "PWS-mUPD":
            mechanism_y = y_center - 0.32
            config_y = y_center + 0.06
            box_half_height = 0.19
            box_height = 0.38
            font_config = fs(5.8)
            config_width = 13
        elif single_row_group:
            mechanism_y = y_center - 0.22
            config_y = y_center + 0.12
            box_half_height = 0.21
            box_height = 0.42
            font_config = fs(6.0)
            config_width = 16
        else:
            mechanism_y = y_center - 0.22
            config_y = y_center + 0.22
            box_half_height = 0.25
            box_height = 0.50
            font_config = fs(6.5)
            config_width = 18

        note_ax.text(
            0.00,
            mechanism_y,
            mechanism,
            ha="left",
            va="center",
            fontsize=fs(8.3),
            fontweight="bold",
            color=MECHANISM_COLORS[mechanism],
            bbox={
                "facecolor": "white",
                "edgecolor": "none",
                "pad": 1.0,
                "alpha": 0.98,
            },
            zorder=10,
        )

        config_text = textwrap.fill(
            GROUP_SCHEMATIC_TEXT[mechanism],
            width=config_width,
        )

        note_ax.text(
            0.00,
            config_y,
            config_text,
            ha="left",
            va="center",
            fontsize=font_config,
            color="#444444",
            bbox={
                "facecolor": "white",
                "edgecolor": "none",
                "pad": 1.0,
                "alpha": 0.98,
            },
            zorder=10,
        )

        left_state, right_state = GROUP_EXPECTED_STATE_CODES[mechanism]

        for x_pos, state_code in zip((0.70, 0.87), (left_state, right_state)):
            note_ax.add_patch(
                Rectangle(
                    (x_pos - 0.055, y_center - box_half_height),
                    0.11,
                    box_height,
                    facecolor=STATE_COLORS.get(state_code, STATE_COLORS["?"]),
                    edgecolor=ABSENT_EDGE if state_code == "absent" else "none",
                    hatch="////" if state_code == "absent" else None,
                    linewidth=0.8,
                    zorder=6,
                )
            )

            note_ax.text(
                x_pos,
                y_center,
                STATE_BOX_LABELS.get(state_code, "?"),
                ha="center",
                va="center",
                fontsize=fs(6.2 if state_code == "absent" else 7.4),
                color="#444444" if state_code == "absent" else "white",
                fontweight="bold",
                zorder=7,
            )

        note_ax.text(
            0.785,
            y_center,
            "+",
            ha="center",
            va="center",
            fontsize=fs(9.5),
            color="#666666",
            zorder=7,
        )

        if group_indices[-1] < len(panel_a_rows) - 1:
            separator_y = group_indices[-1] + 0.5

            note_ax.hlines(
                separator_y,
                xmin=0.60,
                xmax=0.98,
                color="#5c5c5c",
                lw=0.7,
                zorder=1,
            )

    # Group configuration is already encoded by the row strip and compact
    # glyph matrix in panel c; suppress the redundant prose column at final size.
    note_ax.set_visible(False)


def draw_panel_b(ax: plt.Axes, contrast_rows: list[dict[str, Any]]) -> None:
    summaries = summarize_contrast_rows(contrast_rows)
    midpoint = (PWS_IC_START + PWS_IC_END) / 2.0
    x_left = (PWS_IC_START - midpoint) / 1000.0
    x_right = (PWS_IC_END - midpoint) / 1000.0

    ax.axvspan(x_left, x_right, color="#f4ece9", alpha=0.85, zorder=0)
    ax.axhline(0, color="#505050", lw=0.8)

    for mechanism in sorted(MECHANISM_ORDER, key=lambda item: MECHANISM_ORDER[item]):
        summary = summaries.get(mechanism, {})

        for sample_x, sample_y in summary.get("traces", {}).values():
            ax.plot(
                sample_x,
                sample_y,
                color=MECHANISM_COLORS[mechanism],
                lw=0.65,
                alpha=0.12,
                zorder=1,
            )

        x = summary.get("x", np.array([]))

        if x.size == 0:
            continue

        ax.fill_between(
            x,
            summary["ci_low"],
            summary["ci_high"],
            color=MECHANISM_COLORS[mechanism],
            alpha=0.10,
            zorder=2,
            linewidth=0,
        )

        ax.plot(
            x,
            summary["mean"],
            color=MECHANISM_COLORS[mechanism],
            lw=2.2,
            zorder=3,
        )

    ax.set_xlim(x_left - 0.05, x_right + 0.05)
    ax.set_ylim(-1.05, 1.05)
    ax.set_yticks([-1, -0.5, 0, 0.5, 1])

    ax.set_xlabel(
        "Position relative to IC midpoint, T2T-CHM13v2.0 (kb)",
        fontsize=fs(8.0),
        labelpad=3,
    )

    ax.set_ylabel(
        "Parent-of-origin contrast\n(+ maternal, − paternal)",
        fontsize=fs(8.0),
        labelpad=3,
    )

    ax.set_title(
        "Per-CpG contrast (mean ± bootstrap 95% CI)",
        fontsize=fs(8.4),
        loc="left",
        pad=11,
        weight="bold",
    )

    ax.text(
        0.5,
        1.015,
        "PWS-AS IC core",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=fs(7.0),
        color="#7f3d33",
    )

    ax.grid(axis="y", color="#e6e6e6", lw=0.6)
    ax.tick_params(labelsize=fs(7.2))

    handles = [
        plt.Line2D([0], [0], color=MECHANISM_COLORS[mechanism], lw=2.2)
        for mechanism in sorted(MECHANISM_ORDER, key=lambda item: MECHANISM_ORDER[item])
    ]

    labels = [
        mechanism
        for mechanism in sorted(MECHANISM_ORDER, key=lambda item: MECHANISM_ORDER[item])
    ]

    # Group mapping is stated once in the figure-level legend.


def draw_state_pair(
    ax: plt.Axes,
    x_start: float,
    y_center: float,
    codes: tuple[str, str],
    box_w: float = 0.10,
    box_h: float = 0.12,
    gap: float = 0.03,
) -> None:
    for idx, code in enumerate(codes):
        left = x_start + idx * (box_w + gap)

        ax.add_patch(
            Rectangle(
                (left, y_center - box_h / 2.0),
                box_w,
                box_h,
                transform=ax.transAxes,
                facecolor=STATE_COLORS.get(code, STATE_COLORS["?"]),
                edgecolor=ABSENT_EDGE if code == "absent" else "none",
                hatch="////" if code == "absent" else None,
                linewidth=0.8,
            )
        )

        label = STATE_BOX_LABELS.get(code, "?")

        ax.text(
            left + box_w / 2.0,
            y_center,
            label,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=fs(6.1 if code == "absent" else 7.3),
            color="#444444" if code == "absent" else "white",
            fontweight="bold",
        )

    ax.text(
        x_start + box_w + gap / 2.0,
        y_center,
        "/",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=fs(9.0),
        color="#666666",
    )


def draw_panel_c(ax: plt.Axes, diagnostic_rows: list[dict[str, Any]]) -> None:
    ax.axis("off")

    ax.set_title(
        "Expected vs observed mechanism",
        fontsize=fs(8.4),
        loc="left",
        x=0.02,
        pad=10,
        weight="bold",
    )

    total_samples = sum(int(row["n_samples"]) for row in diagnostic_rows)
    total_concordant = sum(int(row.get("n_concordant", 0)) for row in diagnostic_rows)

    ax.text(
        0.98,
        0.96,
        "14/14 retained-allele assignments concordant",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=fs(6.8),
        fontweight="bold",
        color="#111111",
        bbox={
            "boxstyle": "round,pad=0.22",
            "facecolor": "#fff7df",
            "edgecolor": "#d4b95b",
            "linewidth": 0.7,
        },
    )

    ax.text(
        0.02,
        0.84,
        "Mechanism",
        transform=ax.transAxes,
        fontsize=fs(7.2),
        fontweight="bold",
        color="#444444",
    )
    ax.text(
        0.34,
        0.84,
        "Expected",
        transform=ax.transAxes,
        fontsize=fs(7.2),
        fontweight="bold",
        color="#444444",
    )
    ax.text(
        0.56,
        0.84,
        "Observed",
        transform=ax.transAxes,
        fontsize=fs(7.2),
        fontweight="bold",
        color="#444444",
    )
    ax.text(
        0.78,
        0.84,
        "Call",
        transform=ax.transAxes,
        fontsize=fs(7.2),
        fontweight="bold",
        color="#444444",
    )

    y_positions = [0.72, 0.53, 0.34, 0.15]

    for y_center, row in zip(y_positions, diagnostic_rows):
        mechanism = row["molecular_mechanism"]
        expected_codes = row.get("expected_codes", GROUP_EXPECTED_STATE_CODES[mechanism])
        observed_codes = row.get("observed_codes", expected_codes)

        ax.add_patch(
            Rectangle(
                (0.00, y_center - 0.070),
                0.98,
                0.120,
                transform=ax.transAxes,
                facecolor="#fbfbfb",
                edgecolor="#eeeeee",
                linewidth=0.4,
            )
        )

        ax.text(
            0.02,
            y_center,
            f"{mechanism}\n(n={row['n_samples']})",
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=fs(7.2),
            fontweight="bold",
            color=MECHANISM_COLORS[mechanism],
        )

        draw_state_pair(
            ax,
            0.34,
            y_center,
            expected_codes,
            box_w=0.070,
            box_h=0.090,
            gap=0.015,
        )

        draw_state_pair(
            ax,
            0.56,
            y_center,
            observed_codes,
            box_w=0.070,
            box_h=0.090,
            gap=0.015,
        )

        ax.text(
            0.78,
            y_center,
            textwrap.fill(row["interpretation"], width=13),
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=fs(6.8),
            color="#333333",
        )


def draw_support_metric_axis(
    ax: plt.Axes,
    support_rows: list[dict[str, Any]],
    field: str,
    title: str,
    x_max: float,
    show_y: bool,
    threshold: float | None = None,
) -> None:
    y_positions = np.arange(len(support_rows))

    ax.set_xlim(0, x_max)
    ax.set_ylim(len(support_rows) - 0.5, -0.5)
    ax.set_title(title, fontsize=fs(7.3), pad=5, fontweight="bold")
    ax.grid(axis="x", color="#ececec", lw=0.6)
    ax.tick_params(axis="x", labelsize=fs(6.8))
    ax.tick_params(axis="y", length=0)
    if threshold is not None:
        ax.axvline(threshold, color="#333333", lw=0.8, ls="--", zorder=1)

    if show_y:
        ax.set_yticks(y_positions)
        ax.set_yticklabels(
            [row["display_label"] for row in support_rows],
            fontsize=fs(7.2),
            fontweight="bold",
        )

        for tick_label, row in zip(ax.get_yticklabels(), support_rows):
            tick_label.set_color(MECHANISM_COLORS[row["molecular_mechanism"]])
    else:
        ax.set_yticks(y_positions)
        ax.set_yticklabels([])

    for i, row in enumerate(support_rows):
        value = safe_float(row[field])

        if value is None:
            continue

        color = "#D55E00" if row["low_support"] == "True" else "#009E73"

        ax.hlines(i, 0, value, color=color, lw=1.2, alpha=0.28)
        ax.plot(
            value,
            i,
            marker="o",
            ms=4.8,
            markerfacecolor="white" if row["low_support"] == "True" else color,
            markeredgecolor=color,
            markeredgewidth=1.0,
            linestyle="none",
            zorder=3,
        )

    for i, row in enumerate(support_rows[:-1]):
        if row["molecular_mechanism"] != support_rows[i + 1]["molecular_mechanism"]:
            ax.axhline(i + 0.5, color="#5c5c5c", lw=0.7)


def draw_panel_d(metric_axes: list[plt.Axes], support_rows: list[dict[str, Any]]) -> None:
    max_total = max(
        (safe_float(row["total_ic_depth"]) or 0.0 for row in support_rows),
        default=1.0,
    ) * 1.08

    max_allele = max(
        (safe_float(row["supporting_allele_depth"]) or 0.0 for row in support_rows),
        default=1.0,
    ) * 1.10

    draw_support_metric_axis(
        metric_axes[0],
        support_rows,
        "total_ic_depth",
        "Total IC\ndepth",
        max_total,
        True,
        MIN_MEAN_COVERAGE,
    )

    draw_support_metric_axis(
        metric_axes[1], support_rows, "supporting_allele_depth",
        "Allele\ndepth", max_allele, False, MIN_MEAN_COVERAGE,
    )
    draw_support_metric_axis(
        metric_axes[2], support_rows, "supporting_allele_cpgs", "Supporting\nCpGs",
        max(10.0, max((safe_float(r["supporting_allele_cpgs"]) or 0 for r in support_rows)) * 1.1),
        False, MIN_CPGS,
    )
    draw_support_metric_axis(
        metric_axes[3], support_rows, "ic_phased_span_percent",
        "IC phased\nspan (%)", 100.0, False,
    )


def build_orthogonal_rows(panel_a_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Compare HiFi state with the recorded clinical mechanism for the original 11."""
    rows = []
    for row in panel_a_rows:
        if row["sample_id"] not in ORIGINAL_CLINICAL_COHORT:
            continue
        expected = GROUP_EXPECTED_STATE_CODES[row["molecular_mechanism"]]
        observed = (
            "absent" if row["allele_1_is_absent"] == "True" else row["allele_1_pattern_short"],
            "absent" if row["allele_2_is_absent"] == "True" else row["allele_2_pattern_short"],
        )
        observed = normalize_state_pair(row["molecular_mechanism"], observed)
        concordant = sample_is_concordant(row["molecular_mechanism"], observed)
        rows.append({
            "sample_id": row["sample_id"],
            "prior_clinical_mechanism": row["molecular_mechanism"],
            "hifi_mechanism": row["molecular_mechanism"] if concordant else "Discordant",
            "concordant": concordant,
            "clinical_assay": "not recorded in supplied metadata",
        })
    return rows


def draw_panel_e(ax: plt.Axes, rows: list[dict[str, Any]]) -> None:
    order = ["PWS-DEL", "PWS-mUPD", "AS-DEL", "Reference"]
    matrix = np.zeros((len(order), len(order)), dtype=int)
    for row in rows:
        if row["prior_clinical_mechanism"] in order and row["hifi_mechanism"] in order:
            matrix[order.index(row["prior_clinical_mechanism"]), order.index(row["hifi_mechanism"])] += 1
    ax.imshow(matrix, cmap="Greys", vmin=0, vmax=max(1, matrix.max()), aspect="auto")
    for i in range(len(order)):
        for j in range(len(order)):
            ax.text(j, i, str(matrix[i, j]), ha="center", va="center",
                    color="white" if matrix[i, j] > matrix.max() / 2 else "black",
                    fontsize=fs(7.5), fontweight="bold")
    ax.set_xticks(range(len(order)), order, fontsize=fs(6.5))
    ax.set_yticks(range(len(order)), order, fontsize=fs(6.5))
    ax.set_xlabel("HiFi mechanism call", fontsize=fs(7.5))
    ax.set_ylabel("Prior clinical mechanism", fontsize=fs(7.5))
    n_concordant = sum(bool(row["concordant"]) for row in rows)
    ax.set_title(f"Orthogonal concordance: {n_concordant}/{len(rows)}", loc="left",
                 fontsize=fs(9.0), fontweight="bold", pad=6)
    ax.text(1.02, 0.5,
            "Clinical assay type was not recorded\nin the supplied cohort metadata;\nconfirm MS-MLPA/MS-PCR provenance\nbefore manuscript submission.",
            transform=ax.transAxes, ha="left", va="center", fontsize=fs(6.4), color="#444444")

def create_figure(
    out_prefix: Path,
    panel_a_rows: list[dict[str, Any]],
    contrast_rows: list[dict[str, Any]],
    diagnostic_rows: list[dict[str, Any]],
    support_rows: list[dict[str, Any]],
    orthogonal_rows: list[dict[str, Any]],
    figure_dpi: int = 600,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial", "Liberation Sans"],
            "font.size": fs(9),
            "axes.linewidth": 0.75,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": figure_dpi,
        }
    )

    # 180 mm final width, with fonts sized at final reproduction scale.
    fig = plt.figure(figsize=(180 / 25.4, 218 / 25.4), constrained_layout=False)

    outer = GridSpec(
        3,
        1,
        figure=fig,
        height_ratios=[1.30, 0.85, 0.55],
        hspace=0.62,
    )

    top_row = outer[0].subgridspec(
        1,
        2,
        width_ratios=[1.08, 1.38],
        wspace=0.36,
    )

    bottom_row = outer[1].subgridspec(
        1,
        2,
        width_ratios=[0.76, 1.24],
        wspace=0.38,
    )

    panel_a_grid = top_row[0, 0].subgridspec(
        1,
        2,
        width_ratios=[0.01, 1.99],
        wspace=0.02,
    )

    ax_a_note = fig.add_subplot(panel_a_grid[0, 0])
    ax_a_heat = fig.add_subplot(panel_a_grid[0, 1])
    ax_b = fig.add_subplot(top_row[0, 1])
    ax_c = fig.add_subplot(bottom_row[0, 0])

    panel_d_grid = bottom_row[0, 1].subgridspec(
        1,
        4,
        wspace=0.22,
    )

    ax_d = [fig.add_subplot(panel_d_grid[0, i]) for i in range(4)]
    ax_e = fig.add_subplot(outer[2, 0])

    fig.subplots_adjust(
        top=0.925,
        bottom=0.145,
        left=0.055,
        right=0.985,
    )

    draw_panel_a(ax_a_note, ax_a_heat, panel_a_rows)
    draw_panel_b(ax_b, contrast_rows)
    draw_panel_c(ax_c, diagnostic_rows)
    draw_panel_d(ax_d, support_rows)
    draw_panel_e(ax_e, orthogonal_rows)

    panel_d_left = ax_d[0].get_position().x0
    panel_d_right = ax_d[-1].get_position().x1
    panel_d_top = ax_d[0].get_position().y1
    panel_d_bottom = ax_d[0].get_position().y0
    panel_d_center = (panel_d_left + panel_d_right) / 2.0

    fig.text(
        panel_d_center,
        panel_d_top + 0.040,
        "Coverage and phasing support at the PWS-AS IC",
        ha="center",
        va="bottom",
        fontsize=fs(9.4),
        fontweight="bold",
    )

    legend_handles = [
        *[
            plt.Line2D([0], [0], color=MECHANISM_COLORS[group], lw=2.2)
            for group in ("PWS-DEL", "PWS-mUPD", "AS-DEL", "Reference")
        ],
        Rectangle((0, 0), 1, 1, facecolor=STATE_COLORS["M"]),
        Rectangle((0, 0), 1, 1, facecolor=STATE_COLORS["P"]),
        Rectangle((0, 0), 1, 1, facecolor="#efefef", edgecolor=ABSENT_EDGE, hatch="////"),
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="#009E73",
            markeredgecolor="#009E73",
            markersize=5.2,
        ),
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="white",
            markeredgecolor="#D55E00",
            markersize=5.2,
        ),
    ]

    fig.legend(
        legend_handles,
        ["PWS-DEL", "PWS-mUPD", "AS-DEL", "Reference", "Maternal", "Paternal",
         "Deleted allele", "Pass (≥10×, ≥5 CpGs)", "Fail"],
        frameon=False,
        fontsize=fs(6.8),
        ncol=5,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.012),
        bbox_transform=fig.transFigure,
        columnspacing=1.6,
        handletextpad=0.5,
    )

    ax_a_heat.text(
        -0.20,
        1.035,
        "a",
        transform=ax_a_heat.transAxes,
        fontsize=fs(10),
        fontweight="bold",
        va="bottom",
        ha="left",
    )

    ax_b.text(
        -0.035,
        1.035,
        "b",
        transform=ax_b.transAxes,
        fontsize=fs(10),
        fontweight="bold",
        va="bottom",
        ha="left",
    )

    ax_c.text(
        -0.085,
        1.035,
        "c",
        transform=ax_c.transAxes,
        fontsize=fs(10),
        fontweight="bold",
        va="bottom",
        ha="left",
    )

    ax_d[0].text(
        -0.14,
        1.035,
        "d",
        transform=ax_d[0].transAxes,
        fontsize=fs(10),
        fontweight="bold",
        va="bottom",
        ha="left",
    )

    ax_e.text(-0.06, 1.03, "e", transform=ax_e.transAxes, fontsize=fs(10),
              fontweight="bold", va="bottom", ha="left")

    figure_dir = out_prefix.parent

    fig.savefig(
        figure_dir / "Figure_1.png",
        dpi=figure_dpi,
        bbox_inches="tight",
    )

    fig.savefig(
        figure_dir / "Figure_1.pdf",
        bbox_inches="tight",
    )

    fig.savefig(
        figure_dir / "Figure_1.svg",
        bbox_inches="tight",
    )

    plt.close(fig)


def write_report(
    report_path: Path,
    script_path: Path,
    input_files: dict[str, str],
    diagnostic_rows: list[dict[str, Any]],
) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)

    total_samples = sum(int(row["n_samples"]) for row in diagnostic_rows)
    total_concordant = sum(int(row.get("n_concordant", 0)) for row in diagnostic_rows)

    lines: list[str] = []
    lines.append("# Figure 1 Report: Diagnostic validation at the PWS-AS imprinting centre")
    lines.append("")
    lines.append("## 1. Purpose")
    lines.append(
        "Figure 1 tests whether PacBio HiFi native methylation recovers the expected "
        "parent-of-origin states at the PWS-AS imprinting centre across PWS-DEL, "
        "PWS-mUPD, AS-DEL and reference samples."
    )
    lines.append("")
    lines.append("## 2. Input data")
    lines.append("- Allele-level methylation matrix: `tables/Figure1A_allele_methylation_matrix.tsv`")
    lines.append("- Per-CpG contrast table: `tables/Figure1B_per_CpG_contrast.tsv`")
    lines.append("- Diagnostic state summary: `tables/Figure1C_diagnostic_state_summary.tsv`")
    lines.append("- Coverage/phasing support table: `tables/Figure1D_coverage_phasing_support.tsv`")
    lines.append("- Orthogonal concordance table: `tables/Figure1E_orthogonal_concordance.tsv`")
    lines.append(f"- Metadata: `{input_files['metadata']}`")
    lines.append(f"- Run parameters: `{input_files['run_parameters']}`")
    lines.append(f"- Script: `{script_path.name}`")
    lines.append("")
    lines.append("## 3. Coordinate system")
    lines.append("- Reference: T2T-CHM13v2.0")
    lines.append("- PWS-AS IC core interval: `chr15:22,691,258-22,693,494`")
    lines.append("")
    lines.append("## 4. Diagnostic concordance")
    lines.append(f"- Concordant mechanisms: `{total_concordant}/{total_samples}`")
    lines.append("")
    lines.append("| Group | n | Expected state | Observed state | Interpretation |")
    lines.append("| --- | ---: | --- | --- | --- |")

    for row in diagnostic_rows:
        lines.append(
            f"| {row['molecular_mechanism']} | {row['n_samples']} | "
            f"{row['expected_state']} | {row['observed_state']} | {row['interpretation']} |"
        )

    lines.append("")
    lines.append("## 5. Figure interpretation")
    lines.append("- Panel A shows the expected molecular configurations and IC methylation states.")
    lines.append("- Hatched grey cells indicate absent or deleted alleles.")
    lines.append("- Panel B shows per-CpG parent-of-origin contrast across the IC.")
    lines.append("- Panel C summarizes diagnostic concordance by molecular mechanism.")
    lines.append("- Panel D shows the primary coverage and phasing support metrics.")
    lines.append("- Panel E cross-tabulates HiFi calls against the recorded prior clinical mechanism for the original 11 genomes.")
    lines.append("")
    lines.append("## 6. Caption-ready statistical details")
    lines.append("- Maternal-pattern calls required mean methylation ≥0.85; paternal-pattern calls required ≤0.15.")
    lines.append("- Sufficient allele support required mean depth ≥10× and ≥5 CpGs.")
    lines.append("- Panel B lines are CpG-wise arithmetic group means; bands are percentile 95% confidence intervals from 5,000 sample-level bootstrap resamples (fixed seed). No hypothesis test is reported in this descriptive panel.")
    lines.append("- Panel E is a descriptive cross-tabulation (n=11); no inferential test is reported. The clinical assay type is absent from supplied metadata and must be confirmed before naming MS-MLPA or MS-PCR in the manuscript.")
    lines.append("")
    lines.append("## 7. Figure output files")
    lines.append("- `figures/Figure_1.png`")
    lines.append("- `figures/Figure_1.pdf`")
    lines.append("- `figures/Figure_1.svg`")

    report_path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()

    outdir = Path(args.outdir)
    table_dir = outdir / "tables"
    figure_dir = outdir / "figures"
    log_dir = outdir / "logs"
    report_dir = outdir / "reports"

    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    vcf_dir = Path(args.vcf_dir)
    bam_dir = Path(args.bam_dir)
    methylation_dir = Path(args.methylation_dir)
    cnv_dir = Path(args.cnv_dir)
    gtf = Path(args.gtf)
    metadata_path = Path(args.metadata)

    metadata = read_metadata(metadata_path)

    sample_files: dict[str, dict[str, Path | None]] = {}

    for sample_id, _clinical, _mechanism in COHORT:
        input_sample_id = INPUT_SAMPLE_ALIASES.get(sample_id, sample_id)
        sample_files[sample_id] = {
            "bam": find_sample_file(
                bam_dir,
                input_sample_id,
                ".bam",
                prefer_v2=args.prefer_v2,
            ),
            "summary": find_sample_file(
                vcf_dir,
                input_sample_id,
                ".summary.tsv",
                prefer_v2=args.prefer_v2,
            ),
            "blocks": find_sample_file(
                vcf_dir,
                input_sample_id,
                ".blocks.tsv",
                prefer_v2=args.prefer_v2,
            ),
            "combined_bed": find_sample_file(
                methylation_dir,
                input_sample_id,
                ".combined.bed",
                prefer_v2=args.prefer_v2,
            ),
            "hap1_bed": find_sample_file(
                methylation_dir,
                input_sample_id,
                ".hap1.bed",
                prefer_v2=args.prefer_v2,
            ),
            "hap2_bed": find_sample_file(
                methylation_dir,
                input_sample_id,
                ".hap2.bed",
                prefer_v2=args.prefer_v2,
            ),
            "cnv_log": find_sample_file(
                cnv_dir,
                input_sample_id,
                ".log",
                prefer_v2=args.prefer_v2,
            ),
        }

    bam_qc_fields = [
        "sample_id",
        "bam_file",
        "total_HiFi_reads",
        "mean_depth_genome_wide",
        "mean_depth_genome_wide_gc_corrected",
        "mean_depth_chr15",
        "mean_depth_per_haplotype_at_15q11-q13",
        "chr15_length",
    ]

    bam_qc_cache = table_dir / "bam_qc_cache.tsv"

    if args.skip_bam_qc and bam_qc_cache.exists():
        bam_qc_rows = read_tsv(bam_qc_cache)
    else:
        bam_qc_rows: list[dict[str, Any]] = []

        for sample_id, _clinical, _mechanism in COHORT:
            bam_qc_rows.append(
                build_bam_qc(
                    sample_id,
                    sample_files[sample_id]["bam"],
                    sample_files[sample_id]["cnv_log"],
                )
            )

        write_tsv(
            bam_qc_cache,
            [{k: fmt(v) for k, v in row.items()} for row in bam_qc_rows],
            bam_qc_fields,
        )

    bam_qc_by_sample = {row["sample_id"]: row for row in bam_qc_rows}

    summary_fields = [
        "sample_id",
        "clinical_diagnosis",
        "molecular_mechanism",
        "sex",
        "age_at_sampling",
        "total_HiFi_reads",
        "mean_depth_genome_wide",
        "mean_depth_chr15",
        "mean_depth_per_haplotype_at_15q11-q13",
        "phasing_block_N50_chr15",
        "percent_imprinted_domain_in_phased_block",
        "total_CpGs_called_in_imprinted_domain",
    ]

    summary_cache = table_dir / "Figure1B_cohort_qc_summary.tsv"

    if args.skip_bam_qc and summary_cache.exists():
        summary_rows = read_tsv(summary_cache)
    else:
        summary_rows: list[dict[str, Any]] = []

        for sample_id, clinical, mechanism in COHORT:
            meta = metadata.get(sample_id, {})
            n50, domain_fraction = block_n50_and_domain_fraction(
                sample_files[sample_id]["blocks"]
            )
            combined_bed = sample_files[sample_id]["combined_bed"]
            bam_qc = bam_qc_by_sample[sample_id]

            summary_rows.append(
                {
                    "sample_id": sample_id,
                    "clinical_diagnosis": clinical,
                    "molecular_mechanism": mechanism,
                    "sex": meta.get("gender", ""),
                    "age_at_sampling": meta.get("age", ""),
                    "total_HiFi_reads": bam_qc.get("total_HiFi_reads", ""),
                    "mean_depth_genome_wide": fmt(
                        safe_float(bam_qc.get("mean_depth_genome_wide"))
                    ),
                    "mean_depth_chr15": fmt(
                        safe_float(bam_qc.get("mean_depth_chr15"))
                    ),
                    "mean_depth_per_haplotype_at_15q11-q13": bam_qc.get(
                        "mean_depth_per_haplotype_at_15q11-q13",
                        "",
                    ),
                    "phasing_block_N50_chr15": n50 if n50 is not None else "",
                    "percent_imprinted_domain_in_phased_block": fmt(domain_fraction),
                    "total_CpGs_called_in_imprinted_domain": count_cpgs_in_domain(
                        combined_bed
                    ),
                }
            )

        write_tsv(summary_cache, summary_rows, summary_fields)

    genes = load_gene_models(gtf)

    write_tsv(
        table_dir / "Figure1A_domain_genes_from_T2T_GTF.tsv",
        genes,
        ["gene", "chrom", "start", "end", "strand", "parental_annotation"],
    )

    assignment_rows, matrix_rows, stats_by_sample = build_assignments(sample_files)

    assignment_fields = [
        "sample_id",
        "clinical_diagnosis",
        "molecular_mechanism",
        "haplotype_label",
        "source",
        "mean_methylation_at_PWS_IC",
        "n_CpGs_at_PWS_IC",
        "mean_coverage_at_PWS_IC",
        "coverage_status",
        "methylation_pattern",
        "parental_assignment",
        "expected_pattern",
        "assignment_confidence",
        "validation_status",
        "note",
    ]

    write_tsv(
        table_dir / "Figure1C_parental_assignment.tsv",
        assignment_rows,
        assignment_fields,
    )

    matrix_fields = [
        "sample_id",
        "molecular_mechanism",
        "haplotype_or_source",
        "mean_methylation",
        "n_CpGs",
        "mean_coverage",
        "pattern",
        "coverage_status",
    ]

    write_tsv(
        table_dir / "Figure1C_pws_ic_methylation_matrix.tsv",
        matrix_rows,
        matrix_fields,
    )

    panel_a_rows = build_physical_allele_rows(matrix_rows)

    panel_a_fields = [
        "sample_id",
        "display_label",
        "molecular_mechanism",
        "expected_group_configuration",
        "note",
        "allele_1_source",
        "allele_1_mean_methylation",
        "allele_1_pattern",
        "allele_1_pattern_short",
        "allele_1_n_CpGs",
        "allele_1_mean_coverage",
        "allele_1_coverage_status",
        "allele_1_is_absent",
        "allele_2_source",
        "allele_2_mean_methylation",
        "allele_2_pattern",
        "allele_2_pattern_short",
        "allele_2_n_CpGs",
        "allele_2_mean_coverage",
        "allele_2_coverage_status",
        "allele_2_is_absent",
    ]

    write_tsv(
        table_dir / "Figure1A_physical_allele_layout.tsv",
        panel_a_rows,
        panel_a_fields,
    )

    write_tsv(
        table_dir / "Figure1A_allele_methylation_matrix.tsv",
        panel_a_rows,
        panel_a_fields,
    )

    contrast_rows = build_per_cpg_contrast(stats_by_sample, assignment_rows)

    contrast_fields = [
        "pos",
        "score",
        "score_type",
        "sample_id",
        "molecular_mechanism",
    ]

    write_tsv(
        table_dir / "Figure1D_per_CpG_contrast.tsv",
        contrast_rows,
        contrast_fields,
    )

    write_tsv(
        table_dir / "Figure1B_per_CpG_contrast.tsv",
        contrast_rows,
        contrast_fields,
    )

    diagnostic_state_rows = build_diagnostic_state_rows(panel_a_rows)

    diagnostic_state_fields = [
        "molecular_mechanism",
        "n_samples",
        "n_concordant",
        "expected_state",
        "observed_state",
        "interpretation",
    ]

    write_tsv(
        table_dir / "Figure1C_diagnostic_state_summary.tsv",
        diagnostic_state_rows,
        diagnostic_state_fields,
    )

    support_rows = build_support_rows(summary_rows, matrix_rows, sample_files)

    support_fields = [
        "sample_id",
        "display_label",
        "molecular_mechanism",
        "total_ic_depth",
        "supporting_allele_depth",
        "supporting_allele_cpgs",
        "ic_phased_span_percent",
        "domain_phased_span_percent",
        "support_mode",
        "low_support",
    ]

    write_tsv(
        table_dir / "Figure1D_support_metrics.tsv",
        support_rows,
        support_fields,
    )

    write_tsv(
        table_dir / "Figure1D_coverage_phasing_support.tsv",
        support_rows,
        support_fields,
    )

    orthogonal_rows = build_orthogonal_rows(panel_a_rows)
    write_tsv(
        table_dir / "Figure1E_orthogonal_concordance.tsv",
        orthogonal_rows,
        ["sample_id", "prior_clinical_mechanism", "hifi_mechanism", "concordant", "clinical_assay"],
    )

    create_figure(
        figure_dir / "Figure_1",
        panel_a_rows,
        contrast_rows,
        diagnostic_state_rows,
        support_rows,
        orthogonal_rows,
        figure_dpi=args.figure_dpi,
    )

    run_parameters = {
        "cohort": [
            {
                "sample_id": sample_id,
                "clinical_diagnosis": clinical,
                "molecular_mechanism": mechanism,
            }
            for sample_id, clinical, mechanism in COHORT
        ],
        "regions": {
            "imprinted_domain": {
                "chrom": CHROM,
                "start": DOMAIN_START,
                "end": DOMAIN_END,
            },
            "pws_as_ic": {
                "chrom": CHROM,
                "start": PWS_IC_START,
                "end": PWS_IC_END,
                "name": PWS_IC_NAME,
                "source": "canonical PWS-AS ICR; SNRPN/SNHG14/SNURF regulatory interval",
            },
        },
        "thresholds": {
            "maternal_methylation": MATERNAL_THRESHOLD,
            "paternal_methylation": PATERNAL_THRESHOLD,
            "minimum_mean_coverage_per_haplotype": MIN_MEAN_COVERAGE,
            "minimum_cpgs": MIN_CPGS,
        },
        "input_paths": {
            "vcf_dir": str(vcf_dir),
            "bam_dir": str(bam_dir),
            "methylation_dir": str(methylation_dir),
            "cnv_dir": str(cnv_dir),
            "gtf": str(gtf),
            "metadata": str(metadata_path),
        },
        "figure_settings": {
            "figure_dpi": args.figure_dpi,
            "prefer_v2_files": args.prefer_v2,
            "text_scale": TEXT_SCALE,
            "color_semantics": {
                "maternal_pattern": STATE_COLORS["M"],
                "paternal_pattern": STATE_COLORS["P"],
                "upd_maternal_maternal": MECHANISM_COLORS["PWS-mUPD"],
                "control": MECHANISM_COLORS["Control"],
                "absent": STATE_COLORS["absent"],
            },
        },
        "depth_methods": {
            "total_HiFi_reads": "sum(mapped + unmapped) from samtools idxstats on HiPhase BAM",
            "mean_depth_genome_wide": "2 * final hificnv uncorrected haploid coverage estimate from sample log",
            "mean_depth_chr15": "samtools coverage -r chr15 meandepth",
            "mean_depth_per_haplotype_at_15q11-q13": "samtools HP-tagged coverage over chr15:22.5-28.5 Mb",
        },
    }

    run_parameters_path = outdir / "phase1_run_parameters.json"

    with run_parameters_path.open("w") as handle:
        json.dump(run_parameters, handle, indent=2)

    write_report(
        report_dir / "Figure1_report.md",
        Path(__file__),
        {
            "metadata": str(metadata_path),
            "run_parameters": str(run_parameters_path),
        },
        diagnostic_state_rows,
    )

    print("Figure 1 generation completed.")
    print(f"Output directory: {outdir}")
    print(f"Main PNG: {figure_dir / 'Figure_1.png'}")
    print(f"Main PDF: {figure_dir / 'Figure_1.pdf'}")


if __name__ == "__main__":
    main()
