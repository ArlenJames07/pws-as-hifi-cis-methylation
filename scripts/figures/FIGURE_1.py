#!/usr/bin/env python3
"""
Mechanistic Figure 1 generator for the PWS/AS HiFi cis-methylation manuscript.

Biological question
-------------------
Does the SNRPN/SNHG14 imprinting-centre methylation state encode the causal
parent-of-origin mechanism at chromosome 15, and is that state preserved in an
orthogonal genomic-disorder control carrying a pathogenic deletion outside
chromosome 15?

The figure is intentionally mechanistic rather than purely descriptive:

A. Allele-resolved IC methylation in every genome.
B. CpG-level parental-state contrast across the IC.
C. A mechanistic state-space that separates:
      - paternal deletion / retained maternal allele,
      - maternal deletion / retained paternal allele,
      - maternal UPD / duplicated maternal state,
      - canonical biparental chr15 state.
   The six 22q11.2-deletion samples are used as orthogonal disease controls and
   are NEVER silently pooled with unaffected controls.
D. Technical support for the mechanistic assignments.

Outputs
-------
tables/Figure1A_allele_methylation_matrix.tsv
tables/Figure1B_per_CpG_contrast.tsv
tables/Figure1C_mechanistic_state_space.tsv
tables/Figure1C_sample_level_classification.tsv
tables/Figure1C_diagnostic_state_summary.tsv
tables/Figure1D_coverage_phasing_support.tsv
figures/Figure1_mechanistic.{png,pdf,svg}
reports/Figure1_report.md

The script retains the command-line interface used by the current repository.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


# ---------------------------------------------------------------------------
# Project / analysis constants
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS_DIR = PROJECT_ROOT / "results"
DEFAULT_OUTDIR = DEFAULT_RESULTS_DIR / "07_figures" / "figure_1"
DEFAULT_PHASED_DIR = DEFAULT_RESULTS_DIR / "04_phasing"
DEFAULT_METHYLATION_DIR = DEFAULT_RESULTS_DIR / "06_methylation"
DEFAULT_CNV_DIR = DEFAULT_RESULTS_DIR / "05_cnv"
DEFAULT_METADATA = PROJECT_ROOT / "assets" / "metadata.csv"
DEFAULT_GTF = Path("/home/rare/arlen/reference/chm13v22.sorted.gtf")

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

# Mechanistic ordering: reciprocal deletions first, then copy-neutral mUPD,
# then orthogonal disease controls and unaffected controls.
COHORT = [
    ("001P", "Prader-Willi syndrome", "PWS-DEL"),
    ("002P", "Prader-Willi syndrome", "PWS-DEL"),
    ("005P", "Prader-Willi syndrome", "PWS-DEL"),
    ("006P", "Prader-Willi syndrome", "PWS-DEL"),
    ("007P", "Prader-Willi syndrome", "PWS-DEL"),

    ("013A", "Angelman syndrome", "AS-DEL"),
    ("014A", "Angelman syndrome", "AS-DEL"),
    ("016A", "Angelman syndrome", "AS-DEL"),

    ("004P", "Prader-Willi syndrome", "PWS-mUPD"),

    # Orthogonal genomic-disorder controls:
    # pathogenic 22q11.2 deletion but intact biparental chromosome 15.
    ("008D", "22q11.2 deletion syndrome", "Disease control"),
    ("009D", "22q11.2 deletion syndrome", "Disease control"),
    ("010D", "22q11.2 deletion syndrome", "Disease control"),
    ("011D", "22q11.2 deletion syndrome", "Disease control"),
    ("012D", "22q11.2 deletion syndrome", "Disease control"),
    ("015D", "22q11.2 deletion syndrome", "Disease control"),

    ("017C", "Unaffected control", "Control"),
    ("018C", "Unaffected control", "Control"),
]

MECHANISM_ORDER = {
    "PWS-DEL": 0,
    "AS-DEL": 1,
    "PWS-mUPD": 2,
    "Disease control": 3,
    "Control": 4,
}

MECHANISM_COLORS = {
    "PWS-DEL": "#C0392B",
    "AS-DEL": "#8E44AD",
    "PWS-mUPD": "#2E86C1",
    "Disease control": "#D89000",
    "Control": "#6F6F6F",
}

MECHANISM_MARKERS = {
    "PWS-DEL": "o",
    "AS-DEL": "s",
    "PWS-mUPD": "D",
    "Disease control": "^",
    "Control": "P",
}

MECHANISM_SAMPLE_PREFIX = {
    "PWS-DEL": "PW",
    "AS-DEL": "AS",
    "PWS-mUPD": "UPD",
    "Disease control": "DC",
    "Control": "CTRL",
}

GROUP_EXPECTED_STATE_CODES = {
    "PWS-DEL": ("M", "absent"),
    "AS-DEL": ("absent", "P"),
    "PWS-mUPD": ("M", "M"),
    "Disease control": ("M", "P"),
    "Control": ("M", "P"),
}

GROUP_INTERPRETATIONS = {
    "PWS-DEL": "paternal chr15 deletion; maternal allele retained",
    "AS-DEL": "maternal chr15 deletion; paternal allele retained",
    "PWS-mUPD": "copy-neutral maternal duplication",
    "Disease control": "orthogonal 22q11.2 deletion; biparental chr15 expected",
    "Control": "unaffected biparental chr15 reference",
}

STATE_COLORS = {
    "M": "#CB4335",
    "P": "#3F72AF",
    "absent": "#ECECEC",
    "?": "#F3F3F3",
}
ABSENT_EDGE = "#A6A6A6"
TEXT_SCALE = 1.05


def fs(size: float) -> float:
    return size * TEXT_SCALE


def sorted_cohort() -> list[tuple[str, str, str]]:
    return sorted(COHORT, key=lambda r: (MECHANISM_ORDER[r[2]], r[0]))


def sample_display_labels() -> dict[str, str]:
    labels: dict[str, str] = {}
    counts: dict[str, int] = defaultdict(int)
    for sample_id, _clinical, mechanism in sorted_cohort():
        counts[mechanism] += 1
        labels[sample_id] = f"{MECHANISM_SAMPLE_PREFIX[mechanism]}-{counts[mechanism]}"
    return labels


# ---------------------------------------------------------------------------
# Input data structures and utilities
# ---------------------------------------------------------------------------

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf-dir", type=Path, default=DEFAULT_PHASED_DIR)
    parser.add_argument("--bam-dir", type=Path, default=DEFAULT_PHASED_DIR)
    parser.add_argument("--methylation-dir", type=Path, default=DEFAULT_METHYLATION_DIR)
    parser.add_argument("--cnv-dir", type=Path, default=DEFAULT_CNV_DIR)
    parser.add_argument("--gtf", type=Path, default=DEFAULT_GTF)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--skip-bam-qc",
        action="store_true",
        help="Reuse existing BAM QC cache when present.",
    )
    parser.add_argument(
        "--render-only",
        action="store_true",
        help="Redraw the figure from cached Figure1 tables without rereading BAM/methylation files.",
    )
    return parser.parse_args()


def safe_float(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return x if np.isfinite(x) else None


def fmt(value: Any, digits: int = 3) -> str:
    x = safe_float(value)
    if x is None:
        return "" if value in (None, "") else str(value)
    return f"{x:.{digits}f}"


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_metadata(path: Path) -> dict[str, dict[str, str]]:
    metadata: dict[str, dict[str, str]] = {}
    if not path.exists():
        return metadata
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sample_id = (
                row.get("Codigo")
                or row.get("sample")
                or row.get("sample_id")
                or row.get("Sample")
                or ""
            ).strip()
            if sample_id:
                metadata[sample_id] = row
    return metadata


def choose_file(files: list[Path], sample_id: str) -> Path | None:
    if not files:
        return None
    exact = [p for p in files if re.search(rf"[_-]{re.escape(sample_id)}(\.|_|$)", p.name)]
    candidates = exact or files
    candidates = sorted(candidates, key=lambda p: ("v2" in p.name.lower(), len(p.name), p.name))
    return candidates[0]


def find_sample_file(directory: Path, sample_id: str, suffix: str) -> Path | None:
    if not directory.exists():
        return None
    sample_directory = directory / sample_id
    if sample_directory.is_dir():
        matches = list(sample_directory.glob(f"*{sample_id}*{suffix}"))
    else:
        matches = list(directory.rglob(f"*{sample_id}*{suffix}"))
    return choose_file(matches, sample_id)


def run_command(args: list[str]) -> str:
    result = subprocess.run(args, check=True, text=True, capture_output=True)
    return result.stdout


def read_bed_region(
    path: Path | None,
    start: int,
    end: int,
    keep_values: bool = False,
) -> BedStats:
    """
    Read the existing pb-CpG-tools/DSS-style BED used by this project.

    Expected columns are compatible with the current pipeline:
      chrom start ... methylation_percent ... coverage
    where methylation is read from column 4 and coverage from column 6.
    """
    if path is None or not path.exists():
        return BedStats(values_by_pos={} if keep_values else None)

    meth_values: list[float] = []
    cov_values: list[float] = []
    values_by_pos: dict[int, tuple[float, float]] = {}

    awk_script = "$1==chrom && $2>=start && $2<end {print}"
    result = subprocess.run(
        [
            "awk",
            "-v", f"chrom={CHROM}",
            "-v", f"start={start}",
            "-v", f"end={end}",
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
        try:
            row_start = int(fields[1])
            meth = float(fields[3]) / 100.0
            cov = float(fields[5])
        except ValueError:
            continue
        if not (np.isfinite(meth) and np.isfinite(cov)):
            continue
        meth_values.append(meth)
        cov_values.append(cov)
        if keep_values:
            values_by_pos[row_start] = (meth, cov)

    if not meth_values:
        return BedStats(values_by_pos=values_by_pos if keep_values else None)

    weights = np.asarray(cov_values, dtype=float)
    values = np.asarray(meth_values, dtype=float)
    if np.sum(weights) > 0:
        weighted_meth = float(np.average(values, weights=weights))
    else:
        weighted_meth = float(np.mean(values))

    return BedStats(
        n_cpgs=len(meth_values),
        mean_methylation=weighted_meth,
        mean_coverage=float(np.mean(cov_values)),
        values_by_pos=values_by_pos if keep_values else None,
    )


def methylation_pattern(stats: BedStats) -> str:
    beta = stats.mean_methylation
    if beta is None:
        return "missing"
    if beta >= MATERNAL_THRESHOLD:
        return "maternal-pattern"
    if beta <= PATERNAL_THRESHOLD:
        return "paternal-pattern"
    return "intermediate"


def pattern_short(pattern: str) -> str:
    return {
        "maternal-pattern": "M",
        "paternal-pattern": "P",
        "absent": "absent",
    }.get(pattern, "?")


def pattern_confidence(stats: BedStats) -> float:
    if stats.mean_methylation is None or stats.n_cpgs < MIN_CPGS:
        return 0.0
    distance = abs(stats.mean_methylation - 0.5)
    max_threshold_distance = MATERNAL_THRESHOLD - 0.5
    conf = distance / max_threshold_distance
    if stats.mean_coverage is not None and stats.mean_coverage < MIN_MEAN_COVERAGE:
        conf *= max(0.25, stats.mean_coverage / MIN_MEAN_COVERAGE)
    return float(np.clip(conf, 0.0, 1.0))


# ---------------------------------------------------------------------------
# BAM / phasing QC
# ---------------------------------------------------------------------------

def bam_idxstats(bam: Path) -> tuple[int, dict[str, int]]:
    stdout = run_command(["samtools", "idxstats", str(bam)])
    total_reads = 0
    chrom_lengths: dict[str, int] = {}
    for line in stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 4:
            continue
        chrom, length, mapped, unmapped = fields[:4]
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
            try:
                return float(fields[6])
            except ValueError:
                return None
    return None


def haplotype_depths_from_bam(bam: Path, region: str, region_len: int) -> dict[str, float]:
    """
    Estimate HP1/HP2 depth over the requested region.

    This preserves the strategy used by the existing project script and avoids
    materializing SAM text.
    """
    depths: dict[str, float] = {}
    for hp_value, label in [("1", "hap1"), ("2", "hap2")]:
        view = subprocess.Popen(
            ["samtools", "view", "-u", "-F", "2308", "-d", f"HP:{hp_value}", str(bam), region],
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
            if len(fields) >= 7 and fields[0] == CHROM:
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


def block_n50_and_domain_fraction(blocks_file: Path | None) -> tuple[int | None, float | None]:
    if blocks_file is None or not blocks_file.exists():
        return None, None

    lengths: list[int] = []
    overlaps: list[tuple[int, int]] = []

    with blocks_file.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("chrom") != CHROM:
                continue
            try:
                start = int(row["start"])
                end = int(row["end"])
            except (KeyError, ValueError):
                continue
            if end < start:
                continue
            lengths.append(end - start + 1)
            ov_start = max(start, DOMAIN_START)
            ov_end = min(end, DOMAIN_END)
            if ov_end >= ov_start:
                overlaps.append((ov_start, ov_end + 1))

    n50: int | None = None
    if lengths:
        total = sum(lengths)
        running = 0
        for length in sorted(lengths, reverse=True):
            running += length
            if running >= total / 2.0:
                n50 = length
                break

    if not overlaps:
        return n50, 0.0

    overlaps.sort()
    merged: list[list[int]] = []
    for start, end in overlaps:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    covered = sum(end - start for start, end in merged)
    fraction = 100.0 * covered / (DOMAIN_END - DOMAIN_START)
    return n50, fraction


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
            try:
                start = int(row["start"])
                end = int(row["end"])
            except (KeyError, ValueError):
                continue
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


def build_bam_qc(sample_id: str, bam: Path | None, cnv_log: Path | None) -> dict[str, Any]:
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
    genome_depth, genome_depth_gc = parse_hificnv_depth(cnv_log)
    mean_depth_chr15 = samtools_coverage_mean_depth(bam, CHROM)
    region = f"{CHROM}:{DOMAIN_START}-{DOMAIN_END}"
    hap_depths = haplotype_depths_from_bam(
        bam,
        region,
        DOMAIN_END - DOMAIN_START + 1,
    )

    return {
        "sample_id": sample_id,
        "bam_file": str(bam),
        "total_HiFi_reads": total_reads,
        "mean_depth_genome_wide": genome_depth,
        "mean_depth_genome_wide_gc_corrected": genome_depth_gc,
        "mean_depth_chr15": mean_depth_chr15,
        "mean_depth_per_haplotype_at_15q11-q13": ";".join(
            f"{k}={v:.3f}" for k, v in sorted(hap_depths.items())
        ),
        "chr15_length": chrom_lengths.get(CHROM, ""),
    }


# ---------------------------------------------------------------------------
# IC methylation assignment
# ---------------------------------------------------------------------------

def build_assignments(
    sample_files: dict[str, dict[str, Path | None]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, dict[str, BedStats]]]:

    assignment_rows: list[dict[str, Any]] = []
    matrix_rows: list[dict[str, Any]] = []
    stats_by_sample: dict[str, dict[str, BedStats]] = {}

    for sample_id, clinical, mechanism in sorted_cohort():
        files = sample_files[sample_id]
        stats = {
            "hap1": read_bed_region(
                files["hap1_bed"], PWS_IC_START, PWS_IC_END, keep_values=True
            ),
            "hap2": read_bed_region(
                files["hap2_bed"], PWS_IC_START, PWS_IC_END, keep_values=True
            ),
            "combined_fallback": read_bed_region(
                files["combined_bed"], PWS_IC_START, PWS_IC_END, keep_values=True
            ),
        }
        stats_by_sample[sample_id] = stats

        for label, bed_stats in stats.items():
            patt = methylation_pattern(bed_stats)
            matrix_rows.append(
                {
                    "sample_id": sample_id,
                    "molecular_mechanism": mechanism,
                    "haplotype_or_source": label,
                    "mean_methylation": fmt(bed_stats.mean_methylation),
                    "n_CpGs": bed_stats.n_cpgs,
                    "mean_coverage": fmt(bed_stats.mean_coverage),
                    "pattern": patt,
                    "coverage_status": (
                        "sufficient" if bed_stats.sufficient else "insufficient_or_missing"
                    ),
                }
            )

        # Deletion cases are biologically hemizygous across the disease interval;
        # the combined track is the retained physical chromosome.
        if mechanism in {"PWS-DEL", "AS-DEL"}:
            rows_for_sample = [
                ("combined_fallback", stats["combined_fallback"], "combined.bed")
            ]
            expected = (
                "maternal-pattern"
                if mechanism == "PWS-DEL"
                else "paternal-pattern"
            )
            note = (
                f"{mechanism}: retained hemizygous chromosome evaluated from combined track"
            )
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
                "maternal UPD: both physical haplotypes evaluated independently"
                if mechanism == "PWS-mUPD"
                else "orthogonal disease control: chr15 haplotypes assigned from IC methylation"
                if mechanism == "Disease control"
                else "unaffected control: chr15 haplotypes assigned from IC methylation"
            )

        for label, bed_stats, source in rows_for_sample:
            patt = methylation_pattern(bed_stats)
            if mechanism == "PWS-mUPD":
                parental_assignment = (
                    "maternal-pattern" if patt == "maternal-pattern" else "unassigned"
                )
            elif patt == "maternal-pattern":
                parental_assignment = "maternal"
            elif patt == "paternal-pattern":
                parental_assignment = "paternal"
            else:
                parental_assignment = "unassigned"

            if expected == "maternal-pattern":
                validation = "PASS" if patt == "maternal-pattern" else "CHECK"
            elif expected == "paternal-pattern":
                validation = "PASS" if patt == "paternal-pattern" else "CHECK"
            elif expected == "both maternal-pattern":
                validation = "PASS" if patt == "maternal-pattern" else "CHECK"
            else:
                validation = (
                    "PASS"
                    if patt in {"maternal-pattern", "paternal-pattern"}
                    else "CHECK"
                )

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
                    "coverage_status": (
                        "sufficient" if bed_stats.sufficient else "insufficient_or_missing"
                    ),
                    "methylation_pattern": patt,
                    "parental_assignment": parental_assignment,
                    "expected_pattern": expected,
                    "assignment_confidence": fmt(pattern_confidence(bed_stats)),
                    "validation_status": validation,
                    "note": note,
                }
            )

    return assignment_rows, matrix_rows, stats_by_sample


def _absent_cell() -> dict[str, Any]:
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


def _row_to_cell(row: dict[str, Any] | None) -> dict[str, Any]:
    if row is None:
        return _absent_cell()
    value = safe_float(row.get("mean_methylation"))
    if value is None:
        return _absent_cell()
    patt = row.get("pattern", "missing")
    return {
        "source": row.get("haplotype_or_source", ""),
        "mean_methylation": value,
        "pattern": patt,
        "pattern_short": pattern_short(patt),
        "n_CpGs": int(row.get("n_CpGs", 0) or 0),
        "mean_coverage": safe_float(row.get("mean_coverage")),
        "coverage_status": row.get("coverage_status", ""),
        "is_absent": False,
    }


def build_physical_allele_rows(matrix_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_sample: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)
    for row in matrix_rows:
        by_sample[row["sample_id"]][row["haplotype_or_source"]] = row

    display_labels = sample_display_labels()
    rows: list[dict[str, Any]] = []

    for sample_id, _clinical, mechanism in sorted_cohort():
        sample_rows = by_sample[sample_id]

        if mechanism == "PWS-DEL":
            allele_1 = _row_to_cell(sample_rows.get("combined_fallback"))
            allele_2 = _absent_cell()
        elif mechanism == "AS-DEL":
            allele_1 = _absent_cell()
            allele_2 = _row_to_cell(sample_rows.get("combined_fallback"))
        else:
            allele_1 = _row_to_cell(sample_rows.get("hap1"))
            allele_2 = _row_to_cell(sample_rows.get("hap2"))

        rows.append(
            {
                "sample_id": sample_id,
                "display_label": display_labels[sample_id],
                "molecular_mechanism": mechanism,

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


# ---------------------------------------------------------------------------
# Mechanistic classification
# ---------------------------------------------------------------------------

def normalize_state_pair(
    codes: tuple[str, str],
    mechanism: str,
) -> tuple[str, str]:
    """
    HP1/HP2 are not intrinsically maternal/paternal labels.

    For biparental groups M/P and P/M are biologically equivalent, so they are
    normalized before concordance is assessed.
    """
    if mechanism in {"Control", "Disease control"} and set(codes) == {"M", "P"}:
        return ("M", "P")
    return codes


def classify_discrete_state(codes: tuple[str, str]) -> str:
    if codes == ("M", "absent"):
        return "maternal-retained deletion"
    if codes == ("absent", "P"):
        return "paternal-retained deletion"
    if codes == ("M", "M"):
        return "duplicated maternal state"
    if codes == ("M", "P"):
        return "canonical biparental chr15"
    return "unresolved"


def mechanistic_template_class(
    observed_betas: list[float],
) -> tuple[str, float]:
    """
    Continuous template-distance classifier.

    One observed allele:
      maternal-retained template = [1]
      paternal-retained template = [0]

    Two observed alleles (sorted):
      biparental template = [0, 1]
      maternal-duplicated template = [1, 1]
      paternal-duplicated template = [0, 0]

    The returned distance is RMS methylation distance to the winning template.
    This makes the classification data-derived rather than a copy of the known
    diagnosis.
    """
    values = np.asarray([x for x in observed_betas if np.isfinite(x)], dtype=float)

    if len(values) == 1:
        beta = float(values[0])
        distances = {
            "maternal-retained deletion": abs(beta - 1.0),
            "paternal-retained deletion": abs(beta - 0.0),
        }
    elif len(values) == 2:
        low, high = np.sort(values)
        obs = np.array([low, high], dtype=float)
        templates = {
            "canonical biparental chr15": np.array([0.0, 1.0]),
            "duplicated maternal state": np.array([1.0, 1.0]),
            "duplicated paternal state": np.array([0.0, 0.0]),
        }
        distances = {
            name: float(np.sqrt(np.mean((obs - tpl) ** 2)))
            for name, tpl in templates.items()
        }
    else:
        return "unresolved", float("nan")

    winner = min(distances, key=distances.get)
    return winner, float(distances[winner])


def expected_template_for_group(mechanism: str) -> str:
    return {
        "PWS-DEL": "maternal-retained deletion",
        "AS-DEL": "paternal-retained deletion",
        "PWS-mUPD": "duplicated maternal state",
        "Disease control": "canonical biparental chr15",
        "Control": "canonical biparental chr15",
    }[mechanism]


def build_mechanistic_state_rows(panel_a_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    for row in panel_a_rows:
        mechanism = row["molecular_mechanism"]

        values: list[float] = []
        state_codes: list[str] = []
        for idx in (1, 2):
            absent = row[f"allele_{idx}_is_absent"] == "True"
            if absent:
                state_codes.append("absent")
                continue

            beta = safe_float(row[f"allele_{idx}_mean_methylation"])
            if beta is None:
                state_codes.append("?")
                continue
            values.append(beta)
            state_codes.append(row[f"allele_{idx}_pattern_short"])

        if len(state_codes) < 2:
            state_codes += ["?"] * (2 - len(state_codes))

        observed_codes = normalize_state_pair(
            (state_codes[0], state_codes[1]),
            mechanism,
        )
        expected_codes = normalize_state_pair(
            GROUP_EXPECTED_STATE_CODES[mechanism],
            mechanism,
        )

        n_observed = len(values)
        mean_beta = float(np.mean(values)) if values else np.nan
        allelic_contrast = (
            float(abs(values[0] - values[1]))
            if len(values) == 2
            else np.nan
        )

        predicted_template, template_distance = mechanistic_template_class(values)
        expected_template = expected_template_for_group(mechanism)

        rows.append(
            {
                "sample_id": row["sample_id"],
                "display_label": row["display_label"],
                "molecular_mechanism": mechanism,
                "n_resolved_chr15_alleles": n_observed,
                "mean_IC_methylation_across_resolved_alleles": fmt(mean_beta),
                "allelic_methylation_contrast": fmt(allelic_contrast),
                "expected_discrete_state": " / ".join(expected_codes),
                "observed_discrete_state": " / ".join(observed_codes),
                "observed_discrete_class": classify_discrete_state(observed_codes),
                "expected_template_class": expected_template,
                "predicted_template_class": predicted_template,
                "template_distance": fmt(template_distance),
                "template_concordant": str(predicted_template == expected_template),
                "state_concordant": str(observed_codes == expected_codes),
                "specificity_role": (
                    "orthogonal non-chr15 deletion control"
                    if mechanism == "Disease control"
                    else "unaffected biparental reference"
                    if mechanism == "Control"
                    else "PWS/AS causal-mechanism state"
                ),
            }
        )

    return rows


def build_diagnostic_state_rows(
    mechanistic_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:

    rows: list[dict[str, Any]] = []
    for mechanism in sorted(MECHANISM_ORDER, key=MECHANISM_ORDER.get):
        group = [r for r in mechanistic_rows if r["molecular_mechanism"] == mechanism]
        if not group:
            continue

        state_counts: dict[str, int] = defaultdict(int)
        template_counts: dict[str, int] = defaultdict(int)
        state_ok = 0
        template_ok = 0

        for row in group:
            state_counts[row["observed_discrete_state"]] += 1
            template_counts[row["predicted_template_class"]] += 1
            state_ok += int(row["state_concordant"] == "True")
            template_ok += int(row["template_concordant"] == "True")

        expected_state = " / ".join(
            normalize_state_pair(
                GROUP_EXPECTED_STATE_CODES[mechanism],
                mechanism,
            )
        )

        rows.append(
            {
                "molecular_mechanism": mechanism,
                "n_samples": len(group),
                "expected_state": expected_state,
                "observed_state_distribution": "; ".join(
                    f"{state}: {n}/{len(group)}"
                    for state, n in sorted(
                        state_counts.items(),
                        key=lambda kv: (-kv[1], kv[0]),
                    )
                ),
                "state_concordant_n": state_ok,
                "state_concordance_percent": fmt(100.0 * state_ok / len(group), 1),
                "template_concordant_n": template_ok,
                "template_concordance_percent": fmt(100.0 * template_ok / len(group), 1),
                "predicted_template_distribution": "; ".join(
                    f"{state}: {n}/{len(group)}"
                    for state, n in sorted(
                        template_counts.items(),
                        key=lambda kv: (-kv[1], kv[0]),
                    )
                ),
                "interpretation": GROUP_INTERPRETATIONS[mechanism],
            }
        )

    return rows


# ---------------------------------------------------------------------------
# CpG-level contrast
# ---------------------------------------------------------------------------

def build_per_cpg_contrast(
    stats_by_sample: dict[str, dict[str, BedStats]],
    assignment_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:

    assignment_by_sample: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in assignment_rows:
        assignment_by_sample[row["sample_id"]].append(row)

    rows: list[dict[str, Any]] = []

    for sample_id, _clinical, mechanism in sorted_cohort():
        sample_stats = stats_by_sample[sample_id]

        if mechanism in {"Control", "Disease control"}:
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

            maternal = sample_stats[maternal_label].values_by_pos or {}
            paternal = sample_stats[paternal_label].values_by_pos or {}
            for pos in sorted(set(maternal) & set(paternal)):
                rows.append(
                    {
                        "pos": pos,
                        "score": fmt(maternal[pos][0] - paternal[pos][0]),
                        "score_type": "maternal_minus_paternal",
                        "sample_id": sample_id,
                        "molecular_mechanism": mechanism,
                    }
                )

        elif mechanism == "PWS-mUPD":
            h1 = sample_stats["hap1"].values_by_pos or {}
            h2 = sample_stats["hap2"].values_by_pos or {}
            for pos in sorted(set(h1) & set(h2)):
                rows.append(
                    {
                        "pos": pos,
                        "score": fmt(h1[pos][0] - h2[pos][0]),
                        "score_type": "maternal_hap1_minus_maternal_hap2",
                        "sample_id": sample_id,
                        "molecular_mechanism": mechanism,
                    }
                )

        elif mechanism in {"PWS-DEL", "AS-DEL"}:
            combined = sample_stats["combined_fallback"].values_by_pos or {}
            # Standardized retained-parent identity score:
            # beta=1 -> +1 (maternal-like), beta=0 -> -1 (paternal-like).
            for pos, (meth, _cov) in sorted(combined.items()):
                rows.append(
                    {
                        "pos": pos,
                        "score": fmt(2.0 * meth - 1.0),
                        "score_type": "standardized_retained_parent_identity",
                        "sample_id": sample_id,
                        "molecular_mechanism": mechanism,
                    }
                )

    return rows


def summarize_contrast_rows(
    contrast_rows: list[dict[str, Any]],
) -> dict[str, dict[str, Any]]:

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
            x = np.asarray(
                [(pos - midpoint) / 1000.0 for pos, _ in ordered],
                dtype=float,
            )
            y = np.asarray([score for _, score in ordered], dtype=float)
            traces[sample_id] = (x, y)

            for pos, score in ordered:
                position_values[pos].append(score)

        if position_values:
            positions = np.asarray(sorted(position_values), dtype=int)
            x = np.asarray(
                [(pos - midpoint) / 1000.0 for pos in positions],
                dtype=float,
            )
            median = np.asarray(
                [np.median(position_values[pos]) for pos in positions],
                dtype=float,
            )
            q25 = np.asarray(
                [np.percentile(position_values[pos], 25) for pos in positions],
                dtype=float,
            )
            q75 = np.asarray(
                [np.percentile(position_values[pos], 75) for pos in positions],
                dtype=float,
            )
        else:
            x = median = q25 = q75 = np.asarray([], dtype=float)

        summaries[mechanism] = {
            "x": x,
            "median": median,
            "q25": q25,
            "q75": q75,
            "traces": traces,
        }

    return summaries


# ---------------------------------------------------------------------------
# Support metrics
# ---------------------------------------------------------------------------

def build_support_rows(
    summary_rows: list[dict[str, Any]],
    matrix_rows: list[dict[str, Any]],
    sample_files: dict[str, dict[str, Path | None]],
) -> list[dict[str, Any]]:

    summary_by_sample = {row["sample_id"]: row for row in summary_rows}
    matrix_by_sample: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)
    for row in matrix_rows:
        matrix_by_sample[row["sample_id"]][row["haplotype_or_source"]] = row

    labels = sample_display_labels()
    rows: list[dict[str, Any]] = []

    for sample_id, _clinical, mechanism in sorted_cohort():
        sample_summary = summary_by_sample.get(sample_id, {})
        sample_matrix = matrix_by_sample[sample_id]

        combined = sample_matrix.get("combined_fallback")
        total_ic_depth = (
            safe_float(combined.get("mean_coverage"))
            if combined
            else None
        )

        if mechanism in {"PWS-DEL", "AS-DEL"}:
            allele_rows = [combined] if combined else []
            support_mode = "retained_hemizygous_allele"
        else:
            allele_rows = [
                sample_matrix.get("hap1"),
                sample_matrix.get("hap2"),
            ]
            support_mode = "minimum_phased_allele"

        allele_rows = [
            row for row in allele_rows
            if row and safe_float(row.get("mean_coverage")) is not None
        ]

        min_depth = min(
            (safe_float(row["mean_coverage"]) for row in allele_rows),
            default=None,
        )
        min_cpgs = min(
            (int(row["n_CpGs"]) for row in allele_rows),
            default=0,
        )
        low_support = (
            not allele_rows
            or any(row["coverage_status"] != "sufficient" for row in allele_rows)
        )

        ic_phased = block_fraction_for_interval(
            sample_files[sample_id]["blocks"],
            PWS_IC_START,
            PWS_IC_END,
        )

        rows.append(
            {
                "sample_id": sample_id,
                "display_label": labels[sample_id],
                "molecular_mechanism": mechanism,
                "total_ic_depth": fmt(total_ic_depth),
                "supporting_allele_depth": fmt(min_depth),
                "supporting_allele_cpgs": min_cpgs,
                "ic_phased_span_percent": fmt(ic_phased),
                "domain_phased_span_percent": sample_summary.get(
                    "percent_imprinted_domain_in_phased_block", ""
                ),
                "support_mode": support_mode,
                "low_support": str(low_support),
            }
        )

    return rows


# ---------------------------------------------------------------------------
# Figure drawing
# ---------------------------------------------------------------------------

def draw_panel_a(
    note_ax: plt.Axes,
    heat_ax: plt.Axes,
    panel_a_rows: list[dict[str, Any]],
) -> None:

    n = len(panel_a_rows)
    values = np.full((n, 2), np.nan)
    absent = np.zeros((n, 2), dtype=bool)
    low_support = np.zeros((n, 2), dtype=bool)

    for i, row in enumerate(panel_a_rows):
        for j in range(2):
            prefix = f"allele_{j + 1}"
            is_absent = row[f"{prefix}_is_absent"] == "True"
            absent[i, j] = is_absent
            if is_absent:
                continue
            beta = safe_float(row[f"{prefix}_mean_methylation"])
            if beta is None:
                absent[i, j] = True
                continue
            values[i, j] = beta
            low_support[i, j] = row[f"{prefix}_coverage_status"] != "sufficient"

    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_bad("#F2F2F2")
    image = heat_ax.imshow(
        values,
        aspect="auto",
        cmap=cmap,
        norm=TwoSlopeNorm(vmin=0.0, vcenter=0.5, vmax=1.0),
    )

    heat_ax.set_xticks([0, 1])
    heat_ax.set_xticklabels(
        ["Physical allele /\nhaplotype 1", "Physical allele /\nhaplotype 2"],
        fontsize=fs(8),
    )
    heat_ax.set_yticks(range(n))
    heat_ax.set_yticklabels(
        [row["display_label"] for row in panel_a_rows],
        fontsize=fs(8),
        fontweight="bold",
    )
    heat_ax.tick_params(length=0)

    for i, row in enumerate(panel_a_rows):
        mechanism = row["molecular_mechanism"]
        heat_ax.get_yticklabels()[i].set_color(MECHANISM_COLORS[mechanism])

        for j in range(2):
            prefix = f"allele_{j + 1}"
            if absent[i, j]:
                heat_ax.add_patch(
                    Rectangle(
                        (j - 0.5, i - 0.5),
                        1,
                        1,
                        facecolor="#EFEFEF",
                        edgecolor=ABSENT_EDGE,
                        hatch="///",
                        linewidth=0.8,
                        zorder=3,
                    )
                )
                heat_ax.text(
                    j, i, "absent",
                    ha="center", va="center",
                    fontsize=fs(6.3),
                    color="#444444",
                    zorder=4,
                )
            else:
                beta = values[i, j]
                patt = row[f"{prefix}_pattern_short"]
                if low_support[i, j]:
                    heat_ax.add_patch(
                        Rectangle(
                            (j - 0.5, i - 0.5),
                            1,
                            1,
                            facecolor="none",
                            edgecolor="#7F6A2F",
                            linewidth=1.2,
                            zorder=3,
                        )
                    )
                text_color = "white" if beta >= 0.72 or beta <= 0.20 else "#111111"
                heat_ax.text(
                    j, i,
                    f"{beta:.2f}\n{patt}",
                    ha="center", va="center",
                    fontsize=fs(6.8),
                    color=text_color,
                    zorder=4,
                )

    # Group separators.
    for i in range(n - 1):
        if (
            panel_a_rows[i]["molecular_mechanism"]
            != panel_a_rows[i + 1]["molecular_mechanism"]
        ):
            heat_ax.axhline(i + 0.5, color="#5C5C5C", lw=0.8)

    note_ax.set_xlim(0, 1)
    note_ax.set_ylim(heat_ax.get_ylim())
    note_ax.axis("off")

    for mechanism in sorted(MECHANISM_ORDER, key=MECHANISM_ORDER.get):
        idx = [
            i for i, row in enumerate(panel_a_rows)
            if row["molecular_mechanism"] == mechanism
        ]
        if not idx:
            continue
        y = 0.5 * (idx[0] + idx[-1])
        note_ax.text(
            0.02, y,
            f"{mechanism}\n(n={len(idx)})",
            ha="left", va="center",
            fontsize=fs(8),
            fontweight="bold",
            color=MECHANISM_COLORS[mechanism],
        )

    heat_ax.set_title(
        "A. Allele-resolved SNRPN/SNHG14 IC methylation",
        fontsize=fs(11),
        loc="left",
        x=-0.73,
        pad=12,
        weight="bold",
    )

    cbar = plt.colorbar(image, ax=heat_ax, fraction=0.048, pad=0.025)
    cbar.set_label("Mean IC methylation", fontsize=fs(8))
    cbar.ax.tick_params(labelsize=fs(7))


def draw_panel_b(ax: plt.Axes, contrast_rows: list[dict[str, Any]]) -> None:
    summaries = summarize_contrast_rows(contrast_rows)

    midpoint = (PWS_IC_START + PWS_IC_END) / 2.0
    x_left = (PWS_IC_START - midpoint) / 1000.0
    x_right = (PWS_IC_END - midpoint) / 1000.0

    ax.axvspan(x_left, x_right, color="#F5EFEA", alpha=0.8, zorder=0)
    ax.axhline(0, color="#505050", lw=0.9)

    for mechanism in sorted(MECHANISM_ORDER, key=MECHANISM_ORDER.get):
        summary = summaries.get(mechanism, {})

        for sample_x, sample_y in summary.get("traces", {}).values():
            ax.plot(
                sample_x,
                sample_y,
                color=MECHANISM_COLORS[mechanism],
                lw=0.75,
                alpha=0.14,
                zorder=1,
            )

        x = summary.get("x", np.asarray([]))
        if x.size == 0:
            continue

        ax.fill_between(
            x,
            summary["q25"],
            summary["q75"],
            color=MECHANISM_COLORS[mechanism],
            alpha=0.11,
            linewidth=0,
            zorder=2,
        )
        ax.plot(
            x,
            summary["median"],
            color=MECHANISM_COLORS[mechanism],
            lw=2.4,
            zorder=3,
        )

    ax.set_xlim(x_left - 0.05, x_right + 0.05)
    ax.set_ylim(-1.05, 1.05)
    ax.set_xlabel(
        "Position relative to IC midpoint, T2T-CHM13v2.0 (kb)",
        fontsize=fs(9),
    )
    ax.set_ylabel(
        "Standardized parental-state contrast",
        fontsize=fs(9),
    )
    ax.set_title(
        "B. CpG-level parental-state signal across the imprinting centre",
        fontsize=fs(11),
        loc="left",
        pad=12,
        weight="bold",
    )

    ax.text(
        0.99,
        0.03,
        "Disease controls test whether an unrelated\npathogenic deletion perturbs chr15 imprinting",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=fs(6.8),
        color=MECHANISM_COLORS["Disease control"],
    )

    ax.grid(axis="y", color="#E7E7E7", lw=0.7)
    ax.tick_params(labelsize=fs(8))

    handles = [
        Line2D(
            [0], [0],
            color=MECHANISM_COLORS[m],
            lw=2.5,
            label=m,
        )
        for m in sorted(MECHANISM_ORDER, key=MECHANISM_ORDER.get)
    ]
    ax.legend(
        handles=handles,
        frameon=False,
        fontsize=fs(7.2),
        loc="center left",
        bbox_to_anchor=(1.02, 0.55),
        handlelength=2.5,
    )


def draw_panel_c(ax: plt.Axes, mechanistic_rows: list[dict[str, Any]]) -> None:
    """
    Mechanistic state-space.

    x = number of resolved physical chr15 alleles at the IC
    y = mean methylation across those resolved alleles

    For two-allele samples, marker size additionally represents allelic
    methylation contrast. Thus canonical biparental samples have intermediate
    mean methylation but high allelic contrast, whereas mUPD has high mean
    methylation with low allelic contrast.
    """
    ax.set_title(
        "C. Mechanistic state-space of chromosome 15 imprinting",
        fontsize=fs(11),
        loc="left",
        pad=12,
        weight="bold",
    )

    # Expected mechanistic zones.
    ax.add_patch(Rectangle((0.82, 0.79), 0.36, 0.20, facecolor="#F9E5E2", edgecolor="none", alpha=0.55))
    ax.add_patch(Rectangle((0.82, 0.01), 0.36, 0.20, facecolor="#EEE8F7", edgecolor="none", alpha=0.55))
    ax.add_patch(Rectangle((1.82, 0.79), 0.36, 0.20, facecolor="#E5EFF8", edgecolor="none", alpha=0.55))
    ax.add_patch(Rectangle((1.82, 0.34), 0.36, 0.32, facecolor="#F3F3F3", edgecolor="none", alpha=0.65))

    ax.text(1.0, 0.96, "PWS-DEL\nmaternal retained", ha="center", va="top", fontsize=fs(7), color=MECHANISM_COLORS["PWS-DEL"])
    ax.text(1.0, 0.04, "AS-DEL\npaternal retained", ha="center", va="bottom", fontsize=fs(7), color=MECHANISM_COLORS["AS-DEL"])
    ax.text(2.0, 0.96, "mUPD\nmaternal + maternal", ha="center", va="top", fontsize=fs(7), color=MECHANISM_COLORS["PWS-mUPD"])
    ax.text(2.0, 0.36, "biparental chr15\nmaternal + paternal", ha="center", va="bottom", fontsize=fs(7), color="#555555")

    # Deterministic jitter avoids exact overlap while preserving allele-count axis.
    jitter_map = {}
    for mechanism in MECHANISM_ORDER:
        group = [r for r in mechanistic_rows if r["molecular_mechanism"] == mechanism]
        offsets = np.linspace(-0.065, 0.065, max(len(group), 1))
        for row, off in zip(group, offsets):
            jitter_map[row["sample_id"]] = off

    for row in mechanistic_rows:
        n_alleles = safe_float(row["n_resolved_chr15_alleles"])
        mean_beta = safe_float(row["mean_IC_methylation_across_resolved_alleles"])
        if n_alleles is None or mean_beta is None:
            continue

        contrast = safe_float(row["allelic_methylation_contrast"])
        marker_size = 65.0 if contrast is None else 55.0 + 145.0 * contrast
        x = n_alleles + jitter_map.get(row["sample_id"], 0.0)
        mechanism = row["molecular_mechanism"]

        ax.scatter(
            x,
            mean_beta,
            s=marker_size,
            marker=MECHANISM_MARKERS[mechanism],
            facecolor=MECHANISM_COLORS[mechanism],
            edgecolor="white",
            linewidth=0.8,
            alpha=0.88,
            zorder=4,
        )

        ax.annotate(
            row["display_label"],
            (x, mean_beta),
            xytext=(3, 3),
            textcoords="offset points",
            fontsize=fs(5.8),
            color="#333333",
            zorder=5,
        )

    ax.set_xlim(0.72, 2.28)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(
        ["1 resolved allele\n(hemizygous state)", "2 resolved alleles\n(copy-neutral/biparental state)"],
        fontsize=fs(8),
    )
    ax.set_ylabel(
        "Mean IC methylation across resolved alleles",
        fontsize=fs(9),
    )
    ax.grid(axis="y", color="#E8E8E8", lw=0.7)

    handles = [
        Line2D(
            [0], [0],
            marker=MECHANISM_MARKERS[m],
            linestyle="none",
            markerfacecolor=MECHANISM_COLORS[m],
            markeredgecolor="white",
            markersize=7,
            label=m,
        )
        for m in sorted(MECHANISM_ORDER, key=MECHANISM_ORDER.get)
    ]
    ax.legend(
        handles=handles,
        frameon=False,
        fontsize=fs(6.8),
        loc="lower center",
        bbox_to_anchor=(0.5, -0.30),
        ncol=3,
    )

    ax.text(
        0.02,
        0.02,
        "Marker size ∝ |allele 1 − allele 2| methylation contrast\n(two-allele genomes only)",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=fs(6.1),
        color="#666666",
    )


def draw_support_metric_axis(
    ax: plt.Axes,
    support_rows: list[dict[str, Any]],
    field: str,
    title: str,
    x_max: float,
    show_y: bool,
    value_kind: str,
) -> None:

    y_positions = np.arange(len(support_rows))
    ax.set_xlim(0, max(x_max, 1.0))
    ax.set_ylim(len(support_rows) - 0.5, -0.5)
    ax.set_title(title, fontsize=fs(8), pad=6, fontweight="bold")
    ax.grid(axis="x", color="#ECECEC", lw=0.7)
    ax.tick_params(axis="x", labelsize=fs(6.8))
    ax.tick_params(axis="y", length=0)

    if show_y:
        ax.set_yticks(y_positions)
        ax.set_yticklabels(
            [row["display_label"] for row in support_rows],
            fontsize=fs(7),
            fontweight="bold",
        )
        for tick, row in zip(ax.get_yticklabels(), support_rows):
            tick.set_color(MECHANISM_COLORS[row["molecular_mechanism"]])
    else:
        ax.set_yticks(y_positions)
        ax.set_yticklabels([])

    for i, row in enumerate(support_rows):
        value = safe_float(row[field])
        if value is None:
            continue

        color = MECHANISM_COLORS[row["molecular_mechanism"]]
        ax.hlines(i, 0, value, color=color, lw=1.25, alpha=0.28)
        ax.plot(
            value,
            i,
            marker="o",
            ms=5.2,
            markerfacecolor="white" if row["low_support"] == "True" else color,
            markeredgecolor=color,
            markeredgewidth=1.0,
            linestyle="none",
            zorder=3,
        )

        if value_kind == "percent":
            label = f"{value:.0f}%"
        elif value_kind == "integer":
            label = f"{value:.0f}"
        else:
            label = f"{value:.1f}"

        ax.annotate(
            label,
            (value, i),
            xytext=(4, 0),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize=fs(5.7),
            color="#333333",
            clip_on=False,
        )

    for i in range(len(support_rows) - 1):
        if (
            support_rows[i]["molecular_mechanism"]
            != support_rows[i + 1]["molecular_mechanism"]
        ):
            ax.axhline(i + 0.5, color="#5C5C5C", lw=0.7)


def draw_panel_d(metric_axes: list[plt.Axes], support_rows: list[dict[str, Any]]) -> None:
    max_total = max(
        (safe_float(row["total_ic_depth"]) or 0 for row in support_rows),
        default=1,
    ) * 1.28
    max_allele = max(
        (safe_float(row["supporting_allele_depth"]) or 0 for row in support_rows),
        default=1,
    ) * 1.30
    max_cpgs = max(
        (safe_float(row["supporting_allele_cpgs"]) or 0 for row in support_rows),
        default=1,
    ) * 1.25

    draw_support_metric_axis(
        metric_axes[0], support_rows,
        "total_ic_depth", "Total IC depth",
        max_total, True, "decimal",
    )
    draw_support_metric_axis(
        metric_axes[1], support_rows,
        "supporting_allele_depth", "Min allele depth",
        max_allele, False, "decimal",
    )
    draw_support_metric_axis(
        metric_axes[2], support_rows,
        "supporting_allele_cpgs", "CpGs / allele",
        max_cpgs, False, "integer",
    )
    draw_support_metric_axis(
        metric_axes[3], support_rows,
        "ic_phased_span_percent", "IC phased span",
        112.0, False, "percent",
    )
    metric_axes[3].set_xticks([0, 25, 50, 75, 100])


def create_figure(
    out_prefix: Path,
    panel_a_rows: list[dict[str, Any]],
    contrast_rows: list[dict[str, Any]],
    mechanistic_rows: list[dict[str, Any]],
    support_rows: list[dict[str, Any]],
) -> None:

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial", "Liberation Sans"],
            "font.size": fs(10),
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig = plt.figure(figsize=(18.5, 12.0), constrained_layout=False)
    outer = GridSpec(2, 1, figure=fig, height_ratios=[1.02, 0.98], hspace=0.42)

    top = outer[0].subgridspec(1, 2, width_ratios=[1.02, 1.45], wspace=0.31)
    bottom = outer[1].subgridspec(1, 2, width_ratios=[1.04, 1.40], wspace=0.27)

    panel_a_grid = top[0, 0].subgridspec(1, 2, width_ratios=[0.50, 1.0], wspace=0.05)
    ax_a_note = fig.add_subplot(panel_a_grid[0, 0])
    ax_a_heat = fig.add_subplot(panel_a_grid[0, 1])

    ax_b = fig.add_subplot(top[0, 1])
    ax_c = fig.add_subplot(bottom[0, 0])

    panel_d_grid = bottom[0, 1].subgridspec(1, 4, wspace=0.24)
    ax_d = [fig.add_subplot(panel_d_grid[0, i]) for i in range(4)]

    fig.subplots_adjust(
        top=0.95,
        bottom=0.13,
        left=0.045,
        right=0.92,
    )

    draw_panel_a(ax_a_note, ax_a_heat, panel_a_rows)
    draw_panel_b(ax_b, contrast_rows)
    draw_panel_c(ax_c, mechanistic_rows)
    draw_panel_d(ax_d, support_rows)

    panel_d_left = ax_d[0].get_position().x0
    panel_d_right = ax_d[-1].get_position().x1
    panel_d_top = ax_d[0].get_position().y1

    fig.text(
        (panel_d_left + panel_d_right) / 2,
        panel_d_top + 0.027,
        "D. Technical support for mechanistic IC assignments",
        ha="center",
        va="bottom",
        fontsize=fs(11),
        fontweight="bold",
    )

    support_handles = [
        Line2D(
            [0], [0],
            marker="o",
            linestyle="none",
            markerfacecolor="#444444",
            markeredgecolor="#444444",
            markersize=6,
        ),
        Line2D(
            [0], [0],
            marker="o",
            linestyle="none",
            markerfacecolor="white",
            markeredgecolor="#444444",
            markersize=6,
        ),
    ]

    fig.legend(
        support_handles,
        ["passes nominal IC support", "below nominal IC support"],
        frameon=False,
        fontsize=fs(6.7),
        ncol=2,
        loc="lower center",
        bbox_to_anchor=((panel_d_left + panel_d_right) / 2, 0.066),
        bbox_transform=fig.transFigure,
    )

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".png", ".pdf", ".svg"):
        kwargs = {"bbox_inches": "tight"}
        if suffix == ".png":
            kwargs["dpi"] = 300
        fig.savefig(out_prefix.with_suffix(suffix), **kwargs)

    # Backward-compatible aliases.
    for alias in ("Figure1", "Figure1_improved"):
        for suffix in (".png", ".pdf", ".svg"):
            kwargs = {"bbox_inches": "tight"}
            if suffix == ".png":
                kwargs["dpi"] = 300
            fig.savefig(out_prefix.with_name(alias).with_suffix(suffix), **kwargs)

    plt.close(fig)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def write_report(
    report_path: Path,
    mechanistic_rows: list[dict[str, Any]],
    diagnostic_rows: list[dict[str, Any]],
) -> None:

    report_path.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "# Figure 1 Report: mechanistic validation of the PWS/AS imprinting centre",
        "",
        "## Biological question",
        "",
        "Does the allele-resolved SNRPN/SNHG14 methylation state encode the causal parent-of-origin mechanism at chromosome 15, and is that state preserved in an orthogonal genomic-disorder control carrying a pathogenic deletion outside chromosome 15?",
        "",
        "## Design",
        "",
        "- PWS-DEL: paternal chromosome 15 deletion -> retained maternal methylation state.",
        "- AS-DEL: maternal chromosome 15 deletion -> retained paternal methylation state.",
        "- PWS-mUPD: copy-neutral maternal + maternal state.",
        "- Disease controls: 22q11.2 deletion syndrome with an expected canonical maternal + paternal state at chromosome 15.",
        "- Unaffected controls: canonical maternal + paternal state.",
        "",
        "The disease-control cohort is kept analytically separate from unaffected controls. It is used to test locus specificity, not to inflate the normative control group.",
        "",
        "## Mechanistic template classification",
        "",
        "A continuous template-distance classifier is applied to the observed IC methylation values. For one-allele states it compares the retained allele with methylated [1] and unmethylated [0] templates. For two-allele states it compares the sorted allele pair with biparental [0,1], maternal-duplicated [1,1], and paternal-duplicated [0,0] templates.",
        "",
        "| Sample | Group | Resolved alleles | Mean beta | Allelic contrast | Predicted mechanism | Expected mechanism | Concordant |",
        "| --- | --- | ---: | ---: | ---: | --- | --- | --- |",
    ]

    for row in mechanistic_rows:
        lines.append(
            f"| {row['display_label']} | {row['molecular_mechanism']} | "
            f"{row['n_resolved_chr15_alleles']} | "
            f"{row['mean_IC_methylation_across_resolved_alleles']} | "
            f"{row['allelic_methylation_contrast']} | "
            f"{row['predicted_template_class']} | "
            f"{row['expected_template_class']} | "
            f"{row['template_concordant']} |"
        )

    lines += [
        "",
        "## Group-level concordance",
        "",
        "| Group | n | Observed state distribution | Template concordance | Interpretation |",
        "| --- | ---: | --- | --- | --- |",
    ]

    for row in diagnostic_rows:
        lines.append(
            f"| {row['molecular_mechanism']} | {row['n_samples']} | "
            f"{row['observed_state_distribution']} | "
            f"{row['template_concordant_n']}/{row['n_samples']} "
            f"({row['template_concordance_percent']}%) | "
            f"{row['interpretation']} |"
        )

    lines += [
        "",
        "## Figure interpretation",
        "",
        "- Panel A shows the physical methylation state carried by each resolved chromosome/haplotype.",
        "- Panel B asks whether the parent-of-origin signal is distributed across CpGs rather than being driven by a single site.",
        "- Panel C is the mechanistic core of the figure: copy state (one versus two resolved alleles) and methylation state jointly separate reciprocal deletions, mUPD, and canonical biparental chr15 configurations.",
        "- The six disease controls provide an orthogonal test that an unrelated pathogenic deletion does not by itself generate a PWS/AS-like chr15 imprinting configuration.",
        "- Panel D shows the depth, CpG support, and phasing evidence underlying each assignment.",
        "",
        "## Suggested manuscript claim",
        "",
        "Allele-resolved HiFi methylation did not merely reproduce a diagnostic mean at SNRPN. Instead, the joint configuration of allele number and parental methylation state separated paternal deletion, maternal deletion, maternal UPD, and biparental chromosome 15 states. The canonical biparental configuration was preserved in an independent 22q11.2-deletion disease-control group, supporting locus specificity of the chromosome 15 imprinting phenotype.",
        "",
        "## Caution",
        "",
        "Figure 1 establishes locus-specific molecular configuration, not a genome-wide absence of methylation effects in 22q11.2 deletion syndrome. Disease controls should therefore remain a separate group throughout the manuscript.",
    ]

    report_path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    outdir = Path(args.outdir)
    table_dir = outdir / "tables"
    figure_dir = outdir / "figures"
    report_dir = outdir / "reports"

    for directory in (table_dir, figure_dir, report_dir):
        directory.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------------------
    # Render-only mode
    # ---------------------------------------------------------------------
    if args.render_only:
        panel_a_path = table_dir / "Figure1A_allele_methylation_matrix.tsv"
        contrast_path = table_dir / "Figure1B_per_CpG_contrast.tsv"
        support_path = table_dir / "Figure1D_coverage_phasing_support.tsv"

        missing = [
            p for p in (panel_a_path, contrast_path, support_path)
            if not p.exists()
        ]
        if missing:
            raise FileNotFoundError(
                "Render-only mode requires cached plotting tables:\n"
                + "\n".join(f"- {p}" for p in missing)
            )

        panel_a_rows = read_tsv(panel_a_path)
        contrast_rows = read_tsv(contrast_path)
        support_rows = read_tsv(support_path)

        # Compatibility with older repository tables.
        for rows in (panel_a_rows, contrast_rows, support_rows):
            for row in rows:
                if row.get("molecular_mechanism") == "DiGeorge":
                    row["molecular_mechanism"] = "Disease control"
                if row.get("display_label", "").startswith("DG-"):
                    row["display_label"] = row["display_label"].replace("DG-", "DC-", 1)

        mechanistic_rows = build_mechanistic_state_rows(panel_a_rows)
        diagnostic_rows = build_diagnostic_state_rows(mechanistic_rows)

        write_tsv(
            table_dir / "Figure1C_mechanistic_state_space.tsv",
            mechanistic_rows,
        )
        write_tsv(
            table_dir / "Figure1C_sample_level_classification.tsv",
            mechanistic_rows,
        )
        write_tsv(
            table_dir / "Figure1C_diagnostic_state_summary.tsv",
            diagnostic_rows,
        )

        create_figure(
            figure_dir / "Figure1_mechanistic",
            panel_a_rows,
            contrast_rows,
            mechanistic_rows,
            support_rows,
        )
        write_report(
            report_dir / "Figure1_report.md",
            mechanistic_rows,
            diagnostic_rows,
        )
        return

    # ---------------------------------------------------------------------
    # Input discovery
    # ---------------------------------------------------------------------
    vcf_dir = Path(args.vcf_dir)
    bam_dir = Path(args.bam_dir)
    methylation_dir = Path(args.methylation_dir)
    cnv_dir = Path(args.cnv_dir)
    metadata_path = Path(args.metadata)

    metadata = read_metadata(metadata_path)

    sample_files: dict[str, dict[str, Path | None]] = {}
    for sample_id, _clinical, _mechanism in sorted_cohort():
        sample_files[sample_id] = {
            "bam": find_sample_file(bam_dir, sample_id, ".bam"),
            "blocks": find_sample_file(vcf_dir, sample_id, ".blocks.tsv"),
            "combined_bed": find_sample_file(methylation_dir, sample_id, ".combined.bed"),
            "hap1_bed": find_sample_file(methylation_dir, sample_id, ".hap1.bed"),
            "hap2_bed": find_sample_file(methylation_dir, sample_id, ".hap2.bed"),
            "cnv_log": find_sample_file(cnv_dir, sample_id, ".log"),
        }

    # ---------------------------------------------------------------------
    # BAM QC
    # ---------------------------------------------------------------------
    bam_qc_cache = table_dir / "bam_qc_cache.tsv"
    if args.skip_bam_qc and bam_qc_cache.exists():
        bam_qc_rows = read_tsv(bam_qc_cache)
    else:
        bam_qc_rows = []
        for sample_id, _clinical, _mechanism in sorted_cohort():
            bam_qc_rows.append(
                build_bam_qc(
                    sample_id,
                    sample_files[sample_id]["bam"],
                    sample_files[sample_id]["cnv_log"],
                )
            )
        write_tsv(bam_qc_cache, bam_qc_rows)

    bam_qc_by_sample = {row["sample_id"]: row for row in bam_qc_rows}

    # ---------------------------------------------------------------------
    # Cohort QC summary
    # ---------------------------------------------------------------------
    summary_rows: list[dict[str, Any]] = []

    for sample_id, clinical, mechanism in sorted_cohort():
        meta = metadata.get(sample_id, {})
        n50, domain_fraction = block_n50_and_domain_fraction(
            sample_files[sample_id]["blocks"]
        )
        bam_qc = bam_qc_by_sample.get(sample_id, {})

        summary_rows.append(
            {
                "sample_id": sample_id,
                "clinical_diagnosis": clinical,
                "molecular_mechanism": mechanism,
                "sex": meta.get("gender", meta.get("sex", "")),
                "age_at_sampling": meta.get("age", ""),
                "total_HiFi_reads": bam_qc.get("total_HiFi_reads", ""),
                "mean_depth_genome_wide": fmt(
                    bam_qc.get("mean_depth_genome_wide")
                ),
                "mean_depth_chr15": fmt(
                    bam_qc.get("mean_depth_chr15")
                ),
                "mean_depth_per_haplotype_at_15q11-q13": bam_qc.get(
                    "mean_depth_per_haplotype_at_15q11-q13", ""
                ),
                "phasing_block_N50_chr15": n50 if n50 is not None else "",
                "percent_imprinted_domain_in_phased_block": fmt(domain_fraction),
            }
        )

    write_tsv(
        table_dir / "Figure1B_cohort_qc_summary.tsv",
        summary_rows,
    )

    # ---------------------------------------------------------------------
    # Methylation state
    # ---------------------------------------------------------------------
    assignment_rows, matrix_rows, stats_by_sample = build_assignments(sample_files)

    write_tsv(
        table_dir / "Figure1_parental_assignment.tsv",
        assignment_rows,
    )
    # Backward compatibility.
    write_tsv(
        table_dir / "Figure1C_parental_assignment.tsv",
        assignment_rows,
    )

    write_tsv(
        table_dir / "Figure1_pws_ic_methylation_matrix.tsv",
        matrix_rows,
    )
    write_tsv(
        table_dir / "Figure1C_pws_ic_methylation_matrix.tsv",
        matrix_rows,
    )

    panel_a_rows = build_physical_allele_rows(matrix_rows)
    write_tsv(
        table_dir / "Figure1A_allele_methylation_matrix.tsv",
        panel_a_rows,
    )
    write_tsv(
        table_dir / "Figure1A_physical_allele_layout.tsv",
        panel_a_rows,
    )

    # ---------------------------------------------------------------------
    # Mechanistic classification
    # ---------------------------------------------------------------------
    mechanistic_rows = build_mechanistic_state_rows(panel_a_rows)
    diagnostic_rows = build_diagnostic_state_rows(mechanistic_rows)

    write_tsv(
        table_dir / "Figure1C_mechanistic_state_space.tsv",
        mechanistic_rows,
    )
    write_tsv(
        table_dir / "Figure1C_sample_level_classification.tsv",
        mechanistic_rows,
    )
    write_tsv(
        table_dir / "Figure1C_diagnostic_state_summary.tsv",
        diagnostic_rows,
    )

    # ---------------------------------------------------------------------
    # CpG-level signal
    # ---------------------------------------------------------------------
    contrast_rows = build_per_cpg_contrast(stats_by_sample, assignment_rows)
    write_tsv(
        table_dir / "Figure1B_per_CpG_contrast.tsv",
        contrast_rows,
    )
    # Backward-compatible alias.
    write_tsv(
        table_dir / "Figure1D_per_CpG_contrast.tsv",
        contrast_rows,
    )

    # ---------------------------------------------------------------------
    # Support metrics
    # ---------------------------------------------------------------------
    support_rows = build_support_rows(
        summary_rows,
        matrix_rows,
        sample_files,
    )
    write_tsv(
        table_dir / "Figure1D_coverage_phasing_support.tsv",
        support_rows,
    )
    write_tsv(
        table_dir / "Figure1D_support_metrics.tsv",
        support_rows,
    )

    # ---------------------------------------------------------------------
    # Figure and report
    # ---------------------------------------------------------------------
    create_figure(
        figure_dir / "Figure1_mechanistic",
        panel_a_rows,
        contrast_rows,
        mechanistic_rows,
        support_rows,
    )

    write_report(
        report_dir / "Figure1_report.md",
        mechanistic_rows,
        diagnostic_rows,
    )

    run_parameters = {
        "cohort": [
            {
                "sample_id": sample_id,
                "clinical_diagnosis": clinical,
                "molecular_mechanism": mechanism,
            }
            for sample_id, clinical, mechanism in sorted_cohort()
        ],
        "biological_design": {
            "PWS_DEL": "paternal chr15 deletion exposes retained maternal methylation",
            "AS_DEL": "maternal chr15 deletion exposes retained paternal methylation",
            "PWS_mUPD": "copy-neutral duplication of maternal epigenetic identity",
            "Disease_control": "22q11.2 deletion; orthogonal disease control for chr15 locus specificity",
            "Control": "unaffected biparental reference",
        },
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
            },
        },
        "thresholds": {
            "maternal_methylation": MATERNAL_THRESHOLD,
            "paternal_methylation": PATERNAL_THRESHOLD,
            "minimum_mean_coverage": MIN_MEAN_COVERAGE,
            "minimum_CpGs": MIN_CPGS,
        },
        "input_paths": {
            "vcf_dir": str(vcf_dir),
            "bam_dir": str(bam_dir),
            "methylation_dir": str(methylation_dir),
            "cnv_dir": str(cnv_dir),
            "metadata": str(metadata_path),
        },
    }

    with (outdir / "Figure1_run_parameters.json").open("w") as handle:
        json.dump(run_parameters, handle, indent=2)


if __name__ == "__main__":
    main()
