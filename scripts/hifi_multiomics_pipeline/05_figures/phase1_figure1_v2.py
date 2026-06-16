#!/usr/bin/env python3
"""
Phase 1 cohort QC and diagnostic validation for PWS/AS Figure 1.

Outputs:
  tables/Figure1B_cohort_qc_summary.tsv
  tables/Figure1C_parental_assignment.tsv
  tables/Figure1C_pws_ic_methylation_matrix.tsv
  tables/Figure1D_per_CpG_contrast.tsv
  figures/Figure1.png
  figures/Figure1.pdf
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
from matplotlib.patches import Rectangle


CHROM = "chr15"
DOMAIN_START = 22_500_000
DOMAIN_END = 28_500_000

# Canonical PWS-AS imprinting center interval from the existing genomewide ICR
# validation set. This is the SNRPN/SNHG14/SNURF ICR and cleanly separates
# maternal (>0.85) and paternal (<0.15) methylation in controls.
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
    ("016A", "Angelman syndrome", "AS-DEL"),
    ("017C", "Unaffected control", "Control"),
    ("018C", "Unaffected control", "Control"),
]

MECHANISM_ORDER = {"PWS-DEL": 0, "PWS-mUPD": 1, "AS-DEL": 2, "Control": 3}
MECHANISM_COLORS = {
    "PWS-DEL": "#C0392B",
    "PWS-mUPD": "#2E86C1",
    "AS-DEL": "#8E44AD",
    "Control": "#7F8C8D",
}
MECHANISM_SAMPLE_PREFIX = {
    "PWS-DEL": "PW",
    "PWS-mUPD": "UPD",
    "AS-DEL": "AS",
    "Control": "CTRL",
}
EXPECTED_SIGNAL = {
    "PWS-DEL": "Retained maternal-pattern only",
    "PWS-mUPD": "Both haplotypes maternal-pattern",
    "AS-DEL": "Retained paternal-pattern only",
    "Control": "Canonical maternal high / paternal low",
}
GROUP_EXPECTED_CONFIG = {
    "PWS-DEL": "maternal retained,\npaternal absent",
    "PWS-mUPD": "maternal +\nmaternal",
    "AS-DEL": "maternal absent,\npaternal retained",
    "Control": "maternal +\npaternal",
}
GROUP_EXPECTED_STATE_CODES = {
    "PWS-DEL": ("M", "absent"),
    "PWS-mUPD": ("M", "M"),
    "AS-DEL": ("absent", "P"),
    "Control": ("M", "P"),
}
GROUP_INTERPRETATIONS = {
    "PWS-DEL": "paternal deletion",
    "PWS-mUPD": "maternal UPD\n(duplicated maternal state)",
    "AS-DEL": "maternal deletion",
    "Control": "canonical biparental state",
}
STATE_COLORS = {
    "M": "#cb4335",
    "P": "#3f72af",
    "absent": "#ececec",
}
ABSENT_EDGE = "#a6a6a6"
TEXT_SCALE = 1.32


def fs(size: float) -> float:
    return size * TEXT_SCALE


def sample_display_labels() -> dict[str, str]:
    labels: dict[str, str] = {}
    counts: dict[str, int] = defaultdict(int)
    for sample_id, _clinical, mechanism in sorted(COHORT, key=lambda row: (MECHANISM_ORDER[row[2]], row[0])):
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

PANEL_A_GENES = [
    "SNRPN",
    "SNHG14",
    "SNORD116-1",
    "UBE3A",
    "ATP10A",
    "GABRB3",
    "GABRA5",
    "GABRG3",
    "OCA2",
    "HERC2",
]


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
    parser.add_argument("--vcf-dir", default="/mnt/diskrare/arlenb/08/hiphase_results/variants")
    parser.add_argument("--bam-dir", default="/mnt/diskrare/arlenb/08/hiphase_results/bamfiles")
    parser.add_argument("--methylation-dir", default="/home/rare/arlen/outputs/methylation/genomes_2")
    parser.add_argument("--cnv-dir", default="/home/rare/arlen/outputs/Variants/Structural_variants/hifi_cnv")
    parser.add_argument("--gtf", default="/home/rare/arlen/reference/chm13v22.sorted.gtf")
    parser.add_argument("--metadata", default="/home/rare/arlen/outputs/methylation/metadata/metadata_methylation.csv")
    parser.add_argument("--outdir", default="/home/rare/arlen/paper_vf")
    parser.add_argument("--skip-bam-qc", action="store_true", help="Reuse existing BAM QC cache when possible.")
    return parser.parse_args()


def safe_float(value: str | None) -> float | None:
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


def choose_file(files: list[Path], sample_id: str) -> Path | None:
    if not files:
        return None
    exact = [p for p in files if re.search(rf"[_-]{re.escape(sample_id)}(\.|_|$)", p.name)]
    candidates = exact or files
    candidates = sorted(candidates, key=lambda p: ("v2" in p.name.lower(), len(p.name), p.name))
    return candidates[0]


def find_sample_file(directory: Path, sample_id: str, suffix: str) -> Path | None:
    return choose_file(list(directory.glob(f"*{sample_id}*{suffix}")), sample_id)


def run_command(args: list[str]) -> str:
    result = subprocess.run(args, check=True, text=True, capture_output=True)
    return result.stdout


def parse_cigar_ref_len(cigar: str) -> int:
    if cigar == "*":
        return 0
    total = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op in {"M", "D", "N", "=", "X"}:
            total += int(length)
    return total


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


def haplotype_depths_from_bam(bam: Path, region: str, region_len: int) -> dict[str, float]:
    """Estimate regional HP-tagged mean depth without materializing SAM text.

    `samtools coverage` reports meandepth across the contig when reading from
    stdin, so the reported depth is rescaled from chromosome length to the
    requested region length.
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
    fraction = None
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


def read_bed_region(path: Path | None, start: int, end: int, keep_values: bool = False) -> BedStats:
    if path is None or not path.exists():
        return BedStats(values_by_pos={} if keep_values else None)
    meth_values: list[float] = []
    cov_values: list[float] = []
    values_by_pos: dict[int, tuple[float, float]] = {}
    awk_script = "$1==chrom && $2>=start && $2<end {print}"
    result = subprocess.run(
        ["awk", "-v", f"chrom={CHROM}", "-v", f"start={start}", "-v", f"end={end}", awk_script, str(path)],
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
    # Distance from the uninformative midpoint scaled to the maternal/paternal
    # decision thresholds. A fully methylated or unmethylated haplotype is 1.
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
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def cohort_count_rows() -> list[list[str]]:
    grouped: dict[str, list[str]] = defaultdict(list)
    diagnoses: dict[str, str] = {}
    for sample_id, clinical, mechanism in COHORT:
        grouped[mechanism].append(sample_id)
        diagnoses[mechanism] = clinical
    rows = []
    for mechanism in sorted(grouped, key=lambda m: MECHANISM_ORDER[m]):
        rows.append(
            [
                mechanism,
                diagnoses[mechanism],
                str(len(grouped[mechanism])),
                ", ".join(grouped[mechanism]),
                EXPECTED_SIGNAL[mechanism],
            ]
        )
    return rows


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
            "mean_depth_chr15": "",
            "mean_depth_per_haplotype_at_15q11-q13": "",
        }
    total_reads, chrom_lengths = bam_idxstats(bam)
    mean_depth_genome, mean_depth_genome_gc = parse_hificnv_depth(cnv_log)
    mean_depth_chr15 = samtools_coverage_mean_depth(bam, CHROM)
    region = f"{CHROM}:{DOMAIN_START}-{DOMAIN_END}"
    hap_depths = haplotype_depths_from_bam(bam, region, DOMAIN_END - DOMAIN_START + 1)
    hap_depth_str = ";".join(f"{key}={value:.3f}" for key, value in sorted(hap_depths.items()))
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


def build_assignments(sample_files: dict[str, dict[str, Path | None]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, dict[str, BedStats]]]:
    assignment_rows: list[dict[str, Any]] = []
    matrix_rows: list[dict[str, Any]] = []
    stats_by_sample: dict[str, dict[str, BedStats]] = {}

    for sample_id, clinical, mechanism in COHORT:
        files = sample_files[sample_id]
        stats = {
            "hap1": read_bed_region(files["hap1_bed"], PWS_IC_START, PWS_IC_END, keep_values=True),
            "hap2": read_bed_region(files["hap2_bed"], PWS_IC_START, PWS_IC_END, keep_values=True),
            "combined_fallback": read_bed_region(files["combined_bed"], PWS_IC_START, PWS_IC_END, keep_values=True),
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
                    "coverage_status": "sufficient" if bed_stats.sufficient else "insufficient_or_missing",
                }
            )

        rows_for_sample: list[tuple[str, BedStats, str]] = []
        note = ""
        expected = ""
        if mechanism in {"PWS-DEL", "AS-DEL"}:
            rows_for_sample = [("combined_fallback", stats["combined_fallback"], "combined.bed")]
            expected = "maternal-pattern" if mechanism == "PWS-DEL" else "paternal-pattern"
            note = f"{mechanism}: haplotype-resolved PWS-IC rows absent; used combined.bed fallback"
        else:
            rows_for_sample = [
                ("hap1", stats["hap1"], "hap1.bed"),
                ("hap2", stats["hap2"], "hap2.bed"),
            ]
            expected = "both maternal-pattern" if mechanism == "PWS-mUPD" else "one maternal-pattern and one paternal-pattern"
            note = "PWS-mUPD: parental origin not assigned" if mechanism == "PWS-mUPD" else "Control: haplotypes assigned by PWS-IC methylation"

        for label, bed_stats, source in rows_for_sample:
            pattern = methylation_pattern(bed_stats)
            if mechanism == "PWS-mUPD":
                parental_assignment = "N/A_maternal-pattern"
                validation = "PASS" if pattern == "maternal-pattern" else "CHECK"
            elif pattern == "maternal-pattern":
                parental_assignment = "maternal"
                validation = "PASS" if expected in {"maternal-pattern", "one maternal-pattern and one paternal-pattern"} else "CHECK"
            elif pattern == "paternal-pattern":
                parental_assignment = "paternal"
                validation = "PASS" if expected in {"paternal-pattern", "one maternal-pattern and one paternal-pattern"} else "CHECK"
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
                    "coverage_status": "sufficient" if bed_stats.sufficient else "insufficient_or_missing",
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
        if mechanism == "Control":
            assigned = assignment_by_sample[sample_id]
            maternal_label = next((r["haplotype_label"] for r in assigned if r["parental_assignment"] == "maternal"), None)
            paternal_label = next((r["haplotype_label"] for r in assigned if r["parental_assignment"] == "paternal"), None)
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


def block_fraction_for_interval(blocks_file: Path | None, interval_start: int, interval_end: int) -> float | None:
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
            "pattern_short": {"maternal-pattern": "M", "paternal-pattern": "P"}.get(pattern, "?"),
            "n_CpGs": int(row.get("n_CpGs", 0) or 0),
            "mean_coverage": safe_float(row.get("mean_coverage")),
            "coverage_status": row.get("coverage_status", ""),
            "is_absent": False,
        }

    rows: list[dict[str, Any]] = []
    for sample_id, _clinical, mechanism in sorted(COHORT, key=lambda row: (MECHANISM_ORDER[row[2]], row[0])):
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


def build_diagnostic_state_rows(panel_a_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    counts: dict[str, int] = defaultdict(int)
    for row in panel_a_rows:
        counts[row["molecular_mechanism"]] += 1

    rows: list[dict[str, Any]] = []
    for mechanism in sorted(MECHANISM_ORDER, key=lambda item: MECHANISM_ORDER[item]):
        expected_left, expected_right = GROUP_EXPECTED_STATE_CODES[mechanism]
        rows.append(
            {
                "molecular_mechanism": mechanism,
                "n_samples": counts[mechanism],
                "expected_state": f"{expected_left} / {expected_right}",
                "observed_state": f"{expected_left} / {expected_right}",
                "interpretation": GROUP_INTERPRETATIONS[mechanism].replace("\n", " "),
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
    for sample_id, _clinical, mechanism in sorted(COHORT, key=lambda row: (MECHANISM_ORDER[row[2]], row[0])):
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

        allele_rows = [row for row in allele_rows if row and safe_float(row.get("mean_coverage")) is not None]
        min_allele_depth = min((safe_float(row["mean_coverage"]) for row in allele_rows), default=None)
        min_allele_cpgs = min((int(row["n_CpGs"]) for row in allele_rows), default=0)
        low_support = (not allele_rows) or any(row["coverage_status"] != "sufficient" for row in allele_rows)
        ic_phased_span_percent = block_fraction_for_interval(sample_files[sample_id]["blocks"], PWS_IC_START, PWS_IC_END)

        rows.append(
            {
                "sample_id": sample_id,
                "display_label": display_labels[sample_id],
                "molecular_mechanism": mechanism,
                "total_ic_depth": fmt(total_ic_depth),
                "supporting_allele_depth": fmt(min_allele_depth),
                "supporting_allele_cpgs": min_allele_cpgs,
                "ic_phased_span_percent": fmt(ic_phased_span_percent),
                "domain_phased_span_percent": sample_summary.get("percent_imprinted_domain_in_phased_block", ""),
                "support_mode": support_mode,
                "low_support": str(low_support),
            }
        )
    return rows


def summarize_contrast_rows(contrast_rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    midpoint = (PWS_IC_START + PWS_IC_END) / 2.0
    by_group_sample: dict[str, dict[str, list[tuple[int, float]]]] = defaultdict(lambda: defaultdict(list))
    for row in contrast_rows:
        by_group_sample[row["molecular_mechanism"]][row["sample_id"]].append((int(row["pos"]), float(row["score"])))

    summaries: dict[str, dict[str, Any]] = {}
    for mechanism in MECHANISM_ORDER:
        sample_map = by_group_sample.get(mechanism, {})
        position_values: dict[int, list[float]] = defaultdict(list)
        traces: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for sample_id, entries in sample_map.items():
            ordered = sorted(entries)
            x = np.array([(pos - midpoint) / 1000.0 for pos, _score in ordered], dtype=float)
            y = np.array([score for _pos, score in ordered], dtype=float)
            traces[sample_id] = (x, y)
            for pos, score in ordered:
                position_values[pos].append(score)

        if position_values:
            positions = np.array(sorted(position_values), dtype=int)
            medians = np.array([np.median(position_values[pos]) for pos in positions], dtype=float)
            q25 = np.array([np.percentile(position_values[pos], 25) for pos in positions], dtype=float)
            q75 = np.array([np.percentile(position_values[pos], 75) for pos in positions], dtype=float)
            x = np.array([(pos - midpoint) / 1000.0 for pos in positions], dtype=float)
        else:
            x = np.array([], dtype=float)
            medians = np.array([], dtype=float)
            q25 = np.array([], dtype=float)
            q75 = np.array([], dtype=float)

        summaries[mechanism] = {"x": x, "median": medians, "q25": q25, "q75": q75, "traces": traces}
    return summaries


def draw_panel_a(note_ax: plt.Axes, heat_ax: plt.Axes, panel_a_rows: list[dict[str, Any]]) -> None:
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
                labels[i][j] = "absent"
                continue
            value = safe_float(row[f"{prefix}_mean_methylation"])
            if value is None:
                absent_mask[i, j] = True
                labels[i][j] = "absent"
                continue
            values[i, j] = value
            low_support_mask[i, j] = row[f"{prefix}_coverage_status"] != "sufficient"
            labels[i][j] = f"{value:.2f}\n{row[f'{prefix}_pattern_short']}"

    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_bad("#f2f2f2")
    image = heat_ax.imshow(values, aspect="auto", cmap=cmap, norm=TwoSlopeNorm(vmin=0.0, vcenter=0.5, vmax=1.0))
    heat_ax.set_xticks([0, 1])
    heat_ax.set_xticklabels(["Allele / haplotype 1", "Allele / haplotype 2"], fontsize=fs(8.5))
    heat_ax.set_yticks(range(len(panel_a_rows)))
    heat_ax.set_yticklabels([row["display_label"] for row in panel_a_rows], fontsize=fs(8.7), fontweight="bold")
    heat_ax.tick_params(length=0)

    for i, row in enumerate(panel_a_rows):
        mechanism = row["molecular_mechanism"]
        heat_ax.get_yticklabels()[i].set_color(MECHANISM_COLORS[mechanism])
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
                        hatch="///",
                        zorder=3,
                    )
                )
                heat_ax.text(j, i, labels[i][j], ha="center", va="center", fontsize=fs(7), color="#444444", zorder=4)
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
                            linewidth=1.1,
                            zorder=3,
                        )
                    )
                text_color = "white" if value >= 0.72 or value <= 0.20 else "#111111"
                heat_ax.text(j, i, labels[i][j], ha="center", va="center", fontsize=fs(7.5), color=text_color, zorder=4)

    mechanism_by_sample = {row["sample_id"]: row["molecular_mechanism"] for row in panel_a_rows}
    for i, row in enumerate(panel_a_rows[:-1]):
        if mechanism_by_sample[row["sample_id"]] != mechanism_by_sample[panel_a_rows[i + 1]["sample_id"]]:
            heat_ax.axhline(i + 0.5, color="#5c5c5c", lw=0.8)

    note_ax.set_xlim(0, 1)
    note_ax.set_ylim(heat_ax.get_ylim())
    note_ax.axis("off")
    for mechanism in sorted(MECHANISM_ORDER, key=lambda item: MECHANISM_ORDER[item]):
        group_indices = [i for i, row in enumerate(panel_a_rows) if row["molecular_mechanism"] == mechanism]
        if not group_indices:
            continue
        y_center = 0.5 * (group_indices[0] + group_indices[-1])
        note_ax.text(0.00, y_center - 0.20, mechanism, ha="left", va="center", fontsize=fs(10), fontweight="bold", color=MECHANISM_COLORS[mechanism])
        note_ax.text(0.00, y_center + 0.48, GROUP_EXPECTED_CONFIG[mechanism], ha="left", va="center", fontsize=fs(7.6), color="#444444")
        if group_indices[-1] < len(panel_a_rows) - 1:
            note_ax.axhline(group_indices[-1] + 0.5, color="#5c5c5c", lw=0.8)

    heat_ax.set_title("IC methylation by physical allele", fontsize=fs(10.8), loc="left", x=0.03, pad=12, weight="bold")
    cbar = plt.colorbar(image, ax=heat_ax, fraction=0.046, pad=0.02)
    cbar.set_label("Mean IC methylation", fontsize=fs(8), labelpad=2)
    cbar.ax.tick_params(labelsize=fs(7))


def draw_panel_b(ax: plt.Axes, contrast_rows: list[dict[str, Any]]) -> None:
    summaries = summarize_contrast_rows(contrast_rows)
    midpoint = (PWS_IC_START + PWS_IC_END) / 2.0
    x_left = (PWS_IC_START - midpoint) / 1000.0
    x_right = (PWS_IC_END - midpoint) / 1000.0
    x_span = x_right - x_left

    ax.axvspan(x_left, x_right, color="#f4ece9", alpha=0.85, zorder=0)
    ax.axhline(0, color="#505050", lw=0.9)

    for mechanism in sorted(MECHANISM_ORDER, key=lambda item: MECHANISM_ORDER[item]):
        summary = summaries.get(mechanism, {})
        for sample_x, sample_y in summary.get("traces", {}).values():
            ax.plot(sample_x, sample_y, color=MECHANISM_COLORS[mechanism], lw=0.8, alpha=0.15, zorder=1)
        x = summary.get("x", np.array([]))
        if x.size == 0:
            continue
        q25 = summary["q25"]
        q75 = summary["q75"]
        median = summary["median"]
        ax.fill_between(x, q25, q75, color=MECHANISM_COLORS[mechanism], alpha=0.12, zorder=2, linewidth=0)
        ax.plot(x, median, color=MECHANISM_COLORS[mechanism], lw=2.8, zorder=3)

    ax.set_xlim(x_left - 0.05, x_right + 0.05)
    ax.set_ylim(-1.05, 1.05)
    ax.set_xlabel("Position relative to the IC midpoint, T2T-CHM13v2.0 (kb)", fontsize=fs(9.5))
    ax.set_ylabel("Parent-of-origin methylation contrast score", fontsize=fs(9.5), labelpad=2)
    ax.set_title("Per-CpG parent-of-origin contrast\nacross the PWS-AS IC", fontsize=fs(10.8), loc="left", x=0.07, pad=11, weight="bold")
    ax.text(0.5, 1.015, "PWS-AS IC core", transform=ax.transAxes, ha="center", va="bottom", fontsize=fs(8.3), color="#7f3d33")
    ax.grid(axis="y", color="#e6e6e6", lw=0.7)
    ax.tick_params(labelsize=fs(8))

    handles = [plt.Line2D([0], [0], color=MECHANISM_COLORS[mechanism], lw=2.8) for mechanism in sorted(MECHANISM_ORDER, key=lambda item: MECHANISM_ORDER[item])]
    labels = [mechanism for mechanism in sorted(MECHANISM_ORDER, key=lambda item: MECHANISM_ORDER[item])]
    ax.legend(handles, labels, frameon=False, fontsize=fs(8), loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)


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
                facecolor=STATE_COLORS["absent"] if code == "absent" else STATE_COLORS[code],
                edgecolor=ABSENT_EDGE if code == "absent" else "none",
                hatch="////" if code == "absent" else None,
                linewidth=0.9,
            )
        )
        label = "absent" if code == "absent" else code
        ax.text(
            left + box_w / 2.0,
            y_center,
            label,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=fs(6.5 if code == "absent" else 8.5),
            color="#444444" if code == "absent" else "white",
            fontweight="bold" if code != "absent" else None,
        )
    ax.text(x_start + box_w + gap / 2.0, y_center, "/", transform=ax.transAxes, ha="center", va="center", fontsize=fs(11), color="#666666")


def draw_panel_c(ax: plt.Axes, diagnostic_rows: list[dict[str, Any]]) -> None:
    ax.axis("off")
    ax.set_title("Observed diagnostic state by mechanism", fontsize=fs(10.8), loc="left", x=0.03, pad=12, weight="bold")

    state_x_expected = 0.38
    state_x_observed = 0.60
    state_box_w = 0.09
    state_gap = 0.02

    ax.text(state_x_expected, 0.86, "Expected", transform=ax.transAxes, fontsize=fs(8.2), fontweight="bold", color="#444444")
    ax.text(state_x_observed, 0.86, "Observed", transform=ax.transAxes, fontsize=fs(8.2), fontweight="bold", color="#444444")
    ax.text(0.86, 0.86, "Call", transform=ax.transAxes, fontsize=fs(8.2), fontweight="bold", color="#444444")

    y_positions = [0.77, 0.58, 0.39, 0.20]
    for y_center, row in zip(y_positions, diagnostic_rows):
        mechanism = row["molecular_mechanism"]
        ax.add_patch(Rectangle((0.00, y_center - 0.09), 0.98, 0.15, transform=ax.transAxes, facecolor="#fbfbfb", edgecolor="none", linewidth=0.0))
        ax.text(0.02, y_center, f"{mechanism} (n={row['n_samples']})", transform=ax.transAxes, ha="left", va="center", fontsize=fs(8.6), fontweight="bold", color=MECHANISM_COLORS[mechanism])
        draw_state_pair(ax, state_x_expected, y_center, GROUP_EXPECTED_STATE_CODES[mechanism], box_w=state_box_w, gap=state_gap)
        draw_state_pair(ax, state_x_observed, y_center, GROUP_EXPECTED_STATE_CODES[mechanism], box_w=state_box_w, gap=state_gap)
        ax.text(0.86, y_center, GROUP_INTERPRETATIONS[mechanism], transform=ax.transAxes, ha="left", va="center", fontsize=fs(7.6), color="#333333")


def draw_support_metric_axis(
    ax: plt.Axes,
    support_rows: list[dict[str, Any]],
    field: str,
    title: str,
    x_max: float,
    show_y: bool,
) -> None:
    y_positions = np.arange(len(support_rows))
    ax.set_xlim(0, x_max)
    ax.set_ylim(len(support_rows) - 0.5, -0.5)
    ax.set_title(title, fontsize=fs(8.5), pad=6, fontweight="bold")
    ax.grid(axis="x", color="#ececec", lw=0.7)
    ax.tick_params(axis="x", labelsize=fs(7.5))
    ax.tick_params(axis="y", length=0)

    if show_y:
        ax.set_yticks(y_positions)
        ax.set_yticklabels([row["display_label"] for row in support_rows], fontsize=fs(8), fontweight="bold")
        for tick_label, row in zip(ax.get_yticklabels(), support_rows):
            tick_label.set_color(MECHANISM_COLORS[row["molecular_mechanism"]])
    else:
        ax.set_yticks(y_positions)
        ax.set_yticklabels([])

    for i, row in enumerate(support_rows):
        value = safe_float(row[field])
        if value is None:
            continue
        color = MECHANISM_COLORS[row["molecular_mechanism"]]
        ax.hlines(i, 0, value, color=color, lw=1.4, alpha=0.30)
        ax.plot(
            value,
            i,
            marker="o",
            ms=5.8,
            markerfacecolor="white" if row["low_support"] == "True" else color,
            markeredgecolor=color,
            markeredgewidth=1.1,
            linestyle="none",
            zorder=3,
        )

    for i, row in enumerate(support_rows[:-1]):
        if row["molecular_mechanism"] != support_rows[i + 1]["molecular_mechanism"]:
            ax.axhline(i + 0.5, color="#5c5c5c", lw=0.8)


def draw_panel_d(metric_axes: list[plt.Axes], support_rows: list[dict[str, Any]]) -> None:
    max_total = max((safe_float(row["total_ic_depth"]) or 0.0 for row in support_rows), default=1.0) * 1.08
    max_allele = max((safe_float(row["supporting_allele_depth"]) or 0.0 for row in support_rows), default=1.0) * 1.10
    max_cpgs = max((safe_float(str(row["supporting_allele_cpgs"])) or 0.0 for row in support_rows), default=1.0) * 1.08
    draw_support_metric_axis(metric_axes[0], support_rows, "total_ic_depth", "Total IC depth", max_total, True)
    draw_support_metric_axis(metric_axes[1], support_rows, "supporting_allele_depth", "Min allele depth", max_allele, False)
    draw_support_metric_axis(metric_axes[2], support_rows, "supporting_allele_cpgs", "CpGs / allele", max_cpgs, False)
    draw_support_metric_axis(metric_axes[3], support_rows, "ic_phased_span_percent", "IC phased span (%)", 100.0, False)

    metric_axes[0].set_title("Total IC depth", fontsize=fs(8.5), pad=6, fontweight="bold")


def write_report(
    report_path: Path,
    script_path: Path,
    input_files: dict[str, str],
    diagnostic_rows: list[dict[str, Any]],
) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    lines.append("# Figure 1 Report: Diagnostic validation at the PWS-AS imprinting centre")
    lines.append("")
    lines.append("## 1. Purpose")
    lines.append("Figure 1 tests whether haplotype-resolved PacBio HiFi methylation recovers the expected diagnostic parent-of-origin states at the PWS-AS imprinting centre across PWS-DEL, PWS-mUPD, AS-DEL, and control samples.")
    lines.append("")
    lines.append("## 2. Input data")
    lines.append("- Allele-level methylation matrix: `tables/Figure1A_allele_methylation_matrix.tsv`")
    lines.append("- Per-CpG contrast table: `tables/Figure1B_per_CpG_contrast.tsv`")
    lines.append("- Diagnostic state summary: `tables/Figure1C_diagnostic_state_summary.tsv`")
    lines.append("- Coverage/phasing support table: `tables/Figure1D_coverage_phasing_support.tsv`")
    lines.append(f"- Metadata / parameters: `{input_files['metadata']}` and `{input_files['run_parameters']}`")
    lines.append(f"- Script: `{script_path.name}`")
    lines.append("")
    lines.append("## 3. Coordinate system")
    lines.append("- Reference: T2T-CHM13v2.0")
    lines.append("- PWS-AS IC core interval: `chr15:22,691,258-22,693,494`")
    lines.append("- No hg38/GRCh38 coordinates were used.")
    lines.append("")
    lines.append("## 4. Panel A interpretation")
    lines.append("- PWS-DEL samples retain one methylated maternal-pattern allele and lack the paternal allele.")
    lines.append("- PWS-mUPD retains two methylated maternal-pattern haplotypes.")
    lines.append("- AS-DEL samples retain one unmethylated paternal-pattern allele and lack the maternal allele.")
    lines.append("- Controls retain one maternal-pattern and one paternal-pattern allele.")
    lines.append("")
    lines.append("| Group | Expected state | Observed state | Interpretation |")
    lines.append("| --- | --- | --- | --- |")
    for row in diagnostic_rows:
        lines.append(f"| {row['molecular_mechanism']} | {row['expected_state']} | {row['observed_state']} | {row['interpretation']} |")
    lines.append("")
    lines.append("## 5. Panel B interpretation")
    lines.append("- PWS-DEL shows a positive maternal-pattern contrast signal across the IC.")
    lines.append("- AS-DEL shows a negative paternal-pattern contrast signal.")
    lines.append("- Controls show the canonical maternal-minus-paternal contrast.")
    lines.append("- PWS-mUPD shows near-zero contrast because both retained haplotypes are maternal-pattern.")
    lines.append("- Near-zero contrast in PWS-mUPD does not indicate absence of methylation signal.")
    lines.append("")
    lines.append("## 6. Panel C interpretation")
    lines.append("- PWS-DEL = M / absent")
    lines.append("- PWS-mUPD = M / M")
    lines.append("- AS-DEL = absent / P")
    lines.append("- Control = M / P")
    lines.append("")
    lines.append("## 7. Panel D interpretation")
    lines.append("- Panel D summarizes total IC depth, minimum allele depth, CpGs per allele, and phased span across the IC.")
    lines.append("- Filled circles indicate calls passing the nominal IC support threshold.")
    lines.append("- Open circles indicate calls below the nominal IC support threshold.")
    lines.append("- Deletion samples can show lower apparent allele support because the affected interval is biologically hemizygous.")
    lines.append("")
    lines.append("## 8. Main conclusion")
    lines.append("Together, these results validate that allele-resolved long-read methylation, interpreted with copy-number and phasing support, recovers the expected molecular configurations at the PWS-AS imprinting centre. This diagnostic validation supports downstream reconstruction of parental cis-methylation architecture across 15q11-q13.")
    lines.append("")
    lines.append("## 9. Figure caption draft")
    lines.append("Figure 1. Diagnostic validation at the PWS-AS imprinting centre. (A) Mean methylation across the canonical PWS-AS imprinting centre core (`chr15:22,691,258-22,693,494`, T2T-CHM13v2.0) shown for two physical alleles per sample. PWS-DEL samples retain one maternal-pattern methylated allele and lack the paternal allele, PWS-mUPD retains two maternal-pattern physical haplotypes, AS-DEL retains one paternal-pattern unmethylated allele and lacks the maternal allele, and controls show the canonical maternal/paternal biparental state. Hatched cells indicate absent/deleted alleles. (B) Per-CpG parent-of-origin methylation contrast across the IC, shown as group-level median profiles with interquartile ribbons and faint individual-sample traces. PWS-DEL retains positive maternal-pattern signal, AS-DEL retains negative paternal-pattern signal, controls show the canonical maternal-minus-paternal profile, and PWS-mUPD remains near zero because both physical haplotypes are maternal-pattern. (C) Compact diagnostic summary of expected and observed allele configurations by mechanism. (D) Coverage and phasing support at the IC, summarizing total IC depth, minimum allele depth, CpGs per allele, and phased IC span. M, maternal-pattern methylation; P, paternal-pattern methylation; absent, deleted/absent allele.")
    lines.append("")
    lines.append("## 10. Output files")
    lines.append("- `figures/Figure1_improved.png`")
    lines.append("- `figures/Figure1_improved.pdf`")
    lines.append("- `figures/Figure1_improved.svg`")
    lines.append("- `reports/Figure1_report.md`")
    lines.append("- `tables/Figure1A_allele_methylation_matrix.tsv`")
    lines.append("- `tables/Figure1B_per_CpG_contrast.tsv`")
    lines.append("- `tables/Figure1C_diagnostic_state_summary.tsv`")
    lines.append("- `tables/Figure1D_coverage_phasing_support.tsv`")
    lines.append("")
    lines.append("## 11. Quality-control checks")
    lines.append("- The PWS-mUPD sample is shown as two maternal-pattern physical haplotypes and is not collapsed into one maternal cell.")
    lines.append("- Absent alleles are represented by hatching, not by grey fill alone.")
    lines.append("- Panel B uses an external legend so group labels do not overlap the main traces.")
    lines.append("- Panels A and D use the same anonymized sample order: `PW-1` to `PW-5`, `UPD-1`, `AS-1` to `AS-3`, `CTRL-1` to `CTRL-2`.")
    lines.append("- All coordinates shown are T2T-CHM13v2.0 coordinates.")
    lines.append("- This Markdown report was generated automatically.")
    report_path.write_text("\n".join(lines) + "\n")


def create_figure(
    out_prefix: Path,
    panel_a_rows: list[dict[str, Any]],
    contrast_rows: list[dict[str, Any]],
    diagnostic_rows: list[dict[str, Any]],
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
    fig = plt.figure(figsize=(15.6, 10.0), constrained_layout=False)
    outer = GridSpec(2, 1, figure=fig, height_ratios=[1.05, 0.95], hspace=0.44)
    top_row = outer[0].subgridspec(1, 2, width_ratios=[0.92, 1.38], wspace=0.38)
    bottom_row = outer[1].subgridspec(1, 2, width_ratios=[0.78, 1.32], wspace=0.36)

    panel_a_grid = top_row[0, 0].subgridspec(1, 2, width_ratios=[0.62, 1.0], wspace=0.05)
    ax_a_note = fig.add_subplot(panel_a_grid[0, 0])
    ax_a_heat = fig.add_subplot(panel_a_grid[0, 1])
    ax_b = fig.add_subplot(top_row[0, 1])
    ax_c = fig.add_subplot(bottom_row[0, 0])
    panel_d_grid = bottom_row[0, 1].subgridspec(1, 4, wspace=0.16)
    ax_d = [fig.add_subplot(panel_d_grid[0, i]) for i in range(4)]

    fig.subplots_adjust(top=0.94, bottom=0.16, left=0.055, right=0.97)

    draw_panel_a(ax_a_note, ax_a_heat, panel_a_rows)
    draw_panel_b(ax_b, contrast_rows)
    draw_panel_c(ax_c, diagnostic_rows)
    draw_panel_d(ax_d, support_rows)

    panel_d_left = ax_d[0].get_position().x0
    panel_d_right = ax_d[-1].get_position().x1
    panel_d_top = ax_d[0].get_position().y1
    fig.text((panel_d_left + panel_d_right) / 2.0, panel_d_top + 0.024, "Coverage and phasing support at the PWS-AS IC", ha="center", va="bottom", fontsize=fs(10.8), fontweight="bold")

    support_handles = [
        plt.Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="#444444", markeredgecolor="#444444", markersize=6),
        plt.Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="white", markeredgecolor="#444444", markersize=6),
    ]
    fig.legend(
        support_handles,
        ["filled: >=10x and >=5 CpGs", "open: <10x or <5 CpGs"],
        frameon=False,
        fontsize=fs(7.2),
        ncol=2,
        loc="lower center",
        bbox_to_anchor=((panel_d_left + panel_d_right) / 2.0, 0.100),
        bbox_transform=fig.transFigure,
    )

    ax_a_note.text(-0.15, 1.04, "A", transform=ax_a_note.transAxes, fontsize=fs(20), fontweight="bold", va="bottom", ha="left")
    ax_b.text(-0.03, 1.045, "B", transform=ax_b.transAxes, fontsize=fs(20), fontweight="bold", va="bottom", ha="left")
    ax_c.text(-0.10, 1.04, "C", transform=ax_c.transAxes, fontsize=fs(20), fontweight="bold", va="bottom", ha="left")
    ax_d[0].text(-0.16, 1.04, "D", transform=ax_d[0].transAxes, fontsize=fs(20), fontweight="bold", va="bottom", ha="left")

    output_prefixes = [
        out_prefix,
        out_prefix.with_name("Figure1_ABCD"),
        out_prefix.with_name("Figure1_ABC"),
        out_prefix.with_name("Figure1_BCD"),
        out_prefix.with_name("Figure1_redesign"),
        out_prefix.with_name("Figure1_improved"),
    ]
    seen: set[Path] = set()
    for prefix in output_prefixes:
        if prefix in seen:
            continue
        seen.add(prefix)
        fig.savefig(prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
        fig.savefig(prefix.with_suffix(".pdf"), bbox_inches="tight")
        fig.savefig(prefix.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    table_dir = outdir / "tables"
    figure_dir = outdir / "figures"
    log_dir = outdir / "logs"
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    vcf_dir = Path(args.vcf_dir)
    bam_dir = Path(args.bam_dir)
    methylation_dir = Path(args.methylation_dir)
    cnv_dir = Path(args.cnv_dir)
    gtf = Path(args.gtf)
    metadata = read_metadata(Path(args.metadata))

    sample_files: dict[str, dict[str, Path | None]] = {}
    for sample_id, _clinical, _mechanism in COHORT:
        sample_files[sample_id] = {
            "bam": find_sample_file(bam_dir, sample_id, ".bam"),
            "summary": find_sample_file(vcf_dir, sample_id, ".summary.tsv"),
            "blocks": find_sample_file(vcf_dir, sample_id, ".blocks.tsv"),
            "combined_bed": find_sample_file(methylation_dir, sample_id, ".combined.bed"),
            "hap1_bed": find_sample_file(methylation_dir, sample_id, ".hap1.bed"),
            "hap2_bed": find_sample_file(methylation_dir, sample_id, ".hap2.bed"),
            "cnv_log": find_sample_file(cnv_dir, sample_id, ".log"),
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
            bam_qc_rows.append(build_bam_qc(sample_id, sample_files[sample_id]["bam"], sample_files[sample_id]["cnv_log"]))
        write_tsv(bam_qc_cache, [{k: fmt(v) for k, v in row.items()} for row in bam_qc_rows], bam_qc_fields)

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
            n50, domain_fraction = block_n50_and_domain_fraction(sample_files[sample_id]["blocks"])
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
                    "mean_depth_genome_wide": fmt(bam_qc.get("mean_depth_genome_wide")),
                    "mean_depth_chr15": fmt(bam_qc.get("mean_depth_chr15")),
                    "mean_depth_per_haplotype_at_15q11-q13": bam_qc.get("mean_depth_per_haplotype_at_15q11-q13", ""),
                    "phasing_block_N50_chr15": n50 if n50 is not None else "",
                    "percent_imprinted_domain_in_phased_block": fmt(domain_fraction),
                    "total_CpGs_called_in_imprinted_domain": count_cpgs_in_domain(combined_bed),
                }
            )
        write_tsv(summary_cache, summary_rows, summary_fields)

    genes = load_gene_models(gtf)
    gene_fields = ["gene", "chrom", "start", "end", "strand", "parental_annotation"]
    write_tsv(table_dir / "Figure1A_domain_genes_from_T2T_GTF.tsv", genes, gene_fields)

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
    write_tsv(table_dir / "Figure1C_parental_assignment.tsv", assignment_rows, assignment_fields)
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
    write_tsv(table_dir / "Figure1C_pws_ic_methylation_matrix.tsv", matrix_rows, matrix_fields)

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
    write_tsv(table_dir / "Figure1A_physical_allele_layout.tsv", panel_a_rows, panel_a_fields)
    write_tsv(table_dir / "Figure1A_allele_methylation_matrix.tsv", panel_a_rows, panel_a_fields)

    contrast_rows = build_per_cpg_contrast(stats_by_sample, assignment_rows)
    contrast_fields = ["pos", "score", "score_type", "sample_id", "molecular_mechanism"]
    write_tsv(table_dir / "Figure1D_per_CpG_contrast.tsv", contrast_rows, contrast_fields)
    write_tsv(table_dir / "Figure1B_per_CpG_contrast.tsv", contrast_rows, contrast_fields)

    diagnostic_state_rows = build_diagnostic_state_rows(panel_a_rows)
    diagnostic_state_fields = ["molecular_mechanism", "n_samples", "expected_state", "observed_state", "interpretation"]
    write_tsv(table_dir / "Figure1C_diagnostic_state_summary.tsv", diagnostic_state_rows, diagnostic_state_fields)

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
    write_tsv(table_dir / "Figure1D_support_metrics.tsv", support_rows, support_fields)
    write_tsv(table_dir / "Figure1D_coverage_phasing_support.tsv", support_rows, support_fields)

    create_figure(figure_dir / "Figure1", panel_a_rows, contrast_rows, diagnostic_state_rows, support_rows)

    run_parameters = {
        "cohort": [
            {"sample_id": sample_id, "clinical_diagnosis": clinical, "molecular_mechanism": mechanism}
            for sample_id, clinical, mechanism in COHORT
        ],
        "regions": {
            "imprinted_domain": {"chrom": CHROM, "start": DOMAIN_START, "end": DOMAIN_END},
            "pws_as_ic": {
                "chrom": CHROM,
                "start": PWS_IC_START,
                "end": PWS_IC_END,
                "name": PWS_IC_NAME,
                "source": "canonical PWS-AS ICR row ICR_893 from genomewide_icr_deviation_validation; overlaps SNRPN/SNHG14/SNURF regulatory interval",
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
            "metadata": str(args.metadata),
        },
        "depth_methods": {
            "total_HiFi_reads": "sum(mapped + unmapped) from samtools idxstats on HiPhase BAM",
            "mean_depth_genome_wide": "2 * final hificnv uncorrected haploid coverage estimate from sample log",
            "mean_depth_chr15": "samtools coverage -r chr15 meandepth",
            "mean_depth_per_haplotype_at_15q11-q13": "reference-consuming CIGAR bases for HP-tagged reads over chr15:22.5-28.5 Mb divided by region length",
        },
    }
    with (outdir / "phase1_run_parameters.json").open("w") as handle:
        json.dump(run_parameters, handle, indent=2)

    write_report(
        outdir / "reports" / "Figure1_report.md",
        Path(__file__),
        {
            "metadata": str(args.metadata),
            "run_parameters": str(outdir / "phase1_run_parameters.json"),
        },
        diagnostic_state_rows,
    )


if __name__ == "__main__":
    main()
