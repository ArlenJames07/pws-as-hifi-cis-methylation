

from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.metadata
import itertools
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm, to_rgba
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, Polygon, Rectangle


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
DEFAULT_ANALYSIS_DIR = DEFAULT_RESULTS_DIR / "analysis" / "01_evidence_matrix"
DEFAULT_STRUCTURAL_EVIDENCE = DEFAULT_ANALYSIS_DIR / "chr15_structural_evidence.tsv"
DEFAULT_GTF = Path("/home/rare/arlen/reference/chm13v22.sorted.gtf")

CHROM = "chr15"
DOMAIN_START = 22_500_000
DOMAIN_END = 28_500_000

PWS_IC_START = 22_691_258
PWS_IC_END = 22_693_494
PWS_IC_NAME = "ICR_893_SNHG14_SNRPN_SNURF"

# Raw single-molecule window for Figure 1A. This is intentionally a little
# broader than the diagnostic IC interval so that the local transition is visible.
MODBAM_REGION_START = 22_690_500
MODBAM_REGION_END = 22_695_000
MODBAM_REGION = f"{CHROM}:{MODBAM_REGION_START}-{MODBAM_REGION_END}"

# Boundary result shown only as genomic context in Figure 1. Figure 3 remains
# the formal boundary-mapping analysis.
SHARED_CORE_START = 22_691_000
SHARED_CORE_END = 22_694_700


# ---------------------------------------------------------------------------
# USER CONFIGURATION — EDIT ONLY THIS BLOCK
# ---------------------------------------------------------------------------
# If this file remains at scripts/figures/FIGURE_1.py inside the repository,
# the PROJECT_ROOT-derived defaults below should work without CLI arguments.
# Run simply with:
#
#     python3 scripts/figures/FIGURE_1.py
#

# HiPhase outputs: phased BAMs, BAM indexes and *.blocks.tsv files.
VCF_DIR = DEFAULT_PHASED_DIR
BAM_DIR = DEFAULT_PHASED_DIR

# BAMs used by ModBAMtools. They must retain MM/ML modified-base tags.
# Control and PWS-mUPD should also retain HP tags for haplotype separation.
# If the HiPhase BAMs retain MM/ML, leave MODBAM_DIR = BAM_DIR.
MODBAM_DIR = BAM_DIR

# Current pb-CpG-tools outputs: *.cpg.combined.*, *.cpg.hap1.* and
# *.cpg.hap2.*.  Keep the tag explicit so legacy symlinks with shorter names
# can never be selected accidentally.
METHYLATION_DIR = DEFAULT_METHYLATION_DIR
PBCPG_OUTPUT_TAG = ".cpg"

# HiFiCNV outputs used for depth/QC support.
CNV_DIR = DEFAULT_CNV_DIR

# Non-identifiable cohort metadata.
METADATA_PATH = DEFAULT_METADATA

# Canonical structural/deletion evidence produced by scripts/analysis.
# Figure 1 is a consumer of this table and must not re-infer deletion status.
STRUCTURAL_EVIDENCE_PATH = DEFAULT_STRUCTURAL_EVIDENCE

# Output directory for tables, ModBAMtools panels, reports and Figure 1.
OUTDIR = DEFAULT_OUTDIR

# ModBAMtools executable. Replace with a full path if it is not on PATH, e.g.
# MODBAMTOOLS_BIN = "/home/rare/miniforge3/envs/modbamtools/bin/modbamtools"
MODBAMTOOLS_BIN = "/home/rare/miniconda3/envs/modbam38/bin/modbamtools"
USE_EXTERNAL_MODBAMTOOLS = False

# Optional ModBAMtools gene annotation. Leave None unless you have prepared
# the indexed annotation required by ModBAMtools.
MODBAM_GTF: Path | None = None

# Standard GTF used for the compact gene context above the Panel B CNV tracks.
# Only a curated set of established 15q11-q13 genes is shown to avoid turning
# the overview into an unreadable transcript annotation panel.
PANEL_B_GTF: Path | None = DEFAULT_GTF
PANEL_B_GENE_NAMES = (
    "MKRN3", "MAGEL2", "NDN", "NPAP1", "SNHG14", "SNRPN", "UBE3A",
    "ATP10A", "GABRB3", "GABRA5", "GABRG3", "OCA2", "HERC2",
)

# Raw single-molecule interval plotted in Figure 1A.
MODBAM_PLOT_REGION = MODBAM_REGION

# Execution switches. Normally all three remain False for a complete run.
SKIP_BAM_QC = False
SKIP_MODBAMTOOLS = False
RENDER_ONLY = False

# Publication-grade QC / sensitivity controls.
STRICT_SCIENTIFIC_QC = True
GENERATE_EXTENDED_MODBAM = True
GENERATE_READ_LEVEL_SOURCE_DATA = True

# ModBAM display is standardized across every sample. We filter the plotted
# interval to primary alignments with MAPQ >= this threshold before calling
# ModBAMtools, rather than relying on renderer-specific defaults.
MODBAM_MIN_MAPQ = 20
MODBAM_PLOT_WIDTH = 2400
MODBAM_FREQUENCY_BIN_BP = 100

# Raster output is intended for final publication assembly. PDF and SVG remain
# the preferred editable/vector deliverables. The large 900-dpi PNG is intended
# to keep single-molecule marks and compact annotations legible at ordinary zoom.
PUBLICATION_DPI = 900

# Descriptive uncertainty for allele-level beta values is estimated by
# coverage-weighted bootstrap resampling of CpG sites. Formal group inference
# below uses participants, not CpGs or molecules, as the independent unit.
CPG_BOOTSTRAP_REPLICATES = 5000
SAMPLE_BOOTSTRAP_REPLICATES = 10000
BOOTSTRAP_SEED = 20260908

# Prespecified threshold sensitivity analysis.
THRESHOLD_SENSITIVITY = [
    (0.80, 0.20),
    (0.85, 0.15),
    (0.90, 0.10),
]

# HiFiCNV logs in this project historically report a haploid-normalized depth.
# Keep the raw value AND an explicitly scaled diploid-equivalent value; never
# hide the conversion behind an unexplained '* 2'. Verify this constant against
# the exact HiFiCNV release used for the submitted analysis.
HIFICNV_DEPTH_SCALE = 2.0

# Reproducibility manifest. Large BAMs can be hashed too, but that can add many
# minutes/hours. The code supports it; default is to hash files <= 100 MB and
# record size/mtime for larger inputs. Set HASH_LARGE_INPUTS=True for the final
# archival release if desired.
HASH_LARGE_INPUTS = False
HASH_MAX_BYTES = 100 * 1024 * 1024


# One representative raw-ModBAM panel is shown for every heatmap cohort block.
# This order matches the vertical order used by sorted_cohort().
MODBAM_GROUP_ORDER = [
    "PWS-DEL",
    "AS-DEL",
    "PWS-mUPD",
    "Disease control",
    "Control",
]

# Primary parental-state assignment is NOT based on fixed beta cutoffs.
# Maternal-like and paternal-like reference centroids are estimated automatically
# from the unaffected controls in each run. For each control, the more methylated
# haplotype anchors the maternal-like state and the less methylated haplotype
# anchors the paternal-like state; HP1/HP2 labels themselves are not parental.
# The decision boundary is the midpoint between the two empirical centroids.
PARENTAL_REFERENCE_METHOD = "unaffected-control median centroids"
CONTROL_REFERENCE_LEAVE_ONE_OUT = True

# These extreme-state thresholds are retained ONLY as a supplementary sensitivity
# analysis. They are not used for the primary parental-state calls in Figure 1.
EXTREME_MATERNAL_THRESHOLD = 0.85
EXTREME_PATERNAL_THRESHOLD = 0.15

# Depth-aware methylation QC.
# We do NOT convert a lower-depth haplotype into "missing". State assignment
# requires the minimum evidence below. The higher-coverage reference is reported
# descriptively in Figure 1C; falling below it is not a preflight
# warning and does not invalidate an otherwise estimable parental-state call.
MIN_STATE_CPGS = 3
MIN_STATE_MEAN_COVERAGE = 2.0
HIGHER_COVERAGE_CPGS = 5
HIGHER_COVERAGE_MEAN_COVERAGE = 10.0
MIN_MODBAM_FALLBACK_MOLECULES = 3

# Deletion-span parental-profile classifier. Unaffected-control haplotypes are
# first labelled M-like/P-like at the IC. Across each CN-defined deletion span,
# CpGs are retained only when the two control M-like profiles and the two
# control P-like profiles are internally consistent and well separated. The
# retained chromosome is then assigned by RMSE to the two empirical profiles.
DELETION_PROFILE_MIN_REFERENCE_DELTA = 0.50
DELETION_PROFILE_MAX_WITHIN_STATE_RANGE = 0.15
DELETION_PROFILE_MIN_SHARED_CPGS = 50
DELETION_PROFILE_MIN_SCORE_MAGNITUDE = 0.05

# Backward-compatible aliases retained for downstream configuration readers.
NOMINAL_CPGS = HIGHER_COVERAGE_CPGS
NOMINAL_MEAN_COVERAGE = HIGHER_COVERAGE_MEAN_COVERAGE
MIN_MEAN_COVERAGE = HIGHER_COVERAGE_MEAN_COVERAGE
MIN_CPGS = HIGHER_COVERAGE_CPGS

# HiFiCNV / deletion classification on T2T-CHM13v2.0.
# Coordinates are the same chr15 landmarks used by the manuscript Figure 5.
CN_PLOT_START = 18_000_000
CN_PLOT_END = 32_500_000
CN_DELETION_THRESHOLD = 1.35
CN_NORMAL_EXPECTED = 2.0
CN_MERGE_GAP_BP = 150_000
BP_NEAREST_MAX_DISTANCE = 500_000
BREAKPOINT_LANDMARKS = {
    "BP1": 20_940_000,
    "BP2": 21_070_000,
    "BP3": 26_050_000,
    "BP4": 26_460_000,
    "BP5": 31_840_000,
}

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
    # Okabe-Ito-inspired, color-vision-deficiency-friendly palette.
    "PWS-DEL": "#D55E00",
    "AS-DEL": "#0072B2",
    "PWS-mUPD": "#CC79A7",
    "Disease control": "#E69F00",
    "Control": "#4D4D4D",
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
    "M": "#D55E00",
    "P": "#0072B2",
    "deleted": "#E6E6E6",
    "missing": "#F5F5F5",
    "structural_unresolved": "#FFF2CC",
    "?": "#F3F3F3",
}
ABSENT_EDGE = "#A6A6A6"
TEXT_SCALE = 1.18


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
    ci_low: float | None = None
    ci_high: float | None = None
    values_by_pos: dict[int, tuple[float, float]] | None = None
    n_molecules: int | None = None
    data_source: str = "pb-CpG-tools BED"

    @property
    def estimable(self) -> bool:
        """Enough information for a descriptive parental-like state call."""
        return (
            self.mean_methylation is not None
            and self.n_cpgs >= MIN_STATE_CPGS
            and self.mean_coverage is not None
            and self.mean_coverage >= MIN_STATE_MEAN_COVERAGE
        )

    @property
    def sufficient(self) -> bool:
        """Higher-coverage reference category used for descriptive QC only."""
        return (
            self.mean_methylation is not None
            and self.n_cpgs >= HIGHER_COVERAGE_CPGS
            and self.mean_coverage is not None
            and self.mean_coverage >= HIGHER_COVERAGE_MEAN_COVERAGE
        )

@dataclass(frozen=True)
class ParentalReferenceModel:
    """Empirical parental-like methylation references learned from controls.

    The model is intentionally simple and transparent. Each unaffected control
    contributes its lower- and higher-methylated IC haplotype. The median of the
    high values defines the maternal-like reference, the median of the low values
    defines the paternal-like reference, and their midpoint is the equal-distance
    decision boundary. This removes the fixed 0.85/0.15 cutoffs from the primary
    analysis while preserving them as a Supplementary sensitivity analysis.
    """
    maternal_reference: float
    paternal_reference: float
    decision_boundary: float
    control_sample_ids: tuple[str, ...]
    method: str = PARENTAL_REFERENCE_METHOD

    @property
    def separation(self) -> float:
        return self.maternal_reference - self.paternal_reference


def validate_configuration() -> None:
    """Fail early when the hard-coded configuration is inconsistent."""
    if RENDER_ONLY:
        return

    required_dirs = {
        "VCF_DIR": Path(VCF_DIR),
        "BAM_DIR": Path(BAM_DIR),
        "MODBAM_DIR": Path(MODBAM_DIR),
        "METHYLATION_DIR": Path(METHYLATION_DIR),
        "CNV_DIR": Path(CNV_DIR),
    }
    missing_dirs = [
        f"{name}={path}"
        for name, path in required_dirs.items()
        if not path.exists()
    ]
    if missing_dirs:
        raise FileNotFoundError(
            "Configured input directories do not exist:\n"
            + "\n".join(f"- {item}" for item in missing_dirs)
        )

    if not Path(METADATA_PATH).exists():
        raise FileNotFoundError(
            f"METADATA_PATH does not exist: {METADATA_PATH}"
        )

    if not Path(STRUCTURAL_EVIDENCE_PATH).exists():
        raise FileNotFoundError(
            f"STRUCTURAL_EVIDENCE_PATH does not exist: {STRUCTURAL_EVIDENCE_PATH}\n"
            "Run `python3 scripts/analysis/run_analysis.py` before Figure 1."
        )

    if not SKIP_MODBAMTOOLS and USE_EXTERNAL_MODBAMTOOLS:
        executable = str(MODBAMTOOLS_BIN)
        if shutil.which(executable) is None and not Path(executable).exists():
            raise FileNotFoundError(
                "ModBAMtools executable was not found. Edit MODBAMTOOLS_BIN "
                "in the USER CONFIGURATION block. Current value: "
                f"{MODBAMTOOLS_BIN}"
            )


def safe_float(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return x if np.isfinite(x) else None


def safe_int(value: Any) -> int | None:
    x = safe_float(value)
    if x is None:
        return None
    rounded = round(x)
    if not np.isclose(x, rounded, rtol=0.0, atol=1e-6):
        return None
    return int(rounded)


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


def _parse_optional_bool(value: Any) -> bool | None:
    if isinstance(value, bool):
        return value
    if value is None:
        return None
    token = str(value).strip().lower()
    if token in {"true", "1", "yes"}:
        return True
    if token in {"false", "0", "no"}:
        return False
    return None


def load_analysis_structural_evidence(
    path: Path,
) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    """Load canonical chr15 structural evidence produced by scripts/analysis.

    Figure 1 deliberately does not call HiFiCNV/pbsv classifiers. Structural
    inference is performed once in the analysis layer and consumed here.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"Canonical structural evidence not found: {path}\n"
            "Run `python3 scripts/analysis/run_analysis.py` before Figure 1."
        )

    rows: list[dict[str, Any]] = [dict(row) for row in read_tsv(path)]
    expected_mechanism = {
        sample_id: mechanism
        for sample_id, _clinical, mechanism in sorted_cohort()
    }
    expected_samples = set(expected_mechanism)
    observed_samples = {str(row.get("sample_id", "")).strip() for row in rows}

    missing = sorted(expected_samples - observed_samples)
    extra = sorted(observed_samples - expected_samples)
    if missing or extra:
        parts = []
        if missing:
            parts.append("missing: " + ", ".join(missing))
        if extra:
            parts.append("unexpected: " + ", ".join(extra))
        raise RuntimeError(
            "Structural analysis table does not match the Figure 1 cohort ("
            + "; ".join(parts)
            + ")"
        )

    boolean_fields = (
        "expected_ic_deletion",
        "sv_support",
        "cnv_support",
        "pbsv_support",
        "hificnv_support",
        "cn1_eligible",
    )
    normalized_rows: list[dict[str, Any]] = []
    labels = sample_display_labels()
    clinical_by_sample = {
        sample_id: clinical
        for sample_id, clinical, _mechanism in sorted_cohort()
    }
    for row in rows:
        sample_id = str(row["sample_id"]).strip()
        normalized = dict(row)
        for field in boolean_fields:
            if field in normalized:
                normalized[field] = _parse_optional_bool(normalized[field])

        # Keep Figure 1's display vocabulary stable. The analysis layer may call
        # the orthogonal group "DiGeorge" whereas the figure uses "Disease control".
        normalized["molecular_mechanism"] = expected_mechanism[sample_id]
        normalized["display_label"] = labels[sample_id]
        normalized["clinical_diagnosis"] = clinical_by_sample[sample_id]

        # Compatibility aliases expected by existing Figure 1 report/QC code.
        if "pbsv_support" not in normalized:
            normalized["pbsv_support"] = normalized.get("sv_support")
        if "pbsv_note" not in normalized:
            normalized["pbsv_note"] = normalized.get("sv_note", "")
        if "hificnv_support" not in normalized:
            normalized["hificnv_support"] = normalized.get("cnv_support")
        if "hificnv_note" not in normalized:
            normalized["hificnv_note"] = normalized.get("cnv_note", "")

        normalized_rows.append(normalized)

    by_sample = {row["sample_id"]: row for row in normalized_rows}
    return by_sample, normalized_rows


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


def find_modbam_file(directory: Path, sample_id: str) -> Path | None:
    """
    Prefer phased BAMs because Figure 1A uses HP tags when biologically
    meaningful. Fall back to any sample BAM only if no phased BAM is present.
    """
    if not directory.exists():
        return None

    sample_directory = directory / sample_id
    if sample_directory.is_dir():
        matches = list(sample_directory.glob(f"*{sample_id}*.bam"))
    else:
        matches = list(directory.rglob(f"*{sample_id}*.bam"))

    if not matches:
        return None

    phased = [p for p in matches if "phased" in p.name.lower()]
    return choose_file(phased or matches, sample_id)


def run_command(args: list[str]) -> str:
    result = subprocess.run(args, check=True, text=True, capture_output=True)
    return result.stdout


def read_bed_region(
    path: Path | None,
    start: int,
    end: int,
    keep_values: bool = False,
) -> BedStats:
    """Read pb-CpG-tools BED values and estimate beta plus descriptive CI.

    The beta estimate is coverage-weighted across CpG sites. The 95% interval is
    a coverage-weighted CpG bootstrap and is explicitly *descriptive*; CpGs are
    not treated as independent biological replicates in group-level inference.
    """
    if path is None or not path.exists():
        return BedStats(
            values_by_pos={} if keep_values else None,
            data_source="pb-CpG-tools BED missing",
        )

    meth_values: list[float] = []
    cov_values: list[float] = []
    values_by_pos: dict[int, tuple[float, float]] = {}

    awk_script = "$1==chrom && $2>=start && $2<end {print}"
    result = subprocess.run(
        ["awk", "-v", f"chrom={CHROM}", "-v", f"start={start}", "-v", f"end={end}",
         awk_script, str(path)],
        check=True, text=True, capture_output=True,
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
        if not (np.isfinite(meth) and np.isfinite(cov)) or cov < 0:
            continue
        meth_values.append(meth)
        cov_values.append(cov)
        if keep_values:
            values_by_pos[row_start] = (meth, cov)

    if not meth_values:
        return BedStats(
            values_by_pos=values_by_pos if keep_values else None,
            data_source=f"pb-CpG-tools BED: {path.name} (no IC values)",
        )

    values = np.asarray(meth_values, dtype=float)
    weights = np.asarray(cov_values, dtype=float)
    if np.sum(weights) > 0:
        weighted_meth = float(np.average(values, weights=weights))
    else:
        weighted_meth = float(np.mean(values))

    ci_low = ci_high = None
    if len(values) >= 2 and CPG_BOOTSTRAP_REPLICATES > 0:
        seed_material = f"{path.name}|{start}|{end}".encode()
        seed_offset = int(hashlib.sha256(seed_material).hexdigest()[:8], 16)
        rng = np.random.default_rng((BOOTSTRAP_SEED + seed_offset) % (2**32 - 1))
        boot = np.empty(CPG_BOOTSTRAP_REPLICATES, dtype=float)
        n = len(values)
        for i in range(CPG_BOOTSTRAP_REPLICATES):
            idx = rng.integers(0, n, n)
            w = weights[idx]
            v = values[idx]
            boot[i] = float(np.average(v, weights=w)) if np.sum(w) > 0 else float(np.mean(v))
        ci_low, ci_high = np.quantile(boot, [0.025, 0.975])
        ci_low, ci_high = float(ci_low), float(ci_high)

    return BedStats(
        n_cpgs=len(values),
        mean_methylation=weighted_meth,
        mean_coverage=float(np.mean(weights)),
        ci_low=ci_low,
        ci_high=ci_high,
        values_by_pos=values_by_pos if keep_values else None,
        data_source=f"pb-CpG-tools BED: {path.name}",
    )



def read_modbam_haplotype_region_stats(
    bam_path: Path | None,
    hp: int | None,
    start: int,
    end: int,
) -> BedStats:
    """Quantify haplotype-specific methylation directly from MM/ML-tagged reads.

    This is used only when the corresponding pb-CpG-tools haplotype BED has no
    estimable values. It prevents a lower-depth sample from being labelled
    "missing" merely because pb-CpG-tools did not emit a haplotype BED estimate.

    The summary is calculated per reference position from MM/ML probabilities,
    then averaged across observed positions. The data source is explicitly
    recorded so BED-derived and direct ModBAM estimates remain auditable.
    """
    if bam_path is None or not Path(bam_path).exists():
        return BedStats(data_source="direct ModBAM unavailable")

    try:
        import pysam
    except ImportError:
        return BedStats(data_source="direct ModBAM unavailable: pysam missing")

    bam_path = Path(bam_path)
    position_values: dict[int, list[float]] = defaultdict(list)
    molecules: set[str] = set()

    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam.fetch(CHROM, start, end):
                if (
                    read.is_unmapped
                    or read.is_secondary
                    or read.is_supplementary
                    or read.is_duplicate
                    or read.mapping_quality < MODBAM_MIN_MAPQ
                ):
                    continue
                if hp is not None:
                    if not read.has_tag("HP"):
                        continue
                    try:
                        if int(read.get_tag("HP")) != hp:
                            continue
                    except (ValueError, TypeError):
                        continue

                mods = read.modified_bases
                if not mods:
                    continue

                q_to_r = {
                    q: r
                    for q, r in read.get_aligned_pairs(matches_only=False)
                    if q is not None and r is not None
                }

                contributed = False
                for key, calls in mods.items():
                    canonical, _strand, modification = key
                    if str(canonical).upper() != "C":
                        continue
                    # PacBio 5mC ModBAMs typically use C+m. Keep only methyl-C
                    # style codes; unknown non-methyl-C modifications are skipped.
                    mod_code = str(modification).lower()
                    if mod_code not in {"m", "5mc", "c+m"}:
                        continue
                    for qpos, qual in calls:
                        rpos = q_to_r.get(qpos)
                        if rpos is None or not (start <= rpos < end):
                            continue
                        if qual is None or qual < 0:
                            continue
                        position_values[rpos].append(float(qual) / 255.0)
                        contributed = True

                if contributed:
                    molecules.add(read.query_name)
    except Exception as exc:
        return BedStats(data_source=f"direct ModBAM quantification failed: {exc}")

    if not position_values:
        return BedStats(
            n_molecules=len(molecules),
            data_source="direct ModBAM MM/ML",
        )

    positions = sorted(position_values)
    site_beta = np.asarray(
        [np.mean(position_values[pos]) for pos in positions],
        dtype=float,
    )
    site_cov = np.asarray(
        [len(position_values[pos]) for pos in positions],
        dtype=float,
    )
    beta = float(np.average(site_beta, weights=site_cov)) if site_cov.sum() > 0 else float(site_beta.mean())

    ci_low = ci_high = None
    if len(site_beta) >= 2 and CPG_BOOTSTRAP_REPLICATES > 0:
        seed_material = f"{bam_path.name}|HP={hp}|{start}|{end}".encode()
        seed_offset = int(hashlib.sha256(seed_material).hexdigest()[:8], 16)
        rng = np.random.default_rng((BOOTSTRAP_SEED + seed_offset) % (2**32 - 1))
        boot = np.empty(CPG_BOOTSTRAP_REPLICATES, dtype=float)
        n = len(site_beta)
        for i in range(CPG_BOOTSTRAP_REPLICATES):
            idx = rng.integers(0, n, n)
            w = site_cov[idx]
            v = site_beta[idx]
            boot[i] = float(np.average(v, weights=w)) if w.sum() > 0 else float(v.mean())
        ci_low, ci_high = map(float, np.quantile(boot, [0.025, 0.975]))

    return BedStats(
        n_cpgs=len(positions),
        mean_methylation=beta,
        mean_coverage=float(site_cov.mean()),
        ci_low=ci_low,
        ci_high=ci_high,
        values_by_pos={pos: (float(np.mean(position_values[pos])), float(len(position_values[pos]))) for pos in positions},
        n_molecules=len(molecules),
        data_source="direct ModBAM MM/ML",
    )


def complete_haplotype_stats_from_modbam(
    files: dict[str, Path | None],
    stats: dict[str, BedStats],
    mechanism: str,
) -> dict[str, BedStats]:
    """Use ModBAM MM/ML+HP evidence only when BED evidence is absent/non-estimable."""
    modbam = files.get("modbam")
    if modbam is None or not Path(modbam).exists():
        return stats

    if mechanism in {"Control", "Disease control", "PWS-mUPD"}:
        for label, hp in (("hap1", 1), ("hap2", 2)):
            current = stats[label]
            if current.estimable:
                continue
            direct_stats = read_modbam_haplotype_region_stats(
                Path(modbam), hp, PWS_IC_START, PWS_IC_END
            )
            # Prefer direct ModBAM evidence if it supplies more informative CpGs.
            if (
                direct_stats.mean_methylation is not None
                and direct_stats.n_molecules is not None
                and direct_stats.n_molecules >= MIN_MODBAM_FALLBACK_MOLECULES
                and direct_stats.n_cpgs > current.n_cpgs
            ):
                stats[label] = direct_stats

    elif mechanism in {"PWS-DEL", "AS-DEL"}:
        current = stats["combined_fallback"]
        if not current.estimable:
            direct_stats = read_modbam_haplotype_region_stats(
                Path(modbam), None, PWS_IC_START, PWS_IC_END
            )
            if (
                direct_stats.mean_methylation is not None
                and direct_stats.n_molecules is not None
                and direct_stats.n_molecules >= MIN_MODBAM_FALLBACK_MOLECULES
                and direct_stats.n_cpgs > current.n_cpgs
            ):
                stats["combined_fallback"] = direct_stats

    return stats


def _control_sample_ids() -> list[str]:
    return [
        sample_id
        for sample_id, _clinical, mechanism in sorted_cohort()
        if mechanism == "Control"
    ]


def estimate_parental_reference_model(
    stats_by_sample: dict[str, dict[str, BedStats]],
    exclude_sample_id: str | None = None,
) -> ParentalReferenceModel:
    """Estimate maternal-like and paternal-like IC centroids from controls.

    No fixed methylation cutoff is used. For each eligible unaffected control,
    the two observed haplotype beta values are sorted: the higher value anchors
    the maternal-like state and the lower value anchors the paternal-like state.
    The cohort reference is the median of those per-control anchors.

    When a control itself is classified and CONTROL_REFERENCE_LEAVE_ONE_OUT is
    enabled, that sample is excluded from reference estimation to avoid using a
    control to define and validate its own state. With two controls this leaves
    one independent control pair for the held-out sample.
    """
    high_values: list[float] = []
    low_values: list[float] = []
    used_ids: list[str] = []

    for sample_id in _control_sample_ids():
        if exclude_sample_id is not None and sample_id == exclude_sample_id:
            continue
        sample_stats = stats_by_sample.get(sample_id, {})
        h1 = sample_stats.get("hap1")
        h2 = sample_stats.get("hap2")
        if h1 is None or h2 is None:
            continue
        b1 = safe_float(h1.mean_methylation)
        b2 = safe_float(h2.mean_methylation)
        if b1 is None or b2 is None:
            continue
        if not (h1.estimable and h2.estimable):
            continue
        low, high = sorted((float(b1), float(b2)))
        low_values.append(low)
        high_values.append(high)
        used_ids.append(sample_id)

    if not high_values or not low_values:
        raise RuntimeError(
            "Cannot estimate parental methylation references: no unaffected "
            "control has two estimable IC haplotypes."
        )

    maternal = float(np.median(high_values))
    paternal = float(np.median(low_values))
    if not np.isfinite(maternal) or not np.isfinite(paternal) or maternal <= paternal:
        raise RuntimeError(
            "Invalid control-derived parental reference model: maternal-like "
            f"centroid={maternal}, paternal-like centroid={paternal}."
        )

    return ParentalReferenceModel(
        maternal_reference=maternal,
        paternal_reference=paternal,
        decision_boundary=(maternal + paternal) / 2.0,
        control_sample_ids=tuple(used_ids),
    )


def parental_reference_for_sample(
    sample_id: str,
    mechanism: str,
    stats_by_sample: dict[str, dict[str, BedStats]],
    global_reference: ParentalReferenceModel,
) -> ParentalReferenceModel:
    """Return global or leave-one-control-out reference for one sample."""
    if mechanism == "Control" and CONTROL_REFERENCE_LEAVE_ONE_OUT:
        try:
            loo = estimate_parental_reference_model(
                stats_by_sample, exclude_sample_id=sample_id
            )
            return loo
        except RuntimeError:
            pass
    return global_reference


def reference_relative_metrics(
    stats: BedStats,
    reference: ParentalReferenceModel,
) -> dict[str, Any]:
    """Continuous distance metrics relative to empirical parental centroids."""
    beta = safe_float(stats.mean_methylation)
    if beta is None:
        return {
            "distance_to_maternal": None,
            "distance_to_paternal": None,
            "nearest_reference": "unresolved",
            "reference_margin": None,
            "parental_axis_score": None,
        }
    dm = abs(beta - reference.maternal_reference)
    dp = abs(beta - reference.paternal_reference)
    nearest = "maternal" if dm < dp else "paternal" if dp < dm else "equidistant"
    denom = dm + dp
    margin = abs(dp - dm) / denom if denom > 0 else 0.0
    separation = reference.separation
    axis = (2.0 * (beta - reference.decision_boundary) / separation) if separation > 0 else None
    return {
        "distance_to_maternal": dm,
        "distance_to_paternal": dp,
        "nearest_reference": nearest,
        "reference_margin": margin,
        # -1 at the paternal centroid, +1 at the maternal centroid. Values may
        # extend beyond that range if an observation is more extreme than a reference.
        "parental_axis_score": axis,
    }


def methylation_pattern(
    stats: BedStats,
    reference: ParentalReferenceModel,
) -> str:
    """Depth-aware, control-calibrated parental-state assignment.

    Primary state calls no longer use fixed beta thresholds. An estimable allele
    is maternal-like when its descriptive CpG-bootstrap interval lies entirely
    above the empirical equal-distance boundary and paternal-like when the
    interval lies entirely below it. Intervals overlapping the boundary are
    labelled uncertain rather than being forced into either parental state.

    If a descriptive interval is unavailable, the nearest empirical reference
    centroid is used as a fallback. The CpG interval is a measurement-support
    device only; participant-level biological inference remains sample-based.
    """
    beta = safe_float(stats.mean_methylation)
    if beta is None:
        return "missing"
    if not stats.estimable:
        return "low-support"

    lo = safe_float(stats.ci_low)
    hi = safe_float(stats.ci_high)
    boundary = reference.decision_boundary
    if lo is not None and hi is not None:
        if lo > boundary:
            return "maternal-pattern"
        if hi < boundary:
            return "paternal-pattern"
        return "uncertain"

    metrics = reference_relative_metrics(stats, reference)
    if metrics["nearest_reference"] == "maternal":
        return "maternal-pattern"
    if metrics["nearest_reference"] == "paternal":
        return "paternal-pattern"
    return "uncertain"


def pattern_short(pattern: str) -> str:
    return {
        "maternal-pattern": "M",
        "paternal-pattern": "P",
        "uncertain": "U",
        "low-support": "low",
        "absent": "absent",
    }.get(pattern, "?")


def build_parental_reference_rows(
    stats_by_sample: dict[str, dict[str, BedStats]],
    global_reference: ParentalReferenceModel,
) -> list[dict[str, Any]]:
    """Auditable table of global and leave-one-control-out references."""
    rows: list[dict[str, Any]] = [{
        "scope": "global",
        "excluded_control": "",
        "controls_used": ";".join(global_reference.control_sample_ids),
        "maternal_reference_beta": fmt(global_reference.maternal_reference),
        "paternal_reference_beta": fmt(global_reference.paternal_reference),
        "decision_boundary_beta": fmt(global_reference.decision_boundary),
        "reference_separation": fmt(global_reference.separation),
        "method": global_reference.method,
    }]

    if CONTROL_REFERENCE_LEAVE_ONE_OUT:
        for sample_id in _control_sample_ids():
            try:
                ref = estimate_parental_reference_model(
                    stats_by_sample, exclude_sample_id=sample_id
                )
            except RuntimeError:
                continue
            rows.append({
                "scope": "leave-one-control-out",
                "excluded_control": sample_id,
                "controls_used": ";".join(ref.control_sample_ids),
                "maternal_reference_beta": fmt(ref.maternal_reference),
                "paternal_reference_beta": fmt(ref.paternal_reference),
                "decision_boundary_beta": fmt(ref.decision_boundary),
                "reference_separation": fmt(ref.separation),
                "method": ref.method,
            })
    return rows


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
    """Exact average depth over every base in the requested interval."""
    chrom, start, end = parse_region(region) if ":" in region else (region, 1, None)
    if end is None:
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

    # samtools depth uses 1-based closed region strings. -aa retains zero-depth bases.
    stdout = run_command(["samtools", "depth", "-aa", "-r", f"{chrom}:{start}-{end}", str(bam)])
    depths = []
    for line in stdout.splitlines():
        fields = line.split("\t")
        if len(fields) >= 3:
            try:
                depths.append(float(fields[2]))
            except ValueError:
                pass
    region_len = end - start + 1
    if region_len <= 0:
        return None
    # Some samtools builds may omit trailing zero-depth bases despite -aa on a region;
    # summing observed positions and dividing by the explicit interval length is robust.
    return float(sum(depths) / region_len)


def haplotype_depths_from_bam(bam: Path, region: str, region_len: int) -> dict[str, float]:
    """Exact mean HP1/HP2 depth from aligned reference blocks using pysam.

    Each aligned reference base contributed by an HP-tagged primary read counts
    once. Mean depth is total overlapping aligned bases / explicit region length.
    This avoids rescaling contig-level `samtools coverage` summaries.
    """
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("pysam is required for haplotype-specific depth QC") from exc

    chrom, start, end = parse_region(region)
    overlap_bases = {1: 0, 2: 0}
    with pysam.AlignmentFile(str(bam), "rb") as handle:
        for read in handle.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
                continue
            if not read.has_tag("HP"):
                continue
            try:
                hp = int(read.get_tag("HP"))
            except (ValueError, TypeError):
                continue
            if hp not in overlap_bases:
                continue
            for bstart, bend in read.get_blocks():
                ov = max(0, min(bend, end) - max(bstart, start))
                overlap_bases[hp] += ov
    denom = max(region_len, 1)
    return {"hap1": overlap_bases[1] / denom, "hap2": overlap_bases[2] / denom}


def parse_hificnv_depth(cnv_log: Path | None) -> tuple[float | None, float | None, float | None, float | None]:
    """Return raw and explicitly scaled HiFiCNV depth values.

    Returns: raw_uncorrected, raw_gc_corrected, scaled_uncorrected,
    scaled_gc_corrected. Keeping both prevents a hidden normalization factor.
    """
    if cnv_log is None or not cnv_log.exists():
        return None, None, None, None
    matches: list[tuple[str, str]] = []
    pattern = re.compile(r"Uncorrected:\s*([0-9.]+)\s+GC-Corrected:\s*([0-9.]+)")
    for line in cnv_log.read_text(errors="replace").splitlines():
        match = pattern.search(line)
        if match:
            matches.append((match.group(1), match.group(2)))
    if not matches:
        return None, None, None, None
    raw_unc, raw_gc = map(float, matches[-1])
    return raw_unc, raw_gc, raw_unc * HIFICNV_DEPTH_SCALE, raw_gc * HIFICNV_DEPTH_SCALE


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
            "sample_id": sample_id, "bam_file": str(bam or ""), "total_HiFi_reads": "",
            "mean_depth_genome_wide_raw_hificnv": "", "mean_depth_genome_wide_gc_raw_hificnv": "",
            "mean_depth_genome_wide": "", "mean_depth_genome_wide_gc_corrected": "",
            "mean_depth_chr15": "", "mean_depth_per_haplotype_at_15q11-q13": "", "chr15_length": "",
        }

    total_reads, chrom_lengths = bam_idxstats(bam)
    raw_depth, raw_gc, genome_depth, genome_depth_gc = parse_hificnv_depth(cnv_log)
    mean_depth_chr15 = samtools_coverage_mean_depth(bam, CHROM)
    region = f"{CHROM}:{DOMAIN_START}-{DOMAIN_END}"
    hap_depths = haplotype_depths_from_bam(bam, region, DOMAIN_END - DOMAIN_START + 1)

    return {
        "sample_id": sample_id,
        "bam_file": str(bam),
        "total_HiFi_reads": total_reads,
        "mean_depth_genome_wide_raw_hificnv": raw_depth,
        "mean_depth_genome_wide_gc_raw_hificnv": raw_gc,
        "mean_depth_genome_wide": genome_depth,
        "mean_depth_genome_wide_gc_corrected": genome_depth_gc,
        "mean_depth_chr15": mean_depth_chr15,
        "mean_depth_per_haplotype_at_15q11-q13": ";".join(f"{k}={v:.3f}" for k, v in sorted(hap_depths.items())),
        "chr15_length": chrom_lengths.get(CHROM, ""),
    }


# ---------------------------------------------------------------------------
# Structural evidence, scientific QC, uncertainty and reproducibility
# ---------------------------------------------------------------------------

def _open_text_maybe_gzip(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open("rt")


def parse_info_field(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for token in info.split(";"):
        if not token:
            continue
        if "=" in token:
            key, value = token.split("=", 1)
            out[key] = value
        else:
            out[token] = "True"
    return out


def sv_vcf_deletion_overlaps_ic(path: Path | None) -> tuple[bool | None, str]:
    """Conservatively test whether a DEL call spans the complete PWS/AS IC."""
    if path is None or not path.exists():
        return None, "SV VCF unavailable"
    try:
        with _open_text_maybe_gzip(path) as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                f = line.rstrip().split("\t")
                if len(f) < 8 or f[0] != CHROM:
                    continue
                try:
                    pos = int(f[1])
                except ValueError:
                    continue
                info = parse_info_field(f[7])
                svtype = info.get("SVTYPE", "").upper()
                if svtype != "DEL":
                    continue
                try:
                    end = int(info.get("END", pos + abs(int(info.get("SVLEN", "0").split(",")[0]))))
                except ValueError:
                    continue
                left, right = min(pos, end), max(pos, end)
                if left <= PWS_IC_START and right >= PWS_IC_END:
                    return True, f"SV DEL {CHROM}:{left}-{right}"
        return False, "No SV DEL spanning complete IC"
    except Exception as exc:
        return None, f"SV VCF parse error: {exc}"



@dataclass
class CNSegment:
    chrom: str
    start: int
    end: int
    copy_number: float
    source: str


def _extract_cn_from_fields(fields: list[str], path: Path) -> float | None:
    """Extract copy number from HiFiCNV BED/BEDGraph records.

    Supports:
      * modern *.copynum.bedgraph (4th column = copy number),
      * records containing CN=1 / copy_number=1 / copynum=1,
      * legacy HiFiCNV *.cnv.bed with an integer-like CN field.

    Numeric fallback is intentionally restricted to HiFiCNV/CNV-named files.
    """
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
            if 0 <= value <= 10:
                return value
        except ValueError:
            pass

    # Legacy project output can be *.hificnv.cnv.bed with CN as a bare
    # numeric tail field. Restrict the heuristic to these files and prefer
    # integer-like values in the biologically plausible range.
    if "hificnv" in lower_name or "cnv" in lower_name:
        candidates: list[float] = []
        for token in tail:
            clean = token.strip(",;")
            try:
                value = float(clean)
            except ValueError:
                continue
            if 0 <= value <= 6 and abs(value - round(value)) <= 0.05:
                candidates.append(value)
        if candidates:
            # CN tends to be the first integer-like dosage field in legacy BEDs.
            return float(candidates[0])

    return None


def read_hificnv_cn_segments(path: Path | None) -> list[CNSegment]:
    """Read a HiFiCNV CN track without assuming one release-specific layout."""
    if path is None or not Path(path).exists():
        return []

    path = Path(path)
    segments: list[CNSegment] = []
    opener = gzip.open if str(path).endswith(".gz") else open

    try:
        with opener(path, "rt", errors="replace") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip().split()
                if len(fields) < 4:
                    continue
                if fields[0] != CHROM:
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
    except Exception:
        return []

    return sorted(segments, key=lambda s: (s.start, s.end))


def _merge_deleted_cn_segments(
    segments: list[CNSegment],
    threshold: float = CN_DELETION_THRESHOLD,
) -> list[dict[str, Any]]:
    """Merge neighboring CN-loss bins into candidate deletion intervals."""
    loss = [s for s in segments if s.copy_number <= threshold]
    if not loss:
        return []

    events: list[dict[str, Any]] = []
    current = {
        "start": loss[0].start,
        "end": loss[0].end,
        "segments": [loss[0]],
    }

    for seg in loss[1:]:
        if seg.start <= current["end"] + CN_MERGE_GAP_BP:
            current["end"] = max(current["end"], seg.end)
            current["segments"].append(seg)
        else:
            events.append(current)
            current = {"start": seg.start, "end": seg.end, "segments": [seg]}
    events.append(current)

    for event in events:
        segs = event["segments"]
        lengths = np.asarray([s.end - s.start for s in segs], dtype=float)
        cns = np.asarray([s.copy_number for s in segs], dtype=float)
        event["mean_cn"] = float(np.average(cns, weights=lengths)) if lengths.sum() > 0 else float(np.mean(cns))
        event["min_cn"] = float(np.min(cns))
        event["n_segments"] = len(segs)

    return events


def _nearest_breakpoint(position: int, candidates: list[str]) -> tuple[str | None, int | None]:
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
    """Classify recurrent PWS/AS deletion type from CN transition coordinates.

    The label is deliberately '-like': HiFiCNV segmentation estimates dosage
    transitions, whereas sequence-resolved breakpoints inside segmental
    duplications may be more precise when pbsv/assembly evidence is available.
    """
    if event is None:
        return {
            "deletion_type": "no chr15 deletion",
            "left_landmark": "",
            "right_landmark": "",
            "left_distance_bp": "",
            "right_distance_bp": "",
        }

    left_name, left_dist = _nearest_breakpoint(event["start"], ["BP1", "BP2"])
    right_name, right_dist = _nearest_breakpoint(event["end"], ["BP3", "BP4", "BP5"])

    if left_name == "BP1" and right_name == "BP3":
        dtype = "Type I-like (BP1-BP3)"
    elif left_name == "BP2" and right_name == "BP3":
        dtype = "Type II-like (BP2-BP3)"
    elif left_name and right_name:
        dtype = f"atypical/extended {left_name}-{right_name}-like"
    else:
        dtype = "atypical/extended deletion"

    return {
        "deletion_type": dtype,
        "left_landmark": left_name or "unassigned",
        "right_landmark": right_name or "unassigned",
        "left_distance_bp": left_dist if left_dist is not None else "",
        "right_distance_bp": right_dist if right_dist is not None else "",
    }


def cnv_copy_number_evidence(
    path: Path | None,
) -> tuple[bool | None, str, dict[str, Any] | None, list[CNSegment]]:
    """Use CN dosage itself as primary evidence of hemizygous deletion.

    CN≈1 spanning the complete IC is sufficient copy-number evidence even when
    pbsv does not emit a sequence-resolved DEL. This is exactly the distinction
    between a CNV dosage call and an SV breakpoint call.
    """
    if path is None or not Path(path).exists():
        return None, "HiFiCNV copy-number track unavailable", None, []

    segments = read_hificnv_cn_segments(Path(path))
    if not segments:
        return None, "No parseable chr15 copy-number records", None, []

    events = _merge_deleted_cn_segments(segments)
    spanning = [
        event
        for event in events
        if event["start"] <= PWS_IC_START and event["end"] >= PWS_IC_END
    ]

    if spanning:
        # Prefer the event with the greatest overlap around the IC.
        event = max(spanning, key=lambda e: e["end"] - e["start"])
        return (
            True,
            (
                f"HiFiCNV CN-loss {CHROM}:{event['start']}-{event['end']} "
                f"(mean CN={event['mean_cn']:.2f}, threshold<={CN_DELETION_THRESHOLD:.2f})"
            ),
            event,
            segments,
        )

    # If a full bedGraph is available, explicitly summarize dosage across IC.
    overlapping = [
        s for s in segments
        if s.start < PWS_IC_END and s.end > PWS_IC_START
    ]
    if overlapping:
        lengths = np.asarray(
            [max(0, min(s.end, PWS_IC_END) - max(s.start, PWS_IC_START)) for s in overlapping],
            dtype=float,
        )
        cns = np.asarray([s.copy_number for s in overlapping], dtype=float)
        mean_cn = float(np.average(cns, weights=lengths)) if lengths.sum() > 0 else float(np.mean(cns))
        return False, f"IC mean CN={mean_cn:.2f}; no CN<=1 deletion spanning complete IC", None, segments

    return None, "CN track has chr15 records but no bins overlapping the IC", None, segments


def structural_ic_status(
    sample_id: str,
    mechanism: str,
    files: dict[str, Path | None],
) -> dict[str, Any]:
    """Integrate orthogonal dosage (HiFiCNV) and breakpoint (pbsv) evidence.

    HiFiCNV CN=1 is treated as direct evidence of hemizygous dosage. pbsv is
    retained as orthogonal sequence-resolved support, not a mandatory gate.
    """
    sv_ok, sv_note = sv_vcf_deletion_overlaps_ic(files.get("sv_vcf"))
    cn_ok, cn_note, cn_event, _segments = cnv_copy_number_evidence(files.get("cn_track"))

    expected_deletion = mechanism in {"PWS-DEL", "AS-DEL"}
    confirmed = (cn_ok is True) or (sv_ok is True)
    evidence_available = (cn_ok is not None) or (sv_ok is not None)

    if expected_deletion:
        if confirmed:
            status = "confirmed"
        elif evidence_available:
            status = "not_confirmed"
        else:
            status = "unavailable"
    else:
        status = "not_expected"

    deletion_class = classify_deletion_from_cn_event(cn_event)

    # Keep the two evidence modalities distinct in the output.
    if cn_ok is True and sv_ok is not True:
        evidence_basis = "HiFiCNV CN=1 dosage support; pbsv breakpoint DEL not required"
    elif cn_ok is True and sv_ok is True:
        evidence_basis = "concordant HiFiCNV CN-loss + pbsv DEL"
    elif sv_ok is True:
        evidence_basis = "pbsv DEL support; HiFiCNV CN track unavailable/non-confirmatory"
    else:
        evidence_basis = "no confirmed IC deletion"

    return {
        "sample_id": sample_id,
        "molecular_mechanism": mechanism,
        "expected_ic_deletion": expected_deletion,
        "ic_deletion_status": status,
        "evidence_basis": evidence_basis,
        "sv_support": sv_ok,
        "sv_note": sv_note,
        "cnv_support": cn_ok,
        "cnv_note": cn_note,
        "cn_event_start": cn_event["start"] if cn_event else "",
        "cn_event_end": cn_event["end"] if cn_event else "",
        "cn_event_size_mb": fmt((cn_event["end"] - cn_event["start"]) / 1e6, 3) if cn_event else "",
        "cn_event_mean_cn": fmt(cn_event["mean_cn"], 3) if cn_event else "",
        **deletion_class,
    }


def find_hificnv_cn_track(directory: Path, sample_id: str) -> Path | None:
    """Prefer the continuous HiFiCNV copy-number bedGraph, then legacy CNV BED."""
    candidates = [
        find_sample_file(directory, sample_id, ".copynum.bedgraph"),
        find_sample_file(directory, sample_id, ".copynum.bedgraph.gz"),
        find_sample_file(directory, sample_id, ".cnv.bed"),
        find_sample_file(directory, sample_id, ".cnv.bed.gz"),
        find_sample_file(directory, sample_id, ".bedgraph"),
    ]
    for path in candidates:
        if path is not None and Path(path).exists():
            return Path(path)
    return None


def build_cn_classification_rows(
    sample_files: dict[str, dict[str, Path | None]],
    structural_by_sample: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    """Table used by Supplementary CN figures and deletion-type reporting."""
    labels = sample_display_labels()
    rows: list[dict[str, Any]] = []
    for sample_id, clinical, mechanism in sorted_cohort():
        s = structural_by_sample[sample_id]
        rows.append({
            "sample_id": sample_id,
            "display_label": labels[sample_id],
            "clinical_diagnosis": clinical,
            "molecular_mechanism": mechanism,
            "cn_track": str(sample_files[sample_id].get("cn_track") or ""),
            "ic_deletion_status": s["ic_deletion_status"],
            "evidence_basis": s["evidence_basis"],
            "cn_event_start": s["cn_event_start"],
            "cn_event_end": s["cn_event_end"],
            "cn_event_size_mb": s["cn_event_size_mb"],
            "cn_event_mean_cn": s["cn_event_mean_cn"],
            "deletion_type": s["deletion_type"],
            "left_landmark": s["left_landmark"],
            "right_landmark": s["right_landmark"],
            "left_distance_bp": s["left_distance_bp"],
            "right_distance_bp": s["right_distance_bp"],
            "pbsv_support": s["sv_support"],
            "pbsv_note": s["sv_note"],
            "hificnv_support": s["cnv_support"],
            "hificnv_note": s["cnv_note"],
        })
    return rows



def build_cn_segment_source_rows(
    sample_files: dict[str, dict[str, Path | None]],
) -> list[dict[str, Any]]:
    """Long-format CN segments underlying supplementary deletion plots."""
    labels = sample_display_labels()
    rows: list[dict[str, Any]] = []
    for sample_id, _clinical, mechanism in sorted_cohort():
        path = sample_files[sample_id].get("cn_track")
        for seg in read_hificnv_cn_segments(path):
            if seg.end <= CN_PLOT_START or seg.start >= CN_PLOT_END:
                continue
            rows.append({
                "sample_id": sample_id,
                "display_label": labels[sample_id],
                "molecular_mechanism": mechanism,
                "chrom": seg.chrom,
                "start": seg.start,
                "end": seg.end,
                "copy_number": fmt(seg.copy_number, 3),
                "source_file": seg.source,
            })
    return rows


def _plot_cn_profile_axis(
    ax: plt.Axes,
    sample_id: str,
    display_label: str,
    mechanism: str,
    cn_track: Path | None,
    structural_row: dict[str, Any],
    compact: bool = False,
) -> None:
    """Draw a standardized chr15 copy-number profile for one participant."""
    segments = read_hificnv_cn_segments(cn_track)
    in_window = [
        s for s in segments
        if s.start < CN_PLOT_END and s.end > CN_PLOT_START
    ]

    ax.set_xlim(CN_PLOT_START / 1e6, CN_PLOT_END / 1e6)
    ax.set_ylim(-0.1, 3.6)
    ax.axhline(2.0, color="#777777", lw=0.8, ls="--", alpha=0.8)
    ax.axhline(1.0, color="#999999", lw=0.6, ls=":", alpha=0.7)

    # Recurrent breakpoint landmarks.
    for name, pos in BREAKPOINT_LANDMARKS.items():
        ax.axvline(pos / 1e6, color="#B0B0B0", lw=0.6, ls="--", zorder=0)
        if not compact:
            ax.text(
                pos / 1e6, 3.38, name,
                ha="center", va="top", fontsize=6, color="#777777"
            )

    # PWS/AS imprinting centre.
    ax.axvspan(
        PWS_IC_START / 1e6,
        PWS_IC_END / 1e6,
        color="#CC79A7",
        alpha=0.18,
        linewidth=0,
    )

    if in_window:
        # Full bedGraph: draw every segment. Variant-only BEDs simply overlay
        # their called segments on an expected CN=2 baseline.
        if "copynum" not in str(cn_track).lower() and "bedgraph" not in str(cn_track).lower():
            ax.plot(
                [CN_PLOT_START / 1e6, CN_PLOT_END / 1e6],
                [2.0, 2.0],
                color="#BDBDBD", lw=1.0, zorder=1,
            )

        for seg in in_window:
            left = max(seg.start, CN_PLOT_START) / 1e6
            right = min(seg.end, CN_PLOT_END) / 1e6
            color = MECHANISM_COLORS[mechanism]
            ax.plot([left, right], [seg.copy_number, seg.copy_number], color=color, lw=1.8)
            ax.vlines([left, right], max(0, seg.copy_number - 0.06), seg.copy_number + 0.06, color=color, lw=0.5)
    else:
        ax.text(
            0.5, 0.5, "No parseable chr15 CN track",
            transform=ax.transAxes, ha="center", va="center",
            fontsize=8, color="#777777"
        )

    event_start = safe_float(structural_row.get("cn_event_start"))
    event_end = safe_float(structural_row.get("cn_event_end"))
    if event_start is not None and event_end is not None:
        ax.axvspan(
            event_start / 1e6,
            event_end / 1e6,
            color=MECHANISM_COLORS[mechanism],
            alpha=0.08,
            linewidth=0,
        )

    dtype = structural_row.get("deletion_type", "unresolved")
    title = f"{display_label} | {mechanism} | {dtype}"
    ax.set_title(title, fontsize=7.2 if compact else 9, loc="left", fontweight="bold",
                 color=MECHANISM_COLORS[mechanism])
    ax.set_ylabel("Copy number", fontsize=7 if compact else 9)
    ax.set_xlabel("T2T-CHM13v2.0 chr15 position (Mb)", fontsize=7 if compact else 9)
    ax.tick_params(labelsize=6 if compact else 8)
    ax.grid(axis="y", color="#EEEEEE", lw=0.6)

    if not compact:
        note = structural_row.get("evidence_basis", "")
        ax.text(
            0.995, 0.02, note,
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=6.5, color="#555555"
        )


def render_supplementary_cn_profiles(
    sample_files: dict[str, dict[str, Path | None]],
    structural_by_sample: dict[str, dict[str, Any]],
    outdir: Path,
) -> None:
    """Write one CN profile per participant plus an all-cohort contact sheet."""
    profile_dir = outdir / "supplementary" / "chr15_copy_number_profiles"
    profile_dir.mkdir(parents=True, exist_ok=True)
    labels = sample_display_labels()

    # One publication-quality PDF/PNG per participant.
    for sample_id, _clinical, mechanism in sorted_cohort():
        fig, ax = plt.subplots(figsize=(10.5, 3.4))
        _plot_cn_profile_axis(
            ax,
            sample_id,
            labels[sample_id],
            mechanism,
            sample_files[sample_id].get("cn_track"),
            structural_by_sample[sample_id],
            compact=False,
        )
        fig.tight_layout()
        base = profile_dir / f"{labels[sample_id]}_{sample_id}_chr15_CN"
        fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
        fig.savefig(base.with_suffix(".png"), dpi=PUBLICATION_DPI, bbox_inches="tight")
        plt.close(fig)

    # Cohort overview with identical axes for direct depth/CN comparison.
    cohort = sorted_cohort()
    ncol = 2
    nrow = math.ceil(len(cohort) / ncol)
    fig, axes = plt.subplots(
        nrow, ncol,
        figsize=(15, 2.55 * nrow),
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    for ax in axes.flat:
        ax.axis("off")

    for ax, (sample_id, _clinical, mechanism) in zip(axes.flat, cohort):
        ax.axis("on")
        _plot_cn_profile_axis(
            ax,
            sample_id,
            labels[sample_id],
            mechanism,
            sample_files[sample_id].get("cn_track"),
            structural_by_sample[sample_id],
            compact=True,
        )

    fig.suptitle(
        "Supplementary Figure — chr15 copy-number profiles and deletion classes",
        fontsize=14, fontweight="bold", y=0.997
    )
    fig.tight_layout(rect=[0, 0, 1, 0.992])
    base = outdir / "supplementary" / "Supplementary_chr15_CN_all_participants"
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=PUBLICATION_DPI, bbox_inches="tight")
    plt.close(fig)


def classify_beta(beta: float | None, maternal_threshold: float, paternal_threshold: float) -> str:
    if beta is None or not np.isfinite(beta):
        return "?"
    if beta >= maternal_threshold:
        return "M"
    if beta <= paternal_threshold:
        return "P"
    return "I"


def exact_permutation_median_difference(a: list[float], b: list[float]) -> tuple[float, float]:
    """Exact two-sided permutation test using participants as independent units."""
    a = [float(x) for x in a if np.isfinite(x)]
    b = [float(x) for x in b if np.isfinite(x)]
    observed = float(np.median(a) - np.median(b))
    pooled = np.asarray(a + b, dtype=float)
    n_a = len(a)
    if n_a == 0 or len(b) == 0:
        return observed, float("nan")
    deltas = []
    indices = range(len(pooled))
    for chosen in itertools.combinations(indices, n_a):
        mask = np.zeros(len(pooled), dtype=bool)
        mask[list(chosen)] = True
        deltas.append(float(np.median(pooled[mask]) - np.median(pooled[~mask])))
    deltas = np.asarray(deltas)
    p = (np.sum(np.abs(deltas) >= abs(observed) - 1e-12)) / len(deltas)
    return observed, float(p)


def bootstrap_sample_median_difference(a: list[float], b: list[float]) -> tuple[float, float]:
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    if len(a) == 0 or len(b) == 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    boot = np.empty(SAMPLE_BOOTSTRAP_REPLICATES, dtype=float)
    for i in range(SAMPLE_BOOTSTRAP_REPLICATES):
        aa = rng.choice(a, size=len(a), replace=True)
        bb = rng.choice(b, size=len(b), replace=True)
        boot[i] = np.median(aa) - np.median(bb)
    lo, hi = np.quantile(boot, [0.025, 0.975])
    return float(lo), float(hi)


def build_sample_level_inference(panel_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    pws, asd = [], []
    for row in panel_rows:
        mechanism = row["molecular_mechanism"]
        if mechanism == "PWS-DEL" and row.get("allele_1_status") == "observed":
            beta = safe_float(row.get("allele_1_mean_methylation"))
            if beta is not None: pws.append(beta)
        if mechanism == "AS-DEL" and row.get("allele_2_status") == "observed":
            beta = safe_float(row.get("allele_2_mean_methylation"))
            if beta is not None: asd.append(beta)
    delta, p = exact_permutation_median_difference(pws, asd)
    lo, hi = bootstrap_sample_median_difference(pws, asd)
    return [{
        "comparison": "PWS-DEL retained maternal-like vs AS-DEL retained paternal-like",
        "independent_unit": "participant",
        "n_PWS_DEL": len(pws),
        "n_AS_DEL": len(asd),
        "median_PWS_DEL": fmt(np.median(pws) if pws else None),
        "median_AS_DEL": fmt(np.median(asd) if asd else None),
        "delta_median": fmt(delta),
        "bootstrap_95CI_low": fmt(lo),
        "bootstrap_95CI_high": fmt(hi),
        "exact_permutation_p_two_sided": fmt(p, 6),
        "note": "CpGs/molecules are nested technical observations; participants define n.",
    }]


def build_threshold_sensitivity(panel_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for mt, pt in THRESHOLD_SENSITIVITY:
        concordant = total = 0
        group_counts = defaultdict(lambda: [0, 0])
        for row in panel_rows:
            mech = row["molecular_mechanism"]
            codes = []
            for idx in (1, 2):
                status = row.get(f"allele_{idx}_status", "missing")
                if status == "deleted": codes.append("absent")
                elif status != "observed": codes.append("?")
                else: codes.append(classify_beta(safe_float(row.get(f"allele_{idx}_mean_methylation")), mt, pt))
            observed = normalize_state_pair((codes[0], codes[1]), mech)
            expected = normalize_state_pair(GROUP_EXPECTED_STATE_CODES[mech], mech)
            ok = observed == expected
            total += 1; concordant += int(ok)
            group_counts[mech][1] += 1; group_counts[mech][0] += int(ok)
        rows.append({
            "maternal_threshold": mt,
            "paternal_threshold": pt,
            "cohort_concordant_n": concordant,
            "cohort_n": total,
            "cohort_concordance_percent": fmt(100 * concordant / total if total else None, 1),
            **{f"{g}_concordance": f"{v[0]}/{v[1]}" for g, v in group_counts.items()},
        })
    return rows


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def build_input_manifest(sample_files: dict[str, dict[str, Path | None]], metadata_path: Path) -> list[dict[str, Any]]:
    rows = []
    seen = set()
    candidates = [("metadata", metadata_path)]
    for sample_id, files in sample_files.items():
        for role, path in files.items():
            if path is not None:
                candidates.append((f"{sample_id}:{role}", Path(path)))
    for role, path in candidates:
        if not path.exists():
            continue
        key = str(path.resolve())
        if (role, key) in seen: continue
        seen.add((role, key))
        size = path.stat().st_size
        do_hash = HASH_LARGE_INPUTS or size <= HASH_MAX_BYTES
        rows.append({
            "role": role, "path": str(path.resolve()), "size_bytes": size,
            "mtime_epoch": fmt(path.stat().st_mtime, 3),
            "sha256": sha256_file(path) if do_hash else "",
            "hash_status": "sha256" if do_hash else "skipped_large_file",
        })
    return rows


def command_version(command: list[str]) -> str:
    try:
        result = subprocess.run(command, text=True, capture_output=True, timeout=15)
        text = (result.stdout or result.stderr).strip().splitlines()
        return text[0] if text else "unknown"
    except Exception as exc:
        return f"unavailable: {exc}"


def software_versions() -> dict[str, Any]:
    packages = {}
    for package in ["numpy", "matplotlib", "pysam", "pandas", "scipy", "seaborn", "modbamtools"]:
        try: packages[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError: packages[package] = "not installed"
    git_commit = command_version(["git", "-C", str(PROJECT_ROOT), "rev-parse", "HEAD"])
    return {
        "python": sys.version.split()[0], "platform": platform.platform(), "packages": packages,
        "samtools": command_version(["samtools", "--version"]),
        "modbamtools": command_version([str(MODBAMTOOLS_BIN), "--version"]),
        "git_commit": git_commit,
        "reference_build": "T2T-CHM13v2.0",
    }


def write_tsv_gz(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows: writer.writerow({k: row.get(k, "") for k in fields})


def extract_modbam_source_data(sample_id: str, mechanism: str, bam_path: Path, region: str) -> list[dict[str, Any]]:
    """Long-format MM/ML source data at the displayed IC window."""
    try:
        import pysam
    except ImportError:
        return []
    chrom, start, end = parse_region(region)
    rows = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate or read.mapping_quality < MODBAM_MIN_MAPQ:
                continue
            mods = read.modified_bases
            if not mods:
                continue
            q_to_r = {q: r for q, r in read.get_aligned_pairs(matches_only=False) if q is not None and r is not None}
            hp = read.get_tag("HP") if read.has_tag("HP") else ""
            for key, calls in mods.items():
                canonical, strand, modification = key
                for qpos, qual in calls:
                    rpos = q_to_r.get(qpos)
                    if rpos is None or not (start <= rpos < end):
                        continue
                    rows.append({
                        "sample_id": sample_id, "molecular_mechanism": mechanism,
                        "read_id": read.query_name, "HP": hp, "MAPQ": read.mapping_quality,
                        "read_strand": "-" if read.is_reverse else "+",
                        "reference_position_0based": rpos,
                        "reference_position_1based": rpos + 1,
                        "canonical_base": canonical, "modification": modification, "mod_strand": strand,
                        "methylation_probability": qual / 255.0 if qual is not None and qual >= 0 else "",
                    })
    return rows


def build_preflight_qc(
    sample_files: dict[str, dict[str, Path | None]],
    structural: dict[str, dict[str, Any]],
    stats_by_sample: dict[str, dict[str, BedStats]] | None = None,
) -> list[dict[str, Any]]:
    """Scientific preflight that distinguishes usable lower depth from failure.

    Only missing or non-estimable methylation evidence can produce a support
    warning. An estimate below the higher-coverage reference remains a PASS
    when it meets the prespecified state-assignment minimum.
    """
    rows = []
    for sample_id, _clinical, mechanism in sorted_cohort():
        files = sample_files[sample_id]

        # Core files. Haplotype BEDs are no longer hard-required because an
        # MM/ML+HP ModBAM can directly quantify haplotypes when BED estimates
        # are absent or non-estimable.
        required = ["bam", "combined_bed", "modbam"]
        missing = [
            k for k in required
            if files.get(k) is None or not Path(files[k]).exists()
        ]

        mod_tags = hp_tags = None
        if files.get("modbam") and Path(files["modbam"]).exists():
            try:
                tag_status = inspect_modbam_tags(
                    Path(files["modbam"]), MODBAM_PLOT_REGION
                )
                mod_tags = tag_status["modified_base_tags"]
                hp_tags = tag_status["hp_tags"]
            except Exception:
                mod_tags, hp_tags = False, False

        s = structural[sample_id]
        issues: list[str] = []
        support_notes: list[str] = []
        status = "PASS"

        if missing:
            issues.append("missing required: " + ",".join(missing))
            status = "FAIL"

        if mod_tags is not True:
            issues.append("MM/ML not confirmed")
            status = "FAIL"

        if mechanism in {"Control", "Disease control", "PWS-mUPD"} and hp_tags is not True:
            issues.append("HP tags not confirmed")
            status = "FAIL"

        if mechanism in {"PWS-DEL", "AS-DEL"}:
            if s["ic_deletion_status"] == "not_confirmed":
                issues.append(
                    "expected IC deletion not confirmed by HiFiCNV dosage or pbsv"
                )
                status = "FAIL"
            elif s["ic_deletion_status"] == "unavailable":
                issues.append("structural/copy-number confirmation unavailable")
                if status == "PASS":
                    status = "WARN"
        else:
            if s.get("cnv_support") is True or s.get("sv_support") is True:
                issues.append("unexpected chr15 IC deletion evidence in non-deletion sample")
                status = "FAIL"

        # Depth-aware methylation support. Lower coverage is descriptive when
        # the estimate remains state-estimable; it is not a failed sample.
        hap_support = {}
        if stats_by_sample and sample_id in stats_by_sample:
            st = stats_by_sample[sample_id]
            labels = ["combined_fallback"] if mechanism in {"PWS-DEL", "AS-DEL"} else ["hap1", "hap2"]
            for label in labels:
                bst = st[label]
                hap_support[label] = _support_category(bst)
                if bst.mean_methylation is None:
                    issues.append(f"{label}: no methylation estimate")
                    if status == "PASS":
                        status = "WARN"
                elif not bst.estimable:
                    issues.append(
                        f"{label}: insufficient evidence for state assignment "
                        f"(CpGs={bst.n_cpgs}, mean_cov={fmt(bst.mean_coverage,1)}, "
                        f"source={bst.data_source})"
                    )
                    if status == "PASS":
                        status = "WARN"
                elif not bst.sufficient:
                    support_notes.append(
                        f"{label}: state-estimable below higher-coverage reference "
                        f"(CpGs={bst.n_cpgs}, mean_cov={fmt(bst.mean_coverage,1)}, "
                        f"source={bst.data_source})"
                    )

        rows.append({
            "sample_id": sample_id,
            "molecular_mechanism": mechanism,
            "status": status,
            "missing_required_files": ";".join(missing),
            "MM_ML_tags": mod_tags,
            "HP_tags": hp_tags,
            "ic_deletion_status": s["ic_deletion_status"],
            "structural_evidence_basis": s.get("evidence_basis", ""),
            "deletion_type": s.get("deletion_type", ""),
            "sv_note": s["sv_note"],
            "cnv_note": s["cnv_note"],
            "haplotype_support": ";".join(f"{k}={v}" for k, v in hap_support.items()),
            "support_notes": "; ".join(support_notes),
            "issues": "; ".join(issues),
        })
    return rows


def enforce_preflight(rows: list[dict[str, Any]]) -> None:
    fails = [r for r in rows if r["status"] == "FAIL"]
    if STRICT_SCIENTIFIC_QC and fails:
        details = "\n".join(f"- {r['sample_id']}: {r['issues']}" for r in fails)
        raise RuntimeError("Scientific preflight QC failed. Inspect Figure1_preflight_QC.tsv:\n" + details)

# ---------------------------------------------------------------------------
# IC methylation assignment
# ---------------------------------------------------------------------------

def _support_category(stats: BedStats) -> str:
    if stats.mean_methylation is None:
        return "missing"
    if stats.sufficient:
        return "higher_coverage"
    if stats.estimable:
        return "state_estimable"
    return "insufficient"


def build_assignments(
    sample_files: dict[str, dict[str, Path | None]],
    structural_by_sample: dict[str, dict[str, Any]],
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, dict[str, BedStats]],
    ParentalReferenceModel,
    list[dict[str, Any]],
]:
    """Build methylation assignments using control-calibrated references.

    The function is intentionally two-pass. First, all raw/bed-rescued IC
    methylation statistics are collected without imposing parental labels.
    Second, unaffected controls define empirical maternal-like and paternal-like
    centroids, after which every allele/haplotype is classified relative to those
    centroids. Controls themselves use a leave-one-control-out reference when
    possible.
    """
    assignment_rows: list[dict[str, Any]] = []
    matrix_rows: list[dict[str, Any]] = []
    stats_by_sample: dict[str, dict[str, BedStats]] = {}

    # Pass 1: collect quantitative methylation summaries only.
    for sample_id, _clinical, mechanism in sorted_cohort():
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
        stats_by_sample[sample_id] = complete_haplotype_stats_from_modbam(
            files, stats, mechanism
        )

    global_reference = estimate_parental_reference_model(stats_by_sample)
    reference_rows = build_parental_reference_rows(
        stats_by_sample, global_reference
    )

    # Pass 2: classify each quantitative estimate against empirical references.
    for sample_id, clinical, mechanism in sorted_cohort():
        stats = stats_by_sample[sample_id]
        sample_reference = parental_reference_for_sample(
            sample_id, mechanism, stats_by_sample, global_reference
        )
        reference_scope = (
            "leave-one-control-out"
            if mechanism == "Control"
            and CONTROL_REFERENCE_LEAVE_ONE_OUT
            and sample_id not in sample_reference.control_sample_ids
            else "global"
        )

        for label, bed_stats in stats.items():
            patt = methylation_pattern(bed_stats, sample_reference)
            metrics = reference_relative_metrics(bed_stats, sample_reference)
            matrix_rows.append({
                "sample_id": sample_id,
                "molecular_mechanism": mechanism,
                "haplotype_or_source": label,
                "mean_methylation": fmt(bed_stats.mean_methylation),
                "ci_low": fmt(bed_stats.ci_low),
                "ci_high": fmt(bed_stats.ci_high),
                "n_CpGs": bed_stats.n_cpgs,
                "n_molecules": bed_stats.n_molecules if bed_stats.n_molecules is not None else "",
                "mean_coverage": fmt(bed_stats.mean_coverage),
                "pattern": patt,
                "pattern_short": pattern_short(patt),
                "support_category": _support_category(bed_stats),
                "coverage_status": (
                    "higher_coverage"
                    if bed_stats.sufficient
                    else "state_estimable"
                    if bed_stats.estimable
                    else "insufficient"
                ),
                "data_source": bed_stats.data_source,
                "classification_method": "control-calibrated reference + descriptive CI",
                "reference_scope": reference_scope,
                "maternal_reference_beta": fmt(sample_reference.maternal_reference),
                "paternal_reference_beta": fmt(sample_reference.paternal_reference),
                "decision_boundary_beta": fmt(sample_reference.decision_boundary),
                "distance_to_maternal_reference": fmt(metrics["distance_to_maternal"]),
                "distance_to_paternal_reference": fmt(metrics["distance_to_paternal"]),
                "nearest_parental_reference": metrics["nearest_reference"],
                "reference_margin": fmt(metrics["reference_margin"]),
                "parental_axis_score": fmt(metrics["parental_axis_score"]),
            })

        if mechanism in {"PWS-DEL", "AS-DEL"}:
            rows_for_sample = [
                ("combined_fallback", stats["combined_fallback"], stats["combined_fallback"].data_source)
            ]
            expected = "maternal-pattern" if mechanism == "PWS-DEL" else "paternal-pattern"
        else:
            rows_for_sample = [
                ("hap1", stats["hap1"], stats["hap1"].data_source),
                ("hap2", stats["hap2"], stats["hap2"].data_source),
            ]
            expected = (
                "both maternal-pattern"
                if mechanism == "PWS-mUPD"
                else "one maternal-pattern and one paternal-pattern"
            )

        for label, bed_stats, source in rows_for_sample:
            patt = methylation_pattern(bed_stats, sample_reference)
            metrics = reference_relative_metrics(bed_stats, sample_reference)
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

            assignment_rows.append({
                "sample_id": sample_id,
                "clinical_diagnosis": clinical,
                "molecular_mechanism": mechanism,
                "haplotype_label": label,
                "source": source,
                "mean_methylation_at_PWS_IC": fmt(bed_stats.mean_methylation),
                "bootstrap_CI_low": fmt(bed_stats.ci_low),
                "bootstrap_CI_high": fmt(bed_stats.ci_high),
                "n_CpGs_at_PWS_IC": bed_stats.n_cpgs,
                "n_molecules_at_PWS_IC": bed_stats.n_molecules if bed_stats.n_molecules is not None else "",
                "mean_coverage_at_PWS_IC": fmt(bed_stats.mean_coverage),
                "support_category": _support_category(bed_stats),
                "methylation_pattern": patt,
                "parental_assignment": parental_assignment,
                "nearest_parental_reference": metrics["nearest_reference"],
                "maternal_reference_beta": fmt(sample_reference.maternal_reference),
                "paternal_reference_beta": fmt(sample_reference.paternal_reference),
                "decision_boundary_beta": fmt(sample_reference.decision_boundary),
                "distance_to_maternal_reference": fmt(metrics["distance_to_maternal"]),
                "distance_to_paternal_reference": fmt(metrics["distance_to_paternal"]),
                "reference_margin": fmt(metrics["reference_margin"]),
                "parental_axis_score": fmt(metrics["parental_axis_score"]),
                "expected_pattern": expected,
                "structural_ic_deletion_status": structural_by_sample[sample_id]["ic_deletion_status"],
                "note": (
                    "Primary parental-like state is learned from unaffected-control "
                    "methylation centroids; fixed 0.85/0.15 thresholds are used "
                    "only for Supplementary sensitivity analysis. State-estimable "
                    "lower-coverage haplotypes are retained with continuous support "
                    "metrics and are not treated as failed samples."
                ),
            })

    return (
        assignment_rows,
        matrix_rows,
        stats_by_sample,
        global_reference,
        reference_rows,
    )


def _state_cell(status: str, row: dict[str, Any] | None = None) -> dict[str, Any]:
    if status != "observed":
        label = {"deleted": "absent", "missing": "?", "structural_unresolved": "?"}.get(status, "?")
        return {
            "source": status, "mean_methylation": None, "ci_low": None, "ci_high": None,
            "pattern": status, "pattern_short": label, "n_CpGs": 0, "n_molecules": None, "mean_coverage": None,
            "coverage_status": status,
            "nearest_parental_reference": "", "maternal_reference_beta": None, "paternal_reference_beta": None,
            "decision_boundary_beta": None, "distance_to_maternal_reference": None,
            "distance_to_paternal_reference": None, "reference_margin": None, "parental_axis_score": None,
            "status": status, "is_absent": status == "deleted",
        }
    if row is None:
        return _state_cell("missing")
    value = safe_float(row.get("mean_methylation"))
    if value is None:
        return _state_cell("missing")
    patt = row.get("pattern", "missing")
    return {
        "source": row.get("data_source") or row.get("haplotype_or_source", ""), "mean_methylation": value,
        "ci_low": safe_float(row.get("ci_low")), "ci_high": safe_float(row.get("ci_high")),
        "pattern": patt, "pattern_short": row.get("pattern_short") or pattern_short(patt),
        "n_CpGs": int(row.get("n_CpGs", 0) or 0),
        "n_molecules": int(row.get("n_molecules", 0) or 0) if str(row.get("n_molecules", "")).strip() else None,
        "mean_coverage": safe_float(row.get("mean_coverage")), "coverage_status": row.get("coverage_status", ""),
        "nearest_parental_reference": row.get("nearest_parental_reference", ""),
        "maternal_reference_beta": safe_float(row.get("maternal_reference_beta")),
        "paternal_reference_beta": safe_float(row.get("paternal_reference_beta")),
        "decision_boundary_beta": safe_float(row.get("decision_boundary_beta")),
        "distance_to_maternal_reference": safe_float(row.get("distance_to_maternal_reference")),
        "distance_to_paternal_reference": safe_float(row.get("distance_to_paternal_reference")),
        "reference_margin": safe_float(row.get("reference_margin")),
        "parental_axis_score": safe_float(row.get("parental_axis_score")),
        "status": "observed", "is_absent": False,
    }


def _row_to_cell(row: dict[str, Any] | None) -> dict[str, Any]:
    """Backward-compatible wrapper: missing observations are never called deleted."""
    return _state_cell("observed", row) if row is not None else _state_cell("missing")


def build_physical_allele_rows(
    matrix_rows: list[dict[str, Any]],
    structural_by_sample: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    by_sample: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)
    for row in matrix_rows:
        by_sample[row["sample_id"]][row["haplotype_or_source"]] = row
    labels = sample_display_labels(); rows = []
    for sample_id, _clinical, mechanism in sorted_cohort():
        sample_rows = by_sample[sample_id]
        structural_status = structural_by_sample[sample_id]["ic_deletion_status"]
        deleted_status = "deleted" if structural_status == "confirmed" else "structural_unresolved"
        if mechanism == "PWS-DEL":
            allele_1 = _row_to_cell(sample_rows.get("combined_fallback")); allele_2 = _state_cell(deleted_status)
        elif mechanism == "AS-DEL":
            allele_1 = _state_cell(deleted_status); allele_2 = _row_to_cell(sample_rows.get("combined_fallback"))
        else:
            allele_1 = _row_to_cell(sample_rows.get("hap1")); allele_2 = _row_to_cell(sample_rows.get("hap2"))
        out = {"sample_id": sample_id, "display_label": labels[sample_id], "molecular_mechanism": mechanism,
               "structural_ic_deletion_status": structural_status}
        for idx, cell in [(1, allele_1), (2, allele_2)]:
            out.update({
                f"allele_{idx}_source": cell["source"], f"allele_{idx}_status": cell["status"],
                f"allele_{idx}_mean_methylation": fmt(cell["mean_methylation"]),
                f"allele_{idx}_ci_low": fmt(cell["ci_low"]), f"allele_{idx}_ci_high": fmt(cell["ci_high"]),
                f"allele_{idx}_pattern": cell["pattern"], f"allele_{idx}_pattern_short": cell["pattern_short"],
                f"allele_{idx}_n_CpGs": cell["n_CpGs"], f"allele_{idx}_n_molecules": cell.get("n_molecules", ""),
                f"allele_{idx}_mean_coverage": fmt(cell["mean_coverage"]),
                f"allele_{idx}_coverage_status": cell["coverage_status"],
                f"allele_{idx}_nearest_parental_reference": cell.get("nearest_parental_reference", ""),
                f"allele_{idx}_maternal_reference_beta": fmt(cell.get("maternal_reference_beta")),
                f"allele_{idx}_paternal_reference_beta": fmt(cell.get("paternal_reference_beta")),
                f"allele_{idx}_decision_boundary_beta": fmt(cell.get("decision_boundary_beta")),
                f"allele_{idx}_distance_to_maternal_reference": fmt(cell.get("distance_to_maternal_reference")),
                f"allele_{idx}_distance_to_paternal_reference": fmt(cell.get("distance_to_paternal_reference")),
                f"allele_{idx}_reference_margin": fmt(cell.get("reference_margin")),
                f"allele_{idx}_parental_axis_score": fmt(cell.get("parental_axis_score")),
                f"allele_{idx}_is_absent": str(cell["is_absent"]),
            })
        rows.append(out)
    return rows


# ---------------------------------------------------------------------------
# Deletion-span methylation-profile classification
# ---------------------------------------------------------------------------

def read_bigwig_region_values(
    path: Path | None,
    start: int,
    end: int,
) -> dict[int, float]:
    """Return position-level methylation beta values from a pb-CpG BigWig."""
    if path is None or not Path(path).exists() or start >= end:
        return {}
    try:
        import pyBigWig
    except ImportError as exc:
        raise RuntimeError(
            "Deletion-span classification requires pyBigWig."
        ) from exc

    values: dict[int, float] = {}
    with pyBigWig.open(str(path)) as bw:
        intervals = bw.intervals(CHROM, int(start), int(end)) or []
    for interval_start, _interval_end, raw_value in intervals:
        value = safe_float(raw_value)
        if value is None:
            continue
        # pb-CpG-tools BigWigs in this project store percent methylation,
        # including exact 0% and 1% values that must not be mistaken for beta.
        beta = value / 100.0
        if 0.0 <= beta <= 1.0:
            values[int(interval_start)] = float(beta)
    return values


def _slice_profile(
    values: dict[int, float],
    start: int,
    end: int,
) -> dict[int, float]:
    return {pos: value for pos, value in values.items() if start <= pos < end}


def _profile_rmse_score(
    observed: dict[int, float],
    maternal_reference: dict[int, float],
    paternal_reference: dict[int, float],
    informative_positions: set[int],
) -> dict[str, Any]:
    """Classify one profile on a continuous -1 (P) to +1 (M) axis."""
    shared = sorted(set(observed).intersection(informative_positions))
    if len(shared) < DELETION_PROFILE_MIN_SHARED_CPGS:
        return {
            "n_shared_informative_CpGs": len(shared),
            "rmse_to_maternal_profile": None,
            "rmse_to_paternal_profile": None,
            "parental_profile_score": None,
            "profile_parental_class": "insufficient",
        }

    obs = np.asarray([observed[pos] for pos in shared], dtype=float)
    maternal = np.asarray([maternal_reference[pos] for pos in shared], dtype=float)
    paternal = np.asarray([paternal_reference[pos] for pos in shared], dtype=float)
    rmse_m = float(np.sqrt(np.mean(np.square(obs - maternal))))
    rmse_p = float(np.sqrt(np.mean(np.square(obs - paternal))))
    denominator = rmse_m + rmse_p
    score = float((rmse_p - rmse_m) / denominator) if denominator > 0 else 0.0
    if score >= DELETION_PROFILE_MIN_SCORE_MAGNITUDE:
        profile_class = "maternal-like"
    elif score <= -DELETION_PROFILE_MIN_SCORE_MAGNITUDE:
        profile_class = "paternal-like"
    else:
        profile_class = "uncertain"
    return {
        "n_shared_informative_CpGs": len(shared),
        "rmse_to_maternal_profile": rmse_m,
        "rmse_to_paternal_profile": rmse_p,
        "parental_profile_score": score,
        "profile_parental_class": profile_class,
    }


def build_deletion_profile_classification_rows(
    sample_files: dict[str, dict[str, Path | None]],
    structural_by_sample: dict[str, dict[str, Any]],
    panel_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Classify methylation tracks across each CN-defined deletion interval.

    Unaffected controls alone define the M-like and P-like position-level
    references. The control haplotype identity is anchored by its IC state;
    disease-control haplotypes are held out and used only for validation.
    """
    panel_by_sample = {row["sample_id"]: row for row in panel_rows}
    labels = sample_display_labels()
    deletion_samples = [
        (sample_id, mechanism, structural_by_sample[sample_id])
        for sample_id, _clinical, mechanism in sorted_cohort()
        if mechanism in {"PWS-DEL", "AS-DEL"}
        and structural_by_sample[sample_id].get("ic_deletion_status") == "confirmed"
    ]
    interval_by_sample: dict[str, tuple[int, int]] = {}
    for sample_id, _mechanism, row in deletion_samples:
        start = safe_int(row.get("cn_event_start"))
        end = safe_int(row.get("cn_event_end"))
        if start is None or end is None:
            raise ValueError(
                f"Confirmed chr15 deletion for {sample_id} has invalid coordinates: "
                f"start={row.get('cn_event_start')!r}, end={row.get('cn_event_end')!r}"
            )
        if start < 0 or end <= start:
            raise ValueError(
                f"Confirmed chr15 deletion for {sample_id} has a non-positive interval: "
                f"start={start}, end={end}"
            )
        interval_by_sample[sample_id] = (start, end)
    intervals = list(interval_by_sample.values())
    if not intervals:
        return []
    union_start = min(start for start, _end in intervals)
    union_end = max(end for _start, end in intervals)
    common_start = max(start for start, _end in intervals)
    common_end = min(end for _start, end in intervals)

    control_profiles: dict[str, list[dict[int, float]]] = {"M": [], "P": []}
    for sample_id, _clinical, mechanism in sorted_cohort():
        if mechanism != "Control":
            continue
        panel_row = panel_by_sample.get(sample_id, {})
        for allele_index, key in ((1, "hap1_bw"), (2, "hap2_bw")):
            state = str(panel_row.get(f"allele_{allele_index}_pattern_short", ""))
            if state not in {"M", "P"}:
                continue
            profile = read_bigwig_region_values(
                sample_files[sample_id].get(key), union_start, union_end
            )
            if profile:
                control_profiles[state].append(profile)

    if len(control_profiles["M"]) < 2 or len(control_profiles["P"]) < 2:
        raise RuntimeError(
            "Deletion-span classification requires two unaffected-control "
            "M-like and two P-like haplotype profiles."
        )

    def make_reference(start: int, end: int) -> tuple[dict[int, float], dict[int, float], set[int]]:
        maternal_profiles = [_slice_profile(profile, start, end) for profile in control_profiles["M"]]
        paternal_profiles = [_slice_profile(profile, start, end) for profile in control_profiles["P"]]
        shared = set(maternal_profiles[0])
        for profile in maternal_profiles[1:] + paternal_profiles:
            shared.intersection_update(profile)
        maternal_reference = {
            pos: float(np.median([profile[pos] for profile in maternal_profiles]))
            for pos in shared
        }
        paternal_reference = {
            pos: float(np.median([profile[pos] for profile in paternal_profiles]))
            for pos in shared
        }
        informative = {
            pos
            for pos in shared
            if abs(maternal_reference[pos] - paternal_reference[pos])
            >= DELETION_PROFILE_MIN_REFERENCE_DELTA
            and (max(profile[pos] for profile in maternal_profiles)
                 - min(profile[pos] for profile in maternal_profiles))
            <= DELETION_PROFILE_MAX_WITHIN_STATE_RANGE
            and (max(profile[pos] for profile in paternal_profiles)
                 - min(profile[pos] for profile in paternal_profiles))
            <= DELETION_PROFILE_MAX_WITHIN_STATE_RANGE
        }
        return maternal_reference, paternal_reference, informative

    rows: list[dict[str, Any]] = []
    for sample_id, mechanism, structural_row in deletion_samples:
        start, end = interval_by_sample[sample_id]
        maternal_reference, paternal_reference, informative = make_reference(start, end)
        observed = read_bigwig_region_values(
            sample_files[sample_id].get("combined_bw"), start, end
        )
        metrics = _profile_rmse_score(
            observed, maternal_reference, paternal_reference, informative
        )
        panel_row = panel_by_sample.get(sample_id, {})
        observed_index = 1 if mechanism == "PWS-DEL" else 2
        ic_short = str(panel_row.get(f"allele_{observed_index}_pattern_short", "?"))
        ic_class = {"M": "maternal-like", "P": "paternal-like"}.get(ic_short, "uncertain")
        expected_class = "maternal-like" if mechanism == "PWS-DEL" else "paternal-like"
        profile_class = str(metrics["profile_parental_class"])
        if profile_class in {"maternal-like", "paternal-like"}:
            integrated = (
                f"concordant {profile_class}"
                if profile_class == ic_class
                else f"discordant: profile {profile_class}, IC {ic_class}"
            )
        else:
            integrated = f"profile {profile_class}; IC {ic_class}"
        rows.append({
            "sample_id": sample_id,
            "display_label": labels[sample_id],
            "molecular_mechanism": mechanism,
            "track_role": "retained deletion chromosome",
            "haplotype_label": "retained",
            "evaluation_interval": "sample-specific CN deletion",
            "evaluation_start": start,
            "evaluation_end": end,
            "deletion_size_mb": fmt((end - start) / 1e6),
            "n_control_informative_CpGs": len(informative),
            **{key: fmt(value) if isinstance(value, float) else value for key, value in metrics.items()},
            "IC_parental_class": ic_class,
            "expected_parental_class": expected_class,
            "integrated_classification": integrated,
            "profile_matches_expected": str(profile_class == expected_class),
        })

    # Controls and DiGeorge disease controls are evaluated on the interval
    # shared by every deletion. They provide separated haplotype benchmarks but
    # do not contribute to the empirical reference fit unless they are controls.
    common_maternal, common_paternal, common_informative = make_reference(
        common_start, common_end
    )
    for sample_id, _clinical, mechanism in sorted_cohort():
        if mechanism not in {"Control", "Disease control"}:
            continue
        panel_row = panel_by_sample.get(sample_id, {})
        for allele_index, key in ((1, "hap1_bw"), (2, "hap2_bw")):
            observed = read_bigwig_region_values(
                sample_files[sample_id].get(key), common_start, common_end
            )
            metrics = _profile_rmse_score(
                observed, common_maternal, common_paternal, common_informative
            )
            ic_short = str(panel_row.get(f"allele_{allele_index}_pattern_short", "?"))
            ic_class = {"M": "maternal-like", "P": "paternal-like"}.get(ic_short, "uncertain")
            profile_class = str(metrics["profile_parental_class"])
            rows.append({
                "sample_id": sample_id,
                "display_label": labels[sample_id],
                "molecular_mechanism": mechanism,
                "track_role": "reference" if mechanism == "Control" else "held-out disease control",
                "haplotype_label": f"H{allele_index}",
                "evaluation_interval": "shared deletion core",
                "evaluation_start": common_start,
                "evaluation_end": common_end,
                "deletion_size_mb": "",
                "n_control_informative_CpGs": len(common_informative),
                **{key: fmt(value) if isinstance(value, float) else value for key, value in metrics.items()},
                "IC_parental_class": ic_class,
                "expected_parental_class": ic_class,
                "integrated_classification": (
                    f"concordant {profile_class}"
                    if profile_class == ic_class
                    else f"profile {profile_class}; IC {ic_class}"
                ),
                "profile_matches_expected": str(profile_class == ic_class),
            })
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


def closest_methylation_template(
    observed_betas: list[float],
    reference: ParentalReferenceModel,
) -> tuple[str, float]:
    """Continuous template-distance classifier calibrated to controls.

    Templates are no longer fixed at idealized 0 and 1. Instead they use the
    empirical paternal-like and maternal-like control centroids estimated in the
    same run. This makes template classification consistent with the primary
    parental-state model while remaining independent of the disease label.
    """
    values = np.asarray([x for x in observed_betas if np.isfinite(x)], dtype=float)
    m = reference.maternal_reference
    p = reference.paternal_reference

    if len(values) == 1:
        beta = float(values[0])
        distances = {
            "maternal-retained deletion": abs(beta - m),
            "paternal-retained deletion": abs(beta - p),
        }
    elif len(values) == 2:
        low, high = np.sort(values)
        obs = np.array([low, high], dtype=float)
        templates = {
            "canonical biparental chr15": np.array([p, m]),
            "duplicated maternal state": np.array([m, m]),
            "duplicated paternal state": np.array([p, p]),
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


def build_mechanistic_state_rows(
    panel_a_rows: list[dict[str, Any]],
    reference: ParentalReferenceModel,
) -> list[dict[str, Any]]:
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

        # Use the same per-sample reference that generated the discrete state.
        # For unaffected controls this is leave-one-control-out when available;
        # all other samples use the global unaffected-control reference.
        sample_reference = reference
        sample_m = safe_float(row.get("allele_1_maternal_reference_beta"))
        sample_p = safe_float(row.get("allele_1_paternal_reference_beta"))
        sample_b = safe_float(row.get("allele_1_decision_boundary_beta"))
        if sample_m is not None and sample_p is not None and sample_b is not None:
            sample_reference = ParentalReferenceModel(
                sample_m, sample_p, sample_b, tuple(), reference.method
            )

        predicted_template, template_distance = closest_methylation_template(values, sample_reference)
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
                "maternal_reference_beta": fmt(sample_reference.maternal_reference),
                "paternal_reference_beta": fmt(sample_reference.paternal_reference),
                "decision_boundary_beta": fmt(sample_reference.decision_boundary),
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
# Raw single-molecule ModBAM / ModBAMtools evidence
# ---------------------------------------------------------------------------

def parse_region(region: str) -> tuple[str, int, int]:
    """Parse chr:start-end into a tuple suitable for pysam."""
    match = re.fullmatch(r"([^:]+):(\d+)-(\d+)", region.replace(",", ""))
    if not match:
        raise ValueError(
            f"Invalid genomic region '{region}'. Expected e.g. "
            "chr15:22690500-22695000."
        )
    chrom, start, end = match.groups()
    start_i = int(start)
    end_i = int(end)
    if end_i <= start_i:
        raise ValueError(f"Invalid genomic interval '{region}'.")
    return chrom, start_i, end_i


def inspect_modbam_tags(
    bam_path: Path,
    region: str,
    max_reads: int = 300,
) -> dict[str, bool]:
    """
    Check whether the candidate BAM contains modified-base (MM/ML) and HP tags
    in the plotted region.

    ModBAMtools requires MM/ML tags. HP tags are additionally required when the
    plot is split by haplotype (controls and mUPD in this figure).
    """
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError(
            "pysam is required to validate ModBAM/HP tags before Figure 1A."
        ) from exc

    chrom, start, end = parse_region(region)
    found_mod = False
    found_hp = False
    observed = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        try:
            iterator = bam.fetch(chrom, start, end)
        except ValueError as exc:
            raise RuntimeError(
                f"Cannot query {region} in {bam_path}. Check the reference "
                "coordinates and BAM index."
            ) from exc

        for read in iterator:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            observed += 1
            tags = {tag for tag, _value in read.get_tags()}

            # SAM permits the historical mixed-case spellings as well.
            has_mm = "MM" in tags or "Mm" in tags
            has_ml = "ML" in tags or "Ml" in tags
            found_mod = found_mod or (has_mm and has_ml)
            found_hp = found_hp or ("HP" in tags)

            if found_mod and found_hp:
                break
            if observed >= max_reads:
                break

    return {
        "modified_base_tags": found_mod,
        "hp_tags": found_hp,
    }


def choose_representative_modbam_samples(
    support_rows: list[dict[str, Any]],
    sample_files: dict[str, dict[str, Path | None]],
) -> dict[str, dict[str, Any]]:
    """
    Select one representative sample for every cohort block in the heatmap.

    Selection is deterministic: choose the individual whose total IC depth is
    closest to the within-group median IC depth. This prevents aesthetic
    cherry-picking of the cleanest raw single-molecule panel.
    """
    representatives: dict[str, dict[str, Any]] = {}

    for mechanism in MODBAM_GROUP_ORDER:
        candidates: list[dict[str, Any]] = []

        for row in support_rows:
            if row["molecular_mechanism"] != mechanism:
                continue

            sample_id = row["sample_id"]
            modbam = sample_files.get(sample_id, {}).get("modbam")
            depth = safe_float(row.get("total_ic_depth"))

            if modbam is None or not Path(modbam).exists() or depth is None:
                continue

            candidates.append(
                {
                    "sample_id": sample_id,
                    "display_label": row["display_label"],
                    "molecular_mechanism": mechanism,
                    "total_ic_depth": depth,
                    "modbam": Path(modbam),
                }
            )

        if not candidates:
            raise RuntimeError(
                f"No usable ModBAM candidate found for {mechanism}. "
                "Check MODBAM_DIR, BAM names and IC coverage in USER CONFIGURATION."
            )

        group_depths = np.asarray(
            [candidate["total_ic_depth"] for candidate in candidates],
            dtype=float,
        )
        target = float(np.median(group_depths))

        winner = min(
            candidates,
            key=lambda candidate: (
                abs(candidate["total_ic_depth"] - target),
                candidate["sample_id"],
            ),
        )
        winner["group_median_depth"] = target
        representatives[mechanism] = winner

    return representatives


def locate_modbamtools_png(output_dir: Path, prefix: str) -> Path:
    """Resolve the PNG emitted by ModBAMtools across minor naming differences."""
    canonical = output_dir / f"{prefix}.png"
    if canonical.exists():
        return canonical

    matches = sorted(output_dir.glob(f"{prefix}*.png"))
    if not matches:
        raise FileNotFoundError(
            f"ModBAMtools completed but no PNG matching '{prefix}*.png' "
            f"was found in {output_dir}."
        )
    return matches[0]


def prepare_standardized_modbam(bam_path: Path, region: str, output_dir: Path, sample_label: str) -> Path:
    """Create a small region-only primary-alignment modBAM with fixed MAPQ filter."""
    output_dir.mkdir(parents=True, exist_ok=True)
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", sample_label)
    filtered = output_dir / f"{safe}.MAPQ{MODBAM_MIN_MAPQ}.modbam.bam"
    if filtered.exists() and Path(str(filtered) + ".bai").exists():
        return filtered
    command = ["samtools", "view", "-bh", "-q", str(MODBAM_MIN_MAPQ), "-F", "2308", "-o", str(filtered), str(bam_path), region]
    subprocess.run(command, check=True)
    subprocess.run(["samtools", "index", str(filtered)], check=True)
    return filtered


def render_internal_modbam_panel(
    bam_path: Path,
    sample_label: str,
    mechanism: str,
    region: str,
    output_path: Path,
) -> Path:
    """Render a compact IGV-like molecule panel directly from MM/ML tags."""
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("Internal ModBAM rendering requires pysam.") from exc

    chrom, start, end = parse_region(region)
    use_hap = mechanism in {"Control", "Disease control", "PWS-mUPD"}
    reads_by_group: dict[str, list[dict[str, Any]]] = defaultdict(list)
    position_probabilities_by_group: dict[str, dict[int, list[float]]] = defaultdict(
        lambda: defaultdict(list)
    )

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
                continue
            if read.mapping_quality < MODBAM_MIN_MAPQ:
                continue
            hp = None
            if read.has_tag("HP"):
                try:
                    hp = int(read.get_tag("HP"))
                except (TypeError, ValueError):
                    hp = None
            if use_hap and hp not in {1, 2}:
                # Biparental panels intentionally contain exactly two traces and
                # two molecule groups: HP1 and HP2. Untagged reads are excluded.
                continue
            group = f"HP{hp}" if use_hap else "retained"
            q_to_r = {
                query_pos: ref_pos
                for query_pos, ref_pos in read.get_aligned_pairs(matches_only=False)
                if query_pos is not None and ref_pos is not None
            }
            calls: list[tuple[int, float]] = []
            for key, modified_calls in (read.modified_bases or {}).items():
                canonical, _strand, modification = key
                if str(canonical).upper() != "C":
                    continue
                if str(modification).lower() not in {"m", "5mc", "c+m"}:
                    continue
                for query_pos, quality in modified_calls:
                    ref_pos = q_to_r.get(query_pos)
                    if ref_pos is None or quality is None or quality < 0:
                        continue
                    if start <= ref_pos < end:
                        probability = float(quality) / 255.0
                        calls.append((ref_pos, probability))
                        position_probabilities_by_group[group][ref_pos].append(probability)
            if not calls:
                continue
            reads_by_group[group].append({
                "start": max(start, int(read.reference_start)),
                "end": min(end, int(read.reference_end or end)),
                "calls": sorted(calls),
            })

    group_order = [group for group in ("HP1", "HP2", "retained") if reads_by_group.get(group)]
    if not group_order:
        raise RuntimeError(f"No plottable MM/ML calls found in {bam_path} at {region}.")

    # Cap only extremely deep tracks; deterministic evenly spaced selection
    # preserves coverage across the complete stack without visual overplotting.
    max_reads_per_group = 32
    for group in group_order:
        group_reads = sorted(reads_by_group[group], key=lambda item: (item["start"], item["end"]))
        if len(group_reads) > max_reads_per_group:
            indices = np.linspace(0, len(group_reads) - 1, max_reads_per_group).astype(int)
            group_reads = [group_reads[index] for index in indices]
        reads_by_group[group] = group_reads

    height = 1.58 if use_hap else 1.42
    fig = plt.figure(figsize=(8.0, height), constrained_layout=False)
    grid = GridSpec(2, 1, figure=fig, height_ratios=[0.72, 0.42], hspace=0.055)
    frequency_ax = fig.add_subplot(grid[0, 0])
    molecule_ax = fig.add_subplot(grid[1, 0], sharex=frequency_ax)

    trace_colors = {"HP1": "#D55E00", "HP2": "#0072B2", "retained": "#375A7F"}
    for group in group_order:
        group_values = position_probabilities_by_group[group]
        binned_probabilities: dict[int, list[float]] = defaultdict(list)
        for position, probabilities in group_values.items():
            bin_index = (position - start) // MODBAM_FREQUENCY_BIN_BP
            binned_probabilities[int(bin_index)].extend(probabilities)
        positions = np.asarray(
            [
                start + (bin_index + 0.5) * MODBAM_FREQUENCY_BIN_BP
                for bin_index in sorted(binned_probabilities)
            ],
            dtype=float,
        )
        frequency = np.asarray(
            [
                100.0 * np.mean(binned_probabilities[bin_index])
                for bin_index in sorted(binned_probabilities)
            ],
            dtype=float,
        )
        if len(positions):
            frequency_ax.plot(
                positions, frequency, color=trace_colors[group], lw=1.05,
                alpha=0.94, label=group,
            )
            frequency_ax.scatter(
                positions, frequency, s=2.5, color=trace_colors[group],
                alpha=0.70, linewidths=0,
            )
    frequency_ax.set_ylim(-4, 104)
    frequency_ax.set_yticks([0, 50, 100])
    frequency_ax.set_yticklabels(["0", "50", "100"], fontsize=4.6)
    frequency_ax.set_ylabel("%", fontsize=4.8, labelpad=2)
    frequency_ax.tick_params(axis="x", bottom=False, labelbottom=False)
    frequency_ax.tick_params(axis="y", length=1.8, width=0.5)
    frequency_ax.spines[["top", "right"]].set_visible(False)
    frequency_ax.spines[["left", "bottom"]].set_linewidth(0.45)
    frequency_ax.text(
        0.5, 1.02, "Methylation frequency", transform=frequency_ax.transAxes,
        ha="center", va="bottom", fontsize=4.8, color="#444444",
    )
    if use_hap:
        frequency_ax.legend(
            frameon=False, loc="upper right", ncol=2, fontsize=4.5,
            handlelength=1.3, columnspacing=0.7, borderaxespad=0.2,
        )

    current_y = 0
    group_centres: list[tuple[str, float]] = []
    for group_index, group in enumerate(group_order):
        group_start = current_y
        for read in reads_by_group[group]:
            molecule_ax.hlines(
                current_y, read["start"], read["end"],
                color="#97A6B2", lw=0.48, alpha=0.86, zorder=1,
            )
            call_x = [call[0] for call in read["calls"]]
            call_p = [call[1] for call in read["calls"]]
            molecule_ax.scatter(
                call_x, [current_y] * len(call_x), c=call_p,
                cmap="coolwarm", vmin=0, vmax=1, s=4.0,
                marker="|", linewidths=0.55, alpha=0.98, zorder=2,
            )
            current_y += 1
        group_centres.append((group, (group_start + current_y - 1) / 2.0))
        if group_index < len(group_order) - 1:
            molecule_ax.axhline(current_y - 0.5, color="#555555", lw=0.55, alpha=0.55)
            current_y += 1

    molecule_ax.set_xlim(start, end)
    molecule_ax.set_ylim(current_y - 0.3, -0.7)
    molecule_ax.set_yticks([centre for _group, centre in group_centres])
    molecule_ax.set_yticklabels([group for group, _centre in group_centres], fontsize=4.7)
    ticks = np.linspace(start, end, 5)
    molecule_ax.set_xticks(ticks)
    molecule_ax.set_xticklabels([f"{tick / 1e6:.3f}" for tick in ticks], fontsize=4.5)
    molecule_ax.set_xlabel(f"{chrom} position (Mb)", fontsize=4.8, labelpad=1.5)
    molecule_ax.tick_params(axis="both", length=1.8, width=0.45, pad=1.2)
    molecule_ax.spines[["top", "right"]].set_visible(False)
    molecule_ax.spines[["left", "bottom"]].set_linewidth(0.45)

    fig.subplots_adjust(left=0.060, right=0.995, top=0.955, bottom=0.21)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, facecolor="white")
    plt.close(fig)
    return output_path


def run_modbamtools_plot(
    modbamtools_bin: str,
    bam_path: Path,
    sample_label: str,
    mechanism: str,
    region: str,
    output_dir: Path,
    gtf_path: Path | None = None,
) -> Path:
    """
    Generate one publication-facing raw single-molecule methylation panel.

    Haplotype grouping is used for biparental controls and PWS-mUPD. It is not
    forced inside reciprocal deletion intervals, where the IC is biologically
    hemizygous and lack of a second HP class is expected.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    executable = shutil.which(modbamtools_bin)
    if executable is None:
        explicit = Path(modbamtools_bin).expanduser()
        if explicit.exists():
            executable = str(explicit)
        else:
            raise FileNotFoundError(
                f"Cannot find ModBAMtools executable '{modbamtools_bin}'. "
                "Edit MODBAMTOOLS_BIN in USER CONFIGURATION."
            )

    tags = inspect_modbam_tags(bam_path, region)
    if not tags["modified_base_tags"]:
        raise RuntimeError(
            f"{bam_path} does not contain MM/ML modified-base tags in {region}. "
            "Point MODBAM_DIR to phased BAMs that retain native modification "
            "tags; pb-CpG-tools BED files alone are not ModBAMtools input."
        )

    use_hap = mechanism in {"Control", "Disease control", "PWS-mUPD"}
    if use_hap and not tags["hp_tags"]:
        raise RuntimeError(
            f"{bam_path} has MM/ML tags but no HP tags were observed in {region}. "
            f"Figure 1A requires haplotype grouping for {mechanism}. Use a "
            "HiPhase-tagged modBAM."
        )

    filtered_bam = prepare_standardized_modbam(
        bam_path, region, output_dir / "filtered_inputs", sample_label
    )

    safe_mechanism = mechanism.lower().replace("-", "_").replace(" ", "_")
    prefix = f"Figure1A_{safe_mechanism}_{sample_label}"

    if not USE_EXTERNAL_MODBAMTOOLS:
        return render_internal_modbam_panel(
            filtered_bam,
            sample_label,
            mechanism,
            region,
            output_dir / f"{prefix}.png",
        )

    command = [
        executable,
        "plot",
        "--region", region,
        "--out", str(output_dir),
        "--prefix", prefix,
        "--samples", sample_label,
        "--fmt", "png",
        "--width", str(MODBAM_PLOT_WIDTH),
    ]

    if use_hap:
        command.append("--hap")

    if gtf_path is not None:
        gtf_path = Path(gtf_path).expanduser()
        if not gtf_path.exists():
            raise FileNotFoundError(f"ModBAMtools GTF does not exist: {gtf_path}")
        tbi = Path(str(gtf_path) + ".tbi")
        if not tbi.exists():
            raise FileNotFoundError(
                f"{gtf_path} has no tabix index ({tbi}). ModBAMtools expects a "
                "sorted, bgzip-compressed and tabix-indexed GTF."
            )
        command.extend(["--gtf", str(gtf_path), "--track-titles", "Genes"])

    command.append(str(filtered_bam))

    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            "ModBAMtools failed.\n\n"
            f"Command:\n{' '.join(command)}\n\n"
            f"STDOUT:\n{result.stdout}\n\n"
            f"STDERR:\n{result.stderr}"
        )

    return locate_modbamtools_png(output_dir, prefix)


def build_modbamtools_panels(
    representatives: dict[str, dict[str, Any]],
    outdir: Path,
    modbamtools_bin: str,
    region: str,
    gtf_path: Path | None,
) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    """Generate one raw panel per cohort block and a provenance table."""
    panel_dir = outdir / "modbamtools"
    panel_dir.mkdir(parents=True, exist_ok=True)

    panels: dict[str, dict[str, Any]] = {}
    provenance_rows: list[dict[str, Any]] = []

    for mechanism in MODBAM_GROUP_ORDER:
        representative = representatives[mechanism]

        png = run_modbamtools_plot(
            modbamtools_bin=modbamtools_bin,
            bam_path=representative["modbam"],
            sample_label=representative["display_label"],
            mechanism=mechanism,
            region=region,
            output_dir=panel_dir,
            gtf_path=gtf_path,
        )

        use_hap = mechanism in {"Control", "Disease control", "PWS-mUPD"}
        panels[mechanism] = {
            **representative,
            "image": png,
            "haplotype_grouping": use_hap,
        }

        provenance_rows.append(
            {
                "molecular_mechanism": mechanism,
                "sample_id": representative["sample_id"],
                "display_label": representative["display_label"],
                "selection_rule": "IC depth closest to within-group median",
                "sample_total_ic_depth": fmt(representative["total_ic_depth"]),
                "group_median_ic_depth": fmt(representative["group_median_depth"]),
                "region": region,
                "modbam": str(representative["modbam"]),
                "plot_filename": png.name,
                "renderer": "ModBAMtools" if USE_EXTERNAL_MODBAMTOOLS else "internal MM/ML renderer",
                "haplotype_grouping": "HP tag" if use_hap else "not forced (hemizygous IC)",
            }
        )

    return panels, provenance_rows


def build_all_sample_modbamtools_panels(
    sample_files: dict[str, dict[str, Path | None]],
    outdir: Path,
) -> list[dict[str, Any]]:
    """Generate raw ModBAM evidence for every participant as Extended Data."""
    rows = []
    panel_dir = outdir / "extended_data" / "modbamtools_all_samples"
    labels = sample_display_labels()
    for sample_id, _clinical, mechanism in sorted_cohort():
        bam = sample_files[sample_id].get("modbam")
        if bam is None or not Path(bam).exists():
            rows.append({"sample_id": sample_id, "molecular_mechanism": mechanism, "status": "missing_modbam", "plot_filename": ""})
            continue
        try:
            png = run_modbamtools_plot(
                MODBAMTOOLS_BIN, Path(bam), labels[sample_id], mechanism,
                MODBAM_PLOT_REGION, panel_dir, MODBAM_GTF,
            )
            rows.append({"sample_id": sample_id, "display_label": labels[sample_id], "molecular_mechanism": mechanism,
                         "status": "generated", "plot_filename": png.name, "region": MODBAM_PLOT_REGION,
                         "MAPQ_min": MODBAM_MIN_MAPQ})
        except Exception as exc:
            rows.append({"sample_id": sample_id, "display_label": labels[sample_id], "molecular_mechanism": mechanism,
                         "status": "failed", "plot_filename": "", "error": str(exc)})
    return rows


def render_extended_modbam_contact_sheet(outdir: Path, rows: list[dict[str, Any]]) -> None:
    generated = [r for r in rows if r.get("status") == "generated"]
    if not generated:
        return
    panel_dir = outdir / "extended_data" / "modbamtools_all_samples"
    ncol = 3
    nrow = math.ceil(len(generated) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(16, 4.5 * nrow), squeeze=False)
    for ax in axes.flat: ax.axis("off")
    for ax, row in zip(axes.flat, generated):
        image = panel_dir / row["plot_filename"]
        ax.imshow(plt.imread(image)); ax.axis("off")
        ax.set_title(f"{row['display_label']} | {row['molecular_mechanism']}", fontsize=9,
                     color=MECHANISM_COLORS[row['molecular_mechanism']], fontweight="bold")
    fig.suptitle("Extended Data — raw ModBAM methylation profiles for all participants", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    outbase = outdir / "extended_data" / "ExtendedData_Figure1_all_modbamtools"
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".png"), dpi=PUBLICATION_DPI, bbox_inches="tight")
    plt.close(fig)


def load_modbamtools_panels_from_provenance(
    provenance_path: Path,
    outdir: Path,
) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    """Load cached Figure 1A PNGs without rerunning ModBAMtools."""
    if not provenance_path.exists():
        raise FileNotFoundError(
            f"Cached ModBAMtools provenance table not found: {provenance_path}"
        )

    rows = read_tsv(provenance_path)
    panels: dict[str, dict[str, Any]] = {}

    for row in rows:
        mechanism = row.get("molecular_mechanism", "")
        if mechanism not in MODBAM_GROUP_ORDER:
            continue

        filename = row.get("plot_filename", "")
        image = outdir / "modbamtools" / filename
        if not filename or not image.exists():
            raise FileNotFoundError(
                f"Cached ModBAMtools image for {mechanism} is missing: {image}"
            )

        panels[mechanism] = {
            "sample_id": row.get("sample_id", ""),
            "display_label": row.get("display_label", ""),
            "molecular_mechanism": mechanism,
            "total_ic_depth": safe_float(row.get("sample_total_ic_depth")),
            "group_median_depth": safe_float(row.get("group_median_ic_depth")),
            "modbam": Path(row["modbam"]) if row.get("modbam") else None,
            "image": image,
            "haplotype_grouping": row.get("haplotype_grouping", "").startswith("HP"),
        }

    missing = [mechanism for mechanism in MODBAM_GROUP_ORDER if mechanism not in panels]
    if missing:
        raise RuntimeError(
            "Cached Figure 1A provenance is incomplete; missing mechanisms: "
            + ", ".join(missing)
        )

    return panels, rows


# ---------------------------------------------------------------------------
# Support metrics
# ---------------------------------------------------------------------------

def build_support_rows(
    summary_rows: list[dict[str, Any]],
    matrix_rows: list[dict[str, Any]],
    sample_files: dict[str, dict[str, Path | None]],
    structural_by_sample: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    """Build depth-aware IC support metrics.

    BAM depth and haplotype balance are calculated directly from the same IC
    interval. This avoids treating a 10x sample and a 35x sample as though they
    had equal opportunity to generate two well-supported haplotype BED tracks.
    """
    summary_by_sample = {r["sample_id"]: r for r in summary_rows}
    matrix_by_sample: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)
    for row in matrix_rows:
        matrix_by_sample[row["sample_id"]][row["haplotype_or_source"]] = row

    labels = sample_display_labels()
    rows: list[dict[str, Any]] = []
    ic_region = f"{CHROM}:{PWS_IC_START}-{PWS_IC_END}"
    ic_len = PWS_IC_END - PWS_IC_START + 1

    for sample_id, _clinical, mechanism in sorted_cohort():
        sm = matrix_by_sample[sample_id]
        combined = sm.get("combined_fallback")
        hap1 = sm.get("hap1")
        hap2 = sm.get("hap2")
        bam = sample_files[sample_id].get("bam")

        total_bam_depth = (
            samtools_coverage_mean_depth(Path(bam), ic_region)
            if bam is not None and Path(bam).exists()
            else None
        )

        hp1_depth = hp2_depth = None
        if bam is not None and Path(bam).exists():
            try:
                hp_depths = haplotype_depths_from_bam(
                    Path(bam), ic_region, ic_len
                )
                hp1_depth = safe_float(hp_depths.get("hap1"))
                hp2_depth = safe_float(hp_depths.get("hap2"))
            except Exception:
                hp1_depth = hp2_depth = None

        if mechanism in {"PWS-DEL", "AS-DEL"}:
            # HP balance is not a biologically meaningful QC statistic in a
            # hemizygous interval. Use the combined/retained methylation track.
            hp_balance = None
            cpg1 = safe_float(combined.get("n_CpGs")) if combined else None
            cpg2 = None
            support_mode = "retained_hemizygous_allele"
            retained_bed_depth = safe_float(combined.get("mean_coverage")) if combined else None
            supporting_depth = retained_bed_depth
            supporting_cpgs = cpg1
            insufficient_support = (
                combined is None
                or combined.get("support_category")
                not in {"higher_coverage", "state_estimable"}
            )
            below_higher_coverage_reference = (
                combined is not None
                and combined.get("support_category") == "state_estimable"
            )
        else:
            if hp1_depth is not None and hp2_depth is not None and (hp1_depth + hp2_depth) > 0:
                hp_balance = min(hp1_depth, hp2_depth) / (hp1_depth + hp2_depth)
            else:
                hp_balance = None

            cpg1 = safe_float(hap1.get("n_CpGs")) if hap1 else None
            cpg2 = safe_float(hap2.get("n_CpGs")) if hap2 else None
            supporting_depth = min(
                [x for x in (hp1_depth, hp2_depth) if x is not None],
                default=None,
            )
            supporting_cpgs = min(
                [x for x in (cpg1, cpg2) if x is not None],
                default=None,
            )
            support_mode = "diploid_haplotype_resolved"
            insufficient_support = (
                hap1 is None
                or hap2 is None
                or hap1.get("support_category")
                not in {"higher_coverage", "state_estimable"}
                or hap2.get("support_category")
                not in {"higher_coverage", "state_estimable"}
            )
            below_higher_coverage_reference = (
                not insufficient_support
                and (
                    hap1.get("support_category") == "state_estimable"
                    or hap2.get("support_category") == "state_estimable"
                )
            )

        if insufficient_support:
            support_tier = "insufficient"
        elif below_higher_coverage_reference:
            support_tier = "state_estimable"
        else:
            support_tier = "higher_coverage"

        structural_status = structural_by_sample[sample_id]["ic_deletion_status"]
        if mechanism in {"PWS-DEL", "AS-DEL"} and structural_status == "confirmed":
            ic_phased = None
            interpretation = "NA: CN/SV-confirmed hemizygous IC interval"
        else:
            ic_phased = block_fraction_for_interval(
                sample_files[sample_id]["blocks"],
                PWS_IC_START,
                PWS_IC_END,
            )
            interpretation = "diploid interval; phase-block coverage shown"

        rows.append({
            "sample_id": sample_id,
            "display_label": labels[sample_id],
            "molecular_mechanism": mechanism,
            "total_ic_depth": fmt(total_bam_depth),
            "bam_total_ic_depth": fmt(total_bam_depth),
            "hp1_ic_depth": fmt(hp1_depth),
            "hp2_ic_depth": fmt(hp2_depth),
            "hp_balance": fmt(hp_balance),
            "hap1_cpgs": fmt(cpg1, 0),
            "hap2_cpgs": fmt(cpg2, 0),
            "combined_cpgs": fmt(safe_float(combined.get("n_CpGs")) if combined else None, 0),
            # Backward-compatible fields:
            "supporting_allele_depth": fmt(supporting_depth),
            "supporting_allele_cpgs": fmt(supporting_cpgs, 0),
            "ic_phased_span_percent": fmt(ic_phased),
            "ic_phasing_interpretation": interpretation,
            "structural_ic_deletion_status": structural_status,
            "deletion_type": structural_by_sample[sample_id].get("deletion_type", ""),
            "cn_event_mean_cn": structural_by_sample[sample_id].get("cn_event_mean_cn", ""),
            "domain_phased_span_percent": summary_by_sample.get(sample_id, {}).get(
                "percent_imprinted_domain_in_phased_block", ""
            ),
            "support_mode": support_mode,
            "support_tier": support_tier,
            "state_support_adequate": str(not insufficient_support),
            "below_higher_coverage_reference": str(below_higher_coverage_reference),
            # Backward-compatible field: True now means genuinely insufficient
            # for state assignment, not merely below the 10x descriptive tier.
            "low_support": str(insufficient_support),
            "hap1_source": hap1.get("data_source", "") if hap1 else "",
            "hap2_source": hap2.get("data_source", "") if hap2 else "",
            "combined_source": combined.get("data_source", "") if combined else "",
        })

    return rows


def render_parental_reference_calibration(
    panel_rows: list[dict[str, Any]],
    reference: ParentalReferenceModel,
    outdir: Path,
) -> None:
    """Supplementary visualization of the empirical parental-state scale.

    Every observed allele/haplotype is shown on the continuous beta axis. The
    figure makes the control-derived paternal centroid, maternal centroid and
    equal-distance boundary explicit, so reviewers can see that primary state
    calls are not driven by an arbitrary 0.85/0.15 cutoff.
    """
    rows = panel_rows
    y = np.arange(len(rows), dtype=float)
    fig, ax = plt.subplots(figsize=(10.5, max(5.0, 0.42 * len(rows) + 1.8)))

    # Continuous decision regions.
    ax.axvspan(0.0, reference.decision_boundary, color=STATE_COLORS["P"], alpha=0.05, lw=0)
    ax.axvspan(reference.decision_boundary, 1.0, color=STATE_COLORS["M"], alpha=0.05, lw=0)
    ax.axvline(reference.paternal_reference, color=STATE_COLORS["P"], lw=1.8, ls="-", label=f"Paternal ref β={reference.paternal_reference:.3f}")
    ax.axvline(reference.maternal_reference, color=STATE_COLORS["M"], lw=1.8, ls="-", label=f"Maternal ref β={reference.maternal_reference:.3f}")
    ax.axvline(reference.decision_boundary, color="#333333", lw=1.2, ls="--", label=f"Equal-distance boundary β={reference.decision_boundary:.3f}")

    for i, row in enumerate(rows):
        mechanism = row["molecular_mechanism"]
        color = MECHANISM_COLORS[mechanism]
        for idx, offset, marker_symbol in ((1, -0.10, "o"), (2, 0.10, "s")):
            if row.get(f"allele_{idx}_status") != "observed":
                continue
            beta = safe_float(row.get(f"allele_{idx}_mean_methylation"))
            lo = safe_float(row.get(f"allele_{idx}_ci_low"))
            hi = safe_float(row.get(f"allele_{idx}_ci_high"))
            if beta is None:
                continue
            if lo is not None and hi is not None:
                ax.hlines(i + offset, lo, hi, color=color, lw=1.0, alpha=0.8, zorder=2)
            ax.plot(beta, i + offset, marker=marker_symbol, ms=5.5, linestyle="none",
                    markerfacecolor=color, markeredgecolor="white", markeredgewidth=0.6, zorder=3)

    ax.set_yticks(y)
    ax.set_yticklabels([r["display_label"] for r in rows], fontsize=8)
    for tick, row in zip(ax.get_yticklabels(), rows):
        tick.set_color(MECHANISM_COLORS[row["molecular_mechanism"]])
        tick.set_fontweight("bold")
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Mean IC methylation (β)", fontsize=9)
    ax.set_title(
        "Supplementary Figure — control-calibrated parental methylation reference scale",
        fontsize=11, loc="left", fontweight="bold"
    )
    ax.grid(axis="x", color="#ECECEC", lw=0.7)
    ax.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=3)
    ax.text(
        0.5, -0.16,
        "Circles/squares denote chromosome/haplotype 1/2; horizontal segments are descriptive CpG-bootstrap 95% intervals. "
        "Fixed 0.85/0.15 thresholds are not used for primary state assignment.",
        transform=ax.transAxes, ha="center", va="top", fontsize=7, color="#666666"
    )
    fig.tight_layout(rect=[0, 0.08, 1, 1])
    base = outdir / "supplementary" / "Supplementary_parental_reference_calibration"
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=PUBLICATION_DPI, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure drawing
# ---------------------------------------------------------------------------


def draw_local_ic_gene_track(ax: plt.Axes) -> None:
    """Shared gene/ICR context for the five ModBAM profiles in Panel A."""
    start_mb = MODBAM_REGION_START / 1e6
    end_mb = MODBAM_REGION_END / 1e6
    ax.set_xlim(start_mb, end_mb)
    ax.set_ylim(0.0, 2.15)
    ax.set_yticks([])
    ax.tick_params(axis="x", bottom=False, labelbottom=False)

    for y, name, color in (
        (1.42, "SNRPN / SNURF", "#444444"),
        (0.66, "SNHG14", "#777777"),
    ):
        ax.annotate(
            "", xy=(end_mb, y), xytext=(start_mb, y),
            arrowprops={"arrowstyle": "-|>", "lw": 1.15, "color": color, "mutation_scale": 7},
        )
        ax.text(
            start_mb + 0.00005, y + 0.12, name,
            ha="left", va="bottom", fontsize=fs(5.8),
            fontstyle="italic", color=color,
        )

    icr_start = PWS_IC_START / 1e6
    icr_end = PWS_IC_END / 1e6
    ax.axvspan(icr_start, icr_end, color="#CC79A7", alpha=0.18, lw=0)
    ax.text(
        (icr_start + icr_end) / 2.0, 2.02,
        f"PWS/AS ICR  {PWS_IC_START:,}–{PWS_IC_END:,}",
        ha="center", va="top", fontsize=fs(5.8),
        color="#A64C91", fontweight="bold",
    )
    ax.set_title(
        f"Local gene context · {CHROM}:{MODBAM_REGION_START:,}–{MODBAM_REGION_END:,}",
        loc="left", fontsize=fs(6.2), fontweight="bold", pad=1,
    )
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_single_molecule_panel(
    fig: plt.Figure,
    spec: Any,
    panels: dict[str, dict[str, Any]],
    region: str,
    panel_rows: list[dict[str, Any]],
    connector_ax: plt.Axes | None = None,
    panel_label: str = "A",
) -> list[plt.Axes]:
    """Stack readable representative ModBAM tracks and connect each block.

    The enclosing Figure 1 layout owns the shared Panel A heading and caption.
    Track heights follow their source-image aspect ratios, preserving both the
    methylation-frequency panel and raw-molecule rows without vertical crushing.
    """
    available = [mechanism for mechanism in MODBAM_GROUP_ORDER if mechanism in panels]
    loaded_images = {mechanism: plt.imread(panels[mechanism]["image"]) for mechanism in available}
    image_height_ratios = [
        max(0.10, image.shape[0] / max(image.shape[1], 1))
        for image in loaded_images.values()
    ]
    grid = spec.subgridspec(
        len(available) + 1, 1,
        height_ratios=[0.085, *image_height_ratios],
        hspace=0.18,
    )

    gene_ax = fig.add_subplot(grid[0, 0])
    draw_local_ic_gene_track(gene_ax)

    axes: list[plt.Axes] = []
    titles = {
        "Control": "Control · biparental",
        "PWS-DEL": "PWS-DEL · maternal retained",
        "AS-DEL": "AS-DEL · paternal retained",
        "PWS-mUPD": "PWS-mUPD · maternal + maternal",
        "Disease control": "Disease control · biparental",
    }

    for track_index, mechanism in enumerate(available):
        group_indices = [
            i for i, row in enumerate(panel_rows)
            if row["molecular_mechanism"] == mechanism
        ]
        if not group_indices:
            continue

        group_center = float(np.mean(group_indices))
        ax = fig.add_subplot(grid[track_index + 1, 0])
        axes.append(ax)

        panel = panels[mechanism]
        image = loaded_images[mechanism]
        ax.imshow(image, aspect="auto", interpolation="none")
        ax.axis("off")

        grouping = "HP-grouped" if panel.get("haplotype_grouping") else "hemizygous"
        representative_row = next(
            (row for row in panel_rows if row["sample_id"] == panel.get("sample_id")),
            None,
        )
        observed_percentages: list[str] = []
        if representative_row is not None:
            for allele_index in (1, 2):
                if representative_row.get(f"allele_{allele_index}_status") != "observed":
                    continue
                beta = safe_float(representative_row.get(f"allele_{allele_index}_mean_methylation"))
                if beta is not None:
                    observed_percentages.append(f"{100 * beta:.0f}%")
        methylation_label = "/".join(observed_percentages) if observed_percentages else "n/a"
        ax.set_title(
            f"{titles[mechanism]} | {panel['display_label']} | {grouping}",
            loc="left", pad=1.5,
            fontsize=fs(5.9), weight="bold",
            color=MECHANISM_COLORS[mechanism],
        )
        ax.text(
            0.995, 1.015, f"IC methylation: {methylation_label}",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=fs(7.0), fontweight="bold",
            color=MECHANISM_COLORS[mechanism],
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.92, "pad": 0.8},
        )

        if connector_ax is not None:
            connector = ConnectionPatch(
                xyA=(1.5, group_center), coordsA=connector_ax.transData,
                xyB=(0.0, 0.5), coordsB=ax.transAxes,
                arrowstyle="-|>", mutation_scale=8,
                color=MECHANISM_COLORS[mechanism],
                lw=0.9, alpha=0.75,
                connectionstyle="arc3,rad=0.0",
                clip_on=False, zorder=8,
            )
            fig.add_artist(connector)

    return axes


def draw_cohort_methylation_panel(
    note_ax: plt.Axes, heat_ax: plt.Axes, panel_a_rows: list[dict[str, Any]],
    inference_rows: list[dict[str, Any]] | None = None,
    parental_reference: ParentalReferenceModel | None = None,
    panel_label: str = "A",
    cbar_ax: plt.Axes | None = None,
    cbar_orientation: str = "vertical",
    caption_ax: plt.Axes | None = None,
) -> None:
    n = len(panel_a_rows); values = np.full((n, 2), np.nan); statuses = np.empty((n, 2), dtype=object)
    insufficient_support = np.zeros((n, 2), dtype=bool)
    for i, row in enumerate(panel_a_rows):
        for j in range(2):
            prefix = f"allele_{j+1}"; status = row.get(f"{prefix}_status", "missing"); statuses[i,j] = status
            if status == "observed":
                values[i,j] = safe_float(row.get(f"{prefix}_mean_methylation")) or np.nan
                insufficient_support[i,j] = (
                    row.get(f"{prefix}_coverage_status") == "insufficient"
                )
    cmap = plt.get_cmap("coolwarm").copy(); cmap.set_bad("#FFFFFF")
    image = heat_ax.imshow(values, aspect="auto", cmap=cmap, norm=TwoSlopeNorm(vmin=0, vcenter=.5, vmax=1))
    heat_ax.set_xticks([0,1]); heat_ax.set_xticklabels(["Allele / haplotype 1", "Allele / haplotype 2"], fontsize=fs(7.2))
    heat_ax.set_yticks(range(n)); heat_ax.set_yticklabels([r["display_label"] for r in panel_a_rows], fontsize=fs(8), fontweight="bold")
    heat_ax.tick_params(length=0)
    for i,row in enumerate(panel_a_rows):
        mech=row["molecular_mechanism"]; heat_ax.get_yticklabels()[i].set_color(MECHANISM_COLORS[mech])
        for j in range(2):
            p=f"allele_{j+1}"; status=statuses[i,j]
            if status != "observed":
                style = {
                    "deleted": (STATE_COLORS["deleted"], "///", "deleted"),
                    "missing": (STATE_COLORS["missing"], "...", "no hap-CpG"),
                    "structural_unresolved": (STATE_COLORS["structural_unresolved"], "xx", "CN/SV?")
                }.get(status, ("#F5F5F5", "...", status))
                heat_ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor=style[0],edgecolor="#999",hatch=style[1],linewidth=.8,zorder=3))
                heat_ax.text(j,i,style[2],ha="center",va="center",fontsize=fs(6.0),color="#444",zorder=4)
                continue
            beta=values[i,j]; patt=row.get(f"{p}_pattern_short","?")
            lo=safe_float(row.get(f"{p}_ci_low")); hi=safe_float(row.get(f"{p}_ci_high"))
            if insufficient_support[i,j]:
                heat_ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor="none",edgecolor="#B00020",linewidth=1.4,zorder=3))
            tc="white" if beta>=.72 or beta<=.20 else "#111"
            ci_text = f"\n[{lo:.2f},{hi:.2f}]" if lo is not None and hi is not None else ""
            heat_ax.text(j,i,f"{beta:.2f} {patt}{ci_text}",ha="center",va="center",fontsize=fs(5.8),color=tc,zorder=4)
    for i in range(n-1):
        if panel_a_rows[i]["molecular_mechanism"] != panel_a_rows[i+1]["molecular_mechanism"]:
            heat_ax.axhline(i+.5, color="#2F2F2F", lw=1.35, zorder=6)
    note_ax.set_xlim(0,1); note_ax.set_ylim(heat_ax.get_ylim()); note_ax.axis("off")
    for mech in sorted(MECHANISM_ORDER,key=MECHANISM_ORDER.get):
        idx=[i for i,r in enumerate(panel_a_rows) if r["molecular_mechanism"]==mech]
        if idx:
            note_ax.text(.02,.5*(idx[0]+idx[-1]),f"{mech}\n(n={len(idx)})",ha="left",va="center",fontsize=fs(8),fontweight="bold",color=MECHANISM_COLORS[mech])
            heat_ax.add_patch(Rectangle(
                (-0.5, idx[0] - 0.5), 2.0, idx[-1] - idx[0] + 1.0,
                facecolor="none", edgecolor=MECHANISM_COLORS[mech],
                linewidth=2.6, zorder=7, clip_on=False,
            ))
    heat_ax.set_title(
        "Cohort methylation heatmap",
        fontsize=fs(8.2), loc="left", pad=7, weight="bold"
    )
    cbar = heat_ax.figure.colorbar(
        image,
        ax=None if cbar_ax is not None else heat_ax,
        cax=cbar_ax,
        fraction=.048,
        pad=.025,
        orientation=cbar_orientation,
    )
    cbar.set_label("Mean IC methylation (β)", fontsize=fs(7.0))
    cbar.ax.tick_params(labelsize=fs(6.5), length=2)
    if parental_reference is not None:
        state_note = (
            f"Control-calibrated references: P={parental_reference.paternal_reference:.3f}, "
            f"M={parental_reference.maternal_reference:.3f}, boundary={parental_reference.decision_boundary:.3f}.\n"
            "Brackets: descriptive CpG-bootstrap 95% intervals."
        )
    else:
        state_note = "Control-calibrated parental reference; brackets show descriptive CpG-bootstrap 95% intervals"
    note_lines = [state_note]
    if inference_rows:
        r=inference_rows[0]
        note_lines.append(
            f"Participant-level PWS-DEL vs AS-DEL: Δmedian={r['delta_median']} "
            f"[95% CI {r['bootstrap_95CI_low']}, {r['bootstrap_95CI_high']}], "
            f"exact permutation P={r['exact_permutation_p_two_sided']}"
        )
    if caption_ax is not None:
        caption_ax.axis("off")
        caption_ax.text(
            0.0, 0.72, "\n".join(note_lines),
            ha="left", va="top", fontsize=fs(5.5), color="#555555",
            linespacing=1.35,
        )
    else:
        heat_ax.text(
            .5, -.10, "\n".join(note_lines),
            transform=heat_ax.transAxes, ha="center", va="top",
            fontsize=fs(5.6), color="#555555",
        )


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


def _spread_label_positions(
    values: list[float],
    minimum_gap: float = 0.045,
    lower: float = 0.035,
    upper: float = 0.90,
) -> list[float]:
    """Return vertically separated label positions while preserving order."""
    if not values:
        return []
    order = np.argsort(values)
    placed = np.asarray([np.clip(values[i], lower, upper) for i in order], dtype=float)
    for i in range(1, len(placed)):
        placed[i] = max(placed[i], placed[i - 1] + minimum_gap)
    if placed[-1] > upper:
        placed -= placed[-1] - upper
    for i in range(len(placed) - 2, -1, -1):
        placed[i] = min(placed[i], placed[i + 1] - minimum_gap)
    if placed[0] < lower:
        placed += lower - placed[0]
    result = np.empty(len(placed), dtype=float)
    result[order] = placed
    return result.tolist()


def _convex_hull(points: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """Andrew monotone-chain hull used for the closed methylation clusters."""
    unique = sorted(set(points))
    if len(unique) <= 1:
        return unique

    def cross(
        origin: tuple[float, float],
        a: tuple[float, float],
        b: tuple[float, float],
    ) -> float:
        return (a[0] - origin[0]) * (b[1] - origin[1]) - (a[1] - origin[1]) * (b[0] - origin[0])

    lower: list[tuple[float, float]] = []
    for point in unique:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[float, float]] = []
    for point in reversed(unique):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def draw_panel_c(
    ax: plt.Axes,
    deletion_profile_rows: list[dict[str, Any]],
) -> None:
    """Compare each retained deletion track directly with control chromosomes."""
    ax.set_title(
        "C. Retained methylation tracks versus control chromosomes",
        fontsize=fs(10.2), loc="left", pad=12, weight="bold",
    )
    deletion_rows = [
        row for row in deletion_profile_rows
        if row.get("track_role") == "retained deletion chromosome"
    ]
    control_rows = [
        row for row in deletion_profile_rows
        if row.get("track_role") == "reference"
        and safe_float(row.get("parental_profile_score")) is not None
    ]
    maternal_control_scores = [
        float(row["parental_profile_score"])
        for row in control_rows if row.get("IC_parental_class") == "maternal-like"
    ]
    paternal_control_scores = [
        float(row["parental_profile_score"])
        for row in control_rows if row.get("IC_parental_class") == "paternal-like"
    ]
    maternal_centroid = float(np.median(maternal_control_scores)) if maternal_control_scores else 1.0
    paternal_centroid = float(np.median(paternal_control_scores)) if paternal_control_scores else -1.0

    ax.axvspan(-1.05, -DELETION_PROFILE_MIN_SCORE_MAGNITUDE,
               color=STATE_COLORS["P"], alpha=0.055, lw=0, zorder=0)
    ax.axvspan(DELETION_PROFILE_MIN_SCORE_MAGNITUDE, 1.05,
               color=STATE_COLORS["M"], alpha=0.055, lw=0, zorder=0)
    ax.axvspan(-DELETION_PROFILE_MIN_SCORE_MAGNITUDE,
               DELETION_PROFILE_MIN_SCORE_MAGNITUDE,
               color="#9E9E9E", alpha=0.12, lw=0, zorder=0)
    ax.axvline(0, color="#666666", lw=0.75, ls="--", zorder=1)
    ax.axvline(maternal_centroid, color=STATE_COLORS["M"], lw=1.25, ls="--", zorder=1)
    ax.axvline(paternal_centroid, color=STATE_COLORS["P"], lw=1.25, ls="--", zorder=1)

    ordered_rows = sorted(
        deletion_rows,
        key=lambda row: (
            MECHANISM_ORDER.get(str(row.get("molecular_mechanism")), 99),
            str(row.get("sample_id", "")),
        ),
    )
    y_positions = np.arange(1, len(ordered_rows) + 1, dtype=float)

    # The four unaffected-control chromosomes occupy one reference row. Labels
    # are consolidated by parental class so they cannot collide with markers.
    state_counts: dict[str, int] = defaultdict(int)
    state_labels: dict[str, list[str]] = defaultdict(list)
    for row in control_rows:
        score = float(row["parental_profile_score"])
        state = "M" if row.get("IC_parental_class") == "maternal-like" else "P"
        y = -0.09 + 0.18 * state_counts[state]
        state_counts[state] += 1
        state_labels[state].append(f"{row['display_label']} {row['haplotype_label']}")
        ax.scatter(
            score, y, s=42, marker="D", facecolor=STATE_COLORS[state],
            edgecolor="#333333", linewidth=0.8, zorder=4,
        )

    for y, row in zip(y_positions, ordered_rows):
        score = safe_float(row.get("parental_profile_score"))
        mechanism = str(row.get("molecular_mechanism", ""))
        profile_class = str(row.get("profile_parental_class", "uncertain"))
        if score is None:
            score = 0.0
            profile_class = "insufficient"
        nearest_centroid = (
            maternal_centroid
            if abs(score - maternal_centroid) <= abs(score - paternal_centroid)
            else paternal_centroid
        )
        state = (
            "M" if profile_class == "maternal-like"
            else "P" if profile_class == "paternal-like"
            else "?"
        )
        ax.hlines(
            y, min(score, nearest_centroid), max(score, nearest_centroid),
            color=MECHANISM_COLORS[mechanism], lw=0.9,
            ls=":" if profile_class in {"uncertain", "insufficient"} else "-",
            alpha=0.58, zorder=2,
        )
        ax.scatter(
            score, y, s=68, marker=MECHANISM_MARKERS[mechanism],
            facecolor=STATE_COLORS[state], edgecolor=MECHANISM_COLORS[mechanism],
            linewidth=1.1, zorder=4,
        )
        short_class = {
            "maternal-like": "M-like", "paternal-like": "P-like",
            "uncertain": "uncertain", "insufficient": "insufficient",
        }.get(profile_class, profile_class)
        label_x = 1.085
        ax.annotate(
            f"{short_class}  (score {score:+.2f})",
            xy=(score, y), xytext=(label_x, y), textcoords="data",
            ha="left", va="center", fontsize=fs(5.25),
            color=MECHANISM_COLORS[mechanism], fontweight="bold", zorder=5,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.84, "pad": 0.7},
            arrowprops={
                "arrowstyle": "-", "lw": 0.42,
                "color": MECHANISM_COLORS[mechanism], "alpha": 0.55,
            },
        )

    ax.set_xlim(-1.06, 1.46)
    ax.set_ylim(len(ordered_rows) + 0.65, -0.60)
    ax.set_yticks([0.0, *y_positions])
    ax.set_yticklabels(
        ["Control chromosomes", *[row["display_label"] for row in ordered_rows]],
        fontsize=fs(6.6), fontweight="bold",
    )
    for tick, row in zip(ax.get_yticklabels()[1:], ordered_rows):
        tick.set_color(MECHANISM_COLORS[str(row["molecular_mechanism"])])
    ax.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    ax.set_xlabel(
        "Similarity to control chromosome methylation profiles\n"
        "paternal-like ←   profile score   → maternal-like",
        fontsize=fs(7.6),
    )
    ax.grid(axis="y", color="#E7E7E7", lw=0.65)
    ax.tick_params(axis="x", labelsize=fs(6.7))
    ax.tick_params(axis="y", length=0)
    ax.text(
        paternal_centroid, -0.43,
        "P control\n" + ", ".join(state_labels.get("P", [])),
        ha="center", va="center",
        color=STATE_COLORS["P"], fontsize=fs(5.5), fontweight="bold",
    )
    ax.text(
        maternal_centroid, -0.43,
        "M control\n" + ", ".join(state_labels.get("M", [])),
        ha="center", va="center",
        color=STATE_COLORS["M"], fontsize=fs(5.5), fontweight="bold",
    )
    ax.text(
        0.01, 0.01,
        "Control endpoints are training anchors; movement toward 0 indicates "
        "weaker whole-deletion profile similarity.",
        transform=ax.transAxes, ha="left", va="bottom",
        fontsize=fs(4.9), color="#666666",
    )

def _support_y_axis(ax: plt.Axes, support_rows: list[dict[str, Any]], show_y: bool) -> None:
    y = np.arange(len(support_rows))
    ax.set_ylim(len(support_rows) - 0.5, -0.5)
    ax.set_yticks(y)
    if show_y:
        ax.set_yticklabels(
            [row["display_label"] for row in support_rows],
            fontsize=fs(7),
            fontweight="bold",
        )
        for tick, row in zip(ax.get_yticklabels(), support_rows):
            tick.set_color(MECHANISM_COLORS[row["molecular_mechanism"]])
    else:
        ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)
    for i in range(len(support_rows) - 1):
        if support_rows[i]["molecular_mechanism"] != support_rows[i + 1]["molecular_mechanism"]:
            ax.axhline(i + 0.5, color="#5C5C5C", lw=0.7)


def draw_total_depth_axis(
    ax: plt.Axes,
    support_rows: list[dict[str, Any]],
    panel_label: str = "c1",
) -> None:
    values = [safe_float(r.get("bam_total_ic_depth")) or 0 for r in support_rows]
    xmax = max(values, default=1) * 1.22
    ax.set_xlim(0, max(xmax, 1))
    ax.set_title(
        f"{panel_label}. Total IC depth",
        fontsize=fs(8), pad=6, fontweight="bold", loc="left"
    )
    ax.grid(axis="x", color="#ECECEC", lw=0.7)
    ax.tick_params(axis="x", labelsize=fs(6.8))
    _support_y_axis(ax, support_rows, True)

    for i, row in enumerate(support_rows):
        value = safe_float(row.get("bam_total_ic_depth"))
        if value is None:
            continue
        color = MECHANISM_COLORS[row["molecular_mechanism"]]
        ax.hlines(i, 0, value, color=color, lw=1.2, alpha=0.25)
        ax.plot(value, i, "o", ms=5.0, color=color)
        ax.annotate(f"{value:.1f}x", (value, i), xytext=(4, 0),
                    textcoords="offset points", va="center", fontsize=fs(5.6))


def draw_hp_depth_axis(
    ax: plt.Axes,
    support_rows: list[dict[str, Any]],
    panel_label: str = "c2",
) -> None:
    vals = [
        safe_float(r.get(k)) or 0
        for r in support_rows
        for k in ("hp1_ic_depth", "hp2_ic_depth")
    ]
    xmax = max(vals, default=1) * 1.25
    ax.set_xlim(0, max(xmax, 1))
    ax.set_title(
        f"{panel_label}. HP1 / HP2 IC depth",
        fontsize=fs(8), pad=6, fontweight="bold", loc="left"
    )
    ax.grid(axis="x", color="#ECECEC", lw=0.7)
    ax.tick_params(axis="x", labelsize=fs(6.8))
    _support_y_axis(ax, support_rows, False)

    for i, row in enumerate(support_rows):
        mechanism = row["molecular_mechanism"]
        if mechanism in {"PWS-DEL", "AS-DEL"}:
            ax.text(0.02, i, "hemizygous", transform=ax.get_yaxis_transform(),
                    ha="left", va="center", fontsize=fs(5.5), color="#777777")
            continue

        h1 = safe_float(row.get("hp1_ic_depth"))
        h2 = safe_float(row.get("hp2_ic_depth"))
        if h1 is not None:
            ax.plot(h1, i - 0.10, marker="o", ms=4.5, linestyle="none",
                    markerfacecolor="white", markeredgecolor="#0072B2", markeredgewidth=1.1)
        if h2 is not None:
            ax.plot(h2, i + 0.10, marker="s", ms=4.2, linestyle="none",
                    markerfacecolor="white", markeredgecolor="#D55E00", markeredgewidth=1.1)
        if h1 is not None and h2 is not None:
            ax.hlines(i, min(h1, h2), max(h1, h2), color="#AAAAAA", lw=0.8, zorder=0)

    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="white",
                   markeredgecolor="#0072B2", label="HP1"),
            Line2D([0], [0], marker="s", linestyle="none", markerfacecolor="white",
                   markeredgecolor="#D55E00", label="HP2"),
        ],
        frameon=False, fontsize=fs(5.8), loc="lower right"
    )


def draw_cpg_support_axis(
    ax: plt.Axes,
    support_rows: list[dict[str, Any]],
    panel_label: str = "c3",
) -> None:
    vals = []
    for r in support_rows:
        if r["molecular_mechanism"] in {"PWS-DEL", "AS-DEL"}:
            vals.append(safe_float(r.get("combined_cpgs")) or 0)
        else:
            vals.extend([safe_float(r.get("hap1_cpgs")) or 0, safe_float(r.get("hap2_cpgs")) or 0])
    xmax = max(vals, default=1) * 1.20
    ax.set_xlim(0, max(xmax, 1))
    ax.set_title(
        f"{panel_label}. CpGs by allele / HP",
        fontsize=fs(8), pad=6, fontweight="bold", loc="left"
    )
    ax.grid(axis="x", color="#ECECEC", lw=0.7)
    ax.tick_params(axis="x", labelsize=fs(6.8))
    _support_y_axis(ax, support_rows, False)

    for i, row in enumerate(support_rows):
        mechanism = row["molecular_mechanism"]
        if mechanism in {"PWS-DEL", "AS-DEL"}:
            c = safe_float(row.get("combined_cpgs"))
            if c is not None:
                ax.plot(c, i, marker="D", ms=4.5, linestyle="none",
                        color=MECHANISM_COLORS[mechanism])
                ax.annotate(f"{c:.0f}", (c, i), xytext=(4, 0),
                            textcoords="offset points", va="center", fontsize=fs(5.5))
        else:
            c1 = safe_float(row.get("hap1_cpgs"))
            c2 = safe_float(row.get("hap2_cpgs"))
            if c1 is not None:
                ax.plot(c1, i - 0.10, marker="o", ms=4.2, linestyle="none",
                        markerfacecolor="white", markeredgecolor="#0072B2")
            if c2 is not None:
                ax.plot(c2, i + 0.10, marker="s", ms=4.0, linestyle="none",
                        markerfacecolor="white", markeredgecolor="#D55E00")


def draw_hp_balance_axis(
    ax: plt.Axes,
    support_rows: list[dict[str, Any]],
    panel_label: str = "c4",
) -> None:
    ax.set_xlim(0, 0.52)
    ax.set_xticks([0, 0.25, 0.5])
    ax.set_title(
        f"{panel_label}. HP balance",
        fontsize=fs(8), pad=6, fontweight="bold", loc="left"
    )
    ax.grid(axis="x", color="#ECECEC", lw=0.7)
    ax.tick_params(axis="x", labelsize=fs(6.8))
    ax.axvline(0.5, color="#777777", lw=0.8, ls="--")
    _support_y_axis(ax, support_rows, False)

    for i, row in enumerate(support_rows):
        if row["molecular_mechanism"] in {"PWS-DEL", "AS-DEL"}:
            ax.text(0.02, i, "NA", ha="left", va="center",
                    fontsize=fs(5.8), color="#777777", fontstyle="italic")
            continue
        balance = safe_float(row.get("hp_balance"))
        if balance is None:
            ax.text(0.02, i, "unresolved", ha="left", va="center",
                    fontsize=fs(5.3), color="#777777")
            continue
        color = MECHANISM_COLORS[row["molecular_mechanism"]]
        ax.plot(balance, i, marker="o", ms=5.0, color=color)
        ax.annotate(f"{balance:.2f}", (balance, i), xytext=(4, 0),
                    textcoords="offset points", va="center", fontsize=fs(5.4))

    ax.text(
        0.5, -0.08,
        "0.50 = perfectly balanced",
        transform=ax.transAxes, ha="center", va="top",
        fontsize=fs(5.5), color="#777777"
    )


def draw_sequencing_support_panel(
    metric_axes: list[plt.Axes],
    support_rows: list[dict[str, Any]],
    panel_labels: tuple[str, str, str, str] = ("c1", "c2", "c3", "c4"),
) -> None:
    draw_total_depth_axis(metric_axes[0], support_rows, panel_labels[0])
    draw_hp_depth_axis(metric_axes[1], support_rows, panel_labels[1])
    draw_cpg_support_axis(metric_axes[2], support_rows, panel_labels[2])
    draw_hp_balance_axis(metric_axes[3], support_rows, panel_labels[3])



def _compact_deletion_display(structural_row: dict[str, Any]) -> tuple[str, str]:
    """Return a concise Nature-style deletion label and size label."""
    dtype = str(structural_row.get("deletion_type", "") or "").strip()
    size_raw = structural_row.get("cn_event_size_mb", "")
    try:
        size_val = float(size_raw)
        size_label = f"{size_val:.2f} Mb loss"
    except Exception:
        size_label = "size n/a"

    if not dtype or dtype == "no chr15 deletion":
        return "copy-neutral", size_label

    if "Type I-like" in dtype:
        core = "Type I-like"
    elif "Type II-like" in dtype:
        core = "Type II-like"
    elif "BP1-BP3" in dtype:
        core = "BP1-BP3-like"
    elif "BP2-BP3" in dtype:
        core = "BP2-BP3-like"
    elif "atypical/extended" in dtype:
        left = str(structural_row.get("left_landmark", "") or "").strip()
        right = str(structural_row.get("right_landmark", "") or "").strip()
        if left and right and left != "unassigned" and right != "unassigned":
            core = f"Atypical {left}-{right}"
        else:
            core = "Atypical / extended"
    else:
        core = dtype

    # The mechanism is already encoded by the row label and colour. Repeating
    # it inside every segment wastes space and is the main source of collisions.
    return core, size_label


def read_panel_b_gene_annotations(
    gtf_path: Path | None = PANEL_B_GTF,
) -> list[dict[str, Any]]:
    """Read the curated Panel B gene set from a plain or gzipped GTF."""
    if gtf_path is None or not Path(gtf_path).exists():
        return []

    path = Path(gtf_path)
    opener = gzip.open if path.suffix == ".gz" else open
    wanted = set(PANEL_B_GENE_NAMES)
    genes: dict[str, dict[str, Any]] = {}
    with opener(path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[0] != CHROM or fields[2] != "gene":
                continue
            start = int(fields[3])
            end = int(fields[4])
            if end < CN_PLOT_START or start > CN_PLOT_END:
                continue
            match = re.search(r'(?:^|;\s*)gene "([^"]+)"', fields[8])
            if match is None or match.group(1) not in wanted:
                continue
            name = match.group(1)
            genes[name] = {
                "name": name,
                "start": start,
                "end": end,
                "strand": fields[6],
            }

    return [genes[name] for name in PANEL_B_GENE_NAMES if name in genes]


def draw_panel_b_gene_track(
    ax: plt.Axes,
    gtf_path: Path | None = PANEL_B_GTF,
) -> None:
    """Draw a compact directional gene track and the exact PWS/AS ICR."""
    genes = read_panel_b_gene_annotations(gtf_path)
    ax.set_xlim(CN_PLOT_START / 1e6, CN_PLOT_END / 1e6)
    ax.set_ylim(-0.15, 4.15)
    ax.set_yticks([])
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.set_title(
        "B. Chr15 dosage confirms recurrent deletion classes and sizes",
        loc="left", fontsize=fs(11), fontweight="bold", pad=8,
    )

    # Greedy lane allocation separates neighbouring and overlapping labels.
    lane_ends = [float("-inf")] * 4
    for gene in sorted(genes, key=lambda item: (item["start"], item["end"])):
        start = gene["start"] / 1e6
        end = gene["end"] / 1e6
        midpoint = (start + end) / 2.0
        label_width = max(0.22, 0.075 * len(gene["name"]))
        occupied_start = min(start, midpoint - label_width / 2.0)
        occupied_end = max(end, midpoint + label_width / 2.0)
        eligible = [i for i, lane_end in enumerate(lane_ends) if occupied_start > lane_end]
        lane = eligible[0] if eligible else int(np.argmin(lane_ends))
        lane_ends[lane] = occupied_end
        y = 0.42 + lane * 0.88
        arrow_start, arrow_end = (start, end) if gene["strand"] == "+" else (end, start)
        ax.annotate(
            "", xy=(arrow_end, y), xytext=(arrow_start, y),
            arrowprops={"arrowstyle": "-|>", "lw": 0.75, "color": "#4D4D4D", "mutation_scale": 5},
        )
        ax.text(
            midpoint, y + 0.16, gene["name"],
            ha="center", va="bottom", fontsize=fs(5.2),
            color="#333333", fontstyle="italic",
        )

    icr_mid = (PWS_IC_START + PWS_IC_END) / 2.0 / 1e6
    ax.axvline(icr_mid, color="#CC79A7", lw=1.5, zorder=5)
    ax.annotate(
        f"PWS/AS ICR\n{PWS_IC_START / 1e6:.3f}–{PWS_IC_END / 1e6:.3f} Mb",
        xy=(icr_mid, 0.12), xytext=(icr_mid + 0.72, 3.88),
        ha="left", va="top", fontsize=fs(5.6), fontweight="bold",
        color="#A64C91",
        arrowprops={"arrowstyle": "-", "lw": 0.7, "color": "#CC79A7"},
    )

    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_cohort_cnv_track(
    ax: plt.Axes,
    structural_by_sample: dict[str, dict[str, Any]],
    show_legend: bool = True,
    show_title: bool = True,
) -> None:
    """Main Figure 1B: compact cohort-wide chr15 CNV/deletion track."""
    cohort = sorted_cohort()
    labels = sample_display_labels()
    y = np.arange(len(cohort), dtype=float)

    ax.set_xlim(CN_PLOT_START / 1e6, CN_PLOT_END / 1e6)
    ax.set_ylim(len(cohort) - 0.55, -0.55)
    ax.set_yticks(y)
    ax.set_yticklabels([labels[sid] for sid, _c, _m in cohort], fontsize=fs(7.0))
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("T2T-CHM13v2.0 chr15 coordinate (Mb)", fontsize=fs(8.5))
    if show_title:
        ax.set_title(
            "B. Chr15 dosage confirms recurrent deletion classes and sizes",
            loc="left", fontsize=fs(11), fontweight="bold", pad=8,
        )

    for i, (sample_id, _clinical, mechanism) in enumerate(cohort):
        color = MECHANISM_COLORS[mechanism]
        ax.hlines(
            i, CN_PLOT_START / 1e6, CN_PLOT_END / 1e6,
            color="#D8D8D8", lw=2.6, zorder=1
        )
        ax.get_yticklabels()[i].set_color(color)
        ax.get_yticklabels()[i].set_fontweight("bold")

        s = structural_by_sample.get(sample_id, {})
        event_start = safe_float(s.get("cn_event_start"))
        event_end = safe_float(s.get("cn_event_end"))
        event_cn = safe_float(s.get("cn_event_mean_cn"))

        if event_start is not None and event_end is not None and event_end > event_start:
            left = event_start / 1e6
            width = (event_end - event_start) / 1e6
            right = left + width

            ax.add_patch(
                Rectangle(
                    (left, i - 0.22), width, 0.44,
                    facecolor=color, edgecolor=color,
                    linewidth=0.5, alpha=0.74, zorder=3
                )
            )

            class_label, size_label = _compact_deletion_display(s)
            cn_label = f" · CN {event_cn:.1f}" if event_cn is not None else ""
            if width >= 1.4:
                ax.text(
                    left + width / 2, i, f"{class_label} · {size_label}{cn_label}",
                    ha="center", va="center", fontsize=fs(4.9),
                    color="white", fontweight="bold", zorder=4
                )
            else:
                text_x = min(right + 0.10, CN_PLOT_END / 1e6 - 0.02)
                ax.text(
                    text_x, i - 0.08, class_label,
                    ha="left", va="center", fontsize=fs(5.2),
                    color=color, fontweight="bold", zorder=4,
                )
                ax.text(
                    text_x, i + 0.10, f"{size_label}{cn_label}",
                    ha="left", va="center", fontsize=fs(5.0),
                    color="#444444", zorder=4,
                )

        elif mechanism in {"PWS-DEL", "AS-DEL"}:
            ax.add_patch(
                Rectangle(
                    (CN_PLOT_START / 1e6, i - 0.20),
                    (CN_PLOT_END - CN_PLOT_START) / 1e6,
                    0.40, facecolor="none", edgecolor=color,
                    hatch="xx", linewidth=0.8, alpha=0.65, zorder=2
                )
            )
            ax.text(
                CN_PLOT_END / 1e6 - 0.08, i, "structural unresolved",
                ha="right", va="center", fontsize=fs(5.3), color=color
            )
        else:
            ax.hlines(
                i, CN_PLOT_START / 1e6, CN_PLOT_END / 1e6,
                color=color, lw=1.15, alpha=0.78, zorder=2
            )

    for name, pos in BREAKPOINT_LANDMARKS.items():
        ax.axvline(pos / 1e6, color="#B5B5B5", lw=0.6, ls="--", zorder=0)
        ax.text(
            pos / 1e6, -0.43, name,
            ha="center", va="bottom", rotation=90,
            fontsize=fs(5.1), color="#777777"
        )

    ax.axvspan(
        PWS_IC_START / 1e6, PWS_IC_END / 1e6,
        color="#CC79A7", alpha=0.18, linewidth=0, zorder=0
    )
    ax.axvline(
        (PWS_IC_START + PWS_IC_END) / 2.0 / 1e6,
        color="#CC79A7", lw=1.0, alpha=0.85, zorder=2,
    )

    for i in range(len(cohort) - 1):
        if cohort[i][2] != cohort[i + 1][2]:
            ax.axhline(i + 0.5, color="#555555", lw=0.55)

    ax.grid(axis="x", color="#EFEFEF", lw=0.55)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    handles = [
        Line2D([0], [0], color=MECHANISM_COLORS[m], lw=3, label=m)
        for m in sorted(MECHANISM_ORDER, key=MECHANISM_ORDER.get)
    ]
    if show_legend:
        ax.legend(
            handles=handles,
            frameon=False, ncol=3, fontsize=fs(6.0),
            loc="upper center", bbox_to_anchor=(0.5, -0.16),
            handlelength=2.0, columnspacing=1.1, borderaxespad=0.0,
        )

def _state_color_from_pattern(pattern: str, fallback: str = "#666666") -> str:
    if pattern == "maternal-pattern":
        return STATE_COLORS["M"]
    if pattern == "paternal-pattern":
        return STATE_COLORS["P"]
    if pattern == "uncertain":
        return "#B68B00"
    if pattern == "low-support":
        return "#8C8C8C"
    return fallback


def _plot_bed_cpg_track(
    ax: plt.Axes,
    stats: BedStats,
    label: str,
    color: str,
    x_midpoint: float,
    marker: str = "o",
    linestyle: str = "-",
) -> None:
    values = stats.values_by_pos or {}
    if not values:
        return
    ordered = sorted(values.items())
    x = np.asarray([(p - x_midpoint) / 1000.0 for p, _ in ordered], dtype=float)
    y = np.asarray([v[0] for _p, v in ordered], dtype=float)
    ax.plot(
        x, y, marker=marker, ms=2.6, lw=0.8,
        linestyle=linestyle, color=color, alpha=0.9,
        label=label, zorder=3
    )


def render_supplementary_methylation_tracks(
    sample_files: dict[str, dict[str, Path | None]],
    matrix_rows: list[dict[str, Any]],
    outdir: Path,
) -> None:
    """Supplementary Figure S1: per-sample CpG methylation tracks around the IC.

    Unlike the main heatmap, these tracks expose the underlying CpG-to-CpG
    measurements for every participant. Deletion genomes use combined methylation
    because the interval is hemizygous; diploid/mUPD/disease-control genomes use
    hap1 and hap2 tracks. HP1/HP2 remain numerical phase labels; M-like/P-like
    colors reflect the control-calibrated methylation pattern rather than trio
    parent-of-origin assignment.
    """
    labels = sample_display_labels()
    row_map: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)
    for row in matrix_rows:
        row_map[row["sample_id"]][row["haplotype_or_source"]] = row

    cohort = sorted_cohort()
    ncol = 2
    nrow = math.ceil(len(cohort) / ncol)
    fig, axes = plt.subplots(
        nrow, ncol,
        figsize=(14.5, max(12.0, 2.05 * nrow)),
        sharex=True, sharey=True, squeeze=False
    )
    for ax in axes.flat:
        ax.axis("off")

    midpoint = (PWS_IC_START + PWS_IC_END) / 2.0
    display_start = MODBAM_REGION_START
    display_end = MODBAM_REGION_END
    x_left = (display_start - midpoint) / 1000.0
    x_right = (display_end - midpoint) / 1000.0

    # Collect source rows for transparent supplementary source data.
    source_rows: list[dict[str, Any]] = []

    for ax, (sample_id, _clinical, mechanism) in zip(axes.flat, cohort):
        ax.axis("on")
        ax.set_xlim(x_left, x_right)
        ax.set_ylim(-0.03, 1.03)
        ax.axvspan(
            (PWS_IC_START - midpoint) / 1000.0,
            (PWS_IC_END - midpoint) / 1000.0,
            color="#F3E8F1", alpha=0.45, lw=0, zorder=0
        )
        ax.grid(axis="y", color="#EFEFEF", lw=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        rows = row_map.get(sample_id, {})
        files = sample_files[sample_id]

        if mechanism in {"PWS-DEL", "AS-DEL"}:
            stats = read_bed_region(
                files.get("combined_bed"),
                display_start, display_end, keep_values=True
            )
            row = rows.get("combined_fallback", {})
            pattern = row.get("pattern", "")
            fallback = MECHANISM_COLORS[mechanism]
            color = _state_color_from_pattern(pattern, fallback)
            _plot_bed_cpg_track(
                ax, stats, "combined retained chromosome",
                color, midpoint, marker="o", linestyle="-"
            )
            for pos, (beta, cov) in (stats.values_by_pos or {}).items():
                source_rows.append({
                    "sample_id": sample_id,
                    "display_label": labels[sample_id],
                    "molecular_mechanism": mechanism,
                    "track": "combined",
                    "position_0based": pos,
                    "position_relative_to_IC_midpoint_kb": (pos - midpoint) / 1000.0,
                    "methylation_beta": beta,
                    "coverage": cov,
                    "pattern": pattern,
                })
        else:
            for hap, marker, ls in [
                ("hap1", "o", "-"),
                ("hap2", "s", "--"),
            ]:
                stats = read_bed_region(
                    files.get(f"{hap}_bed"),
                    display_start, display_end, keep_values=True
                )
                row = rows.get(hap, {})
                pattern = row.get("pattern", "")
                color = _state_color_from_pattern(
                    pattern, "#4D4D4D" if hap == "hap1" else "#8A8A8A"
                )
                short = row.get("pattern_short", "?")
                _plot_bed_cpg_track(
                    ax, stats, f"{hap.upper()} ({short})",
                    color, midpoint, marker=marker, linestyle=ls
                )
                for pos, (beta, cov) in (stats.values_by_pos or {}).items():
                    source_rows.append({
                        "sample_id": sample_id,
                        "display_label": labels[sample_id],
                        "molecular_mechanism": mechanism,
                        "track": hap,
                        "position_0based": pos,
                        "position_relative_to_IC_midpoint_kb": (pos - midpoint) / 1000.0,
                        "methylation_beta": beta,
                        "coverage": cov,
                        "pattern": pattern,
                    })

        ax.set_title(
            f"{labels[sample_id]} | {mechanism}",
            loc="left", fontsize=fs(7.2), fontweight="bold",
            color=MECHANISM_COLORS[mechanism], pad=3
        )
        ax.legend(
            frameon=False, fontsize=fs(5.3),
            loc="upper right", handlelength=1.8
        )

    # Shared labels.
    for ax in axes[-1, :]:
        if ax.axison:
            ax.set_xlabel("Position relative to IC midpoint (kb)", fontsize=fs(7.2))
    for row_axes in axes:
        if row_axes[0].axison:
            row_axes[0].set_ylabel("Methylation β", fontsize=fs(7.2))

    fig.suptitle(
        "Supplementary Figure S1 | Per-sample CpG methylation tracks across the PWS/AS imprinting centre",
        fontsize=fs(12), fontweight="bold", y=0.997
    )
    fig.text(
        0.5, 0.006,
        "Deletion genomes: combined retained-chromosome methylation. "
        "Diploid/mUPD genomes: numerical HP1/HP2 tracks; M-like/P-like colors are control-calibrated methylation states.",
        ha="center", va="bottom", fontsize=fs(6.1), color="#666666"
    )
    fig.tight_layout(rect=[0.02, 0.025, 0.99, 0.985])

    base = outdir / "supplementary" / "Supplementary_FigureS1_per_sample_methylation_tracks"
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=PUBLICATION_DPI, bbox_inches="tight")
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

    write_tsv(
        outdir / "source_data" / "Supplementary_FigureS1_per_sample_methylation_tracks_source_data.tsv",
        source_rows,
    )


def render_supplementary_sequencing_qc(
    support_rows: list[dict[str, Any]],
    outdir: Path,
) -> None:
    """Supplementary Figure S2: sequencing depth / haplotype-resolution QC."""
    fig = plt.figure(figsize=(16.0, 8.7))
    grid = GridSpec(
        2, 4, figure=fig,
        height_ratios=[0.10, 1.0],
        width_ratios=[1.0, 1.15, 1.05, 1.0],
        hspace=0.04, wspace=0.30,
    )
    header = fig.add_subplot(grid[0, :])
    header.axis("off")
    header.text(
        0.0, 0.72,
        "Supplementary Figure S2 | Sequencing depth and haplotype-resolution QC at the PWS/AS IC",
        ha="left", va="center", fontsize=fs(11.5), fontweight="bold"
    )
    axes = [fig.add_subplot(grid[1, i]) for i in range(4)]
    draw_sequencing_support_panel(
        axes, support_rows,
        panel_labels=("S2a", "S2b", "S2c", "S2d")
    )
    fig.text(
        0.5, 0.015,
        "Depth is interpreted as measurement support, not biological weight. "
        "Hemizygous PWS/AS deletion intervals are not expected to show balanced HP1/HP2 depth.",
        ha="center", va="bottom", fontsize=fs(6.2), color="#666666"
    )
    fig.subplots_adjust(top=0.96, bottom=0.09, left=0.07, right=0.985)
    base = outdir / "supplementary" / "Supplementary_FigureS2_sequencing_QC"
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=PUBLICATION_DPI, bbox_inches="tight")
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

def _save_figure_formats(
    fig: plt.Figure,
    out_prefix: Path,
    aliases: tuple[str, ...] = (),
) -> None:
    """Save publication figures once, then copy byte-identical aliases."""
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    canonical_paths: list[Path] = []
    for suffix in (".png", ".pdf", ".svg"):
        output_path = out_prefix.with_suffix(suffix)
        kwargs: dict[str, Any] = {
            "bbox_inches": "tight",
            "facecolor": "white",
        }
        if suffix == ".png":
            kwargs["dpi"] = PUBLICATION_DPI
        fig.savefig(output_path, **kwargs)
        canonical_paths.append(output_path)

    # Re-rendering the same large figure for every historical filename is slow
    # and can leave a partially written raster if a job is interrupted.  The
    # aliases are therefore exact copies of the fully written canonical files.
    for alias in aliases:
        alias_prefix = out_prefix.with_name(alias)
        for source_path in canonical_paths:
            shutil.copy2(source_path, alias_prefix.with_suffix(source_path.suffix))


def create_main_figure(
    out_prefix: Path,
    panel_rows: list[dict[str, Any]],
    mechanistic_rows: list[dict[str, Any]],
    structural_by_sample: dict[str, dict[str, Any]],
    deletion_profile_rows: list[dict[str, Any]],
    inference_rows: list[dict[str, Any]] | None = None,
    parental_reference: ParentalReferenceModel | None = None,
    modbam_panels: dict[str, dict[str, Any]] | None = None,
) -> None:
    """Main publication Figure 1.

    A = complete-cohort heatmap plus one representative ModBAM profile per block
    B = cohort-wide chromosome-15 copy-number/deletion track
    C = deletion-span methylation profile classification against empirical
        maternal-like and paternal-like control profiles
    """
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
            "font.size": fs(9),
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": PUBLICATION_DPI,
            "savefig.facecolor": "white",
        }
    )

    fig = plt.figure(figsize=(19.0, 18.2), constrained_layout=False)
    outer = GridSpec(
        3, 1, figure=fig,
        height_ratios=[1.85, 1.0, 0.16],
        hspace=0.25,
    )

    panel_a = outer[0].subgridspec(
        3, 2,
        height_ratios=[0.08, 1.0, 0.14],
        width_ratios=[0.46, 0.54],
        hspace=0.05,
        wspace=0.09,
    )

    ax_a_header = fig.add_subplot(panel_a[0, :])
    ax_a_header.axis("off")
    ax_a_header.text(
        0.0, 0.62,
        "A. IC methylation states and representative single-molecule profiles",
        ha="left", va="center", fontsize=fs(11), fontweight="bold",
    )

    heat_grid = panel_a[1, 0].subgridspec(
        1, 2, width_ratios=[0.34, 1.0], wspace=0.05
    )
    ax_a_note = fig.add_subplot(heat_grid[0, 0])
    ax_a_heat = fig.add_subplot(heat_grid[0, 1])
    caption_grid = panel_a[2, :].subgridspec(
        1, 2, width_ratios=[0.46, 0.54], wspace=0.09
    )
    left_footer = caption_grid[0, 0].subgridspec(
        2, 1, height_ratios=[0.68, 0.32], hspace=0.12
    )
    ax_a_caption = fig.add_subplot(left_footer[0, 0])
    ax_a_cbar = fig.add_subplot(left_footer[1, 0])
    ax_modbam_caption = fig.add_subplot(caption_grid[0, 1])
    ax_modbam_caption.axis("off")
    draw_cohort_methylation_panel(
        ax_a_note,
        ax_a_heat,
        panel_rows,
        inference_rows,
        parental_reference,
        panel_label="A",
        cbar_ax=ax_a_cbar,
        cbar_orientation="horizontal",
        caption_ax=ax_a_caption,
    )
    if modbam_panels:
        draw_single_molecule_panel(
            fig, panel_a[1, 1], modbam_panels, MODBAM_PLOT_REGION, panel_rows,
            connector_ax=ax_a_heat,
            panel_label="A",
        )
    else:
        # A clear placeholder is preferable to silently dropping half of Panel A
        # in render-only workflows with incomplete cached ModBAM provenance.
        missing_ax = fig.add_subplot(panel_a[1, 1])
        missing_ax.axis("off")
        missing_ax.text(
            0.5, 0.5,
            "Representative ModBAM panels unavailable",
            ha="center", va="center", fontsize=fs(8), color="#777777",
        )
    ax_modbam_caption.text(
        0.0, 0.72,
        (
            f"One representative per cohort block; IC depth closest to the within-group median. "
            f"Direct MM/ML rendering, MAPQ ≥ {MODBAM_MIN_MAPQ}; {MODBAM_PLOT_REGION}."
        ),
        ha="left", va="top", fontsize=fs(5.5), color="#555555",
        wrap=True,
    )

    bottom = outer[1].subgridspec(
        1, 2, width_ratios=[0.52, 0.48], wspace=0.16
    )

    panel_b = bottom[0, 0].subgridspec(
        2, 1, height_ratios=[0.24, 1.0], hspace=0.025
    )
    ax_b_genes = fig.add_subplot(panel_b[0, 0])
    draw_panel_b_gene_track(ax_b_genes)
    ax_b = fig.add_subplot(panel_b[1, 0], sharex=ax_b_genes)
    draw_cohort_cnv_track(
        ax_b, structural_by_sample, show_legend=False, show_title=False
    )

    ax_c = fig.add_subplot(bottom[0, 1])
    draw_panel_c(
        ax_c,
        deletion_profile_rows,
    )

    # B and C encode the same molecular mechanisms. A single dedicated legend
    # band keeps both plotting areas unobstructed and preserves panel widths.
    legend_ax = fig.add_subplot(outer[2])
    legend_ax.axis("off")
    legend_handles = [
        Line2D(
            [0], [0],
            color=MECHANISM_COLORS[m], lw=2.2,
            marker=MECHANISM_MARKERS[m], markersize=6.5,
            markerfacecolor=MECHANISM_COLORS[m], markeredgecolor="white",
            label=m,
        )
        for m in sorted(MECHANISM_ORDER, key=MECHANISM_ORDER.get)
    ]
    mechanism_legend = legend_ax.legend(
        handles=legend_handles,
        title="Cohort / mechanism",
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.05, 1.0),
        ncol=len(legend_handles),
        fontsize=fs(7.2),
        title_fontsize=fs(7.4),
        handlelength=2.2,
        columnspacing=1.6,
    )
    legend_ax.add_artist(mechanism_legend)
    template_handles = [
        Line2D(
            [0], [0], marker="o", linestyle="none", markersize=5.5,
            markerfacecolor=STATE_COLORS["M"], markeredgecolor="#222222",
            label="M-like profile",
        ),
        Line2D([0], [0], marker="o", linestyle="none", markersize=5.5,
               markerfacecolor=STATE_COLORS["P"], markeredgecolor="#222222",
               label="P-like profile"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=5.5,
               markerfacecolor=STATE_COLORS["?"], markeredgecolor="#222222",
               label="Uncertain profile"),
    ]
    legend_ax.legend(
        handles=template_handles,
        title="Panel C classification",
        frameon=False,
        loc="upper right",
        bbox_to_anchor=(0.98, 1.0),
        ncol=3,
        fontsize=fs(6.8),
        title_fontsize=fs(7.0),
        columnspacing=1.0,
    )
    legend_ax.text(
        0.5, 0.02,
        "Panel C score = (RMSE to P − RMSE to M) / (RMSE to P + RMSE to M); positive is maternal-like, negative is paternal-like; | score | < 0.05 is uncertain.",
        ha="center", va="bottom", fontsize=fs(5.7), color="#666666",
    )

    fig.subplots_adjust(
        top=0.985,
        bottom=0.035,
        left=0.060,
        right=0.975,
    )

    # Earlier versions emitted byte-identical aliases. Remove only those exact
    # legacy stems so the figures directory contains one canonical Figure 1.
    for legacy_stem in ("Figure1_mechanistic", "Figure1_natural_experiment"):
        for suffix in (".png", ".pdf", ".svg"):
            legacy_path = out_prefix.parent / f"{legacy_stem}{suffix}"
            if legacy_path.exists():
                legacy_path.unlink()

    _save_figure_formats(
        fig,
        out_prefix,
    )
    plt.close(fig)

# ---------------------------------------------------------------------------
# Extensive report
# ---------------------------------------------------------------------------


def _md_escape(value: Any) -> str:
    """Escape values for compact Markdown tables."""
    if value is None:
        return ""
    return str(value).replace("|", "\\|").replace("\n", " ").strip()


def _md_table(headers: list[str], rows: list[list[Any]]) -> list[str]:
    """Return a Markdown table as a list of lines."""
    if not rows:
        return ["_No data available for this section._"]
    lines = [
        "| " + " | ".join(_md_escape(h) for h in headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append(
            "| " + " | ".join(_md_escape(v) for v in row) + " |"
        )
    return lines


def _report_image(
    report_path: Path,
    image_path: Path,
    alt_text: str,
    caption: str | None = None,
) -> list[str]:
    """Embed an image using a path relative to the Markdown report."""
    image_path = Path(image_path)
    if not image_path.exists():
        return [
            f"> **Image not available:** `{image_path}`",
        ]
    rel = Path(os.path.relpath(image_path, report_path.parent)).as_posix()
    lines = [f"![{_md_escape(alt_text)}]({rel})"]
    if caption:
        lines += ["", f"*{caption}*"]
    return lines


def _median_numeric(rows: list[dict[str, Any]], field: str) -> float | None:
    values = [
        value
        for row in rows
        if (value := safe_float(row.get(field))) is not None
    ]
    if not values:
        return None
    return float(np.median(values))


def _format_range(values: list[float], digits: int = 2) -> str:
    values = [v for v in values if np.isfinite(v)]
    if not values:
        return "NA"
    return f"{min(values):.{digits}f}–{max(values):.{digits}f}"


def _observed_state_from_panel_row(row: dict[str, Any]) -> str:
    codes: list[str] = []
    for idx in (1, 2):
        status = row.get(f"allele_{idx}_status", "missing")
        if status == "deleted":
            codes.append("absent")
        elif status != "observed":
            codes.append("?")
        else:
            code = row.get(f"allele_{idx}_pattern_short", "?")
            codes.append(code if code else "?")
    return " / ".join(
        normalize_state_pair((codes[0], codes[1]), row["molecular_mechanism"])
    )


def _report_main_findings(
    panel_rows: list[dict[str, Any]],
    support_rows: list[dict[str, Any]],
    structural_rows: list[dict[str, Any]],
    diagnostic_rows: list[dict[str, Any]],
    inference_rows: list[dict[str, Any]] | None,
    threshold_rows: list[dict[str, Any]] | None,
) -> list[str]:
    """Generate restrained, data-driven headline findings."""
    lines: list[str] = []

    expected_del = [
        r for r in structural_rows
        if str(r.get("expected_ic_deletion", "")).lower() in {"true", "1"}
        or r.get("molecular_mechanism") in {"PWS-DEL", "AS-DEL"}
    ]
    confirmed_del = [
        r for r in expected_del
        if r.get("ic_deletion_status") == "confirmed"
    ]
    if expected_del:
        lines.append(
            f"- **Structural dosage:** {len(confirmed_del)}/{len(expected_del)} "
            "expected PWS/AS deletion genomes had an IC-spanning deletion "
            "confirmed by HiFiCNV dosage and/or pbsv evidence."
        )

    for mechanism in ("PWS-DEL", "AS-DEL", "PWS-mUPD", "Disease control", "Control"):
        group = [
            r for r in diagnostic_rows
            if r.get("molecular_mechanism") == mechanism
        ]
        if not group:
            continue
        r = group[0]
        lines.append(
            f"- **{mechanism}:** observed state concordance "
            f"{r.get('state_concordant_n','')}/{r.get('n_samples','')} "
            f"({r.get('state_concordance_percent','')}%). "
            f"Observed distribution: {r.get('observed_state_distribution','')}."
        )

    if inference_rows:
        r = inference_rows[0]
        lines.append(
            "- **Reciprocal retained-allele contrast:** "
            f"PWS-DEL median β={r.get('median_PWS_DEL','NA')} versus "
            f"AS-DEL median β={r.get('median_AS_DEL','NA')}; "
            f"Δmedian={r.get('delta_median','NA')} "
            f"(participant-bootstrap 95% CI "
            f"{r.get('bootstrap_95CI_low','NA')} to "
            f"{r.get('bootstrap_95CI_high','NA')}); "
            f"exact two-sided permutation P="
            f"{r.get('exact_permutation_p_two_sided','NA')}."
        )

    if threshold_rows:
        concordances = [
            safe_float(r.get("cohort_concordance_percent"))
            for r in threshold_rows
        ]
        concordances = [x for x in concordances if x is not None]
        if concordances:
            lines.append(
                "- **Supplementary extreme-threshold sensitivity:** cohort concordance across "
                f"the prespecified threshold grid ranged from "
                f"{min(concordances):.1f}% to {max(concordances):.1f}%."
            )

    # Depth context by group.
    depth_parts = []
    for mechanism in sorted(MECHANISM_ORDER, key=MECHANISM_ORDER.get):
        group = [r for r in support_rows if r["molecular_mechanism"] == mechanism]
        med = _median_numeric(group, "bam_total_ic_depth")
        if med is not None:
            depth_parts.append(f"{mechanism} {med:.1f}×")
    if depth_parts:
        lines.append(
            "- **Sequencing-depth context:** median total IC depth was "
            + "; ".join(depth_parts)
            + ". Depth is therefore interpreted as a support variable rather "
              "than as a binary missingness criterion."
        )

    direct_modbam_cells = 0
    for row in panel_rows:
        for idx in (1, 2):
            source = str(row.get(f"allele_{idx}_source", ""))
            if "ModBAM MM/ML" in source:
                direct_modbam_cells += 1
    if direct_modbam_cells:
        lines.append(
            f"- **Direct ModBAM quantification:** {direct_modbam_cells} allele/haplotype "
            "estimate(s) shown in the cohort matrix were quantified directly "
            "from MM/ML+HP ModBAM evidence because the corresponding "
            "pb-CpG-tools haplotype BED did not provide an estimate."
        )

    return lines


def write_report(
    report_path: Path,
    mechanistic_rows: list[dict[str, Any]],
    diagnostic_rows: list[dict[str, Any]],
    modbam_provenance: list[dict[str, Any]],
    inference_rows: list[dict[str, Any]] | None = None,
    threshold_rows: list[dict[str, Any]] | None = None,
    panel_rows: list[dict[str, Any]] | None = None,
    support_rows: list[dict[str, Any]] | None = None,
    structural_rows: list[dict[str, Any]] | None = None,
    cn_classification_rows: list[dict[str, Any]] | None = None,
    preflight_rows: list[dict[str, Any]] | None = None,
    parental_reference: ParentalReferenceModel | None = None,
    parental_reference_rows: list[dict[str, Any]] | None = None,
    deletion_profile_rows: list[dict[str, Any]] | None = None,
    outdir: Path | None = None,
) -> None:
    """Write an extensive, figure-linked scientific results report.

    The report is intentionally generated from the result tables rather than
    from hard-coded conclusions. It contains:
      * headline findings,
      * cohort-level and sample-level methylation interpretation,
      * CN/SV structural evidence and deletion classes,
      * sequencing-depth and haplotype-support analysis,
      * participant-level statistical inference,
      * threshold sensitivity,
      * preflight/QC findings,
      * representative ModBAM provenance,
      * main, Extended Data and Supplementary images,
      * limitations and manuscript-ready interpretation.
    """
    report_path = Path(report_path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    outdir = Path(outdir) if outdir is not None else report_path.parent.parent

    panel_rows = panel_rows or []
    support_rows = support_rows or []
    structural_rows = structural_rows or []
    cn_classification_rows = cn_classification_rows or []
    preflight_rows = preflight_rows or []
    inference_rows = inference_rows or []
    threshold_rows = threshold_rows or []
    parental_reference_rows = parental_reference_rows or []
    deletion_profile_rows = deletion_profile_rows or []

    lines: list[str] = [
        "# Extensive Figure 1 results report",
        "",
        "## Haplotype-resolved chromosome-15 methylation, copy number and molecular configuration",
        "",
        f"**Reference:** T2T-CHM13v2.0  ",
        f"**PWS/AS imprinting centre:** {CHROM}:{PWS_IC_START:,}–{PWS_IC_END:,}  ",
        f"**Raw ModBAM display window:** {MODBAM_PLOT_REGION}  ",
        f"**Copy-number display window:** {CHROM}:{CN_PLOT_START:,}–{CN_PLOT_END:,}  ",
        "",
        "This report is generated automatically from the same result objects used "
        "to construct Figure 1 and its Supplementary/Extended Data outputs. "
        "The language below is intentionally constrained to what the current "
        "analysis supports: blood-derived methylation states, copy-number dosage, "
        "haplotype support and the prespecified participant-level comparisons.",
        "",
        "---",
        "",
        "## 1. Executive summary and main findings",
        "",
    ]
    lines += _report_main_findings(
        panel_rows,
        support_rows,
        structural_rows,
        diagnostic_rows,
        inference_rows,
        threshold_rows,
    )
    if parental_reference is not None:
        lines.append(
            "- **Control-calibrated parental-state model:** "
            f"paternal-like reference β={parental_reference.paternal_reference:.3f}, "
            f"maternal-like reference β={parental_reference.maternal_reference:.3f}, "
            f"equal-distance boundary β={parental_reference.decision_boundary:.3f}. "
            "Primary M/P calls use these empirical references and descriptive "
            "uncertainty, not fixed 0.85/0.15 thresholds."
        )
    retained_profile_rows = [
        row for row in deletion_profile_rows
        if row.get("track_role") == "retained deletion chromosome"
    ]
    if retained_profile_rows:
        called = [
            row for row in retained_profile_rows
            if row.get("profile_parental_class") in {"maternal-like", "paternal-like"}
        ]
        matching = [
            row for row in called
            if str(row.get("profile_matches_expected", "")).lower() == "true"
        ]
        uncertain = [
            row.get("display_label", row.get("sample_id", ""))
            for row in retained_profile_rows
            if row.get("profile_parental_class") == "uncertain"
        ]
        lines.append(
            "- **Deletion-span profile classifier:** "
            f"{len(called)}/{len(retained_profile_rows)} retained tracks received "
            f"a profile-level M/P call and {len(matching)}/{len(called)} called "
            "tracks matched the syndrome-expected parental state. "
            + (
                "Uncertain profiles: " + ", ".join(map(str, uncertain)) + "."
                if uncertain else "No retained profile fell in the uncertainty zone."
            )
        )

    # Main Figure.
    lines += [
        "",
        "### Main Figure 1",
        "",
    ]
    lines += _report_image(
        report_path,
        outdir / "figures" / "Figure1.png",
        "Main Figure 1",
        (
            "Main Figure 1 integrates the complete-cohort IC methylation heatmap "
            "with compact, depth-matched representative ModBAM profiles, cohort-wide "
            "chr15 dosage and the deletion-span parental-profile classifier. Lower "
            "per-haplotype depth is not treated as sample failure when the "
            "methylation state remains estimable."
        ),
    )

    # Cohort composition.
    counts = defaultdict(int)
    for row in panel_rows:
        counts[row["molecular_mechanism"]] += 1
    lines += [
        "",
        "## 2. Cohort and analytical design",
        "",
        "The analysis treats reciprocal deletion genomes as natural "
        "hemizygous configurations: PWS-DEL exposes the retained maternal-like "
        "chromosome, whereas AS-DEL exposes the retained paternal-like "
        "chromosome. PWS-mUPD represents a copy-neutral duplicated maternal "
        "state. The 22q11.2 deletion cohort serves as an orthogonal genomic-"
        "disorder control because chromosome 15 is expected to remain "
        "biparental.",
        "",
    ]
    lines += _md_table(
        ["Molecular mechanism", "n", "Expected IC state"],
        [
            [
                mechanism,
                counts.get(mechanism, 0),
                " / ".join(
                    normalize_state_pair(
                        GROUP_EXPECTED_STATE_CODES[mechanism], mechanism
                    )
                ),
            ]
            for mechanism in sorted(MECHANISM_ORDER, key=MECHANISM_ORDER.get)
        ],
    )

    # Structural evidence.
    lines += [
        "",
        "## 3. Structural and copy-number evidence at chromosome 15",
        "",
        "A deleted chromosome is not inferred from missing methylation. "
        "Hemizygous dosage is confirmed when HiFiCNV identifies a CN-loss "
        f"(configured threshold CN ≤ {CN_DELETION_THRESHOLD:.2f}) spanning "
        "the complete IC and/or pbsv provides an IC-spanning DEL. "
        "HiFiCNV and pbsv are therefore complementary: copy number establishes "
        "dosage, while sequence-resolved SV calls can provide breakpoint support.",
        "",
    ]

    class_counts = defaultdict(int)
    for row in cn_classification_rows:
        if row.get("molecular_mechanism") in {"PWS-DEL", "AS-DEL"}:
            class_counts[row.get("deletion_type", "unresolved")] += 1
    if class_counts:
        lines += ["### Deletion-class distribution", ""]
        for label, n in sorted(class_counts.items(), key=lambda x: (-x[1], x[0])):
            lines.append(f"- {label}: **{n}** sample(s)")

    lines += [
        "",
        "### Per-sample structural classification",
        "",
    ]
    lines += _md_table(
        [
            "Sample", "Group", "IC deletion", "Evidence basis",
            "Mean CN", "CN event (Mb)", "Deletion class", "pbsv"
        ],
        [
            [
                r.get("display_label", r.get("sample_id", "")),
                r.get("molecular_mechanism", ""),
                r.get("ic_deletion_status", ""),
                r.get("evidence_basis", ""),
                r.get("cn_event_mean_cn", ""),
                (
                    f"{r.get('cn_event_size_mb','')} Mb"
                    if str(r.get("cn_event_size_mb", "")).strip()
                    else ""
                ),
                r.get("deletion_type", ""),
                r.get("pbsv_support", ""),
            ]
            for r in cn_classification_rows
        ],
    )

    # Supplementary CN overview.
    lines += [
        "",
        "### Supplementary Figure: chr15 copy-number profiles across all participants",
        "",
    ]
    lines += _report_image(
        report_path,
        outdir / "supplementary" / "Supplementary_chr15_CN_all_participants.png",
        "Supplementary chr15 copy-number profiles",
        (
            "All participants are shown on identical chr15 axes. PWS-DEL and "
            "AS-DEL should exhibit CN≈1 across their deletion interval, whereas "
            "PWS-mUPD, disease controls and unaffected controls are expected to "
            "remain approximately diploid at chromosome 15."
        ),
    )

    # Per-sample CN images.
    lines += [
        "",
        "### Supplementary per-participant copy-number profiles",
        "",
        "The following panels are generated individually so that CN transitions, "
        "recurrent breakpoint landmarks and atypical/extended events can be "
        "inspected without compression from a cohort contact sheet.",
        "",
    ]
    labels = sample_display_labels()
    for sample_id, _clinical, mechanism in sorted_cohort():
        image = (
            outdir
            / "supplementary"
            / "chr15_copy_number_profiles"
            / f"{labels[sample_id]}_{sample_id}_chr15_CN.png"
        )
        structural = next(
            (
                r for r in cn_classification_rows
                if r.get("sample_id") == sample_id
            ),
            {},
        )
        caption = (
            f"**{labels[sample_id]} ({sample_id}; {mechanism}).** "
            f"Deletion classification: "
            f"{structural.get('deletion_type', 'not available')}. "
            f"Structural evidence: "
            f"{structural.get('evidence_basis', 'not available')}."
        )
        lines += [f"#### {labels[sample_id]} — {mechanism}", ""]
        lines += _report_image(
            report_path,
            image,
            f"{labels[sample_id]} chr15 copy-number profile",
            caption,
        )
        lines.append("")

    # Methylation states.
    lines += [
        "",
        "## 4. IC methylation architecture across the cohort",
        "",
        "The IC is used as an anchoring locus for maternal-like versus paternal-"
        "like methylation identity. Primary parental-state calls are calibrated "
        "to unaffected controls rather than fixed beta cutoffs: within each "
        "control the higher-methylated haplotype anchors the maternal-like state "
        "and the lower-methylated haplotype anchors the paternal-like state. The "
        "median high and low values define the cohort centroids, and their midpoint "
        "defines the equal-distance decision boundary. An allele is called M-like "
        "or P-like only when its descriptive CpG-bootstrap interval lies wholly "
        "on the corresponding side of that boundary; otherwise it is labelled "
        "uncertain. Controls are evaluated with a leave-one-control-out reference "
        "when possible. This is not presented as an independent blinded diagnostic "
        "assay. Confirmed deletions, missing evidence and low-depth observations "
        "remain separate categories.",
        "",
        "### Empirical parental-reference model",
        "",
        "",
    ]
    if parental_reference_rows:
        lines += _md_table(
            [
                "Scope", "Excluded control", "Controls used",
                "Paternal ref β", "Maternal ref β", "Decision boundary β",
                "Reference separation"
            ],
            [
                [
                    r.get("scope", ""),
                    r.get("excluded_control", ""),
                    r.get("controls_used", ""),
                    r.get("paternal_reference_beta", ""),
                    r.get("maternal_reference_beta", ""),
                    r.get("decision_boundary_beta", ""),
                    r.get("reference_separation", ""),
                ]
                for r in parental_reference_rows
            ],
        )
    lines += [
        "",
        "### Supplementary Figure: empirical parental-reference calibration",
        "",
    ]
    lines += _report_image(
        report_path,
        outdir / "supplementary" / "Supplementary_parental_reference_calibration.png",
        "Control-calibrated parental methylation reference scale",
        (
            "Observed allele/haplotype methylation values are displayed on the "
            "same continuous beta scale as the control-derived paternal-like and "
            "maternal-like centroids and their equal-distance boundary. This plot "
            "makes the primary classification rule visually auditable."
        ),
    )
    lines += [
        "",
        "### Per-sample methylation state",
        "",
    ]

    lines += _md_table(
        [
            "Sample", "Group", "Observed state",
            "Allele/HP1 β", "Allele/HP2 β",
            "Allele/HP1 source", "Allele/HP2 source"
        ],
        [
            [
                r.get("display_label", ""),
                r.get("molecular_mechanism", ""),
                _observed_state_from_panel_row(r),
                r.get("allele_1_mean_methylation", ""),
                r.get("allele_2_mean_methylation", ""),
                r.get("allele_1_source", ""),
                r.get("allele_2_source", ""),
            ]
            for r in panel_rows
        ],
    )

    deletion_rows = [
        row for row in deletion_profile_rows
        if row.get("track_role") == "retained deletion chromosome"
    ]
    validation_rows = [
        row for row in deletion_profile_rows
        if row.get("track_role") in {"reference", "held-out disease control"}
    ]
    lines += [
        "",
        "### Deletion-span parental-profile classification (Panel C)",
        "",
        "Panel C asks whether the methylation track retained inside each "
        "participant-specific CN deletion resembles a maternal or paternal "
        "chromosome. Unaffected-control haplotypes are first anchored as M-like "
        "or P-like by their IC state. At each deletion span, a CpG enters the "
        "profile classifier only when it is observed in all four control "
        "haplotypes, the M-like and P-like control medians differ by at least "
        f"{DELETION_PROFILE_MIN_REFERENCE_DELTA:.2f} beta, and the within-state "
        f"range is no greater than {DELETION_PROFILE_MAX_WITHIN_STATE_RANGE:.2f}. "
        "The retained profile is compared with the two position-matched reference "
        "profiles by RMSE. The plotted score is `(RMSE_P - RMSE_M) / "
        "(RMSE_P + RMSE_M)`: positive values are M-like and negative values are "
        f"P-like; absolute scores below {DELETION_PROFILE_MIN_SCORE_MAGNITUDE:.2f} "
        "are called uncertain. DiGeorge haplotypes are held out from training and "
        "retained in the source table as supplementary validation; Panel C itself "
        "shows only the unaffected-control chromosomes and deletion-sample tracks.",
        "",
        "#### Retained chromosome within each sample-specific deletion",
        "",
    ]
    lines += _md_table(
        [
            "Sample", "Group", "Deletion (Mb)", "Interval", "Informative CpGs",
            "RMSE M", "RMSE P", "Profile score", "Profile class", "IC class",
            "Expected", "Integrated interpretation",
        ],
        [
            [
                row.get("display_label", ""),
                row.get("molecular_mechanism", ""),
                row.get("deletion_size_mb", ""),
                f"{row.get('evaluation_start', '')}-{row.get('evaluation_end', '')}",
                row.get("n_shared_informative_CpGs", ""),
                row.get("rmse_to_maternal_profile", ""),
                row.get("rmse_to_paternal_profile", ""),
                row.get("parental_profile_score", ""),
                row.get("profile_parental_class", ""),
                row.get("IC_parental_class", ""),
                row.get("expected_parental_class", ""),
                row.get("integrated_classification", ""),
            ]
            for row in deletion_rows
        ],
    )
    lines += [
        "",
        "#### Separated control and DiGeorge haplotype validation",
        "",
    ]
    lines += _md_table(
        [
            "Sample", "Group", "Haplotype", "Informative CpGs",
            "Profile score", "Profile class", "IC class", "Concordance",
        ],
        [
            [
                row.get("display_label", ""),
                row.get("molecular_mechanism", ""),
                row.get("haplotype_label", ""),
                row.get("n_shared_informative_CpGs", ""),
                row.get("parental_profile_score", ""),
                row.get("profile_parental_class", ""),
                row.get("IC_parental_class", ""),
                row.get("integrated_classification", ""),
            ]
            for row in validation_rows
        ],
    )
    lines += [
        "",
        "The deletion-span result is a profile-similarity classification, not "
        "proof of biological parent of origin. It is internally calibrated from "
        "two unaffected controls, and long-range HP labels can be affected by "
        "phase-block boundaries. The IC call is therefore reported independently "
        "and disagreements are retained rather than forced into an M/P category.",
        "",
        "The control chromosomes lie close to scores of -1 and +1 because those "
        "same control tracks train the position-matched parental profiles. This "
        "is an in-sample reference property, not evidence that independent tracks "
        "should also reach the endpoints. Deletion-sample scores move toward zero "
        "when both RMSE values are substantial or similar, as can occur from "
        "inter-individual methylation variation, incomplete CpG overlap and "
        "long-range phase-block switching. A point can therefore be correctly "
        "closer to one parental profile while remaining far from its control "
        "centroid; its IC classification is reported separately.",
    ]

    lines += [
        "",
        "### Group-level state concordance",
        "",
    ]
    lines += _md_table(
        [
            "Group", "n", "Expected", "Observed distribution",
            "State concordance", "Template concordance"
        ],
        [
            [
                r.get("molecular_mechanism", ""),
                r.get("n_samples", ""),
                r.get("expected_state", ""),
                r.get("observed_state_distribution", ""),
                f"{r.get('state_concordant_n','')}/{r.get('n_samples','')} "
                f"({r.get('state_concordance_percent','')}%)",
                f"{r.get('template_concordant_n','')}/{r.get('n_samples','')} "
                f"({r.get('template_concordance_percent','')}%)",
            ]
            for r in diagnostic_rows
        ],
    )

    # Inference.
    lines += [
        "",
        "## 5. Participant-level statistical inference",
        "",
        "Formal inference treats **participants**, not CpGs or individual "
        "molecules, as independent biological replicates. CpG-level bootstrap "
        "intervals in the figure are descriptive uncertainty intervals only.",
        "",
    ]
    if inference_rows:
        r = inference_rows[0]
        lines += [
            f"- PWS-DEL participants: **n={r.get('n_PWS_DEL','')}**",
            f"- AS-DEL participants: **n={r.get('n_AS_DEL','')}**",
            f"- Median retained-allele methylation, PWS-DEL: "
            f"**β={r.get('median_PWS_DEL','NA')}**",
            f"- Median retained-allele methylation, AS-DEL: "
            f"**β={r.get('median_AS_DEL','NA')}**",
            f"- Difference in medians: **{r.get('delta_median','NA')}**",
            f"- Participant-bootstrap 95% CI: "
            f"**{r.get('bootstrap_95CI_low','NA')} to "
            f"{r.get('bootstrap_95CI_high','NA')}**",
            f"- Exact two-sided permutation P: "
            f"**{r.get('exact_permutation_p_two_sided','NA')}**",
            "",
            "The effect size and its participant-level consistency should be "
            "given greater interpretive weight than the nominal P value because "
            "the deletion subgroups are small.",
        ]
    else:
        lines.append("_Participant-level inference was not available._")

    # Supplementary sequencing-support context.
    lines += [
        "",
        "## 6. Supplementary sequencing depth and haplotype support",
        "",
        "Sequencing support is reported continuously because total IC depth, "
        "haplotype-specific depth and haplotype balance measure different aspects "
        "of the evidence. A state-estimable observation is retained even when it "
        "falls below the descriptive higher-coverage reference; only genuinely "
        "non-estimable evidence is counted as insufficient.",
        "",
        "### Group-level depth summary",
        "",
    ]
    depth_table = []
    for mechanism in sorted(MECHANISM_ORDER, key=MECHANISM_ORDER.get):
        group = [r for r in support_rows if r["molecular_mechanism"] == mechanism]
        total_values = [
            x for r in group
            if (x := safe_float(r.get("bam_total_ic_depth"))) is not None
        ]
        hpb = [
            x for r in group
            if (x := safe_float(r.get("hp_balance"))) is not None
        ]
        depth_table.append([
            mechanism,
            len(group),
            (
                f"{np.median(total_values):.1f}× "
                f"({_format_range(total_values,1)}×)"
                if total_values else "NA"
            ),
            (
                f"{np.median(hpb):.2f} ({_format_range(hpb,2)})"
                if hpb else "NA/hemizygous"
            ),
            sum(str(r.get("low_support", "")).lower() == "true" for r in group),
        ])
    lines += _md_table(
        ["Group", "n", "Median total IC depth (range)", "Median HP balance (range)", "Non-estimable n"],
        depth_table,
    )

    lines += [
        "",
        "### Per-sample depth and haplotype support",
        "",
    ]
    lines += _md_table(
        [
            "Sample", "Group", "Total depth", "HP1 depth", "HP2 depth",
            "HP balance", "HP1 CpGs", "HP2 CpGs", "Combined CpGs",
            "Support tier"
        ],
        [
            [
                r.get("display_label", ""),
                r.get("molecular_mechanism", ""),
                r.get("bam_total_ic_depth", ""),
                r.get("hp1_ic_depth", ""),
                r.get("hp2_ic_depth", ""),
                r.get("hp_balance", ""),
                r.get("hap1_cpgs", ""),
                r.get("hap2_cpgs", ""),
                r.get("combined_cpgs", ""),
                r.get("support_tier", ""),
            ]
            for r in support_rows
        ],
    )

    # Threshold sensitivity.
    lines += [
        "",
        "## 7. Supplementary extreme-threshold sensitivity",
        "",
        "Fixed beta thresholds are not used for the primary parental-state calls. "
        "This prespecified 0.80/0.20, 0.85/0.15 and 0.90/0.10 grid is retained "
        "only as a secondary sensitivity analysis to show how an extreme-state "
        "classifier behaves relative to the control-calibrated continuous model.",
        "",
    ]
    lines += _md_table(
        ["M threshold", "P threshold", "Concordant", "Cohort n", "Concordance"],
        [
            [
                r.get("maternal_threshold", ""),
                r.get("paternal_threshold", ""),
                r.get("cohort_concordant_n", ""),
                r.get("cohort_n", ""),
                f"{r.get('cohort_concordance_percent','')}%",
            ]
            for r in threshold_rows
        ],
    )

    # Preflight.
    lines += [
        "",
        "## 8. Scientific preflight and QC findings",
        "",
    ]
    status_counts = defaultdict(int)
    for r in preflight_rows:
        status_counts[r.get("status", "UNKNOWN")] += 1
    if status_counts:
        lines.append(
            "- " + "; ".join(
                f"**{status}: {n}**"
                for status, n in sorted(status_counts.items())
            )
        )
        lines.append("")

    lines += _md_table(
        ["Sample", "Group", "Status", "Deletion status", "Deletion type", "Haplotype support", "Support notes", "Issues"],
        [
            [
                r.get("sample_id", ""),
                r.get("molecular_mechanism", ""),
                r.get("status", ""),
                r.get("ic_deletion_status", ""),
                r.get("deletion_type", ""),
                r.get("haplotype_support", ""),
                r.get("support_notes", ""),
                r.get("issues", ""),
            ]
            for r in preflight_rows
        ],
    )

    # Raw molecules and provenance.
    lines += [
        "",
        "## 9. Raw single-molecule methylation evidence",
        "",
        f"All main-panel ModBAM visualizations use {MODBAM_PLOT_REGION}, "
        f"primary alignments and MAPQ ≥ {MODBAM_MIN_MAPQ}. Representative "
        "selection is prespecified as the sample whose total IC depth is closest "
        "to the within-group median, reducing the risk of visual cherry-picking.",
        "",
        "### Main-panel representative provenance",
        "",
    ]
    lines += _md_table(
        ["Group", "Sample", "Selection rule", "Sample depth", "Group median depth", "Grouping"],
        [
            [
                r.get("molecular_mechanism", ""),
                r.get("display_label", ""),
                r.get("selection_rule", ""),
                r.get("sample_total_ic_depth", ""),
                r.get("group_median_ic_depth", ""),
                r.get("haplotype_grouping", ""),
            ]
            for r in modbam_provenance
        ],
    )

    lines += [
        "",
        "### Extended Data: ModBAM profiles for all participants",
        "",
    ]
    lines += _report_image(
        report_path,
        outdir / "extended_data" / "ExtendedData_Figure1_all_modbamtools.png",
        "Extended Data ModBAM profiles for all participants",
        (
            "All generated raw ModBAM profiles are shown together to expose "
            "between-participant heterogeneity and to ensure that the main-panel "
            "representatives are not the only visual evidence."
        ),
    )

    individual_modbam_dir = (
        outdir / "extended_data" / "modbamtools_all_samples"
    )
    individual_modbam_pngs = (
        sorted(individual_modbam_dir.glob("*.png"))
        if individual_modbam_dir.exists()
        else []
    )
    if individual_modbam_pngs:
        lines += [
            "",
            "### Extended Data: individual ModBAM profiles",
            "",
            "Individual raw-molecule panels are embedded below to permit "
            "sample-level inspection without relying only on the cohort contact "
            "sheet.",
            "",
        ]
        for image in individual_modbam_pngs:
            lines += [f"#### {image.stem}", ""]
            lines += _report_image(
                report_path,
                image,
                image.stem,
                (
                    "Single-participant ModBAM profile generated with the "
                    f"standardized {MODBAM_PLOT_REGION} window and MAPQ ≥ "
                    f"{MODBAM_MIN_MAPQ} filtering."
                ),
            )
            lines.append("")

    # Interpretative synthesis.
    lines += [
        "",
        "## 10. Integrated biological interpretation",
        "",
        "### Finding 1 — reciprocal deletions provide direct parental-state exposure",
        "",
        "When the structural deletion is confirmed, PWS-DEL and AS-DEL provide "
        "complementary hemizygous configurations in which one parental chromosome "
        "is physically absent. This makes the retained IC methylation state "
        "directly interpretable as maternal-like in PWS-DEL or paternal-like in "
        "AS-DEL, without requiring a two-haplotype decomposition inside the "
        "deleted interval.",
        "",
        "### Finding 2 — copy number and methylation answer different questions",
        "",
        "HiFiCNV establishes dosage (for example, CN≈1 across the pathogenic "
        "interval), whereas methylation establishes the epigenetic state of the "
        "retained chromosome. pbsv adds sequence-resolved structural support when "
        "a breakpoint-spanning DEL is recovered. Failure of pbsv to emit one "
        "canonical large DEL is therefore not equivalent to absence of a deletion "
        "when the CN dosage signal is clear.",
        "",
        "### Finding 3 — depth differences must not be mistaken for biology",
        "",
        "The ability to resolve both haplotypes depends on total IC depth and HP "
        "balance. Disease-control or other lower-depth samples can therefore have "
        "valid phase blocks but sparser haplotype-specific methylation evidence. "
        "The depth-aware logic retains state-estimable observations, reports their "
        "support continuously and can quantify sparse haplotype BEDs directly from "
        "MM/ML+HP ModBAM evidence.",
        "",
        "### Finding 4 — orthogonal disease controls test locus specificity",
        "",
        "The 22q11.2 deletion cohort carries a pathogenic genomic deletion outside "
        "chromosome 15. A preserved biparental M/P configuration on chromosome 15 "
        "therefore argues that the PWS/AS IC pattern is locus-specific rather than "
        "a generic consequence of carrying a large pathogenic deletion.",
        "",
        "### Finding 5 — deletion classes should be interpreted as CN-transition classes",
        "",
        "Type I-like, Type II-like and atypical/extended labels are derived from "
        "HiFiCNV transition coordinates relative to T2T BP landmarks. These labels "
        "are appropriate for dosage architecture but should not be described as "
        "base-pair-resolved breakpoint assignments unless junction-spanning or "
        "assembly evidence independently resolves the breakpoints.",
    ]

    # Limitations.
    lines += [
        "",
        "## 11. Limitations and claims that should remain bounded",
        "",
        "- **Tissue:** these data are blood-derived and do not establish neuronal "
        "chromatin architecture or neuronal UBE3A/SNORD116 regulation.",
        "- **Sample size:** participant-level PWS-DEL and AS-DEL inference is based "
        "on small groups; effect size, confidence interval and directional "
        "consistency are more informative than a highly significant P value.",
        "- **Parental anchoring:** IC methylation anchors parental-like identity. "
        "Unless parental genotypes or an independent assay are available, the same "
        "IC should not be described as an independent validation of the labels it "
        "helps define.",
        "- **Reference calibration:** the empirical maternal-like and paternal-like "
        "centroids are internally calibrated from only two unaffected controls. "
        "They are appropriate anchors for this cohort but should not be presented "
        "as population-wide diagnostic cutoffs. Leave-one-control-out evaluation "
        "reduces, but does not eliminate, the limitations of the small reference set.",
        "- **CN breakpoint precision:** HiFiCNV defines dosage transitions, not "
        "necessarily nucleotide-resolution breakpoints inside segmental "
        "duplications.",
        "- **CpG-bootstrap intervals:** the displayed CpG bootstrap describes "
        "measurement uncertainty within an allele/haplotype; it does not increase "
        "the biological sample size.",
        "- **Direct ModBAM estimates:** ModBAM-derived methylation estimates are "
        "explicitly tagged and should be sensitivity-checked against BED-derived "
        "estimates where both are available.",
    ]

    # Manuscript-ready take-home.
    lines += [
        "",
        "## 12. Manuscript-ready take-home statement",
        "",
        "> Reciprocal PWS and Angelman deletions provide natural hemizygous "
        "configurations that expose opposite parental methylation states at the "
        "PWS/AS imprinting centre. Integrating HiFiCNV dosage, pbsv structural "
        "support, phased MM/ML-tagged HiFi molecules and depth-aware methylation "
        "quantification separates true chromosome loss from incomplete "
        "haplotype resolution and provides a reproducible framework for "
        "parent-of-origin methylation mapping.",
        "",
        "This statement should be revised if the automatically generated tables "
        "show incomplete structural confirmation, unresolved controls, or "
        "threshold-sensitive state assignments.",
    ]

    # Reproducibility index.
    lines += [
        "",
        "## 13. Reproducibility and source-data index",
        "",
        "- `../../../analysis/01_evidence_matrix/chr15_structural_evidence.tsv` — canonical analysis-derived HiFiCNV/pbsv structural evidence consumed by Figure 1.",
        "- `../tables/Figure1_parental_reference_model.tsv` — control-derived parental centroids, decision boundary and leave-one-control-out references.",
        "- `../tables/Figure1A_modbamtools_representatives.tsv` — representative raw single-molecule panels.",
        "- `../tables/Figure1B_allele_methylation_matrix.tsv` — plotted cohort methylation values and states.",
        "- `../tables/Figure1B_sample_level_inference.tsv` — participant-level PWS-DEL versus AS-DEL inference.",
        "- `../tables/Figure1C_deletion_span_parental_profile.tsv` — sample-specific deletion-span RMSE classification and separated control/DiGeorge haplotype validation.",
        "- `../tables/Figure1C_coverage_phasing_support.tsv` — depth, HP depth, CpG support and HP balance.",
        "- `../supplementary/Supplementary_chr15_deletion_classification.tsv` — figure-local snapshot of the canonical analysis-derived deletion classification.",
        "- `../supplementary/Supplementary_parental_reference_calibration.png` — continuous control-calibrated parental reference visualization.",
        "- `../supplementary/Figure1_threshold_sensitivity.tsv` — sensitivity to M/P thresholds.",
        "- `../source_data/Figure1A_single_molecule_MM_ML_source_data.tsv.gz` — raw long-format MM/ML calls used as source data.",
        "- `../source_data/Figure1B_cohort_methylation_source_data.tsv` — source data for the cohort methylation matrix.",
        "- `../source_data/Figure1C_deletion_span_parental_profile_source_data.tsv` — source data for Panel C profile scores.",
        "- `../source_data/Figure1C_sequencing_support_source_data.tsv` — source data for sequencing/phasing support.",
        "- `../source_data/Supplementary_chr15_copy_number_segments_source_data.tsv` — CN segments underlying supplementary CN plots.",
        "- `../Figure1_configuration.json` — exact thresholds and genomic intervals.",
        "- `../Figure1_software_versions.json` — software versions and Git commit.",
        "- `../Figure1_input_manifest.tsv` — input paths, file sizes and available checksums.",
        "",
        "---",
        "",
        "_End of automatically generated Figure 1 report._",
    ]

    report_path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    validate_configuration()
    outdir = Path(OUTDIR); table_dir = outdir / "tables"; figure_dir = outdir / "figures"; report_dir = outdir / "reports"
    source_dir = outdir / "source_data"; supp_dir = outdir / "supplementary"
    for d in (table_dir, figure_dir, report_dir, source_dir, supp_dir): d.mkdir(parents=True, exist_ok=True)

    if RENDER_ONLY:
        panel_path = table_dir / "Figure1B_allele_methylation_matrix.tsv"
        if not panel_path.exists():
            panel_path = table_dir / "Figure1C_allele_methylation_matrix.tsv"
        if not panel_path.exists():
            panel_path = table_dir / "Figure1A_allele_methylation_matrix.tsv"
        support_path = table_dir / "Figure1C_coverage_phasing_support.tsv"
        if not support_path.exists():
            support_path = table_dir / "SupplementaryFigure1_coverage_phasing_support.tsv"
        if not support_path.exists():
            support_path = table_dir / "Figure1D_coverage_phasing_support.tsv"
        structural_path = Path(STRUCTURAL_EVIDENCE_PATH)
        mechanistic_path = supp_dir / "closest_methylation_template_per_sample.tsv"
        deletion_profile_path = table_dir / "Figure1C_deletion_span_parental_profile.tsv"
        missing=[p for p in (panel_path,support_path,structural_path,mechanistic_path,deletion_profile_path) if not p.exists()]
        if missing: raise FileNotFoundError("Render-only mode missing:\n"+"\n".join(map(str,missing)))
        panel_rows=read_tsv(panel_path); support_rows=read_tsv(support_path)
        mechanistic_rows=read_tsv(mechanistic_path)
        deletion_profile_rows=read_tsv(deletion_profile_path)
        structural_cached, structural_rows_cached = load_analysis_structural_evidence(
            structural_path
        )
        inference_path=table_dir/"Figure1B_sample_level_inference.tsv"
        if not inference_path.exists():
            inference_path=table_dir/"Figure1C_sample_level_inference.tsv"
        inference=read_tsv(inference_path) if inference_path.exists() else []
        reference_path=table_dir/"Figure1_parental_reference_model.tsv"
        render_reference=None
        if reference_path.exists():
            ref_rows=read_tsv(reference_path)
            global_rows=[r for r in ref_rows if r.get("scope")=="global"]
            if global_rows:
                rr=global_rows[0]
                m=safe_float(rr.get("maternal_reference_beta")); p=safe_float(rr.get("paternal_reference_beta")); b=safe_float(rr.get("decision_boundary_beta"))
                if m is not None and p is not None and b is not None:
                    render_reference=ParentalReferenceModel(m,p,b,tuple((rr.get("controls_used") or "").split(";")) if rr.get("controls_used") else tuple())
        render_modbam_panels = None
        provenance_path = table_dir / "Figure1A_modbamtools_representatives.tsv"
        legacy_provenance_path = table_dir / "Figure1B_modbamtools_representatives.tsv"
        cached_provenance_path = (
            provenance_path if provenance_path.exists() else legacy_provenance_path
        )
        if cached_provenance_path.exists():
            try:
                render_modbam_panels, _ = load_modbamtools_panels_from_provenance(
                    cached_provenance_path, outdir
                )
            except (FileNotFoundError, RuntimeError) as exc:
                print(f"[WARN] Could not load cached ModBAM panels: {exc}", file=sys.stderr)
        create_main_figure(
            figure_dir / "Figure1",
            panel_rows,
            mechanistic_rows,
            structural_cached,
            deletion_profile_rows,
            inference,
            render_reference,
            modbam_panels=render_modbam_panels,
        )
        render_supplementary_sequencing_qc(support_rows, outdir)
        return

    vcf_dir=Path(VCF_DIR); bam_dir=Path(BAM_DIR); modbam_dir=Path(MODBAM_DIR); methylation_dir=Path(METHYLATION_DIR); cnv_dir=Path(CNV_DIR); metadata_path=Path(METADATA_PATH)
    metadata=read_metadata(metadata_path)
    sample_files={}
    for sample_id,_clinical,_mechanism in sorted_cohort():
        sv_vcf = (find_sample_file(vcf_dir,sample_id,".sv.phased.vcf.gz") or find_sample_file(vcf_dir,sample_id,".sv.pass.vcf.gz") or find_sample_file(vcf_dir,sample_id,".sv.vcf.gz"))
        sample_files[sample_id]={
            "bam":find_sample_file(bam_dir,sample_id,".bam"), "modbam":find_modbam_file(modbam_dir,sample_id),
            "blocks":find_sample_file(vcf_dir,sample_id,".blocks.tsv"), "sv_vcf":sv_vcf,
            "combined_bed":find_sample_file(methylation_dir,sample_id,f"{PBCPG_OUTPUT_TAG}.combined.bed"),
            "hap1_bed":find_sample_file(methylation_dir,sample_id,f"{PBCPG_OUTPUT_TAG}.hap1.bed"), "hap2_bed":find_sample_file(methylation_dir,sample_id,f"{PBCPG_OUTPUT_TAG}.hap2.bed"),
            "combined_bw":find_sample_file(methylation_dir,sample_id,f"{PBCPG_OUTPUT_TAG}.combined.bw"),
            "hap1_bw":find_sample_file(methylation_dir,sample_id,f"{PBCPG_OUTPUT_TAG}.hap1.bw"), "hap2_bw":find_sample_file(methylation_dir,sample_id,f"{PBCPG_OUTPUT_TAG}.hap2.bw"),
            "cnv_log":find_sample_file(cnv_dir,sample_id,".log"),
            "cnv_bed":find_sample_file(cnv_dir,sample_id,".cnv.bed"),
            "cn_track":find_hificnv_cn_track(cnv_dir,sample_id),
        }

    # Structural/deletion inference is performed once by scripts/analysis.
    # Figure 1 consumes the canonical result and never reclassifies CN/SV status.
    structural, cn_classification_rows = load_analysis_structural_evidence(
        Path(STRUCTURAL_EVIDENCE_PATH)
    )

    # Keep a figure-local snapshot for manuscript packaging/backward compatibility.
    # This is a copy of the analysis result, not a Figure 1-derived classification.
    write_tsv(
        supp_dir/"Supplementary_chr15_deletion_classification.tsv",
        cn_classification_rows,
    )
    write_tsv(
        source_dir/"Supplementary_chr15_copy_number_segments_source_data.tsv",
        build_cn_segment_source_rows(sample_files),
    )
    render_supplementary_cn_profiles(sample_files, structural, outdir)

    bam_qc_cache=table_dir/"bam_qc_cache.tsv"
    if SKIP_BAM_QC and bam_qc_cache.exists(): bam_qc_rows=read_tsv(bam_qc_cache)
    else:
        bam_qc_rows=[build_bam_qc(sid,sample_files[sid]["bam"],sample_files[sid]["cnv_log"]) for sid,_c,_m in sorted_cohort()]
        write_tsv(bam_qc_cache,bam_qc_rows)
    bam_qc_by_sample={r["sample_id"]:r for r in bam_qc_rows}

    summary_rows=[]
    for sid,clinical,mech in sorted_cohort():
        meta=metadata.get(sid,{}); n50,domain_fraction=block_n50_and_domain_fraction(sample_files[sid]["blocks"]); q=bam_qc_by_sample.get(sid,{})
        summary_rows.append({
            "sample_id":sid,"clinical_diagnosis":clinical,"molecular_mechanism":mech,
            "sex":meta.get("gender",meta.get("sex","")),"age_at_sampling":meta.get("age",""),
            "total_HiFi_reads":q.get("total_HiFi_reads",""),
            "mean_depth_genome_wide_raw_hificnv":fmt(q.get("mean_depth_genome_wide_raw_hificnv")),
            "mean_depth_genome_wide":fmt(q.get("mean_depth_genome_wide")),"mean_depth_chr15":fmt(q.get("mean_depth_chr15")),
            "mean_depth_per_haplotype_at_15q11-q13":q.get("mean_depth_per_haplotype_at_15q11-q13",""),
            "phasing_block_N50_chr15":n50 if n50 is not None else "","percent_imprinted_domain_in_phased_block":fmt(domain_fraction),
        })
    write_tsv(table_dir/"Figure1_cohort_QC_summary.tsv",summary_rows)

    assignment_rows,matrix_rows,stats_by_sample,parental_reference,parental_reference_rows=build_assignments(sample_files,structural)
    write_tsv(table_dir/"Figure1_parental_like_assignment.tsv",assignment_rows)
    write_tsv(table_dir/"Figure1_IC_methylation_matrix.tsv",matrix_rows)
    write_tsv(table_dir/"Figure1_parental_reference_model.tsv",parental_reference_rows)
    panel_rows=build_physical_allele_rows(matrix_rows,structural)
    write_tsv(table_dir/"Figure1B_allele_methylation_matrix.tsv",panel_rows)
    deletion_profile_rows=build_deletion_profile_classification_rows(
        sample_files, structural, panel_rows
    )
    write_tsv(
        table_dir/"Figure1C_deletion_span_parental_profile.tsv",
        deletion_profile_rows,
    )
    render_parental_reference_calibration(panel_rows,parental_reference,outdir)

    preflight=build_preflight_qc(sample_files,structural,stats_by_sample)
    write_tsv(table_dir/"Figure1_preflight_QC.tsv",preflight)
    enforce_preflight(preflight)

    mechanistic_rows=build_mechanistic_state_rows(panel_rows,parental_reference); diagnostic_rows=build_diagnostic_state_rows(mechanistic_rows)
    write_tsv(supp_dir/"closest_methylation_template_per_sample.tsv",mechanistic_rows)
    write_tsv(supp_dir/"methylation_template_group_summary.tsv",diagnostic_rows)

    threshold_rows=build_threshold_sensitivity(panel_rows); write_tsv(supp_dir/"Figure1_threshold_sensitivity.tsv",threshold_rows)
    inference_rows=build_sample_level_inference(panel_rows); write_tsv(table_dir/"Figure1B_sample_level_inference.tsv",inference_rows)

    # Legacy mixed contrast is retained only outside the Figure1 panel namespace.
    contrast_rows=build_per_cpg_contrast(stats_by_sample,assignment_rows)
    write_tsv(supp_dir/"legacy_mixed_parental_state_contrast.tsv",contrast_rows)

    support_rows=build_support_rows(summary_rows,matrix_rows,sample_files,structural)
    write_tsv(table_dir/"Figure1C_coverage_phasing_support.tsv",support_rows)

    # Supplementary figures now carry the per-sample measurement detail and
    # technical sequencing/QC evidence, keeping the main Figure 1 biologically focused.
    render_supplementary_methylation_tracks(sample_files,matrix_rows,outdir)
    render_supplementary_sequencing_qc(support_rows,outdir)

    provenance_path=table_dir/"Figure1A_modbamtools_representatives.tsv"
    if SKIP_MODBAMTOOLS:
        cached_provenance_path = provenance_path
        legacy_provenance_path = table_dir/"Figure1B_modbamtools_representatives.tsv"
        if not cached_provenance_path.exists() and legacy_provenance_path.exists():
            cached_provenance_path = legacy_provenance_path
        modbam_panels,modbam_provenance=load_modbamtools_panels_from_provenance(cached_provenance_path,outdir)
        if cached_provenance_path != provenance_path:
            write_tsv(provenance_path,modbam_provenance)
    else:
        representatives=choose_representative_modbam_samples(support_rows,sample_files)
        modbam_panels,modbam_provenance=build_modbamtools_panels(representatives,outdir,MODBAMTOOLS_BIN,MODBAM_PLOT_REGION,MODBAM_GTF)
        # Add explicit standardized filtering settings to provenance.
        for r in modbam_provenance:
            r["MAPQ_min"]=MODBAM_MIN_MAPQ; r["primary_alignment_filter"]="-F 2308"; r["plot_width_px"]=MODBAM_PLOT_WIDTH
        write_tsv(provenance_path,modbam_provenance)

    if GENERATE_EXTENDED_MODBAM and not SKIP_MODBAMTOOLS:
        ext_rows=build_all_sample_modbamtools_panels(sample_files,outdir)
        write_tsv(outdir/"extended_data"/"ExtendedData_Figure1_modbam_provenance.tsv",ext_rows)
        render_extended_modbam_contact_sheet(outdir,ext_rows)

    if GENERATE_READ_LEVEL_SOURCE_DATA:
        read_rows=[]
        for sid,_clinical,mech in sorted_cohort():
            mb=sample_files[sid].get("modbam")
            if mb and Path(mb).exists():
                try: read_rows.extend(extract_modbam_source_data(sid,mech,Path(mb),MODBAM_PLOT_REGION))
                except Exception as exc: print(f"[WARN] source-data extraction failed for {sid}: {exc}",file=sys.stderr)
        write_tsv_gz(source_dir/"Figure1A_single_molecule_MM_ML_source_data.tsv.gz",read_rows)

    # Nature-style minimum underlying source data.
    write_tsv(source_dir/"Figure1B_cohort_methylation_source_data.tsv",panel_rows)
    write_tsv(source_dir/"Figure1C_deletion_span_parental_profile_source_data.tsv",deletion_profile_rows)
    write_tsv(source_dir/"Figure1_parental_reference_model_source_data.tsv",parental_reference_rows)
    write_tsv(source_dir/"Figure1C_sequencing_support_source_data.tsv",support_rows)
    write_tsv(source_dir/"Supplementary_FigureS2_sequencing_QC_source_data.tsv",support_rows)

    create_main_figure(
        figure_dir / "Figure1",
        panel_rows,
        mechanistic_rows,
        structural,
        deletion_profile_rows,
        inference_rows,
        parental_reference,
        modbam_panels=modbam_panels,
    )
    write_report(
        report_dir / "Figure1_report.md",
        mechanistic_rows=mechanistic_rows,
        diagnostic_rows=diagnostic_rows,
        modbam_provenance=modbam_provenance,
        inference_rows=inference_rows,
        threshold_rows=threshold_rows,
        panel_rows=panel_rows,
        support_rows=support_rows,
        structural_rows=list(structural.values()),
        cn_classification_rows=cn_classification_rows,
        preflight_rows=preflight,
        parental_reference=parental_reference,
        parental_reference_rows=parental_reference_rows,
        deletion_profile_rows=deletion_profile_rows,
        outdir=outdir,
    )

    manifest=build_input_manifest(sample_files,metadata_path); write_tsv(outdir/"Figure1_input_manifest.tsv",manifest)
    with (outdir/"Figure1_software_versions.json").open("w") as h: json.dump(software_versions(),h,indent=2)
    run_parameters={
        "reference":"T2T-CHM13v2.0","cohort":[{"sample_id":s,"diagnosis":c,"mechanism":m} for s,c,m in sorted_cohort()],
        "regions":{"domain":[CHROM,DOMAIN_START,DOMAIN_END],"PWS_AS_IC":[CHROM,PWS_IC_START,PWS_IC_END],"modbam_display":MODBAM_PLOT_REGION},
        "figure_layout":{
            "panel_a":"complete-cohort IC methylation heatmap plus one vertically matched ModBAM representative per cohort block",
            "panel_b":"cohort-wide chr15 copy-number/deletion track with curated gene and PWS/AS ICR annotation",
            "panel_c":"each retained deletion-sample methylation track compared directly with unaffected-control maternal-like and paternal-like chromosome profiles",
            "shared_legend":"external legend band below panels B and C",
            "png_dpi":PUBLICATION_DPI,
            "supplementary_S1":"per-sample CpG methylation tracks around the IC",
            "supplementary_S2":"sequencing depth and haplotype-resolution QC",
            "extended_data":"raw ModBAM single-molecule profiles when generated",
        },
        "parental_state_model":{
            "method":PARENTAL_REFERENCE_METHOD,
            "controls_used":list(parental_reference.control_sample_ids),
            "maternal_reference_beta":parental_reference.maternal_reference,
            "paternal_reference_beta":parental_reference.paternal_reference,
            "decision_boundary_beta":parental_reference.decision_boundary,
            "leave_one_control_out_for_controls":CONTROL_REFERENCE_LEAVE_ONE_OUT,
            "classification_rule":"M if descriptive CpG-bootstrap CI is entirely above control-derived boundary; P if entirely below; otherwise uncertain; nearest centroid fallback if CI unavailable",
        },
        "deletion_span_profile_model":{
            "reference_training":"unaffected controls only; haplotypes anchored by IC state",
            "validation":"DiGeorge haplotypes held out and evaluated on the shared deletion core",
            "minimum_control_M_vs_P_delta":DELETION_PROFILE_MIN_REFERENCE_DELTA,
            "maximum_within_state_control_range":DELETION_PROFILE_MAX_WITHIN_STATE_RANGE,
            "minimum_shared_informative_CpGs":DELETION_PROFILE_MIN_SHARED_CPGS,
            "score":"(RMSE_P - RMSE_M) / (RMSE_P + RMSE_M)",
            "uncertain_if_absolute_score_below":DELETION_PROFILE_MIN_SCORE_MAGNITUDE,
        },
        "supplementary_extreme_threshold_sensitivity":{
            "legacy_maternal_threshold":EXTREME_MATERNAL_THRESHOLD,
            "legacy_paternal_threshold":EXTREME_PATERNAL_THRESHOLD,
            "grid":THRESHOLD_SENSITIVITY,
        },
        "technical_thresholds":{
            "state_min_mean_coverage":MIN_STATE_MEAN_COVERAGE,
            "state_min_CpGs":MIN_STATE_CPGS,
            "higher_coverage_reference_mean_coverage":HIGHER_COVERAGE_MEAN_COVERAGE,
            "higher_coverage_reference_CpGs":HIGHER_COVERAGE_CPGS,
            "interpretation":"Values below the higher-coverage reference remain valid when they meet state-assignment minima; they are not sample failures.",
            "modbam_min_MAPQ":MODBAM_MIN_MAPQ,
            "modbam_plot_width_px":MODBAM_PLOT_WIDTH,
            "modbam_fallback_min_molecules":MIN_MODBAM_FALLBACK_MOLECULES,
        },
        "copy_number":{
            "deletion_threshold":CN_DELETION_THRESHOLD,
            "normal_expected":CN_NORMAL_EXPECTED,
            "plot_interval":[CHROM,CN_PLOT_START,CN_PLOT_END],
            "breakpoint_landmarks":BREAKPOINT_LANDMARKS,
            "nearest_landmark_max_distance_bp":BP_NEAREST_MAX_DISTANCE,
        },
        "bootstrap":{"CpG_descriptive_replicates":CPG_BOOTSTRAP_REPLICATES,"participant_replicates":SAMPLE_BOOTSTRAP_REPLICATES,"seed":BOOTSTRAP_SEED},
        "structural_rule":"HiFiCNV CN-loss spanning the IC or pbsv DEL spanning the IC confirms hemizygous dosage; pbsv is orthogonal breakpoint support, not a mandatory gate.",
        "input_paths":{"vcf_dir":str(vcf_dir),"bam_dir":str(bam_dir),"modbam_dir":str(modbam_dir),"methylation_dir":str(methylation_dir),"cnv_dir":str(cnv_dir),"metadata":str(metadata_path),"structural_evidence":str(Path(STRUCTURAL_EVIDENCE_PATH))},
    }
    with (outdir/"Figure1_configuration.json").open("w") as h: json.dump(run_parameters,h,indent=2)



if __name__ == "__main__":
    main()
