"""
Every path, region and threshold used by scripts/duplicons/ lives here.

Reference and program locations are read from the same params.local.yml the
Nextflow workflow uses (keys: reference, cpg_model, samtools_bin, pbcpgtools_bin,
minimap2_bin, meryl_bin, repeatmasker_bin). Missing keys fall back to the
defaults below. Nothing here needs editing if params.local.yml is complete.
"""
from __future__ import annotations

import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
PARAMS_FILE = Path(os.environ.get("PARAMS_FILE", PROJECT_ROOT / "params.local.yml"))


def _read_params(path: Path) -> dict[str, str]:
    """Flat 'key: value' YAML, enough for params.local.yml (no PyYAML needed)."""
    params: dict[str, str] = {}
    if not path.is_file():
        return params
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].rstrip()
        if ":" not in line or line.startswith((" ", "\t")):
            continue
        key, value = line.split(":", 1)
        value = value.strip().strip("'\"")
        if value and value.lower() != "null":
            params[key.strip()] = value
    return params


PARAMS = _read_params(PARAMS_FILE)


def _path(key: str, default: str | Path) -> Path:
    value = Path(PARAMS.get(key, default))
    return value if value.is_absolute() else PROJECT_ROOT / value


def _env_path(variable: str, default: Path) -> Path:
    value = Path(os.environ.get(variable, default))
    return value if value.is_absolute() else PROJECT_ROOT / value


# ---------------------------------------------------------------- inputs
REFERENCE = _path("reference", "/absolute/path/to/chm13v2.0.fa")
CPG_MODEL = _path("cpg_model", "/absolute/path/to/pileup_calling_model.v1.tflite")
METADATA_PATH = _path("metadata", "assets/metadata.csv")
# sample,hap1,hap2,primary -> hifiasm *.bp.hap1/hap2.p_ctg.gfa, *.bp.p_ctg.gfa (contigs, not RagTag scaffolds)
ASSEMBLIES_PATH = Path(os.environ.get("ASSEMBLIES", PROJECT_ROOT / "assemblies.local.csv"))
# Optional sample,bam table: overrides the BAM lookup below (several BAMs per sample: ';').
BAMS_PATH = _env_path("DUPLICON_BAMS", PROJECT_ROOT / "bams.local.csv")

RESULTS = PROJECT_ROOT / "results"
# HiFi reads aligned to the SAME reference as REFERENCE, with MM/ML tags (04 only).
# Found as <sample>.aligned.bam or any *.bam whose name contains the sample ID as a
# token (SMRT Link: 08_1_A01_bc2043_001P.bam); see duplicon_analysis/bams.py.
ALIGNMENT_DIR = _env_path("DUPLICON_ALIGNMENT_DIR", RESULTS / "01_alignment")
PHASING_DIR = RESULTS / "04_phasing"                # <sample>.blocks.tsv
METHYLATION_DIR = RESULTS / "06_methylation"        # <sample>/<sample>.cpg.{combined,hap1,hap2}.bed
STRUCTURAL_EVIDENCE_PATH = RESULTS / "analysis" / "01_evidence_matrix" / "chr15_structural_evidence.tsv"

# --------------------------------------------------------------- outputs
OUT = RESULTS / "08_duplicons"
REFERENCE_DIR = OUT / "reference"                   # 01
SUNK_DIR = OUT / "sunk_counts"                      # 02
ASSEMBLY_DIR = OUT / "assembly"                     # 03, 04
BREAKPOINT_DIR = OUT / "breakpoints"                # 05
HAPLOTYPE_REPEAT_DIR = OUT / "haplotype_repeats"    # 06
REPEAT_METHYLATION_DIR = OUT / "repeat_methylation" # 07

# --------------------------------------------------------------- regions
CHROM = "chr15"
# BP1-BP5 plus unique flank in T2T-CHM13 (Hoeps et al. 2026 place 15q11.2-q13.3 at
# chr15:20.0-30.8 Mb). Must contain every HiFiCNV deletion interval.
WINDOW = os.environ.get("DUPLICON_WINDOW", "chr15:17500000-33000000")
# Diploid, SD-free normalisation region; confirm CN=2 in all 17 genomes with HiFiCNV.
CONTROL_REGION = os.environ.get("DUPLICON_CONTROL_REGION", "chr15:60000000-64000000")
# Same analysis domain, IC and CN=1 buffer as scripts/analysis/01-03.
DOMAIN_START = int(os.environ.get("DUPLICON_DOMAIN_START", 22_000_000))
DOMAIN_END = int(os.environ.get("DUPLICON_DOMAIN_END", 28_000_000))
PWS_IC_START = int(os.environ.get("DUPLICON_IC_START", 22_691_258))
PWS_IC_END = int(os.environ.get("DUPLICON_IC_END", 22_693_494))
CN1_BREAKPOINT_BUFFER = int(os.environ.get("DUPLICON_CN1_BUFFER", 75_000))

# ------------------------------------------------------------- programs
SAMTOOLS = PARAMS.get("samtools_bin", "samtools")
MINIMAP2 = PARAMS.get("minimap2_bin", "minimap2")
MERYL = PARAMS.get("meryl_bin", "meryl")
REPEATMASKER = PARAMS.get("repeatmasker_bin", "RepeatMasker")
PBCPGTOOLS = PARAMS.get("pbcpgtools_bin", "aligned_bam_to_cpg_scores")
THREADS = int(os.environ.get("DUPLICON_THREADS", 16))
MERYL_MEMORY_GB = int(os.environ.get("DUPLICON_MERYL_GB", 64))
METHYLATION_MIN_COVERAGE = 4          # as CALL_METHYLATION in the workflow

# ---------------------------------------------------------------- cohort
MECHANISMS_DIPLOID_CHR15 = ("Control", "DiGeorge", "PWS-mUPD")
MECHANISMS_BIPARENTAL = ("Control", "DiGeorge")
MECHANISMS_DELETION = ("PWS-DEL", "AS-DEL")
