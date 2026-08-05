#!/usr/bin/env python3
"""
Self-contained Figure 2 pipeline for the hifi_multiomics_pipeline layout.

This file vendors the local paper_vf helper code, the Phase 2 reciprocal
cis-architecture analysis, and the improved Figure 2 renderer so it can run
without delegating to sibling source scripts at runtime.
"""

from __future__ import annotations

import sys

# --- Vendored from scripts/paper_vf/paper_vf_q1_pipeline.py ---

"""
Q1 PWS/AS T2T methylation architecture pipeline

Purpose
-------
Transforms five research questions into analysis-ready tables and Q1-style figures:
  RQ1 Diagnostic validation at SNRPN/PWS-IC
  RQ2 Parental cis-architecture across 15q11-q13
  RQ3 Data-driven methylation boundary definition
  RQ4 Molecule-level / phased-block methylation coordination
  RQ5 SV/CNV structural context

Design principles
-----------------
- T2T-CHM13 only.
- No GRCh38/hg38 hard-coded SNRPN coordinates.
- The PWS/SNRPN imprinting-center proxy is inferred from the T2T GTF using SNRPN/SNURF/SNHG14.
- The broader domain is inferred from canonical 15q11-q13 gene annotations in the T2T GTF.
- PWS-DEL and AS-DEL are treated as reciprocal hemizygous models exposing retained maternal-like and paternal-like architectures.

Minimum expected inputs
-----------------------
1) methylation BED-like or DSS-like files with sample codes in filenames, e.g. 001P, 013A, 018C
   layers are inferred from filenames containing: hap1, hap2, combined
2) T2T GTF with gene annotations

Optional inputs
---------------
- BAM directory with MM/ML methylation tags for molecule-level barcode plots
- SV VCF directory
- CNV bedGraph/BED directory
- structural annotation BED, e.g. BP/segmental-duplication/repeat regions in T2T coordinates

Author: generated for PWS/AS long-read methylome analysis
"""


import argparse
import gzip
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle, FancyArrow


# -----------------------------------------------------------------------------
# Cohort definition
# -----------------------------------------------------------------------------

DEFAULT_GROUPS = {
    "PWS_DEL": ["001P", "002P", "005P", "006P", "007P"],
    "PWS_mUPD": ["004P"],
    "AS_DEL": ["013A", "014A", "016A"],
    "CONTROL": ["017C", "018C"],
}

DEFAULT_EXCLUDE = {"003P"}

DEFAULT_METADATA_PATH = Path("/home/rare/arlen/outputs/methylation/metadata/metadata_methylation.csv")
DEFAULT_METHYLATION_DIR = Path("/home/rare/arlen/outputs/methylation/genomes_2")
DEFAULT_VCF_DIR = Path("/mnt/diskrare/arlenb/08/hiphase_results/variants")
DEFAULT_BAM_DIR = Path("/mnt/diskrare/arlenb/08/hiphase_results/bamfiles")
DEFAULT_CNV_DIR = Path("/home/rare/arlen/outputs/Variants/Structural_variants/hifi_cnv")
DEFAULT_GTF = Path("/home/rare/arlen/reference/chm13v22.sorted.gtf")
DEFAULT_OUTDIR = Path("/home/rare/arlen/pws-as-hifi-cis-methylation/scripts/hifi_multiomics_pipeline/06_results")

EXPECTED_STATES = {
    "CONTROL": {
        "retained_parental_allele": "maternal + paternal",
        "expected_SNRPN_IC_state": "one maternal-like methylated haplotype and one paternal-like unmethylated haplotype",
    },
    "PWS_DEL": {
        "retained_parental_allele": "maternal",
        "expected_SNRPN_IC_state": "retained maternal-like methylated allele; paternal allele deleted",
    },
    "AS_DEL": {
        "retained_parental_allele": "paternal",
        "expected_SNRPN_IC_state": "retained paternal-like unmethylated allele; maternal allele deleted",
    },
    "PWS_mUPD": {
        "retained_parental_allele": "maternal + maternal",
        "expected_SNRPN_IC_state": "two maternal-like methylated haplotypes",
    },
    "PWS_UNKNOWN": {
        "retained_parental_allele": "unknown",
        "expected_SNRPN_IC_state": "unknown PWS mechanism; inspect IC methylation and CNV/UPD evidence",
    },
    "UNKNOWN": {
        "retained_parental_allele": "unknown",
        "expected_SNRPN_IC_state": "unknown",
    },
}

CANONICAL_DOMAIN_GENES = [
    "MKRN3", "MAGEL2", "NDN", "NPAP1", "SNRPN", "SNURF", "SNHG14",
    "SNORD116-1", "SNORD116", "SNORD115", "UBE3A", "UBE3A-AS1",
    "ATP10A", "GABRB3", "GABRA5", "GABRG3", "OCA2", "HERC2"
]

IC_PRIORITY_GENES = ["SNRPN", "SNURF", "SNHG14"]


# -----------------------------------------------------------------------------
# Generic utilities
# -----------------------------------------------------------------------------

def log(msg: str) -> None:
    print(f"[Q1-PWS-AS] {msg}", file=sys.stderr, flush=True)


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def open_maybe_gzip(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def detect_sample_code(path: Path) -> Optional[str]:
    """Detect sample code such as 001P, 013A, 018C from filename."""
    name = path.name
    m = re.search(r"(?<!\d)(\d{3}[PACD])(?![A-Za-z0-9])", name)
    if m:
        return m.group(1)
    # fallback: anywhere in filename
    m = re.search(r"(\d{3}[PACD])", name)
    return m.group(1) if m else None


def infer_layer(path: Path) -> str:
    name = path.name.lower()
    if re.search(r"hap[_-]?1|\.h1\.|_h1_|haplotype1", name):
        return "hap1"
    if re.search(r"hap[_-]?2|\.h2\.|_h2_|haplotype2", name):
        return "hap2"
    if "combined" in name:
        return "combined"
    return "unknown"


def sample_group(code: str, groups: Dict[str, List[str]]) -> str:
    for g, codes in groups.items():
        if code in codes:
            return g
    if code.endswith("C"):
        return "CONTROL"
    if code.endswith("A"):
        return "AS_DEL"
    if code.endswith("P"):
        return "PWS_UNKNOWN"
    if code.endswith("D"):
        return "DIGEORGE"
    return "UNKNOWN"


def parse_attrs(attr: str) -> Dict[str, str]:
    out = {}
    for part in attr.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if " " in part:
            k, v = part.split(" ", 1)
            out[k] = v.strip().strip('"')
        elif "=" in part:
            k, v = part.split("=", 1)
            out[k] = v.strip().strip('"')
    return out


def weighted_mean(values: np.ndarray, weights: Optional[np.ndarray] = None) -> float:
    values = np.asarray(values, dtype=float)
    mask = np.isfinite(values)
    if weights is None:
        return float(np.nanmean(values[mask])) if mask.any() else np.nan
    weights = np.asarray(weights, dtype=float)
    mask = mask & np.isfinite(weights) & (weights > 0)
    if not mask.any():
        return np.nan
    return float(np.average(values[mask], weights=weights[mask]))



def normalize_group_from_text(sample: str, text: str) -> str:
    """Infer one of the manuscript groups from free-text metadata plus sample code."""
    sample = str(sample)
    t = str(text).upper().replace("-", "_").replace(" ", "_")
    if any(x in t for x in ["CONTROL", "CTRL", "HEALTHY", "UNAFFECTED", "UNAFF", "NORMAL"]):
        return "CONTROL"
    if "ANGELMAN" in t or re.search(r"(^|_)AS(_|$)", t):
        return "AS_DEL" if "DEL" in t or "DELETION" in t or sample.endswith("A") else "AS_DEL"
    if "MUPD" in t or "MATERNAL_UPD" in t or "UPD" in t:
        return "PWS_mUPD" if "PWS" in t or sample.endswith("P") else "UNKNOWN"
    if "PRADER" in t or re.search(r"(^|_)PWS(_|$)", t):
        if "UPD" in t:
            return "PWS_mUPD"
        if "DEL" in t or "DELETION" in t:
            return "PWS_DEL"
        return "PWS_DEL" if sample != "004P" else "PWS_mUPD"
    # Hard fallback from sample code convention used in this project.
    if sample == "004P":
        return "PWS_mUPD"
    if sample.endswith("P"):
        return "PWS_DEL"
    if sample.endswith("A"):
        return "AS_DEL"
    if sample.endswith("C"):
        return "CONTROL"
    return "UNKNOWN"


def load_metadata_table(metadata_path: Optional[Path], outdir: Optional[Path] = None) -> pd.DataFrame:
    """Load metadata_methylation.csv if present and use it to update DEFAULT_GROUPS in-place.

    The file can have flexible column names. We search for a sample-like column and infer
    group from group/disease/class/molecular text columns, falling back to sample code.
    """
    if metadata_path is None or not Path(metadata_path).exists():
        log(f"Metadata not found or not supplied: {metadata_path}. Falling back to filename/sample-code groups.")
        rows = []
        for g, codes in DEFAULT_GROUPS.items():
            for sample in codes:
                rows.append({"sample": sample, "group": g})
        meta = pd.DataFrame(rows)
    else:
        meta_raw = pd.read_csv(metadata_path, sep=None, engine="python")
        # sample column candidates
        sample_col = None
        candidates = ["sample", "sample_id", "id", "code", "Sample", "Sample_ID", "sampleID", "participant", "individual"]
        lower_map = {c.lower(): c for c in meta_raw.columns}
        for c in candidates:
            if c in meta_raw.columns:
                sample_col = c
                break
            if c.lower() in lower_map:
                sample_col = lower_map[c.lower()]
                break
        if sample_col is None:
            # Try any column containing project-style codes.
            for c in meta_raw.columns:
                vals = meta_raw[c].astype(str).head(100).tolist()
                if any(re.search(r"\d{3}[PACD]", v) for v in vals):
                    sample_col = c
                    break
        if sample_col is None:
            raise ValueError(f"Could not identify sample column in metadata: {metadata_path}")

        preferred_text_cols = [
            c for c in meta_raw.columns
            if any(k in c.lower() for k in ["group", "diagn", "disease", "class", "molecular", "status", "phenotype", "type"])
        ]
        if not preferred_text_cols:
            preferred_text_cols = list(meta_raw.columns)

        rows = []
        for _, r in meta_raw.iterrows():
            raw_sample = str(r[sample_col])
            code = detect_sample_code(Path(raw_sample)) or raw_sample.strip()
            if not code or code.lower() in {"nan", "none"}:
                continue
            text = " ".join(str(r[c]) for c in preferred_text_cols if c in meta_raw.columns)
            group = normalize_group_from_text(code, text)
            state = EXPECTED_STATES.get(group, EXPECTED_STATES["UNKNOWN"])
            rows.append({
                "sample": code,
                "group": group,
                "metadata_sample_value": raw_sample,
                "metadata_group_text": text,
                "expected_retained_parental_allele": state["retained_parental_allele"],
                "expected_SNRPN_IC_state": state["expected_SNRPN_IC_state"],
            })
        meta = pd.DataFrame(rows).drop_duplicates(subset=["sample"], keep="first")

    # Update runtime group dictionary in-place so all existing calls use metadata-aware labels.
    grouped = meta.groupby("group")["sample"].apply(list).to_dict() if not meta.empty else {}
    for g, codes in grouped.items():
        if g in {"CONTROL", "PWS_DEL", "PWS_mUPD", "AS_DEL"}:
            DEFAULT_GROUPS[g] = sorted(set(codes))

    if outdir is not None:
        tbldir = Path(outdir) / "tables"
        tbldir.mkdir(parents=True, exist_ok=True)
        meta.to_csv(tbldir / "metadata_loaded_and_expected_states.tsv", sep="\t", index=False)
        exp = []
        for group, state in EXPECTED_STATES.items():
            exp.append({"group": group, **state})
        pd.DataFrame(exp).to_csv(tbldir / "expected_parental_methylation_states.tsv", sep="\t", index=False)
    return meta

# -----------------------------------------------------------------------------
# Input discovery
# -----------------------------------------------------------------------------

def discover_methylation_files(meth_dir: Path, exclude: set[str]) -> pd.DataFrame:
    patterns = ["*.bed", "*.bed.gz", "*.bedgraph", "*.bedGraph", "*.txt", "*.txt.gz", "*.tsv", "*.tsv.gz"]
    rows = []
    for pat in patterns:
        for p in meth_dir.rglob(pat):
            code = detect_sample_code(p)
            if not code or code in exclude:
                continue
            layer = infer_layer(p)
            rows.append({"sample": code, "layer": layer, "path": str(p)})
    df = pd.DataFrame(rows).drop_duplicates()
    if df.empty:
        raise FileNotFoundError(f"No methylation files with sample codes found in: {meth_dir}")
    return df.sort_values(["sample", "layer", "path"]).reset_index(drop=True)


def discover_files_by_sample(directory: Optional[Path], suffixes: Tuple[str, ...], exclude: set[str]) -> pd.DataFrame:
    if directory is None or not directory.exists():
        return pd.DataFrame(columns=["sample", "path"])
    rows = []
    for p in directory.rglob("*"):
        if not p.is_file():
            continue
        if not any(str(p).endswith(s) for s in suffixes):
            continue
        code = detect_sample_code(p)
        if code and code not in exclude:
            rows.append({"sample": code, "path": str(p)})
    return pd.DataFrame(rows).drop_duplicates().sort_values(["sample", "path"]).reset_index(drop=True)


# -----------------------------------------------------------------------------
# GTF parsing and T2T locus inference
# -----------------------------------------------------------------------------

def load_gtf_genes(gtf_path: Path, chrom: str) -> pd.DataFrame:
    rows = []
    with open_maybe_gzip(gtf_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, source, feature, start, end, score, strand, frame, attrs = fields
            if seqid != chrom or feature != "gene":
                continue
            a = parse_attrs(attrs)
            gene = a.get("gene_name") or a.get("gene") or a.get("Name") or a.get("gene_id")
            gene_id = a.get("gene_id", gene)
            if gene is None:
                continue
            rows.append({
                "chrom": seqid,
                "start": int(start) - 1,
                "end": int(end),
                "strand": strand,
                "gene": gene,
                "gene_id": gene_id,
                "source": source,
            })
    genes = pd.DataFrame(rows)
    if genes.empty:
        raise ValueError(f"No gene records found for {chrom} in {gtf_path}")
    genes = genes.drop_duplicates(subset=["chrom", "start", "end", "gene", "strand"])
    return genes.sort_values(["chrom", "start", "end"]).reset_index(drop=True)


def infer_domain_from_gtf(
    genes: pd.DataFrame,
    domain_genes: List[str],
    padding: int,
    chrom: str,
) -> Tuple[int, int, pd.DataFrame]:
    gene_upper = {g.upper() for g in domain_genes}
    g = genes[genes["gene"].str.upper().isin(gene_upper)].copy()
    if g.empty:
        raise ValueError(
            "Could not infer 15q11-q13 domain from GTF because none of the canonical genes "
            f"were found on {chrom}. Check gene symbols in your T2T GTF."
        )
    start = max(0, int(g["start"].min()) - padding)
    end = int(g["end"].max()) + padding
    return start, end, g.sort_values("start").reset_index(drop=True)


def infer_ic_region_from_gtf(
    genes: pd.DataFrame,
    priority_genes: List[str],
    flank: int,
) -> Tuple[str, int, int, str, pd.Series]:
    for gene in priority_genes:
        hit = genes[genes["gene"].str.upper() == gene.upper()].copy()
        if not hit.empty:
            row = hit.sort_values("start").iloc[0]
            tss = int(row["start"]) if row["strand"] == "+" else int(row["end"])
            start = max(0, tss - flank)
            end = tss + flank
            return row["chrom"], start, end, gene, row
    raise ValueError(
        "Could not infer PWS/SNRPN imprinting-center proxy from GTF. "
        f"None of these genes were found: {priority_genes}"
    )


# -----------------------------------------------------------------------------
# Methylation loaders
# -----------------------------------------------------------------------------

def _normalize_methyl_fraction(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    med = x.dropna().median() if x.notna().any() else np.nan
    # Many tools store methylation as 0-100 or percent.
    if np.isfinite(med) and med > 1.5:
        x = x / 100.0
    return x.clip(lower=0, upper=1)


def _parse_column_arg(value: Optional[Union[int, str]], allow_auto: bool = True):
    """Parse 1-based CLI column arguments.

    Returns an integer 1-based column, "auto", or None.
    """
    if value is None:
        return None
    if isinstance(value, str):
        v = value.strip().lower()
        if v in {"", "none", "null", "na"}:
            return None
        if allow_auto and v == "auto":
            return "auto"
        return int(v)
    return int(value)


def _is_numeric_candidate(vals: pd.Series, min_nonmissing: int = 5) -> bool:
    x = pd.to_numeric(vals, errors="coerce").dropna()
    return len(x) >= min_nonmissing


def _auto_bed_columns(sub: pd.DataFrame, bed_meth_col: Optional[Union[int, str]], bed_cov_col: Optional[Union[int, str]]) -> Tuple[int, Optional[int]]:
    """Infer 0-based methylation and coverage columns for BED-like methylation files.

    Supports common layouts:
    - simple BED: chrom, start, end, methylation_fraction_or_percent
    - modkit/bedMethyl: methylation fraction/percent in column 11, coverage in column 10
    - custom BED with score/coverage columns
    """
    meth_arg = _parse_column_arg(bed_meth_col, allow_auto=True)
    cov_arg = _parse_column_arg(bed_cov_col, allow_auto=True)

    if meth_arg != "auto":
        meth_idx = int(meth_arg) - 1
    else:
        preferred_zero_based = [10, 3, 4, 5, 11, 9]  # 11th modkit; 4th simple BED; then fallbacks.
        candidates = []
        sample = sub.head(min(len(sub), 5000))
        for idx in preferred_zero_based:
            if idx >= sample.shape[1]:
                continue
            vals = pd.to_numeric(sample.iloc[:, idx], errors="coerce").dropna()
            if len(vals) < 5:
                continue
            q01, q99 = np.nanpercentile(vals, [1, 99])
            med = float(np.nanmedian(vals))
            # Methylation columns should usually be 0-1 or 0-100. Allow tiny numerical noise.
            if q01 < -1e-6 or q99 > 100 + 1e-6:
                continue
            score = 0
            if idx == 10:  # canonical modkit/bedMethyl fraction_modified
                score += 20
            if idx == 3:   # simple BED methylation column
                score += 12
            if 0 <= med <= 1:
                score += 8
            elif 1 < med <= 100:
                score += 6
            # Prefer columns with many unique values over binary-only columns when possible.
            score += min(vals.nunique(), 50) / 100.0
            candidates.append((score, idx))
        if not candidates:
            raise ValueError(
                "Could not infer methylation column automatically. Provide --bed-meth-col, "
                "for example --bed-meth-col 11 for modkit bedMethyl or --bed-meth-col 4 for simple BED."
            )
        meth_idx = sorted(candidates, reverse=True)[0][1]

    if cov_arg is None:
        cov_idx = None
    elif cov_arg != "auto":
        cov_idx = int(cov_arg) - 1
    else:
        sample = sub.head(min(len(sub), 5000))
        cov_idx = None
        # modkit/bedMethyl valid_coverage is column 10 (0-based 9).
        for idx in [9, 4, 5, 10]:
            if idx >= sample.shape[1] or idx == meth_idx:
                continue
            vals = pd.to_numeric(sample.iloc[:, idx], errors="coerce").dropna()
            if len(vals) >= 5 and np.nanmedian(vals) > 0 and np.nanpercentile(vals, 99) > 1:
                cov_idx = idx
                break
    return meth_idx, cov_idx


def read_methyl_region(
    path: Path,
    chrom: str,
    start: int,
    end: int,
    fmt: str = "bed",
    bed_meth_col: Union[int, str] = "auto",
    bed_cov_col: Optional[Union[int, str]] = "auto",
    chunksize: int = 1_000_000,
) -> pd.DataFrame:
    """Read methylation records overlapping chrom:start-end.

    Supported formats:
    - bed: no header; columns are chrom, start, end, methylation by default.
           bed_meth_col and bed_cov_col are 1-based column numbers.
    - dss: no header; columns are chr, pos, N, X.
    """
    dfs = []
    fmt = fmt.lower()
    if fmt not in {"bed", "dss"}:
        raise ValueError("fmt must be 'bed' or 'dss'")

    try:
        reader = pd.read_csv(
            path,
            sep="\t",
            header=None,
            comment="#",
            chunksize=chunksize,
            low_memory=False,
        )
        for chunk in reader:
            if chunk.empty or chunk.shape[1] < 4:
                continue
            c0 = chunk.iloc[:, 0].astype(str)
            sub = chunk[c0 == chrom].copy()
            if sub.empty:
                continue

            if fmt == "dss":
                pos = pd.to_numeric(sub.iloc[:, 1], errors="coerce").astype("Int64")
                n = pd.to_numeric(sub.iloc[:, 2], errors="coerce")
                x = pd.to_numeric(sub.iloc[:, 3], errors="coerce")
                tmp = pd.DataFrame({
                    "chrom": sub.iloc[:, 0].astype(str),
                    "start": pos.astype(float) - 1,
                    "end": pos.astype(float),
                    "meth": x / n.replace(0, np.nan),
                    "coverage": n,
                })
            else:
                meth_idx, cov_idx = _auto_bed_columns(sub, bed_meth_col, bed_cov_col)
                if meth_idx >= sub.shape[1]:
                    raise IndexError(f"bed_meth_col={bed_meth_col} exceeds number of columns in {path}")
                tmp = pd.DataFrame({
                    "chrom": sub.iloc[:, 0].astype(str),
                    "start": pd.to_numeric(sub.iloc[:, 1], errors="coerce"),
                    "end": pd.to_numeric(sub.iloc[:, 2], errors="coerce"),
                    "meth": _normalize_methyl_fraction(sub.iloc[:, meth_idx]),
                })
                if cov_idx is not None and cov_idx < sub.shape[1]:
                    tmp["coverage"] = pd.to_numeric(sub.iloc[:, cov_idx], errors="coerce")
                else:
                    tmp["coverage"] = 1.0

            tmp = tmp.dropna(subset=["start", "end", "meth"])
            tmp["start"] = tmp["start"].astype(int)
            tmp["end"] = tmp["end"].astype(int)
            tmp = tmp[(tmp["end"] > start) & (tmp["start"] < end)]
            if not tmp.empty:
                dfs.append(tmp)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=["chrom", "start", "end", "meth", "coverage"])

    if not dfs:
        return pd.DataFrame(columns=["chrom", "start", "end", "meth", "coverage"])
    out = pd.concat(dfs, ignore_index=True)
    out["mid"] = ((out["start"] + out["end"]) / 2).astype(int)
    return out


def summarize_region_for_files(
    file_table: pd.DataFrame,
    chrom: str,
    start: int,
    end: int,
    args,
) -> pd.DataFrame:
    rows = []
    for _, r in file_table.iterrows():
        path = Path(r["path"])
        fmt = args.methyl_format
        if fmt == "auto":
            fmt = "dss" if "dss" in path.name.lower() else "bed"
        m = read_methyl_region(path, chrom, start, end, fmt, args.bed_meth_col, args.bed_cov_col, args.chunksize)
        rows.append({
            "sample": r["sample"],
            "group": sample_group(r["sample"], DEFAULT_GROUPS),
            "layer": r["layer"],
            "path": str(path),
            "chrom": chrom,
            "start": start,
            "end": end,
            "n_cpg": int(len(m)),
            "total_coverage": float(m["coverage"].sum()) if not m.empty else 0.0,
            "mean_methylation": weighted_mean(m["meth"].values, m["coverage"].values) if not m.empty else np.nan,
        })
    return pd.DataFrame(rows)


def load_binned_profiles(
    file_table: pd.DataFrame,
    chrom: str,
    start: int,
    end: int,
    bin_size: int,
    args,
) -> pd.DataFrame:
    rows = []
    for i, r in file_table.iterrows():
        path = Path(r["path"])
        log(f"Binning methylation {i+1}/{len(file_table)}: {path.name}")
        fmt = args.methyl_format
        if fmt == "auto":
            fmt = "dss" if "dss" in path.name.lower() else "bed"
        m = read_methyl_region(path, chrom, start, end, fmt, args.bed_meth_col, args.bed_cov_col, args.chunksize)
        if m.empty:
            continue
        m["bin_start"] = ((m["mid"] - start) // bin_size) * bin_size + start
        m["bin_end"] = m["bin_start"] + bin_size
        grouped = []
        for (bs, be), z in m.groupby(["bin_start", "bin_end"], sort=True):
            grouped.append({
                "sample": r["sample"],
                "group": sample_group(r["sample"], DEFAULT_GROUPS),
                "layer": r["layer"],
                "bin_start": int(bs),
                "bin_end": int(be),
                "bin_mid": int((bs + be) / 2),
                "n_cpg": int(len(z)),
                "total_coverage": float(z["coverage"].sum()),
                "mean_methylation": weighted_mean(z["meth"].values, z["coverage"].values),
            })
        if grouped:
            rows.extend(grouped)
    return pd.DataFrame(rows)


def load_single_cpg_region_for_plot(
    file_table: pd.DataFrame,
    chrom: str,
    start: int,
    end: int,
    args,
    max_points_per_file: int = 20000,
) -> pd.DataFrame:
    rows = []
    for _, r in file_table.iterrows():
        if r["layer"] == "unknown":
            continue
        path = Path(r["path"])
        fmt = args.methyl_format
        if fmt == "auto":
            fmt = "dss" if "dss" in path.name.lower() else "bed"
        m = read_methyl_region(path, chrom, start, end, fmt, args.bed_meth_col, args.bed_cov_col, args.chunksize)
        if m.empty:
            continue
        if len(m) > max_points_per_file:
            m = m.sample(max_points_per_file, random_state=1).sort_values("mid")
        m["sample"] = r["sample"]
        m["group"] = sample_group(r["sample"], DEFAULT_GROUPS)
        m["layer"] = r["layer"]
        rows.append(m[["sample", "group", "layer", "chrom", "start", "end", "mid", "meth", "coverage"]])
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


# -----------------------------------------------------------------------------
# Allele assignment
# -----------------------------------------------------------------------------

def assign_allele_labels(ic_summary: pd.DataFrame, min_cpgs: int = 3) -> pd.DataFrame:
    """Classify hap1/hap2 as maternal-like/paternal-like using IC methylation.

    For controls: higher methylated haplotype is maternal-like, lower is paternal-like.
    For PWS-DEL: retained methylated haplotype is maternal-like.
    For AS-DEL: retained unmethylated haplotype is paternal-like.
    For PWS-mUPD: both informative haplotypes are expected maternal-like.
    """
    rows = []
    x = ic_summary[ic_summary["layer"].isin(["hap1", "hap2"])].copy()
    for sample, z in x.groupby("sample"):
        z = z[z["n_cpg"] >= min_cpgs].sort_values("mean_methylation")
        group = sample_group(sample, DEFAULT_GROUPS)
        if z.empty:
            continue
        if group == "CONTROL":
            if len(z) >= 2:
                low = z.iloc[0]
                high = z.iloc[-1]
                rows.append({"sample": sample, "layer": high["layer"], "allele_label": "control_maternal_like", "ic_mean": high["mean_methylation"]})
                rows.append({"sample": sample, "layer": low["layer"], "allele_label": "control_paternal_like", "ic_mean": low["mean_methylation"]})
            else:
                only = z.iloc[0]
                lab = "control_maternal_like" if only["mean_methylation"] >= 0.5 else "control_paternal_like"
                rows.append({"sample": sample, "layer": only["layer"], "allele_label": lab, "ic_mean": only["mean_methylation"]})
        elif group == "PWS_DEL":
            for _, rr in z.iterrows():
                lab = "pwsdel_retained_maternal_like" if rr["mean_methylation"] >= 0.5 else "pwsdel_low_or_deleted_like"
                rows.append({"sample": sample, "layer": rr["layer"], "allele_label": lab, "ic_mean": rr["mean_methylation"]})
        elif group == "AS_DEL":
            for _, rr in z.iterrows():
                lab = "asdel_retained_paternal_like" if rr["mean_methylation"] < 0.5 else "asdel_high_or_noise_like"
                rows.append({"sample": sample, "layer": rr["layer"], "allele_label": lab, "ic_mean": rr["mean_methylation"]})
        elif group == "PWS_mUPD":
            for _, rr in z.iterrows():
                lab = "upd_maternal_like" if rr["mean_methylation"] >= 0.5 else "upd_unexpected_low"
                rows.append({"sample": sample, "layer": rr["layer"], "allele_label": lab, "ic_mean": rr["mean_methylation"]})
        else:
            for _, rr in z.iterrows():
                lab = "unknown_high" if rr["mean_methylation"] >= 0.5 else "unknown_low"
                rows.append({"sample": sample, "layer": rr["layer"], "allele_label": lab, "ic_mean": rr["mean_methylation"]})
    return pd.DataFrame(rows)


def add_allele_labels(profiles: pd.DataFrame, allele_table: pd.DataFrame) -> pd.DataFrame:
    if profiles.empty or allele_table.empty:
        profiles["allele_label"] = profiles.get("layer", "unknown")
        return profiles
    out = profiles.merge(allele_table[["sample", "layer", "allele_label"]], on=["sample", "layer"], how="left")
    out["allele_label"] = out["allele_label"].fillna(out["group"] + "_" + out["layer"].astype(str))
    return out


# -----------------------------------------------------------------------------
# Figure helpers
# -----------------------------------------------------------------------------

def savefig(fig: plt.Figure, out: Path) -> None:
    fig.savefig(out, dpi=350, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_gene_track(ax, genes: pd.DataFrame, start: int, end: int, title: str = ""):
    ax.set_xlim(start, end)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.hlines(0.5, start, end, color="black", lw=0.8)
    if genes.empty:
        ax.text(0.5, 0.5, "No genes", transform=ax.transAxes, ha="center", va="center")
        return
    levels = np.linspace(0.25, 0.8, 4)
    for i, (_, g) in enumerate(genes.sort_values("start").iterrows()):
        y = levels[i % len(levels)]
        x0, x1 = int(g["start"]), int(g["end"])
        width = max(1, x1 - x0)
        ax.add_patch(Rectangle((x0, y - 0.035), width, 0.07, alpha=0.35, ec="black", lw=0.3))
        ax.text((x0 + x1) / 2, y + 0.07, g["gene"], ha="center", va="bottom", fontsize=6, rotation=45)
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")


def plot_cohort_table(ax, file_table: pd.DataFrame):
    ax.axis("off")
    samples = sorted(file_table["sample"].unique())
    groups = pd.DataFrame({"sample": samples})
    groups["group"] = groups["sample"].apply(lambda s: sample_group(s, DEFAULT_GROUPS))
    counts = groups.groupby("group")["sample"].apply(lambda x: ", ".join(x)).reset_index()
    counts["n"] = counts["sample"].apply(lambda s: len(s.split(", ")) if s else 0)
    table_data = counts[["group", "n", "sample"]].values.tolist()
    tbl = ax.table(
        cellText=table_data,
        colLabels=["Group", "n", "Samples"],
        cellLoc="left",
        loc="center",
        colWidths=[0.25, 0.12, 0.63],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1, 1.4)
    ax.set_title("B. Cohort composition", loc="left", fontsize=10, fontweight="bold")


def plot_ic_heatmap(ax, ic_summary: pd.DataFrame):
    x = ic_summary[ic_summary["layer"].isin(["hap1", "hap2", "combined"])].copy()
    if x.empty:
        ax.text(0.5, 0.5, "No IC methylation data", ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        return
    x["row"] = x["sample"] + "\n" + x["group"]
    order = sorted(x["row"].unique(), key=lambda r: (r.split("\n")[-1], r))
    cols = [c for c in ["hap1", "hap2", "combined"] if c in x["layer"].unique()]
    mat = x.pivot_table(index="row", columns="layer", values="mean_methylation", aggfunc="mean").reindex(order)[cols]
    im = ax.imshow(mat.values, aspect="auto", vmin=0, vmax=1, cmap="viridis")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=45, ha="right")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=7)
    ax.set_title("C. SNRPN/PWS-IC methylation", loc="left", fontsize=10, fontweight="bold")
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("Methylation", fontsize=8)


def plot_ic_profiles(ax, single_cpg: pd.DataFrame, ic_start: int, ic_end: int):
    if single_cpg.empty:
        ax.text(0.5, 0.5, "No CpGs in IC window", transform=ax.transAxes, ha="center", va="center")
        return
    # keep haplotype tracks only for clarity
    x = single_cpg[single_cpg["layer"].isin(["hap1", "hap2"])].copy()
    if x.empty:
        x = single_cpg.copy()
    groups = ["CONTROL", "PWS_DEL", "AS_DEL", "PWS_mUPD"]
    for g in groups:
        z = x[x["group"] == g]
        if z.empty:
            continue
        # binned smooth inside panel
        z = z.copy()
        nb = 120
        bins = np.linspace(ic_start, ic_end, nb + 1)
        z["bin"] = pd.cut(z["mid"], bins=bins, labels=False, include_lowest=True)
        sm = z.groupby("bin").agg(mid=("mid", "mean"), meth=("meth", "mean"), n=("meth", "size")).dropna()
        if not sm.empty:
            ax.plot(sm["mid"], sm["meth"], lw=1.5, label=g)
    ax.axvspan(ic_start, ic_end, color="grey", alpha=0.08)
    ax.set_ylim(-0.05, 1.05)
    ax.set_ylabel("Methylation")
    ax.set_xlabel("T2T coordinate")
    ax.set_title("D. Per-CpG profile near inferred IC", loc="left", fontsize=10, fontweight="bold")
    ax.legend(fontsize=7, frameon=False, ncol=2)


# -----------------------------------------------------------------------------
# Boundary detection
# -----------------------------------------------------------------------------

def mean_profile_by_label(profiles: pd.DataFrame, min_cpg: int = 3) -> pd.DataFrame:
    x = profiles[profiles["n_cpg"] >= min_cpg].copy()
    if x.empty:
        return pd.DataFrame()
    rows = []
    for (lab, bs, be, bm), z in x.groupby(["allele_label", "bin_start", "bin_end", "bin_mid"]):
        rows.append({
            "allele_label": lab,
            "bin_start": bs,
            "bin_end": be,
            "bin_mid": bm,
            "n_samples": z["sample"].nunique(),
            "mean_methylation": weighted_mean(z["mean_methylation"].values, z["total_coverage"].values),
            "sd_methylation": float(np.nanstd(z["mean_methylation"].values)),
        })
    return pd.DataFrame(rows)


def pick_architecture_labels(mean_prof: pd.DataFrame) -> Tuple[Optional[str], Optional[str]]:
    labels = set(mean_prof["allele_label"].unique())
    maternal_candidates = [
        "pwsdel_retained_maternal_like",
        "control_maternal_like",
        "upd_maternal_like",
    ]
    paternal_candidates = [
        "asdel_retained_paternal_like",
        "control_paternal_like",
    ]
    mat = next((x for x in maternal_candidates if x in labels), None)
    pat = next((x for x in paternal_candidates if x in labels), None)
    return mat, pat


def calculate_delta_profile(mean_prof: pd.DataFrame, maternal_label: str, paternal_label: str) -> pd.DataFrame:
    a = mean_prof[mean_prof["allele_label"] == maternal_label][["bin_start", "bin_end", "bin_mid", "mean_methylation", "n_samples"]]
    b = mean_prof[mean_prof["allele_label"] == paternal_label][["bin_start", "bin_end", "bin_mid", "mean_methylation", "n_samples"]]
    m = a.merge(b, on=["bin_start", "bin_end", "bin_mid"], suffixes=("_maternal", "_paternal"))
    m["delta_methylation"] = m["mean_methylation_maternal"] - m["mean_methylation_paternal"]
    m = m.sort_values("bin_mid").reset_index(drop=True)
    # rolling smooth
    m["delta_smooth"] = m["delta_methylation"].rolling(7, center=True, min_periods=1).median()
    return m


def detect_boundaries(delta: pd.DataFrame, threshold: float = 0.20, min_consecutive_bins: int = 3) -> Dict[str, float]:
    if delta.empty:
        return {"left_boundary": np.nan, "right_boundary": np.nan, "n_bins_domain": 0}
    flag = delta["delta_smooth"].abs() >= threshold
    # Find longest consecutive TRUE run.
    best_start = best_end = None
    cur_start = None
    best_len = 0
    for i, val in enumerate(flag.values):
        if val and cur_start is None:
            cur_start = i
        if (not val or i == len(flag) - 1) and cur_start is not None:
            cur_end = i if val and i == len(flag) - 1 else i - 1
            run_len = cur_end - cur_start + 1
            if run_len > best_len:
                best_len = run_len
                best_start, best_end = cur_start, cur_end
            cur_start = None
    if best_len < min_consecutive_bins or best_start is None:
        return {"left_boundary": np.nan, "right_boundary": np.nan, "n_bins_domain": int(best_len)}
    left = int(delta.iloc[best_start]["bin_start"])
    right = int(delta.iloc[best_end]["bin_end"])
    return {"left_boundary": left, "right_boundary": right, "n_bins_domain": int(best_len)}


# -----------------------------------------------------------------------------
# SV/CNV loaders
# -----------------------------------------------------------------------------

def parse_info_field(info: str) -> Dict[str, str]:
    out = {}
    for part in info.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
        elif part:
            out[part] = "TRUE"
    return out


def read_sv_vcf(path: Path, sample: str) -> pd.DataFrame:
    rows = []
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 8:
                    continue
                chrom, pos, vid, ref, alt, qual, flt, info = f[:8]
                pos = int(pos)
                inf = parse_info_field(info)
                svtype = inf.get("SVTYPE")
                if not svtype:
                    if "<DEL>" in alt:
                        svtype = "DEL"
                    elif "<DUP>" in alt:
                        svtype = "DUP"
                    elif "<INV>" in alt:
                        svtype = "INV"
                    elif "<INS>" in alt:
                        svtype = "INS"
                    else:
                        svtype = "UNK"
                end = inf.get("END")
                svlen = inf.get("SVLEN")
                try:
                    end = int(end) if end is not None else None
                except Exception:
                    end = None
                try:
                    svlen_i = abs(int(str(svlen).split(",")[0])) if svlen is not None else None
                except Exception:
                    svlen_i = None
                if end is None:
                    end = pos + (svlen_i if svlen_i else 1)
                rows.append({
                    "sample": sample,
                    "group": sample_group(sample, DEFAULT_GROUPS),
                    "chrom": chrom,
                    "start": pos - 1,
                    "end": max(end, pos),
                    "svtype": svtype,
                    "svlen": svlen_i if svlen_i is not None else max(end - pos, 1),
                    "source": str(path),
                })
    except Exception as e:
        log(f"WARNING: failed to parse VCF {path}: {e}")
    return pd.DataFrame(rows)


def load_all_svs(vcf_table: pd.DataFrame) -> pd.DataFrame:
    all_rows = []
    for _, r in vcf_table.iterrows():
        df = read_sv_vcf(Path(r["path"]), r["sample"])
        if not df.empty:
            all_rows.append(df)
    return pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()


def read_bed_like(path: Path, sample: Optional[str] = None) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep="\t", header=None, comment="#")
    except Exception:
        return pd.DataFrame()
    if df.shape[1] < 3:
        return pd.DataFrame()
    out = pd.DataFrame({
        "chrom": df.iloc[:, 0].astype(str),
        "start": pd.to_numeric(df.iloc[:, 1], errors="coerce"),
        "end": pd.to_numeric(df.iloc[:, 2], errors="coerce"),
    }).dropna()
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)
    if df.shape[1] >= 4:
        out["value"] = df.iloc[:, 3]
    else:
        out["value"] = 1
    if sample:
        out["sample"] = sample
        out["group"] = sample_group(sample, DEFAULT_GROUPS)
    out["source"] = str(path)
    return out


def interval_overlap(a: pd.DataFrame, b: pd.DataFrame) -> pd.DataFrame:
    if a.empty or b.empty:
        return pd.DataFrame()
    rows = []
    for chrom, aa in a.groupby("chrom"):
        bb = b[b["chrom"] == chrom].sort_values("start")
        if bb.empty:
            continue
        for _, r in aa.iterrows():
            hit = bb[(bb["end"] > r["start"]) & (bb["start"] < r["end"])]
            for _, h in hit.iterrows():
                rows.append({
                    "chrom": chrom,
                    "a_start": r["start"], "a_end": r["end"],
                    "b_start": h["start"], "b_end": h["end"],
                    "overlap_bp": int(min(r["end"], h["end"]) - max(r["start"], h["start"])),
                    **{f"a_{k}": v for k, v in r.items() if k not in ["chrom", "start", "end"]},
                    **{f"b_{k}": v for k, v in h.items() if k not in ["chrom", "start", "end"]},
                })
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# RQ4 molecule-level barcode using BAM MM/ML tags
# -----------------------------------------------------------------------------

def extract_molecule_methylation_from_bam(
    bam_path: Path,
    sample: str,
    chrom: str,
    start: int,
    end: int,
    min_sites_per_read: int = 3,
    max_reads: int = 1500,
) -> pd.DataFrame:
    try:
        import pysam
    except ImportError:
        log("pysam is not installed; skipping molecule-level methylation extraction.")
        return pd.DataFrame()

    rows = []
    try:
        bam = pysam.AlignmentFile(str(bam_path), "rb")
    except Exception as e:
        log(f"Could not open BAM {bam_path}: {e}")
        return pd.DataFrame()

    n_reads = 0
    try:
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            mods = None
            try:
                mods = read.modified_bases
            except Exception:
                mods = None
            if not mods:
                continue
            # query-position -> reference-position map
            q2r = {}
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if qpos is not None and rpos is not None:
                    q2r[qpos] = rpos
            hp = "NA"
            try:
                hp = str(read.get_tag("HP"))
            except Exception:
                pass
            read_rows = []
            for key, vals in mods.items():
                # key usually resembles ('C', strand, 'm') for 5mC
                canonical = key[0] if isinstance(key, tuple) and len(key) > 0 else None
                mod_code = key[2] if isinstance(key, tuple) and len(key) > 2 else str(key)
                if canonical != "C":
                    continue
                if mod_code not in {"m", "C+m", "5mC", "mC"} and "m" not in str(mod_code):
                    continue
                for qpos, qual in vals:
                    rpos = q2r.get(qpos)
                    if rpos is None or rpos < start or rpos >= end:
                        continue
                    prob = float(qual) / 255.0 if qual is not None else np.nan
                    read_rows.append({
                        "sample": sample,
                        "group": sample_group(sample, DEFAULT_GROUPS),
                        "read_id": read.query_name,
                        "haplotype_tag": hp,
                        "chrom": chrom,
                        "pos": int(rpos),
                        "meth_prob": prob,
                        "meth_state": int(prob >= 0.5) if np.isfinite(prob) else np.nan,
                    })
            if len(read_rows) >= min_sites_per_read:
                rows.extend(read_rows)
                n_reads += 1
            if n_reads >= max_reads:
                break
    except ValueError as e:
        log(f"Could not fetch {chrom}:{start}-{end} from {bam_path.name}; check contig names and BAM index: {e}")
    finally:
        bam.close()
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Main figures
# -----------------------------------------------------------------------------

def figure1(outdir: Path, file_table: pd.DataFrame, domain_genes: pd.DataFrame, ic_summary: pd.DataFrame,
            ic_single: pd.DataFrame, domain_start: int, domain_end: int, ic_start: int, ic_end: int):
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[0.9, 1.4], width_ratios=[1.2, 1])
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0])
    axD = fig.add_subplot(gs[1, 1])
    plot_gene_track(axA, domain_genes, domain_start, domain_end, "A. T2T-inferred 15q11-q13 gene context")
    plot_cohort_table(axB, file_table)
    plot_ic_heatmap(axC, ic_summary)
    plot_ic_profiles(axD, ic_single, ic_start, ic_end)
    fig.suptitle("Figure 1. Cohort structure and diagnostic methylation validation", y=0.99, fontsize=14, fontweight="bold")
    savefig(fig, outdir / "Figure1_diagnostic_validation.png")


def figure2(outdir: Path, mean_prof: pd.DataFrame, domain_genes: pd.DataFrame, domain_start: int, domain_end: int):
    fig = plt.figure(figsize=(15, 8))
    gs = GridSpec(2, 1, height_ratios=[0.9, 3.0], figure=fig)
    ax0 = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[1, 0], sharex=ax0)
    plot_gene_track(ax0, domain_genes, domain_start, domain_end, "A. Gene context")

    preferred_order = [
        "control_maternal_like", "pwsdel_retained_maternal_like", "upd_maternal_like",
        "control_paternal_like", "asdel_retained_paternal_like",
    ]
    labels = [x for x in preferred_order if x in set(mean_prof["allele_label"])]
    labels += [x for x in sorted(mean_prof["allele_label"].unique()) if x not in labels]

    for lab in labels:
        z = mean_prof[mean_prof["allele_label"] == lab].sort_values("bin_mid")
        if z.empty:
            continue
        ax.plot(z["bin_mid"], z["mean_methylation"], lw=1.8, label=lab.replace("_", " "))
        ax.fill_between(
            z["bin_mid"],
            (z["mean_methylation"] - z["sd_methylation"]).clip(0, 1),
            (z["mean_methylation"] + z["sd_methylation"]).clip(0, 1),
            alpha=0.10,
        )
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(domain_start, domain_end)
    ax.set_ylabel("Mean methylation")
    ax.set_xlabel("T2T-CHM13 coordinate")
    ax.set_title("B. Reconstructed parental cis-methylation architecture", loc="left", fontsize=11, fontweight="bold")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=8)
    fig.suptitle("Figure 2. Maternal and paternal methylation architecture across 15q11-q13", y=0.98, fontsize=14, fontweight="bold")
    savefig(fig, outdir / "Figure2_parental_cis_architecture.png")


def figure3(outdir: Path, delta: pd.DataFrame, boundaries: Dict[str, float], domain_genes: pd.DataFrame,
            structural_bed: pd.DataFrame, domain_start: int, domain_end: int, maternal_label: str, paternal_label: str):
    fig = plt.figure(figsize=(15, 9))
    gs = GridSpec(3, 1, height_ratios=[0.8, 2.2, 0.8], figure=fig)
    ax0 = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[1, 0], sharex=ax0)
    ax2 = fig.add_subplot(gs[2, 0], sharex=ax0)
    plot_gene_track(ax0, domain_genes, domain_start, domain_end, "A. Gene context")

    if not delta.empty:
        ax.plot(delta["bin_mid"], delta["delta_methylation"], lw=0.8, alpha=0.4, label="raw Δ methylation")
        ax.plot(delta["bin_mid"], delta["delta_smooth"], lw=2.0, label="smoothed Δ methylation")
        ax.axhline(0, color="black", lw=0.8)
        ax.axhline(0.20, color="grey", lw=0.8, ls="--")
        ax.axhline(-0.20, color="grey", lw=0.8, ls="--")
        if np.isfinite(boundaries.get("left_boundary", np.nan)):
            ax.axvline(boundaries["left_boundary"], color="red", lw=1.5, ls="--", label="left boundary")
            ax.axvline(boundaries["right_boundary"], color="red", lw=1.5, ls="--", label="right boundary")
            ax.axvspan(boundaries["left_boundary"], boundaries["right_boundary"], color="red", alpha=0.06)
    ax.set_ylabel(f"Δ methylation\n{maternal_label} - {paternal_label}")
    ax.set_title("B. Data-driven boundary signal", loc="left", fontsize=11, fontweight="bold")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=8)

    ax2.set_ylim(0, 1)
    ax2.set_yticks([])
    ax2.set_xlabel("T2T-CHM13 coordinate")
    ax2.set_title("C. Optional structural annotations / breakpoint features", loc="left", fontsize=11, fontweight="bold")
    if structural_bed is not None and not structural_bed.empty:
        z = structural_bed[(structural_bed["chrom"] == domain_genes["chrom"].iloc[0]) &
                           (structural_bed["end"] > domain_start) & (structural_bed["start"] < domain_end)]
        for i, (_, r) in enumerate(z.iterrows()):
            y = 0.2 + (i % 4) * 0.15
            ax2.add_patch(Rectangle((r["start"], y), max(1, r["end"] - r["start"]), 0.09, alpha=0.45, ec="black", lw=0.2))
            if "name" in r:
                ax2.text((r["start"] + r["end"]) / 2, y + 0.12, str(r["name"]), ha="center", fontsize=6, rotation=45)
    else:
        ax2.text(0.5, 0.5, "No structural BED supplied", transform=ax2.transAxes, ha="center", va="center", fontsize=9)
    fig.suptitle("Figure 3. Methylation-domain boundary definition", y=0.98, fontsize=14, fontweight="bold")
    savefig(fig, outdir / "Figure3_boundary_definition.png")


def figure4_upd(outdir: Path, mean_prof: pd.DataFrame, domain_genes: pd.DataFrame, domain_start: int, domain_end: int):
    fig = plt.figure(figsize=(15, 8))
    gs = GridSpec(2, 1, height_ratios=[0.8, 3], figure=fig)
    ax0 = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[1, 0], sharex=ax0)
    plot_gene_track(ax0, domain_genes, domain_start, domain_end, "A. Gene context")
    labels = [x for x in ["control_maternal_like", "pwsdel_retained_maternal_like", "upd_maternal_like", "control_paternal_like"] if x in set(mean_prof["allele_label"])]
    for lab in labels:
        z = mean_prof[mean_prof["allele_label"] == lab].sort_values("bin_mid")
        ax.plot(z["bin_mid"], z["mean_methylation"], lw=1.8, label=lab.replace("_", " "))
    ax.set_ylim(-0.05, 1.05)
    ax.set_ylabel("Mean methylation")
    ax.set_xlabel("T2T-CHM13 coordinate")
    ax.set_title("B. Is PWS-mUPD equivalent to duplicated maternal architecture?", loc="left", fontsize=11, fontweight="bold")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=8)
    fig.suptitle("Figure 4. PWS-mUPD as a maternal-architecture test", y=0.98, fontsize=14, fontweight="bold")
    savefig(fig, outdir / "Figure4_PWS_mUPD_maternal_architecture.png")


def figure5_molecule_barcode(outdir: Path, mol: pd.DataFrame, region_label: str):
    if mol.empty:
        return
    # Create a compact matrix: rows are reads, columns are coarse CpG positions.
    x = mol.dropna(subset=["meth_state"]).copy()
    if x.empty:
        return
    # keep top reads with most sites
    read_counts = x.groupby(["sample", "read_id"]).size().sort_values(ascending=False).head(250)
    keep = set(read_counts.index)
    x["key"] = list(zip(x["sample"], x["read_id"]))
    x = x[x["key"].isin(keep)].copy()
    if x.empty:
        return
    # bin positions to reduce matrix width
    bins = np.linspace(x["pos"].min(), x["pos"].max(), 160)
    x["pos_bin"] = pd.cut(x["pos"], bins=bins, labels=False, include_lowest=True)
    x["row"] = x["sample"] + "|" + x["read_id"].astype(str).str.slice(0, 18)
    mat = x.pivot_table(index="row", columns="pos_bin", values="meth_state", aggfunc="mean")
    mat = mat.loc[mat.notna().sum(axis=1).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(13, max(5, min(12, 0.04 * len(mat) + 3))))
    im = ax.imshow(mat.values, aspect="auto", vmin=0, vmax=1, cmap="viridis", interpolation="nearest")
    ax.set_yticks([])
    ax.set_xlabel("Position bins across region")
    ax.set_ylabel("HiFi molecules")
    ax.set_title(f"Figure 5. Molecule-level methylation barcode: {region_label}", loc="left", fontsize=12, fontweight="bold")
    cb = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("Methylated state", fontsize=8)
    savefig(fig, outdir / "Figure5_molecule_level_barcode.png")


def figure6_structural_context(outdir: Path, svs: pd.DataFrame, cnv: pd.DataFrame, chrom: str, domain_start: int, domain_end: int):
    fig = plt.figure(figsize=(14, 9))
    gs = GridSpec(2, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])

    if not svs.empty:
        burden = svs.groupby(["group", "svtype"]).size().reset_index(name="n")
        piv = burden.pivot(index="group", columns="svtype", values="n").fillna(0)
        piv.plot(kind="bar", stacked=True, ax=ax1, width=0.8)
        ax1.set_ylabel("SV count")
        ax1.set_title("A. Genome-wide SV burden", loc="left", fontsize=10, fontweight="bold")
        ax1.legend(fontsize=7, frameon=False)

        chr15 = svs[(svs["chrom"] == chrom) & (svs["end"] > domain_start) & (svs["start"] < domain_end)].copy()
        if not chr15.empty:
            for i, (sample, z) in enumerate(chr15.groupby("sample")):
                y = i
                for _, r in z.iterrows():
                    ax3.plot([r["start"], r["end"]], [y, y], lw=3, label=r["svtype"] if i == 0 else None)
                ax3.text(domain_start, y, sample, va="center", ha="right", fontsize=7)
            ax3.set_xlim(domain_start, domain_end)
            ax3.set_xlabel("T2T-CHM13 coordinate")
            ax3.set_yticks([])
            ax3.set_title("C. SVs/CNV-relevant intervals across inferred domain", loc="left", fontsize=10, fontweight="bold")
        else:
            ax3.text(0.5, 0.5, "No SVs overlapping domain", transform=ax3.transAxes, ha="center", va="center")
    else:
        ax1.text(0.5, 0.5, "No SV VCF supplied", transform=ax1.transAxes, ha="center", va="center")
        ax1.axis("off")
        ax3.text(0.5, 0.5, "No SV VCF supplied", transform=ax3.transAxes, ha="center", va="center")
        ax3.axis("off")

    if cnv is not None and not cnv.empty:
        z = cnv[(cnv["chrom"] == chrom) & (cnv["end"] > domain_start) & (cnv["start"] < domain_end)].copy()
        # If value is numeric, plot; otherwise show interval counts.
        z["value_num"] = pd.to_numeric(z.get("value", 1), errors="coerce")
        if z["value_num"].notna().any():
            for sample, zz in z.groupby("sample") if "sample" in z.columns else [("CNV", z)]:
                ax2.plot((zz["start"] + zz["end"]) / 2, zz["value_num"], lw=1, label=sample)
            ax2.set_ylabel("CNV value")
            ax2.legend(fontsize=6, frameon=False, ncol=2)
        else:
            counts = z.groupby("group").size() if "group" in z.columns else pd.Series({"CNV": len(z)})
            counts.plot(kind="bar", ax=ax2)
            ax2.set_ylabel("CNV intervals")
        ax2.set_title("B. CNV context", loc="left", fontsize=10, fontweight="bold")
    else:
        ax2.text(0.5, 0.5, "No CNV files supplied", transform=ax2.transAxes, ha="center", va="center")
        ax2.axis("off")

    fig.suptitle("Figure 6. Structural and copy-number context", y=0.98, fontsize=14, fontweight="bold")
    savefig(fig, outdir / "Figure6_structural_context.png")


# -----------------------------------------------------------------------------
# Main pipeline
# -----------------------------------------------------------------------------

def _paper_vf_q1_pipeline_main():
    parser = argparse.ArgumentParser(description="Q1 PWS/AS methylation architecture pipeline for paper_vf")
    parser.add_argument("--methylation-dir", "--meth-dir", dest="meth_dir", default=DEFAULT_METHYLATION_DIR, type=Path,
                        help="Directory with methylation BED/DSS files")
    parser.add_argument("--gtf", default=DEFAULT_GTF, type=Path, help="T2T-CHM13 GTF file")
    parser.add_argument("--metadata", default=DEFAULT_METADATA_PATH, type=Path,
                        help="Metadata CSV with sample/group information")
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR, type=Path, help="Output directory")
    parser.add_argument("--chrom", default="chr15", help="Chromosome name in T2T files, default chr15")
    parser.add_argument("--exclude", default=",".join(sorted(DEFAULT_EXCLUDE)), help="Comma-separated sample codes to exclude")
    parser.add_argument("--methyl-format", default="auto", choices=["auto", "bed", "dss"], help="Methylation input format")
    parser.add_argument("--bed-meth-col", default="auto", help="1-based methylation column for BED-like files, or auto. Use 11 for modkit bedMethyl; 4 for simple BED.")
    parser.add_argument("--bed-cov-col", default="auto", help="1-based coverage column for BED-like files, none, or auto. Use 10 for modkit bedMethyl.")
    parser.add_argument("--chunksize", type=int, default=1_000_000, help="Rows per chunk when reading methylation files")
    parser.add_argument("--bin-size", type=int, default=5000, help="Bin size for domain architecture")
    parser.add_argument("--domain-padding", type=int, default=500000, help="Padding around canonical genes to infer domain")
    parser.add_argument("--ic-flank", type=int, default=5000, help="Flank around inferred SNRPN/SNURF/SNHG14 TSS for IC proxy")
    parser.add_argument("--ic-plot-flank", type=int, default=75000, help="Flank around IC for Figure 1 per-CpG plot")
    parser.add_argument("--boundary-delta", type=float, default=0.20, help="Absolute delta-beta threshold for domain boundary")
    parser.add_argument("--min-cpgs", type=int, default=3, help="Minimum CpGs per region/bin")

    parser.add_argument("--bam-dir", default=DEFAULT_BAM_DIR, type=Path, help="Optional BAM directory for MM/ML molecule-level methylation")
    parser.add_argument("--run-molecule", action="store_true", help="Run molecule-level barcode extraction from BAM files")
    parser.add_argument("--molecule-flank", type=int, default=15000, help="Flank around IC for molecule barcode")
    parser.add_argument("--vcf-dir", default=DEFAULT_VCF_DIR, type=Path, help="Optional SV VCF directory")
    parser.add_argument("--cnv-dir", default=DEFAULT_CNV_DIR, type=Path, help="Optional CNV bedGraph/BED directory")
    parser.add_argument("--structural-bed", default=None, type=Path, help="Optional T2T BED with BP/segdup/repeat annotations")

    args = parser.parse_args()

    outdir = mkdir(args.outdir)
    figdir = mkdir(outdir / "figures")
    tbldir = mkdir(outdir / "tables")
    exclude = {x.strip() for x in args.exclude.split(",") if x.strip()}

    log("Loading metadata and expected group states")
    metadata = load_metadata_table(args.metadata, outdir=outdir)

    log("Discovering methylation files")
    file_table = discover_methylation_files(args.meth_dir, exclude)
    file_table["group"] = file_table["sample"].apply(lambda s: sample_group(s, DEFAULT_GROUPS))
    file_table.to_csv(tbldir / "input_methylation_files.tsv", sep="\t", index=False)

    log("Loading T2T GTF and inferring domain/IC coordinates")
    genes = load_gtf_genes(args.gtf, args.chrom)
    domain_start, domain_end, domain_genes = infer_domain_from_gtf(
        genes, CANONICAL_DOMAIN_GENES, args.domain_padding, args.chrom
    )
    ic_chrom, ic_start, ic_end, ic_gene, ic_row = infer_ic_region_from_gtf(genes, IC_PRIORITY_GENES, args.ic_flank)
    if ic_chrom != args.chrom:
        raise ValueError(f"IC inferred on {ic_chrom}, but chrom argument is {args.chrom}")

    coord_summary = pd.DataFrame([
        {"feature": "inferred_domain", "chrom": args.chrom, "start": domain_start, "end": domain_end, "source": "T2T_GTF_gene_span_plus_padding"},
        {"feature": f"inferred_IC_proxy_{ic_gene}", "chrom": ic_chrom, "start": ic_start, "end": ic_end, "source": "T2T_GTF_TSS_plus_flank"},
    ])
    coord_summary.to_csv(tbldir / "inferred_T2T_coordinates.tsv", sep="\t", index=False)
    domain_genes.to_csv(tbldir / "domain_genes_from_T2T_GTF.tsv", sep="\t", index=False)

    log(f"Inferred domain: {args.chrom}:{domain_start:,}-{domain_end:,}")
    log(f"Inferred IC proxy: {ic_chrom}:{ic_start:,}-{ic_end:,} from {ic_gene}")

    log("RQ1: summarizing diagnostic IC methylation")
    ic_summary = summarize_region_for_files(file_table, args.chrom, ic_start, ic_end, args)
    ic_summary.to_csv(tbldir / "RQ1_IC_methylation_summary.tsv", sep="\t", index=False)
    allele_table = assign_allele_labels(ic_summary, args.min_cpgs)
    allele_table.to_csv(tbldir / "RQ1_haplotype_parental_like_assignment.tsv", sep="\t", index=False)

    ic_plot_start = max(0, ic_start - args.ic_plot_flank)
    ic_plot_end = ic_end + args.ic_plot_flank
    ic_single = load_single_cpg_region_for_plot(file_table, args.chrom, ic_plot_start, ic_plot_end, args)
    if not ic_single.empty:
        ic_single.to_csv(tbldir / "RQ1_IC_single_CpG_plot_data.tsv.gz", sep="\t", index=False, compression="gzip")

    figure1(figdir, file_table, domain_genes, ic_summary, ic_single, domain_start, domain_end, ic_plot_start, ic_plot_end)

    log("RQ2: loading binned methylation architecture across inferred domain")
    profiles = load_binned_profiles(file_table, args.chrom, domain_start, domain_end, args.bin_size, args)
    profiles = add_allele_labels(profiles, allele_table)
    profiles.to_csv(tbldir / "RQ2_binned_haplotype_methylation_profiles.tsv.gz", sep="\t", index=False, compression="gzip")

    mean_prof = mean_profile_by_label(profiles, min_cpg=args.min_cpgs)
    mean_prof.to_csv(tbldir / "RQ2_group_mean_parental_architecture.tsv", sep="\t", index=False)
    figure2(figdir, mean_prof, domain_genes, domain_start, domain_end)
    figure4_upd(figdir, mean_prof, domain_genes, domain_start, domain_end)

    log("RQ3: detecting methylation-domain boundaries")
    maternal_label, paternal_label = pick_architecture_labels(mean_prof)
    if maternal_label is not None and paternal_label is not None:
        delta = calculate_delta_profile(mean_prof, maternal_label, paternal_label)
        delta.to_csv(tbldir / "RQ3_maternal_paternal_delta_profile.tsv", sep="\t", index=False)
        boundaries = detect_boundaries(delta, threshold=args.boundary_delta)
    else:
        log("WARNING: could not identify maternal/paternal architecture labels for boundary detection")
        delta = pd.DataFrame()
        boundaries = {"left_boundary": np.nan, "right_boundary": np.nan, "n_bins_domain": 0}
        maternal_label = maternal_label or "NA"
        paternal_label = paternal_label or "NA"
    pd.DataFrame([{**boundaries, "maternal_label": maternal_label, "paternal_label": paternal_label, "threshold": args.boundary_delta}]).to_csv(
        tbldir / "RQ3_boundary_coordinates.tsv", sep="\t", index=False
    )

    structural_bed = pd.DataFrame()
    if args.structural_bed is not None and args.structural_bed.exists():
        structural_bed = read_bed_like(args.structural_bed)
        if not structural_bed.empty:
            structural_bed = structural_bed.rename(columns={"value": "name"})
            structural_bed.to_csv(tbldir / "RQ3_structural_annotation_bed_loaded.tsv", sep="\t", index=False)
    figure3(figdir, delta, boundaries, domain_genes, structural_bed, domain_start, domain_end, maternal_label, paternal_label)

    log("RQ4: molecule-level methylation coordination")
    mol_all = pd.DataFrame()
    if args.run_molecule:
        bam_table = discover_files_by_sample(args.bam_dir, (".bam",), exclude)
        bam_table.to_csv(tbldir / "RQ4_input_bam_files.tsv", sep="\t", index=False)
        mol_start = max(0, ic_start - args.molecule_flank)
        mol_end = ic_end + args.molecule_flank
        mol_rows = []
        for _, r in bam_table.iterrows():
            log(f"Extracting molecule methylation: {Path(r['path']).name}")
            df = extract_molecule_methylation_from_bam(Path(r["path"]), r["sample"], args.chrom, mol_start, mol_end)
            if not df.empty:
                mol_rows.append(df)
        if mol_rows:
            mol_all = pd.concat(mol_rows, ignore_index=True)
            mol_all.to_csv(tbldir / "RQ4_molecule_level_methylation_calls.tsv.gz", sep="\t", index=False, compression="gzip")
            # coordination score: per read fraction methylated and entropy proxy
            coord = mol_all.groupby(["sample", "group", "read_id", "haplotype_tag"]).agg(
                n_sites=("meth_state", "size"),
                mean_methylation=("meth_state", "mean"),
                sd_methylation=("meth_state", "std"),
            ).reset_index()
            coord["coordination_score"] = (coord["mean_methylation"] - 0.5).abs() * 2
            coord.to_csv(tbldir / "RQ4_read_coordination_scores.tsv", sep="\t", index=False)
            figure5_molecule_barcode(figdir, mol_all, f"{args.chrom}:{mol_start}-{mol_end}")
        else:
            log("No molecule-level methylation calls extracted. Check BAM MM/ML tags and contig names.")
    else:
        log("Skipping molecule-level extraction; use --run-molecule to enable.")

    log("RQ5: loading SV/CNV structural context")
    vcf_table = discover_files_by_sample(args.vcf_dir, (".vcf", ".vcf.gz"), exclude)
    if not vcf_table.empty:
        vcf_table.to_csv(tbldir / "RQ5_input_vcf_files.tsv", sep="\t", index=False)
    svs = load_all_svs(vcf_table) if not vcf_table.empty else pd.DataFrame()
    if not svs.empty:
        svs.to_csv(tbldir / "RQ5_all_SVs_parsed.tsv.gz", sep="\t", index=False, compression="gzip")
        domain_svs = svs[(svs["chrom"] == args.chrom) & (svs["end"] > domain_start) & (svs["start"] < domain_end)]
        domain_svs.to_csv(tbldir / "RQ5_domain_overlapping_SVs.tsv", sep="\t", index=False)

    cnv_table = discover_files_by_sample(args.cnv_dir, (".bed", ".bed.gz", ".bedgraph", ".bedGraph", ".tsv", ".txt"), exclude)
    cnv = pd.DataFrame()
    if not cnv_table.empty:
        cnv_table.to_csv(tbldir / "RQ5_input_cnv_files.tsv", sep="\t", index=False)
        cnv_rows = []
        for _, r in cnv_table.iterrows():
            z = read_bed_like(Path(r["path"]), sample=r["sample"])
            if not z.empty:
                cnv_rows.append(z)
        if cnv_rows:
            cnv = pd.concat(cnv_rows, ignore_index=True)
            cnv.to_csv(tbldir / "RQ5_CNV_intervals_loaded.tsv.gz", sep="\t", index=False, compression="gzip")

    # Overlap altered methylation bins with SVs as contextual table.
    if not svs.empty and not delta.empty:
        altered = delta[delta["delta_smooth"].abs() >= args.boundary_delta][["bin_start", "bin_end", "delta_smooth"]].copy()
        altered["chrom"] = args.chrom
        altered = altered.rename(columns={"bin_start": "start", "bin_end": "end"})
        ov = interval_overlap(altered[["chrom", "start", "end", "delta_smooth"]], svs)
        if not ov.empty:
            ov.to_csv(tbldir / "RQ5_altered_methylation_bins_overlapping_SVs.tsv", sep="\t", index=False)
    figure6_structural_context(figdir, svs, cnv, args.chrom, domain_start, domain_end)

    log("Writing final run summary")
    summary = {
        "n_methylation_files": len(file_table),
        "n_samples": file_table["sample"].nunique(),
        "metadata_path": str(args.metadata),
        "methylation_dir": str(args.meth_dir),
        "bam_dir": str(args.bam_dir),
        "vcf_dir": str(args.vcf_dir),
        "cnv_dir": str(args.cnv_dir),
        "gtf": str(args.gtf),
        "domain": f"{args.chrom}:{domain_start}-{domain_end}",
        "ic_proxy": f"{ic_chrom}:{ic_start}-{ic_end} ({ic_gene})",
        "figure_dir": str(figdir),
        "table_dir": str(tbldir),
        "boundary_left": boundaries.get("left_boundary"),
        "boundary_right": boundaries.get("right_boundary"),
    }
    pd.DataFrame([summary]).to_csv(outdir / "run_summary.tsv", sep="\t", index=False)
    log("Done.")

# --- Vendored from scripts/paper_vf/paper_vf_phase2_reciprocal_cis_architecture.py ---

"""
Phase 2 reciprocal cis-architecture analysis for the PWS/AS manuscript.

This script implements the fixed-coordinate Phase 2 prompts:
  - 1 kb per-haplotype methylation profiles across chr15:22,500,000-28,500,000
  - control maternal/paternal reference architecture
  - PWS-DEL retained maternal and AS-DEL retained paternal natural dissection
  - PWS-mUPD maternal/maternal architecture
  - a Figure 2 composite with a gene/IC track and reciprocal-delta boundary panel

The input methylation BEDs in this project are produced by PacBio
aligned_bam_to_cpg_scores. For those files the methylation percentage is column
9 and per-site read coverage is column 6; these are the defaults here because
the Phase 2 prompt explicitly filters on haplotype read coverage.
"""


import argparse
import ast
import gzip
import math
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter


PHASE2_CHROM = "chr15"
PHASE2_START = 17_600_000
PHASE2_END = 28_000_000
PHASE2_WINDOW = 1_000
AS_IC_CENTER = 24_920_000
DEFAULT_PHASE2_SAMPLES = [
    "001P", "002P", "004P", "005P", "006P", "007P",
    "013A", "014A", "016A", "017C", "018C",
]

MATERNAL = "#c7254e"
PATERNAL = "#2266aa"
MUPD_1 = "#d65f2d"
MUPD_2 = "#9b287b"
GREY = "#626262"
LIGHT_GREY = "#d0d0d0"
DEFAULT_ICR_BED = Path("/home/rare/arlen/reference/ICR_t2t.bed")
DEFAULT_SEG_DUP = Path("/home/rare/arlen/reference/dupseg")
DEFAULT_BP_SOURCE_SCRIPT = Path("/home/rare/arlen/scripts/Daniela/pws_chr15_hifi_deletion_analysis.py")
FALLBACK_PWS_BREAKPOINTS = {"BP1": 20_940_000, "BP2": 21_070_000, "BP3": 26_050_000}


def load_pws_breakpoints_from_daniela(script_path: Path = DEFAULT_BP_SOURCE_SCRIPT) -> Dict[str, int]:
    """Read BP1/BP2/BP3 from Daniela's CHM13 deletion-analysis script."""
    try:
        tree = ast.parse(script_path.read_text())
        for node in tree.body:
            if not isinstance(node, ast.Assign):
                continue
            if not any(isinstance(target, ast.Name) and target.id == "PWS_BREAKPOINTS" for target in node.targets):
                continue
            raw = ast.literal_eval(node.value)
            breakpoints = {str(k): int(v) for k, v in raw.items() if str(k) in {"BP1", "BP2", "BP3"}}
            if set(breakpoints) == {"BP1", "BP2", "BP3"}:
                return breakpoints
    except Exception:
        pass
    return dict(FALLBACK_PWS_BREAKPOINTS)


PWS_BREAKPOINTS_T2T = load_pws_breakpoints_from_daniela()
BP_HOTSPOTS_T2T = [(name, PWS_BREAKPOINTS_T2T[name]) for name in ("BP1", "BP2", "BP3")]
BP_CLUSTER_INTERVALS_T2T = [
    {
        "name": "BP1",
        "start": 17_691_439,
        "end": 20_454_275,
        "approx_size": "2.763 Mb",
        "interpretation": "Proximal SD block / BP1 cluster",
    },
    {
        "name": "BP2",
        "start": 20_753_698,
        "end": 21_183_655,
        "approx_size": "430 kb",
        "interpretation": "BP2 SD block",
    },
    {
        "name": "BP3",
        "start": 25_875_912,
        "end": 26_632_507,
        "approx_size": "757 kb",
        "interpretation": "Distal PWS/AS BP3 SD block",
    },
]


def log(msg: str) -> None:
    print(f"[Phase2] {msg}", file=sys.stderr, flush=True)


@dataclass(frozen=True)
class WindowSpec:
    chrom: str
    start: int
    end: int
    size: int

    def frame(self) -> pd.DataFrame:
        starts = np.arange(self.start, self.end, self.size, dtype=int)
        out = pd.DataFrame({"window_start": starts})
        out["window_end"] = np.minimum(out["window_start"] + self.size, self.end).astype(int)
        out["window_mid"] = ((out["window_start"] + out["window_end"]) / 2).astype(int)
        return out


def phase_samples_from_metadata(file_table: pd.DataFrame, requested: Optional[List[str]] = None) -> List[str]:
    if requested:
        present = set(file_table["sample"].unique())
        return [s for s in requested if s in present]
    wanted_groups = {"CONTROL", "PWS_DEL", "PWS_mUPD", "AS_DEL"}
    samples = []
    for sample in sorted(file_table["sample"].unique()):
        if sample_group(sample, DEFAULT_GROUPS) in wanted_groups:
            samples.append(sample)
    return samples


def one_based_column(value, default: Optional[int] = None) -> Optional[int]:
    if value is None:
        return default
    if isinstance(value, str):
        v = value.strip().lower()
        if v in {"", "none", "null", "na"}:
            return None
        if v == "auto":
            return default
        return int(v)
    return int(value)


def read_bed_region_fast(
    path: Path,
    chrom: str,
    start: int,
    end: int,
    bed_meth_col,
    bed_cov_col,
    timeout: int = 240,
) -> pd.DataFrame:
    """Extract a region from plain-text BED methylation files without loading the genome."""
    if path is None or not path.exists() or str(path).endswith(".gz"):
        return pd.DataFrame(columns=["chrom", "start", "end", "meth", "coverage", "mid"])

    meth_col = one_based_column(bed_meth_col, default=9) or 9
    cov_col = one_based_column(bed_cov_col, default=6)
    meth_i = meth_col - 1
    cov_i = cov_col - 1 if cov_col is not None else None
    max_i = max(meth_i, cov_i if cov_i is not None else 0, 2)

    try:
        proc = subprocess.run(
            [
                "awk",
                "-v",
                f"chrom={chrom}",
                "-v",
                f"start={start}",
                "-v",
                f"end={end}",
                '$1==chrom && $2>=start && $2<end {print}',
                str(path),
            ],
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except Exception as exc:
        log(f"WARNING: fast BED extraction failed for {path.name}: {exc}")
        return pd.DataFrame(columns=["chrom", "start", "end", "meth", "coverage", "mid"])

    rows = []
    for line in proc.stdout.splitlines():
        f = line.split()
        if len(f) <= max_i:
            continue
        try:
            s = int(float(f[1]))
            e = int(float(f[2]))
            meth = float(f[meth_i])
            cov = float(f[cov_i]) if cov_i is not None else 1.0
        except Exception:
            continue
        if meth > 1.5:
            meth /= 100.0
        if not np.isfinite(meth):
            continue
        rows.append(
            {
                "chrom": f[0],
                "start": s,
                "end": e,
                "meth": min(1.0, max(0.0, meth)),
                "coverage": cov if np.isfinite(cov) and cov > 0 else 0.0,
                "mid": int((s + e) / 2),
            }
        )
    if not rows:
        return pd.DataFrame(columns=["chrom", "start", "end", "meth", "coverage", "mid"])
    return pd.DataFrame(rows)


def read_phase2_methyl_region(path: Path, chrom: str, start: int, end: int, args) -> pd.DataFrame:
    fmt = args.methyl_format
    if fmt == "auto":
        fmt = "dss" if "dss" in path.name.lower() else "bed"
    if fmt == "bed" and not str(path).endswith(".gz"):
        return read_bed_region_fast(path, chrom, start, end, args.bed_meth_col, args.bed_cov_col)
    return read_methyl_region(
        path,
        chrom,
        start,
        end,
        fmt=fmt,
        bed_meth_col=args.bed_meth_col,
        bed_cov_col=args.bed_cov_col,
        chunksize=args.chunksize,
    )


def read_layer_windows(path: Optional[Path], spec: WindowSpec, args) -> pd.DataFrame:
    base = spec.frame()
    if path is None:
        return base.assign(n_cpg=0, n_reads=0.0, mean_methylation=np.nan)

    m = read_phase2_methyl_region(path, spec.chrom, spec.start, spec.end, args)
    if m.empty:
        return base.assign(n_cpg=0, n_reads=0.0, mean_methylation=np.nan)

    m["window_start"] = ((m["mid"] - spec.start) // spec.size) * spec.size + spec.start
    m["window_end"] = m["window_start"] + spec.size
    rows = []
    for (ws, we), z in m.groupby(["window_start", "window_end"], sort=True):
        rows.append(
            {
                "window_start": int(ws),
                "window_end": int(min(we, spec.end)),
                "n_cpg": int(len(z)),
                "n_reads": float(z["coverage"].sum()),
                "mean_methylation": weighted_mean(z["meth"].values, z["coverage"].values),
            }
        )
    z = pd.DataFrame(rows)
    out = base.merge(z, on=["window_start", "window_end"], how="left")
    out["n_cpg"] = out["n_cpg"].fillna(0).astype(int)
    out["n_reads"] = out["n_reads"].fillna(0.0)
    return out


def sample_layer_paths(file_table: pd.DataFrame, sample: str) -> Dict[str, Optional[Path]]:
    z = file_table[file_table["sample"] == sample].copy()
    out: Dict[str, Optional[Path]] = {"hap1": None, "hap2": None, "combined": None}
    for layer in out:
        hit = z[z["layer"] == layer]
        if not hit.empty:
            out[layer] = Path(hit.iloc[0]["path"])
    return out


def summarize_sample_windows(
    sample: str,
    file_table: pd.DataFrame,
    spec: WindowSpec,
    labels: Dict[str, Dict[str, str]],
    args,
) -> pd.DataFrame:
    paths = sample_layer_paths(file_table, sample)
    h1 = read_layer_windows(paths["hap1"], spec, args).rename(
        columns={
            "n_cpg": "n_CpGs_haplotype1",
            "n_reads": "n_reads_haplotype1",
            "mean_methylation": "mean_meth_haplotype1_raw",
        }
    )
    h2 = read_layer_windows(paths["hap2"], spec, args).rename(
        columns={
            "n_cpg": "n_CpGs_haplotype2",
            "n_reads": "n_reads_haplotype2",
            "mean_methylation": "mean_meth_haplotype2_raw",
        }
    )
    combined = read_layer_windows(paths["combined"], spec, args).rename(
        columns={
            "n_cpg": "n_CpGs_combined",
            "n_reads": "n_reads_combined",
            "mean_methylation": "mean_meth_combined_raw",
        }
    )

    out = h1.merge(h2, on=["window_start", "window_end", "window_mid"], how="outer")
    out = out.merge(combined, on=["window_start", "window_end", "window_mid"], how="outer")
    for c in ["n_CpGs_haplotype1", "n_CpGs_haplotype2", "n_CpGs_combined"]:
        out[c] = out[c].fillna(0).astype(int)
    for c in ["n_reads_haplotype1", "n_reads_haplotype2", "n_reads_combined"]:
        out[c] = out[c].fillna(0.0)

    out["n_CpGs"] = out["n_CpGs_combined"]
    missing_combined = out["n_CpGs"] == 0
    out.loc[missing_combined, "n_CpGs"] = out.loc[
        missing_combined, ["n_CpGs_haplotype1", "n_CpGs_haplotype2"]
    ].max(axis=1)

    out["mean_meth_haplotype1"] = out["mean_meth_haplotype1_raw"]
    out["mean_meth_haplotype2"] = out["mean_meth_haplotype2_raw"]
    out.loc[
        (out["n_reads_haplotype1"] < args.min_reads) | (out["n_CpGs_haplotype1"] < args.min_cpgs),
        "mean_meth_haplotype1",
    ] = np.nan
    out.loc[
        (out["n_reads_haplotype2"] < args.min_reads) | (out["n_CpGs_haplotype2"] < args.min_cpgs),
        "mean_meth_haplotype2",
    ] = np.nan
    out["mean_meth_combined"] = out["mean_meth_combined_raw"]
    out.loc[
        (out["n_reads_combined"] < args.min_reads) | (out["n_CpGs_combined"] < args.min_cpgs),
        "mean_meth_combined",
    ] = np.nan

    lab = labels.get(sample, {})
    out["sample"] = sample
    out["group"] = sample_group(sample, DEFAULT_GROUPS)
    out["parental_label_haplotype1"] = lab.get("hap1", f"{out['group'].iloc[0]}_hap1_unresolved")
    out["parental_label_haplotype2"] = lab.get("hap2", f"{out['group'].iloc[0]}_hap2_unresolved")

    requested_cols = [
        "window_start",
        "window_end",
        "n_CpGs",
        "mean_meth_haplotype1",
        "mean_meth_haplotype2",
        "parental_label_haplotype1",
        "parental_label_haplotype2",
        "n_reads_haplotype1",
        "n_reads_haplotype2",
    ]
    aux_cols = [
        "sample",
        "group",
        "window_mid",
        "n_CpGs_haplotype1",
        "n_CpGs_haplotype2",
        "n_CpGs_combined",
        "mean_meth_combined",
        "n_reads_combined",
    ]
    return out[aux_cols + requested_cols]


def summarize_ic_for_labels(file_table: pd.DataFrame, chrom: str, ic_start: int, ic_end: int, args) -> pd.DataFrame:
    rows = []
    for _, r in file_table[file_table["layer"].isin(["hap1", "hap2"])].iterrows():
        path = Path(r["path"])
        m = read_phase2_methyl_region(path, chrom, ic_start, ic_end, args)
        rows.append(
            {
                "sample": r["sample"],
                "group": sample_group(r["sample"], DEFAULT_GROUPS),
                "layer": r["layer"],
                "n_cpg": int(len(m)),
                "total_coverage": float(m["coverage"].sum()) if not m.empty else 0.0,
                "mean_methylation": weighted_mean(m["meth"].values, m["coverage"].values) if not m.empty else np.nan,
            }
        )
    return pd.DataFrame(rows)


def labels_from_ic_and_group(ic_summary: pd.DataFrame, min_cpgs: int) -> Dict[str, Dict[str, str]]:
    allele_table = assign_allele_labels(ic_summary, min_cpgs=min_cpgs)
    labels: Dict[str, Dict[str, str]] = {}
    for _, r in allele_table.iterrows():
        labels.setdefault(r["sample"], {})[r["layer"]] = r["allele_label"]

    for sample in ic_summary["sample"].unique():
        group = sample_group(sample, DEFAULT_GROUPS)
        labels.setdefault(sample, {})
        if group == "PWS_mUPD":
            labels[sample]["hap1"] = "upd_maternal_like"
            labels[sample]["hap2"] = "upd_maternal_like"
        elif group == "PWS_DEL":
            labels[sample].setdefault("hap1", "pwsdel_hap1_unresolved")
            labels[sample].setdefault("hap2", "pwsdel_hap2_unresolved")
        elif group == "AS_DEL":
            labels[sample].setdefault("hap1", "asdel_hap1_unresolved")
            labels[sample].setdefault("hap2", "asdel_hap2_unresolved")
    return labels


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def nearby_plot_interval(start: int, end: int, spec: WindowSpec, rank: int, side: str, width: int = 8_000) -> Tuple[int, int]:
    if side == "left":
        x0 = spec.start + rank * 15_000
        return x0, min(x0 + width, spec.end)
    x1 = spec.end - rank * 15_000
    return max(spec.start, x1 - width), x1


def assign_plot_coordinates(df: pd.DataFrame, spec: WindowSpec, width: int = 8_000) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    out["within_phase2_window"] = (out["end"] > spec.start) & (out["start"] < spec.end)
    out["plot_start"] = out["start"].clip(lower=spec.start, upper=spec.end)
    out["plot_end"] = out["end"].clip(lower=spec.start, upper=spec.end)
    left_idx = out[out["end"] <= spec.start].sort_values("start").index
    for rank, idx in enumerate(left_idx):
        out.loc[idx, ["plot_start", "plot_end"]] = nearby_plot_interval(
            int(out.loc[idx, "start"]), int(out.loc[idx, "end"]), spec, rank, "left", width
        )
    right_idx = out[out["start"] >= spec.end].sort_values("start").index
    for rank, idx in enumerate(right_idx):
        out.loc[idx, ["plot_start", "plot_end"]] = nearby_plot_interval(
            int(out.loc[idx, "start"]), int(out.loc[idx, "end"]), spec, rank, "right", width
        )
    return out


def read_icr_annotations(path: Path, spec: WindowSpec, flank: int = 1_500_000) -> pd.DataFrame:
    rows = []
    if path is None or not Path(path).exists():
        return pd.DataFrame()
    with open_text(Path(path)) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] != spec.chrom:
                continue
            try:
                start = int(float(f[1]))
                end = int(float(f[2]))
            except Exception:
                continue
            if end <= spec.start - flank or start >= spec.end + flank:
                continue
            rows.append(
                {
                    "chrom": f[0],
                    "start": start,
                    "end": end,
                    "name": f[3] if len(f) > 3 else "ICR",
                    "parent": f[4] if len(f) > 4 else "",
                    "gene": f[5] if len(f) > 5 else "",
                    "source": str(path),
                }
            )
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows).sort_values(["start", "end"]).reset_index(drop=True)
    return assign_plot_coordinates(out, spec, width=7_000)


def read_segdup_annotations(path: Path, spec: WindowSpec) -> Tuple[pd.DataFrame, pd.DataFrame]:
    raw_rows = []
    if path is None or not Path(path).exists():
        return pd.DataFrame(), pd.DataFrame()
    with open_text(Path(path)) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] != spec.chrom:
                continue
            try:
                start = int(float(f[1]))
                end = int(float(f[2]))
            except Exception:
                continue
            if end <= spec.start or start >= spec.end:
                continue
            raw_rows.append(
                {
                    "chrom": f[0],
                    "start": max(start, spec.start),
                    "end": min(end, spec.end),
                    "name": f[3] if len(f) > 3 else "segdup",
                    "target": f[9] if len(f) > 9 else "",
                    "target_start": f[10] if len(f) > 10 else "",
                    "target_end": f[11] if len(f) > 11 else "",
                    "frac_match": f[23] if len(f) > 23 else "",
                    "aln_len": f[15] if len(f) > 15 else "",
                }
            )
    if not raw_rows:
        return pd.DataFrame(), pd.DataFrame()
    raw = pd.DataFrame(raw_rows).sort_values(["start", "end"]).reset_index(drop=True)

    merged = []
    for _, r in raw.sort_values(["start", "end"]).iterrows():
        s, e = int(r["start"]), int(r["end"])
        if not merged or s > merged[-1]["end"]:
            merged.append({"chrom": spec.chrom, "start": s, "end": e, "n_raw_segments": 1})
        else:
            merged[-1]["end"] = max(merged[-1]["end"], e)
            merged[-1]["n_raw_segments"] += 1
    blocks = pd.DataFrame(merged)
    blocks["size_bp"] = blocks["end"] - blocks["start"]
    blocks = blocks[blocks["size_bp"] >= 1_000].reset_index(drop=True)
    return raw, blocks


def nearest_interval(pos: int, intervals: pd.DataFrame, label_col: str) -> Tuple[str, int, str]:
    if intervals is None or intervals.empty:
        return "", np.nan, ""
    z = intervals.copy()
    starts = pd.to_numeric(z["start"], errors="coerce")
    ends = pd.to_numeric(z["end"], errors="coerce")
    dist = np.where((starts <= pos) & (ends >= pos), 0, np.minimum((starts - pos).abs(), (ends - pos).abs()))
    idx = int(np.nanargmin(dist))
    row = z.iloc[idx]
    label = str(row.get(label_col, ""))
    interval = f"{int(row['start'])}-{int(row['end'])}"
    return label, int(dist[idx]), interval


def annotate_deletion_breakpoints(
    deletion_intervals: pd.DataFrame,
    genes: pd.DataFrame,
    icrs: pd.DataFrame,
    segdup_blocks: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    if deletion_intervals.empty:
        return pd.DataFrame()
    gene_intervals = genes[["gene", "start", "end", "strand"]].copy() if not genes.empty else pd.DataFrame()
    icr_intervals = icrs.rename(columns={"name": "icr_name"}).copy() if icrs is not None and not icrs.empty else pd.DataFrame()
    segdup_intervals = segdup_blocks.copy() if segdup_blocks is not None and not segdup_blocks.empty else pd.DataFrame()
    if not segdup_intervals.empty:
        segdup_intervals["segdup_block"] = segdup_intervals["start"].astype(str) + "-" + segdup_intervals["end"].astype(str)

    for _, r in deletion_intervals.iterrows():
        for kind, pos in [("left_or_start", int(r["start"])), ("right_or_end", int(r["end"]))]:
            gene, gene_dist, gene_interval = nearest_interval(pos, gene_intervals, "gene")
            icr, icr_dist, icr_interval = nearest_interval(pos, icr_intervals, "icr_name")
            segdup, segdup_dist, segdup_interval = nearest_interval(pos, segdup_intervals, "segdup_block")
            rows.append(
                {
                    "sample": r["sample"],
                    "group": r["group"],
                    "chrom": r["chrom"],
                    "breakpoint_position": pos,
                    "breakpoint_side": kind,
                    "copy_number": r.get("copy_number", np.nan),
                    "nearest_gene": gene,
                    "distance_to_nearest_gene_bp": gene_dist,
                    "nearest_gene_interval": gene_interval,
                    "nearest_ICR": icr,
                    "distance_to_nearest_ICR_bp": icr_dist,
                    "nearest_ICR_interval": icr_interval,
                    "nearest_or_overlapping_segdup_block": segdup,
                    "distance_to_segdup_block_bp": segdup_dist,
                    "segdup_block_interval": segdup_interval,
                    "source": r.get("source", ""),
                }
            )
    return pd.DataFrame(rows).sort_values(["group", "sample", "breakpoint_position", "breakpoint_side"]).reset_index(drop=True)


def discover_cnv_files(cnv_dir: Path, samples: Iterable[str]) -> pd.DataFrame:
    suffixes = (".bedgraph", ".bedGraph", ".bed", ".bed.gz", ".tsv", ".txt")
    wanted = set(samples)
    rows = []
    if cnv_dir is None or not cnv_dir.exists():
        return pd.DataFrame(columns=["sample", "path"])
    for p in cnv_dir.rglob("*"):
        if not p.is_file() or not any(str(p).endswith(s) for s in suffixes):
            continue
        sample = detect_sample_code(p)
        if sample in wanted:
            rows.append({"sample": sample, "path": str(p)})
    return pd.DataFrame(rows).drop_duplicates().sort_values(["sample", "path"]).reset_index(drop=True)


def read_cnv_deletion_intervals(
    cnv_table: pd.DataFrame,
    spec: WindowSpec,
    copy_threshold: float,
) -> pd.DataFrame:
    rows = []
    for _, r in cnv_table.iterrows():
        path = Path(r["path"])
        try:
            with open_text(path) as fh:
                for line in fh:
                    if not line.strip() or line.startswith("#"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 4 or f[0] != spec.chrom:
                        continue
                    try:
                        s = int(float(f[1]))
                        e = int(float(f[2]))
                        value = float(f[3])
                    except Exception:
                        continue
                    if e <= spec.start or s >= spec.end or value > copy_threshold:
                        continue
                    rows.append(
                        {
                            "sample": r["sample"],
                            "group": sample_group(r["sample"], DEFAULT_GROUPS),
                            "chrom": spec.chrom,
                            "start": max(s, spec.start),
                            "end": min(e, spec.end),
                            "copy_number": value,
                            "source": str(path),
                        }
                    )
        except Exception as exc:
            log(f"WARNING: failed to read CNV file {path}: {exc}")
    if not rows:
        return pd.DataFrame(columns=["sample", "group", "chrom", "start", "end", "copy_number", "source"])
    out = pd.DataFrame(rows)
    out = out[out["end"] > out["start"]].copy()
    return out.sort_values(["sample", "start", "end"]).reset_index(drop=True)


def windows_in_intervals(windows: pd.DataFrame, intervals: pd.DataFrame, sample: Optional[str] = None) -> pd.Series:
    mask = pd.Series(False, index=windows.index)
    if intervals.empty:
        return mask
    z = intervals
    if sample is not None and "sample" in z.columns:
        z = z[z["sample"] == sample]
    for _, r in z.iterrows():
        mask |= (windows["window_end"] > int(r["start"])) & (windows["window_start"] < int(r["end"]))
    return mask


def infer_deletion_intervals_from_windows(all_windows: pd.DataFrame, spec: WindowSpec, min_reads: int) -> pd.DataFrame:
    rows = []
    for sample, z in all_windows.groupby("sample"):
        group = sample_group(sample, DEFAULT_GROUPS)
        if group not in {"PWS_DEL", "AS_DEL"}:
            continue
        zz = z.sort_values("window_start").copy()
        candidate = (
            (zz["n_reads_combined"] >= min_reads)
            & ((zz["n_reads_haplotype1"] < min_reads) | (zz["n_reads_haplotype2"] < min_reads))
        )
        for start_i, end_i in boolean_runs(candidate.to_numpy()):
            start = int(zz.iloc[start_i]["window_start"])
            end = int(zz.iloc[end_i]["window_end"])
            if end - start >= 50_000:
                rows.append(
                    {
                        "sample": sample,
                        "group": group,
                        "chrom": spec.chrom,
                        "start": start,
                        "end": end,
                        "copy_number": np.nan,
                        "source": "coverage_fallback",
                    }
                )
    return pd.DataFrame(rows)


def boolean_runs(mask: np.ndarray) -> List[Tuple[int, int]]:
    runs: List[Tuple[int, int]] = []
    start: Optional[int] = None
    for i, val in enumerate(mask):
        if bool(val) and start is None:
            start = i
        if start is not None and (not bool(val) or i == len(mask) - 1):
            end = i if bool(val) and i == len(mask) - 1 else i - 1
            runs.append((start, end))
            start = None
    return runs


def update_del_labels_by_coverage(
    labels: Dict[str, Dict[str, str]],
    all_windows: pd.DataFrame,
    deletion_intervals: pd.DataFrame,
    fold: float = 1.5,
) -> Dict[str, Dict[str, str]]:
    labels = {s: dict(v) for s, v in labels.items()}
    for sample in sorted(all_windows["sample"].unique()):
        group = sample_group(sample, DEFAULT_GROUPS)
        if group not in {"PWS_DEL", "AS_DEL"}:
            continue
        z = all_windows[all_windows["sample"] == sample].copy()
        mask = windows_in_intervals(z, deletion_intervals, sample=sample)
        if not mask.any():
            continue
        h1 = float(z.loc[mask, "n_reads_haplotype1"].sum())
        h2 = float(z.loc[mask, "n_reads_haplotype2"].sum())
        labels.setdefault(sample, {})
        if max(h1, h2) < 50:
            continue
        if h1 >= h2 * fold:
            retained = "hap1"
            absent = "hap2"
        elif h2 >= h1 * fold:
            retained = "hap2"
            absent = "hap1"
        else:
            continue
        if group == "PWS_DEL":
            labels[sample][retained] = "pwsdel_retained_maternal_like"
            labels[sample][absent] = "pwsdel_deleted_paternal_like"
        else:
            labels[sample][retained] = "asdel_retained_paternal_like"
            labels[sample][absent] = "asdel_deleted_maternal_like"
    return labels


def apply_labels_to_windows(all_windows: pd.DataFrame, labels: Dict[str, Dict[str, str]]) -> pd.DataFrame:
    out = all_windows.copy()
    out["parental_label_haplotype1"] = out["sample"].map(lambda s: labels.get(s, {}).get("hap1", "hap1_unresolved"))
    out["parental_label_haplotype2"] = out["sample"].map(lambda s: labels.get(s, {}).get("hap2", "hap2_unresolved"))
    return out


def write_per_sample_tables(all_windows: pd.DataFrame, outdir: Path) -> None:
    per_dir = mkdir(outdir / "tables" / "phase2_per_sample_windows")
    requested_cols = [
        "window_start",
        "window_end",
        "n_CpGs",
        "mean_meth_haplotype1",
        "mean_meth_haplotype2",
        "parental_label_haplotype1",
        "parental_label_haplotype2",
        "n_reads_haplotype1",
        "n_reads_haplotype2",
    ]
    for sample, z in all_windows.groupby("sample"):
        z = z.sort_values("window_start")
        z[requested_cols].to_csv(
            per_dir / f"Phase2_{sample}_haplotype_1kb_methylation.tsv",
            sep="\t",
            index=False,
            na_rep="NaN",
        )


def value_for_label(z: pd.DataFrame, label_contains: str) -> pd.Series:
    h1 = z["parental_label_haplotype1"].astype(str).str.contains(label_contains, regex=False)
    h2 = z["parental_label_haplotype2"].astype(str).str.contains(label_contains, regex=False)
    out = pd.Series(np.nan, index=z.index, dtype=float)
    out.loc[h1] = z.loc[h1, "mean_meth_haplotype1"]
    out.loc[h2] = z.loc[h2, "mean_meth_haplotype2"]
    return out


def control_reference(all_windows: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    controls = all_windows[all_windows["group"] == "CONTROL"].copy()
    rows = []
    sample_rows = []
    for sample, z in controls.groupby("sample"):
        zz = z.sort_values("window_start")
        maternal = value_for_label(zz, "maternal")
        paternal = value_for_label(zz, "paternal")
        tmp = zz[["sample", "window_start", "window_end", "window_mid"]].copy()
        tmp["control_maternal"] = maternal
        tmp["control_paternal"] = paternal
        sample_rows.append(tmp)
    sample_ref = pd.concat(sample_rows, ignore_index=True) if sample_rows else pd.DataFrame()
    if sample_ref.empty:
        return pd.DataFrame(), sample_ref
    for (ws, we, wm), z in sample_ref.groupby(["window_start", "window_end", "window_mid"]):
        rows.append(
            {
                "window_start": ws,
                "window_end": we,
                "window_mid": wm,
                "control_maternal_mean": float(np.nanmean(z["control_maternal"])),
                "control_paternal_mean": float(np.nanmean(z["control_paternal"])),
                "n_controls_maternal": int(z["control_maternal"].notna().sum()),
                "n_controls_paternal": int(z["control_paternal"].notna().sum()),
            }
        )
    ref = pd.DataFrame(rows).sort_values("window_start").reset_index(drop=True)
    ref["maternal_paternal_delta"] = ref["control_maternal_mean"] - ref["control_paternal_mean"]
    ref["abs_delta"] = ref["maternal_paternal_delta"].abs()
    return ref, sample_ref


def retained_profiles(all_windows: pd.DataFrame, deletion_intervals: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sample, z in all_windows.groupby("sample"):
        group = sample_group(sample, DEFAULT_GROUPS)
        if group not in {"PWS_DEL", "AS_DEL"}:
            continue
        zz = z.sort_values("window_start").copy()
        in_del = windows_in_intervals(zz, deletion_intervals, sample=sample)
        retained_kind = "maternal" if group == "PWS_DEL" else "paternal"
        outside = value_for_label(zz, retained_kind)
        retained = outside.copy()
        retained.loc[in_del] = zz.loc[in_del, "mean_meth_combined"]
        source = pd.Series("assigned_haplotype", index=zz.index, dtype=object)
        source.loc[in_del] = "combined_hemizygous_deletion"
        absent_cov = np.minimum(zz["n_reads_haplotype1"], zz["n_reads_haplotype2"])
        if group == "PWS_DEL":
            hap2_cov = zz["n_reads_haplotype2"]
        else:
            hap2_cov = absent_cov
        tmp = zz[
            [
                "sample",
                "group",
                "window_start",
                "window_end",
                "window_mid",
                "n_reads_haplotype1",
                "n_reads_haplotype2",
                "n_reads_combined",
            ]
        ].copy()
        tmp["retained_parent"] = retained_kind
        tmp["retained_methylation"] = retained
        tmp["retained_source"] = source
        tmp["within_deletion_footprint"] = in_del.to_numpy()
        tmp["absent_haplotype_coverage"] = absent_cov
        tmp["prompt_haplotype2_coverage"] = hap2_cov
        rows.append(tmp)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def pearson_summary(a: pd.Series, b: pd.Series) -> Dict[str, float]:
    m = pd.DataFrame({"a": a, "b": b}).dropna()
    if len(m) < 3 or m["a"].std() == 0 or m["b"].std() == 0:
        r = np.nan
    else:
        r = float(m["a"].corr(m["b"]))
    return {
        "n_windows": int(len(m)),
        "pearson_r": r,
        "mean_abs_delta": float((m["a"] - m["b"]).abs().mean()) if len(m) else np.nan,
        "mean_delta": float((m["a"] - m["b"]).mean()) if len(m) else np.nan,
    }


def compute_correlations(
    retained: pd.DataFrame,
    all_windows: pd.DataFrame,
    control_ref: pd.DataFrame,
    deletion_intervals: pd.DataFrame,
    deviation_threshold: float,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    corr_rows = []
    for sample, z in retained.groupby("sample"):
        group = sample_group(sample, DEFAULT_GROUPS)
        if group not in {"PWS_DEL", "AS_DEL"}:
            continue
        ref_col = "control_maternal_mean" if group == "PWS_DEL" else "control_paternal_mean"
        zz = z.merge(control_ref[["window_start", ref_col]], on="window_start", how="left")
        zz = zz[zz["within_deletion_footprint"]]
        stats = pearson_summary(zz["retained_methylation"], zz[ref_col])
        intervals = deletion_intervals[deletion_intervals["sample"] == sample]
        corr_rows.append(
            {
                "sample": sample,
                "group": group,
                "comparison": f"{sample}_retained_vs_{ref_col}",
                "reference": ref_col,
                "deletion_start": int(intervals["start"].min()) if not intervals.empty else np.nan,
                "deletion_end": int(intervals["end"].max()) if not intervals.empty else np.nan,
                **stats,
            }
        )

    deviations = []
    upd = all_windows[all_windows["sample"] == "004P"].copy()
    if not upd.empty and not control_ref.empty:
        zz = upd.merge(control_ref[["window_start", "control_maternal_mean"]], on="window_start", how="left")
        for hap, col in [("haplotype1", "mean_meth_haplotype1"), ("haplotype2", "mean_meth_haplotype2")]:
            stats = pearson_summary(zz[col], zz["control_maternal_mean"])
            corr_rows.append(
                {
                    "sample": "004P",
                    "group": "PWS_mUPD",
                    "comparison": f"004P_{hap}_vs_control_maternal",
                    "reference": "control_maternal_mean",
                    "deletion_start": np.nan,
                    "deletion_end": np.nan,
                    **stats,
                }
            )
            diff = zz[col] - zz["control_maternal_mean"]
            hit = zz[diff.abs() >= deviation_threshold].copy()
            if not hit.empty:
                hit["haplotype"] = hap
                hit["methylation"] = hit[col]
                hit["control_maternal_mean"] = hit["control_maternal_mean"]
                hit["delta_from_control_maternal"] = diff.loc[hit.index]
                deviations.append(
                    hit[
                        [
                            "sample",
                            "haplotype",
                            "window_start",
                            "window_end",
                            "methylation",
                            "control_maternal_mean",
                            "delta_from_control_maternal",
                        ]
                    ]
                )

    corr = pd.DataFrame(corr_rows)
    dev = pd.concat(deviations, ignore_index=True) if deviations else pd.DataFrame(
        columns=[
            "sample",
            "haplotype",
            "window_start",
            "window_end",
            "methylation",
            "control_maternal_mean",
            "delta_from_control_maternal",
        ]
    )
    return corr, dev


def reciprocal_delta_and_boundaries(
    retained: pd.DataFrame,
    spec: WindowSpec,
    genes: pd.DataFrame,
    icrs: pd.DataFrame,
    segdup_blocks: pd.DataFrame,
    threshold: float = 0.10,
    smooth_window: int = 31,
    min_bp: int = 10_000,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if retained.empty:
        return pd.DataFrame(), pd.DataFrame()

    base = spec.frame()

    def summarize_group(group: str, mean_col: str, n_col: str) -> pd.DataFrame:
        z = retained[retained["group"] == group]
        if z.empty:
            return base[["window_start", "window_end", "window_mid"]].assign(**{mean_col: np.nan, n_col: 0})
        out = (
            z.groupby(["window_start", "window_end", "window_mid"])["retained_methylation"]
            .agg([("mean_value", "mean"), ("n_value", "count")])
            .reset_index()
            .rename(columns={"mean_value": mean_col, "n_value": n_col})
        )
        return base.merge(out, on=["window_start", "window_end", "window_mid"], how="left")

    pws = summarize_group("PWS_DEL", "pwsdel_retained_maternal_mean", "n_pwsdel_samples")
    asdel = summarize_group("AS_DEL", "asdel_retained_paternal_mean", "n_asdel_samples")
    profile = pws.merge(
        asdel[["window_start", "asdel_retained_paternal_mean", "n_asdel_samples"]],
        on="window_start",
        how="left",
    )
    profile = profile.sort_values("window_start").reset_index(drop=True)
    profile["reciprocal_delta"] = profile["pwsdel_retained_maternal_mean"] - profile["asdel_retained_paternal_mean"]
    profile["reciprocal_delta_smoothed"] = smooth_values(profile["reciprocal_delta"], window=smooth_window)
    profile["abs_reciprocal_delta_smoothed"] = profile["reciprocal_delta_smoothed"].abs()
    profile["boundary_threshold"] = threshold
    profile["delta_state"] = "low_delta"
    profile.loc[profile["reciprocal_delta_smoothed"] >= threshold, "delta_state"] = "maternal_higher"
    profile.loc[profile["reciprocal_delta_smoothed"] <= -threshold, "delta_state"] = "paternal_higher"

    gene_intervals = genes[["gene", "start", "end", "strand"]].copy() if genes is not None and not genes.empty else pd.DataFrame()
    icr_intervals = icrs.rename(columns={"name": "icr_name"}).copy() if icrs is not None and not icrs.empty else pd.DataFrame()
    segdup_intervals = segdup_blocks.copy() if segdup_blocks is not None and not segdup_blocks.empty else pd.DataFrame()
    if not segdup_intervals.empty:
        segdup_intervals["segdup_block"] = segdup_intervals["start"].astype(str) + "-" + segdup_intervals["end"].astype(str)

    rows = []
    for state in ["maternal_higher", "paternal_higher"]:
        mask = (profile["delta_state"] == state).fillna(False).to_numpy()
        for start_i, end_i in boolean_runs(mask):
            start = int(profile.iloc[start_i]["window_start"])
            end = int(profile.iloc[end_i]["window_end"])
            if end - start < min_bp:
                continue
            z = profile.iloc[start_i : end_i + 1]
            row = {
                "chrom": spec.chrom,
                "start": start,
                "end": end,
                "size_bp": end - start,
                "delta_state": state,
                "interpretation": "PWS_DEL_maternal_higher_than_AS_DEL_paternal"
                if state == "maternal_higher"
                else "AS_DEL_paternal_higher_than_PWS_DEL_maternal",
                "n_windows": int(end_i - start_i + 1),
                "mean_delta_smoothed": float(z["reciprocal_delta_smoothed"].mean()),
                "max_abs_delta_smoothed": float(z["abs_reciprocal_delta_smoothed"].max()),
                "threshold_abs_delta": threshold,
                "smoothing_windows": smooth_window,
            }
            for prefix, pos in [("left_boundary", start), ("right_boundary", end)]:
                gene, gene_dist, gene_interval = nearest_interval(pos, gene_intervals, "gene")
                icr, icr_dist, icr_interval = nearest_interval(pos, icr_intervals, "icr_name")
                segdup, segdup_dist, segdup_interval = nearest_interval(pos, segdup_intervals, "segdup_block")
                row.update(
                    {
                        f"{prefix}_position": pos,
                        f"{prefix}_nearest_gene": gene,
                        f"{prefix}_distance_to_nearest_gene_bp": gene_dist,
                        f"{prefix}_nearest_gene_interval": gene_interval,
                        f"{prefix}_nearest_ICR": icr,
                        f"{prefix}_distance_to_nearest_ICR_bp": icr_dist,
                        f"{prefix}_nearest_ICR_interval": icr_interval,
                        f"{prefix}_nearest_or_overlapping_segdup_block": segdup,
                        f"{prefix}_distance_to_segdup_block_bp": segdup_dist,
                        f"{prefix}_segdup_block_interval": segdup_interval,
                    }
                )
            rows.append(row)

    boundaries = pd.DataFrame(rows)
    if not boundaries.empty:
        boundaries = boundaries.sort_values(["start", "end", "delta_state"]).reset_index(drop=True)
    return profile, boundaries


def classify_control_regions(
    control_ref: pd.DataFrame,
    divergent_threshold: float,
    convergent_threshold: float,
    min_bp: int,
) -> pd.DataFrame:
    if control_ref.empty:
        return pd.DataFrame()
    rows = []
    z = control_ref.sort_values("window_start").reset_index(drop=True)
    for label, mask in [
        ("imprinted_divergent", (z["abs_delta"] >= divergent_threshold).fillna(False).to_numpy()),
        ("biallelic_convergent", (z["abs_delta"] <= convergent_threshold).fillna(False).to_numpy()),
    ]:
        for start_i, end_i in boolean_runs(mask):
            start = int(z.iloc[start_i]["window_start"])
            end = int(z.iloc[end_i]["window_end"])
            if end - start < min_bp:
                continue
            rows.append(
                {
                    "region_class": label,
                    "chrom": PHASE2_CHROM,
                    "start": start,
                    "end": end,
                    "size_bp": end - start,
                    "mean_abs_delta": float(z.iloc[start_i : end_i + 1]["abs_delta"].mean()),
                    "n_windows": int(end_i - start_i + 1),
                }
            )
    return pd.DataFrame(rows).sort_values(["region_class", "start"]).reset_index(drop=True) if rows else pd.DataFrame()


def feature_from_genes(genes: pd.DataFrame, name: str, patterns: List[str], label: Optional[str] = None) -> Optional[dict]:
    hits = []
    for pat in patterns:
        if pat.endswith("*"):
            prefix = pat[:-1].upper()
            z = genes[genes["gene"].str.upper().str.startswith(prefix)]
        else:
            z = genes[genes["gene"].str.upper() == pat.upper()]
        if not z.empty:
            hits.append(z)
    if not hits:
        return None
    h = pd.concat(hits, ignore_index=True)
    return {
        "feature": label or name,
        "chrom": h["chrom"].iloc[0],
        "start": int(h["start"].min()),
        "end": int(h["end"].max()),
        "strand": h["strand"].iloc[0] if h["strand"].nunique() == 1 else ".",
        "source": "T2T_GTF",
    }


def build_phase2_gene_features(genes: pd.DataFrame, pws_ic_start: int, pws_ic_end: int, spec: WindowSpec) -> pd.DataFrame:
    specs = [
        ("MKRN3", ["MKRN3"], "MKRN3"),
        ("MAGEL2", ["MAGEL2"], "MAGEL2"),
        ("NDN", ["NDN"], "NDN"),
        ("NPAP1", ["NPAP1"], "NPAP1"),
        ("SNURF-SNRPN", ["SNURF", "SNRPN"], "SNURF-SNRPN"),
        ("SNORD116 cluster", ["SNORD116*"], "SNORD116"),
        ("SNORD115 cluster", ["SNORD115*"], "SNORD115"),
        ("UBE3A-ATS", ["SNHG14", "UBE3A-AS1"], "UBE3A-ATS"),
        ("UBE3A", ["UBE3A"], "UBE3A"),
    ]
    rows = [
        {
            "feature": "PWS-IC",
            "chrom": PHASE2_CHROM,
            "start": pws_ic_start,
            "end": pws_ic_end,
            "strand": ".",
            "source": "SNRPN_TSS_proxy",
        },
        {
            "feature": "AS-IC",
            "chrom": PHASE2_CHROM,
            "start": AS_IC_CENTER - 10_000,
            "end": AS_IC_CENTER + 10_000,
            "strand": ".",
            "source": "prompt_coordinate_approximation",
        },
    ]
    for name, patterns, label in specs:
        f = feature_from_genes(genes, name, patterns, label)
        if f is not None:
            rows.append(f)
    out = pd.DataFrame(rows)
    out = out.drop_duplicates("feature", keep="first")
    out = assign_plot_coordinates(out, spec, width=9_000)
    left = out["end"] <= spec.start
    right = out["start"] >= spec.end
    out = out[out["within_phase2_window"] | left | right].copy()
    return out.sort_values(["start", "end"]).reset_index(drop=True)


def add_region_spans(ax, regions: pd.DataFrame, y_text: float = 1.02) -> None:
    if regions.empty:
        return
    div = regions[regions["region_class"] == "imprinted_divergent"].copy()
    if div.empty:
        return
    div = div.sort_values("size_bp", ascending=False).head(3)
    for i, (_, r) in enumerate(div.iterrows()):
        ax.axvspan(r["start"], r["end"], color="#f7c6c7", alpha=0.24, lw=0)
        if i == 0:
            x = (r["start"] + r["end"]) / 2
            ax.text(x, y_text, "imprinted divergence", ha="center", va="bottom", fontsize=7, color="#8a1c2b")


def add_feature_bands(ax, features: pd.DataFrame) -> None:
    if features.empty:
        return
    for _, r in features.iterrows():
        name = str(r["feature"])
        start = int(r["plot_start"]) if "plot_start" in r and pd.notna(r["plot_start"]) else int(r["start"])
        end = int(r["plot_end"]) if "plot_end" in r and pd.notna(r["plot_end"]) else int(r["end"])
        if end <= start:
            continue
        color = "#e7b3c2" if "IC" in name else "#b8c7d9"
        alpha = 0.09 if "IC" in name else 0.035
        ax.axvspan(start, end, color=color, alpha=alpha, lw=0, zorder=0)


def format_mb_axis(ax) -> None:
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x / 1e6:.1f}"))


def smooth_values(values: pd.Series, window: int = 31) -> pd.Series:
    return pd.to_numeric(values, errors="coerce").rolling(window, center=True, min_periods=max(3, window // 5)).median()


def plot_smoothed_track(
    ax,
    x: pd.Series,
    y: pd.Series,
    color: str,
    label: Optional[str] = None,
    lw: float = 1.4,
    alpha: float = 0.95,
    window: int = 31,
    ls: str = "-",
) -> None:
    yy = smooth_values(y, window=window)
    ax.plot(x, yy, color=color, lw=lw, alpha=alpha, label=label, ls=ls)


def plot_breakpoint_lines(ax, intervals: pd.DataFrame, group: Optional[str], color: str) -> None:
    if intervals.empty:
        return
    z = intervals if group is None else intervals[intervals["group"] == group]
    if z.empty:
        return
    positions = sorted(set(z["start"].astype(int).tolist() + z["end"].astype(int).tolist()))
    for pos in positions:
        ax.axvline(pos, color=color, lw=0.55, ls="--", alpha=0.35, zorder=1)


def unique_breakpoint_positions(intervals: pd.DataFrame, group: str, spec: WindowSpec) -> List[int]:
    if intervals.empty:
        return []
    z = intervals[intervals["group"] == group]
    if z.empty:
        return []
    positions = sorted(set(z["start"].astype(int).tolist() + z["end"].astype(int).tolist()))
    return [p for p in positions if spec.start < p < spec.end]


def bp_hotspots_table(spec: WindowSpec) -> pd.DataFrame:
    rows = []
    for name, pos in BP_HOTSPOTS_T2T:
        rows.append(
            {
                "name": name,
                "chrom": spec.chrom,
                "position": pos,
                "position_mb": pos / 1e6,
                "in_plotted_range": spec.start <= pos <= spec.end,
                "source_script": str(DEFAULT_BP_SOURCE_SCRIPT),
                "source_variable": "PWS_BREAKPOINTS",
            }
        )
    return pd.DataFrame(rows)


def bp_cluster_table(spec: WindowSpec) -> pd.DataFrame:
    rows = []
    for cluster in BP_CLUSTER_INTERVALS_T2T:
        start = int(cluster["start"])
        end = int(cluster["end"])
        rows.append(
            {
                "name": cluster["name"],
                "chrom": spec.chrom,
                "start": start,
                "end": end,
                "start_mb": start / 1e6,
                "end_mb": end / 1e6,
                "approx_size": cluster["approx_size"],
                "interpretation": cluster["interpretation"],
                "overlaps_plotted_range": start < spec.end and end > spec.start,
            }
        )
    return pd.DataFrame(rows)


def snrpn_marker_position(gene_features: pd.DataFrame, spec: WindowSpec) -> int:
    if gene_features is not None and not gene_features.empty:
        for feature in ["SNURF-SNRPN", "PWS-IC"]:
            z = gene_features[gene_features["feature"].astype(str) == feature]
            if not z.empty:
                return int(pd.to_numeric(z.iloc[0]["start"], errors="coerce"))
    return min(max(22_560_000, spec.start), spec.end)


def add_panel_vertical_annotations(
    ax,
    deletion_intervals: pd.DataFrame,
    snrpn_pos: int,
    spec: WindowSpec,
    include_labels: bool = True,
) -> None:
    for cluster in BP_CLUSTER_INTERVALS_T2T:
        cluster_start = int(cluster["start"])
        cluster_end = int(cluster["end"])
        plot_start = max(cluster_start, spec.start)
        plot_end = min(cluster_end, spec.end)
        if plot_start >= plot_end:
            continue
        ax.axvspan(plot_start, plot_end, color="#8B6914", alpha=0.075, lw=0, zorder=0)
        ax.axvline(plot_start, color="#8B6914", lw=0.55, ls=":", alpha=0.55, zorder=1)
        ax.axvline(plot_end, color="#8B6914", lw=0.55, ls=":", alpha=0.55, zorder=1)
        if include_labels:
            ax.text(
                (plot_start + plot_end) / 2,
                -0.135,
                f"{cluster['name']} SD cluster\n{cluster_start / 1e6:.3f}-{cluster_end / 1e6:.3f} Mb",
                transform=ax.get_xaxis_transform(),
                rotation=0,
                ha="center",
                va="top",
                fontsize=5.25,
                color="#6e520f",
                linespacing=0.9,
                clip_on=False,
                bbox=dict(fc="white", ec="none", alpha=0.78, pad=0.5),
            )

    if spec.start <= snrpn_pos <= spec.end:
        ax.axvline(snrpn_pos, color="#1f1f1f", lw=0.8, ls="-", alpha=0.58, zorder=3)
        if include_labels:
            ax.text(
                snrpn_pos,
                -0.135,
                f"SNRPN\n{snrpn_pos / 1e6:.3f} Mb",
                transform=ax.get_xaxis_transform(),
                rotation=0,
                ha="center",
                va="top",
                fontsize=5.25,
                color="#1f1f1f",
                linespacing=0.9,
                clip_on=False,
                bbox=dict(fc="white", ec="none", alpha=0.78, pad=0.5),
            )

    for group, color in [
        ("PWS_DEL", MATERNAL),
        ("AS_DEL", PATERNAL),
    ]:
        positions = unique_breakpoint_positions(deletion_intervals, group, spec)
        if not positions:
            continue
        for pos in positions:
            ax.axvline(pos, color=color, lw=0.6, ls="--", alpha=0.34, zorder=2)


def legend_outside(ax, include_annotations: bool = True, **kwargs) -> None:
    handles, labels = ax.get_legend_handles_labels()
    if include_annotations:
        extra = [
            Line2D([0], [0], color="#1f1f1f", lw=0.8, ls="-", alpha=0.72, label="SNRPN"),
            Line2D([0], [0], color=MATERNAL, lw=0.8, ls="--", alpha=0.72, label="PWS-DEL breakpoints"),
            Line2D([0], [0], color=PATERNAL, lw=0.8, ls="--", alpha=0.72, label="AS-DEL breakpoints"),
            Rectangle((0, 0), 1, 1, facecolor="#8B6914", alpha=0.14, edgecolor="#8B6914", label="BP1/BP2/BP3 T2T SD clusters"),
        ]
        existing = set(labels)
        for h in extra:
            label = h.get_label()
            if label not in existing:
                handles.append(h)
                labels.append(label)
                existing.add(label)
    if not handles:
        return
    defaults = {
        "loc": "upper left",
        "bbox_to_anchor": (1.12, 1.0),
        "borderaxespad": 0.0,
        "frameon": False,
        "fontsize": 6.4,
        "handlelength": 2.1,
    }
    defaults.update(kwargs)
    ax.legend(handles, labels, **defaults)


def add_correlation_box(ax, text: str, y: float = 0.02) -> None:
    if not text:
        return
    ax.text(
        1.12,
        0.02 if y is None else y,
        text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=6.5,
        color="#222222",
        bbox=dict(fc="white", ec="#cfcfcf", lw=0.35, alpha=0.82, pad=2.0),
        clip_on=False,
    )


def correlation_label(a: pd.Series, b: pd.Series, prefix: str) -> str:
    stats = pearson_summary(a, b)
    r = stats["pearson_r"]
    n = stats["n_windows"]
    if not np.isfinite(r):
        return f"{prefix}\nr=NA, n={n}"
    return f"{prefix}\nr={r:.3f}, n={n}"


def retained_group_correlation_label(
    retained: pd.DataFrame,
    control_ref: pd.DataFrame,
    group: str,
    deletion_intervals: pd.DataFrame,
) -> str:
    if retained.empty or control_ref.empty:
        return ""
    ref_col = "control_maternal_mean" if group == "PWS_DEL" else "control_paternal_mean"
    label = "retained vs control maternal" if group == "PWS_DEL" else "retained vs control paternal"
    values = []
    for _, z in retained[retained["group"] == group].groupby("sample"):
        zz = z.merge(control_ref[["window_start", ref_col]], on="window_start", how="left")
        zz = zz[zz["within_deletion_footprint"]]
        stats = pearson_summary(zz["retained_methylation"], zz[ref_col])
        if np.isfinite(stats["pearson_r"]):
            values.append(stats["pearson_r"])
    if not values:
        return ""
    return f"Pearson {label}\nmean r={np.mean(values):.3f}\nrange {min(values):.3f}-{max(values):.3f}"


def mupd_correlation_label(all_windows: pd.DataFrame, control_ref: pd.DataFrame) -> str:
    if all_windows.empty or control_ref.empty:
        return ""
    upd = all_windows[all_windows["sample"] == "004P"].copy()
    if upd.empty:
        return ""
    zz = upd.merge(control_ref[["window_start", "control_maternal_mean"]], on="window_start", how="left")
    parts = []
    for hap, col in [("hap1", "mean_meth_haplotype1"), ("hap2", "mean_meth_haplotype2")]:
        stats = pearson_summary(zz[col], zz["control_maternal_mean"])
        if np.isfinite(stats["pearson_r"]):
            parts.append(f"{hap} r={stats['pearson_r']:.3f}")
    return "Pearson vs control maternal\n" + "\n".join(parts) if parts else ""


def reciprocal_overlay_correlation_label(retained: pd.DataFrame) -> str:
    if retained.empty:
        return ""
    pws = (
        retained[retained["group"] == "PWS_DEL"]
        .groupby("window_start")["retained_methylation"]
        .mean()
        .reset_index(name="pws_maternal")
    )
    asdel = (
        retained[retained["group"] == "AS_DEL"]
        .groupby("window_start")["retained_methylation"]
        .mean()
        .reset_index(name="asdel_paternal")
    )
    z = pws.merge(asdel, on="window_start", how="inner")
    return correlation_label(z["pws_maternal"], z["asdel_paternal"], "Pearson PWS maternal vs AS paternal")


def plot_gene_track(
    ax,
    features: pd.DataFrame,
    spec: WindowSpec,
    icrs: Optional[pd.DataFrame] = None,
    segdup_blocks: Optional[pd.DataFrame] = None,
    deletion_intervals: Optional[pd.DataFrame] = None,
) -> None:
    ax.set_xlim(spec.start, spec.end)
    ax.set_ylim(0, 1.08)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="#eeeeee", lw=0.4)
    for y in [0.79, 0.55, 0.34, 0.13]:
        ax.hlines(y, spec.start, spec.end, color="#d6d6d6", lw=0.6)
    for label, y in [("Genes", 0.84), ("ICRs", 0.60), ("Segmental duplications", 0.38), ("CNV breakpoints", 0.18)]:
        ax.text(-0.012, y, label, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7.0, fontweight="bold")

    gene_lanes = {
        "UBE3A-ATS": 0.88,
        "SNURF-SNRPN": 0.72,
        "SNORD116": 0.80,
        "SNORD115": 0.72,
        "UBE3A": 0.80,
    }
    default_levels = [0.73, 0.81, 0.89]
    upstream_gene_labels = []
    regulatory_features = []
    for i, (_, r) in enumerate(features.iterrows()):
        if not bool(r.get("within_phase2_window", True)):
            upstream_gene_labels.append(str(r["feature"]))
            continue
        if "IC" in str(r["feature"]):
            regulatory_features.append(r)
            continue
        y = gene_lanes.get(str(r["feature"]), default_levels[i % len(default_levels)])
        x0 = int(r["plot_start"]) if "plot_start" in r and pd.notna(r["plot_start"]) else int(r["start"])
        x1 = int(r["plot_end"]) if "plot_end" in r and pd.notna(r["plot_end"]) else int(r["end"])
        width = max(1, x1 - x0)
        ax.add_patch(Rectangle((x0, y - 0.020), width, 0.040, fc="#b8c7d9", ec="#333333", lw=0.3))
        ax.text(
            (x0 + x1) / 2,
            y + 0.030,
            str(r["feature"]),
            ha="center",
            va="bottom",
            fontsize=5.2,
            rotation=0,
            bbox=dict(fc="white", ec="none", alpha=0.72, pad=0.5),
        )

    upstream_icr_labels = []
    for i, r in enumerate(regulatory_features):
        x0 = int(r["plot_start"]) if "plot_start" in r and pd.notna(r["plot_start"]) else int(r["start"])
        x1 = int(r["plot_end"]) if "plot_end" in r and pd.notna(r["plot_end"]) else int(r["end"])
        y = 0.60 if str(r["feature"]) == "AS-IC" else 0.52
        ax.add_patch(Rectangle((x0, y - 0.026), max(1, x1 - x0), 0.052, fc="#f2b8c5", ec="#7d3a44", lw=0.3))
        ax.text(
            (x0 + x1) / 2,
            y + 0.038,
            str(r["feature"]),
            ha="center",
            va="bottom",
            fontsize=5.4,
            rotation=28,
            bbox=dict(fc="white", ec="none", alpha=0.74, pad=0.5),
        )

    if icrs is not None and not icrs.empty:
        for i, (_, r) in enumerate(icrs.iterrows()):
            if not bool(r.get("within_phase2_window", True)):
                upstream_icr_labels.append(f"{r['name']} {r['gene']}".strip())
                continue
            x0 = int(r["plot_start"])
            x1 = int(r["plot_end"])
            y = 0.49 + (i % 2) * 0.060
            ax.add_patch(Rectangle((x0, y - 0.025), max(1, x1 - x0), 0.05, fc="#e89aaa", ec="#7d3a44", lw=0.3))
            label = re.sub(r"[*^#]", "", str(r["name"]))
            ax.text(
                (x0 + x1) / 2,
                y + 0.035,
                label,
                ha="center",
                va="bottom",
                fontsize=5.0,
                rotation=28,
                bbox=dict(fc="white", ec="none", alpha=0.72, pad=0.5),
            )

    upstream_bits = []
    if upstream_gene_labels:
        upstream_bits.append("upstream genes: " + ", ".join(upstream_gene_labels))
    if upstream_icr_labels:
        upstream_bits.append("upstream ICRs: " + ", ".join(upstream_icr_labels))
    if upstream_bits:
        ax.text(
            spec.start + 0.004 * (spec.end - spec.start),
            0.985,
            "Outside fixed x-range: " + "; ".join(upstream_bits),
            ha="left",
            va="top",
            fontsize=5.6,
            color="#333333",
            bbox=dict(fc="white", ec="none", alpha=0.80, pad=0.6),
        )

    if segdup_blocks is not None and not segdup_blocks.empty:
        for _, r in segdup_blocks.iterrows():
            ax.add_patch(
                Rectangle(
                    (int(r["start"]), 0.275),
                    max(1, int(r["end"]) - int(r["start"])),
                    0.05,
                    fc="#8f8f8f",
                    ec="none",
                    alpha=0.45,
                )
            )

    if deletion_intervals is not None and not deletion_intervals.empty:
        colors = {"PWS_DEL": MATERNAL, "AS_DEL": PATERNAL}
        for group, z in deletion_intervals.groupby("group"):
            y = 0.09 if group == "PWS_DEL" else 0.15
            for _, r in z.iterrows():
                for pos in [int(r["start"]), int(r["end"])]:
                    ax.vlines(pos, y - 0.035, y + 0.035, color=colors.get(group, GREY), lw=0.65, alpha=0.75)
        ax.text(spec.end, 0.09, "PWS-DEL", ha="right", va="center", fontsize=5.6, color=MATERNAL)
        ax.text(spec.end, 0.15, "AS-DEL", ha="right", va="center", fontsize=5.6, color=PATERNAL)

    ax.text(
        spec.start,
        1.06,
        "T2T annotation tracks: genes, ICRs, segmental duplications, deletion breakpoints",
        ha="left",
        va="top",
        fontsize=8.2,
        fontweight="bold",
    )


def plot_deletion_shading(ax, intervals: pd.DataFrame, group: str) -> None:
    z = intervals[intervals["group"] == group] if not intervals.empty else intervals
    if z.empty:
        return
    for _, r in z.iterrows():
        ax.axvspan(int(r["start"]), int(r["end"]), color="#bdbdbd", alpha=0.10, lw=0)


def plot_phase2_figure(
    outdir: Path,
    spec: WindowSpec,
    all_windows: pd.DataFrame,
    control_ref: pd.DataFrame,
    control_samples: pd.DataFrame,
    retained: pd.DataFrame,
    deletion_intervals: pd.DataFrame,
    regions: pd.DataFrame,
    gene_features: pd.DataFrame,
    reciprocal_delta: Optional[pd.DataFrame] = None,
    reciprocal_boundaries: Optional[pd.DataFrame] = None,
    icrs: Optional[pd.DataFrame] = None,
    segdup_blocks: Optional[pd.DataFrame] = None,
) -> None:
    fig = plt.figure(figsize=(19, 15.8))
    gs = GridSpec(
        5,
        1,
        height_ratios=[1.12, 1.0, 1.0, 0.92, 1.0],
        hspace=0.56,
        figure=fig,
    )
    axes = [fig.add_subplot(gs[i, 0]) for i in range(5)]
    snrpn_pos = snrpn_marker_position(gene_features, spec)

    for ax in axes:
        ax.set_xlim(spec.start, spec.end)
        ax.set_ylim(-0.05, 1.08)
        ax.set_ylabel("Meth.")
        ax.grid(axis="y", color="#ececec", lw=0.5)
        ax.tick_params(axis="x", labelbottom=False)
        ax.spines[["top", "right"]].set_visible(False)
        add_panel_vertical_annotations(ax, deletion_intervals, snrpn_pos, spec)
        format_mb_axis(ax)

    # Panel A: controls, separate lines per haplotype and sample.
    ax = axes[0]
    add_region_spans(ax, regions)
    if not control_samples.empty:
        for sample, z in control_samples.groupby("sample"):
            z = z.sort_values("window_start")
            plot_smoothed_track(ax, z["window_mid"], z["control_maternal"], color=MATERNAL, lw=0.75, alpha=0.28)
            plot_smoothed_track(ax, z["window_mid"], z["control_paternal"], color=PATERNAL, lw=0.75, alpha=0.28)
    if not control_ref.empty:
        plot_smoothed_track(ax, control_ref["window_mid"], control_ref["control_maternal_mean"], color=MATERNAL, lw=2.0, alpha=0.98, label="maternal mean")
        plot_smoothed_track(ax, control_ref["window_mid"], control_ref["control_paternal_mean"], color=PATERNAL, lw=2.0, alpha=0.98, label="paternal mean")
        add_correlation_box(
            ax,
            correlation_label(
                control_ref["control_maternal_mean"],
                control_ref["control_paternal_mean"],
                "Pearson maternal vs paternal",
            ),
        )
    ax.set_title("A. Controls: reciprocal maternal and paternal cis-methylation architecture", loc="left", fontsize=9, fontweight="bold", pad=7)
    legend_outside(ax, fontsize=6.4)

    # Panel B: PWS-DEL retained maternal signal plus coverage track.
    ax = axes[1]
    plot_deletion_shading(ax, deletion_intervals, "PWS_DEL")
    pws = retained[retained["group"] == "PWS_DEL"]
    for sample, z in pws.groupby("sample"):
        z = z.sort_values("window_start")
        plot_smoothed_track(ax, z["window_mid"], z["retained_methylation"], color=MATERNAL, lw=0.65, alpha=0.18)
    if not control_ref.empty:
        plot_smoothed_track(ax, control_ref["window_mid"], control_ref["control_maternal_mean"], color=MATERNAL, lw=1.2, ls=":", alpha=0.85, label="control maternal mean")
    if not pws.empty:
        mean = pws.groupby("window_mid")["retained_methylation"].mean().reset_index()
        plot_smoothed_track(ax, mean["window_mid"], mean["retained_methylation"], color=MATERNAL, lw=2.0, alpha=0.98, label="PWS-DEL retained mean")
    cov_ax = ax.twinx()
    if not pws.empty:
        cov = pws.groupby("window_mid")["prompt_haplotype2_coverage"].mean().reset_index()
        cov_ax.plot(cov["window_mid"], cov["prompt_haplotype2_coverage"], color=GREY, lw=0.65, alpha=0.60)
    cov_ax.set_ylabel("hap2 reads", fontsize=8, color=GREY)
    cov_ax.tick_params(axis="y", labelsize=7, colors=GREY)
    cov_ax.spines["top"].set_visible(False)
    add_correlation_box(ax, retained_group_correlation_label(retained, control_ref, "PWS_DEL", deletion_intervals))
    ax.set_title("B. PWS-DEL: retained maternal-only architecture across the deletion footprint", loc="left", fontsize=9, fontweight="bold", pad=7)
    legend_outside(ax, fontsize=6.4)

    # Panel C: AS-DEL retained paternal signal plus absent-haplotype coverage.
    ax = axes[2]
    plot_deletion_shading(ax, deletion_intervals, "AS_DEL")
    asdel = retained[retained["group"] == "AS_DEL"]
    for sample, z in asdel.groupby("sample"):
        z = z.sort_values("window_start")
        plot_smoothed_track(ax, z["window_mid"], z["retained_methylation"], color=PATERNAL, lw=0.65, alpha=0.18)
    if not control_ref.empty:
        plot_smoothed_track(ax, control_ref["window_mid"], control_ref["control_paternal_mean"], color=PATERNAL, lw=1.2, ls=":", alpha=0.85, label="control paternal mean")
    if not asdel.empty:
        mean = asdel.groupby("window_mid")["retained_methylation"].mean().reset_index()
        plot_smoothed_track(ax, mean["window_mid"], mean["retained_methylation"], color=PATERNAL, lw=2.0, alpha=0.98, label="AS-DEL retained mean")
    cov_ax = ax.twinx()
    if not asdel.empty:
        cov = asdel.groupby("window_mid")["absent_haplotype_coverage"].mean().reset_index()
        cov_ax.plot(cov["window_mid"], cov["absent_haplotype_coverage"], color=GREY, lw=0.65, alpha=0.60)
    cov_ax.set_ylabel("absent-hap reads", fontsize=8, color=GREY)
    cov_ax.tick_params(axis="y", labelsize=7, colors=GREY)
    cov_ax.spines["top"].set_visible(False)
    add_correlation_box(ax, retained_group_correlation_label(retained, control_ref, "AS_DEL", deletion_intervals))
    ax.set_title("C. AS-DEL: retained paternal-only architecture across the deletion footprint", loc="left", fontsize=9, fontweight="bold", pad=7)
    legend_outside(ax, fontsize=6.4)

    # Panel D: PWS-mUPD both haplotypes.
    ax = axes[3]
    upd = all_windows[all_windows["sample"] == "004P"].sort_values("window_start")
    if not upd.empty:
        plot_smoothed_track(ax, upd["window_mid"], upd["mean_meth_haplotype1"], color=MUPD_1, lw=1.35, label="004P haplotype 1")
        plot_smoothed_track(ax, upd["window_mid"], upd["mean_meth_haplotype2"], color=MUPD_2, lw=1.35, label="004P haplotype 2")
    if not control_ref.empty:
        plot_smoothed_track(ax, control_ref["window_mid"], control_ref["control_maternal_mean"], color=MATERNAL, lw=1.1, ls=":", label="control maternal mean")
    add_correlation_box(ax, mupd_correlation_label(all_windows, control_ref))
    ax.set_title("D. PWS-mUPD: both haplotypes test the duplicated maternal architecture", loc="left", fontsize=9, fontweight="bold", pad=7)
    legend_outside(ax, fontsize=6.4)

    # Panel E: reciprocal natural dissection overlay.
    ax = axes[4]
    plot_deletion_shading(ax, deletion_intervals, "PWS_DEL")
    plot_deletion_shading(ax, deletion_intervals, "AS_DEL")
    overlay_rows = []
    for group, color, label in [
        ("PWS_DEL", MATERNAL, "PWS-DEL retained maternal"),
        ("AS_DEL", PATERNAL, "AS-DEL retained paternal"),
    ]:
        z = retained[retained["group"] == group]
        if z.empty:
            continue
        mean = z.groupby("window_mid")["retained_methylation"].mean().reset_index()
        plot_smoothed_track(ax, mean["window_mid"], mean["retained_methylation"], color=color, lw=2.0, label=label)
    ax.set_title("E. Reciprocal natural dissection: maternal-only and paternal-only signals overlaid", loc="left", fontsize=9, fontweight="bold", pad=7)
    add_correlation_box(ax, reciprocal_overlay_correlation_label(retained))
    legend_outside(ax, fontsize=6.4)
    ax.tick_params(axis="x", labelbottom=True)
    ax.set_xlabel(f"{spec.chrom} coordinate (Mb)", labelpad=34)

    fig.suptitle("Figure 2. Reciprocal cis-architecture across the chr15 imprinted domain", y=0.986, fontsize=13, fontweight="bold")
    fig.text(
        0.07,
        0.955,
        "T2T breakpoint SD clusters shaded in all panels: BP1=17.691-20.454 Mb, BP2=20.754-21.184 Mb, BP3=25.876-26.633 Mb.",
        ha="left",
        va="top",
        fontsize=7.2,
        color="#333333",
    )
    fig.subplots_adjust(top=0.925, bottom=0.125, left=0.07, right=0.70)
    fig.savefig(outdir / "figures" / "Figure2_reciprocal_cis_architecture.png", dpi=350, bbox_inches="tight")
    fig.savefig(outdir / "figures" / "Figure2_reciprocal_cis_architecture.pdf", bbox_inches="tight")
    plt.close(fig)


def _phase2_main() -> None:
    parser = argparse.ArgumentParser(description="Phase 2 reciprocal cis-architecture analysis")
    parser.add_argument("--methylation-dir", "--meth-dir", dest="meth_dir", default=DEFAULT_METHYLATION_DIR, type=Path)
    parser.add_argument("--metadata", default=DEFAULT_METADATA_PATH, type=Path)
    parser.add_argument("--gtf", default=DEFAULT_GTF, type=Path)
    parser.add_argument("--cnv-dir", default=DEFAULT_CNV_DIR, type=Path)
    parser.add_argument("--icr-bed", default=DEFAULT_ICR_BED, type=Path, help="T2T imprinting-control-region BED annotations.")
    parser.add_argument("--segdup-bed", default=DEFAULT_SEG_DUP, type=Path, help="T2T segmental-duplication BED-like annotation table.")
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR, type=Path)
    parser.add_argument("--chrom", default=PHASE2_CHROM)
    parser.add_argument("--start", default=PHASE2_START, type=int)
    parser.add_argument("--end", default=PHASE2_END, type=int)
    parser.add_argument("--window-size", default=PHASE2_WINDOW, type=int)
    parser.add_argument("--samples", default=",".join(DEFAULT_PHASE2_SAMPLES),
                        help="Comma-separated sample list for Phase 2; defaults to PWS-DEL, PWS-mUPD, AS-DEL, and controls only.")
    parser.add_argument("--exclude", default=",".join(sorted(DEFAULT_EXCLUDE)))
    parser.add_argument("--methyl-format", default="auto", choices=["auto", "bed", "dss"])
    parser.add_argument("--bed-meth-col", default="9", help="1-based methylation column; pb-CpG BED default is 9.")
    parser.add_argument("--bed-cov-col", default="6", help="1-based read coverage column; pb-CpG BED default is 6.")
    parser.add_argument("--chunksize", default=1_000_000, type=int)
    parser.add_argument("--min-reads", default=10, type=int)
    parser.add_argument("--min-cpgs", default=3, type=int)
    parser.add_argument("--deletion-copy-threshold", default=1.5, type=float)
    parser.add_argument("--divergent-delta", default=0.20, type=float)
    parser.add_argument("--convergent-delta", default=0.10, type=float)
    parser.add_argument("--min-region-bp", default=25_000, type=int)
    parser.add_argument("--deviation-threshold", default=0.25, type=float)
    parser.add_argument("--reciprocal-delta-threshold", default=0.10, type=float)
    parser.add_argument("--reciprocal-smooth-windows", default=31, type=int)
    parser.add_argument("--reciprocal-min-boundary-bp", default=10_000, type=int)
    args = parser.parse_args()

    outdir = mkdir(args.outdir)
    tbldir = mkdir(outdir / "tables")
    figdir = mkdir(outdir / "figures")
    spec = WindowSpec(args.chrom, args.start, args.end, args.window_size)

    exclude = {x.strip() for x in args.exclude.split(",") if x.strip()}

    log("Loading metadata and methylation file inventory")
    load_metadata_table(args.metadata, outdir=outdir)
    file_table = discover_methylation_files(args.meth_dir, exclude)
    file_table["group"] = file_table["sample"].apply(lambda s: sample_group(s, DEFAULT_GROUPS))
    requested_samples = [x.strip() for x in args.samples.split(",") if x.strip()]
    samples = phase_samples_from_metadata(file_table, requested=requested_samples)
    file_table = file_table[file_table["sample"].isin(samples)].copy()
    file_table.to_csv(tbldir / "Phase2_input_methylation_files.tsv", sep="\t", index=False)

    log("Loading T2T gene annotations and IC proxy")
    genes = load_gtf_genes(args.gtf, args.chrom)
    ic_chrom, ic_start, ic_end, ic_gene, _ = infer_ic_region_from_gtf(genes, IC_PRIORITY_GENES, flank=5_000)
    if ic_chrom != args.chrom:
        raise ValueError(f"IC inferred on {ic_chrom}, but Phase 2 is running on {args.chrom}")

    gene_features = build_phase2_gene_features(genes, ic_start, ic_end, spec)
    gene_features.to_csv(tbldir / "Phase2_gene_track_features.tsv", sep="\t", index=False)
    icr_annotations = read_icr_annotations(args.icr_bed, spec)
    if not icr_annotations.empty:
        icr_annotations.to_csv(tbldir / "Phase2_T2T_ICR_annotations.tsv", sep="\t", index=False)
    segdup_raw, segdup_blocks = read_segdup_annotations(args.segdup_bed, spec)
    if not segdup_raw.empty:
        segdup_raw.to_csv(tbldir / "Phase2_T2T_segmental_duplications_raw.tsv.gz", sep="\t", index=False, compression="gzip")
    if not segdup_blocks.empty:
        segdup_blocks.to_csv(tbldir / "Phase2_T2T_segmental_duplication_blocks.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "feature": "phase2_fixed_domain",
                "chrom": args.chrom,
                "start": args.start,
                "end": args.end,
                "window_size": args.window_size,
            },
            {
                "feature": f"PWS_IC_proxy_{ic_gene}",
                "chrom": ic_chrom,
                "start": ic_start,
                "end": ic_end,
                "window_size": np.nan,
            },
            {
                "feature": "AS_IC_prompt_approx",
                "chrom": args.chrom,
                "start": AS_IC_CENTER - 10_000,
                "end": AS_IC_CENTER + 10_000,
                "window_size": np.nan,
            },
        ]
    ).to_csv(tbldir / "Phase2_fixed_coordinates.tsv", sep="\t", index=False)
    bp_hotspots_table(spec).to_csv(tbldir / "Phase2_T2T_BP_hotspots.tsv", sep="\t", index=False)
    bp_cluster_table(spec).to_csv(tbldir / "Phase2_T2T_BP_cluster_annotations.tsv", sep="\t", index=False)

    log("Assigning parental labels from IC methylation and group expectations")
    ic_summary = summarize_ic_for_labels(file_table, args.chrom, ic_start, ic_end, args)
    ic_summary.to_csv(tbldir / "Phase2_IC_haplotype_methylation_for_labeling.tsv", sep="\t", index=False)
    labels = labels_from_ic_and_group(ic_summary, min_cpgs=args.min_cpgs)

    log("Computing 1 kb per-haplotype windows")
    window_tables = []
    for i, sample in enumerate(samples, start=1):
        log(f"Windowing sample {i}/{len(samples)}: {sample}")
        window_tables.append(summarize_sample_windows(sample, file_table, spec, labels, args))
    all_windows = pd.concat(window_tables, ignore_index=True) if window_tables else pd.DataFrame()

    log("Loading CNV deletion footprints")
    cnv_table = discover_cnv_files(args.cnv_dir, samples)
    cnv_table.to_csv(tbldir / "Phase2_input_cnv_files.tsv", sep="\t", index=False)
    deletion_intervals = read_cnv_deletion_intervals(cnv_table, spec, args.deletion_copy_threshold)
    deletion_intervals = deletion_intervals[deletion_intervals["group"].isin(["PWS_DEL", "AS_DEL"])].copy()
    if deletion_intervals.empty:
        log("No CNV deletion intervals found in fixed domain; using haplotype coverage fallback.")
        deletion_intervals = infer_deletion_intervals_from_windows(all_windows, spec, min_reads=args.min_reads)
    deletion_intervals.to_csv(tbldir / "Phase2_deletion_footprints.tsv", sep="\t", index=False)
    breakpoint_annotations = annotate_deletion_breakpoints(deletion_intervals, genes, icr_annotations, segdup_blocks)
    if not breakpoint_annotations.empty:
        breakpoint_annotations.to_csv(tbldir / "Phase2_deletion_breakpoints_T2T_annotated.tsv", sep="\t", index=False)

    labels = update_del_labels_by_coverage(labels, all_windows, deletion_intervals)
    all_windows = apply_labels_to_windows(all_windows, labels)
    write_per_sample_tables(all_windows, outdir)
    all_windows.to_csv(
        tbldir / "Phase2_all_samples_haplotype_1kb_methylation.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
        na_rep="NaN",
    )
    pd.DataFrame(
        [
            {"sample": sample, "haplotype1_label": lab.get("hap1"), "haplotype2_label": lab.get("hap2")}
            for sample, lab in sorted(labels.items())
            if sample in samples
        ]
    ).to_csv(tbldir / "Phase2_haplotype_parental_labels.tsv", sep="\t", index=False)

    log("Building control reference and retained-haplotype profiles")
    ref, control_samples = control_reference(all_windows)
    ref.to_csv(tbldir / "Phase2_control_reference_architecture.tsv", sep="\t", index=False, na_rep="NaN")
    control_samples.to_csv(tbldir / "Phase2_control_sample_architecture.tsv", sep="\t", index=False, na_rep="NaN")
    regions = classify_control_regions(ref, args.divergent_delta, args.convergent_delta, args.min_region_bp)
    regions.to_csv(tbldir / "Phase2_control_divergence_convergence_regions.tsv", sep="\t", index=False)

    retained = retained_profiles(all_windows, deletion_intervals)
    retained.to_csv(tbldir / "Phase2_retained_haplotype_profiles.tsv", sep="\t", index=False, na_rep="NaN")
    corr, dev = compute_correlations(retained, all_windows, ref, deletion_intervals, args.deviation_threshold)
    corr.to_csv(tbldir / "Phase2_retained_and_mUPD_correlations.tsv", sep="\t", index=False, na_rep="NaN")
    dev.to_csv(tbldir / "Phase2_004P_candidate_postzygotic_deviation_windows.tsv", sep="\t", index=False, na_rep="NaN")

    reciprocal_delta, reciprocal_boundaries = reciprocal_delta_and_boundaries(
        retained,
        spec,
        genes,
        icr_annotations,
        segdup_blocks,
        threshold=args.reciprocal_delta_threshold,
        smooth_window=args.reciprocal_smooth_windows,
        min_bp=args.reciprocal_min_boundary_bp,
    )
    reciprocal_delta.to_csv(
        tbldir / "Phase2_reciprocal_delta_profile.tsv",
        sep="\t",
        index=False,
        na_rep="NaN",
    )
    reciprocal_boundaries.to_csv(
        tbldir / "Phase2_reciprocal_delta_boundary_candidates.tsv",
        sep="\t",
        index=False,
        na_rep="NaN",
    )

    log("Rendering Figure 2 composite")
    plot_phase2_figure(
        outdir,
        spec,
        all_windows,
        ref,
        control_samples,
        retained,
        deletion_intervals,
        regions,
        gene_features,
        reciprocal_delta=reciprocal_delta,
        reciprocal_boundaries=reciprocal_boundaries,
        icrs=icr_annotations,
        segdup_blocks=segdup_blocks,
    )

    pd.DataFrame(
        [
            {
                "phase": "Phase2_reciprocal_cis_architecture",
                "domain": f"{args.chrom}:{args.start}-{args.end}",
                "window_size": args.window_size,
                "n_samples": len(samples),
                "n_windows_per_sample": int((args.end - args.start) / args.window_size),
                "min_reads_per_haplotype": args.min_reads,
                "min_cpgs_per_haplotype": args.min_cpgs,
                "figure": str(outdir / "figures" / "Figure2_reciprocal_cis_architecture.png"),
                "all_windows_table": str(tbldir / "Phase2_all_samples_haplotype_1kb_methylation.tsv.gz"),
                "reciprocal_delta_profile_table": str(tbldir / "Phase2_reciprocal_delta_profile.tsv"),
                "reciprocal_boundary_candidate_table": str(tbldir / "Phase2_reciprocal_delta_boundary_candidates.tsv"),
            }
        ]
    ).to_csv(outdir / "phase2_run_summary.tsv", sep="\t", index=False)
    log("Done.")

# --- Vendored from scripts/paper_vf/create_figure2_reciprocal_cis_architecture_improved.py ---

"""
Create an improved Figure 2 for the Q1 genomics paper from existing Phase 2 tables.

This second-version figure keeps the original Phase 2 computations intact and
re-renders the interpretation as a cleaner reciprocal parental-architecture
figure with:
  - a schematic logic panel,
  - a shared T2T annotation track,
  - control, PWS-DEL, AS-DEL, and PWS-mUPD architecture panels,
  - a delta-architecture panel,
  - panel-specific plotting tables,
  - and a Markdown report.
"""


import argparse
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import re
import sys
import textwrap
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.ticker import FuncFormatter


DEFAULT_OUTDIR = Path("/home/rare/arlen/pws-as-hifi-cis-methylation/scripts/hifi_multiomics_pipeline/06_results")
DEFAULT_GTF = Path("/home/rare/arlen/reference/chm13v22.sorted.gtf")
DEFAULT_ICR_BED = Path("/home/rare/arlen/reference/ICR_t2t.bed")
DISPLAY_START = 18_000_000
DISPLAY_END = 28_000_000
DISPLAY_WINDOW = 1_000
STRUCT_GREY = "#8a8a8a"
DARK_GREY = "#3a3a3a"
SOFT_GREY = "#efefef"
BP_REGION_COLORS = {
    "BP1": "#c6a96e",
    "BP2": "#8faa92",
    "BP3": "#92a8c9",
}

ANNOTATION_GENE_SPECS = [
    ("NIPA1", "NIPA1", "top"),
    ("NIPA2", "NIPA2", "bottom"),
    ("CYFIP1", "CYFIP1", "top"),
    ("TUBGCP5", "TUBGCP5", "bottom"),
    ("MKRN3", "MKRN3", "top"),
    ("MAGEL2", "MAGEL2", "bottom"),
    ("NDN", "NDN", "top"),
    ("NPAP1", "NPAP1", "bottom"),
    ("SNORD116", "SNORD116", "top"),
    ("SNORD115", "SNORD115", "bottom"),
    ("UBE3A-ATS", "SNHG14", "top"),
    ("UBE3A", "UBE3A", "bottom"),
    ("ATP10A", "ATP10A", "top"),
    ("GABRB3", "GABRB3", "bottom"),
    ("GABRA5", "GABRA5", "top"),
    ("GABRG3", "GABRG3", "bottom"),
    ("OCA2", "OCA2", "top"),
    ("HERC2", "HERC2", "bottom"),
    ("APBA2", "APBA2", "top"),
    ("TJP1", "TJP1", "bottom"),
]


@dataclass(frozen=True)
class PanelStats:
    title: str
    stats_box: str
    badge: str


def log(message: str) -> None:
    print(f"[Figure2-improved] {message}", file=sys.stderr, flush=True)


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


def parse_gtf_attrs(attrs: str) -> Dict[str, str]:
    return {m.group(1): m.group(2) for m in re.finditer(r'(\S+) "([^"]*)";', attrs)}


def mb_formatter(value: float, _pos: float) -> str:
    return f"{value / 1e6:.1f}"


def overlap_filter(df: pd.DataFrame, start_col: str, end_col: str, spec: WindowSpec) -> pd.DataFrame:
    mask = (pd.to_numeric(df[start_col], errors="coerce") < spec.end) & (
        pd.to_numeric(df[end_col], errors="coerce") > spec.start
    )
    return df.loc[mask].copy()


def clipped_interval(start: float, end: float, spec: WindowSpec) -> Tuple[int, int]:
    return max(int(start), spec.start), min(int(end), spec.end)


def add_smoothed_columns(df: pd.DataFrame, columns: List[str], window: int) -> pd.DataFrame:
    out = df.copy()
    for col in columns:
        out[f"{col}_smoothed"] = smooth_values(out[col], window=window)
    return out


def structural_landmark_color(name: str) -> str:
    return BP_REGION_COLORS.get(name, STRUCT_GREY)


def logic_schematic_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "state": "Control",
                "allele_1": "maternal",
                "allele_2": "paternal",
                "interpretation": "Biparental reference",
            },
            {
                "state": "PWS-DEL",
                "allele_1": "retained_maternal",
                "allele_2": "paternal_deleted",
                "interpretation": "Retained maternal, paternal deleted",
            },
            {
                "state": "AS-DEL",
                "allele_1": "maternal_deleted",
                "allele_2": "retained_paternal",
                "interpretation": "Retained paternal, maternal deleted",
            },
            {
                "state": "PWS-mUPD",
                "allele_1": "maternal_haplotype_1",
                "allele_2": "maternal_haplotype_2",
                "interpretation": "Copy-neutral duplicated maternal state",
            },
        ]
    )


def load_annotation_genes(gtf_path: Path, spec: WindowSpec) -> pd.DataFrame:
    rows = []
    with open(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[0] != spec.chrom or fields[2] != "gene":
                continue
            start = int(fields[3])
            end = int(fields[4])
            if end < spec.start or start > spec.end:
                continue
            attrs = parse_gtf_attrs(fields[8])
            gene = attrs.get("gene") or attrs.get("gene_name") or attrs.get("gene_id", "")
            biotype = attrs.get("gene_biotype", "")
            description = attrs.get("description", "")
            rows.append(
                {
                    "gene": gene,
                    "start": start,
                    "end": end,
                    "strand": fields[6],
                    "gene_biotype": biotype,
                    "description": description,
                }
            )
    genes = pd.DataFrame(rows)
    if genes.empty:
        return genes
    genes = genes[~genes["gene"].str.startswith("LOC")].copy()
    genes = genes[~genes["gene_biotype"].str.contains("pseudogene", case=False, na=False)].copy()
    return genes.sort_values(["start", "end", "gene"]).reset_index(drop=True)


def load_pws_ic_anchor(icr_path: Path, spec: WindowSpec) -> pd.DataFrame:
    rows = []
    with open(icr_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6 or fields[0] != spec.chrom:
                continue
            start = int(fields[1])
            end = int(fields[2])
            if end < spec.start or start > spec.end:
                continue
            rows.append(
                {
                    "chrom": fields[0],
                    "start": start,
                    "end": end,
                    "name": fields[3],
                    "parent": fields[4],
                    "gene": fields[5],
                }
            )
    if not rows:
        return pd.DataFrame()
    icrs = pd.DataFrame(rows)
    icrs = icrs[icrs["gene"].astype(str).str.contains("SNRPN", na=False)].copy()
    return icrs.sort_values(["start", "end"]).reset_index(drop=True)


def select_gene_interval(genes: pd.DataFrame, query: str) -> pd.DataFrame:
    if genes.empty:
        return pd.DataFrame()
    if query in {"SNORD116", "SNORD115"}:
        return genes[genes["gene"].astype(str).str.startswith(f"{query}-")].copy()
    return genes[genes["gene"].astype(str) == query].copy()


def estimate_label_width_bp(label: str) -> int:
    return max(220_000, 24_000 * len(label) + 110_000)


def assign_non_overlapping_label_positions(
    feature_rows: List[Dict[str, object]],
    spec: WindowSpec,
    side: str,
    n_lanes: int = 4,
    gap_bp: int = 80_000,
    center_gap_bp: int = 450_000,
) -> List[Dict[str, object]]:
    if not feature_rows:
        return []
    lane_right = [spec.start - gap_bp] * n_lanes
    lane_last_center = [spec.start - center_gap_bp * 2] * n_lanes
    placed: List[Dict[str, object]] = []
    for row in sorted(feature_rows, key=lambda x: x["center"]):
        width = estimate_label_width_bp(str(row["label"]))
        chosen_lane = None
        chosen_left = None
        for lane in range(n_lanes):
            center_spacing_ok = (float(row["center"]) - lane_last_center[lane]) >= center_gap_bp
            left = max(int(row["center"] - width / 2), int(lane_right[lane] + gap_bp))
            if left + width <= spec.end and center_spacing_ok:
                chosen_lane = lane
                chosen_left = left
                break
        if chosen_lane is None:
            chosen_lane = int(np.argmin(lane_right))
            chosen_left = min(max(int(row["center"] - width / 2), spec.start), spec.end - width)
        lane_right[chosen_lane] = chosen_left + width
        lane_last_center[chosen_lane] = float(row["center"])
        out = dict(row)
        out["label_side"] = side
        out["lane_index"] = chosen_lane
        out["label_x"] = chosen_left + width / 2
        out["label_width_bp"] = width
        placed.append(out)
    return placed


def build_annotation_track_table(
    gene_features: pd.DataFrame,
    gtf_path: Path,
    icr_path: Path,
    spec: WindowSpec,
) -> pd.DataFrame:
    genes = load_annotation_genes(gtf_path, spec)
    pws_ic = load_pws_ic_anchor(icr_path, spec)
    feature_rows = []
    for cluster in BP_CLUSTER_INTERVALS_T2T:
        start, end = clipped_interval(cluster["start"], cluster["end"], spec)
        if start >= end:
            continue
        feature_rows.append(
            {
                "label": cluster["name"],
                "kind": "structural_landmark",
                "start": start,
                "end": end,
                "marker_position": np.nan,
                "lane": "structural",
                "color_role": "grey",
                "center": (start + end) / 2,
                "display_color": structural_landmark_color(cluster["name"]),
            }
        )

    snrpn_pos = snrpn_marker_position(gene_features, spec)
    if not pws_ic.empty:
        icr_row = pws_ic.iloc[0]
        start, end = clipped_interval(icr_row["start"], icr_row["end"], spec)
        if start < end:
            feature_rows.append(
                {
                    "label": "PWS-IC / SNRPN",
                    "kind": "regulatory_anchor",
                    "start": start,
                    "end": end,
                    "marker_position": snrpn_pos,
                    "lane": "anchor",
                    "color_role": "grey",
                    "center": (start + end) / 2,
                }
            )

    top_rows: List[Dict[str, object]] = []
    bottom_rows: List[Dict[str, object]] = []
    for label, query, side in ANNOTATION_GENE_SPECS:
        hit = select_gene_interval(genes, query)
        if hit.empty:
            continue
        start, end = clipped_interval(hit["start"].min(), hit["end"].max(), spec)
        if start >= end:
            continue
        row = {
            "label": label,
            "kind": "gene_context",
            "start": start,
            "end": end,
            "marker_position": np.nan,
            "lane": "gene_track",
            "color_role": "grey",
            "center": (start + end) / 2,
            "query_gene": query,
        }
        if side == "top":
            top_rows.append(row)
        else:
            bottom_rows.append(row)

    feature_rows.extend(assign_non_overlapping_label_positions(top_rows, spec, "top"))
    feature_rows.extend(assign_non_overlapping_label_positions(bottom_rows, spec, "bottom"))
    out = pd.DataFrame(feature_rows)
    if not out.empty:
        out = out.sort_values(["kind", "start", "end", "label"]).reset_index(drop=True)
    return out


def build_panel_tables(
    display_spec: WindowSpec,
    control_ref: pd.DataFrame,
    retained: pd.DataFrame,
    all_windows: pd.DataFrame,
    reciprocal_delta: pd.DataFrame,
    smooth_window: int,
) -> Dict[str, pd.DataFrame]:
    tables: Dict[str, pd.DataFrame] = {}

    panel_b = overlap_filter(control_ref, "window_start", "window_end", display_spec)
    panel_b = add_smoothed_columns(
        panel_b,
        ["control_maternal_mean", "control_paternal_mean"],
        smooth_window,
    )
    tables["Figure2B_control_reference_plot.tsv"] = panel_b

    pws_mean = (
        retained[retained["group"] == "PWS_DEL"]
        .groupby(["window_start", "window_end", "window_mid"], as_index=False)
        .agg(
            pwsdel_retained_maternal_mean=("retained_methylation", "mean"),
            n_pwsdel_samples=("sample", "nunique"),
        )
    )
    panel_c = overlap_filter(
        control_ref.merge(pws_mean, on=["window_start", "window_end", "window_mid"], how="left"),
        "window_start",
        "window_end",
        display_spec,
    )
    panel_c = add_smoothed_columns(
        panel_c,
        ["control_maternal_mean", "pwsdel_retained_maternal_mean"],
        smooth_window,
    )
    tables["Figure2C_pwsdel_retained_maternal_plot.tsv"] = panel_c

    as_mean = (
        retained[retained["group"] == "AS_DEL"]
        .groupby(["window_start", "window_end", "window_mid"], as_index=False)
        .agg(
            asdel_retained_paternal_mean=("retained_methylation", "mean"),
            n_asdel_samples=("sample", "nunique"),
        )
    )
    panel_d = overlap_filter(
        control_ref.merge(as_mean, on=["window_start", "window_end", "window_mid"], how="left"),
        "window_start",
        "window_end",
        display_spec,
    )
    panel_d = add_smoothed_columns(
        panel_d,
        ["control_paternal_mean", "asdel_retained_paternal_mean"],
        smooth_window,
    )
    tables["Figure2D_asdel_retained_paternal_plot.tsv"] = panel_d

    panel_e = overlap_filter(
        all_windows[all_windows["sample"] == "004P"][
            [
                "window_start",
                "window_end",
                "window_mid",
                "mean_meth_haplotype1",
                "mean_meth_haplotype2",
            ]
        ].merge(
            control_ref[["window_start", "control_maternal_mean"]],
            on="window_start",
            how="left",
        ),
        "window_start",
        "window_end",
        display_spec,
    )
    panel_e = add_smoothed_columns(
        panel_e,
        ["mean_meth_haplotype1", "mean_meth_haplotype2", "control_maternal_mean"],
        smooth_window,
    )
    tables["Figure2E_pwsmupd_maternal_plot.tsv"] = panel_e

    panel_f = overlap_filter(reciprocal_delta, "window_start", "window_end", display_spec)
    panel_f = add_smoothed_columns(
        panel_f,
        ["pwsdel_retained_maternal_mean", "asdel_retained_paternal_mean", "reciprocal_delta"],
        smooth_window,
    )
    panel_f["panel_label"] = "shared methylation scaffold"
    panel_f.loc[panel_f["delta_state"] == "paternal_higher", "panel_label"] = "localized_parental_divergence"
    tables["Figure2F_reciprocal_delta_architecture_plot.tsv"] = panel_f
    return tables


def panel_stats(
    control_ref: pd.DataFrame,
    corr: pd.DataFrame,
    reciprocal_delta: pd.DataFrame,
) -> Dict[str, PanelStats]:
    control_stats = pearson_summary(
        control_ref["control_maternal_mean"],
        control_ref["control_paternal_mean"],
    )
    pws = corr[corr["group"] == "PWS_DEL"].copy()
    asdel = corr[corr["group"] == "AS_DEL"].copy()
    upd = corr[corr["group"] == "PWS_mUPD"].copy().sort_values("comparison")
    reciprocal_nonnull = reciprocal_delta.dropna(
        subset=["pwsdel_retained_maternal_mean", "asdel_retained_paternal_mean"]
    )
    reciprocal_stats = pearson_summary(
        reciprocal_nonnull["pwsdel_retained_maternal_mean"],
        reciprocal_nonnull["asdel_retained_paternal_mean"],
    )

    return {
        "A": PanelStats(
            title="A. Controls",
            stats_box=f"r = {control_stats['pearson_r']:.3f}\nn = {control_stats['n_windows']:,} informative windows",
            badge="Biparental reference",
        ),
        "B": PanelStats(
            title="B. PWS-DEL",
            stats_box=(
                f"mean r = {pws['pearson_r'].mean():.3f}\n"
                f"range = {pws['pearson_r'].min():.3f}-{pws['pearson_r'].max():.3f}"
            ),
            badge="Maternal-only architecture exposed by paternal deletion",
        ),
        "C": PanelStats(
            title="C. AS-DEL",
            stats_box=(
                f"mean r = {asdel['pearson_r'].mean():.3f}\n"
                f"range = {asdel['pearson_r'].min():.3f}-{asdel['pearson_r'].max():.3f}"
            ),
            badge="Paternal-only architecture exposed by maternal deletion",
        ),
        "D": PanelStats(
            title="D. PWS-mUPD",
            stats_box=(
                f"hap1 r = {upd.iloc[0]['pearson_r']:.3f}\n"
                f"hap2 r = {upd.iloc[1]['pearson_r']:.3f}"
            ),
            badge="Copy-neutral duplicated maternal architecture",
        ),
        "E": PanelStats(
            title="E. Delta architecture",
            stats_box=(
                f"r = {reciprocal_stats['pearson_r']:.3f}\n"
                f"n = {reciprocal_stats['n_windows']:,} shared windows"
            ),
            badge="Shared methylation scaffold",
        ),
    }


def add_shared_axis_style(ax, spec: WindowSpec, ylabel: str, show_x: bool) -> None:
    ax.set_xlim(spec.start, spec.end)
    ax.grid(axis="y", color="#e8e8e8", lw=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.tick_params(axis="both", labelsize=10, width=0.9)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.xaxis.set_major_formatter(FuncFormatter(mb_formatter))
    ax.tick_params(axis="x", labelbottom=show_x)


def add_stats_box(ax, text: str) -> None:
    ax.text(
        1.015,
        0.98,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        color="#222222",
        clip_on=False,
        bbox=dict(fc="white", ec="#d0d0d0", lw=0.7, alpha=0.92, boxstyle="round,pad=0.35"),
    )


def add_badge(ax, text: str, edgecolor: str) -> None:
    ax.text(
        1.015,
        0.64,
        textwrap.fill(text, width=34),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        color="#1f1f1f",
        clip_on=False,
        bbox=dict(fc="white", ec=edgecolor, lw=0.9, alpha=0.92, boxstyle="round,pad=0.35"),
    )


def add_major_landmarks(ax, spec: WindowSpec, snrpn_pos: int) -> None:
    for cluster in BP_CLUSTER_INTERVALS_T2T:
        start, end = clipped_interval(cluster["start"], cluster["end"], spec)
        if start >= end:
            continue
        ax.axvspan(start, end, color=structural_landmark_color(cluster["name"]), alpha=0.10, lw=0, zorder=0)
    if spec.start <= snrpn_pos <= spec.end:
        ax.axvline(snrpn_pos, color=DARK_GREY, lw=1.0, alpha=0.65, zorder=2)


def add_logic_block(ax, x: float, title: str, top_label: str, bottom_label: str, top_style: str, bottom_style: str) -> None:
    width = 0.18
    box_h = 0.18
    ax.text(x + width / 2, 0.93, title, ha="center", va="top", fontsize=13, fontweight="bold")
    for y0, label, style in [
        (0.60, top_label, top_style),
        (0.35, bottom_label, bottom_style),
    ]:
        fc = "white"
        ec = STRUCT_GREY
        hatch = None
        text_color = "#222222"
        if style in {"maternal", "retained_maternal", "maternal_haplotype_1", "maternal_haplotype_2"}:
            fc = MATERNAL if style in {"maternal", "retained_maternal"} else MUPD_1 if style == "maternal_haplotype_1" else MUPD_2
            ec = fc
            text_color = "white"
        elif style in {"paternal", "retained_paternal"}:
            fc = PATERNAL
            ec = fc
            text_color = "white"
        elif "deleted" in style:
            fc = "#f4f4f4"
            ec = STRUCT_GREY
            hatch = "////"
        rect = FancyBboxPatch(
            (x, y0),
            width,
            box_h,
            boxstyle="round,pad=0.02,rounding_size=0.025",
            fc=fc,
            ec=ec,
            lw=1.0,
            hatch=hatch,
        )
        ax.add_patch(rect)
        ax.text(x + width / 2, y0 + box_h / 2, label, ha="center", va="center", fontsize=11, color=text_color)


def draw_logic_panel(ax) -> None:
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(
        0.0,
        1.02,
        "A. Reciprocal study logic",
        ha="left",
        va="bottom",
        fontsize=16,
        fontweight="bold",
    )
    blocks = [
        ("Control", "maternal", "paternal", "maternal", "paternal"),
        ("PWS-DEL", "retained maternal", "paternal deleted", "retained_maternal", "paternal_deleted"),
        ("AS-DEL", "maternal deleted", "retained paternal", "maternal_deleted", "retained_paternal"),
        ("PWS-mUPD", "maternal haplotype 1", "maternal haplotype 2", "maternal_haplotype_1", "maternal_haplotype_2"),
    ]
    xs = [0.03, 0.285, 0.54, 0.795]
    for x, (title, top_label, bottom_label, top_style, bottom_style) in zip(xs, blocks):
        add_logic_block(ax, x, title, top_label, bottom_label, top_style, bottom_style)
    ax.text(
        0.5,
        0.07,
        "Reciprocal deletions expose parental cis-methylation architecture, while PWS-mUPD confirms duplicated maternal identity.",
        ha="center",
        va="center",
        fontsize=12,
        color="#333333",
    )


def draw_annotation_track(ax, annotation_track: pd.DataFrame, spec: WindowSpec) -> None:
    ax.set_xlim(spec.start, spec.end)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.tick_params(axis="x", labelbottom=False)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.xaxis.set_major_formatter(FuncFormatter(mb_formatter))
    ax.hlines(0.73, spec.start, spec.end, color="#d8d8d8", lw=0.8)
    ax.hlines(0.58, spec.start, spec.end, color="#d8d8d8", lw=0.8)
    ax.hlines(0.43, spec.start, spec.end, color="#d8d8d8", lw=0.8)
    ax.text(
        0.0,
        1.06,
        "Shared T2T annotation track",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11.5,
        fontweight="bold",
    )

    lane_y = {
        "structural": 0.73,
        "anchor": 0.58,
        "gene_track": 0.43,
    }
    top_label_y = {0: 0.92, 1: 1.02, 2: 1.12, 3: 1.22}
    bottom_label_y = {0: 0.28, 1: 0.18, 2: 0.08, 3: -0.02}
    xtrans = ax.get_xaxis_transform()
    for _, row in annotation_track.iterrows():
        y = lane_y[row["lane"]]
        start, end = int(row["start"]), int(row["end"])
        if row["kind"] == "structural_landmark":
            rect_fc = row.get("display_color", STRUCT_GREY)
            rect_ec = row.get("display_color", STRUCT_GREY)
            alpha = 0.30
            height = 0.10
        elif row["kind"] == "regulatory_anchor":
            rect_fc = DARK_GREY
            rect_ec = DARK_GREY
            alpha = 0.25
            height = 0.10
        else:
            rect_fc = "#bdbdbd"
            rect_ec = "#7f7f7f"
            alpha = 0.45
            height = 0.09
        ax.add_patch(
            Rectangle(
                (start, y - height / 2),
                max(1, end - start),
                height,
                fc=rect_fc,
                ec=rect_ec,
                alpha=alpha,
                lw=0.8,
            )
        )
        center = float(row.get("center", (start + end) / 2))
        if row["kind"] in {"structural_landmark", "regulatory_anchor"}:
            label_y = 0.84 if row["kind"] == "structural_landmark" else 0.69
            ax.plot([center, center], [y + height / 2, label_y - 0.03], color="#757575", lw=0.7, alpha=0.8, transform=xtrans)
            ax.text(
                center,
                label_y,
                row["label"],
                transform=xtrans,
                ha="center",
                va="bottom",
                fontsize=10,
                color="#222222",
                clip_on=False,
            )
        else:
            lane_index = int(row.get("lane_index", 0))
            if row.get("label_side") == "top":
                label_y = top_label_y.get(lane_index, 1.20)
                line_end = label_y - 0.035
                va = "bottom"
            else:
                label_y = bottom_label_y.get(lane_index, 0.04)
                line_end = label_y + 0.035
                va = "top"
            label_x = float(row.get("label_x", center))
            ax.plot([center, center], [y + height / 2, 0.50], color="#969696", lw=0.55, alpha=0.75, transform=xtrans, clip_on=False)
            ax.plot([center, label_x], [0.50, line_end], color="#969696", lw=0.55, alpha=0.75, transform=xtrans, clip_on=False)
            ax.text(
                label_x,
                label_y,
                row["label"],
                transform=xtrans,
                ha="center",
                va=va,
                fontsize=8.6,
                color="#222222",
                clip_on=False,
            )
        if pd.notna(row["marker_position"]):
            marker = int(row["marker_position"])
            if spec.start <= marker <= spec.end:
                ax.axvline(marker, color=DARK_GREY, lw=1.1, alpha=0.75)

    ax.text(spec.start, 0.765, "BP landmarks", ha="left", va="bottom", fontsize=9.5, color="#5a5a5a")
    ax.text(spec.start, 0.615, "PWS imprinting center", ha="left", va="bottom", fontsize=9.5, color="#5a5a5a")
    ax.text(spec.start, 0.465, "Broader gene context", ha="left", va="bottom", fontsize=9.5, color="#5a5a5a")


def plot_series(ax, x: pd.Series, y: pd.Series, color: str, label: str, lw: float = 2.0, ls: str = "-") -> None:
    ax.plot(x, y, color=color, lw=lw, ls=ls, label=label)


def boundary_highlights(boundaries: pd.DataFrame, gene_features: pd.DataFrame) -> List[Tuple[int, int, str]]:
    if boundaries.empty or gene_features.empty:
        return []
    features = gene_features.set_index("feature")
    regions: List[Tuple[int, int, str]] = []
    for label, feature_names in [
        ("localized parental divergence", ["SNURF-SNRPN", "UBE3A-ATS"]),
        ("localized parental divergence", ["SNORD116", "UBE3A-ATS"]),
    ]:
        feature_intervals = []
        for name in feature_names:
            if name in features.index:
                row = features.loc[name]
                feature_intervals.append((int(row["start"]), int(row["end"])))
        if not feature_intervals:
            continue
        fstart = min(x[0] for x in feature_intervals)
        fend = max(x[1] for x in feature_intervals)
        overlap = boundaries[(boundaries["start"] < fend) & (boundaries["end"] > fstart)].copy()
        if overlap.empty:
            continue
        best = overlap.sort_values("max_abs_delta_smoothed", ascending=False).iloc[0]
        regions.append((int(best["start"]), int(best["end"]), label))
    unique = []
    seen = set()
    for start, end, label in regions:
        key = (start, end, label)
        if key not in seen:
            unique.append((start, end, label))
            seen.add(key)
    return unique


def render_figure(
    outdir: Path,
    display_spec: WindowSpec,
    annotation_track: pd.DataFrame,
    panel_tables: Dict[str, pd.DataFrame],
    boundaries: pd.DataFrame,
    gene_features: pd.DataFrame,
    stats: Dict[str, PanelStats],
) -> None:
    fig = plt.figure(figsize=(19.5, 15.9))
    gs = GridSpec(
        6,
        1,
        height_ratios=[1.10, 1.14, 1.14, 1.14, 1.14, 1.22],
        hspace=0.34,
        figure=fig,
    )

    anno_ax = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0], sharex=anno_ax)
    ax_c = fig.add_subplot(gs[2, 0], sharex=anno_ax)
    ax_d = fig.add_subplot(gs[3, 0], sharex=anno_ax)
    ax_e = fig.add_subplot(gs[4, 0], sharex=anno_ax)
    ax_f = fig.add_subplot(gs[5, 0], sharex=anno_ax)

    draw_annotation_track(anno_ax, annotation_track, display_spec)
    snrpn_markers = annotation_track.loc[annotation_track["kind"] == "regulatory_anchor", "marker_position"].dropna()
    if snrpn_markers.empty:
        raise ValueError("Annotation track is missing the PWS-IC/SNRPN anchor row.")
    snrpn_pos = int(snrpn_markers.iloc[0])

    # Panel A
    table_b = panel_tables["Figure2B_control_reference_plot.tsv"]
    add_shared_axis_style(ax_b, display_spec, "Mean methylation", show_x=False)
    ax_b.set_ylim(-0.02, 1.04)
    add_major_landmarks(ax_b, display_spec, snrpn_pos)
    plot_series(ax_b, table_b["window_mid"], table_b["control_maternal_mean_smoothed"], MATERNAL, "control maternal")
    plot_series(ax_b, table_b["window_mid"], table_b["control_paternal_mean_smoothed"], PATERNAL, "control paternal")
    ax_b.set_title(stats["A"].title, loc="left", fontsize=15, fontweight="bold", pad=8)
    add_badge(ax_b, stats["A"].badge, MATERNAL)
    add_stats_box(ax_b, stats["A"].stats_box)
    ax_b.legend(loc="lower left", bbox_to_anchor=(1.015, 0.03), fontsize=10, frameon=False, borderaxespad=0.0)

    # Panel B
    table_c = panel_tables["Figure2C_pwsdel_retained_maternal_plot.tsv"]
    add_shared_axis_style(ax_c, display_spec, "Mean methylation", show_x=False)
    ax_c.set_ylim(-0.02, 1.04)
    add_major_landmarks(ax_c, display_spec, snrpn_pos)
    plot_series(ax_c, table_c["window_mid"], table_c["control_maternal_mean_smoothed"], MATERNAL, "control maternal reference", lw=1.6, ls="--")
    plot_series(ax_c, table_c["window_mid"], table_c["pwsdel_retained_maternal_mean_smoothed"], MATERNAL, "PWS-DEL retained maternal", lw=2.2)
    ax_c.set_title(stats["B"].title, loc="left", fontsize=15, fontweight="bold", pad=8)
    add_badge(ax_c, stats["B"].badge, MATERNAL)
    add_stats_box(ax_c, stats["B"].stats_box)
    ax_c.legend(loc="lower left", bbox_to_anchor=(1.015, 0.03), fontsize=10, frameon=False, borderaxespad=0.0)

    # Panel C
    table_d = panel_tables["Figure2D_asdel_retained_paternal_plot.tsv"]
    add_shared_axis_style(ax_d, display_spec, "Mean methylation", show_x=False)
    ax_d.set_ylim(-0.02, 1.04)
    add_major_landmarks(ax_d, display_spec, snrpn_pos)
    plot_series(ax_d, table_d["window_mid"], table_d["control_paternal_mean_smoothed"], PATERNAL, "control paternal reference", lw=1.6, ls="--")
    plot_series(ax_d, table_d["window_mid"], table_d["asdel_retained_paternal_mean_smoothed"], PATERNAL, "AS-DEL retained paternal", lw=2.2)
    ax_d.set_title(stats["C"].title, loc="left", fontsize=15, fontweight="bold", pad=8)
    add_badge(ax_d, stats["C"].badge, PATERNAL)
    add_stats_box(ax_d, stats["C"].stats_box)
    ax_d.legend(loc="lower left", bbox_to_anchor=(1.015, 0.03), fontsize=10, frameon=False, borderaxespad=0.0)

    # Panel D
    table_e = panel_tables["Figure2E_pwsmupd_maternal_plot.tsv"]
    add_shared_axis_style(ax_e, display_spec, "Mean methylation", show_x=False)
    ax_e.set_ylim(-0.02, 1.04)
    add_major_landmarks(ax_e, display_spec, snrpn_pos)
    plot_series(ax_e, table_e["window_mid"], table_e["control_maternal_mean_smoothed"], MATERNAL, "control maternal reference", lw=1.5, ls="--")
    plot_series(ax_e, table_e["window_mid"], table_e["mean_meth_haplotype1_smoothed"], MUPD_1, "004P haplotype 1", lw=2.0)
    plot_series(ax_e, table_e["window_mid"], table_e["mean_meth_haplotype2_smoothed"], MUPD_2, "004P haplotype 2", lw=2.0)
    ax_e.set_title(stats["D"].title, loc="left", fontsize=15, fontweight="bold", pad=8)
    add_badge(ax_e, stats["D"].badge, MUPD_2)
    add_stats_box(ax_e, stats["D"].stats_box)
    ax_e.legend(loc="lower left", bbox_to_anchor=(1.015, 0.03), fontsize=10, frameon=False, borderaxespad=0.0)

    # Panel E
    table_f = panel_tables["Figure2F_reciprocal_delta_architecture_plot.tsv"]
    add_shared_axis_style(ax_f, display_spec, "Delta methylation\n(PWS-DEL maternal - AS-DEL paternal)", show_x=True)
    ax_f.set_ylim(-0.33, 0.18)
    add_major_landmarks(ax_f, display_spec, snrpn_pos)
    ax_f.axhline(0.0, color=DARK_GREY, lw=1.0, alpha=0.8)
    delta = table_f["reciprocal_delta_smoothed"]
    ax_f.fill_between(table_f["window_mid"], 0, delta, where=delta >= 0, color=MATERNAL, alpha=0.16, lw=0)
    ax_f.fill_between(table_f["window_mid"], 0, delta, where=delta < 0, color=PATERNAL, alpha=0.18, lw=0)
    plot_series(ax_f, table_f["window_mid"], delta, DARK_GREY, "PWS-DEL maternal minus AS-DEL paternal", lw=2.1)

    for start, end, _label in boundary_highlights(boundaries, gene_features):
        ax_f.axvspan(start, end, color="#ffd9a8", alpha=0.22, lw=0)

    ax_f.annotate(
        "shared methylation scaffold",
        xy=(24_250_000, 0.01),
        xytext=(24_250_000, 0.12),
        ha="center",
        va="bottom",
        fontsize=10,
        color="#333333",
        arrowprops=dict(arrowstyle="-", color="#666666", lw=0.8),
    )
    ax_f.set_title(stats["E"].title, loc="left", fontsize=15, fontweight="bold", pad=8)
    add_badge(ax_f, stats["E"].badge, PATERNAL)
    add_stats_box(ax_f, stats["E"].stats_box)
    ax_f.legend(loc="lower left", bbox_to_anchor=(1.015, 0.03), fontsize=10, frameon=False, borderaxespad=0.0)
    ax_f.set_xlabel("chr15 coordinate, T2T-CHM13v2.0 (Mb)", fontsize=14)

    fig.subplots_adjust(top=0.985, bottom=0.06, left=0.08, right=0.80)

    figdir = mkdir(outdir / "figures")
    fig.savefig(figdir / "Figure2_reciprocal_cis_architecture_improved.png", dpi=400, bbox_inches="tight")
    fig.savefig(figdir / "Figure2_reciprocal_cis_architecture_improved.pdf", bbox_inches="tight")
    fig.savefig(figdir / "Figure2_reciprocal_cis_architecture_improved.svg", bbox_inches="tight")
    plt.close(fig)


def report_text(
    outdir: Path,
    report_path: Path,
    input_files: List[str],
    output_files: List[str],
    control_ref: pd.DataFrame,
    corr: pd.DataFrame,
    reciprocal_delta: pd.DataFrame,
    boundaries: pd.DataFrame,
    smooth_window: int,
) -> str:
    timestamp = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %z")
    control_stats = pearson_summary(
        control_ref["control_maternal_mean"],
        control_ref["control_paternal_mean"],
    )
    pws = corr[corr["group"] == "PWS_DEL"].copy()
    asdel = corr[corr["group"] == "AS_DEL"].copy()
    upd = corr[corr["group"] == "PWS_mUPD"].copy().sort_values("comparison")
    reciprocal_nonnull = reciprocal_delta.dropna(
        subset=["pwsdel_retained_maternal_mean", "asdel_retained_paternal_mean"]
    )
    reciprocal_stats = pearson_summary(
        reciprocal_nonnull["pwsdel_retained_maternal_mean"],
        reciprocal_nonnull["asdel_retained_paternal_mean"],
    )

    rel = lambda path: path.relative_to(report_path.parent.parent).as_posix()
    lines = [
        "# Figure 2 Report: Reciprocal cis-methylation architecture across the 15q11-q13 imprinted domain",
        "",
        f"Generated: {timestamp}",
        "",
        "## 1. Purpose of Figure 2",
        "Figure 2 tests whether reciprocal deletion genomes reconstruct the parental cis-methylation architecture across the T2T chr15 imprinted domain: controls provide the biparental reference, PWS-DEL exposes the retained maternal-only profile, AS-DEL exposes the retained paternal-only profile, and PWS-mUPD confirms duplicated maternal identity.",
        "",
        "## 2. Input files",
    ]
    lines.extend([f"- `{path}`" for path in input_files])
    lines.extend(
        [
            "",
            "## 3. Coordinate system",
            "- Reference: T2T-CHM13v2.0",
            "- Analysis domain: `chr15:17,600,000-28,000,000`",
            "- Display domain in the improved figure: `chr15:18,000,000-28,000,000`",
            "- Major landmarks retained in the main figure: BP1, BP2, BP3, and PWS-IC/SNRPN",
            "",
            "## 4. Window size and smoothing parameters",
            f"- Base window size: 1 kb",
            f"- Smoothed display tracks: centered rolling median, {smooth_window} windows",
            f"- Approximate smoothing span: {smooth_window} kb",
            "",
            "## 5. Number of informative windows per panel",
            f"- Shared annotation track: structural and gene context only, no methylation windows",
            f"- Panel A (Controls): {control_stats['n_windows']:,} informative windows with both control maternal and control paternal means",
            f"- Panel B (PWS-DEL): {int(pws['n_windows'].sum()):,} within-footprint sample-windows across 5 deletion genomes",
            f"- Panel C (AS-DEL): {int(asdel['n_windows'].sum()):,} within-footprint sample-windows across 3 deletion genomes",
            f"- Panel D (PWS-mUPD): {int(upd.iloc[0]['n_windows']):,} windows for haplotype 1 and {int(upd.iloc[1]['n_windows']):,} windows for haplotype 2",
            f"- Panel E (Delta architecture): {reciprocal_stats['n_windows']:,} shared windows with both PWS-DEL retained maternal and AS-DEL retained paternal means",
            "",
            "## 6. Correlation statistics",
            f"- Controls: `r = {control_stats['pearson_r']:.3f}` across {control_stats['n_windows']:,} informative windows",
            f"- PWS-DEL retained maternal versus control maternal: mean `r = {pws['pearson_r'].mean():.3f}`, range `{pws['pearson_r'].min():.3f}-{pws['pearson_r'].max():.3f}`",
            f"- AS-DEL retained paternal versus control paternal: mean `r = {asdel['pearson_r'].mean():.3f}`, range `{asdel['pearson_r'].min():.3f}-{asdel['pearson_r'].max():.3f}`",
            f"- PWS-mUPD haplotype 1 versus control maternal: `r = {upd.iloc[0]['pearson_r']:.3f}`",
            f"- PWS-mUPD haplotype 2 versus control maternal: `r = {upd.iloc[1]['pearson_r']:.3f}`",
            f"- Reciprocal overlay (PWS-DEL retained maternal versus AS-DEL retained paternal): `r = {reciprocal_stats['pearson_r']:.3f}` across {reciprocal_stats['n_windows']:,} shared windows",
            "",
            "## 7. Interpretation of each panel",
            "- Shared annotation track: provides T2T structural and gene context for BP1, BP2, BP3, the PWS imprinting center, and broader 15q11-q13 genes.",
            "- Panel A: controls provide the biparental reference and show the reciprocal maternal and paternal methylation architecture used throughout the rest of the figure.",
            "- Panel B: PWS-DEL retained methylation closely follows the control maternal reference, supporting exposure of a maternal-only architecture by paternal deletion.",
            "- Panel C: AS-DEL retained methylation closely follows the control paternal reference, supporting exposure of a paternal-only architecture by maternal deletion.",
            "- Panel D: both 004P haplotypes track the control maternal reference, consistent with copy-neutral duplicated maternal architecture rather than biparental inheritance.",
            "- Panel E: the delta track is near zero across broad intervals, indicating a shared methylation scaffold, but shows localized parental divergence centered on the SNRPN/SNHG14 and SNORD116 interval.",
            "",
            "## 8. BP1/BP2/BP3 are structural landmarks, not methylation-boundary calls",
            "BP1, BP2, and BP3 are shown only as T2T structural breakpoint-cluster landmarks. They are not interpreted here as methylation boundaries. The localized methylation divergences highlighted in Panel E are data-driven parental differences within the imprinted gene domain, not calls anchored to BP1/BP2/BP3.",
            "",
            "## 9. Main conclusion",
            "Reciprocal PWS and AS deletion genomes reconstruct the parental cis-methylation architecture across the chr15 imprinted domain. PWS-DEL reveals the retained maternal scaffold, AS-DEL reveals the retained paternal scaffold, PWS-mUPD independently confirms duplicated maternal identity, and controls supply the biparental reference needed to interpret localized parental divergence near SNRPN/SNHG14 and SNORD116.",
            "",
            "## 10. Output file list",
        ]
    )
    lines.extend([f"- `{path}`" for path in output_files])
    return "\n".join(lines) + "\n"


def write_panel_tables(outdir: Path, annotation_track: pd.DataFrame, panel_tables: Dict[str, pd.DataFrame]) -> List[str]:
    tbldir = mkdir(outdir / "tables")
    outputs = []
    annotation_name = "Figure2_shared_annotation_track.tsv"
    annotation_track.to_csv(tbldir / annotation_name, sep="\t", index=False, na_rep="NaN")
    outputs.append(f"tables/{annotation_name}")
    for name, table in panel_tables.items():
        table.to_csv(tbldir / name, sep="\t", index=False, na_rep="NaN")
        outputs.append(f"tables/{name}")
    return outputs


def _figure2_render_main() -> None:
    parser = argparse.ArgumentParser(description="Create the improved Figure 2 from existing Phase 2 outputs")
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR, type=Path)
    parser.add_argument("--display-start", default=DISPLAY_START, type=int)
    parser.add_argument("--display-end", default=DISPLAY_END, type=int)
    parser.add_argument("--smooth-window", default=31, type=int)
    parser.add_argument("--gtf", default=DEFAULT_GTF, type=Path)
    parser.add_argument("--icr-bed", default=DEFAULT_ICR_BED, type=Path)
    args = parser.parse_args()

    outdir = args.outdir
    tbldir = outdir / "tables"
    reports_dir = mkdir(outdir / "reports")

    log("Loading existing Phase 2 tables")
    run_summary = read_tsv(outdir / "phase2_run_summary.tsv")
    window_size = int(run_summary.iloc[0]["window_size"])
    analysis_spec = WindowSpec(
        chrom=str(run_summary.iloc[0]["domain"]).split(":")[0],
        start=int(str(run_summary.iloc[0]["domain"]).split(":")[1].split("-")[0]),
        end=int(str(run_summary.iloc[0]["domain"]).split(":")[1].split("-")[1]),
        size=window_size,
    )
    display_spec = WindowSpec(analysis_spec.chrom, args.display_start, args.display_end, DISPLAY_WINDOW)

    control_ref = read_tsv(tbldir / "Phase2_control_reference_architecture.tsv")
    retained = read_tsv(tbldir / "Phase2_retained_haplotype_profiles.tsv")
    all_windows = read_tsv(tbldir / "Phase2_all_samples_haplotype_1kb_methylation.tsv.gz")
    corr = read_tsv(tbldir / "Phase2_retained_and_mUPD_correlations.tsv")
    reciprocal_delta = read_tsv(tbldir / "Phase2_reciprocal_delta_profile.tsv")
    boundaries = read_tsv(tbldir / "Phase2_reciprocal_delta_boundary_candidates.tsv")
    gene_features = read_tsv(tbldir / "Phase2_gene_track_features.tsv")

    annotation_track = build_annotation_track_table(gene_features, args.gtf, args.icr_bed, display_spec)
    panel_tables = build_panel_tables(
        display_spec,
        control_ref,
        retained,
        all_windows,
        reciprocal_delta,
        args.smooth_window,
    )
    stats = panel_stats(control_ref, corr, reciprocal_delta)

    log("Writing panel plotting tables")
    output_table_paths = write_panel_tables(outdir, annotation_track, panel_tables)

    log("Rendering improved Figure 2")
    render_figure(
        outdir,
        display_spec,
        annotation_track,
        panel_tables,
        boundaries,
        gene_features,
        stats,
    )

    report_path = reports_dir / "Figure2_reciprocal_cis_architecture_report.md"
    input_files = [
        "phase2_run_summary.tsv",
        "tables/Phase2_control_reference_architecture.tsv",
        "tables/Phase2_retained_haplotype_profiles.tsv",
        "tables/Phase2_all_samples_haplotype_1kb_methylation.tsv.gz",
        "tables/Phase2_retained_and_mUPD_correlations.tsv",
        "tables/Phase2_reciprocal_delta_profile.tsv",
        "tables/Phase2_reciprocal_delta_boundary_candidates.tsv",
        "tables/Phase2_gene_track_features.tsv",
        str(args.gtf),
        str(args.icr_bed),
        "scripts/paper_vf/create_figure2_reciprocal_cis_architecture_improved.py",
    ]
    output_files = [
        "figures/Figure2_reciprocal_cis_architecture_improved.png",
        "figures/Figure2_reciprocal_cis_architecture_improved.pdf",
        "figures/Figure2_reciprocal_cis_architecture_improved.svg",
        "reports/Figure2_reciprocal_cis_architecture_report.md",
    ] + output_table_paths
    report = report_text(
        outdir,
        report_path,
        input_files,
        output_files,
        control_ref,
        corr,
        reciprocal_delta,
        boundaries,
        args.smooth_window,
    )
    report_path.write_text(report)
    log("Done")


def _forward_args(argv: list[str], allowed_with_values: set[str], allowed_flags: set[str]) -> list[str]:
    forwarded: list[str] = []
    idx = 0
    while idx < len(argv):
        token = argv[idx]
        matched_key = next((key for key in allowed_with_values if token.startswith(key + '=')), None)
        if matched_key is not None:
            forwarded.append(token)
            idx += 1
            continue
        if token in allowed_with_values:
            forwarded.append(token)
            if idx + 1 < len(argv):
                forwarded.append(argv[idx + 1])
                idx += 2
            else:
                idx += 1
            continue
        if token in allowed_flags:
            forwarded.append(token)
        idx += 1
    return forwarded


def main() -> None:
    original_argv = sys.argv[:]
    try:
        sys.argv = original_argv[:]
        _phase2_main()
        render_args = _forward_args(
            original_argv[1:],
            allowed_with_values={
                '--outdir',
                '--display-start',
                '--display-end',
                '--smooth-window',
                '--gtf',
                '--icr-bed',
            },
            allowed_flags=set(),
        )
        sys.argv = [original_argv[0], *render_args]
        _figure2_render_main()
    finally:
        sys.argv = original_argv


if __name__ == "__main__":
    main()
