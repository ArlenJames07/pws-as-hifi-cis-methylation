"""Readers for pipeline outputs: sample file lookup, VCF (pbsv, HiFiCNV), HiFiCNV
copy-number bedGraph, FASTA index."""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd


def sample_files(root: Path, sample: str, patterns: tuple[str, ...]) -> list[Path]:
    """Files under root (first root/<sample>/, then root, then recursively) matching one of
    the glob patterns, whose name contains the sample ID as a whole token (so SMRT Link
    names like 08_1_C01_bc2049_007P.bam are found and 007P never matches 1007P)."""
    token = re.compile(rf"(?:^|[._-]){re.escape(sample)}(?:[._-]|$)")
    for base, recursive in ((root / sample, False), (root, False), (root, True)):
        if not base.is_dir():
            continue
        hits = []
        for pat in patterns:
            it = base.rglob(pat) if recursive else base.glob(pat)
            hits += [p for p in it if p.is_file() and token.search(p.name)]
        if hits:
            return sorted(set(hits))
    return []


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith((".gz", ".bgz")) else open(path)


def _info(text: str) -> dict[str, str]:
    out = {}
    for item in text.split(";"):
        k, _, v = item.partition("=")
        out[k] = v if _ else "1"
    return out


def read_vcf(path: Path, pass_only: bool = True) -> pd.DataFrame:
    """chrom, pos (1-based), end, svtype, svlen, alt, filter, info (dict), format fields of
    the first sample (dict)."""
    rows = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            if pass_only and f[6] not in ("PASS", "."):
                continue
            info = _info(f[7])
            fmt = dict(zip(f[8].split(":"), f[9].split(":"))) if len(f) > 9 else {}
            pos = int(f[1])
            svtype = info.get("SVTYPE") or (f[4].strip("<>").split(":")[0] if f[4].startswith("<") else "")
            try:
                end = int(info.get("END", pos))
            except ValueError:
                end = pos
            try:
                svlen = abs(int(str(info.get("SVLEN", "0")).split(",")[0]))
            except ValueError:
                svlen = 0
            if not svlen and svtype in ("DEL", "DUP", "INV", "CNV") and end > pos:
                svlen = end - pos
            if not svlen and svtype == "" and len(f[3]) != len(f[4]):
                svlen = abs(len(f[3]) - len(f[4]))
            rows.append((f[0], pos, end, svtype, svlen, f[4], f[6], info, fmt))
    return pd.DataFrame(rows, columns=["chrom", "pos", "end", "svtype", "svlen", "alt", "filter", "info", "format"])


def hificnv_calls(path: Path) -> pd.DataFrame:
    """Copy-number calls of a HiFiCNV VCF: chrom, start, end, size_bp, svtype (DEL/DUP),
    copy_number (FORMAT CN when present)."""
    v = read_vcf(path, pass_only=True)
    v = v[v["svtype"].isin(["DEL", "DUP"])].copy()
    if v.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "size_bp", "svtype", "copy_number"])
    v["start"] = v["pos"] - 1
    v["size_bp"] = (v["end"] - v["start"]).clip(lower=0)
    v["copy_number"] = [pd.to_numeric(f.get("CN"), errors="coerce") for f in v["format"]]
    return v[["chrom", "start", "end", "size_bp", "svtype", "copy_number"]].reset_index(drop=True)


def hificnv_copynum(path: Path, chrom: str | None = None) -> pd.DataFrame:
    """HiFiCNV copy-number bedGraph (chrom, start, end, copy_number)."""
    rows = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(("#", "track")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or (chrom and f[0] != chrom):
                continue
            try:
                rows.append((f[0], int(f[1]), int(f[2]), float(f[3])))
            except ValueError:
                continue
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "copy_number"])


def sv_counts(path: Path, types=("DEL", "INS", "DUP", "INV", "BND")) -> dict:
    """pbsv PASS calls: count per type and total span (|SVLEN|, BND excluded) in Mb."""
    v = read_vcf(path, pass_only=True)
    out = {t: int((v["svtype"] == t).sum()) for t in types}
    out["TOTAL_COUNT"] = int(sum(out[t] for t in types))
    out["TOTAL_SPAN_MB"] = float(v.loc[v["svtype"] != "BND", "svlen"].sum()) / 1e6
    return out


def fai_lengths(path: Path | None, keep=None) -> dict[str, int]:
    if path is None or not Path(path).is_file():
        return {}
    out = {}
    with open(path) as handle:
        for line in handle:
            f = line.split("\t")
            if len(f) >= 2 and (keep is None or f[0] in keep):
                out[f[0]] = int(f[1])
    return out


CHM13_LENGTHS = {   # fallback when the reference index is not available
    "chr1": 248387328, "chr2": 242696752, "chr3": 201105948, "chr4": 193574945, "chr5": 182045439,
    "chr6": 172126628, "chr7": 160567428, "chr8": 146259331, "chr9": 150617247, "chr10": 134758134,
    "chr11": 135127769, "chr12": 133324548, "chr13": 113566686, "chr14": 101161492, "chr15": 99753195,
    "chr16": 96330374, "chr17": 84276897, "chr18": 80542538, "chr19": 61707364, "chr20": 66210255,
    "chr21": 45090682, "chr22": 51324926, "chrX": 154259566, "chrY": 62460029}


def interval_overlap(a0, a1, b0, b1) -> float:
    return max(0.0, min(a1, b1) - max(a0, b0))


def weighted_mean(values, weights) -> float:
    v = np.asarray(values, float)
    w = np.asarray(weights, float)
    ok = np.isfinite(v) & np.isfinite(w) & (w > 0)
    return float(np.sum(v[ok] * w[ok]) / np.sum(w[ok])) if ok.any() else np.nan
