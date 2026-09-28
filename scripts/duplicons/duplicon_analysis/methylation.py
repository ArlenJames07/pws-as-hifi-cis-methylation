"""
CpG methylation readers for the duplicon analyses.

Same interface as scripts/analysis/cis_analysis.methylation (find_track,
read_track -> MethylationTrack), with one difference that matters: the 0-100 vs
0-1 scale is decided once per FILE, never per value. pb-CpG-tools writes
mod_score on a 0-100 scale, so a site scored 0.5 is 0.5 %, not beta 0.5.

Formats (auto-detected from the header and the column count):
  pb-CpG-tools   chrom begin end mod_score type cov est_mod est_unmod [...]   -> percent
  modkit         bedMethyl, percent modified in column 11                   -> percent
  plain 5-col    chrom start end beta cov (tracks written by these scripts) -> fraction
"""
from __future__ import annotations

import gzip
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np


@dataclass(frozen=True)
class MethylationTrack:
    position: np.ndarray
    beta: np.ndarray
    coverage: np.ndarray
    source: Path


def find_track(root: Path, sample_id: str, track_kind: str) -> Path | None:
    """<root>/<sample>/<sample>.cpg.<kind>.bed[.gz], as written by the workflow; otherwise
    one pb-CpG-tools file under <root> named *.<kind>.bed[.gz] that contains the sample ID
    as a token (e.g. 08_1_A01_bc2043_001P.combined.bed.gz). Several matches: ValueError."""
    for name in (f"{sample_id}.cpg.{track_kind}.bed.gz", f"{sample_id}.cpg.{track_kind}.bed",
                 f"{sample_id}.{track_kind}.bed.gz", f"{sample_id}.{track_kind}.bed"):
        for candidate in (root / sample_id / name, root / name):
            if candidate.is_file():
                return candidate
    if not root.is_dir():
        return None
    token = re.compile(rf"(?:^|[._-]){re.escape(sample_id)}(?:[._-]|$)")
    hits = sorted({p.resolve() for suffix in (f".{track_kind}.bed", f".{track_kind}.bed.gz")
                   for p in root.rglob(f"*{suffix}") if p.is_file()
                   and token.search(p.name[: -len(suffix)])})
    if len(hits) > 1:
        raise ValueError(f"ambiguous {track_kind} methylation tracks for {sample_id}: "
                         + ", ".join(map(str, hits)))
    return hits[0] if hits else None


def _lines(path: Path) -> Iterator[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as handle:
        yield from handle


def track_format(path: Path) -> str:
    """'pbcpg', 'bedmethyl' or 'plain', from the header and the first data line."""
    for line in _lines(path):
        if line.startswith("##") and "pb-cpg-tools" in line.lower():
            return "pbcpg"
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 11 and fields[3] in ("m", "h", "a", "C", "5mC"):
            return "bedmethyl"
        if len(fields) >= 8:
            return "pbcpg"
        return "plain"
    return "plain"


def _region_lines(path: Path, chrom: str, start: int, end: int) -> Iterable[str]:
    if str(path).endswith(".gz") and Path(f"{path}.tbi").is_file():
        try:
            import pysam
            with pysam.TabixFile(str(path)) as tabix:
                yield from tabix.fetch(chrom, max(0, start), end)
            return
        except (ImportError, OSError, ValueError):
            pass
    if not str(path).endswith(".gz"):
        program = "$1==chrom && $2>=start && $2<end {print}"
        with subprocess.Popen(["awk", "-v", f"chrom={chrom}", "-v", f"start={start}", "-v", f"end={end}",
                               program, str(path)], stdout=subprocess.PIPE, text=True) as process:
            yield from process.stdout
        return
    yield from _lines(path)


def read_track(path: Path, chrom: str, start: int, end: int) -> MethylationTrack:
    fmt = track_format(path)
    pos, beta, cov = [], [], []
    for line in _region_lines(path, chrom, start, end):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[0] != chrom:
            continue
        try:
            p = int(f[1])
            if fmt == "bedmethyl":
                b, c = float(f[10]) / 100.0, float(f[9])
            elif fmt == "pbcpg":
                b, c = float(f[3]) / 100.0, float(f[5])
            else:
                b, c = float(f[3]), float(f[4]) if len(f) > 4 else 1.0
        except (ValueError, IndexError):
            continue
        if p < start or p >= end or not (0.0 <= b <= 1.0) or c <= 0:
            continue
        pos.append(p)
        beta.append(b)
        cov.append(c)
    if not pos:
        empty = np.array([], dtype=float)
        return MethylationTrack(np.array([], dtype=np.int64), empty, empty, path)
    order = np.argsort(pos, kind="stable")
    p = np.asarray(pos, dtype=np.int64)[order]
    b = np.asarray(beta, dtype=float)[order]
    c = np.asarray(cov, dtype=float)[order]
    if np.any(np.diff(p) == 0):          # same CpG listed twice: coverage-weighted mean
        u, inv = np.unique(p, return_inverse=True)
        num = np.bincount(inv, weights=b * c)
        den = np.bincount(inv, weights=c)
        p, b, c = u, num / den, den
    return MethylationTrack(p, b, c, path)


def element_means(position, beta, starts, ends, lo, hi, min_cpgs):
    """Mean beta of the CpGs inside each [start, end); NaN unless the element lies
    fully inside [lo, hi) and has at least min_cpgs CpGs. position must be sorted."""
    position = np.asarray(position)
    cs = np.concatenate(([0.0], np.cumsum(beta)))
    i0 = np.searchsorted(position, starts, side="left")
    i1 = np.searchsorted(position, ends, side="left")
    n = i1 - i0
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = (cs[i1] - cs[i0]) / np.maximum(n, 1)
    ok = (np.asarray(starts) >= lo) & (np.asarray(ends) <= hi) & (n >= min_cpgs)
    return np.where(ok, mean, np.nan), np.where(ok, n, 0)
