from __future__ import annotations

import gzip
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

    def __post_init__(self) -> None:
        if not (len(self.position) == len(self.beta) == len(self.coverage)):
            raise ValueError("Track arrays have unequal lengths")
        if len(self.position) and np.any(np.diff(self.position) <= 0):
            raise ValueError("CpG coordinates must be unique and sorted")


def find_track(root: Path, sample_id: str, track_kind: str) -> Path | None:
    suffixes = (
        f"{sample_id}.cpg.{track_kind}.bed.gz",
        f"{sample_id}.cpg.{track_kind}.bed",
        f"{sample_id}.{track_kind}.bed.gz",
        f"{sample_id}.{track_kind}.bed",
    )
    candidates: list[Path] = []
    sample_dir = root / sample_id
    for suffix in suffixes:
        direct = sample_dir / suffix
        if direct.is_file():
            candidates.append(direct)
        flat = root / suffix
        if flat.is_file():
            candidates.append(flat)
    if not candidates:
        for suffix in suffixes:
            candidates.extend(path for path in root.rglob(suffix) if path.is_file())
    unique = sorted(set(path.resolve() for path in candidates))
    if not unique:
        return None
    if len(unique) > 1:
        raise ValueError(f"Ambiguous {track_kind} tracks for {sample_id}: {unique}")
    return unique[0]


def _plain_lines(path: Path) -> Iterator[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        yield from handle


def read_track_metadata(path: Path) -> dict[str, str]:
    metadata: dict[str, str] = {}
    for line in _plain_lines(path):
        if line.startswith("##") and "=" in line:
            key, value = line[2:].rstrip("\n").split("=", 1)
            metadata[key.strip().replace("-", "_")] = value.strip()
            continue
        if line.startswith("#"):
            continue
        break
    return metadata


def _region_lines(path: Path, chrom: str, start: int, end: int) -> Iterable[str]:
    if path.suffix == ".gz" and Path(f"{path}.tbi").is_file():
        try:
            import pysam

            with pysam.TabixFile(str(path)) as tabix:
                yield from tabix.fetch(chrom, start, end)
            return
        except (ImportError, OSError, ValueError):
            pass
    if path.suffix != ".gz":
        program = "$1==chrom && $2>=start && $2<end {print}"
        with subprocess.Popen(
            [
                "awk",
                "-v",
                f"chrom={chrom}",
                "-v",
                f"start={start}",
                "-v",
                f"end={end}",
                program,
                str(path),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        ) as process:
            if process.stdout is None:
                raise RuntimeError("Failed to stream methylation BED")
            yield from process.stdout
            stderr = process.stderr.read() if process.stderr is not None else ""
            return_code = process.wait()
            if return_code != 0:
                raise RuntimeError(f"awk failed for {path}: {stderr.strip()}")
        return
    yield from _plain_lines(path)


def read_track(path: Path, chrom: str, start: int, end: int) -> MethylationTrack:
    positions: list[int] = []
    beta: list[float] = []
    coverage: list[float] = []
    for line in _region_lines(path, chrom, start, end):
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            fields = line.split()
        if len(fields) < 6 or fields[0] != chrom:
            continue
        try:
            position = int(fields[1])
            score = float(fields[3])
            depth = float(fields[5])
        except ValueError:
            continue
        if position < start or position >= end:
            continue
        value = score / 100.0 if score > 1.0 else score
        if not np.isfinite(value) or not np.isfinite(depth):
            continue
        if value < 0.0 or value > 1.0 or depth <= 0.0:
            continue
        positions.append(position)
        beta.append(value)
        coverage.append(depth)
    if not positions:
        return MethylationTrack(
            np.array([], dtype=np.int64),
            np.array([], dtype=float),
            np.array([], dtype=float),
            path,
        )
    order = np.argsort(np.asarray(positions), kind="stable")
    pos = np.asarray(positions, dtype=np.int64)[order]
    values = np.asarray(beta, dtype=float)[order]
    depths = np.asarray(coverage, dtype=float)[order]
    if np.any(np.diff(pos) == 0):
        unique, inverse = np.unique(pos, return_inverse=True)
        numerator = np.bincount(inverse, weights=values * depths)
        denominator = np.bincount(inverse, weights=depths)
        pos = unique
        values = numerator / denominator
        depths = denominator
    return MethylationTrack(pos, values, depths, path)
