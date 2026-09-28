"""
Find each participant's HiFi BAM(s) and check that they can be used.

Lookup order for one participant:
  1. bams.local.csv (columns sample,bam; several BAMs of one participant
     separated by ';') -- always wins, use it to resolve ambiguous names;
  2. <alignment_dir>/<sample>.aligned.bam        (the Nextflow workflow's name);
  3. any *.bam under <alignment_dir> whose name contains the sample ID as a
     separate token, e.g. 08_1_A01_bc2043_001P.bam (SMRT Link / demultiplexed
     names). Two or more matches (e.g. ..._bc2044_002P.bam and
     ..._bc2044v2_002P.bam) are refused: say which to use in bams.local.csv.
"""
from __future__ import annotations

import csv
import re
import subprocess
from pathlib import Path


class BamLookupError(Exception):
    """The BAM of a participant cannot be decided or does not exist."""


def read_bam_table(path: Path) -> dict[str, list[Path]]:
    """sample -> [bam, ...] from bams.local.csv; {} if the file does not exist."""
    if not path.is_file():
        return {}
    table: dict[str, list[Path]] = {}
    with open(path) as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not {"sample", "bam"} <= {f.strip() for f in reader.fieldnames}:
            raise SystemExit(f"{path}: needs the columns sample,bam")
        for row in reader:
            row = {k.strip(): (v or "").strip() for k, v in row.items() if k}
            if not row.get("sample") or not row.get("bam"):
                continue
            paths = []
            for value in row["bam"].split(";"):
                value = value.strip()
                if value:
                    p = Path(value)
                    paths.append(p if p.is_absolute() else path.parent / p)
            table[row["sample"]] = paths
    return table


def _token_pattern(sample: str) -> re.Pattern:
    return re.compile(rf"(?:^|[._-]){re.escape(sample)}(?:[._-]|$)")


def find_bams(sample: str, alignment_dir: Path, table: dict[str, list[Path]]) -> list[Path]:
    """The BAM(s) of one participant; [] if none exists. Raises BamLookupError
    when a listed BAM is missing, a link is broken, or the name is ambiguous."""
    if sample in table:
        missing = [p for p in table[sample] if not p.is_file()]
        if missing:
            raise BamLookupError(f"{sample}: listed in bams.local.csv but not found: "
                                 + ", ".join(map(str, missing)))
        return list(table[sample])
    exact = alignment_dir / f"{sample}.aligned.bam"
    if exact.is_file():
        return [exact]
    if not alignment_dir.is_dir():
        return []
    pattern = _token_pattern(sample)
    hits = sorted(p for p in alignment_dir.rglob("*.bam") if pattern.search(p.name[:-4]))
    broken = [p for p in hits if p.is_symlink() and not p.exists()]
    if broken:
        raise BamLookupError(f"{sample}: broken link(s): " + ", ".join(map(str, broken)))
    if len(hits) > 1:
        raise BamLookupError(
            f"{sample}: {len(hits)} BAMs match -- " + ", ".join(p.name for p in hits)
            + ". Put the one to use (or several, separated by ';', if they are different "
              "sequencing runs of this participant) in bams.local.csv."
        )
    return hits


def bam_index(bam: Path) -> Path | None:
    for candidate in (Path(f"{bam}.bai"), Path(f"{bam}.csi"), bam.with_suffix(".bai")):
        if candidate.is_file():
            return candidate
    return None


def header_contig_lengths(samtools: str, bam: Path) -> dict[str, int]:
    """@SQ name -> length from the BAM header."""
    text = subprocess.run([samtools, "view", "-H", str(bam)], check=True,
                          capture_output=True, text=True).stdout
    lengths = {}
    for line in text.splitlines():
        if line.startswith("@SQ"):
            fields = dict(f.split(":", 1) for f in line.split("\t")[1:] if ":" in f)
            if "SN" in fields and "LN" in fields:
                lengths[fields["SN"]] = int(fields["LN"])
    return lengths


def reference_contig_length(reference: Path, contig: str) -> int | None:
    fai = Path(f"{reference}.fai")
    if not fai.is_file():
        return None
    with open(fai) as handle:
        for line in handle:
            f = line.split("\t")
            if f[0] == contig:
                return int(f[1])
    return None


def reference_check(samtools: str, bam: Path, reference: Path, contig: str) -> tuple[bool, str]:
    """(usable for coordinate lookups, message). A BAM aligned to another assembly
    (e.g. GRCh38: chr15 = 101,991,189 bp; CHM13: 99,753,195 bp) is not usable."""
    lengths = header_contig_lengths(samtools, bam)
    if not lengths:
        return False, "no @SQ lines: the BAM is unaligned"
    if contig not in lengths:
        return False, f"no contig named {contig} in the BAM header (first: {next(iter(lengths))})"
    expected = reference_contig_length(reference, contig)
    if expected is None:
        return True, f"{contig} length {lengths[contig]:,} (reference .fai not found; not compared)"
    if lengths[contig] != expected:
        return False, (f"{contig} is {lengths[contig]:,} bp in the BAM but {expected:,} bp in the "
                       f"reference: aligned to a different assembly")
    return True, "same reference"


def has_methylation_tags(samtools: str, bam: Path, region: str, n_reads: int = 200) -> bool | None:
    """True if reads in the region carry MM/ML tags; None if it cannot be checked."""
    if bam_index(bam) is None:
        return None
    proc = subprocess.Popen([samtools, "view", str(bam), region], stdout=subprocess.PIPE, text=True)
    seen = found = 0
    try:
        for line in proc.stdout:
            seen += 1
            if "\tMM:Z:" in line or "\tMm:Z:" in line:
                found += 1
            if seen >= n_reads:
                break
    finally:
        proc.kill()
        proc.wait()
    return None if seen == 0 else found > 0


def fastq_command(samtools: str, bams: list[Path], threads: int, tags: str | None = None) -> str:
    """Shell fragment writing the primary reads of one or several BAMs to stdout."""
    tag = f" -T {tags}" if tags else ""
    parts = [f"{samtools} fastq -@ {threads}{tag} {b}" for b in bams]
    return parts[0] if len(parts) == 1 else "{ " + "; ".join(parts) + "; }"
