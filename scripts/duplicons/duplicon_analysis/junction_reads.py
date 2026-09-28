"""
Junction reads: HiFi reads that carry singly-unique k-mers (SUNKs) of BOTH copies
of one directly oriented SD pair, in the order expected for an NAHR fusion.

In a deletion carrier the deleted chromosome holds a hybrid duplicon: copy A up to
the crossover, copy B after it. A read anywhere across the crossover contains
A-specific SUNKs followed by B-specific SUNKs at increasing homologous offsets, so
the crossover lies between the last A-SUNK and the first B-SUNK on that read. The
precision is set by the spacing of paralog-specific variants, not by read depth.
Genomes without the fusion (the diploid panel) should yield no such reads, which
makes them a built-in negative control.
"""
from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .kmers import header_offset, kmer_strings_to_codes, read_fasta, seq_codes


@dataclass(frozen=True)
class SdPair:
    name: str
    chrom: str
    a0: int
    a1: int
    b0: int
    b1: int
    orientation: str


def read_direct_pairs(bedpe: Path, min_length: int = 5_000) -> list[SdPair]:
    """Unique, directly oriented SD pairs (copy A upstream of copy B)."""
    seen, pairs = set(), []
    with open(bedpe) as handle:
        lines = handle.readlines()
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        a = (int(f[1]), int(f[2]))
        b = (int(f[4]), int(f[5]))
        orient = "+" if (len(f) < 10 or f[8] == f[9]) else "-"
        if orient != "+" or min(a[1] - a[0], b[1] - b[0]) < min_length:
            continue
        a, b = sorted([a, b])
        if a[1] > b[0] or (a, b) in seen:          # overlapping copies are not a pair here
            continue
        seen.add((a, b))
        pairs.append(SdPair(f[6] if len(f) > 6 else f"pair{len(pairs)}", f[0], a[0], a[1], b[0], b[1], orient))
    return pairs


def sunk_index(window_fasta: Path, sunks_txt: Path, pairs: list[SdPair], k: int):
    """Sorted canonical codes of SUNKs inside the pair blocks with their
    (pair index, copy 0/1, homologous offset, reference position)."""
    name, seq = read_fasta(window_fasta)[0]
    chrom, offset = header_offset(name)
    codes, valid = seq_codes(seq, k)
    sunk_codes = np.unique(kmer_strings_to_codes(
        pd.read_csv(sunks_txt, sep=r"\s+", header=None, usecols=[0], names=["kmer"], dtype=str)["kmer"], k))
    rows = []
    for pi, pr in enumerate(pairs):
        for copy, (s, e) in enumerate(((pr.a0, pr.a1), (pr.b0, pr.b1))):
            lo, hi = s - offset, e - offset - k + 1
            idx = np.arange(max(lo, 0), max(min(hi, codes.size), 0))
            keep = valid[idx] & np.isin(codes[idx], sunk_codes)
            idx = idx[keep]
            ref = idx + offset
            rows.append(pd.DataFrame({"code": codes[idx], "pair": pi, "copy": copy,
                                      "offset": ref - s, "ref": ref}))
    table = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=["code", "pair", "copy", "offset", "ref"])
    table = table.drop_duplicates("code", keep=False).sort_values("code").reset_index(drop=True)
    return table


def fetch_reads(samtools: str, bam: Path, regions_bed: Path) -> dict[str, str]:
    """Primary alignments overlapping the SD blocks (any MAPQ): name -> sequence."""
    out = subprocess.run([samtools, "view", "-F", "0x904", "-L", str(regions_bed), str(bam)],
                         check=True, capture_output=True, text=True).stdout
    reads = {}
    for line in out.splitlines():
        f = line.split("\t", 10)
        if len(f) > 9 and f[9] != "*":
            reads[f[0]] = f[9]
    return reads


def classify_read(seq: str, index: pd.DataFrame, k: int, min_hits: int = 2):
    """Return per-pair junction evidence for one read, or []."""
    codes, valid = seq_codes(seq, k)
    pos = np.nonzero(valid)[0]
    q = codes[pos]
    ic = index["code"].to_numpy()
    j = np.searchsorted(ic, q)
    jc = np.minimum(j, ic.size - 1)
    hit = (j < ic.size) & (ic[jc] == q)
    if hit.sum() < 2 * min_hits:
        return []
    h = index.iloc[jc[hit]].assign(read_pos=pos[hit])
    calls = []
    for pair, g in h.groupby("pair"):
        a, b = g[g["copy"] == 0], g[g["copy"] == 1]
        if len(a) < min_hits or len(b) < min_hits:
            continue
        # orient the read along the reference using the larger copy's hits
        ref_side = a if len(a) >= len(b) else b
        forward = np.corrcoef(ref_side["read_pos"], ref_side["ref"])[0, 1] >= 0 if len(ref_side) > 2 else True
        rp = g["read_pos"] if forward else len(seq) - g["read_pos"]
        g = g.assign(rp=rp.to_numpy())
        a, b = g[g["copy"] == 0], g[g["copy"] == 1]
        ordered = a["rp"].max() < b["rp"].min()
        fused = a["offset"].max() < b["offset"].min()
        calls.append({"pair": int(pair), "a_hits": len(a), "b_hits": len(b), "ordered_A_then_B": bool(ordered),
                      "offsets_consistent": bool(fused),
                      "crossover_offset_low": int(a["offset"].max()), "crossover_offset_high": int(b["offset"].min())})
    return calls


def prepare(window_fasta: Path, sunks_txt: Path, sd_pairs_bedpe: Path, k: int = 31):
    """Build once per cohort: the direct SD pairs and their SUNK index."""
    pairs = read_direct_pairs(sd_pairs_bedpe)
    index = sunk_index(window_fasta, sunks_txt, pairs, k) if pairs else None
    return pairs, index


def junction_reads(sample: str, bam: Path | list[Path], samtools: str, pairs: list[SdPair],
                   index: pd.DataFrame, workdir: Path, k: int = 31):
    """bam: one BAM or several (e.g. two sequencing runs); reads are merged by name."""
    if not pairs or index is None or index.empty:
        return pd.DataFrame(), pd.DataFrame()
    regions = workdir / f"{sample}.sd_pair_blocks.bed"
    with open(regions, "w") as out:
        for pr in pairs:
            out.write(f"{pr.chrom}\t{pr.a0}\t{pr.a1}\n{pr.chrom}\t{pr.b0}\t{pr.b1}\n")
    reads = {}
    for one in ([bam] if isinstance(bam, (str, Path)) else bam):
        reads.update(fetch_reads(samtools, Path(one), regions))
    rows = []
    for name, seq in reads.items():
        for call in classify_read(seq, index, k):
            rows.append({"sample": sample, "read": name, **call})
    per_read = pd.DataFrame(rows)
    summary = []
    for pi, pr in enumerate(pairs):
        sub = per_read[per_read["pair"] == pi] if len(per_read) else per_read
        good = sub[sub["ordered_A_then_B"] & sub["offsets_consistent"]] if len(sub) else sub
        row = {"sample": sample, "pair": pr.name, "copy_A": f"{pr.chrom}:{pr.a0}-{pr.a1}",
               "copy_B": f"{pr.chrom}:{pr.b0}-{pr.b1}", "reads_in_blocks": len(reads),
               "junction_reads": len(good), "reads_both_copies_other_order": len(sub) - len(good)}
        if len(good):
            lo, hi = int(good["crossover_offset_low"].max()), int(good["crossover_offset_high"].min())
            if lo >= hi:        # reads disagree: report their spread instead of an empty intersection
                lo, hi = int(good["crossover_offset_low"].median()), int(good["crossover_offset_high"].median())
            row.update({"crossover_offset_low": lo, "crossover_offset_high": hi,
                        "crossover_in_A": f"{pr.chrom}:{pr.a0 + lo}-{pr.a0 + hi}",
                        "crossover_in_B": f"{pr.chrom}:{pr.b0 + lo}-{pr.b0 + hi}"})
        summary.append(row)
    regions.unlink(missing_ok=True)
    summary = pd.DataFrame(summary)
    for col in ("crossover_offset_low", "crossover_offset_high"):
        if col in summary:
            summary[col] = summary[col].astype("Int64")
    return per_read, summary
