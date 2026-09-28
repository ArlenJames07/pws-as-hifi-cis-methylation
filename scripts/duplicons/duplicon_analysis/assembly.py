"""
Assembly helpers: minimap2 PAF blocks (contigs or haplotype pieces aligned to
T2T-CHM13) with exact coordinate projection through the cg:Z: CIGAR, plus small
interval utilities. Used by scripts/duplicons/03-07.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

CG_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_region(text: str) -> tuple[str, int, int] | None:
    """'chr15:22,000,000-28,000,000' -> ('chr15', 22000000, 28000000); '.'/'' -> None."""
    if text in (None, "", ".", "NA"):
        return None
    match = re.match(r"^([^:]+):([\d,]+)-([\d,]+)$", str(text))
    if not match:
        raise ValueError(f"bad region '{text}', expected chrom:start-end")
    return match.group(1), int(match.group(2).replace(",", "")), int(match.group(3).replace(",", ""))


def window_offset(name: str) -> tuple[str, int]:
    """Record name written by `samtools faidx ref chr15:S-E` -> ('chr15', S-1)."""
    match = re.match(r"^(.+):(\d+)-(\d+)$", name)
    return (match.group(1), int(match.group(2)) - 1) if match else (name, 0)


def read_bed_intervals(path: Path | str | None, chrom: str) -> list[tuple[int, int, str]]:
    out: list[tuple[int, int, str]] = []
    if not path or not Path(path).is_file():
        return out
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[0] == chrom:
                name = fields[3] if len(fields) > 3 and fields[3] else f"{fields[1]}-{fields[2]}"
                out.append((int(fields[1]), int(fields[2]), name))
    return sorted(out)


def names_at(intervals: list[tuple[int, int, str]], pos: int) -> list[str]:
    return [name for start, end, name in intervals if start <= pos < end]


def names_overlapping(intervals: list[tuple[int, int, str]], a: int, b: int) -> list[str]:
    return sorted({name for start, end, name in intervals if start < b and end > a})


def merge_intervals(intervals) -> list[list[int]]:
    out: list[list[int]] = []
    for start, end in sorted(intervals):
        if out and start <= out[-1][1]:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([start, end])
    return out


def covered_bp(merged: list[list[int]], a: int, b: int) -> int:
    return sum(max(0, min(end, b) - max(start, a)) for start, end in merged)


class Block:
    __slots__ = ("hap", "q", "qlen", "qs", "qe", "strand", "t", "tlen", "ts", "te",
                 "nmatch", "alen", "mapq", "tp", "cg", "de")

    def __init__(self, hap, f):
        self.hap = hap
        self.q, self.qlen, self.qs, self.qe = f[0], int(f[1]), int(f[2]), int(f[3])
        self.strand = f[4]
        self.t, self.tlen, self.ts, self.te = f[5], int(f[6]), int(f[7]), int(f[8])
        self.nmatch, self.alen, self.mapq = int(f[9]), int(f[10]), int(f[11])
        tags = {}
        for x in f[12:]:
            p = x.split(":", 2)
            if len(p) == 3:
                tags[p[0]] = p[2]
        self.tp = tags.get("tp", "P")
        self.cg = tags.get("cg")
        self.de = float(tags["de"]) if "de" in tags else None

    @property
    def identity(self):
        return self.nmatch / self.alen if self.alen else 0.0

    def t_to_q(self, tpos: int) -> int:
        """Project a reference coordinate inside [ts, te] onto contig coordinates.
        Uses the cg:Z: CIGAR when present, linear interpolation otherwise."""
        tpos = min(max(tpos, self.ts), self.te)
        if not self.cg:
            frac = (tpos - self.ts) / max(1, self.te - self.ts)
            off = round(frac * (self.qe - self.qs))
        else:
            t, off = self.ts, 0
            for n, op in CG_RE.findall(self.cg):
                n = int(n)
                if op in "M=X":
                    if t + n >= tpos:
                        off += tpos - t
                        t = tpos
                        break
                    t += n
                    off += n
                elif op == "I":
                    off += n
                elif op in "DN":
                    if t + n >= tpos:
                        t = tpos
                        break
                    t += n
            else:
                off = self.qe - self.qs
        return self.qs + off if self.strand == "+" else self.qe - off

    # -- helpers used by contig_origin.py and repeat_content.py -----------------
    def _segments(self):
        """Matched segments as arrays (t_start, q_offset_start, length); q_offset is
        measured along the alignment in reference orientation (0 = qs for '+',
        0 = qe-1 for '-')."""
        ts, qo, ln = [], [], []
        t, q = self.ts, 0
        ops = CG_RE.findall(self.cg) if self.cg else [(str(self.te - self.ts), "M")]
        for n, op in ops:
            n = int(n)
            if op in "M=X":
                ts.append(t); qo.append(q); ln.append(n)
                t += n; q += n
            elif op == "I":
                q += n
            elif op in "DN":
                t += n
        return ts, qo, ln

    def q_to_t_array(self, qpos, cpg=False):
        """Vectorized contig -> reference projection. Returns -1 where the base is
        outside the block or inside an insertion. With cpg=True on a '-' block, the
        CpG starting at contig q is reported at the reference C (t of q+1)."""
        qpos = np.asarray(qpos, dtype=np.int64)
        if self.strand == "-" and cpg:
            qpos = qpos + 1
        off = qpos - self.qs if self.strand == "+" else (self.qe - 1) - qpos
        ts, qo, ln = (np.asarray(x, dtype=np.int64) for x in self._segments())
        out = np.full(qpos.shape, -1, dtype=np.int64)
        if ts.size == 0:
            return out
        i = np.searchsorted(qo, off, side="right") - 1
        ok = (i >= 0) & (off >= 0)
        i_c = np.clip(i, 0, ts.size - 1)
        ok &= off < qo[i_c] + ln[i_c]
        out[ok] = ts[i_c[ok]] + (off[ok] - qo[i_c[ok]])
        return out

    def indels(self, min_len=50):
        """Yield ('INS'|'DEL', ref_pos, contig_start, contig_end, length) from cg."""
        if not self.cg:
            return
        t, q = self.ts, 0
        for n, op in CG_RE.findall(self.cg):
            n = int(n)
            if op in "M=X":
                t += n; q += n
            elif op == "I":
                if n >= min_len:
                    a, b = (self.qs + q, self.qs + q + n) if self.strand == "+" else \
                           (self.qe - q - n, self.qe - q)
                    yield ("INS", t, a, b, n)
                q += n
            elif op in "DN":
                if n >= min_len:
                    c = self.qs + q if self.strand == "+" else self.qe - q
                    yield ("DEL", t, c, c, n)
                t += n


def load_paf(path, hap, min_len):
    """Keep minimap2 tp:A:P (primary, which includes supplementary chains in PAF) and
    tp:A:I (primary inversion); drop tp:A:S / tp:A:i (secondary)."""
    blocks = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            b = Block(hap, f)
            if b.tp not in ("P", "I"):
                continue
            if b.alen < min_len:
                continue
            blocks.append(b)
    return blocks
