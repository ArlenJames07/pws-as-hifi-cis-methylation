"""Per-molecule CpG methylation from MM/ML-tagged HiFi BAMs (HiPhase output keeps
HP/PS). Needs pysam.

Every call is placed on the C of the CpG on the forward reference strand: a call on
a reverse-strand read sits on the G of the CpG (ref position + 1), so it is moved one
base left. Probabilities are ML/255."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .readers import sample_files

MOD_CODES = {"m", "5mC", "c+m", 27551}      # 5mC (27551 = ChEBI code)


def pysam_or_exit():
    try:
        import pysam
    except ImportError as error:
        raise SystemExit("pysam is required for reading modified-base BAMs (conda env envs/figures.yml "
                         "or `pip install pysam`)") from error
    return pysam


@dataclass
class Molecule:
    read_id: str
    sample: str
    hp: int            # HiPhase haplotype (0 = unphased)
    ps: int            # phase block (-1 = none)
    start: int         # alignment start (0-based)
    end: int
    reverse: bool
    positions: np.ndarray   # CpG C positions, forward strand, sorted
    probs: np.ndarray       # methylation probability 0-1

    def within(self, start: int, end: int) -> tuple[np.ndarray, np.ndarray]:
        m = (self.positions >= start) & (self.positions < end)
        return self.positions[m], self.probs[m]


def find_bam(paths, sample: str) -> tuple[Path | None, bool]:
    """(BAM, phased). HiPhase BAM in 04_phasing first; else the alignment BAM (no HP/PS)."""
    for pattern in (f"{sample}.phased.bam", "*.phased.bam", "*.haplotagged.bam"):
        hits = [p for p in sample_files(paths.phasing, sample, (pattern,)) if Path(f"{p}.bai").is_file()
                or Path(str(p)[:-4] + ".bai").is_file() or Path(f"{p}.csi").is_file()]
        if hits:
            return hits[0], True
    hits = [p for p in sample_files(paths.alignment, sample, ("*.bam",))
            if Path(f"{p}.bai").is_file() or Path(f"{p}.csi").is_file() or Path(str(p)[:-4] + ".bai").is_file()]
    return (hits[0], False) if hits else (None, False)


def molecules(bam_path: Path, sample: str, chrom: str, start: int, end: int, min_mapq: int = 20,
              min_cpgs: int = 1) -> list[Molecule]:
    pysam = pysam_or_exit()
    out = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(chrom, max(0, start), end):
            if (read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate
                    or read.is_qcfail or read.mapping_quality < min_mapq):
                continue
            mods = read.modified_bases
            if not mods:
                continue
            calls = []
            for (base, _strand, code), items in mods.items():
                if base not in ("C", "G") or code not in MOD_CODES:
                    continue
                calls.extend(items)
            if not calls:
                continue
            q2r = np.full(read.query_length or len(read.query_sequence or ""), -1, dtype=np.int64)
            pairs = np.array(read.get_aligned_pairs(matches_only=True), dtype=np.int64)
            if pairs.size == 0:
                continue
            q2r[pairs[:, 0]] = pairs[:, 1]
            q = np.array([c[0] for c in calls], dtype=np.int64)
            prob = np.array([c[1] for c in calls], dtype=float) / 255.0
            keep = (q >= 0) & (q < len(q2r))
            q, prob = q[keep], prob[keep]
            ref = q2r[q]
            ok = ref >= 0
            ref, prob = ref[ok], prob[ok]
            if read.is_reverse:
                ref = ref - 1
            order = np.argsort(ref)
            ref, prob = ref[order], prob[order]
            if len(ref) > 1:            # one call per CpG
                uniq, idx = np.unique(ref, return_index=True)
                ref, prob = uniq, prob[idx]
            if ((ref >= start) & (ref < end)).sum() < min_cpgs:
                continue
            hp = read.get_tag("HP") if read.has_tag("HP") else 0
            ps = read.get_tag("PS") if read.has_tag("PS") else -1
            out.append(Molecule(read.query_name, sample, int(hp), int(ps), read.reference_start,
                                read.reference_end, read.is_reverse, ref, prob.astype(np.float32)))
    return out


def binary_entropy(fraction: np.ndarray) -> np.ndarray:
    """Shannon entropy (bits) of a molecule's methylated fraction: 0 for a molecule that is
    fully methylated or unmethylated, 1 for half-and-half (maximal discordance)."""
    p = np.clip(np.asarray(fraction, float), 1e-12, 1 - 1e-12)
    h = -(p * np.log2(p) + (1 - p) * np.log2(1 - p))
    return np.where((np.asarray(fraction) <= 0) | (np.asarray(fraction) >= 1), 0.0, h)


def switch_rate(probs: np.ndarray) -> float:
    """Fraction of adjacent CpG pairs whose binary state differs (along the molecule)."""
    s = np.asarray(probs) >= 0.5
    return float(np.mean(s[1:] != s[:-1])) if len(s) > 1 else np.nan
