"""
k-mer utilities for paralog-specific copy number at singly-unique k-mers (SUNKs;
Sudmant et al. 2010). k-mers are 2-bit encoded in uint64 (k <= 32) and compared
in a canonical form computed here, so the canonical convention of the counting
tool (meryl, jellyfish) does not matter.
"""
from __future__ import annotations

import math
import re

import numpy as np
import pandas as pd

LUT = np.full(256, 4, dtype=np.uint8)
for ch, v in zip("ACGTacgt", (0, 1, 2, 3, 0, 1, 2, 3)):
    LUT[ord(ch)] = v


def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def header_offset(name):
    """chr15:19500001-31500000 (samtools, 1-based) -> ('chr15', 19500000)"""
    m = re.match(r"^(.+):(\d+)-(\d+)$", name)
    if not m:
        return name, 0
    return m.group(1), int(m.group(2)) - 1


def seq_codes(seq: str, k: int):
    """Canonical 2-bit codes of every k-mer in seq (uint64) and a validity mask."""
    a = LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]
    n = len(a) - k + 1
    if n <= 0:
        return np.zeros(0, np.uint64), np.zeros(0, bool)
    bad = (a == 4).astype(np.int32)
    cs = np.concatenate(([0], np.cumsum(bad)))
    valid = (cs[k:k + n] - cs[:n]) == 0
    a = np.where(a == 4, 0, a).astype(np.uint64)
    comp = (3 - a).astype(np.uint64)
    fwd = np.zeros(n, np.uint64)
    rev = np.zeros(n, np.uint64)
    two = np.uint64(2)
    for j in range(k):
        fwd = (fwd << two) | a[j:j + n]
        rev = (rev << two) | comp[k - 1 - j:k - 1 - j + n]
    return np.minimum(fwd, rev), valid


def kmer_strings_to_codes(kmers: pd.Series, k: int):
    s = kmers.str.upper().to_numpy(dtype=f"S{k}")
    if s.size == 0:
        return np.zeros(0, np.uint64)
    m = LUT[np.frombuffer(s.tobytes(), dtype=np.uint8).reshape(-1, k)]
    if (m == 4).any():
        raise ValueError("non-ACGT k-mer in counts file")
    m = m.astype(np.uint64)
    fwd = np.zeros(len(s), np.uint64)
    rev = np.zeros(len(s), np.uint64)
    two = np.uint64(2)
    for j in range(k):
        fwd = (fwd << two) | m[:, j]
        rev = (rev << two) | (np.uint64(3) - m[:, k - 1 - j])
    return np.minimum(fwd, rev)


def read_counts(path, k):
    df = pd.read_csv(path, sep=r"\s+", header=None, usecols=[0, 1],
                     names=["kmer", "count"], dtype={"kmer": str, "count": np.int64},
                     engine="c")
    df = df[df["kmer"].str.len() == k]
    codes = kmer_strings_to_codes(df["kmer"], k)
    order = np.argsort(codes, kind="stable")
    return codes[order], df["count"].to_numpy()[order]


def lookup(sorted_codes, sorted_vals, query):
    if sorted_codes.size == 0:
        return np.zeros(query.size, np.int64)
    idx = np.searchsorted(sorted_codes, query)
    idx_c = np.minimum(idx, sorted_codes.size - 1)
    found = (idx < sorted_codes.size) & (sorted_codes[idx_c] == query)
    return np.where(found, sorted_vals[idx_c], 0)


def viterbi(x, sd, states, p_switch, prior_state=2):
    """x: observed CN per bin; sd: per-bin sd for a CN-2 bin. Gaussian emissions,
    sd scaled by sqrt(state/2) (Poisson-like), floored."""
    S, n = len(states), len(x)
    logA = np.full((S, S), math.log(p_switch / (S - 1)))
    np.fill_diagonal(logA, math.log(1 - p_switch))
    init = np.full(S, math.log(0.1 / (S - 1)))
    init[states.index(prior_state)] = math.log(0.9)
    mu = np.array(states, float)
    scale = np.sqrt(np.maximum(mu, 0.5) / 2.0)
    sdm = np.maximum(np.outer(sd, scale), 0.05)                      # n x S
    E = -0.5 * ((x[:, None] - mu[None, :]) / sdm) ** 2 - np.log(sdm)  # n x S
    V = np.zeros((n, S))
    B = np.zeros((n, S), int)
    V[0] = init + E[0]
    for t in range(1, n):
        cand = V[t - 1][:, None] + logA
        B[t] = cand.argmax(0)
        V[t] = cand.max(0) + E[t]
    path = np.zeros(n, int)
    path[-1] = V[-1].argmax()
    for t in range(n - 1, 0, -1):
        path[t - 1] = B[t, path[t]]
    return np.array(states)[path]
