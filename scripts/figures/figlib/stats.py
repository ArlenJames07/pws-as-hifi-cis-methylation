"""Tests used by the figures, in numpy only (the duplicon environment has no scipy).

With 3-5 genomes per group, every test that decides a claim works on per-genome
values (exact label permutations); molecule- or CpG-level tests are reported as
descriptive because molecules of one genome are not independent."""
from __future__ import annotations

import itertools
import math

import numpy as np


def bh(pvalues) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    q = p[ok]
    order = np.argsort(q)
    ranked = q[order] * len(q) / (np.arange(len(q)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty_like(q)
    res[order] = np.clip(ranked, 0, 1)
    out[ok] = res
    return out


def rankdata(x) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x))
    sx = x[order]
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and sx[j + 1] == sx[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2 + 1
        i = j + 1
    return ranks


def _norm_sf(z: float) -> float:
    return 0.5 * math.erfc(z / math.sqrt(2))


def mann_whitney(a, b) -> tuple[float, float]:
    """U of `a` and a two-sided p: exact over all splits when feasible, else normal with
    tie and continuity corrections."""
    a = np.asarray(a, float)[np.isfinite(a)]
    b = np.asarray(b, float)[np.isfinite(b)]
    na, nb = len(a), len(b)
    if na == 0 or nb == 0:
        return np.nan, np.nan
    r = rankdata(np.concatenate([a, b]))
    u = r[:na].sum() - na * (na + 1) / 2
    n = na + nb
    if math.comb(n, na) <= 100_000:
        mu = na * nb / 2
        count = total = 0
        for idx in itertools.combinations(range(n), na):
            uu = r[list(idx)].sum() - na * (na + 1) / 2
            total += 1
            count += abs(uu - mu) >= abs(u - mu) - 1e-9
        return float(u), count / total
    _, counts = np.unique(r, return_counts=True)
    tie = (counts ** 3 - counts).sum() / (n * (n - 1))
    sigma = math.sqrt(na * nb / 12 * ((n + 1) - tie))
    if sigma == 0:
        return float(u), 1.0
    z = (abs(u - na * nb / 2) - 0.5) / sigma
    return float(u), float(min(1.0, 2 * _norm_sf(max(z, 0))))


def rank_biserial(u: float, na: int, nb: int) -> float:
    """> 0 when values of the first group tend to be larger."""
    return 2 * u / (na * nb) - 1 if na and nb else np.nan


def _gammq(s: float, x: float) -> float:
    """Regularized upper incomplete gamma Q(s, x) (series / continued fraction)."""
    if x <= 0:
        return 1.0
    gln = math.lgamma(s)
    if x < s + 1:
        ap, total, delta = s, 1 / s, 1 / s
        for _ in range(1000):
            ap += 1
            delta *= x / ap
            total += delta
            if abs(delta) < abs(total) * 1e-14:
                break
        return max(0.0, 1 - total * math.exp(-x + s * math.log(x) - gln))
    b = x + 1 - s
    c, d = 1 / 1e-300, 1 / b
    h = d
    for i in range(1, 1000):
        an = -i * (i - s)
        b += 2
        d = an * d + b
        d = 1e-300 if abs(d) < 1e-300 else d
        c = b + an / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1 / d
        h *= d * c
        if abs(d * c - 1) < 1e-14:
            break
    return math.exp(-x + s * math.log(x) - gln) * h


def chi2_sf(x: float, df: int) -> float:
    return _gammq(df / 2, x / 2) if np.isfinite(x) else np.nan


def kruskal(groups) -> tuple[float, float, float]:
    """H (tie-corrected), chi-square p, and eta-squared_H = (H - k + 1) / (n - k)."""
    groups = [np.asarray(g, float)[np.isfinite(g)] for g in groups]
    groups = [g for g in groups if len(g)]
    k = len(groups)
    n = sum(len(g) for g in groups)
    if k < 2 or n <= k:
        return np.nan, np.nan, np.nan
    r = rankdata(np.concatenate(groups))
    h, i = 0.0, 0
    for g in groups:
        h += r[i:i + len(g)].sum() ** 2 / len(g)
        i += len(g)
    h = 12 / (n * (n + 1)) * h - 3 * (n + 1)
    _, counts = np.unique(r, return_counts=True)
    corr = 1 - (counts ** 3 - counts).sum() / (n ** 3 - n)
    h = h / corr if corr > 0 else np.nan
    return float(h), chi2_sf(h, k - 1), max(0.0, (h - k + 1) / (n - k))


def label_permutation(a, b, stat: str = "mean", max_exact: int = 200_000, n_mc: int = 20_000,
                      seed: int = 1) -> tuple[float, float, int]:
    """Difference (a - b) of means or medians and its two-sided permutation p over all
    relabellings of the pooled values (exact when feasible, else Monte Carlo)."""
    a = np.asarray(a, float)[np.isfinite(a)]
    b = np.asarray(b, float)[np.isfinite(b)]
    if len(a) < 1 or len(b) < 1:
        return np.nan, np.nan, 0
    f = np.mean if stat == "mean" else np.median
    vals = np.concatenate([a, b])
    n, na = len(vals), len(a)
    obs = f(a) - f(b)
    if math.comb(n, na) <= max_exact:
        stats = []
        for idx in itertools.combinations(range(n), na):
            m = np.zeros(n, bool)
            m[list(idx)] = True
            stats.append(f(vals[m]) - f(vals[~m]))
        stats = np.abs(np.asarray(stats))
        return float(obs), float(np.mean(stats >= abs(obs) - 1e-12)), len(stats)
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(n_mc):
        p = rng.permutation(vals)
        ge += abs(f(p[:na]) - f(p[na:])) >= abs(obs) - 1e-12
    return float(obs), (ge + 1) / (n_mc + 1), n_mc


def kruskal_permutation(groups, n_mc: int = 20_000, seed: int = 1729) -> float:
    """Permutation p of Kruskal-Wallis H (Monte Carlo over relabellings)."""
    groups = [np.asarray(g, float)[np.isfinite(g)] for g in groups]
    groups = [g for g in groups if len(g)]
    h_obs = kruskal(groups)[0]
    if not np.isfinite(h_obs):
        return np.nan
    vals = np.concatenate(groups)
    sizes = np.cumsum([len(g) for g in groups])[:-1]
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(n_mc):
        ge += kruskal(np.split(rng.permutation(vals), sizes))[0] >= h_obs - 1e-12
    return (ge + 1) / (n_mc + 1)


def bootstrap_ci(values, func=np.mean, n: int = 5_000, seed: int = 13, alpha: float = 0.05):
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if len(v) < 2:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    boots = np.array([func(v[rng.integers(0, len(v), len(v))]) for _ in range(n)])
    return float(np.quantile(boots, alpha / 2)), float(np.quantile(boots, 1 - alpha / 2))


def bootstrap_diff_ci(a, b, func=np.median, n: int = 5_000, seed: int = 20260525, alpha: float = 0.05):
    a = np.asarray(a, float)[np.isfinite(a)]
    b = np.asarray(b, float)[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    d = [func(a[rng.integers(0, len(a), len(a))]) - func(b[rng.integers(0, len(b), len(b))]) for _ in range(n)]
    return float(np.quantile(d, alpha / 2)), float(np.quantile(d, 1 - alpha / 2))


def sign_flip_p(values, max_exact: int = 16, n_mc: int = 20_000, seed: int = 1) -> float:
    """Two-sided sign-flip p for mean(values) = 0."""
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return np.nan
    obs = abs(v.mean())
    if len(v) <= max_exact:
        signs = np.array(list(itertools.product((-1, 1), repeat=len(v))))
        return float(np.mean(np.abs((signs * v).mean(axis=1)) >= obs - 1e-12))
    rng = np.random.default_rng(seed)
    flips = rng.choice((-1, 1), size=(n_mc, len(v)))
    return float((np.sum(np.abs((flips * v).mean(axis=1)) >= obs - 1e-12) + 1) / (n_mc + 1))
