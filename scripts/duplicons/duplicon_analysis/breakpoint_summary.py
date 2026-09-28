"""
Per-carrier breakpoint summary from the outputs of step 05 (no reads are re-read).

Why this exists (lessons from the real cohort)
  * The SUNK HMM produces many copy-number switches inside the segmental duplications
    of every genome, controls included: duplicon copy number and paralog-specific
    variants are polymorphic in the population. Matching "the nearest 2->1 switch" to a
    HiFiCNV edge therefore picks unrelated switches. Here the deletion is taken as the
    CN=1 run that overlaps the HiFiCNV interval, runs of CN=1 separated only by bins
    lying inside the HiFiCNV interval are joined, and its two outer edges are the
    breakpoints.
  * Reads carrying copy-A SUNKs followed by copy-B SUNKs of a pair ("junction-like"
    reads) occur in normal chromosomes too (gene conversion, paralog variants that
    differ from CHM13). Their presence is not evidence of a deletion; an EXCESS over the
    rate seen in genomes with two copies of chr15, at a pair whose copies lie at the
    carrier's two edges, is.

Outputs (in the breakpoints folder)
  breakpoints_vs_hificnv.tsv  one row per carrier and edge: HiFiCNV edge, SUNK edge,
                              support interval, shift, SD block and BP cluster
  nahr_test.tsv               one row per carrier: do the two edges fall at homologous
                              positions of one directly oriented SD pair?
  junction_enrichment.tsv     carrier x pair: junction-like reads observed vs expected
                              from the panel rate, Poisson p, whether the pair spans
                              the carrier's edges
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from .junction_reads import read_direct_pairs


def _poisson_sf(k: int, mu: float) -> float:
    """P(X >= k) for X ~ Poisson(mu)."""
    if k <= 0:
        return 1.0
    if mu <= 0:
        return 0.0
    term, cdf = math.exp(-mu), 0.0
    for i in range(k):
        cdf += term
        term *= mu / (i + 1)
    return max(0.0, 1.0 - cdf)


def _runs(states: np.ndarray, value: int) -> list[tuple[int, int]]:
    """Index ranges [i, j] (inclusive) of consecutive entries equal to value."""
    runs, start = [], None
    for i, s in enumerate(states):
        if s == value and start is None:
            start = i
        elif s != value and start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, len(states) - 1))
    return runs


def deletion_run(bins: pd.DataFrame, h0: int, h1: int):
    """Outer edges of the CN=1 run that carries the deletion. Returns
    (first_bin_row, last_bin_row, at_left_window_edge, at_right_window_edge) or None."""
    b = bins[bins["usable"]].sort_values("start").reset_index(drop=True)
    if b.empty:
        return None
    runs = _runs(b["state"].to_numpy(), 1)
    if not runs:
        return None

    def overlap(r):
        s, e = b.loc[r[0], "start"], b.loc[r[1], "end"]
        return max(0, min(e, h1) - max(s, h0))

    best = max(runs, key=overlap)
    if overlap(best) == 0:
        return None
    i, j = best
    # join CN=1 runs whose separating bins all lie inside the HiFiCNV interval
    changed = True
    while changed:
        changed = False
        for r in runs:
            if r[1] < i and b.loc[r[1] + 1, "start"] >= h0 and b.loc[i - 1, "end"] <= h1:
                i, changed = r[0], True
            if r[0] > j and b.loc[j + 1, "start"] >= h0 and b.loc[r[0] - 1, "end"] <= h1:
                j, changed = r[1], True
    return b.loc[i], b.loc[j], i == 0, j == len(b) - 1


def _nearest_transition(tr: pd.DataFrame, pos: int, into_one: bool, max_bp: int):
    if tr.empty:
        return None
    cand = tr[tr["to_cn"] == 1] if into_one else tr[tr["from_cn"] == 1]
    if cand.empty:
        return None
    d = (cand["bin_boundary"] - pos).abs()
    k = d.idxmin()
    return cand.loc[k] if d.loc[k] <= max_bp else None


def _block(sd_blocks: list[tuple[int, int]], pos: int, slack: int = 5_000) -> str:
    for s, e in sd_blocks:
        if s - slack <= pos < e + slack:
            return f"{s}-{e}"
    return "."


def _cluster(clusters: dict[str, tuple[int, int]], pos: int, slack: int = 50_000) -> str:
    for name, (s, e) in clusters.items():
        if s - slack <= pos <= e + slack:
            return name
    return "."


def add_position_specific_junctions(outdir: Path, nahr_rows: list[dict], pairs, deletion_mechanisms,
                                    panel_mechanisms, min_window: int = 10_000) -> None:
    """For each carrier whose edges sit in one direct pair, count junction-like reads of
    that pair whose crossover interval overlaps the carrier's own SUNK crossover
    (+/- max(10 kb, support)), and compare with the same position in the panel genomes.
    A real fusion adds reads exactly there; paralog variation spreads them along the pair."""
    per_read, crossovers = outdir / "junction_reads.tsv", outdir / "junction_crossovers.tsv"
    if not per_read.is_file() or per_read.stat().st_size == 0 or not crossovers.is_file():
        return
    reads = pd.read_csv(per_read, sep="\t")
    summary = pd.read_csv(crossovers, sep="\t")
    if reads.empty or summary.empty:
        return
    reads = reads[reads["ordered_A_then_B"].astype(bool) & reads["offsets_consistent"].astype(bool)]
    depth = summary.groupby("sample")["reads_in_blocks"].first()
    mech = summary.groupby("sample")["mechanism"].first()
    panel = [s for s in depth.index if mech[s] in panel_mechanisms]
    index = {pr.name: k for k, pr in enumerate(pairs)}
    for row in nahr_rows:
        if row.get("pair", ".") in (".", None) or row.get("pair") not in index:
            continue
        k, pr = index[row["pair"]], pairs[index[row["pair"]]]
        offset = row["proximal"] - pr.a0
        w = max(min_window, int(row.get("combined_support_bp", 0) or 0))
        lo, hi = offset - w, offset + w

        def count(sample):
            r = reads[(reads["sample"] == sample) & (reads["pair"] == k)]
            return int(((r["crossover_offset_high"] >= lo) & (r["crossover_offset_low"] <= hi)).sum())

        obs = count(row["sample_id"])
        panel_reads = sum(count(s) for s in panel)
        panel_depth = sum(depth[s] for s in panel)
        mu = panel_reads / panel_depth * depth.get(row["sample_id"], 0) if panel_depth else np.nan
        row.update({"junction_reads_at_crossover": obs, "expected_from_panel_at_crossover": mu,
                    "panel_genomes_with_reads_at_crossover": int(sum(count(s) > 0 for s in panel)),
                    "poisson_p_at_crossover": _poisson_sf(obs, mu) if mu == mu else np.nan})


def summarise(outdir: Path, structural_path: Path, sd_bed: Path, sd_pairs_bedpe: Path,
              deletion_mechanisms: tuple[str, ...], panel_mechanisms: tuple[str, ...],
              bp_clusters: dict[str, tuple[int, int]], bin_bp: int = 5_000):
    bins = pd.read_csv(outdir / "sunk15q.bins.tsv", sep="\t")
    try:
        tr = pd.read_csv(outdir / "sunk15q.transitions.tsv", sep="\t")
    except (pd.errors.EmptyDataError, FileNotFoundError):
        tr = pd.DataFrame(columns=["sample", "from_cn", "to_cn", "bin_boundary", "refined_pos",
                                   "support_start", "support_end"])
    structural = pd.read_csv(structural_path, sep="\t")
    structural = structural[structural["molecular_mechanism"].isin(deletion_mechanisms)]
    sd_blocks = []
    with open(sd_bed) as handle:
        for line in handle:
            f = line.split("\t")
            if len(f) >= 3:
                sd_blocks.append((int(f[1]), int(f[2])))
    pairs = read_direct_pairs(sd_pairs_bedpe)

    edge_rows, nahr_rows, edges = [], [], {}
    for rec in structural.itertuples(index=False):
        s = rec.sample_id
        if pd.isna(rec.cn_event_start):
            continue
        h0, h1 = int(rec.cn_event_start), int(rec.cn_event_end)
        found = deletion_run(bins[bins["sample"] == s], h0, h1)
        t = tr[tr["sample"] == s] if len(tr) else tr
        if found is None:
            for edge, h in (("proximal", h0), ("distal", h1)):
                edge_rows.append({"sample_id": s, "mechanism": rec.molecular_mechanism, "edge": edge,
                                  "hificnv_pos": h, "note": "no CN=1 run overlapping the HiFiCNV interval"})
            continue
        first, last, left_edge, right_edge = found
        info = {}
        for edge, h, bin_pos, into_one, at_window in (
                ("proximal", h0, int(first["start"]), True, left_edge),
                ("distal", h1, int(last["end"]), False, right_edge)):
            row = {"sample_id": s, "mechanism": rec.molecular_mechanism, "edge": edge, "hificnv_pos": h}
            if at_window:
                row.update({"sunk_pos": np.nan, "note": "CN=1 reaches the end of WINDOW: breakpoint outside"})
                edge_rows.append(row)
                continue
            best = _nearest_transition(t, bin_pos, into_one, 2 * bin_bp)
            if best is not None:
                pos, lo, hi = int(best["refined_pos"]), int(best["support_start"]), int(best["support_end"])
                other = int(best["from_cn"] if into_one else best["to_cn"])
            else:
                pos, lo, hi, other = bin_pos, bin_pos - bin_bp, bin_pos + bin_bp, np.nan
            row.update({"sunk_pos": pos, "support_start": lo, "support_end": hi, "shift_bp": pos - h,
                        "cn_outside": other, "sd_block": _block(sd_blocks, pos),
                        "bp_cluster": _cluster(bp_clusters, pos), "note": ""})
            if other == other and other != 2:
                row["note"] = f"flank reads CN {int(other)}, not 2 (duplicon polymorphism?)"
            info[edge] = (pos, lo, hi)
            edge_rows.append(row)
        edges[s] = info
        # NAHR: proximal edge in copy A, distal edge in copy B of one direct pair
        nrow = {"sample_id": s, "mechanism": rec.molecular_mechanism}
        if len(info) == 2:
            (p, plo, phi), (d, dlo, dhi) = info["proximal"], info["distal"]
            tol = 50_000
            best_pair = None
            for pr in pairs:
                if pr.a0 - tol <= p <= pr.a1 + tol and pr.b0 - tol <= d <= pr.b1 + tol:
                    frac = (min(max(p, pr.a0), pr.a1) - pr.a0) / max(1, pr.a1 - pr.a0)
                    hom = int(round(pr.b0 + frac * (pr.b1 - pr.b0)))
                    dist = abs(hom - d)
                    if best_pair is None or dist < best_pair[1]:
                        best_pair = (pr, dist, hom)
            if best_pair is not None:
                pr, dist, hom = best_pair
                support = (phi - plo) + (dhi - dlo)
                nrow.update({"proximal": p, "distal": d, "pair": pr.name,
                             "copy_A": f"{pr.chrom}:{pr.a0}-{pr.a1}", "copy_B": f"{pr.chrom}:{pr.b0}-{pr.b1}",
                             "homolog_of_proximal": hom, "distance_bp": dist, "combined_support_bp": support,
                             "consistent_with_NAHR": bool(dist <= max(20_000, support))})
            else:
                nrow.update({"proximal": p, "distal": d, "pair": ".",
                             "note": "no directly oriented SD pair has one copy at each edge"})
        else:
            nrow["note"] = "an edge lies outside WINDOW"
        nahr_rows.append(nrow)

    add_position_specific_junctions(outdir, nahr_rows, pairs, deletion_mechanisms, panel_mechanisms)
    pd.DataFrame(edge_rows).to_csv(outdir / "breakpoints_vs_hificnv.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(nahr_rows).to_csv(outdir / "nahr_test.tsv", sep="\t", index=False, na_rep="NA")

    # junction-like reads: carrier excess over the panel rate
    jpath = outdir / "junction_crossovers.tsv"
    if not jpath.is_file() or jpath.stat().st_size == 0:
        return
    j = pd.read_csv(jpath, sep="\t")
    if j.empty:
        return
    j["is_panel"] = j["mechanism"].isin(panel_mechanisms)
    rows = []
    coords = {pr.name: pr for pr in pairs}
    for pair, d in j.groupby("pair"):
        P = d[d["is_panel"]]
        denom = P["reads_in_blocks"].sum()
        rate = P["junction_reads"].sum() / denom if denom else np.nan
        prate = (P["junction_reads"] / P["reads_in_blocks"].replace(0, np.nan))
        for r in d[d["mechanism"].isin(deletion_mechanisms)].itertuples(index=False):
            if r.junction_reads == 0:
                continue
            mu = rate * r.reads_in_blocks if rate == rate else np.nan
            pr = coords.get(pair)
            at_edges = False
            if pr is not None and r.sample in edges and len(edges[r.sample]) == 2:
                p, d_ = edges[r.sample]["proximal"][0], edges[r.sample]["distal"][0]
                at_edges = (pr.a0 - 100_000 <= p <= pr.a1 + 100_000) and (pr.b0 - 100_000 <= d_ <= pr.b1 + 100_000)
            rows.append({"sample": r.sample, "mechanism": r.mechanism, "pair": pair,
                         "copy_A": r.copy_A, "copy_B": r.copy_B, "junction_reads": int(r.junction_reads),
                         "expected_from_panel": mu, "fold": r.junction_reads / mu if mu and mu > 0 else np.inf,
                         "poisson_p": _poisson_sf(int(r.junction_reads), mu) if mu == mu else np.nan,
                         "carrier_rate_per_1000": 1000 * r.junction_reads / r.reads_in_blocks,
                         "panel_max_rate_per_1000": 1000 * prate.max() if len(prate) else np.nan,
                         "panel_genomes_with_junction_reads": int((P["junction_reads"] > 0).sum()),
                         "pair_spans_carrier_edges": at_edges,
                         "crossover_in_A": r.crossover_in_A, "crossover_in_B": r.crossover_in_B})
    E = pd.DataFrame(rows)
    if not E.empty:
        E["excess"] = (E["poisson_p"] < 1e-3) & (E["carrier_rate_per_1000"] > E["panel_max_rate_per_1000"])
        E = E.sort_values(["sample", "poisson_p"])
    E.to_csv(outdir / "junction_enrichment.tsv", sep="\t", index=False, float_format="%.4g", na_rep="NA")
