#!/usr/bin/env python3
"""
Figure 5 -- Structural definition of 15q11-q13 deletion classes and breakpoint-proximal context.

Built from the new outputs of scripts/analysis and scripts/duplicons (plus HiFiCNV/pbsv VCFs):
  a  paralog-specific copy number (SUNK k-mers, 5-kb bins; scripts/duplicons 05) across
     chr15:17.5-33 Mb for every deletion carrier and PWS-mUPD, with the best available edges
     (figlib/deletions.py: assembled contig > contig junction > SUNK > HiFiCNV), the HiFiCNV
     edges (analysis 01) and the class and breakpoint status (duplicons 08)
  b  how far each edge moved from the HiFiCNV call, by the evidence that placed it
  c  pbsv PASS structural variants per genome, by group
  d  genome-wide HiFiCNV calls >= 2 Mb; the defining lesion of each genome (chr15 deletion over
     the imprinting centre, largest chr22 deletion of a 22q11.2DS genome) is outlined and left
     out of the burden test
  e  CpG methylation (pb-CpG-tools) in 10-100 kb bins on both sides of each edge, measured from
     the edge of unique sequence (an edge inside a segmental duplication is anchored at the
     duplication boundary; duplication CpGs excluded), minus the biparental mean
  f  near (<= 25 kb) minus far (50-100 kb) |difference| per group and side; bootstrap interval
     over carriers, exact sign-flip test

Usage: python3 scripts/figures/FIGURE_5.py [--results DIR] [--render-only]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figlib import annotation as ann  # noqa: E402
from figlib import cohort as cohort_lib  # noqa: E402
from figlib import deletions  # noqa: E402
from figlib import paths as paths_lib  # noqa: E402
from figlib import readers  # noqa: E402
from figlib import stats  # noqa: E402
from figlib import style  # noqa: E402
from figlib.style import mb  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

CHROM = ann.CHROM
PLOT = (17_500_000, 33_000_000)
CNV_MIN = 2_000_000
FLANK_BINS = ((0, 10_000), (10_000, 25_000), (25_000, 50_000), (50_000, 100_000))
NEAR_MAX, FAR_MIN = 25_000, 50_000
MIN_CPG_BIN = 5


# ------------------------------------------------------------------ a, b: copy number and edges
def cn_strips(paths, samples: list[str]) -> tuple[pd.DataFrame, str]:
    """5-kb paralog-specific copy number (duplicons 05), else HiFiCNV copy number in 5-kb bins."""
    path = paths.duplicons / "breakpoints" / "sunk15q.bins.tsv"
    if path.is_file():
        b = pd.read_csv(path, sep="\t", usecols=["sample", "start", "end", "cn_norm", "usable", "state"], low_memory=False)
        b = b[b["sample"].isin(samples)].rename(columns={"sample": "sample_id"})
        b["usable"] = b["usable"].astype(str).isin(["True", "true", "1"])
        return b, "SUNK copy number (scripts/duplicons 05)"
    rows = []
    bins = np.arange(PLOT[0], PLOT[1], 5_000)
    for s in samples:
        cn = readers.sample_files(paths.cnv, s, ("*.copynum.bedgraph", "*.copynum.bedgraph.gz"))
        if not cn:
            continue
        seg = readers.hificnv_copynum(cn[0], CHROM)
        val = np.full(len(bins), np.nan)
        for r in seg.itertuples(index=False):
            val[(bins + 2_500 >= r.start) & (bins + 2_500 < r.end)] = r.copy_number
        rows.append(pd.DataFrame({"sample_id": s, "start": bins, "end": bins + 5_000, "cn_norm": val,
                                  "usable": np.isfinite(val), "state": val}))
    return (pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()), "HiFiCNV copy number (no SUNK table)"


def edge_shifts(dels: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for r in dels.itertuples(index=False):
        for edge, pos, h, src in (("proximal", r.start, r.hificnv_start, r.start_source),
                                  ("distal", r.end, r.hificnv_end, r.end_source)):
            kind = ("assembled contig" if str(src).startswith("assembled") else
                    "contig junction" if str(src).startswith("contig junction") else
                    "SUNK copy number" if str(src).startswith("SUNK") else "HiFiCNV")
            rows.append({"sample_id": r.sample_id, "label": r.label, "group": r.group, "edge": edge,
                         "refined": pos, "hificnv": h, "shift_kb": (pos - h) / 1e3 if np.isfinite(pos) and np.isfinite(h) else np.nan,
                         "evidence": kind, "source": src})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ b, c: CNV and SV burden
def cnv_calls(paths, cohort) -> tuple[pd.DataFrame, list[str]]:
    rows, notes = [], []
    for s in cohort.samples:
        files = [f for f in readers.sample_files(paths.cnv, s, ("*.vcf.gz", "*.vcf")) if "hificnv" in f.name.lower()] \
            or readers.sample_files(paths.cnv, s, ("*.vcf.gz", "*.vcf"))
        if not files:
            notes.append(f"{s}: no HiFiCNV VCF")
            continue
        c = readers.hificnv_calls(files[0])
        c.insert(0, "sample_id", s)
        c.insert(1, "group", cohort.group(s))
        rows.append(c)
    calls = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=["sample_id", "group", "chrom", "start", "end", "size_bp", "svtype", "copy_number"])
    calls["defining_lesion"] = False
    for s, d in calls.groupby("sample_id"):
        g = cohort.group(s)
        if g in cohort_lib.DELETIONS:
            hit = d[(d["chrom"] == CHROM) & (d["svtype"] == "DEL") & (d["start"] < ann.IC[1]) & (d["end"] > ann.IC[0])]
        elif g == "DiGeorge":
            hit = d[(d["chrom"] == "chr22") & (d["svtype"] == "DEL") & (d["size_bp"] >= 1_500_000)]
            hit = hit.nlargest(1, "size_bp")
        else:
            hit = d.iloc[0:0]
        calls.loc[hit.index, "defining_lesion"] = True
    return calls, notes


def sv_burden(paths, cohort) -> tuple[pd.DataFrame, list[str]]:
    rows, notes = [], []
    for s in cohort.samples:
        f = readers.sample_files(paths.sv_calls, s, ("*.sv.pass.vcf.gz", "*.vcf.gz")) or \
            readers.sample_files(paths.phasing, s, ("*.sv.phased.vcf.gz",))
        if not f:
            notes.append(f"{s}: no pbsv VCF")
            continue
        rows.append({"sample_id": s, "group": cohort.group(s), "source": f[0].name, **readers.sv_counts(f[0])})
    return pd.DataFrame(rows), notes


def group_tests(df: pd.DataFrame, metrics: list[str], groups: list[str]) -> pd.DataFrame:
    rows = []
    for m in metrics:
        vals = [df.loc[df["group"] == g, m].to_numpy(float) for g in groups]
        h, p, eta = stats.kruskal(vals)
        rows.append({"metric": m, "kruskal_H": h, "eta2": eta, "permutation_p": stats.kruskal_permutation(vals, 10_000)})
    t = pd.DataFrame(rows)
    t["q"] = stats.bh(t["permutation_p"])
    return t


# ------------------------------------------------------------------ d, e: flanking methylation
def flank_table(paths, cohort, dels: pd.DataFrame, cis) -> pd.DataFrame:
    blocks, _ = ann.sd_blocks(paths)
    lo = int(np.nanmin(dels[["start", "end"]].to_numpy()) - 150_000)
    hi = int(np.nanmax(dels[["start", "end"]].to_numpy()) + 150_000)

    def load(s):
        p = cis.find_track(paths.methylation, s, "combined")
        if p is None:
            return None
        t = cis.read_track(p, CHROM, lo, hi)
        keep = np.ones(len(t.position), bool)
        for b in blocks.itertuples(index=False):
            keep &= ~((t.position >= b.start) & (t.position < b.end))
        return pd.Series(t.beta[keep], index=t.position[keep])

    ref = [load(s) for s in cohort.of(*cohort_lib.BIPARENTAL)]
    ref = [r for r in ref if r is not None]
    if not ref:
        return pd.DataFrame()
    control = pd.concat(ref, axis=1).mean(axis=1)
    tracks = {s: load(s) for s in list(dels["sample_id"]) + list(cohort.of("PWS-mUPD"))}
    def anchor(pos: float, direction: int) -> float:
        """The edge itself, or the boundary of the duplication containing it in `direction`."""
        for b in blocks.itertuples(index=False):
            if b.start <= pos < b.end:
                return float(b.end if direction > 0 else b.start)
        return float(pos)

    rows = []
    for r in dels.itertuples(index=False):
        for edge, pos, inward in (("proximal", r.start, 1), ("distal", r.end, -1)):
            if not np.isfinite(pos):
                continue
            for side, direction in (("inside (CN 1)", inward), ("outside (CN 2)", -inward)):
                a0 = anchor(pos, direction)
                for genome in [r.sample_id] + list(cohort.of("PWS-mUPD")):
                    tr = tracks.get(genome)
                    if tr is None:
                        continue
                    j = pd.concat([tr.rename("x"), control.rename("c")], axis=1, join="inner")
                    rel = (j.index.to_numpy() - a0) * direction
                    for b0, b1 in FLANK_BINS:
                        m = (rel >= b0) & (rel < b1)
                        n = int(m.sum())
                        rows.append({"carrier": r.sample_id, "genome": genome, "group": cohort.group(genome),
                                     "edge": edge, "edge_pos": pos, "side": side, "anchor": a0,
                                     "anchor_offset_bp": abs(a0 - pos), "bin": f"{b0 // 1000}-{b1 // 1000} kb",
                                     "distance_start": b0, "distance_end": b1, "n_cpg": n,
                                     "methylation": float(j["x"].to_numpy()[m].mean()) if n >= MIN_CPG_BIN else np.nan,
                                     "control": float(j["c"].to_numpy()[m].mean()) if n >= MIN_CPG_BIN else np.nan})
    t = pd.DataFrame(rows)
    t["delta"] = t["methylation"] - t["control"]
    return t


def flank_effects(flank: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if flank.empty:
        return pd.DataFrame()
    unit = flank.groupby(["carrier", "genome", "group", "edge", "side"])
    for (carrier, genome, group, edge, side), d in unit:
        near = d.loc[d["distance_end"] <= NEAR_MAX, "delta"].abs().mean()
        far = d.loc[d["distance_start"] >= FAR_MIN, "delta"].abs().mean()
        rows.append({"carrier": carrier, "genome": genome, "group": group, "edge": edge, "side": side,
                     "near_abs_delta": near, "far_abs_delta": far, "near_minus_far": near - far})
    u = pd.DataFrame(rows)
    out = []
    for (group, side), d in u.groupby(["group", "side"]):
        per_carrier = d.groupby("carrier")["near_minus_far"].mean().dropna().to_numpy()
        lo, hi = stats.bootstrap_ci(per_carrier) if group != "PWS-mUPD" else (np.nan, np.nan)
        out.append({"group": group, "side": side, "carriers": len(per_carrier),
                    "mean_near_minus_far": float(np.mean(per_carrier)) if len(per_carrier) else np.nan,
                    "ci_low": lo, "ci_high": hi,
                    "sign_flip_p": stats.sign_flip_p(per_carrier) if group != "PWS-mUPD" else np.nan})
    e = pd.DataFrame(out)
    e["q"] = stats.bh(e["sign_flip_p"])
    return e


# ------------------------------------------------------------------ figure
CN_CMAP = None


def cn_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(
        "cn", [(0.0, "#B2182B"), (0.25, "#EF8A62"), (0.5, "#F4F4F4"), (0.75, "#67A9CF"), (1.0, "#2166AC")])
    cmap.set_bad("#FFFFFF")
    return cmap


EVIDENCE_MARKER = {"assembled contig": "D", "contig junction": "*", "SUNK copy number": "o", "HiFiCNV": "s"}


def render(outdir: Path, t: dict, meta: dict) -> list[Path]:
    style.setup()
    dels, strips = t["deletions"], t["cn"]
    rows = list(dels.itertuples(index=False))
    extra = meta.get("extra_rows", [])
    n = len(rows) + len(extra)
    fig = plt.figure(figsize=(style.WIDTH, 10.8))
    gs = GridSpec(7, 1, figure=fig, height_ratios=[0.75 + 0.24 * n, 0.95, 1.35, 1.2, 1.25, 0.95, 1.55], hspace=0.0,
                  left=0.13, right=0.9, top=0.965, bottom=0.05)
    cmap = cn_cmap()

    # a -- copy-number strips
    sub = GridSpecFromSubplotSpec(n + 2, 1, subplot_spec=gs[0], height_ratios=[0.3, 0.55] + [1] * n, hspace=0.18)
    xl = (mb(PLOT[0]), mb(PLOT[1]))
    ax_sd = fig.add_subplot(sub[0])
    ann.draw_sd_track(ax_sd, t["sd_blocks"], xl)
    ax_sd.tick_params(labelbottom=False, bottom=False)
    ax_sd.spines["bottom"].set_visible(False)
    ax_g = fig.add_subplot(sub[1], sharex=ax_sd)
    ann.draw_gene_track(ax_g, t["genes"], xl, fontsize=style.BASE_FONT - 2)
    ax_g.tick_params(labelbottom=False, bottom=False)
    ax_g.spines["bottom"].set_visible(False)
    items = [(r.sample_id, r.label, r.group, r) for r in rows] + [(e["sample_id"], e["label"], e["group"], None) for e in extra]
    image = None
    for k, (sid, label, group, r) in enumerate(items):
        ax = fig.add_subplot(sub[k + 2], sharex=ax_sd)
        d = strips[strips["sample_id"] == sid].sort_values("start") if len(strips) else strips
        if len(d):
            v = np.where(d["usable"].to_numpy(bool), d["cn_norm"].to_numpy(float), np.nan)
            image = ax.imshow(v[None, :], aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=4,
                              extent=(mb(d["start"].min()), mb(d["end"].max()), 0, 1))
        if r is not None:
            for pos, h in ((r.start, r.hificnv_start), (r.end, r.hificnv_end)):
                if np.isfinite(h):
                    ax.plot([mb(h)], [1.12], marker="v", ms=3, mfc="white", mec=style.INK, mew=0.6, clip_on=False)
                if np.isfinite(pos):
                    ax.plot([mb(pos), mb(pos)], [-0.05, 1.05], color=style.INK, lw=0.9, clip_on=False)
            ax.text(1.005, 0.5, f"{r.deletion_class.split(' (')[0]} · {r.size_mb:.1f} Mb", transform=ax.transAxes,
                    va="center", fontsize=style.BASE_FONT - 2, color=style.INK)
        else:
            ax.text(1.005, 0.5, "copy-neutral", transform=ax.transAxes, va="center", fontsize=style.BASE_FONT - 2,
                    color=style.MUTED)
        ax.set_yticks([])
        ax.set_ylim(0, 1)
        ax.text(-0.005, 0.5, label, transform=ax.transAxes, ha="right", va="center",
                fontsize=style.BASE_FONT - 1.5, color=cohort_lib.COLOR[group], fontweight="bold")
        ax.tick_params(labelbottom=k == n - 1, left=False)
        for sp in ("left", "top", "right"):
            ax.spines[sp].set_visible(False)
    ax.set_xlim(*xl)
    ax.set_xlabel("chr15 (Mb, T2T-CHM13)")
    if image is not None:
        top = ax_sd.get_position().y1
        cax = fig.add_axes([0.76, top + 0.03, 0.14, 0.006])
        bar = fig.colorbar(image, cax=cax, orientation="horizontal", ticks=[0, 1, 2, 3, 4])
        bar.ax.tick_params(labelsize=style.BASE_FONT - 2, length=1.5, pad=1)
        cax.text(-0.04, 0.5, "copy number", transform=cax.transAxes, ha="right", va="center",
                 fontsize=style.BASE_FONT - 1.5)
    style.panel_label(fig, ax_sd, "a", f"Deletion architecture ({meta['cn_source'].split(' (')[0]}); "
                      "line, best edge; open triangle, HiFiCNV", dx=-0.12, dy=0.035)

    # b -- edge refinement | c -- SV burden
    sub = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[2], width_ratios=[1.35, 1, 1], wspace=0.55)
    ax = fig.add_subplot(sub[0])
    sh = t["edges"]
    labels_b = []
    for k, r in enumerate(dels.itertuples(index=False)):
        labels_b.append(r.label)
        for e in sh[sh["sample_id"] == r.sample_id].itertuples(index=False):
            if not np.isfinite(e.shift_kb):
                continue
            ax.plot(e.shift_kb, k + (-0.15 if e.edge == "proximal" else 0.15), marker=EVIDENCE_MARKER[e.evidence],
                    ls="", ms=5 if e.evidence != "contig junction" else 7, color=cohort_lib.COLOR[r.group],
                    mfc=cohort_lib.COLOR[r.group] if e.edge == "proximal" else "white", mew=0.9)
    ax.axvline(0, color=style.MUTED, lw=0.6)
    ax.set_yticks(range(len(labels_b)))
    ax.set_yticklabels(labels_b, fontsize=style.BASE_FONT - 1.5)
    ax.set_ylim(len(labels_b) - 0.5, -0.5)
    ax.set_xlabel("edge − HiFiCNV edge (kb)")
    short = {"assembled contig": "assembled", "contig junction": "junction", "SUNK copy number": "SUNK",
             "HiFiCNV": "HiFiCNV"}
    handles = [Line2D([], [], marker=m, ls="", color=style.MUTED, ms=4.5, label=short[k])
               for k, m in EVIDENCE_MARKER.items() if (sh["evidence"] == k).any()]
    handles += [Line2D([], [], marker="o", ls="", color=style.MUTED, ms=4.5, label="proximal"),
                Line2D([], [], marker="o", ls="", color=style.MUTED, mfc="white", ms=4.5, label="distal")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.4, -0.26), ncol=3, fontsize=style.BASE_FONT - 2,
              handletextpad=0.1, borderaxespad=0, columnspacing=0.6, frameon=False)
    style.panel_label(fig, ax, "b", "Edge refinement", dx=-0.1, dy=0.035)
    sv = t["sv"]
    rng = np.random.default_rng(7)
    for j, (metric, lab) in enumerate((("TOTAL_COUNT", "SVs per genome"), ("TOTAL_SPAN_MB", "SV span (Mb)"))):
        ax = fig.add_subplot(sub[j + 1])
        for i, g in enumerate(meta["groups"]):
            v = sv.loc[sv["group"] == g, metric].to_numpy(float) if len(sv) else []
            ax.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=14, marker=cohort_lib.MARKER[g],
                       color=cohort_lib.COLOR[g], edgecolor="white", lw=0.4)
            if len(v) > 1:
                ax.hlines(np.median(v), i - 0.3, i + 0.3, color=style.INK, lw=1)
        ax.set_xticks(range(len(meta["groups"])))
        ax.set_xticklabels([cohort_lib.PREFIX[g] for g in meta["groups"]], fontsize=style.BASE_FONT - 2)
        ax.set_ylabel(lab, fontsize=style.BASE_FONT - 1)
        r = t["sv_stats"][t["sv_stats"]["metric"] == metric] if len(t["sv_stats"]) else pd.DataFrame()
        if len(r):
            ax.set_title(f"p = {style.fmt_p(r['permutation_p'].iloc[0])}", fontsize=style.BASE_FONT - 1.5, color=style.MUTED)
        if j == 0:
            style.panel_label(fig, ax, "c", "pbsv structural variants", dx=-0.07, dy=0.035)

    # d -- genome-wide CNVs
    ax = fig.add_subplot(gs[4])
    calls, lens = t["cnv"], meta["chrom_lengths"]
    chroms = [c for c in lens if c not in ("chrX", "chrY", "chrM")]
    offs, x0 = {}, 0
    for c in chroms:
        offs[c] = x0
        x0 += lens[c]
    for i, c in enumerate(chroms):
        if i % 2 == 0:
            ax.axvspan(offs[c] / 1e6, (offs[c] + lens[c]) / 1e6, color="#F5F5F5", lw=0)
    big = calls[(calls["size_bp"] >= CNV_MIN) & calls["chrom"].isin(chroms)]
    lesion = big["defining_lesion"].astype(str).isin(["True", "true", "1"]) if len(big) else pd.Series(dtype=bool)
    for g in cohort_lib.GROUPS:
        sel = (big["group"] == g) if len(big) else pd.Series(dtype=bool)
        d, les = big[sel], lesion[sel]
        if d.empty:
            continue
        x = (d["chrom"].map(offs) + (d["start"] + d["end"]) / 2) / 1e6
        ax.scatter(x, d["size_bp"] / 1e6, s=np.where(les, 46, 16), marker=cohort_lib.MARKER[g], color=cohort_lib.COLOR[g],
                   edgecolor=np.where(les, style.INK, "white"), lw=0.5, alpha=0.9, zorder=3)
    ax.set_xticks([(offs[c] + lens[c] / 2) / 1e6 for c in chroms])
    ax.set_xticklabels([c.replace("chr", "") for c in chroms], fontsize=style.BASE_FONT - 2)
    ax.set_xlim(0, x0 / 1e6)
    ax.set_ylabel("CNV size (Mb)")
    ax.set_ylim(0, max(4, big["size_bp"].max() / 1e6 * 1.12 if len(big) else 4))
    st = t["cnv_stats"]
    txt = ", ".join(f"{'count' if 'count' in r.metric else 'total Mb'} p = {style.fmt_p(r.permutation_p)}"
                    for r in st.itertuples()) if len(st) else ""
    ax.text(1.0, 1.035, f"other CNVs ≥ 2 Mb by group: {txt}", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=style.BASE_FONT - 1.5, color=style.MUTED, clip_on=False)
    ax.legend(handles=style.group_legend_handles(meta["groups"]) +
              [Line2D([], [], marker="o", ls="", ms=6, mfc="white", mec=style.INK, label="defining lesion")],
              loc="upper center", ncol=6, fontsize=style.BASE_FONT - 2,
              bbox_to_anchor=(0.5, -0.15), columnspacing=0.75, handletextpad=0.25,
              borderaxespad=0)
    style.panel_label(fig, ax, "d", "Genome-wide HiFiCNV calls ≥ 2 Mb (autosomes)", dx=-0.12, dy=0.02)

    # e -- flank methylation | f -- effects
    sub = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[6], width_ratios=[1.25, 1], wspace=0.55)
    ax = fig.add_subplot(sub[0])
    flank = t["flank"]
    groups_d = [g for g in ("PWS-DEL", "AS-DEL", "PWS-mUPD") if len(flank) and (flank["group"] == g).any()]
    ends = []
    for g in groups_d:
        d = flank[flank["group"] == g]
        for side, sgn in (("inside (CN 1)", 1), ("outside (CN 2)", -1)):
            e = d[d["side"] == side].groupby(["distance_start", "distance_end", "carrier"])["delta"].mean().reset_index()
            m = e.groupby(["distance_start", "distance_end"])["delta"].agg(["mean", "count"]).reset_index()
            x = sgn * (m["distance_start"] + m["distance_end"]) / 2 / 1e3
            ax.plot(x, m["mean"], marker=cohort_lib.MARKER[g], color=cohort_lib.COLOR[g], lw=1, ms=3.5,
                    ls="--" if g == "PWS-mUPD" else "-")
            if sgn == 1 and len(m):
                ends.append((float(x.iloc[-1]), float(m["mean"].iloc[-1]), g))
    ax.axhline(0, color=style.LIGHT, lw=0.6)
    ax.axvline(0, color=style.INK, lw=0.6)
    ax.text(0.02, 0.97, "outside (CN 2)", transform=ax.transAxes, fontsize=style.BASE_FONT - 1.5, color=style.MUTED, va="top")
    ax.text(0.98, 0.97, "inside (CN 1)", transform=ax.transAxes, fontsize=style.BASE_FONT - 1.5, color=style.MUTED,
            va="top", ha="right")
    ax.set_xlabel("distance from the edge of unique sequence (kb)")
    ax.set_ylabel("Δ methylation vs\nbiparental genomes")
    if ends:                                   # direct labels at the right end of each line
        x_end = max(e[0] for e in ends)
        ys = sorted(ends, key=lambda e: e[1])
        lo_, hi_ = ax.get_ylim()
        gap = (hi_ - lo_) * 0.07
        placed = []
        for xe, ye, g in ys:
            y = ye if not placed or ye - placed[-1] >= gap else placed[-1] + gap
            placed.append(y)
            ax.text(x_end + 6, y, cohort_lib.DISPLAY[g], color=cohort_lib.COLOR[g], fontsize=style.BASE_FONT - 1.5,
                    va="center", ha="left")
        ax.set_xlim(right=x_end + 45)
    style.panel_label(fig, ax, "e", "Methylation flanking the breakpoints", dx=-0.07, dy=0.03)
    ax = fig.add_subplot(sub[1])
    eff = t["effects"]
    labels = []
    for k, r in enumerate(eff.itertuples(index=False)):
        c = cohort_lib.COLOR[r.group]
        if np.isfinite(r.ci_low):
            ax.plot([r.ci_low, r.ci_high], [k, k], color=c, lw=1.4)
        ax.plot(r.mean_near_minus_far, k, marker=cohort_lib.MARKER[r.group], color=c, ms=5,
                mfc="white" if r.group == "PWS-mUPD" else c)
        labels.append(f"{cohort_lib.DISPLAY[r.group]} {r.side.split(' (')[0]} (n={r.carriers}"
                      + (f"; q={style.fmt_p(r.q)})" if np.isfinite(r.q) else ")"))
    ax.axvline(0, color=style.MUTED, lw=0.6)
    ax.set_yticks(range(len(eff)))
    ax.set_yticklabels(labels, fontsize=style.BASE_FONT - 2)
    ax.set_ylim(len(eff) - 0.5, -0.5)
    ax.set_xlabel("near − far |Δ|")
    style.panel_label(fig, ax, "f", "Near vs far", dx=-0.07, dy=0.03)
    return style.save(fig, outdir, "Figure5")


def report(outdir, t, meta, notes, paths_out) -> Path:
    r = style.Report("Figure 5 report: deletion structure and breakpoint context")
    r.p("Generated by `scripts/figures/FIGURE_5.py` from scripts/analysis 01, scripts/duplicons 03/05/08, "
        "HiFiCNV, pbsv and pb-CpG-tools outputs.")
    r.h("Deletions (panel a)")
    r.p(f"Copy number drawn from: {meta['cn_source']}.")
    r.table(t["deletions"][["sample_id", "label", "group", "start", "end", "size_mb", "start_source", "end_source",
                            "deletion_class", "status", "analysis01_type"]])
    if len(t["junctions"]):
        r.p("Contig junctions used to refine an edge (hifiasm alignments of scripts/duplicons 03):")
        r.table(t["junctions"])
    r.h("Edge refinement (panel b)")
    r.table(t["edges"])
    r.h("Genome-wide CNVs >= 2 Mb (panel d)")
    r.table(t["cnv_burden"])
    r.table(t["cnv_stats"])
    r.h("pbsv structural variants (panel c)")
    r.table(t["sv"])
    r.table(t["sv_stats"])
    r.h("Breakpoint-flanking methylation (panels e, f)")
    fl = t["flank"]
    if len(fl):
        anchors = fl.drop_duplicates(["carrier", "edge", "side"])[["carrier", "edge", "edge_pos", "side", "anchor",
                                                                    "anchor_offset_bp"]]
        r.p("Anchors (an edge inside a segmental duplication is measured from the duplication boundary):")
        r.table(anchors)
    r.table(t["effects"])
    if notes:
        r.h("Notes")
        for n_ in notes:
            r.p(f"- {n_}")
    cls = t["deletions"]["deletion_class"].value_counts().to_dict()
    r.h("Caption draft")
    r.p(f"**Structural definition of 15q11–q13 deletion classes and breakpoint-proximal context.** **a,** Paralog-specific "
        f"copy number (singly-unique k-mers, 5-kb bins) across chr15:17.5–33 Mb in the {len(t['deletions'])} deletion "
        f"genomes and PWS-mUPD; vertical lines, best-resolved deletion edges; open triangles, HiFiCNV edges; classes "
        f"{'; '.join(f'{v} {k}' for k, v in cls.items())}. Top, segmental duplications (BP1–BP3) and genes. "
        f"**b,** Displacement of each edge from the HiFiCNV call, by the evidence that placed it (assembled contig, "
        f"contig junction, SUNK copy number). **c,** pbsv structural variants per genome. **d,** Genome-wide HiFiCNV "
        f"calls ≥ 2 Mb; outlined, the defining lesion of each genome, excluded from the burden test (Kruskal–Wallis "
        f"permutation). **e,** CpG methylation in 10–100-kb bins on both sides of each edge, measured from the edge of unique "
        f"sequence, minus the biparental mean; dashed, PWS-mUPD at the carriers' coordinates. **f,** Near-minus-far |Δ| "
        f"with carrier bootstrap 95% CI; q, sign-flip test with Benjamini–Hochberg correction.")
    r.p("")
    r.p("Files: " + ", ".join(p_.name for p_ in paths_out))
    return r.write(outdir, "Figure5_report.md")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    paths_lib.add_common_arguments(ap)
    ap.add_argument("--render-only", action="store_true", help="redraw from the tables of a previous run")
    a = ap.parse_args(argv)
    paths = paths_lib.resolve(a)
    outdir = paths.figure_dir(5, a.outdir)
    names = {"deletions": "Figure5a_deletions.tsv", "junctions": "Figure5a_contig_junctions.tsv",
             "cn": "Figure5a_copy_number_5kb.tsv.gz", "sd_blocks": "Figure5a_sd_blocks.tsv", "genes": "Figure5a_genes.tsv",
             "edges": "Figure5b_edge_refinement.tsv", "cnv": "Figure5d_cnv_calls.tsv",
             "cnv_burden": "Figure5d_cnv_burden.tsv", "cnv_stats": "Figure5d_statistics.tsv",
             "sv": "Figure5c_sv_burden.tsv", "sv_stats": "Figure5c_statistics.tsv",
             "flank": "Figure5e_flanking_methylation.tsv", "effects": "Figure5f_effects.tsv"}
    if a.render_only:
        t = {k: style.read_table(outdir, v) for k, v in names.items()}
        meta = json.loads((outdir / "tables" / "Figure5_meta.json").read_text())
        paths_out = render(outdir, t, meta)
        print(report(outdir, t, meta, meta.get("notes", []), paths_out))
        return
    cis = paths_lib.import_analysis_package()
    cohort = cohort_lib.load(paths.metadata)
    groups = cohort.groups_present()
    dels, junctions = deletions.deletion_table(paths, cohort)
    extra = [{"sample_id": s, "label": cohort.label(s), "group": cohort.group(s)} for s in cohort.of("PWS-mUPD")]
    cn, cn_source = cn_strips(paths, list(dels["sample_id"]) + [e["sample_id"] for e in extra])
    edges = edge_shifts(dels)
    cnv, n2 = cnv_calls(paths, cohort)
    other = cnv[(cnv["size_bp"] >= CNV_MIN) & ~cnv["defining_lesion"] & ~cnv["chrom"].isin(["chrX", "chrY"])]
    burden = pd.DataFrame({"sample_id": list(cohort.samples), "group": [cohort.group(s) for s in cohort.samples]})
    agg = other.groupby("sample_id").agg(other_cnv_ge2mb_count=("size_bp", "size"),
                                         other_cnv_ge2mb_total_mb=("size_bp", lambda x: x.sum() / 1e6))
    burden = burden.merge(agg, on="sample_id", how="left").fillna({"other_cnv_ge2mb_count": 0, "other_cnv_ge2mb_total_mb": 0})
    cnv_stats = group_tests(burden, ["other_cnv_ge2mb_count", "other_cnv_ge2mb_total_mb"], groups)
    sv, n3 = sv_burden(paths, cohort)
    sv_stats = group_tests(sv, ["TOTAL_COUNT", "TOTAL_SPAN_MB", "DEL", "INS", "DUP", "INV", "BND"], groups) \
        if len(sv) else pd.DataFrame()
    flank = flank_table(paths, cohort, dels, cis)
    effects = flank_effects(flank)
    blocks, _ = ann.sd_blocks(paths)
    genes = ann.landmark_genes(paths, *PLOT)
    lens = readers.fai_lengths(paths.fai()) or readers.CHM13_LENGTHS
    meta = {"groups": groups, "chrom_lengths": {c: lens[c] for c in readers.CHM13_LENGTHS if c in lens},
            "notes": n2 + n3, "cn_source": cn_source, "extra_rows": extra}
    t = {"deletions": dels, "junctions": junctions, "cn": cn, "sd_blocks": blocks, "genes": genes, "edges": edges,
         "cnv": cnv, "cnv_burden": burden, "cnv_stats": cnv_stats, "sv": sv, "sv_stats": sv_stats, "flank": flank,
         "effects": effects}
    for k, v in names.items():
        style.write_table(t[k], outdir, v)
    (outdir / "tables" / "Figure5_meta.json").write_text(json.dumps(meta, indent=1))
    paths_out = render(outdir, t, meta)
    print(report(outdir, t, meta, meta["notes"], paths_out))


if __name__ == "__main__":
    main()
