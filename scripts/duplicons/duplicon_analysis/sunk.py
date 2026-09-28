"""
sunk.py -- paralog-specific copy number across the 15q11-q13 BP blocks from
HiFi read k-mer counts at singly-unique k-mers (SUNKs), to localize deletion
breakpoints INSIDE the segmental duplications where read-depth callers
(HiFiCNV) go blind.

Principle (Sudmant et al. 2010; Nuttle et al. 2014): a SUNK is a k-mer that
occurs once in the reference genome, so it tags one paralog of a duplicon. Its
count in a sample's reads, divided by the diploid k-mer coverage, is the copy
number of THAT paralog. In a hemizygous NAHR deletion the SUNKs of the deleted
part of each duplicon drop to CN 1 and the switch point is the crossover. For a
direct-orientation NAHR event the two switches (one in the proximal BP, one in
BP3) must fall at homologous positions of the paralog pair -- the sharpest
internal test that the event is NAHR between those two copies.

Inputs (made by 01_prepare_reference.py and 02_count_sunks.py)
  --window-fasta   CHM13 window (samtools faidx, header chr15:START-END)
  --control-fasta  a diploid, SD-free CHM13 region used for normalization
  --sunks          meryl print of the SUNK set (window+control k-mers with
                   genome-wide CHM13 count == 1): "KMER<TAB>COUNT"
  --samples        TSV: sample, group, counts[, kcov]; counts = meryl print of
                   (sample read k-mers  intersect  SUNKs)
  --sd-bed, --landmarks-bed, --sd-pairs (BEDPE; optional)

Outputs
  <prefix>.bins.tsv          per sample x bin: SUNK n, raw CN, panel-normalized CN, HMM state
  <prefix>.transitions.tsv   state changes, bin- and SUNK-level position + support
                             interval, SD/landmark annotation, homologous positions
  <prefix>.nahr_pairs.tsv    for each deletion sample: distance between one
                             transition's homologous position and the other transition
  <prefix>.sunk_qc.tsv       SUNK counts, fraction behaving diploid in the panel
  <prefix>.overview.png      CN along the window for every sample
  <prefix>.zoom_<sample>.png per-SUNK view around each transition
"""
from __future__ import annotations

import argparse
import math
import os
import re
import sys

from pathlib import Path

import numpy as np
import pandas as pd

from .kmers import (  # noqa: E402
    header_offset,
    kmer_strings_to_codes,
    lookup,
    read_counts,
    read_fasta,
    seq_codes,
    viterbi,
)

# plotting colours (validated reference categorical slots 1-3; CN0/CN>=4 folded to grey)
C_CN2, C_CN1, C_CN3, C_OTHER = "#2a78d6", "#eb6834", "#1baf7a", "#8a8984"
INK, INK2, GRID, SDFILL = "#0b0b0b", "#52514e", "#e4e3df", "#ecebe7"


# ------------------------------------------------------------------ intervals
def read_bed(path, chrom):
    out = []
    if not path:
        return out
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == chrom:
                out.append((int(f[1]), int(f[2]), f[3] if len(f) > 3 and f[3] else f"{f[1]}-{f[2]}"))
    return sorted(out)


def read_pairs(path, chrom):
    """BEDPE-like: chr1 s1 e1 chr2 s2 e2 [name score strand1 strand2 ...].
    Returns both directions so any side can be projected."""
    out = []
    if not path:
        return out
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            c1, s1, e1, c2, s2, e2 = f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5])
            st1 = f[8] if len(f) > 9 else "+"
            st2 = f[9] if len(f) > 9 else "+"
            same = "+" if st1 == st2 else "-"
            name = f[6] if len(f) > 6 else "."
            if c1 == chrom and c2 == chrom:
                out.append((s1, e1, s2, e2, same, name))
                out.append((s2, e2, s1, e1, same, name))
    return out


def hits(iv, pos):
    return [n for s, e, n in iv if s <= pos < e]


def homologs(pairs, pos):
    res = []
    for s1, e1, s2, e2, same, name in pairs:
        if s1 <= pos < e1:
            f = (pos - s1) / max(1, e1 - s1)
            h = s2 + f * (e2 - s2) if same == "+" else e2 - f * (e2 - s2)
            res.append((int(round(h)), same, name))
    return res


# ----------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--window-fasta", required=True)
    ap.add_argument("--control-fasta", required=True)
    ap.add_argument("--sunks", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--panel-groups", default="Control,DiGeorge,PWS-mUPD",
                    help="groups with diploid chr15 used as the per-SUNK reference")
    ap.add_argument("--sd-bed")
    ap.add_argument("--landmarks-bed")
    ap.add_argument("--sd-pairs", help="BEDPE of SD pairs (BISER/SEDEF) for the NAHR check")
    ap.add_argument("--bin-bp", type=int, default=5000)
    ap.add_argument("--min-sunks", type=int, default=8)
    ap.add_argument("--p-switch", type=float, default=1e-4)
    ap.add_argument("--zoom-bp", type=int, default=20_000)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--prefix", default="sunk15q")
    a = ap.parse_args(argv)
    k = a.k
    os.makedirs(a.outdir, exist_ok=True)
    pre = os.path.join(a.outdir, a.prefix)
    log = lambda m: print(f"[05 copy number] {m}", file=sys.stderr, flush=True)  # noqa: E731

    # --- reference sequences
    win = read_fasta(a.window_fasta)
    if len(win) != 1:
        sys.exit("--window-fasta must contain exactly one record")
    chrom, W0 = header_offset(win[0][0])
    wseq = win[0][1]
    W1 = W0 + len(wseq)
    log(f"window {chrom}:{W0}-{W1} ({len(wseq) / 1e6:.2f} Mb)")
    wcodes, wvalid = seq_codes(wseq, k)
    ccodes_all = []
    for _, s in read_fasta(a.control_fasta):
        c, v = seq_codes(s, k)
        ccodes_all.append(c[v])
    ccodes = np.concatenate(ccodes_all) if ccodes_all else np.zeros(0, np.uint64)

    # --- SUNK set, restricted to k-mers seen exactly once in window+control
    sunk_codes = np.unique(kmer_strings_to_codes(
        pd.read_csv(a.sunks, sep=r"\s+", header=None, usecols=[0], names=["kmer"], dtype=str)["kmer"], k))
    allpos = np.concatenate([wcodes[wvalid], ccodes])
    u, cnt = np.unique(allpos[np.isin(allpos, sunk_codes)], return_counts=True)
    sunk_codes = u[cnt == 1]
    w_is = wvalid & np.isin(wcodes, sunk_codes)
    w_pos = np.nonzero(w_is)[0]
    w_q = wcodes[w_pos]
    gpos = W0 + w_pos + k // 2                                     # genomic midpoint
    c_q = ccodes[np.isin(ccodes, sunk_codes)]
    log(f"{w_q.size:,} SUNKs in window, {c_q.size:,} in control region(s)")
    if c_q.size < 1000:
        log("WARNING: <1000 control SUNKs; normalization will be noisy")

    sds = read_bed(a.sd_bed, chrom)
    lms = read_bed(a.landmarks_bed, chrom)
    pairs = read_pairs(a.sd_pairs, chrom)
    in_sd = np.zeros(gpos.size, bool)
    for s, e, _ in sds:
        in_sd |= (gpos >= s) & (gpos < e)

    # --- samples
    sheet = pd.read_csv(a.samples, sep="\t", dtype=str).fillna("")
    for col in ("sample", "group", "counts"):
        if col not in sheet.columns:
            sys.exit(f"--samples needs a '{col}' column")
    panel_groups = {g.strip() for g in a.panel_groups.split(",")}
    cn = {}
    kcov = {}
    for _, r in sheet.iterrows():
        codes, vals = read_counts(r["counts"], k)
        wc = lookup(codes, vals, w_q).astype(np.float32)
        if r.get("kcov", ""):
            kc = float(r["kcov"])
        else:
            cc = lookup(codes, vals, c_q)
            kc = float(np.median(cc))           # diploid k-mer coverage
        if kc <= 0:
            sys.exit(f"{r['sample']}: control k-mer coverage is 0 -- wrong counts file?")
        kcov[r["sample"]] = kc
        cn[r["sample"]] = 2.0 * wc / kc
        log(f"{r['sample']} ({r['group']}): diploid k-mer coverage {kc:.1f}")

    groups = dict(zip(sheet["sample"], sheet["group"]))
    panel = [s for s in cn if groups[s] in panel_groups]
    log(f"panel ({len(panel)}): {', '.join(panel) or 'NONE'}")

    def panel_median(exclude=None):
        mem = [s for s in panel if s != exclude]
        if len(mem) < 3:
            return None
        return np.median(np.vstack([cn[s] for s in mem]), axis=0)

    pmed_all = panel_median()
    if pmed_all is None:
        log("WARNING: <3 panel samples; using raw CN without per-SUNK normalization")
        ok_all = np.ones(w_q.size, bool)
    else:
        ok_all = (pmed_all >= 1.5) & (pmed_all <= 2.5)
    pd.DataFrame({
        "set": ["window", "window_in_SD", "window_outside_SD"],
        "n_sunks": [w_q.size, int(in_sd.sum()), int((~in_sd).sum())],
        "fraction_diploid_in_panel": [ok_all.mean() if ok_all.size else np.nan,
                                      ok_all[in_sd].mean() if in_sd.any() else np.nan,
                                      ok_all[~in_sd].mean() if (~in_sd).any() else np.nan],
    }).to_csv(f"{pre}.sunk_qc.tsv", sep="\t", index=False)

    # --- per-SUNK normalized CN, bins
    bin_id = (gpos - W0) // a.bin_bp
    rel = {}
    for s in cn:
        pm = panel_median(exclude=s) if groups[s] in panel_groups else pmed_all
        if pm is None:
            r_ = cn[s].copy()
            ok = np.ones_like(r_, bool)
        else:
            ok = (pm >= 1.5) & (pm <= 2.5)
            r_ = np.where(ok, 2.0 * cn[s] / np.where(pm > 0, pm, np.nan), np.nan)
        rel[s] = r_

    rows = []
    for s in cn:
        d = pd.DataFrame({"bin": bin_id, "raw": cn[s], "rel": np.clip(rel[s], 0, 6)})
        g = d.dropna(subset=["rel"]).groupby("bin")
        # clipped mean rather than median: k-mer counts are integers and bin medians
        # of integers tie, which would understate the noise the HMM sees
        b = pd.DataFrame({"n_sunks": g.size(), "cn_raw": g["raw"].median(),
                          "cn_norm": g["rel"].mean()}).reset_index()
        b.insert(0, "sample", s)
        b.insert(1, "group", groups[s])
        rows.append(b)
    bins = pd.concat(rows, ignore_index=True)
    bins["start"] = W0 + bins["bin"] * a.bin_bp
    bins["end"] = bins["start"] + a.bin_bp
    bins["usable"] = bins["n_sunks"] >= a.min_sunks

    # per-SUNK noise from the panel's own per-SUNK values: sd(bin) = sd1 / sqrt(n)
    pv = np.concatenate([np.clip(rel[s], 0, 6) for s in cn if groups[s] in panel_groups] or [np.array([])])
    pv = pv[np.isfinite(pv)]
    sd1 = 1.4826 * float(np.median(np.abs(pv - np.median(pv)))) if pv.size > 1000 else 0.5
    sd1 = max(sd1, 0.15)
    # Read depth also varies along the genome in a correlated way (read placement,
    # GC), which averaging many SUNKs in a bin does not remove. Estimate that
    # bin-level component from the panel: total bin variance minus the part
    # explained by independent per-SUNK noise.
    pb = bins[bins["group"].isin(panel_groups) & bins["usable"]]
    if len(pb) > 20:
        z = pb["cn_norm"].to_numpy() - 2.0
        total = 1.4826 * float(np.median(np.abs(z - np.median(z))))
        indep = float(np.median(sd1 ** 2 / pb["n_sunks"].to_numpy()))
        sd_bin = float(np.sqrt(max(0.0, total ** 2 - indep)))
    else:
        sd_bin = 0.15
    sd_bin = max(sd_bin, 0.05)
    log(f"per-SUNK CN noise sd ~ {sd1:.2f}; correlated bin-level noise sd ~ {sd_bin:.2f}")

    states = [0, 1, 2, 3, 4]
    bins["state"] = np.nan
    for s in cn:
        m = (bins["sample"] == s) & bins["usable"]
        sub = bins[m].sort_values("bin")
        if sub.empty:
            continue
        sd = np.sqrt(sd_bin ** 2 + sd1 ** 2 / sub["n_sunks"].to_numpy())
        bins.loc[sub.index, "state"] = viterbi(sub["cn_norm"].to_numpy(), sd, states, a.p_switch)
    bins["in_sd"] = [bool(hits(sds, int((s + e) // 2))) for s, e in zip(bins["start"], bins["end"])]
    bins["landmark"] = [";".join(hits(lms, int((s + e) // 2))) or "." for s, e in zip(bins["start"], bins["end"])]
    bins.to_csv(f"{pre}.bins.tsv", sep="\t", index=False, float_format="%.4f")

    # --- transitions with SUNK-level refinement
    # The HMM places a switch between two bins, but with correlated depth noise the
    # position is only approximate: the SUNK-level search covers the two bins at
    # the switch, and the support interval comes from a 1-kb block bootstrap of the
    # SUNK values (keeps local correlation). Expect several kb of error at 20-30x;
    # junction reads (05) give the crossover itself.
    trows = []
    rng_boot = np.random.default_rng(2026)

    def best_split(x, s1, s2):
        c1 = np.concatenate(([0.0], np.cumsum((x - s1) ** 2)))
        c2 = np.concatenate((np.cumsum(((x - s2) ** 2)[::-1])[::-1], [0.0]))
        return int(np.argmin(c1 + c2))                      # split before index j

    for s in cn:
        sub = bins[(bins["sample"] == s) & bins["usable"]].sort_values("bin").reset_index(drop=True)
        st = sub["state"].to_numpy()
        change = [i for i in range(1, len(sub)) if st[i] != st[i - 1]]
        for k, i in enumerate(change):
            s1, s2 = float(st[i - 1]), float(st[i])
            lo_i = max(i - 1, change[k - 1] if k > 0 else 0)
            hi_i = min(i, (change[k + 1] - 1) if k + 1 < len(change) else len(sub) - 1)
            lo, hi = int(sub.loc[lo_i, "start"]), int(sub.loc[hi_i, "end"])
            m = (gpos >= lo) & (gpos < hi) & np.isfinite(rel[s])
            x, p = np.clip(rel[s][m], 0, 6), gpos[m]
            if x.size >= 2:
                j = best_split(x, s1, s2)
                refined = int(p[min(j, x.size - 1)])
                block = (p - lo) // 1000
                ids = np.unique(block)
                boots = []
                for _ in range(200):
                    # resample residuals by 1-kb block around the fitted step
                    fit = np.where(np.arange(x.size) < j, s1, s2)
                    res = x - fit
                    pick = rng_boot.choice(ids, ids.size)
                    shuffled = np.concatenate([res[block == b] for b in pick])[: x.size]
                    if shuffled.size < x.size:
                        shuffled = np.resize(shuffled, x.size)
                    jb = best_split(fit + shuffled, s1, s2)
                    boots.append(p[min(jb, x.size - 1)])
                ci_lo, ci_hi = (int(v) for v in np.percentile(boots, [2.5, 97.5]))
                ci_lo, ci_hi = min(ci_lo, refined), max(ci_hi, refined)
            else:
                refined, ci_lo, ci_hi = int(sub.loc[i - 1, "end"]), lo, hi
            hom = homologs(pairs, refined)
            trows.append({
                "sample": s, "group": groups[s], "from_cn": int(s1), "to_cn": int(s2),
                "bin_boundary": int(sub.loc[i - 1, "end"]), "refined_pos": refined,
                "support_start": ci_lo, "support_end": ci_hi,
                "n_sunks_used": int(x.size),
                "sd_at_refined": ";".join(hits(sds, refined)) or ".",
                "landmark": ";".join(hits(lms, refined)) or ".",
                "homologous_positions": ";".join(f"{h}({o},{n})" for h, o, n in hom) or ".",
            })
    tr = pd.DataFrame(trows)
    tr.to_csv(f"{pre}.transitions.tsv", sep="\t", index=False)

    # --- NAHR consistency: does one switch project onto the other?
    nrows = []
    if not tr.empty and pairs:
        for s, g in tr.groupby("sample"):
            g = g.reset_index(drop=True)
            for i in range(len(g)):
                for j in range(len(g)):
                    if i == j:
                        continue
                    a_, b_ = g.loc[i], g.loc[j]
                    if not (a_["from_cn"] == 2 and a_["to_cn"] == 1 and b_["from_cn"] == 1 and b_["to_cn"] == 2):
                        continue
                    hs = homologs(pairs, int(a_["refined_pos"]))
                    if not hs:
                        continue
                    h, o, n = min(hs, key=lambda t: abs(t[0] - b_["refined_pos"]))
                    width = (a_["support_end"] - a_["support_start"]) + (b_["support_end"] - b_["support_start"])
                    nrows.append({"sample": s, "enter_deletion": int(a_["refined_pos"]),
                                  "exit_deletion": int(b_["refined_pos"]),
                                  "homolog_of_enter": h, "pair_orientation": o, "pair": n,
                                  "distance_bp": int(abs(h - b_["refined_pos"])),
                                  "combined_support_bp": int(width),
                                  "consistent_with_NAHR": bool(abs(h - b_["refined_pos"]) <= max(width, 2 * a.bin_bp))})
    pd.DataFrame(nrows).to_csv(f"{pre}.nahr_pairs.tsv", sep="\t", index=False)

    # --- plots
    plot_overview(bins, tr, sds, lms, W0, W1, pre, panel_groups)
    for s in tr["sample"].unique() if not tr.empty else []:
        plot_zoom(s, tr[tr["sample"] == s], gpos, rel[s], sds, a.zoom_bp, pre)
    log(f"done: {len(tr)} transitions; outputs {pre}.*")


# ------------------------------------------------------------------- plotting
def _state_color(v):
    return {1: C_CN1, 2: C_CN2, 3: C_CN3}.get(int(v), C_OTHER) if np.isfinite(v) else C_OTHER


def plot_overview(bins, tr, sds, lms, W0, W1, pre, panel_groups):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    samples = list(dict.fromkeys(
        list(bins.loc[~bins["group"].isin(panel_groups), "sample"]) +
        list(bins.loc[bins["group"].isin(panel_groups), "sample"])))
    n = len(samples)
    fig, axes = plt.subplots(n, 1, figsize=(11, 0.9 * n + 1.2), sharex=True, squeeze=False)
    for ax, s in zip(axes[:, 0], samples):
        b = bins[(bins["sample"] == s) & bins["usable"]]
        for x0, x1, _ in sds:
            if x1 > W0 and x0 < W1:
                ax.axvspan(x0 / 1e6, x1 / 1e6, color=SDFILL, lw=0, zorder=0)
        for y in (1, 2, 3):
            ax.axhline(y, color=GRID, lw=0.6, zorder=1)
        ax.scatter((b["start"] + b["end"]) / 2e6, b["cn_norm"].clip(0, 4.5), s=4, lw=0,
                   c=[_state_color(v) for v in b["state"]], zorder=2)
        if not tr.empty:
            for p in tr.loc[tr["sample"] == s, "refined_pos"]:
                ax.axvline(p / 1e6, color=INK2, lw=0.6, ls=(0, (2, 2)), zorder=3)
        ax.set_ylim(-0.2, 4.6)
        ax.set_yticks([0, 1, 2, 3, 4])
        ax.tick_params(labelsize=6, colors=INK2, length=2)
        ax.set_ylabel(f"{s}\n{b['group'].iloc[0] if len(b) else ''}", fontsize=6.5,
                      color=INK, rotation=0, ha="right", va="center")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
    top = axes[0, 0]
    for x0, x1, nme in lms:
        if x1 > W0 and x0 < W1:
            top.text((x0 + x1) / 2e6, 4.9, nme, ha="center", va="bottom", fontsize=7, color=INK)
    axes[-1, 0].set_xlabel(f"CHM13 position (Mb)", fontsize=8, color=INK)
    handles = [Line2D([], [], marker="o", ls="", color=c, label=l, markersize=4)
               for c, l in ((C_CN2, "CN 2"), (C_CN1, "CN 1"), (C_CN3, "CN 3"), (C_OTHER, "CN 0 / ≥4"))]
    handles.append(Line2D([], [], color=SDFILL, lw=6, label="segmental duplication"))
    fig.legend(handles=handles, loc="upper right", ncol=5, fontsize=7, frameon=False)
    fig.suptitle("Paralog-specific copy number at SUNKs (panel-normalized, 5-state HMM)",
                 x=0.01, ha="left", fontsize=9, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(f"{pre}.overview.png", dpi=200)
    plt.close(fig)


def plot_zoom(sample, trs, gpos, rel, sds, zoom, pre):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(trs)
    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 2.8), squeeze=False)
    for ax, (_, t) in zip(axes[0], trs.iterrows()):
        c = t["refined_pos"]
        m = (gpos >= c - zoom) & (gpos < c + zoom) & np.isfinite(rel)
        for x0, x1, _ in sds:
            if x1 > c - zoom and x0 < c + zoom:
                ax.axvspan(max(x0, c - zoom) / 1e3, min(x1, c + zoom) / 1e3, color=SDFILL, lw=0)
        ax.axvspan(t["support_start"] / 1e3, t["support_end"] / 1e3, color=C_CN1, alpha=0.15, lw=0)
        xs, ys = gpos[m], np.clip(rel[m], 0, 6)
        # per-SUNK values are Poisson-noisy; show 500-bp bin means as the readable signal
        if xs.size:
            bb = (xs - (c - zoom)) // 500
            dfz = pd.DataFrame({"b": bb, "y": ys}).groupby("b")["y"].agg(["mean", "size"])
            dfz = dfz[dfz["size"] >= 3]
            bx = (c - zoom + dfz.index.to_numpy() * 500 + 250) / 1e3
            ax.plot(bx, dfz["mean"].to_numpy(), color=C_CN2, lw=1.2, marker="o", ms=2.5,
                    label="500-bp mean")
        ax.axvline(c / 1e3, color=INK, lw=0.8)
        for y in (1, 2):
            ax.axhline(y, color=GRID, lw=0.6)
        ax.set_ylim(-0.2, 3.5)
        ax.set_title(f"CN {t['from_cn']}→{t['to_cn']} at {c:,}\n{t['sd_at_refined']}", fontsize=7, color=INK)
        ax.set_xlabel("CHM13 position (kb)", fontsize=7, color=INK)
        ax.tick_params(labelsize=6, colors=INK2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0, 0].set_ylabel("paralog-specific CN", fontsize=7, color=INK)
    fig.suptitle(sample, x=0.01, ha="left", fontsize=9, color=INK)
    fig.tight_layout()
    fig.savefig(f"{pre}.zoom_{sample}.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
