#!/usr/bin/env python3
"""
Figure 3 -- Boundary mapping of the parent-of-origin methylation transition across 15q11-q13.

Computed from the pb-CpG-tools tracks (results/06_methylation):
  per CpG   controls   |maternal - paternal| (hap1/hap2 oriented by Figure 1's assignment,
                       else by imprinting-centre methylation: the more methylated haplotype
                       is maternal)
            PWS-DEL    combined - control baseline      (maternal copy retained)
            AS-DEL     control baseline - combined      (paternal copy retained)
            baseline = mean of the control combined tracks
  windows   1 kb every 100 bp across chr15:22-29 Mb (mean of per-CpG signals)
  boundary  a run of >= 5 windows above 0.4 enters the transition; the first later run of
            >= 5 windows below 0.1 exits it. Called per genome and on each group's mean profile.
  core      the control consensus segment that agrees best with the PWS-DEL and AS-DEL calls;
            shared core = [latest entry, earliest exit] of the three groups.
  sensitivity  250/500-bp and 2-kb windows, CpG-count windows (25, 50 CpGs), the 0.4/0.15
            threshold, change-point on the smoothed profiles, and a CpG bootstrap (100 replicates).

Panels  a locus-wide contrast (22-29 Mb) with duplication blocks and genes
        b boundary intervals per group and the shared core
        c CpG-level signals and smoothed profiles over the core
        d genomic context: distance to BP1-BP3 duplications; genes and ICR around the core
        e boundary calls under every sensitivity setting
Every coordinate and number in the figure and report is computed; none is hard-coded.

Usage: python3 scripts/figures/FIGURE_3.py [--results DIR] [--include-22q] [--signed]
       [--render-only]   (redraw from results/07_figures/figure_3/tables)
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figlib import annotation as ann  # noqa: E402
from figlib import cohort as cohort_lib  # noqa: E402
from figlib import paths as paths_lib  # noqa: E402
from figlib import style  # noqa: E402
from figlib.style import mb  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec  # noqa: E402

CHROM = ann.CHROM
REGION = (22_000_000, 29_000_000)
WINDOW, STEP, RUN = 1_000, 100, 5
ENTER, EXIT = 0.4, 0.1
PAIR_BP = 10_000          # control segments pair when their entries are this close
CONVERGE_BP = 5_000       # PWS/AS call within this of the control call = convergent
LOCAL_FLANK = 30_000      # sensitivity analyses rerun on core +/- this
BOOT_FLANK, BOOT_N, BOOT_SEED = 20_000, 100, 1729
GROUPS = ("Control", "PWS-DEL", "AS-DEL")
GROUP_LABEL = {"Control": "Controls |mat − pat|", "PWS-DEL": "PWS-DEL retained maternal",
               "AS-DEL": "AS-DEL retained paternal"}


def color(group: str) -> str:
    return cohort_lib.COLOR[group]


# ------------------------------------------------------------------ inputs
def orientation(paths, samples, cis) -> dict[str, dict[str, str]]:
    """{sample: {'maternal': 'hap1'|'hap2', 'paternal': ..., 'source': ...}}."""
    table = paths.figures / "figure_1" / "tables" / "Figure1_parental_like_assignment.tsv"
    out = {}
    if table.is_file():
        t = pd.read_csv(table, sep="\t")
        for s in samples:
            rows = t[(t["sample_id"] == s) & t["haplotype_label"].isin(["hap1", "hap2"])]
            mat = rows.loc[rows["parental_assignment"] == "maternal", "haplotype_label"].tolist()
            pat = rows.loc[rows["parental_assignment"] == "paternal", "haplotype_label"].tolist()
            if len(mat) == 1 and len(pat) == 1:
                out[s] = {"maternal": mat[0], "paternal": pat[0], "source": "Figure 1 assignment"}
    for s in samples:
        if s in out:
            continue
        means = {}
        for kind in ("hap1", "hap2"):
            path = cis.find_track(paths.methylation, s, kind)
            t = cis.read_track(path, CHROM, *ann.IC) if path else None
            means[kind] = float(np.average(t.beta, weights=t.coverage)) if t is not None and len(t.beta) else np.nan
        if np.all(np.isfinite(list(means.values()))):
            mat = max(means, key=means.get)
            out[s] = {"maternal": mat, "paternal": "hap2" if mat == "hap1" else "hap1",
                      "source": "IC methylation (higher = maternal)"}
    return out


def load_signals(paths, cohort, include_22q: bool, signed: bool, cis):
    """Per-CpG signals of every genome used: DataFrame(sample, group, pos, signal)."""
    lo, hi = REGION[0] - 50_000, REGION[1] + 50_000
    controls = list(cohort.of("Control")) + (list(cohort.of("DiGeorge")) if include_22q else [])
    orient = orientation(paths, controls, cis)
    rows, files = [], []

    def track(sample, kind):
        path = cis.find_track(paths.methylation, sample, kind)
        if path is None:
            return None
        files.append({"sample_id": sample, "track": kind, "path": str(path)})
        t = cis.read_track(path, CHROM, lo, hi)
        return pd.Series(t.beta, index=t.position)

    combined = []
    for s in controls:
        if s not in orient:
            print(f"[figure 3] {s}: no parental orientation for hap1/hap2; skipped")
            continue
        m, p = track(s, orient[s]["maternal"]), track(s, orient[s]["paternal"])
        c = track(s, "combined")
        if m is None or p is None:
            continue
        j = pd.concat([m.rename("m"), p.rename("p")], axis=1, join="inner")
        sig = (j["m"] - j["p"]) if signed else (j["m"] - j["p"]).abs()
        rows.append(pd.DataFrame({"sample": s, "group": "Control", "pos": sig.index, "signal": sig.values}))
        if c is not None:
            combined.append(c.rename(s))
    if not combined:
        raise SystemExit("no control combined tracks: the baseline cannot be built")
    baseline = pd.concat(combined, axis=1).mean(axis=1, skipna=True)
    for group, sign in (("PWS-DEL", 1), ("AS-DEL", -1)):
        for s in cohort.of(group):
            c = track(s, "combined")
            if c is None:
                continue
            j = pd.concat([c.rename("c"), baseline.rename("b")], axis=1, join="inner")
            rows.append(pd.DataFrame({"sample": s, "group": group, "pos": j.index,
                                      "signal": sign * (j["c"] - j["b"]).values}))
    sig = pd.concat(rows, ignore_index=True)
    sig["pos"] = sig["pos"].astype(np.int64)
    orient_table = pd.DataFrame([{"sample_id": s, **o} for s, o in orient.items()])
    return sig.sort_values(["sample", "pos"]).reset_index(drop=True), pd.DataFrame(files), orient_table


# ------------------------------------------------------------------ windows and calls
def windows_for(pos: np.ndarray, val: np.ndarray, starts: np.ndarray, ends: np.ndarray, min_cpgs: int = 1):
    i0 = np.searchsorted(pos, starts, "left")
    i1 = np.searchsorted(pos, ends, "left")
    cs = np.concatenate([[0.0], np.cumsum(val)])
    n = i1 - i0
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = (cs[i1] - cs[i0]) / n
    mean[n < min_cpgs] = np.nan
    return mean, n


@dataclass
class Spec:
    name: str
    window: int = WINDOW
    step: int = STEP
    enter: float = ENTER
    exit: float = EXIT
    cpg_window: int = 0          # > 0: windows of this many consecutive CpGs


def grid(spec: Spec, region: tuple[int, int], all_pos: np.ndarray):
    if spec.cpg_window:
        p = all_pos[(all_pos >= region[0]) & (all_pos < region[1])]
        idx = np.arange(0, max(0, len(p) - spec.cpg_window), max(1, spec.cpg_window // 5))
        return p[idx], p[idx + spec.cpg_window - 1] + 1
    starts = np.arange(region[0], region[1] - spec.window + 1, spec.step)
    return starts, starts + spec.window


def profiles(sig: pd.DataFrame, spec: Spec, region: tuple[int, int]) -> pd.DataFrame:
    """Window means per genome: sample, group, start, end, mid, n_cpg, value."""
    all_pos = np.unique(sig["pos"].to_numpy())
    starts, ends = grid(spec, region, all_pos)
    out = []
    for (s, g), d in sig.groupby(["sample", "group"], sort=False):
        v, n = windows_for(d["pos"].to_numpy(), d["signal"].to_numpy(float), starts, ends)
        out.append(pd.DataFrame({"sample": s, "group": g, "start": starts, "end": ends,
                                 "mid": (starts + ends) // 2, "n_cpg": n, "value": v}))
    return pd.concat(out, ignore_index=True)


def group_profile(prof: pd.DataFrame) -> pd.DataFrame:
    return (prof.groupby(["group", "start", "end", "mid"], sort=False)["value"]
            .agg(value="mean", sd="std", n_samples="count").reset_index())


def call(start: np.ndarray, end: np.ndarray, value: np.ndarray, enter: float, exit_: float) -> list[dict]:
    """Segments: a run of >= RUN windows above `enter`; exit = end of the first window of the
    first later run of >= RUN windows below `exit_`."""
    high = np.nan_to_num(value, nan=-np.inf) > enter
    low = np.nan_to_num(value, nan=np.inf) < exit_
    segs, i, n = [], 0, len(value)
    while i < n:
        if not high[i]:
            i += 1
            continue
        j = i
        while j < n and high[j]:
            j += 1
        if j - i >= RUN:
            k, exit_pos = j, np.nan
            while k < n:
                if low[k]:
                    m = k
                    while m < n and low[m]:
                        m += 1
                    if m - k >= RUN:
                        exit_pos = float(end[k])
                        break
                    k = m
                else:
                    k += 1
            segs.append({"entry": float(start[i]), "exit": exit_pos, "max_signal": float(np.nanmax(value[i:j]))})
            i = max(j, k if np.isfinite(exit_pos) else j)
        else:
            i = j
    return segs


def consensus(sample_calls: dict[str, list[dict]]) -> list[dict]:
    """Control segments grouped by entry (within PAIR_BP of the first genome's segment)."""
    names = [s for s in sample_calls if sample_calls[s]]
    if not names:
        return []
    out = []
    for seg in sample_calls[names[0]]:
        members = [seg]
        for other in names[1:]:
            near = [x for x in sample_calls[other] if abs(x["entry"] - seg["entry"]) <= PAIR_BP]
            if near:
                members.append(min(near, key=lambda x: abs(x["entry"] - seg["entry"])))
        entries = np.array([m["entry"] for m in members])
        exits = np.array([m["exit"] for m in members], float)
        out.append({"entry": float(entries.mean()), "exit": float(np.nanmean(exits)) if np.isfinite(exits).any() else np.nan,
                    "entry_low": float(entries.min()), "entry_high": float(entries.max()),
                    "exit_low": float(np.nanmin(exits)) if np.isfinite(exits).any() else np.nan,
                    "exit_high": float(np.nanmax(exits)) if np.isfinite(exits).any() else np.nan,
                    "n_genomes": len(members), "high_confidence":
                        len(members) == len(names) and len(names) > 1 and (entries.max() - entries.min()) <= 2_000})
    return out


def _d(a, b) -> float:
    return abs(a - b) if np.isfinite(a) and np.isfinite(b) else 1e9


def primary(ctrl: list[dict], group_calls: dict[str, list[dict]]) -> dict | None:
    """Control segment agreeing best with the PWS-DEL and AS-DEL group calls, and those calls."""
    if not ctrl:
        return None
    def nearest(seg, calls):
        return min(calls, key=lambda x: _d(x["entry"], seg["entry"]) + 0.5 * _d(x["exit"], seg["exit"])) if calls else None
    def score(seg):
        sc = 0.0
        for g in ("PWS-DEL", "AS-DEL"):
            n = nearest(seg, group_calls.get(g, []))
            sc += (_d(n["entry"], seg["entry"]) + 0.5 * _d(n["exit"], seg["exit"])) if n else 1e7
        return sc - (1_000 if seg["high_confidence"] else 0)
    best = min(ctrl, key=score)
    out = {"Control": best}
    for g in ("PWS-DEL", "AS-DEL"):
        out[g] = nearest(best, group_calls.get(g, []))
    return out


def shared_core(prim: dict | None) -> tuple[float, float]:
    if not prim or any(prim.get(g) is None for g in GROUPS):
        return np.nan, np.nan
    s = max(prim[g]["entry"] for g in GROUPS)
    e = min(prim[g]["exit"] for g in GROUPS)
    return (float(s), float(e)) if np.isfinite(e) and e > s else (np.nan, np.nan)


def run_calls(sig: pd.DataFrame, spec: Spec, region: tuple[int, int]):
    prof = profiles(sig, spec, region)
    gp = group_profile(prof)
    sample_calls = {s: call(d["start"].to_numpy(), d["end"].to_numpy(), d["value"].to_numpy(), spec.enter, spec.exit)
                    for s, d in prof.groupby("sample", sort=False)}
    group_calls = {g: call(d["start"].to_numpy(), d["end"].to_numpy(), d["value"].to_numpy(), spec.enter, spec.exit)
                   for g, d in gp.groupby("group", sort=False)}
    ctrl_samples = prof.loc[prof["group"] == "Control", "sample"].unique()
    ctrl = consensus({s: sample_calls[s] for s in ctrl_samples})
    prim = primary(ctrl, group_calls)
    return prof, gp, sample_calls, group_calls, ctrl, prim, shared_core(prim)


def change_point(gp: pd.DataFrame, core: tuple[float, float]) -> tuple[float, float]:
    ent, ex = [], []
    for g in GROUPS:
        d = gp[(gp["group"] == g) & (gp["mid"] >= core[0] - 10_000) & (gp["mid"] <= core[1] + 10_000)]
        if len(d) < 8:
            return np.nan, np.nan
        y = d["value"].rolling(11, center=True, min_periods=3).mean().interpolate().to_numpy()
        grad = np.gradient(y, d["mid"].to_numpy())
        i = int(np.nanargmax(grad))
        j = i + int(np.nanargmin(grad[i:]))
        ent.append(float(d["mid"].iloc[i]))
        ex.append(float(d["mid"].iloc[j]))
    return max(ent), min(ex)


def bootstrap(sig: pd.DataFrame, core: tuple[float, float]) -> pd.DataFrame:
    region = (int(core[0] - BOOT_FLANK), int(core[1] + BOOT_FLANK))
    local = sig[(sig["pos"] >= region[0] - WINDOW) & (sig["pos"] < region[1] + WINDOW)]
    rng = np.random.default_rng(BOOT_SEED)
    rows = []
    for b in range(BOOT_N):
        parts = []
        for _, d in local.groupby("sample", sort=False):
            idx = np.sort(rng.integers(0, len(d), len(d)))
            parts.append(d.iloc[idx])
        res = run_calls(pd.concat(parts, ignore_index=True), Spec("bootstrap"), region)[-1]
        rows.append({"replicate": b, "start": res[0], "end": res[1]})
    return pd.DataFrame(rows)


def sensitivity(sig: pd.DataFrame, core: tuple[float, float], gp_primary: pd.DataFrame) -> pd.DataFrame:
    region = (int(core[0] - LOCAL_FLANK), int(core[1] + LOCAL_FLANK))
    local = sig[(sig["pos"] >= region[0] - 3_000) & (sig["pos"] < region[1] + 3_000)]
    specs = [Spec("250-bp windows", 250), Spec("500-bp windows", 500), Spec("1-kb windows (primary)"),
             Spec("2-kb windows", 2_000), Spec("threshold 0.4/0.15", exit=0.15), Spec("threshold 0.3/0.1", enter=0.3),
             Spec("threshold 0.5/0.1", enter=0.5), Spec("25-CpG windows", cpg_window=25), Spec("50-CpG windows", cpg_window=50)]
    rows = []
    for spec in specs:
        s, e = run_calls(local, spec, region)[-1]
        rows.append({"method": spec.name, "start": s, "end": e})
    s, e = change_point(gp_primary, core)
    rows.append({"method": "change-point on smoothed profiles", "start": s, "end": e})
    boot = bootstrap(sig, core)
    ok = boot.dropna()
    rows.append({"method": f"CpG bootstrap ({len(ok)}/{len(boot)} replicates)",
                 "start": float(np.quantile(ok["start"], 0.975)) if len(ok) else np.nan,
                 "end": float(np.quantile(ok["end"], 0.025)) if len(ok) else np.nan})
    t = pd.DataFrame(rows)
    t["width_bp"] = t["end"] - t["start"]
    inter = (np.minimum(t["end"], core[1]) - np.maximum(t["start"], core[0])).clip(lower=0)
    t["core_covered_pct"] = 100 * inter / (core[1] - core[0])
    t["jaccard"] = inter / (np.maximum(t["end"], core[1]) - np.minimum(t["start"], core[0]))
    return t, boot


# ------------------------------------------------------------------ figure
def smooth_broad(gp: pd.DataFrame, group: str) -> pd.DataFrame:
    d = gp[gp["group"] == group].sort_values("mid").copy()
    win = 901 if group == "Control" else 701
    d["smooth"] = (d["value"].rolling(win, center=True, min_periods=1).median()
                   .rolling(121, center=True, min_periods=1).mean().interpolate())
    return d


def render(outdir: Path, t: dict, meta: dict) -> list[Path]:
    style.setup()
    core = (meta["core_start"], meta["core_end"])
    cmid = (core[0] + core[1]) / 2
    zoom = (core[0] - 4_000, core[1] + 4_000)
    fig = plt.figure(figsize=(style.WIDTH, 9.4))
    gs = GridSpec(5, 1, figure=fig, height_ratios=[1.55, 0.72, 1.25, 1.05, 0.95], hspace=0.76,
                  left=0.2, right=0.86, top=0.955, bottom=0.05)

    # a -- locus-wide contrast
    sub = GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[0], height_ratios=[1, 0.13, 0.26], hspace=0.08)
    ax = fig.add_subplot(sub[0])
    gp = t["group_profiles"]
    for g in GROUPS:
        d = smooth_broad(gp, g)
        ax.plot(mb(d["mid"]), d["smooth"], color=color(g), lw=1.4 if g == "Control" else 1.1, label=GROUP_LABEL[g])
    ax.axhline(0, color=style.LIGHT, lw=0.6, ls=(0, (3, 2)))
    ax.axvline(mb(cmid), color=style.CORE, lw=1.0)
    ax.text(mb(cmid) + 0.03, 0.5, "shared core", transform=ax.get_xaxis_transform(), color=style.CORE,
            fontsize=style.BASE_FONT - 1, va="center", zorder=6,
            bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.85))
    ax.set_xlim(mb(REGION[0]), mb(REGION[1]))
    ax.set_ylabel("parent-of-origin\nsignal (smoothed)")
    ax.tick_params(labelbottom=False)
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.0), fontsize=style.BASE_FONT - 1, ncol=3, borderaxespad=0.2)
    sd = fig.add_subplot(sub[1], sharex=ax)
    ann.draw_sd_track(sd, t["sd_blocks"], (mb(REGION[0]), mb(REGION[1])))
    sd.tick_params(labelbottom=False, bottom=False)
    sd.spines["bottom"].set_visible(False)
    sd.text(-0.01, 0.5, "SD", transform=sd.transAxes, ha="right", va="center", fontsize=style.BASE_FONT - 1, color=style.MUTED)
    gt = fig.add_subplot(sub[2], sharex=ax)
    ann.draw_gene_track(gt, t["genes"], (mb(REGION[0]), mb(REGION[1])))
    gt.set_xlabel("chr15 (Mb, T2T-CHM13)")
    style.panel_label(fig, ax, "a", "Parent-of-origin methylation contrast across chr15:22–29 Mb", dx=-0.15, dy=0.03)

    # b -- boundary intervals
    ax = fig.add_subplot(gs[1])
    iv = t["intervals"]
    rows = list(GROUPS) + ["shared core"]
    for k, lab in enumerate(rows):
        r = iv[iv["label"] == lab]
        if r.empty or not np.isfinite(r["start"].iloc[0]):
            ax.text(mb(zoom[0]), k, "  no call", va="center", fontsize=style.BASE_FONT - 1, color=style.MUTED)
            continue
        s, e = float(r["start"].iloc[0]), float(r["end"].iloc[0])
        c = style.CORE if lab == "shared core" else color(lab)
        ax.barh(k, mb(e) - mb(s), left=mb(s), height=0.58, color=c, alpha=0.9 if lab != "shared core" else 1)
        ax.text(mb(e), k, f" {(e - s) / 1e3:.1f} kb", va="center",
                fontsize=style.BASE_FONT - 1, color=style.INK)
    ax.axvspan(mb(core[0]), mb(core[1]), color=style.CORE, alpha=0.1, lw=0)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(["Controls", "PWS-DEL", "AS-DEL", "Shared core"])
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xlim(mb(zoom[0]), mb(zoom[1]))
    ax.ticklabel_format(axis="x", useOffset=False, style="plain")
    ax.set_xlabel("chr15 (Mb)")
    style.panel_label(fig, ax, "b", f"Boundary calls converge on a {(core[1] - core[0]) / 1e3:.1f}-kb core "
                      f"(chr15:{int(core[0]):,}–{int(core[1]):,})", dx=-0.15)

    # c -- CpG-level evidence
    ax = fig.add_subplot(gs[2])
    pts = t["cpg_points"]
    for g in GROUPS:
        d = pts[pts["group"] == g]
        ax.scatter(mb(d["pos"]), d["signal"], s=2.5, color=color(g), alpha=0.25, lw=0, rasterized=True)
        z = gp[(gp["group"] == g) & (gp["mid"] >= zoom[0]) & (gp["mid"] <= zoom[1])].sort_values("mid")
        ax.plot(mb(z["mid"]), z["value"].rolling(7, center=True, min_periods=3).mean(), color=color(g), lw=1.4,
                label=GROUP_LABEL[g])
    ax.axvspan(mb(core[0]), mb(core[1]), color=style.CORE, alpha=0.1, lw=0)
    for v in (ENTER, EXIT):
        ax.axhline(v, color=style.LIGHT, lw=0.6, ls=(0, (2, 2)))
    ax.text(1.0, ENTER, " enter", transform=ax.get_yaxis_transform(), fontsize=style.BASE_FONT - 1.5, va="center",
            color=style.MUTED)
    ax.text(1.0, EXIT, " exit", transform=ax.get_yaxis_transform(), fontsize=style.BASE_FONT - 1.5, va="center",
            color=style.MUTED)
    ax.set_xlim(mb(zoom[0]), mb(zoom[1]))
    ax.ticklabel_format(axis="x", useOffset=False, style="plain")
    ax.set_ylabel("per-CpG signal")
    ax.set_xlabel("chr15 (Mb)")
    # Reserve the inter-panel whitespace for the legend so it cannot cover the
    # CpG points or the grey enter/exit threshold labels.
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25),
              fontsize=style.BASE_FONT - 1.5, ncol=3, columnspacing=0.8,
              handlelength=1.8, borderaxespad=0)
    style.panel_label(fig, ax, "c", "CpG-level evidence (points) and 1-kb profiles (lines)", dx=-0.15)

    # d -- genomic context
    sub = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[3], width_ratios=[1.25, 1], wspace=0.25)
    ax = fig.add_subplot(sub[0])
    span = (min(ann.BP_CLUSTERS["BP1"][0], core[0]) - 200_000, max(ann.BP_CLUSTERS["BP3"][1], core[1]) + 200_000)
    blocks = t["sd_blocks"]
    for r in blocks.itertuples(index=False):
        if r.end > span[0] and r.start < span[1]:
            ax.add_patch(ann.matplotlib_rect((mb(r.start), 0.35), mb(r.end) - mb(r.start), 0.3, style.SD))
    for name, (s, e) in ann.BP_CLUSTERS.items():
        ax.text(mb((s + e) / 2), 0.72, name, ha="center", fontsize=style.BASE_FONT - 1, fontweight="bold")
    ax.plot([mb(cmid)], [0.5], marker="v", color=style.CORE, ms=6)
    near = meta["nearest_sd_edge"]
    ax.annotate("", xy=(mb(near), 0.2), xytext=(mb(cmid), 0.2),
                arrowprops=dict(arrowstyle="<->", color=style.MUTED, lw=0.7))
    ax.text(mb((near + cmid) / 2), 0.08, f"{abs(cmid - near) / 1e6:.2f} Mb to the nearest duplication",
            ha="center", fontsize=style.BASE_FONT - 1, color=style.MUTED)
    ax.set_xlim(mb(span[0]), mb(span[1]))
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("chr15 (Mb)")
    style.panel_label(fig, ax, "d", "Genomic context of the shared core", dx=-0.15)
    ax = fig.add_subplot(sub[1])
    near_genes = t["core_genes"]
    y = 0.9
    for r in near_genes.itertuples(index=False):
        s, e = max(r.start, zoom[0] - 6_000), min(r.end, zoom[1] + 6_000)
        ax.add_patch(ann.matplotlib_rect((mb(s), y - 0.05), mb(e) - mb(s), 0.1, "#A8B5C7"))
        ax.text(mb(s), y + 0.06, r.gene, fontsize=style.BASE_FONT - 1.5, style="italic", va="bottom")
        y -= 0.2
    reps = t["core_repeats"]
    rep_color = {"maternal_higher": style.MATERNAL, "paternal_higher": style.PATERNAL}
    for r in reps.itertuples(index=False):
        c = rep_color.get(str(r.direction), "#9A9A9A" if str(r.testable) in ("True", "true", "1") else "white")
        ax.add_patch(ann.matplotlib_rect((mb(r.start), 0.25), max(mb(r.end) - mb(r.start), 1e-5), 0.08, c))
        ax.add_patch(plt.Rectangle((mb(r.start), 0.25), max(mb(r.end) - mb(r.start), 1e-5), 0.08, fill=False,
                                   edgecolor="#6F6F6F", lw=0.4))
    if len(reps):
        ax.text(mb(zoom[0] - 6_000), 0.36, "repeat elements (duplicons 07; filled = testable)",
                fontsize=style.BASE_FONT - 2, color=style.MUTED, va="bottom")
    for r in t["icr"].itertuples(index=False):
        ax.add_patch(ann.matplotlib_rect((mb(r.start), 0.04), mb(r.end) - mb(r.start), 0.08, "#7A9EC8"))
        ax.text(mb(r.end), 0.08, f" {r.name}", fontsize=style.BASE_FONT - 1.5, va="center", color=style.MUTED)
    ax.axvspan(mb(core[0]), mb(core[1]), color=style.CORE, alpha=0.18, lw=0)
    ax.set_xlim(mb(zoom[0] - 6_000), mb(zoom[1] + 6_000))
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.ticklabel_format(axis="x", useOffset=False, style="plain")
    ax.xaxis.set_major_locator(plt.MaxNLocator(3))
    ax.set_xlabel("chr15 (Mb)")

    # e -- sensitivity
    ax = fig.add_subplot(gs[4])
    sens = t["sensitivity"]
    for k, r in enumerate(sens.itertuples(index=False)):
        if np.isfinite(r.start) and np.isfinite(r.end):
            ax.barh(k, mb(r.end) - mb(r.start), left=mb(r.start), height=0.55,
                    color=style.CORE if "primary" in r.method else "#6F6F6F")
            ax.text(1.005, k, f"{r.width_bp / 1e3:.1f} kb · {r.core_covered_pct:.0f}%", transform=ax.get_yaxis_transform(),
                    va="center", fontsize=style.BASE_FONT - 1.5, color=style.INK)
        else:
            ax.text(mb(zoom[0]), k, "  no call", va="center", fontsize=style.BASE_FONT - 1.5, color=style.MUTED)
    ax.axvspan(mb(core[0]), mb(core[1]), color=style.CORE, alpha=0.1, lw=0)
    ax.set_yticks(range(len(sens)))
    ax.set_yticklabels(sens["method"], fontsize=style.BASE_FONT - 1)
    ax.set_ylim(len(sens) - 0.5, -0.5)
    ax.set_xlim(mb(zoom[0] - 10_000), mb(zoom[1] + 10_000))
    ax.ticklabel_format(axis="x", useOffset=False, style="plain")
    ax.set_xlabel("chr15 (Mb)")
    ax.text(1.005, -0.9, "width · core\ncovered", transform=ax.get_yaxis_transform(), fontsize=style.BASE_FONT - 1.5,
            color=style.MUTED, va="bottom")
    style.panel_label(fig, ax, "e", "Boundary calls under alternative settings", dx=-0.15)
    return style.save(fig, outdir, "Figure3")


def report(outdir: Path, t: dict, meta: dict, paths_out) -> Path:
    core = (meta["core_start"], meta["core_end"])
    r = style.Report("Figure 3 report: boundary mapping")
    r.p(f"Generated by `scripts/figures/FIGURE_3.py`. Signal mode: {meta['signal_mode']}; controls: "
        f"{', '.join(meta['controls'])} ({meta['orientation_sources']}).")
    r.h("Shared core")
    r.p(f"chr15:{int(core[0]):,}–{int(core[1]):,} ({(core[1] - core[0]) / 1e3:.2f} kb), "
        f"{abs((core[0] + core[1]) / 2 - meta['nearest_sd_edge']) / 1e6:.2f} Mb from the nearest segmental duplication "
        f"({meta['sd_source']}). Overlap with the imprinting centre (chr15:{ann.IC[0]:,}–{ann.IC[1]:,}): "
        f"{meta['ic_overlap_bp']:,} bp.")
    r.table(t["intervals"].assign(width_kb=lambda d: ((d["end"] - d["start"]) / 1e3).round(2)))
    r.p(f"PWS-DEL and AS-DEL calls within {CONVERGE_BP / 1e3:.0f} kb of the control call (entry and exit): "
        f"{meta['convergent']}.")
    r.h("Outside the core")
    r.p("Mean smoothed signal per group, chr15:22–29 Mb excluding core ± 50 kb (a regional parent-of-origin "
        "offset elsewhere would show here):")
    r.table(t["outside"])
    r.h("Repeat elements around the core (scripts/duplicons 07)")
    reps = t["core_repeats"]
    if len(reps):
        r.p("Elements within the imprinting centre +/- 5 kb are not testable in 07 (IC excluded), so their parental "
            "state comes from the CpG analysis of this figure, not from 07.")
        r.table(reps)
    else:
        r.p("No element table from scripts/duplicons 07.")
    r.h("Sensitivity")
    r.table(t["sensitivity"].round({"start": 0, "end": 0, "width_bp": 0, "core_covered_pct": 1, "jaccard": 3}))
    r.h("Caption draft")
    r.p(f"**Boundary mapping of the parent-of-origin methylation transition across 15q11–q13.** "
        f"**a,** Smoothed parent-of-origin signal across chr15:22–29 Mb: |maternal − paternal| in controls "
        f"(n = {len(meta['controls'])}), retained maternal (PWS-DEL) and paternal (AS-DEL) copies relative to the control "
        f"baseline; segmental duplications (BP1–BP3) and genes below. **b,** Boundary intervals called per group converge on "
        f"a {(core[1] - core[0]) / 1e3:.1f}-kb shared core at chr15:{int(core[0]):,}–{int(core[1]):,}. **c,** CpG-level "
        f"signals and 1-kb profiles; dashed lines, entry (0.4) and exit (0.1) thresholds. **d,** The core lies "
        f"{abs((core[0] + core[1]) / 2 - meta['nearest_sd_edge']) / 1e6:.2f} Mb from the nearest duplication and within the "
        f"SNRPN/SNHG14 imprinting-centre interval. **e,** Boundary calls under alternative window sizes, thresholds, "
        f"CpG-count windows, change-point detection and CpG bootstrap; labels give call width and the percentage of the "
        f"core covered.")
    r.p("")
    r.p("Files: " + ", ".join(p.name for p in paths_out))
    return r.write(outdir, "Figure3_report.md")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    paths_lib.add_common_arguments(ap)
    ap.add_argument("--include-22q", action="store_true", help="add the 22q11.2DS genomes to the control references")
    ap.add_argument("--signed", action="store_true", help="controls: maternal - paternal instead of |difference|")
    ap.add_argument("--render-only", action="store_true", help="redraw from the tables of a previous run")
    a = ap.parse_args(argv)
    paths = paths_lib.resolve(a)
    outdir = paths.figure_dir(3, a.outdir)
    names = ["group_profiles", "intervals", "cpg_points", "sensitivity", "sd_blocks", "genes", "core_genes", "icr",
             "outside", "core_repeats"]
    if a.render_only:
        t = {n: style.read_table(outdir, f"Figure3_{n}.tsv") for n in names}
        meta = json.loads((outdir / "tables" / "Figure3_meta.json").read_text())
        paths_out = render(outdir, t, meta)
        print(report(outdir, t, meta, paths_out))
        return
    cis = paths_lib.import_analysis_package()
    cohort = cohort_lib.load(paths.metadata)
    sig, files, orient = load_signals(paths, cohort, a.include_22q, a.signed, cis)
    style.write_table(files, outdir, "Figure3_input_files.tsv")
    style.write_table(orient, outdir, "Figure3_control_orientation.tsv")
    style.write_table(sig, outdir, "Figure3_cpg_signals.tsv.gz")
    prof, gp, sample_calls, group_calls, ctrl, prim, core = run_calls(sig, Spec("primary"), REGION)
    if not np.isfinite(core[0]):
        raise SystemExit("no shared boundary core: see Figure3_sample_calls.tsv and Figure3_group_calls.tsv")
    style.write_table(prof, outdir, "Figure3_window_profiles.tsv.gz")
    style.write_table(pd.DataFrame([{"sample": s, **c} for s, cs in sample_calls.items() for c in cs]), outdir,
                      "Figure3_sample_calls.tsv")
    style.write_table(pd.DataFrame([{"group": g, **c} for g, cs in group_calls.items() for c in cs]), outdir,
                      "Figure3_group_calls.tsv")
    style.write_table(pd.DataFrame(ctrl), outdir, "Figure3_control_consensus.tsv")
    intervals = pd.DataFrame([{"label": g, "start": prim[g]["entry"] if prim[g] else np.nan,
                               "end": prim[g]["exit"] if prim[g] else np.nan} for g in GROUPS]
                             + [{"label": "shared core", "start": core[0], "end": core[1]}])
    sens, boot = sensitivity(sig, core, gp)
    style.write_table(boot, outdir, "Figure3_bootstrap_replicates.tsv")
    blocks, sd_source = ann.sd_blocks(paths)
    cmid = (core[0] + core[1]) / 2
    edges = np.concatenate([blocks["start"].to_numpy(), blocks["end"].to_numpy()])
    nearest_edge = float(edges[np.argmin(np.abs(edges - cmid))]) if len(edges) else np.nan
    genes = ann.landmark_genes(paths, REGION[0], REGION[1])
    core_genes = ann.read_gtf_genes(paths.gtf, CHROM, int(core[0]) - 10_000, int(core[1]) + 10_000)
    core_genes = core_genes[core_genes["gene"].isin(["SNURF", "SNRPN", "SNHG14", "SNURF-SNRPN"])] \
        if len(core_genes) else core_genes
    icr = ann.read_bed(paths.icr_bed)
    icr = icr[(icr["end"] > core[0] - 10_000) & (icr["start"] < core[1] + 10_000)] if len(icr) else icr
    zoom = (core[0] - 4_000, core[1] + 4_000)
    pts = sig[(sig["pos"] >= zoom[0]) & (sig["pos"] <= zoom[1])]
    out_mask = (gp["mid"] < core[0] - 50_000) | (gp["mid"] > core[1] + 50_000)
    outside = gp[out_mask].groupby("group")["value"].agg(mean="mean", median="median").reset_index()
    prim_ok = all(prim.get(g) for g in GROUPS)
    convergent = prim_ok and all(_d(prim[g]["entry"], prim["Control"]["entry"]) <= CONVERGE_BP
                                 and _d(prim[g]["exit"], prim["Control"]["exit"]) <= CONVERGE_BP
                                 for g in ("PWS-DEL", "AS-DEL"))
    core_repeats = pd.DataFrame(columns=["start", "end", "name", "class", "family", "testable", "direction", "delta_direct"])
    es = paths.duplicons / "repeat_methylation" / "element_summary.tsv"
    if es.is_file():
        e = pd.read_csv(es, sep="\t", low_memory=False)
        e = e[(e["end"] > zoom[0] - 6_000) & (e["start"] < zoom[1] + 6_000)]
        core_repeats = e.reindex(columns=core_repeats.columns)
    t = {"group_profiles": gp, "intervals": intervals, "cpg_points": pts, "sensitivity": sens, "sd_blocks": blocks,
         "genes": genes, "core_genes": core_genes, "icr": icr, "outside": outside, "core_repeats": core_repeats}
    meta = {"core_start": core[0], "core_end": core[1], "nearest_sd_edge": nearest_edge, "sd_source": sd_source,
            "signal_mode": "signed" if a.signed else "absolute", "controls": sorted(orient["sample_id"].tolist()),
            "orientation_sources": "; ".join(sorted(orient["source"].unique())) if len(orient) else "none",
            "ic_overlap_bp": int(max(0, min(core[1], ann.IC[1]) - max(core[0], ann.IC[0]))),
            "convergent": "yes" if convergent else "no"}
    for n, df in t.items():
        style.write_table(df, outdir, f"Figure3_{n}.tsv")
    (outdir / "tables" / "Figure3_meta.json").write_text(json.dumps(meta, indent=1))
    paths_out = render(outdir, t, meta)
    print(report(outdir, t, meta, paths_out))


if __name__ == "__main__":
    main()
