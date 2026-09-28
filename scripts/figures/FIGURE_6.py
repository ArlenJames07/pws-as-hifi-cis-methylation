#!/usr/bin/env python3
"""
Figure 6 -- Parent-of-origin methylation of repeat elements across 15q11-q13.

Built from scripts/duplicons 06-08 (nothing is recomputed from raw data):
  a  every testable repeat element of the common CN=1 core: maternal minus paternal
     methylation (PWS-DEL retained minus AS-DEL retained, duplicons 07); direct
     parent-of-origin calls coloured by the higher parent; 250-kb bin medians; the
     largest methylation domain (duplicons 08) shaded; duplications and genes below
  b  number of direct calls against the other 55 labellings of the deletion carriers
  c  agreement of the two designs: direct |difference| vs biparental |H1 - H2| per element
  d  repeat content of haplotype insertions and deletions (>= 50 bp) relative to CHM13 on the
     assembled haplotypes (duplicons 06) -- structural polymorphism, no parental label
  e  repeat types by group: methylation level and the shift inside the domain relative to
     biparental genomes (Control + 22q11.2DS), PWS-mUPD included as an independent genome

Usage: python3 scripts/figures/FIGURE_6.py [--results DIR] [--domain START-END]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figlib import annotation as ann  # noqa: E402
from figlib import cohort as cohort_lib  # noqa: E402
from figlib import paths as paths_lib  # noqa: E402
from figlib import stats  # noqa: E402
from figlib import style  # noqa: E402
from figlib.style import mb  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

GROUPS = ["Biparental", "PWS-mUPD", "PWS-DEL", "AS-DEL"]
GROUP_OF = {"Control": "Biparental", "DiGeorge": "Biparental", "PWS-mUPD": "PWS-mUPD", "PWS-DEL": "PWS-DEL",
            "AS-DEL": "AS-DEL"}
G_COLOR = {"Biparental": style.BIPARENTAL, "PWS-mUPD": cohort_lib.COLOR["PWS-mUPD"],
           "PWS-DEL": cohort_lib.COLOR["PWS-DEL"], "AS-DEL": cohort_lib.COLOR["AS-DEL"]}
G_MARKER = {"Biparental": "o", "PWS-mUPD": cohort_lib.MARKER["PWS-mUPD"], "PWS-DEL": cohort_lib.MARKER["PWS-DEL"],
            "AS-DEL": cohort_lib.MARKER["AS-DEL"]}
TYPE_ORDER = ["Alu", "MIR", "L1", "L2", "LTR", "DNA", "Simple/low compl.", "Satellite"]
BIN = 250_000
EFFECT = 0.2          # scripts/duplicons 07: |direct delta| and biparental |H1 - H2| threshold
SV_HAPLOTYPES = ("hap1", "hap2")   # duplicons 06 also lists 'primary', which repeats hap1/hap2 sequence


def repeat_type(cls: str, fam: str) -> str:
    cls, fam = str(cls), str(fam)
    if cls == "SINE":
        return "Alu" if fam.startswith("Alu") else "MIR" if fam.startswith("MIR") else "other SINE"
    if cls == "LINE":
        return "L1" if fam.startswith("L1") else "L2" if fam.startswith("L2") else "other LINE"
    if cls in ("LTR", "DNA", "Satellite"):
        return cls
    if cls in ("Simple_repeat", "Low_complexity"):
        return "Simple/low compl."
    return "other"


def truthy(series: pd.Series) -> pd.Series:
    return series.astype(str).isin(["True", "true", "1"])


def validation_value(val: pd.DataFrame, prefix: str) -> float:
    hit = val[val["quantity"].astype(str).str.startswith(prefix)] if len(val) else val
    return float(hit["value"].iloc[0]) if len(hit) else np.nan


def domain_of(paths, arg: str | None):
    """(start, end, direction) of the domain: --domain, else the largest of duplicons 08."""
    if arg:
        a, b = arg.replace(",", "").split("-")
        return int(a), int(b), ""
    path = paths.duplicons / "followup" / "methylation_domains.tsv"
    if path.is_file() and path.stat().st_size:
        d = pd.read_csv(path, sep="\t")
        if len(d):
            d = d.assign(size=d["end"] - d["start"]).sort_values("size", ascending=False)
            direction = str(d["direction"].iloc[0]) if "direction" in d else ""
            return int(d["start"].iloc[0]), int(d["end"].iloc[0]), "" if direction == "nan" else direction
    return None


def domain_name(domain) -> str:
    return f"{domain[2].replace('_', '-')} domain" if domain and domain[2] else "methylation domain"


def per_type(ebp: pd.DataFrame, summ: pd.DataFrame, domain) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    s = summ[truthy(summ["testable"]) & truthy(summ["in_common_core"])].copy()
    s["type"] = [repeat_type(c, f) for c, f in zip(s["class"], s["family"])]
    inside = (s["start"] >= domain[0]) & (s["end"] <= domain[1])
    outside = (s["end"] <= domain[0]) | (s["start"] >= domain[1])
    s["region"] = np.where(inside, "inside", np.where(outside, "outside", "edge"))
    s = s[s["region"] != "edge"].set_index("element_id")[["type", "region"]]
    keep = (ebp["mechanism"].isin(["PWS-DEL", "AS-DEL"]) & (ebp["measure"] == "retained_combined")) | \
           (~ebp["mechanism"].isin(["PWS-DEL", "AS-DEL"]) & (ebp["measure"] == "scaffold_combined"))
    e = ebp[keep & ebp["element_id"].isin(s.index)].dropna(subset=["value"]).join(s, on="element_id")
    e["group"] = e["mechanism"].map(GROUP_OF)
    med = (e.groupby(["sample_id", "group", "type", "region"])["value"].agg(median="median", elements="size")
           .reset_index())
    med = med[med["elements"] >= 10]
    w = med.pivot_table(index=["sample_id", "group", "type"], columns="region", values="median").reset_index()
    if not {"inside", "outside"} <= set(w.columns):
        return med, pd.DataFrame(), pd.DataFrame()
    w["in_minus_out"] = w["inside"] - w["outside"]
    ref = w[w["group"] == "Biparental"].groupby("type")["in_minus_out"].median()
    w["shift"] = w["in_minus_out"] - w["type"].map(ref)
    w = w.dropna(subset=["shift"])
    n_in = s[s["region"] == "inside"]["type"].value_counts()
    rows = []
    for t in TYPE_ORDER:
        if n_in.get(t, 0) < 20 or t not in set(w["type"]):
            continue
        x = w[w["type"] == t]
        a = x.loc[x["group"] == "PWS-DEL", "shift"]
        b = x.loc[x["group"] == "AS-DEL", "shift"]
        diff, p, n = stats.label_permutation(a, b)
        u = x.loc[x["group"] == "PWS-mUPD", "shift"]
        bip = x.loc[x["group"] == "Biparental", "shift"]
        rows.append({"type": t, "elements_in_domain": int(n_in.get(t, 0)),
                     "median_shift_PWS-DEL": a.median(), "median_shift_AS-DEL": b.median(),
                     "median_shift_PWS-mUPD": u.median() if len(u) else np.nan,
                     "mupd_below_n_biparental": int((bip > u.median()).sum()) if len(u) else np.nan,
                     "n_biparental": len(bip), "pws_minus_as_mean": diff, "label_permutation_p": p,
                     "labellings": n})
    tests = pd.DataFrame(rows)
    if len(tests):
        tests["q"] = stats.bh(tests["label_permutation_p"])
    return med, w, tests


def load(paths, domain) -> dict:
    d = paths.duplicons / "repeat_methylation"
    need = [d / "element_summary.tsv", d / "element_by_participant.tsv.gz", d / "validation.tsv"]
    missing = [str(x) for x in need if not x.is_file()]
    if missing:
        raise SystemExit(f"missing scripts/duplicons 07 outputs: {missing}")
    summ = pd.read_csv(need[0], sep="\t", low_memory=False)
    ebp = pd.read_csv(need[1], sep="\t", low_memory=False)
    val = pd.read_csv(need[2], sep="\t")
    reg = pd.read_csv(d / "direct_delta_by_region.tsv", sep="\t") if (d / "direct_delta_by_region.tsv").is_file() \
        else pd.DataFrame()
    hr = paths.duplicons / "haplotype_repeats" / "cohort_sv_repeats.tsv"
    sv = pd.read_csv(hr, sep="\t", low_memory=False) if hr.is_file() and hr.stat().st_size else pd.DataFrame()
    med, shifts, tests = per_type(ebp, summ, domain) if domain else (pd.DataFrame(),) * 3
    return {"summary": summ, "validation": val, "regions": reg, "sv": sv, "medians": med, "shifts": shifts,
            "tests": tests}


# ------------------------------------------------------------------ figure
def render(outdir: Path, t: dict, paths, domain) -> list[Path]:
    style.setup()
    summ, val = t["summary"], t["validation"]
    offset = validation_value(val, "parental offset removed")
    core = summ[truthy(summ["in_common_core"]) & truthy(summ["testable"])]
    lo, hi = int(summ["start"].min()), int(summ["end"].max())
    fig = plt.figure(figsize=(style.WIDTH, 10.4))
    gs = GridSpec(7, 1, figure=fig, height_ratios=[1.6, 0.28, 0.42, 1.05, 1.45, 1.2, 2.6], hspace=0.0,
                  left=0.1, right=0.97, top=0.955, bottom=0.045)

    # a -- element map
    ax = fig.add_subplot(gs[0])
    raw = core["delta_direct_raw"] if "delta_direct_raw" in core else core["delta_direct"] + offset
    called = truthy(core["direct_parent_of_origin"])
    ax.scatter(mb((core["start"] + core["end"]) / 2)[~called.to_numpy()], raw[~called], s=1.5, color=style.LIGHT,
               lw=0, rasterized=True, zorder=2)
    for direction, col in (("maternal_higher", style.MATERNAL), ("paternal_higher", style.PATERNAL)):
        sel = called & (core["direction"].astype(str) == direction)
        ax.scatter(mb((core.loc[sel, "start"] + core.loc[sel, "end"]) / 2), raw[sel], s=9, color=col,
                   edgecolor="white", lw=0.3, zorder=4, label=f"direct call, {direction.replace('_', ' ')} "
                                                               f"(n = {int(sel.sum())})")
    reg = t["regions"]
    if len(reg):
        x = mb((reg["bin_start"] + reg["bin_end"]) / 2)
        ax.step(x, reg["observed_median_delta"] + offset, where="mid", color=style.INK, lw=1.1, zorder=5,
                label="250-kb median")
    if domain:
        ax.axvspan(mb(domain[0]), mb(domain[1]), color=style.DOMAIN, alpha=0.8, lw=0, zorder=0)
        ax.text(mb(domain[0]) + 0.02, 0.97, domain_name(domain), transform=ax.get_xaxis_transform(),
                fontsize=style.BASE_FONT - 1.5, color=style.MUTED, va="top")
    c0, c1 = core["start"].min(), core["end"].max()
    for a0, a1 in ((lo, c0), (c1, hi)):
        if a1 > a0:
            ax.axvspan(mb(a0), mb(a1), color="#F2F2F2", lw=0, zorder=0)
    ax.axhline(0, color=style.MUTED, lw=0.6, zorder=1)
    ax.axvline(mb(sum(ann.IC) / 2), color=style.INK, lw=0.6, ls=(0, (2, 2)))
    lim = float(np.nanpercentile(np.abs(raw), 99.5)) if raw.notna().any() else 0.5
    if called.any():
        lim = max(lim, float(np.nanmax(np.abs(raw[called]))))
    lim *= 1.08
    ax.set_ylim(-lim, lim)
    ax.set_xlim(mb(lo), mb(hi))
    ax.set_ylabel("maternal − paternal\nmethylation")
    ax.tick_params(labelbottom=False)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=style.BASE_FONT - 1.5, borderaxespad=0.2)
    sd = fig.add_subplot(gs[1], sharex=ax)
    blocks, _ = ann.sd_blocks(paths)
    ann.draw_sd_track(sd, blocks, (mb(lo), mb(hi)))
    sd.tick_params(labelbottom=False, bottom=False)
    sd.spines["bottom"].set_visible(False)
    gt = fig.add_subplot(gs[2], sharex=ax)
    ann.draw_gene_track(gt, ann.landmark_genes(paths, lo, hi), (mb(lo), mb(hi)), fontsize=style.BASE_FONT - 2)
    gt.set_xlabel("chr15 (Mb, T2T-CHM13); grey background, outside the common CN=1 core (not tested)")
    style.panel_label(fig, ax, "a", "Repeat elements: maternal (PWS-DEL) − paternal (AS-DEL) methylation", dx=-0.08,
                      dy=0.035)

    # b, c
    sub = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[4], width_ratios=[0.8, 1.2, 1.0], wspace=0.55)
    ax = fig.add_subplot(sub[0])
    obs = float(truthy(summ["direct_parent_of_origin"]).sum())    # the count behind p and FDR (07)
    exp = validation_value(val, "direct calls expected under relabelling")
    pcalls = validation_value(val, "label-permutation p")
    nlab = validation_value(val, "label permutations of the deletion carriers")
    ax.bar([0, 1], [obs, exp], color=[style.TRACE, style.LIGHT], width=0.6)
    for x, v in ((0, obs), (1, exp)):
        ax.text(x, v, f"{v:.0f}" if x == 0 else f"{v:.1f}", ha="center", va="bottom", fontsize=style.BASE_FONT - 1)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["observed", f"other {int(nlab) - 1 if np.isfinite(nlab) else ''}\nlabellings"],
                       fontsize=style.BASE_FONT - 1.5)
    ax.set_ylabel("direct calls")
    ax.set_title(f"p = {style.fmt_p(pcalls)}; FDR ≈ {style.fmt(validation_value(val, 'false discovery estimate'), 2)}",
                 fontsize=style.BASE_FONT - 1.5, color=style.MUTED)
    style.panel_label(fig, ax, "b", "Direct calls vs relabelling", dx=-0.07, dy=0.045)
    ax = fig.add_subplot(sub[1])
    v = core[truthy(core["in_common_core"])]
    if "n_biparental" in v:
        v = v[pd.to_numeric(v["n_biparental"], errors="coerce") >= 3]
    x_ = v["delta_direct"].abs()
    y_ = v["asm_mean"]
    ax.scatter(x_[~truthy(v["direct_parent_of_origin"])], y_[~truthy(v["direct_parent_of_origin"])], s=2,
               color=style.LIGHT, lw=0, rasterized=True)
    for direction, col in (("maternal_higher", style.MATERNAL), ("paternal_higher", style.PATERNAL)):
        sel = truthy(v["direct_parent_of_origin"]) & (v["direction"].astype(str) == direction)
        ax.scatter(x_[sel], y_[sel], s=9, color=col, edgecolor="white", lw=0.3, zorder=3)
    ax.axhline(EFFECT, color=style.MUTED, lw=0.6, ls=(0, (2, 2)))
    ax.axvline(EFFECT, color=style.MUTED, lw=0.6, ls=(0, (2, 2)))
    ax.text(0.02, EFFECT, f"effect {EFFECT}", ha="left", va="bottom", fontsize=style.BASE_FONT - 2,
            color=style.MUTED, transform=ax.get_yaxis_transform())
    ax.set_xlabel("direct |maternal − paternal|")
    ax.set_ylabel("biparental |H1 − H2|")
    ax.set_title(f"agreement {style.fmt(validation_value(val, 'agreement'), 2)} "
                 f"(circular-shift p = {style.fmt_p(validation_value(val, 'circular-shift p'))}); "
                 f"ρ = {style.fmt(validation_value(val, 'Spearman rho'), 2)}", fontsize=style.BASE_FONT - 1.5,
                 color=style.MUTED)
    style.panel_label(fig, ax, "c", "Direct vs allelic design", dx=-0.07, dy=0.045)
    ax = fig.add_subplot(sub[2])
    sv = t["sv"]
    if len(sv):
        sv = sv[sv["hap"].astype(str).isin(SV_HAPLOTYPES)].copy()
        sv["dominant"] = sv["repeat_classes"].fillna(".").astype(str).str.split(":").str[0]
        sv.loc[sv["repeat_fraction"] < 0.8, "dominant"] = "< 80% repeat"
        order = sv["dominant"].value_counts().index[:7].tolist()
        counts = sv[sv["dominant"].isin(order)].groupby(["dominant", "type"]).size().unstack(fill_value=0).reindex(order)
        y = np.arange(len(order))
        left = np.zeros(len(order))
        for typ, col in (("INS", "#4D4D4D"), ("DEL", "#BDBDBD")):
            if typ in counts:
                ax.barh(y, counts[typ], left=left, color=col, height=0.6, label=typ)
                left += counts[typ].to_numpy()
        ax.set_yticks(y)
        ax.set_yticklabels(order, fontsize=style.BASE_FONT - 1.5)
        ax.set_ylim(len(order) - 0.5, -0.5)
        ax.set_xlabel("events (hap1 + hap2)")
        ax.legend(fontsize=style.BASE_FONT - 1.5, loc="lower right")
        ax.set_title(f"{sv['sample'].nunique()} genomes, {len(sv):,} indels ≥ 50 bp (hap1 + hap2)",
                     fontsize=style.BASE_FONT - 1.5,
                     color=style.MUTED)
    else:
        ax.text(0.5, 0.5, "no haplotype_repeats table\n(scripts/duplicons 06)", ha="center", va="center",
                transform=ax.transAxes, color=style.MUTED)
        ax.set_axis_off()
    style.panel_label(fig, ax, "d", "Haplotype indels vs CHM13 by repeat", dx=-0.1, dy=0.045)

    # e -- repeat types by group
    tests, shifts, med = t["tests"], t["shifts"], t["medians"]
    types = tests["type"].tolist() if len(tests) else []
    if types:
        sub = GridSpecFromSubplotSpec(2, len(types), subplot_spec=gs[6], hspace=0.55, wspace=0.12)
        rng = np.random.default_rng(1)
        level = med.assign(w=med["median"] * med["elements"]).groupby(["sample_id", "group", "type"])[["w", "elements"]] \
            .sum().reset_index()
        level["level"] = level["w"] / level["elements"]
        lim = max(0.05, float(np.nanmax(np.abs(shifts["shift"]))) * 1.15)
        first = None
        for j, ty in enumerate(types):
            a = fig.add_subplot(sub[0, j], sharey=first[0] if first else None)
            b = fig.add_subplot(sub[1, j], sharey=first[1] if first else None)
            first = first or (a, b)
            for i, g in enumerate(GROUPS):
                for axx, vals in ((a, level.loc[(level["type"] == ty) & (level["group"] == g), "level"]),
                                  (b, shifts.loc[(shifts["type"] == ty) & (shifts["group"] == g), "shift"])):
                    vv = vals.to_numpy()
                    if not len(vv):
                        continue
                    axx.scatter(i + rng.uniform(-0.15, 0.15, len(vv)), vv, s=10, marker=G_MARKER[g], color=G_COLOR[g],
                                edgecolor="white", lw=0.3, zorder=3)
                    if len(vv) > 1:
                        axx.hlines(np.median(vv), i - 0.3, i + 0.3, color=style.INK, lw=1, zorder=4)
            r = tests[tests["type"] == ty].iloc[0]
            a.set_title(f"{ty}\nn = {int(r['elements_in_domain'])}", fontsize=style.BASE_FONT - 1.5)
            b.set_title(f"q = {style.fmt_p(r['q'])}", fontsize=style.BASE_FONT - 1.5, color=style.MUTED)
            b.axhline(0, color=style.MUTED, lw=0.6)
            for axx in (a, b):
                axx.set_xticks(range(len(GROUPS)))
                axx.set_xticklabels(["Bip", "UPD", "PW", "AS"] if axx is b else [], fontsize=style.BASE_FONT - 2.5)
                axx.set_xlim(-0.6, len(GROUPS) - 0.4)
                if j:
                    axx.tick_params(labelleft=False)
            a.set_ylim(0, 1.02)
            b.set_ylim(-lim, lim)
        first[0].set_ylabel("methylation", fontsize=style.BASE_FONT - 1)
        first[1].set_ylabel("domain shift", fontsize=style.BASE_FONT - 1)
        handles = [Line2D([], [], marker=G_MARKER[g], ls="", ms=4.5, color=G_COLOR[g],
                          label={"Biparental": "biparental", "PWS-mUPD": "PWS-mUPD (2 × mat.)",
                                 "PWS-DEL": "PWS-DEL (mat.)", "AS-DEL": "AS-DEL (pat.)"}[g]) for g in GROUPS]
        first[0].legend(handles=handles, loc="lower left", bbox_to_anchor=(0, 1.22), ncol=4, fontsize=style.BASE_FONT - 1.5)
        style.panel_label(fig, first[0], "e", "Repeat types by group: level (top) and shift in the domain (bottom)",
                          dx=-0.08, dy=0.06)
    return style.save(fig, outdir, "Figure6")


def report(outdir, t, domain, paths_out) -> Path:
    r = style.Report("Figure 6 report: repeat-element methylation")
    r.p("Generated by `scripts/figures/FIGURE_6.py` from scripts/duplicons 06, 07 and 08.")
    r.h("Validation (duplicons 07 validation.tsv)")
    r.table(t["validation"])
    r.h("Domain")
    r.p(f"chr15:{domain[0]:,}–{domain[1]:,} ({domain_name(domain)})" if domain
        else "no domain (duplicons 08 methylation_domains.tsv missing)")
    r.h("Repeat types (panel e)")
    r.p("Shift = (inside − outside) − biparental median, per genome; p = exact label permutation of PWS-DEL vs AS-DEL "
        "shifts (difference of means); q = Benjamini–Hochberg across types. A regional effect shifts every type; a "
        "type-specific one only some. The domain was found from the same PWS-DEL and AS-DEL genomes, so these q values "
        "describe how the effect is spread across repeat types rather than test it; PWS-mUPD, which took no part in "
        "finding the domain, is the independent check (mupd_below_n_biparental).")
    r.table(t["tests"])
    if len(t["sv"]):
        r.h("Haplotype indels (panel d)")
        sv = t["sv"][t["sv"]["hap"].astype(str).isin(SV_HAPLOTYPES)].copy()
        sv["dominant"] = sv["repeat_classes"].fillna(".").astype(str).str.split(":").str[0]
        sv.loc[sv["repeat_fraction"] < 0.8, "dominant"] = "< 80% repeat"
        r.p("hap1 and hap2 assemblies only (the primary assembly repeats their sequence); dominant class of "
            "events with >= 80% repeat sequence.")
        r.table(sv.groupby(["type", "dominant"]).size().rename("events").reset_index())
    val = t["validation"]
    r.h("Caption draft")
    n_calls = int(truthy(t["summary"]["direct_parent_of_origin"]).sum())
    offset = validation_value(val, "parental offset removed")
    r.p(f"**Parent-of-origin methylation of repeat elements across 15q11–q13.** **a,** Maternal (PWS-DEL retained) minus "
        f"paternal (AS-DEL retained) methylation for each testable repeat element of the common CN=1 core; coloured, "
        f"elements with complete separation of the two groups and |Δ| ≥ {EFFECT} after removing the genome-wide "
        f"parental offset ({style.fmt(offset, 3, signed=True)}; direct calls); line, 250-kb medians; shaded, the "
        f"{domain_name(domain)}. **b,** Direct calls ({n_calls}) "
        f"compared with the mean of the other labellings of the deletion carriers "
        f"({validation_value(val, 'direct calls expected under relabelling'):.1f}; p = "
        f"{style.fmt_p(validation_value(val, 'label-permutation p'))}). **c,** Agreement between the direct design and "
        f"allelic methylation in biparental genomes (agreement {style.fmt(validation_value(val, 'agreement'), 2)}, "
        f"circular-shift p = {style.fmt_p(validation_value(val, 'circular-shift p'))}). **d,** Repeat content "
        f"of insertions and deletions (≥ 50 bp) of the assembled haplotypes relative to T2T-CHM13; sequence carries no "
        f"parental label, so these are structural polymorphisms. **e,** Methylation level and domain shift by repeat "
        f"type and group; PWS-mUPD did not take part in finding the domain.")
    r.p("")
    r.p("Files: " + ", ".join(p.name for p in paths_out))
    return r.write(outdir, "Figure6_report.md")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    paths_lib.add_common_arguments(ap)
    ap.add_argument("--domain", help="START-END (default: largest domain of scripts/duplicons 08)")
    a = ap.parse_args(argv)
    paths = paths_lib.resolve(a)
    outdir = paths.figure_dir(6, a.outdir)
    domain = domain_of(paths, a.domain)
    t = load(paths, domain)
    style.write_table(t["tests"], outdir, "Figure6e_repeat_type_tests.tsv")
    style.write_table(t["shifts"], outdir, "Figure6e_repeat_type_shifts.tsv")
    style.write_table(t["medians"], outdir, "Figure6e_repeat_type_levels.tsv")
    paths_out = render(outdir, t, paths, domain)
    print(report(outdir, t, domain, paths_out))


if __name__ == "__main__":
    main()
