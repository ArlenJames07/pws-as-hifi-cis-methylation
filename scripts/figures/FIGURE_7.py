#!/usr/bin/env python3
"""
Figure 7 -- How the 15q11-q13 deletions formed.

Built from scripts/duplicons 01, 03, 05 and 08 and the phased BAMs:
  a  directly oriented duplication pairs of the window (arcs; the pairs that link the two
     edges of a deletion in black) above the deletion interval of every carrier, coloured by
     mechanism
  b  homologous-position test: copy number along copy A (proximal) and along copy B
     (distal, projected onto copy A) of the pair at each carrier's edges. Non-allelic
     homologous recombination leaves the two transitions at the same homologous position
  c  evidence for each deletion (Table 2): edges at homologous positions, junction-read
     excess at the crossover, an assembled contig across the deletion, a contig joining an
     edge to another locus, split reads in the carrier and none in the other genomes
  d  assembled junctions: alignment blocks of each contig that crosses a deletion, drawn
     from the contig (top) to the reference (bottom, the two edges side by side)
  e  non-recurrent junction: copy number at the edge, the contig junction to the partner
     locus (microhomology / insertion) and split reads in every genome of the cohort

Mechanism (figure level), from the evidence, in this order:
  contig junction to another locus at an edge      -> non-recurrent junction
  duplicons 08 status 'junction resolved at bp level' -> junction at bp (split reads)
  'confirmed NAHR' / 'compatible with NAHR'         -> NAHR, confirmed / compatible
  'confirmed by assembly'                           -> assembled, switch not homologous
  'edges near one SD pair, not homologous'          -> SD pair, edges not homologous
  otherwise                                         -> unresolved

Usage: python3 scripts/figures/FIGURE_7.py [--results DIR] [--render-only] [--no-split-reads]
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
from figlib import style  # noqa: E402
from figlib.style import mb  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Arc, Polygon, Rectangle  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

CHROM = ann.CHROM
PLOT = (17_500_000, 33_000_000)
FLANK = 100_000                  # bins drawn on each side of a duplication copy (panel b)
EDGE_VIEW = 200_000              # half-width of the copy-number view at a non-recurrent edge (panel e)
MIN_PAIR_BP = 5_000              # as scripts/duplicons read_direct_pairs
MIN_BLOCK_BP = 5_000             # contig blocks drawn in panel d
HOMOLOGOUS_BP = 20_000           # scripts/duplicons limit for 'homologous positions'
MAX_B_PANELS = 6

MECHANISMS = ["NAHR, confirmed", "NAHR, compatible", "junction at bp (split reads)", "non-recurrent junction",
              "assembled, switch not homologous", "SD pair, edges not homologous", "unresolved"]
MECH_COLOR = {"NAHR, confirmed": "#1B7837", "NAHR, compatible": "#A6DBA0",
              "junction at bp (split reads)": "#00441B", "non-recurrent junction": style.TRACE,
              "assembled, switch not homologous": "#C2A5CF", "SD pair, edges not homologous": "#80CDC1",
              "unresolved": style.LIGHT}
COPY_A, COPY_B = style.INK, style.TRACE


def truthy(v) -> bool:
    return str(v) in ("True", "true", "1", "1.0")


def num(v) -> float:
    try:
        x = float(v)
    except (TypeError, ValueError):
        return np.nan
    return x


def display_pair(value) -> str:
    """Turn internal SD-pair codes into reader-facing labels."""
    text = str(value)
    if text.startswith("SDpair_"):
        return f"SD pair {text.removeprefix('SDpair_')}"
    if text.startswith("SDpair"):
        return text.replace("SDpair", "SD pair ", 1).replace("_", " ").strip()
    return text


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False) if path.is_file() and path.stat().st_size else pd.DataFrame()


# ------------------------------------------------------------------ inputs
def direct_pairs(path: Path) -> pd.DataFrame:
    """Unique, directly oriented duplication pairs (copy A upstream of copy B), as
    scripts/duplicons junction_reads.read_direct_pairs."""
    cols = ["pair", "a0", "a1", "b0", "b1", "identity"]
    if not path.is_file():
        return pd.DataFrame(columns=cols)
    seen, rows = set(), []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 6 or line.startswith("#") or f[0] != CHROM or f[3] != CHROM:
            continue
        a, b = (int(f[1]), int(f[2])), (int(f[4]), int(f[5]))
        if len(f) >= 10 and f[8] != f[9]:
            continue
        if min(a[1] - a[0], b[1] - b[0]) < MIN_PAIR_BP:
            continue
        a, b = sorted([a, b])
        if a[1] > b[0] or (a, b) in seen:
            continue
        seen.add((a, b))
        rows.append((f[6] if len(f) > 6 else f"pair{len(rows)}", a[0], a[1], b[0], b[1], num(f[7]) if len(f) > 7 else np.nan))
    return pd.DataFrame(rows, columns=cols)


def to_copy_a(pos: float, p: pd.Series) -> float:
    """Position in copy B -> homologous position in copy A (linear, direct orientation)."""
    return p.a0 + (pos - p.b0) * (p.a1 - p.a0) / max(1, p.b1 - p.b0)


def sunk_bins(paths) -> pd.DataFrame:
    path = paths.duplicons / "breakpoints" / "sunk15q.bins.tsv"
    if not path.is_file():
        return pd.DataFrame(columns=["sample_id", "start", "end", "cn_norm", "usable"])
    b = pd.read_csv(path, sep="\t", usecols=["sample", "start", "end", "cn_norm", "usable"], low_memory=False)
    b["usable"] = b["usable"].astype(str).isin(["True", "true", "1"])
    return b.rename(columns={"sample": "sample_id"})


def split_counts(paths, cohort, junctions: pd.DataFrame, notes: list[str]) -> pd.DataFrame:
    """Split reads supporting each non-recurrent junction in every genome of the cohort."""
    if junctions.empty:
        return pd.DataFrame()
    import importlib.util
    if importlib.util.find_spec("pysam") is None:
        notes.append("pysam not available: split reads at non-recurrent junctions not counted "
                     "(conda env envs/figures.yml)")
        return pd.DataFrame()
    from figlib.modbam import find_bam
    rows = []
    for j in junctions.itertuples(index=False):
        for s in cohort.samples:
            bam, _ = find_bam(paths, s)
            n = deletions.split_read_support(bam, CHROM, int(j.chr15_pos), j.partner_chrom, int(j.partner_pos)) \
                if bam is not None else np.nan
            rows.append({"junction": j.junction, "carrier": j.sample_id, "sample_id": s, "label": cohort.label(s),
                         "group": cohort.group(s), "is_carrier": s == j.sample_id, "split_reads": n})
    return pd.DataFrame(rows)


def join_text(microhomology, overlap, inserted) -> str:
    mh, ov, ins = num(microhomology), num(overlap), num(inserted)
    if np.isfinite(mh) and mh > 0:
        return f"{mh:.0f}-bp microhomology"
    if np.isfinite(ov) and ov > 0:
        return f"{ov:,.0f}-bp alignment overlap"
    if np.isfinite(ins) and ins > 0:
        return f"{ins:,.0f}-bp insertion"
    return "blunt join"


def mechanism(status: str, junction: bool) -> str:
    if junction:
        return "non-recurrent junction"
    return {"junction resolved at bp level": "junction at bp (split reads)", "confirmed NAHR": "NAHR, confirmed",
            "compatible with NAHR": "NAHR, compatible",
            "confirmed by assembly": "assembled, switch not homologous",
            "edges near one SD pair, not homologous": "SD pair, edges not homologous"}.get(status, "unresolved")


def build(paths, cohort, count_split: bool = True) -> tuple[dict, dict]:
    notes: list[str] = []
    dels, junctions = deletions.deletion_table(paths, cohort)
    if dels.empty:
        raise SystemExit("no PWS-DEL / AS-DEL carriers in the metadata")
    d = paths.duplicons
    nahr = read(d / "breakpoints" / "nahr_test.tsv")
    status = read(d / "followup" / "breakpoint_status.tsv")
    fusions = read(d / "followup" / "assembled_fusions.tsv")
    if len(fusions):            # reference order: on a minus-strand contig ref_left is the distal edge
        lo_, hi_ = fusions[["ref_left", "ref_right"]].min(axis=1), fusions[["ref_left", "ref_right"]].max(axis=1)
        fusions["ref_left"], fusions["ref_right"] = lo_, hi_
        if "fraction_of_cnv_deletion" in fusions:        # 08 uses the fusion covering most of the deletion
            fusions = fusions.sort_values("fraction_of_cnv_deletion", ascending=False, kind="stable")
    split08 = read(d / "followup" / "split_reads.tsv")
    pairs = direct_pairs(d / "reference" / "window.sd_pairs.bedpe")
    bins = sunk_bins(paths)
    for name, t in (("breakpoints/nahr_test.tsv", nahr), ("followup/breakpoint_status.tsv", status),
                    ("reference/window.sd_pairs.bedpe", pairs), ("breakpoints/sunk15q.bins.tsv", bins)):
        if t.empty:
            notes.append(f"scripts/duplicons {name} missing or empty")

    # non-recurrent junctions: contig junctions that set a deletion edge
    used = []
    for r in dels.itertuples(index=False):
        for side in ("start", "end"):
            if str(getattr(r, f"{side}_source")).startswith("contig junction") and len(junctions):
                j = junctions[(junctions["sample_id"] == r.sample_id) & (junctions["chr15_pos"] == getattr(r, side))]
                if len(j):
                    used.append(j.iloc[[0]].assign(edge="proximal" if side == "start" else "distal"))
    nrj = pd.concat(used, ignore_index=True) if used else pd.DataFrame()
    if len(nrj):
        nrj["junction"] = [f"{r.sample_id} {CHROM}:{int(r.chr15_pos):,}–{r.partner_chrom}:{int(r.partner_pos):,}"
                           for r in nrj.itertuples(index=False)]
    split = split_counts(paths, cohort, nrj, notes) if count_split else pd.DataFrame()
    if not count_split and len(nrj):
        notes.append("--no-split-reads: split reads at non-recurrent junctions not counted")

    # evidence per carrier (Table 2)
    ev = []
    for r in dels.itertuples(index=False):
        s = r.sample_id
        n = nahr[nahr["sample_id"] == s].iloc[0] if len(nahr) and (nahr["sample_id"] == s).any() else {}
        st = status[status["sample"] == s].iloc[0] if len(status) and (status["sample"] == s).any() else {}
        fu = fusions[fusions["sample"] == s] if len(fusions) else pd.DataFrame()
        jn = nrj[nrj["sample_id"] == s] if len(nrj) else pd.DataFrame()
        sp_c = sp_p = np.nan
        if len(jn) and len(split):
            k = split[split["junction"] == jn["junction"].iloc[0]]
            sp_c = float(k.loc[k["is_carrier"], "split_reads"].max())
            sp_p = float(k.loc[~k["is_carrier"], "split_reads"].max()) if (~k["is_carrier"]).any() else np.nan
        if not len(jn):     # deletion-type split reads of 08 (carrier; maximum over the other genomes)
            if len(split08) and {"carrier", "genome", "split_reads_deletion_type"} <= set(split08.columns):
                k = split08[split08["carrier"] == s]
                v = pd.to_numeric(k["split_reads_deletion_type"], errors="coerce").fillna(0)
                if len(k):
                    sp_c = float(v[k["genome"] == s].max()) if (k["genome"] == s).any() else np.nan
                    sp_p = float(v[k["genome"] != s].max()) if (k["genome"] != s).any() else np.nan
        pair = str(n.get("pair", ".")) if len(n) else "."
        fd = num(fu["fusion_homolog_distance_bp"].iloc[0]) if len(fu) and "fusion_homolog_distance_bp" in fu else np.nan
        stat = str(st.get("status", ".")) if len(st) else "."
        ev.append({
            "sample_id": s, "label": r.label, "group": r.group, "deletion_class": r.deletion_class,
            "start": r.start, "start_source": r.start_source, "end": r.end, "end_source": r.end_source,
            "size_mb": r.size_mb, "sd_pair": pair,
            "homologous_distance_bp": num(n.get("distance_bp")) if len(n) else np.nan,
            "homologous_edges": truthy(n.get("consistent_with_NAHR")) if len(n) else False,
            "junction_reads": num(n.get("junction_reads_at_crossover")) if len(n) else np.nan,
            "junction_reads_expected": num(n.get("expected_from_panel_at_crossover")) if len(n) else np.nan,
            "junction_reads_p": num(n.get("poisson_p_at_crossover")) if len(n) else np.nan,
            "assembled_contig": f"{fu['hap'].iloc[0]}:{fu['contig'].iloc[0]}" if len(fu) else ".",
            "assembled_contigs_all": ";".join(f"{a}:{b}" for a, b in zip(fu["hap"], fu["contig"])) if len(fu) else ".",
            "assembled_fusion_pair": str(fu["fusion_pair"].iloc[0]) if len(fu) and "fusion_pair" in fu else ".",
            "assembled_homologous_distance_bp": fd,
            "contig_junction": (f"{jn['partner_chrom'].iloc[0]}:{int(jn['partner_pos'].iloc[0]):,} "
                                f"({jn['assembly'].iloc[0]} {jn['contig'].iloc[0]})") if len(jn) else ".",
            "microhomology_bp": num(jn["microhomology_bp"].iloc[0]) if len(jn) else np.nan,
            "alignment_overlap_bp": num(jn["alignment_overlap_bp"].iloc[0]) if len(jn) and "alignment_overlap_bp" in jn
            else np.nan,
            "inserted_bp": num(jn["inserted_bp"].iloc[0]) if len(jn) else np.nan,
            "split_reads_carrier": sp_c, "split_reads_panel_max": sp_p,
            "duplicons08_status": stat, "crossover_estimate": str(st.get("crossover_estimate", ".")) if len(st) else ".",
            "mechanism": mechanism(stat, len(jn) > 0)})
    evidence = pd.DataFrame(ev)

    # panel b: copy number along the two copies of each carrier's pair
    hom, hom_meta = [], []
    for e in evidence.itertuples(index=False):
        if e.sd_pair in (".", "nan") or pairs.empty or not (pairs["pair"] == e.sd_pair).any() or bins.empty:
            continue
        p = pairs[pairs["pair"] == e.sd_pair].iloc[0]
        b = bins[(bins["sample_id"] == e.sample_id)]
        mid = (b["start"] + b["end"]) / 2
        for copy, lo, hi in (("A", p.a0, p.a1), ("B", p.b0, p.b1)):
            x = b[(mid >= lo - FLANK) & (mid <= hi + FLANK)].copy()
            x["copy"] = copy
            xm = (x["start"] + x["end"]) / 2
            x["x_copy_a"] = xm if copy == "A" else [to_copy_a(v, p) for v in xm]
            hom.append(x.assign(pair=e.sd_pair)[["sample_id", "pair", "copy", "x_copy_a", "cn_norm", "usable"]])
        n = nahr[nahr["sample_id"] == e.sample_id].iloc[0]
        fu = fusions[fusions["sample"] == e.sample_id] if len(fusions) else pd.DataFrame()
        fu = fu[fu.get("fusion_pair", pd.Series(dtype=str)).astype(str) == e.sd_pair] if len(fu) else fu
        hom_meta.append({"sample_id": e.sample_id, "label": e.label, "group": e.group, "pair": e.sd_pair,
                         "a0": int(p.a0), "a1": int(p.a1), "b0": int(p.b0), "b1": int(p.b1),
                         "proximal_edge": num(n.get("proximal")), "distal_edge": num(n.get("distal")),
                         "distal_edge_in_copy_a": to_copy_a(num(n.get("distal")), p),
                         "distance_bp": e.homologous_distance_bp, "mechanism": e.mechanism,
                         "fusion_left": num(fu["ref_left"].iloc[0]) if len(fu) else np.nan,
                         "fusion_right_in_copy_a": to_copy_a(num(fu["ref_right"].iloc[0]), p) if len(fu) else np.nan,
                         "fusion_distance_bp": num(fu["fusion_homolog_distance_bp"].iloc[0]) if len(fu) else np.nan})
    hom_meta = pd.DataFrame(hom_meta)
    if len(hom_meta):
        order = {m: i for i, m in enumerate(MECHANISMS)}
        hom_meta = hom_meta.assign(o=hom_meta["mechanism"].map(order)).sort_values(["o", "distance_bp"]) \
            .drop(columns="o").head(MAX_B_PANELS)

    # panel d: blocks of the assembled contigs that cross a deletion
    blocks = []
    for f in fusions.itertuples(index=False) if len(fusions) else []:
        path = d / "assembly" / f.sample / f"{f.sample}.blocks.tsv"
        b = read(path)
        if b.empty:
            continue
        b = b[(b["hap"] == f.hap) & (b["contig"] == f.contig) & (b["chrom"] == CHROM) & (b["aln_len"] >= MIN_BLOCK_BP)]
        near_l = (b["t_end"] - f.ref_left).abs().clip(lower=0)
        near_r = (b["t_start"] - f.ref_right).abs()
        b = b[(np.minimum(near_l, near_r) <= 1_000_000)].copy()
        if b.empty:
            continue
        mid = (b["t_start"] + b["t_end"]) / 2
        b["side"] = np.where((mid - f.ref_left).abs() <= (mid - f.ref_right).abs(), "proximal", "distal")
        b["ref_left"], b["ref_right"] = f.ref_left, f.ref_right
        b["fusion_pair"] = getattr(f, "fusion_pair", ".")
        b["fusion_homolog_distance_bp"] = num(getattr(f, "fusion_homolog_distance_bp", np.nan))
        blocks.append(b)
    blocks = pd.concat(blocks, ignore_index=True) if blocks else pd.DataFrame()

    # panel e: copy number around each non-recurrent edge (carrier and the other genomes)
    edge_cn = []
    for j in nrj.itertuples(index=False) if len(nrj) else []:
        m = (bins["end"] >= j.chr15_pos - EDGE_VIEW) & (bins["start"] <= j.chr15_pos + EDGE_VIEW)
        x = bins[m].copy()
        x["junction"] = j.junction
        x["is_carrier"] = x["sample_id"] == j.sample_id
        edge_cn.append(x)
    edge_cn = pd.concat(edge_cn, ignore_index=True) if edge_cn else pd.DataFrame()

    # panel a: pairs to draw
    used_pairs = set(evidence["sd_pair"]) | set(evidence["assembled_fusion_pair"])
    arcs = pairs.copy()
    arcs["used"] = arcs["pair"].isin(used_pairs)
    arcs["length"] = np.minimum(arcs["a1"] - arcs["a0"], arcs["b1"] - arcs["b0"])
    arcs = arcs[(arcs["a0"] < PLOT[1]) & (arcs["b1"] > PLOT[0])]
    background = arcs[~arcs["used"]].sort_values("length", ascending=False).head(300)
    arcs = pd.concat([background, arcs[arcs["used"]]], ignore_index=True)
    sd, sd_source = ann.sd_blocks(paths)
    tables = {"evidence": evidence, "arcs": arcs, "sd_blocks": sd, "hom": pd.concat(hom, ignore_index=True) if hom
              else pd.DataFrame(), "hom_meta": hom_meta, "blocks": blocks, "junctions": nrj, "split": split,
              "edge_cn": edge_cn}
    meta = {"sd_source": sd_source, "notes": notes,
            "groups": cohort.groups_present(), "n_genomes": len(cohort.samples)}
    return tables, meta


# ------------------------------------------------------------------ drawing
def draw_map(fig, spec, t) -> None:
    ev, arcs = t["evidence"], t["arcs"]
    n = len(ev)
    sub = GridSpecFromSubplotSpec(3, 1, subplot_spec=spec, height_ratios=[1.0, 0.28, 0.2 * n + 0.25], hspace=0.05)
    xl = (mb(PLOT[0]), mb(PLOT[1]))
    ax = fig.add_subplot(sub[0])
    span = float((arcs["b0"] + arcs["b1"] - arcs["a0"] - arcs["a1"]).max() / 2) if len(arcs) else 1.0
    placed: list[tuple[float, float]] = []
    for r in arcs.sort_values(["used", "length"]).itertuples(index=False):
        ma, mb_ = (r.a0 + r.a1) / 2e6, (r.b0 + r.b1) / 2e6
        w = mb_ - ma
        h = 1.8 * (w * 1e6) / span
        ax.add_patch(Arc(((ma + mb_) / 2, 0), w, h, theta1=0, theta2=180,
                         color=style.INK if r.used else style.LIGHT, lw=1.1 if r.used else 0.35,
                         alpha=1 if r.used else 0.6, zorder=3 if r.used else 1))
        if r.used:                                  # stack labels of pairs with nearby apices
            x, y = (ma + mb_) / 2, h / 2 + 0.02
            while any(abs(x - px) < 2.0 and abs(y - py) < 0.09 for px, py in placed):
                y += 0.09
            placed.append((x, y))
            ax.text(x, y, display_pair(r.pair), ha="center", va="bottom", fontsize=style.BASE_FONT - 2,
                    color=style.INK, zorder=6,
                    bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none", alpha=0.8))
    ax.set_xlim(*xl)
    ax.set_ylim(0, max([1.12] + [py + 0.12 for _, py in placed]))
    ax.axis("off")
    ax.text(1.0, 1.0, "arcs, directly oriented duplication pairs\n(black, one copy at each edge of a deletion)\n"
            "bars, deletion intervals coloured by mechanism (c)", transform=ax.transAxes, ha="right",
            fontsize=style.BASE_FONT - 2, color=style.MUTED, va="top")
    sd = fig.add_subplot(sub[1], sharex=ax)
    ann.draw_sd_track(sd, t["sd_blocks"], xl, fontsize=style.BASE_FONT - 1.5)
    for r in arcs[arcs["used"]].itertuples(index=False) if len(arcs) else []:
        for s0, s1 in ((r.a0, r.a1), (r.b0, r.b1)):
            sd.add_patch(Rectangle((mb(s0), 0.2), mb(s1) - mb(s0), 0.6, facecolor="#8C7A4E", lw=0, zorder=1.5))
    sd.tick_params(labelbottom=False, bottom=False)
    sd.spines["bottom"].set_visible(False)
    ax_d = fig.add_subplot(sub[2], sharex=ax)
    for k, r in enumerate(ev.itertuples(index=False)):
        c = MECH_COLOR[r.mechanism]
        if np.isfinite(r.start) and np.isfinite(r.end):
            ax_d.add_patch(Rectangle((mb(r.start), k - 0.3), mb(r.end) - mb(r.start), 0.6, facecolor=c, lw=0))
        ax_d.text(mb(PLOT[0]) - 0.1, k, r.label, ha="right", va="center",
                  fontsize=style.BASE_FONT - 1.5, color=cohort_lib.COLOR[r.group], fontweight="bold")
        ax_d.text(mb(PLOT[1]) + 0.1, k, f"{r.deletion_class.split(' (')[0]} · {r.size_mb:.1f} Mb", ha="left",
                  va="center", fontsize=style.BASE_FONT - 2, color=style.INK)
    ax_d.set_ylim(n - 0.5, -0.6)
    ax_d.set_yticks([])
    ax_d.spines["left"].set_visible(False)
    ax_d.set_xlim(*xl)
    ax_d.set_xlabel("chr15 (Mb, T2T-CHM13)")
    style.panel_label(fig, ax, "a", "Duplication pairs and deletion intervals, coloured by mechanism", dx=-0.1, dy=0.02)


def draw_homologous(fig, spec, t) -> None:
    meta, hom = t["hom_meta"], t["hom"]
    if meta.empty:
        ax = fig.add_subplot(spec)
        ax.text(0.5, 0.5, "no carrier with a direct duplication pair at both edges", ha="center", va="center",
                transform=ax.transAxes, color=style.MUTED)
        ax.axis("off")
        style.panel_label(fig, ax, "b", "Homologous-position test", dx=-0.1, dy=0.02)
        return
    k = len(meta)
    sub = GridSpecFromSubplotSpec(1, max(k, 4), subplot_spec=spec, wspace=0.22)
    first = None
    for i, m in enumerate(meta.itertuples(index=False)):
        ax = fig.add_subplot(sub[i], sharey=first)
        first = first or ax
        h = hom[hom["sample_id"] == m.sample_id]
        for copy, col, lab in (("A", COPY_A, "copy A (proximal)"), ("B", COPY_B, "copy B (distal), projected")):
            x = h[(h["copy"] == copy) & h["usable"].astype(str).isin(["True", "true", "1"])].sort_values("x_copy_a")
            xs = (x["x_copy_a"] - m.a0) / 1e3
            ax.scatter(xs, x["cn_norm"], s=1.5, color=col, alpha=0.35, lw=0, rasterized=True)
            ax.plot(xs, x["cn_norm"].rolling(5, center=True, min_periods=2).median(), color=col, lw=1.0, label=lab)
        ax.axvspan(0, (m.a1 - m.a0) / 1e3, color=style.SD, alpha=0.35, lw=0, zorder=0)
        for pos, col, ls in ((m.proximal_edge, COPY_A, "-"), (m.distal_edge_in_copy_a, COPY_B, (0, (2, 1.5)))):
            if np.isfinite(pos):
                ax.axvline((pos - m.a0) / 1e3, color=col, lw=0.8, ls=ls)
        if np.isfinite(m.fusion_left):
            ax.plot([(m.fusion_left - m.a0) / 1e3], [3.35], marker="v", ms=4, color=MECH_COLOR["NAHR, confirmed"],
                    clip_on=False)
        ax.set_ylim(-0.1, 3.5)
        ax.set_yticks([0, 1, 2, 3])
        ax.set_xlim(-FLANK / 1e3, (m.a1 - m.a0 + FLANK) / 1e3)
        ax.axhline(2, color=style.GRID, lw=0.6, zorder=0)
        ax.axhline(1, color=style.GRID, lw=0.6, zorder=0)
        if i:
            ax.tick_params(labelleft=False)
        else:
            ax.set_ylabel("copy number\n(SUNK)", fontsize=style.BASE_FONT - 1)
        ax.tick_params(labelsize=style.BASE_FONT - 2)
        ax.xaxis.set_major_locator(MaxNLocator(3))
        dist = f" · Δ {m.distance_bp / 1e3:.1f} kb" if np.isfinite(m.distance_bp) else ""
        fus = f" · contig {m.fusion_distance_bp:.0f} bp" if np.isfinite(m.fusion_distance_bp) else ""
        ax.set_title(m.label, fontsize=style.BASE_FONT - 1, color=cohort_lib.COLOR[m.group],
                     fontweight="bold", pad=9)
        ax.text(0.5, 1.02, f"{display_pair(m.pair)}{dist}{fus}", transform=ax.transAxes, ha="center", va="bottom",
                fontsize=style.BASE_FONT - 2.2, color=MECH_COLOR[m.mechanism] if m.mechanism != "unresolved"
                else style.MUTED)
    fig.text(0.5, first.get_position().y0 - 0.03, "position in copy A (kb from its start; copy B projected by "
             "homology); shaded, the duplication copy", ha="center", fontsize=style.BASE_FONT - 1)
    first.legend(handles=[Line2D([], [], color=COPY_A, lw=1, label="copy A (proximal)"),
                          Line2D([], [], color=COPY_B, lw=1, label="copy B (distal), projected"),
                          Line2D([], [], color=style.INK, lw=0.8, label="SUNK edges"),
                          Line2D([], [], color=MECH_COLOR["NAHR, confirmed"], marker="v", ls="", ms=4,
                                 label="assembled switch")],
                 loc="lower left", bbox_to_anchor=(0, 1.2), ncol=4, fontsize=style.BASE_FONT - 1.5, frameon=False)
    style.panel_label(fig, first, "b", "Homologous-position test (Δ, distance between the two edges in homologous "
                      "coordinates)", dx=-0.1, dy=0.05)


EVIDENCE_COLS = [("homologous\nedges", 0.0), ("junction reads\n(obs / exp)", 1.12),
                 ("assembled\ncontig", 2.27), ("junction to\nother locus", 3.55),
                 ("split reads (deletion\ntype; carrier / others)", 4.92),
                 ("mechanism", 6.18)]


def draw_evidence(fig, spec, t) -> None:
    ev = t["evidence"]
    ax = fig.add_subplot(spec)
    n = len(ev)
    on = style.INK

    def cell(j, k, state, text):
        x = EVIDENCE_COLS[j][1]
        if state == "yes":
            ax.scatter(x + 0.08, k, s=26, color=on, zorder=3)
        elif state == "partial":
            ax.scatter(x + 0.08, k, s=26, facecolor="white", edgecolor=on, lw=0.9, zorder=3)
        else:
            ax.scatter(x + 0.08, k, s=10, color=style.GRID, zorder=2)
        if text:
            ax.text(x + 0.2, k, text, va="center", ha="left", fontsize=style.BASE_FONT - 2.2, color=style.INK)

    for k, r in enumerate(ev.itertuples(index=False)):
        if k % 2 == 0:
            ax.axhspan(k - 0.5, k + 0.5, color="#F7F7F7", lw=0, zorder=0)
        d = r.homologous_distance_bp
        cell(0, k, "yes" if r.homologous_edges else "partial" if np.isfinite(d) else "no",
             f"{display_pair(r.sd_pair)}\n{d / 1e3:.1f} kb" if np.isfinite(d) else "")
        jr, je, jp = r.junction_reads, r.junction_reads_expected, r.junction_reads_p
        cell(1, k, "yes" if np.isfinite(jp) and jp < 1e-3 else "partial" if np.isfinite(jr) and jr > (je if np.isfinite(je)
             else np.inf) else "no", f"{jr:.0f} / {je:.1f}" + (f"\np = {style.fmt_p(jp)}" if np.isfinite(jp) else "")
             if np.isfinite(jr) and np.isfinite(je) else "")
        fd = r.assembled_homologous_distance_bp
        cell(2, k, "yes" if np.isfinite(fd) and fd <= HOMOLOGOUS_BP else "partial" if r.assembled_contig != "." else "no",
             (f"{fd:,.0f} bp from\nhomologous pos." if np.isfinite(fd) else "switch not\nhomologous")
             if r.assembled_contig != "." else "")
        if r.contig_junction != ".":
            mh = join_text(r.microhomology_bp, getattr(r, "alignment_overlap_bp", np.nan), r.inserted_bp)
            mh = mh.replace("microhomology", "microhom.").replace("alignment overlap", "aln overlap")
            cell(3, k, "yes", f"{r.contig_junction.split(' (')[0]}\n{mh}")
        else:
            cell(3, k, "no", "")
        sc, sp = r.split_reads_carrier, r.split_reads_panel_max
        ok = np.isfinite(sc) and sc >= 2 and (not np.isfinite(sp) or sp == 0)
        cell(4, k, "yes" if ok else "partial" if np.isfinite(sc) and sc > 0 else "no",
             f"{sc:.0f} / {sp:.0f}" if np.isfinite(sc) and np.isfinite(sp) else "not counted")
        x = EVIDENCE_COLS[5][1]
        ax.add_patch(Rectangle((x, k - 0.3), 0.14, 0.6, facecolor=MECH_COLOR[r.mechanism], lw=0))
        ax.text(x + 0.22, k, r.mechanism, va="center", fontsize=style.BASE_FONT - 2, color=style.INK)
        ax.text(-0.12, k, r.label, ha="right", va="center", fontsize=style.BASE_FONT - 1.5,
                color=cohort_lib.COLOR[r.group], fontweight="bold")
    for name, x in EVIDENCE_COLS:
        ax.text(x, -0.75, name, ha="left", va="bottom", fontsize=style.BASE_FONT - 1.8, color=style.MUTED)
    ax.set_xlim(-0.05, 7.55)
    ax.set_ylim(n - 0.5, -0.5)
    ax.axis("off")
    ax.text(7.55, -1.35, "● supports   ○ partial   · absent", va="bottom", ha="right", fontsize=style.BASE_FONT - 2,
            color=style.MUTED)
    style.panel_label(fig, ax, "c", "Evidence for each deletion (Table 2)", dx=-0.1, dy=0.035)


def draw_fusions(fig, spec, t) -> None:
    blocks, ev = t["blocks"], t["evidence"]
    keys = blocks[["sample", "hap", "contig"]].drop_duplicates().values.tolist() if len(blocks) else []
    if not keys:
        ax = fig.add_subplot(spec)
        ax.text(0.5, 0.5, "no assembled contig across a deletion\n(scripts/duplicons 08 assembled_fusions.tsv)",
                ha="center", va="center", transform=ax.transAxes, color=style.MUTED, fontsize=style.BASE_FONT - 1)
        ax.axis("off")
        style.panel_label(fig, ax, "d", "Assembled junctions", dx=-0.05, dy=0.045)
        return
    keys = keys[:3]
    sub = GridSpecFromSubplotSpec(1, max(2, len(keys)), subplot_spec=spec, wspace=0.12)
    arcs = t["arcs"]
    for i, (s, hap, contig) in enumerate(keys):
        ax = fig.add_subplot(sub[i])
        b = blocks[(blocks["sample"] == s) & (blocks["hap"] == hap) & (blocks["contig"] == contig)]
        f = b.iloc[0]
        q0, q1 = b["q_start"].min(), b["q_end"].max()
        segs = {}
        for side, (x0, x1) in (("proximal", (0.0, 0.47)), ("distal", (0.53, 1.0))):
            bb = b[b["side"] == side]
            lo, hi = (bb["t_start"].min(), bb["t_end"].max()) if len(bb) else (np.nan, np.nan)
            pr = arcs[arcs["pair"] == f.fusion_pair] if len(arcs) and "pair" in arcs else pd.DataFrame()
            if len(pr):
                c0, c1 = (pr["a0"].iloc[0], pr["a1"].iloc[0]) if side == "proximal" else (pr["b0"].iloc[0], pr["b1"].iloc[0])
                lo, hi = np.nanmin([lo, c0]), np.nanmax([hi, c1])
            segs[side] = (lo, hi, x0, x1)

        def tx(pos, side):
            lo, hi, x0, x1 = segs[side]
            return x0 + (pos - lo) / max(1, hi - lo) * (x1 - x0)

        # a contig aligned mostly to the minus strand is drawn reverse-complemented, so that
        # its blocks run left to right like the reference
        flip = b.loc[b["strand"] == "-", "aln_len"].sum() > b.loc[b["strand"] == "+", "aln_len"].sum()

        def qx(q):
            v = (q - q0) / max(1, q1 - q0)
            return 1 - v if flip else v

        for r in b.itertuples(index=False):
            a, c = tx(r.t_start, r.side), tx(r.t_end, r.side)
            top0, top1 = sorted((qx(r.q_start), qx(r.q_end)))
            if (r.strand == "-") != flip:
                a, c = c, a
            col = "#6B6B6B" if r.side == "proximal" else "#A8A8A8"
            ax.add_patch(Polygon([(top0, 0.78), (top1, 0.78), (c, 0.22), (a, 0.22)], closed=True,
                                 facecolor=col, alpha=0.55 if r.identity >= 0.995 else 0.25, lw=0))
            ax.add_patch(Rectangle((top0, 0.78), top1 - top0, 0.1, facecolor=col, lw=0))
        for side in ("proximal", "distal"):
            lo, hi, x0, x1 = segs[side]
            if not np.isfinite(lo):
                continue
            ax.plot([x0, x1], [0.17, 0.17], color=style.INK, lw=1.2, solid_capstyle="butt")
            ax.text(x0, 0.05, f"{lo / 1e6:.2f}", ha="left", va="top", fontsize=style.BASE_FONT - 2.2, color=style.MUTED)
            ax.text(x1, 0.05, f"{hi / 1e6:.2f} Mb", ha="right", va="top", fontsize=style.BASE_FONT - 2.2,
                    color=style.MUTED)
            pr = arcs[arcs["pair"] == f.fusion_pair] if len(arcs) and "pair" in arcs else pd.DataFrame()
            if len(pr):
                c0, c1 = (pr["a0"].iloc[0], pr["a1"].iloc[0]) if side == "proximal" else (pr["b0"].iloc[0], pr["b1"].iloc[0])
                ax.add_patch(Rectangle((tx(c0, side), 0.12), tx(c1, side) - tx(c0, side), 0.1, facecolor=style.SD,
                                       lw=0, zorder=0))
            pos = f.ref_left if side == "proximal" else f.ref_right
            ax.plot([tx(pos, side)] * 2, [0.1, 0.27], color=MECH_COLOR["NAHR, confirmed"], lw=1.1)
        e = ev[ev["sample_id"] == s]
        label = e["label"].iloc[0] if len(e) else s
        fd = f.fusion_homolog_distance_bp
        txt = (f"switch {fd:,.0f} bp from homologous positions of {display_pair(f.fusion_pair)}" if np.isfinite(fd)
               else "switch not at homologous positions of a direct pair")
        rc = ", reverse complement" if flip else ""
        ax.text(0, 1.13, f"{label} · {hap} {contig} ({b['contig_len'].iloc[0] / 1e3:,.0f} kb{rc})", transform=ax.transAxes,
                fontsize=style.BASE_FONT - 1.5, color=cohort_lib.COLOR[e["group"].iloc[0]] if len(e) else style.INK,
                fontweight="bold", va="bottom")
        ax.text(0, 1.0, txt, transform=ax.transAxes, ha="left", va="bottom", fontsize=style.BASE_FONT - 2,
                color=style.MUTED)
        ax.set_xlim(-0.01, 1.01)
        ax.set_ylim(0, 1)
        ax.axis("off")
        if i == 0:
            style.panel_label(fig, ax, "d", "Assembled contigs across a deletion (top, contig; bottom, reference at "
                              "the proximal | distal edge)", dx=-0.1, dy=0.05)


def draw_junction(fig, spec, t) -> None:
    nrj, split, cn = t["junctions"], t["split"], t["edge_cn"]
    if nrj.empty:
        ax = fig.add_subplot(spec)
        ax.text(0.5, 0.5, "no contig junction to another locus at a deletion edge", ha="center", va="center",
                transform=ax.transAxes, color=style.MUTED)
        ax.axis("off")
        style.panel_label(fig, ax, "e", "Non-recurrent junctions", dx=-0.1, dy=0.03)
        return
    j = nrj.iloc[0]
    # Show the non-recurrent junction with the strongest carrier-specific
    # split-read support instead of relying on input-table row order.
    if len(split) and {"junction", "is_carrier", "split_reads"} <= set(split.columns):
        carrier_split = split[split["is_carrier"].astype(str).isin(["True", "true", "1"])].copy()
        carrier_split["split_reads"] = pd.to_numeric(carrier_split["split_reads"], errors="coerce")
        carrier_split = carrier_split.sort_values("split_reads", ascending=False, na_position="last")
        if len(carrier_split):
            supported = nrj[nrj["junction"] == carrier_split.iloc[0]["junction"]]
            if len(supported):
                j = supported.iloc[0]
    ev = t["evidence"]
    e = ev[ev["sample_id"] == j["sample_id"]].iloc[0]
    sub = GridSpecFromSubplotSpec(1, 3, subplot_spec=spec, width_ratios=[1.1, 1.0, 1.25], wspace=0.35)
    # copy number around the edge
    ax = fig.add_subplot(sub[0])
    x = cn[cn["junction"] == j["junction"]] if len(cn) else cn
    if len(x):
        x = x[x["usable"].astype(str).isin(["True", "true", "1"])]
        x = x.assign(mid=(x["start"] + x["end"]) / 2)
        other = x[~x["is_carrier"].astype(str).isin(["True", "true", "1"])].groupby("mid")["cn_norm"]
        med = other.median()
        q1, q3 = other.quantile(0.1), other.quantile(0.9)
        ax.fill_between(mb(med.index), q1, q3, color=style.LIGHT, alpha=0.5, lw=0, label="other genomes (10–90%)")
        ax.plot(mb(med.index), med, color=style.MUTED, lw=0.9, label="other genomes (median)")
        c = x[x["is_carrier"].astype(str).isin(["True", "true", "1"])].sort_values("mid")
        col = cohort_lib.COLOR[e["group"]]
        ax.scatter(mb(c["mid"]), c["cn_norm"], s=3, color=col, alpha=0.5, lw=0)
        ax.plot(mb(c["mid"]), c["cn_norm"].rolling(5, center=True, min_periods=2).median(), color=col, lw=1.1,
                label=e["label"])
        ax.legend(fontsize=style.BASE_FONT - 2.2, loc="upper left", frameon=False, ncol=1, borderaxespad=0.2)
    else:
        ax.text(0.5, 0.5, "edge outside the SUNK window", ha="center", va="center", transform=ax.transAxes,
                color=style.MUTED, fontsize=style.BASE_FONT - 1.5)
    ax.axvline(mb(j["chr15_pos"]), color=style.INK, lw=0.7, ls=(0, (2, 2)))
    ax.set_ylim(-0.1, 3.9)
    ax.set_xlabel("chr15 (Mb)")
    ax.set_ylabel("copy number (SUNK)", fontsize=style.BASE_FONT - 1)
    ax.tick_params(labelsize=style.BASE_FONT - 2)
    style.panel_label(fig, ax, "e", f"Non-recurrent junction of {e['label']} "
                      f"({j['edge']} edge)", dx=-0.1, dy=0.035)
    # junction schematic
    ax = fig.add_subplot(sub[1])
    col_a, col_b = style.INK, style.TRACE
    ax.add_patch(Rectangle((0.02, 0.28), 0.46, 0.1, facecolor=col_a, lw=0))
    ax.add_patch(Rectangle((0.52, 0.28), 0.46, 0.1, facecolor=col_b, lw=0))
    ax.text(0.25, 0.2, f"{CHROM}\n…{int(j['chr15_pos']):,}", ha="center", va="top", fontsize=style.BASE_FONT - 2,
            color=col_a)
    idt = num(j.get("partner_identity"))
    ax.text(0.75, 0.2, f"{j['partner_chrom']}\n{int(j['partner_pos']):,}…"
            + (f"\nidentity {idt:.3f}" if np.isfinite(idt) else ""), ha="center", va="top",
            fontsize=style.BASE_FONT - 2, color=col_b)
    ax.add_patch(Rectangle((0.12, 0.66), 0.38, 0.08, facecolor=col_a, lw=0))
    ax.add_patch(Rectangle((0.50, 0.66), 0.38, 0.08, facecolor=col_b, lw=0))
    ax.annotate("", xy=(0.9, 0.7), xytext=(0.88, 0.7), arrowprops=dict(arrowstyle="-|>", color=col_b, lw=0.8))
    ax.plot([0.5, 0.5], [0.4, 0.64], color=style.MUTED, lw=0.6, ls=(0, (2, 2)))
    joint = join_text(j.get("microhomology_bp"), j.get("alignment_overlap_bp"), j.get("inserted_bp"))
    nb = num(j.get("filtered_blocks_between"))
    if np.isfinite(nb) and nb > 0:
        joint += f"\n({nb:.0f} short/ambiguous block{'s' if nb > 1 else ''} between)"
    ax.text(0.5, 0.8, f"{j['assembly']} {j['contig']}", ha="center", va="bottom", fontsize=style.BASE_FONT - 2,
            color=style.INK)
    ax.text(0.5, 0.5, joint, ha="center", va="center", fontsize=style.BASE_FONT - 1.5, color=style.INK,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=style.LIGHT, lw=0.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("assembled junction", fontsize=style.BASE_FONT - 1, color=style.MUTED)
    # split reads in every genome
    ax = fig.add_subplot(sub[2])
    s = split[split["junction"] == j["junction"]] if len(split) else split
    if len(s):
        order = {g: i for i, g in enumerate(cohort_lib.GROUPS)}
        s = s.assign(o=s["group"].map(order)).sort_values(["o", "sample_id"]).reset_index(drop=True)
        carrier = s["is_carrier"].astype(str).isin(["True", "true", "1"])
        colors = [cohort_lib.COLOR[g] if c else style.LIGHT for g, c in zip(s["group"], carrier)]
        values = pd.to_numeric(s["split_reads"], errors="coerce").fillna(0)
        ax.bar(range(len(s)), values, color=colors, width=0.7)
        ax.set_xticks(range(len(s)))
        ax.set_xticklabels(s["label"], rotation=90, fontsize=style.BASE_FONT - 2.5)
        ax.set_ylabel(f"split reads {CHROM}–{j['partner_chrom']}", fontsize=style.BASE_FONT - 1)
        carrier_max = values.loc[carrier].max()
        others = values.loc[~carrier]
        ax.set_title(f"carrier {carrier_max:.0f}; other {len(others)} genomes, max {others.max():.0f}",
                     fontsize=style.BASE_FONT - 1.5, color=style.MUTED)
        ax.tick_params(axis="y", labelsize=style.BASE_FONT - 2)
        ax.set_xlim(-0.6, len(s) - 0.4)
        if values.max() <= 0:
            ax.scatter(range(len(s)), np.zeros(len(s)), s=9, color=colors, zorder=3, clip_on=False)
            ax.set_ylim(-0.08, 1.0)
            ax.text(0.5, 0.62, f"No supporting split reads\ncarrier 0; all {len(others)} other genomes 0",
                    transform=ax.transAxes, ha="center", va="center", fontsize=style.BASE_FONT - 1.5,
                    color=style.INK, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=style.LIGHT, lw=0.5))
    else:
        ax.text(0.5, 0.5, "split reads not counted\n(pysam or BAMs missing)", ha="center", va="center",
                transform=ax.transAxes, color=style.MUTED, fontsize=style.BASE_FONT - 1.5)
        ax.axis("off")


def render(outdir: Path, t: dict, meta: dict) -> list[Path]:
    style.setup()
    n = len(t["evidence"])
    fig = plt.figure(figsize=(style.WIDTH, 10.8 + 0.25 * max(0, n - 8)))
    gs = GridSpec(9, 1, figure=fig, height_ratios=[1.35 + 0.2 * n, 1.55, 1.05, 1.3, 0.45 + 0.26 * n, 1.0, 0.8, 0.95,
                                                   1.35], hspace=0.0, left=0.12, right=0.9, top=0.975, bottom=0.045)
    draw_map(fig, gs[0], t)
    draw_homologous(fig, gs[2], t)
    draw_evidence(fig, gs[4], t)
    draw_fusions(fig, gs[6], t)
    draw_junction(fig, gs[8], t)
    return style.save(fig, outdir, "Figure7")


# ------------------------------------------------------------------ Table 2
TABLE2 = [("label", "Patient"), ("sample_id", "Sample"), ("group", "Group"), ("deletion_class", "Class"),
          ("start", "Proximal edge"), ("start_source", "Proximal source"), ("end", "Distal edge"),
          ("end_source", "Distal source"), ("size_mb", "Size (Mb)"), ("sd_pair", "SD pair"),
          ("homologous_distance_bp", "Edges from homologous positions (bp)"),
          ("junction_reads", "Junction reads"), ("junction_reads_expected", "Expected"),
          ("junction_reads_p", "Junction-read p"), ("assembled_contig", "Assembled contig"),
          ("assembled_homologous_distance_bp", "Contig switch from homologous positions (bp)"),
          ("contig_junction", "Contig junction to other locus"), ("microhomology_bp", "Microhomology (bp)"),
          ("alignment_overlap_bp", "Alignment overlap at junction (bp)"),
          ("split_reads_carrier", "Split reads (carrier)"), ("split_reads_panel_max", "Split reads (others, max)"),
          ("duplicons08_status", "Status (duplicons 08)"), ("mechanism", "Mechanism")]


def table2(ev: pd.DataFrame) -> pd.DataFrame:
    return ev[[c for c, _ in TABLE2 if c in ev]].rename(columns=dict(TABLE2))


def latex(df: pd.DataFrame) -> str:
    keep = ["Patient", "Group", "Class", "Proximal edge", "Distal edge", "Size (Mb)", "SD pair",
            "Edges from homologous positions (bp)", "Junction reads", "Contig switch from homologous positions (bp)",
            "Contig junction to other locus", "Split reads (carrier)", "Mechanism"]
    d = df[[c for c in keep if c in df]].copy()

    def esc(v) -> str:
        if isinstance(v, float):
            if not np.isfinite(v):
                return "--"
            v = f"{v:,.0f}" if v >= 100 or float(v).is_integer() else f"{v:.2f}"
        s = str(v)
        if s in (".", "nan", "NA"):
            return "--"
        for a, b in (("\\", "\\textbackslash{}"), ("_", "\\_"), ("%", "\\%"), ("&", "\\&"), ("#", "\\#"),
                     ("–", "--"), ("≥", "$\\geq$")):
            s = s.replace(a, b)
        return s

    head = " & ".join(esc(c) for c in d.columns)
    body = "\n".join(" & ".join(esc(v) for v in row) + r" \\" for row in d.itertuples(index=False))
    return ("% Table 2 -- generated by scripts/figures/FIGURE_7.py\n\\begin{table}[ht]\n\\centering\\scriptsize\n"
            f"\\begin{{tabular}}{{{'l' * len(d.columns)}}}\n\\toprule\n{head} \\\\\n\\midrule\n{body}\n"
            "\\bottomrule\n\\end{tabular}\n\\caption{Breakpoint evidence and mechanism of each deletion.}\n"
            "\\label{tab:mechanisms}\n\\end{table}\n")


def report(outdir, t, meta, paths_out) -> Path:
    r = style.Report("Figure 7 report: deletion mechanisms")
    r.p("Generated by `scripts/figures/FIGURE_7.py` from scripts/duplicons 01/03/05/08 and the phased BAMs.")
    for note in meta.get("notes", []):
        r.p(f"- {note}")
    ev = t["evidence"]
    r.h("Table 2: evidence per deletion (panel c)")
    r.table(table2(ev))
    r.h("Mechanisms")
    r.table(ev.groupby(["mechanism", "group"]).size().rename("carriers").reset_index())
    if len(t["hom_meta"]):
        r.h("Homologous-position test (panel b)")
        r.p("Proximal SUNK edge in copy A and distal SUNK edge in copy B projected onto copy A by linear homology; "
            f"scripts/duplicons calls them homologous within {HOMOLOGOUS_BP // 1000} kb (or the combined SUNK support).")
        r.table(t["hom_meta"])
    if len(t["blocks"]):
        r.h("Assembled contigs across a deletion (panel d)")
        r.table(t["blocks"][["sample", "hap", "contig", "q_start", "q_end", "strand", "t_start", "t_end", "identity",
                             "mapq", "side"]])
    if len(t["junctions"]):
        r.h("Non-recurrent junctions (panel e)")
        r.table(t["junctions"])
        if len(t["split"]):
            r.table(t["split"])
    r.h("Caption draft")
    counts = ev["mechanism"].value_counts()
    summary = "; ".join(f"{m}: {counts[m]}" for m in MECHANISMS if m in counts)
    r.p(f"**Mechanisms of the 15q11–q13 deletions.** **a,** Directly oriented segmental-duplication pairs of the "
        f"region (arcs; black, pairs with one copy at each edge of a deletion) and the deletion interval of each "
        f"carrier coloured by mechanism ({summary}). **b,** Paralog-specific (SUNK) copy number along copy A "
        f"(proximal) and along copy B (distal, projected onto copy A by homology) of the pair at each carrier's "
        f"edges; after non-allelic homologous recombination both transitions fall at the same homologous position. "
        f"**c,** Evidence for each deletion (Table 2). **d,** Alignment blocks of assembled contigs that cross a "
        f"deletion, from the contig (top) to the proximal and distal edges of the reference (bottom); tan, the two "
        f"copies of the duplication pair; green, the switch between them. **e,** The non-recurrent junction with the "
        f"strongest carrier-specific split-read support: copy number at the edge, the assembled junction to the partner "
        f"locus and split reads in every genome of the cohort.")
    r.p("")
    r.p("Files: " + ", ".join(p.name for p in paths_out))
    return r.write(outdir, "Figure7_report.md")


NAMES = {"evidence": "Figure7c_evidence.tsv", "arcs": "Figure7a_sd_pairs.tsv", "sd_blocks": "Figure7a_sd_blocks.tsv",
         "hom": "Figure7b_homologous_copy_number.tsv.gz", "hom_meta": "Figure7b_homologous_edges.tsv",
         "blocks": "Figure7d_contig_blocks.tsv", "junctions": "Figure7e_junctions.tsv",
         "split": "Figure7e_split_reads.tsv", "edge_cn": "Figure7e_copy_number.tsv.gz"}


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    paths_lib.add_common_arguments(ap)
    ap.add_argument("--render-only", action="store_true", help="redraw from the tables of a previous run")
    ap.add_argument("--no-split-reads", action="store_true", help="do not count split reads in the BAMs")
    a = ap.parse_args(argv)
    paths = paths_lib.resolve(a)
    outdir = paths.figure_dir(7, a.outdir)
    if a.render_only:
        t = {}
        for k, v in NAMES.items():
            p = outdir / "tables" / v
            t[k] = pd.read_csv(p, sep="\t", low_memory=False) if p.is_file() and p.stat().st_size > 1 else pd.DataFrame()
        meta = json.loads((outdir / "tables" / "Figure7_meta.json").read_text())
    else:
        cohort = cohort_lib.load(paths.metadata)
        t, meta = build(paths, cohort, count_split=not a.no_split_reads)
        for k, v in NAMES.items():
            style.write_table(t[k], outdir, v)
        (outdir / "tables" / "Figure7_meta.json").write_text(json.dumps(meta, indent=1))
        t2 = table2(t["evidence"])
        style.write_table(t2, outdir, "Table2_mechanisms.tsv")
        (outdir / "tables" / "Table2_mechanisms.tex").write_text(latex(t2))
    if "used" not in t["arcs"]:
        t["arcs"] = pd.DataFrame(columns=["pair", "a0", "a1", "b0", "b1", "identity", "used", "length"])
    t["arcs"]["used"] = t["arcs"]["used"].astype(str).isin(["True", "true", "1"])
    paths_out = render(outdir, t, meta)
    print(report(outdir, t, meta, paths_out))


if __name__ == "__main__":
    main()
