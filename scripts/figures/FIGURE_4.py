#!/usr/bin/env python3
"""
Figure 4 -- Single-molecule methylation at SNORD116.

Computed from the HiPhase BAMs (results/04_phasing, MM/ML + HP/PS tags):
  molecules  primary reads with MAPQ >= 20 over three regions of analysis 00's
             prespecified_regions.tsv: the imprinting centre, the SNORD116 cluster and the
             OCA2 non-imprinted control. Per molecule and region: CpGs, mean methylation
             probability, methylated fraction (probability >= 0.5), entropy of that fraction
             (bits; 0 = all-or-none, 1 = half methylated) and the switch rate between
             adjacent CpGs.
  parents    deletion carriers: every molecule is the retained copy (PWS-DEL maternal,
             AS-DEL paternal). Controls: HP haplotypes are oriented by Figure 1's assignment
             (or IC methylation); a SNORD116 molecule inherits that orientation only if it is
             in the same phase block (PS) as the imprinting centre. Genomes without such a
             block are oriented by SNORD116 methylation itself (the haplotype closer to the
             AS-DEL molecules is paternal-like) and flagged.
Panels  a,b  paternal-like and maternal-like control molecules across SNORD116 (heatmap of
             representative molecules, mean profile of all, SNORD116 copies)
        c    per-molecule entropy by region and group; tests on per-genome medians
        d    SNORD116 molecule methylation, AS-DEL vs PWS-DEL; per-genome test (exact label
             permutation) and molecule-level effect sizes (descriptive)

Usage: python3 scripts/figures/FIGURE_4.py [--results DIR] [--include-22q] [--render-only]
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
from figlib import modbam  # noqa: E402
from figlib import paths as paths_lib  # noqa: E402
from figlib import stats  # noqa: E402
from figlib import style  # noqa: E402
from figlib.style import mb  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

CHROM = ann.CHROM
REGION_IDS = {"IC": "PWS/AS imprinting centre", "SNORD116": "SNORD116", "control": "OCA2 downstream control"}
REGION_LABEL = {"IC": "Imprinting centre", "SNORD116": "SNORD116 cluster", "control": "OCA2 control"}
FALLBACK_REGIONS = {"IC": ann.IC}
MIN_MAPQ = 20
MIN_CPGS = {"IC": 10, "SNORD116": 20, "control": 20}
MIN_SPAN_D = 10_000            # panel d: molecules spanning >= 10 kb of the cluster
DISPLAY_ROWS = 18
MAX_COLUMNS = 240
UNMETH, METH, MISSING = "#3A74B7", "#D6544A", "#F0F0F0"


# ------------------------------------------------------------------ regions and orientation
def regions(paths) -> dict[str, tuple[int, int]]:
    pre = ann.prespecified_regions(paths)
    out = {}
    for key, rid in REGION_IDS.items():
        row = pre[pre["region_id"] == rid] if len(pre) else pre
        if len(row):
            out[key] = (int(row["start"].iloc[0]), int(row["end"].iloc[0]))
        elif key in FALLBACK_REGIONS:
            out[key] = FALLBACK_REGIONS[key]
    if "SNORD116" not in out:
        g = ann.collapse_families(ann.read_gtf_genes(paths.gtf, CHROM, 22_000_000, 24_000_000))
        g = g[g["gene"] == "SNORD116 cluster"]
        if len(g):
            out["SNORD116"] = (int(g["start"].iloc[0]) - 5_000, int(g["end"].iloc[0]) + 5_000)
    missing = [k for k in REGION_IDS if k not in out]
    if missing:
        raise SystemExit(f"regions not found: {missing} (run scripts/analysis/00_build_prespecified_regions.py)")
    return out


def hap_orientation(paths, sample: str, cis) -> tuple[dict[int, str], str]:
    """{HP: 'maternal'|'paternal'} from Figure 1, else IC methylation of the hap tracks."""
    table = paths.figures / "figure_1" / "tables" / "Figure1_parental_like_assignment.tsv"
    if table.is_file():
        t = pd.read_csv(table, sep="\t")
        t = t[(t["sample_id"] == sample) & t["haplotype_label"].isin(["hap1", "hap2"])]
        mp = {int(h[-1]): a for h, a in zip(t["haplotype_label"], t["parental_assignment"])
              if a in ("maternal", "paternal")}
        if len(mp) == 2 and set(mp.values()) == {"maternal", "paternal"}:
            return mp, "Figure 1 assignment"
    means = {}
    for hp in (1, 2):
        path = cis.find_track(paths.methylation, sample, f"hap{hp}")
        tr = cis.read_track(path, CHROM, *ann.IC) if path else None
        means[hp] = float(np.average(tr.beta, weights=tr.coverage)) if tr is not None and len(tr.beta) else np.nan
    if np.all(np.isfinite(list(means.values()))):
        mat = max(means, key=means.get)
        return {mat: "maternal", 3 - mat: "paternal"}, "IC methylation of hap tracks"
    return {}, "none"


# ------------------------------------------------------------------ extraction
def extract(paths, cohort, regs, include_22q: bool, cis):
    mol_rows, calls = [], []
    notes = []
    for s in cohort.samples:
        g = cohort.group(s)
        bam, phased = modbam.find_bam(paths, s)
        if bam is None:
            notes.append(f"{s}: no indexed BAM in 04_phasing or 01_alignment")
            continue
        if not phased:
            notes.append(f"{s}: {bam.name} has no HiPhase tags; haplotypes unavailable")
        per_region = {}
        for key, (a, b) in regs.items():
            mols = modbam.molecules(bam, s, CHROM, a, b, MIN_MAPQ)
            per_region[key] = mols
            for m in mols:
                pos, pr = m.within(a, b)
                if len(pos) < MIN_CPGS[key]:
                    continue
                frac = float(np.mean(pr >= 0.5))
                mol_rows.append({"sample_id": s, "group": g, "region": key, "read_id": m.read_id, "hp": m.hp,
                                 "ps": m.ps, "start": int(pos[0]), "end": int(pos[-1]) + 1, "n_cpg": len(pos),
                                 "mean_methylation": float(np.mean(pr)), "methylated_fraction": frac,
                                 "entropy_bits": float(modbam.binary_entropy(np.array([frac]))[0]),
                                 "switch_rate": modbam.switch_rate(pr)})
                if key == "SNORD116" and g in (("Control", "DiGeorge") if include_22q else ("Control",)):
                    calls.append(pd.DataFrame({"read_id": m.read_id, "sample_id": s, "pos": pos, "prob": pr}))
        # IC phase block per HP (for controls)
        ic_ps = {}
        for hp in (1, 2):
            ps = [m.ps for m in per_region.get("IC", []) if m.hp == hp and m.ps >= 0]
            if ps:
                ic_ps[hp] = max(set(ps), key=ps.count)
        mol_rows.append({"sample_id": s, "group": "_meta", "region": "_ic_ps", "read_id": json.dumps(ic_ps)})
    mol = pd.DataFrame(mol_rows)
    meta = mol[mol["group"] == "_meta"].set_index("sample_id")["read_id"].map(json.loads).to_dict()
    mol = mol[mol["group"] != "_meta"].reset_index(drop=True)
    mol["hp"] = mol["hp"].astype(int)
    mol["ps"] = mol["ps"].astype(int)
    return mol, (pd.concat(calls, ignore_index=True) if calls else pd.DataFrame(columns=["read_id", "sample_id", "pos", "prob"])), meta, notes


def assign_parents(mol: pd.DataFrame, ic_ps: dict, paths, cohort, cis) -> tuple[pd.DataFrame, pd.DataFrame]:
    """parent_label and assignment_source for every molecule."""
    mol = mol.copy()
    mol["parent_label"] = "unassigned"
    mol["assignment_source"] = "."
    mol.loc[mol["group"] == "PWS-DEL", ["parent_label", "assignment_source"]] = ["maternal", "retained copy (PWS-DEL)"]
    mol.loc[mol["group"] == "AS-DEL", ["parent_label", "assignment_source"]] = ["paternal", "retained copy (AS-DEL)"]
    mol.loc[mol["group"] == "PWS-mUPD", ["parent_label", "assignment_source"]] = ["maternal", "two maternal copies (mUPD)"]
    sn = mol[mol["region"] == "SNORD116"]
    as_med = sn.loc[sn["group"] == "AS-DEL", "mean_methylation"].median()
    pws_med = sn.loc[sn["group"] == "PWS-DEL", "mean_methylation"].median()
    paternal_higher = bool(as_med > pws_med) if np.isfinite(as_med) and np.isfinite(pws_med) else True
    deletion_refs_available = np.isfinite(as_med) and np.isfinite(pws_med)
    fallback_boundary = float((as_med + pws_med) / 2) if deletion_refs_available else np.nan
    summary = []
    for s in cohort.of("Control", "DiGeorge"):
        orient, source = hap_orientation(paths, s, cis)
        blocks = {int(k): v for k, v in ic_ps.get(s, {}).items()}
        for key in mol.loc[mol["sample_id"] == s, "region"].unique():
            region_idx = (mol["sample_id"] == s) & (mol["region"] == key)
            phased_idx = region_idx & (mol["hp"] > 0)
            same = phased_idx & mol.apply(lambda r: blocks.get(int(r["hp"])) == int(r["ps"]), axis=1)
            if orient:
                mol.loc[same, "parent_label"] = mol.loc[same, "hp"].map(orient)
                mol.loc[same, "assignment_source"] = f"HP in the IC phase block ({source})"
            rest = phased_idx & ~same
            if key == "SNORD116" and rest.any():
                # different phase block: orient the two HP groups of each block by methylation
                for ps, d in mol[rest].groupby("ps"):
                    means = d.groupby("hp")["mean_methylation"].median()
                    if len(means) == 2:
                        hi = int(means.idxmax())
                        lab = {hi: "paternal" if paternal_higher else "maternal",
                               3 - hi: "maternal" if paternal_higher else "paternal"}
                        sel = rest & (mol["ps"] == ps)
                        mol.loc[sel, "parent_label"] = mol.loc[sel, "hp"].map(lab)
                        mol.loc[sel, "assignment_source"] = "SNORD116 methylation (phase block without the IC)"
            # Some BAMs contain MM/ML calls but no HP/PS tags. The old code
            # silently left every such control molecule unassigned, making
            # panels a and b empty. Apply the documented SNORD116 methylation
            # fallback at molecule level and retain explicit provenance.
            unresolved = region_idx & (mol["parent_label"] == "unassigned")
            if key == "SNORD116" and deletion_refs_available and unresolved.any():
                values = mol.loc[unresolved, "mean_methylation"]
                high_label = "paternal" if paternal_higher else "maternal"
                low_label = "maternal" if paternal_higher else "paternal"
                mol.loc[unresolved, "parent_label"] = np.where(values >= fallback_boundary, high_label, low_label)
                mol.loc[unresolved, "assignment_source"] = "SNORD116 molecule methylation (unphased fallback)"
            summary.append({"sample_id": s, "region": key,
                            "molecules": int(region_idx.sum()),
                            "hp_in_ic_block": int(same.sum()),
                            "other_block_or_unphased": int((region_idx & ~same).sum()),
                            "orientation_source": source})
    return mol, pd.DataFrame(summary)


# ------------------------------------------------------------------ statistics
def entropy_stats(mol: pd.DataFrame, groups: list[str]) -> pd.DataFrame:
    rows = []
    for key in REGION_IDS:
        d = mol[mol["region"] == key]
        per = d.groupby(["sample_id", "group"])["entropy_bits"].median().reset_index()
        h, p, eta = stats.kruskal([per.loc[per["group"] == g, "entropy_bits"] for g in groups])
        pp = stats.kruskal_permutation([per.loc[per["group"] == g, "entropy_bits"] for g in groups], n_mc=5_000)
        diff, pl, nperm = stats.label_permutation(per.loc[per["group"] == "AS-DEL", "entropy_bits"],
                                                  per.loc[per["group"] == "PWS-DEL", "entropy_bits"])
        hm, pm, etam = stats.kruskal([d.loc[d["group"] == g, "entropy_bits"] for g in groups])
        rows.append({"region": key, "genomes": len(per), "molecules": len(d), "kruskal_H_genomes": h,
                     "kruskal_p_genomes": p, "permutation_p_genomes": pp, "eta2_genomes": eta,
                     "as_minus_pws_genome_medians": diff, "as_vs_pws_label_permutation_p": pl, "labellings": nperm,
                     "kruskal_H_molecules": hm, "kruskal_p_molecules_descriptive": pm, "eta2_molecules": etam})
    t = pd.DataFrame(rows)
    t["permutation_q_genomes"] = stats.bh(t["permutation_p_genomes"])
    return t


def snord116_stats(mol: pd.DataFrame) -> dict:
    d = mol[(mol["region"] == "SNORD116") & mol["group"].isin(["AS-DEL", "PWS-DEL"])
            & (mol["end"] - mol["start"] >= MIN_SPAN_D)]
    a = d.loc[d["group"] == "AS-DEL", "mean_methylation"].to_numpy()
    b = d.loc[d["group"] == "PWS-DEL", "mean_methylation"].to_numpy()
    per = d.groupby(["sample_id", "group"])["mean_methylation"].median().reset_index()
    diff, p_perm, nperm = stats.label_permutation(per.loc[per["group"] == "AS-DEL", "mean_methylation"],
                                                  per.loc[per["group"] == "PWS-DEL", "mean_methylation"])
    u, p_mw = stats.mann_whitney(a, b)
    lo, hi = stats.bootstrap_diff_ci(a, b)
    return {"molecules_AS": len(a), "molecules_PWS": len(b), "genomes_AS": int((per["group"] == "AS-DEL").sum()),
            "genomes_PWS": int((per["group"] == "PWS-DEL").sum()),
            "median_AS": float(np.median(a)) if len(a) else np.nan, "median_PWS": float(np.median(b)) if len(b) else np.nan,
            "delta_median_molecules": float(np.median(a) - np.median(b)) if len(a) and len(b) else np.nan,
            "delta_ci_low": lo, "delta_ci_high": hi, "mann_whitney_U": u, "mann_whitney_p_descriptive": p_mw,
            "rank_biserial": stats.rank_biserial(u, len(a), len(b)),
            "delta_mean_of_genome_medians": diff, "genome_label_permutation_p": p_perm, "labellings": nperm}


# ------------------------------------------------------------------ figure
def heatmap_panel(fig, spec, calls: pd.DataFrame, mol: pd.DataFrame, parent: str, window, genes, letter, title):
    sub = GridSpecFromSubplotSpec(3, 1, subplot_spec=spec, height_ratios=[0.32, 1, 0.14], hspace=0.06)
    ax_p = fig.add_subplot(sub[0])
    ax_h = fig.add_subplot(sub[1], sharex=ax_p)
    ax_g = fig.add_subplot(sub[2], sharex=ax_p)
    m = mol[(mol["region"] == "SNORD116") & (mol["parent_label"] == parent) & mol["group"].isin(["Control", "DiGeorge"])]
    keys = m[["sample_id", "read_id"]].drop_duplicates()
    c = calls.merge(keys, on=["sample_id", "read_id"], how="inner")
    edges = np.linspace(window[0], window[1], MAX_COLUMNS + 1)
    m = m.merge(c[["sample_id", "read_id"]].drop_duplicates(), on=["sample_id", "read_id"], how="inner")
    if len(c):
        c = c.assign(col=np.clip(np.searchsorted(edges, c["pos"], side="right") - 1, 0, MAX_COLUMNS - 1))
        mat = c.pivot_table(index=["sample_id", "read_id"], columns="col", values="prob", aggfunc="mean").reindex(columns=range(MAX_COLUMNS))
        order = m.set_index(["sample_id", "read_id"]).loc[mat.index].sort_values("mean_methylation", ascending=False)
        pick = order.index[np.unique(np.linspace(0, len(order) - 1, min(DISPLAY_ROWS, len(order))).round().astype(int))]
        shown = mat.loc[pick].to_numpy()
        img = np.where(np.isnan(shown), 0, np.where(shown >= 0.5, 2, 1))
        ax_h.imshow(img, aspect="auto", interpolation="nearest", cmap=ListedColormap([MISSING, UNMETH, METH]),
                    vmin=0, vmax=2, extent=(mb(window[0]), mb(window[1]), len(pick), 0))
        prof = mat.mean(axis=0).rolling(5, center=True, min_periods=1).mean()
        centres = (edges[:-1] + edges[1:]) / 2
        ax_p.plot(mb(centres), prof.to_numpy(), color=style.MATERNAL if parent == "maternal" else style.PATERNAL, lw=1.1)
        n_genomes = m["sample_id"].nunique()
        ax_p.text(0.0, 1.02, f"heatmap: {len(pick)} representative rows; line: all {len(m)} molecules "
                  f"({n_genomes} genomes)", transform=ax_p.transAxes,
                  fontsize=style.BASE_FONT - 1.5, color=style.MUTED, va="bottom")
    else:
        ax_h.text(0.5, 0.5, "no molecules", transform=ax_h.transAxes, ha="center", color=style.MUTED)
    ax_p.set_ylim(0, 1)
    ax_p.set_yticks([0, 1])
    ax_p.set_ylabel("mean", fontsize=style.BASE_FONT - 1)
    ax_p.tick_params(labelbottom=False)
    ax_h.set_yticks([])
    ax_h.set_ylabel("representative\nmolecules", fontsize=style.BASE_FONT - 1)
    ax_h.tick_params(labelbottom=False)
    for r in genes.itertuples(index=False):
        ax_g.add_patch(ann.matplotlib_rect((mb(r.start), 0.2), max(mb(r.end) - mb(r.start), 0.0002), 0.6, "#8E8E8E"))
    ax_g.set_ylim(0, 1)
    ax_g.set_yticks([])
    ax_g.spines["left"].set_visible(False)
    ax_g.set_xlim(mb(window[0]), mb(window[1]))
    ax_g.ticklabel_format(axis="x", useOffset=False, style="plain")
    ax_g.xaxis.set_major_locator(plt.MaxNLocator(4))
    ax_g.set_xlabel(f"chr15 (Mb); SNORD116 copies (n = {len(genes)})", fontsize=style.BASE_FONT - 0.5)
    style.panel_label(fig, ax_p, letter, title, dx=-0.05, dy=0.03)
    return ax_h


def render(outdir: Path, mol: pd.DataFrame, calls: pd.DataFrame, ent: pd.DataFrame, sn: dict, meta: dict) -> list[Path]:
    style.setup()
    window = tuple(meta["snord116_window"])
    genes = pd.DataFrame(meta["snord116_copies"], columns=["start", "end"])
    fig = plt.figure(figsize=(style.WIDTH, 7.4))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.15, 1], hspace=0.55, wspace=0.28,
                  left=0.08, right=0.98, top=0.95, bottom=0.11)
    heatmap_panel(fig, gs[0, 0], calls, mol, "paternal", window, genes, "a", "Paternal-like control molecules")
    heatmap_panel(fig, gs[0, 1], calls, mol, "maternal", window, genes, "b", "Maternal-like control molecules")
    fig.legend(handles=[Line2D([], [], marker="s", ls="", ms=5, color=METH, label="methylated CpG"),
                        Line2D([], [], marker="s", ls="", ms=5, color=UNMETH, label="unmethylated"),
                        Line2D([], [], marker="s", ls="", ms=5, color=MISSING, mec=style.LIGHT, label="not observed")],
               loc="upper center", ncol=3, fontsize=style.BASE_FONT - 1, bbox_to_anchor=(0.53, 0.5))

    # c -- entropy
    ax = fig.add_subplot(gs[1, 0])
    groups = meta["groups"]
    rng = np.random.default_rng(3)
    width = 0.8 / len(groups)
    for i, key in enumerate(REGION_IDS):
        if key == "SNORD116":
            ax.axvspan(i - 0.45, i + 0.45, color="#F4F4F4", lw=0, zorder=0)
        for j, g in enumerate(groups):
            x = i - 0.4 + width * (j + 0.5)
            v = mol.loc[(mol["region"] == key) & (mol["group"] == g), "entropy_bits"].to_numpy()
            if len(v):
                q = np.percentile(v, [25, 50, 75])
                ax.plot([x, x], [q[0], q[2]], color=cohort_lib.COLOR[g], lw=3.2, solid_capstyle="butt", alpha=0.45)
                per = mol[(mol["region"] == key) & (mol["group"] == g)].groupby("sample_id")["entropy_bits"].median()
                ax.scatter(x + rng.uniform(-width * 0.2, width * 0.2, len(per)), per, s=11,
                           marker=cohort_lib.MARKER[g], color=cohort_lib.COLOR[g], edgecolor="white", lw=0.4, zorder=3)
        r = ent[ent["region"] == key]
        if len(r):
            ax.text(i, 1.02, f"q = {style.fmt_p(r['permutation_q_genomes'].iloc[0])}", ha="center",
                    fontsize=style.BASE_FONT - 1.5, color=style.MUTED, transform=ax.get_xaxis_transform())
    ax.set_xticks(range(len(REGION_IDS)))
    ax.set_xticklabels([REGION_LABEL[k] for k in REGION_IDS])
    ax.set_xlim(-0.5, len(REGION_IDS) - 0.5)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("per-molecule entropy (bits)")
    ax.legend(handles=style.group_legend_handles(groups), loc="upper center", fontsize=style.BASE_FONT - 1.5,
              ncol=len(groups), handletextpad=0.1, columnspacing=0.6, bbox_to_anchor=(0.5, -0.1))
    style.panel_label(fig, ax, "c", "Per-molecule entropy by region", dx=-0.07, dy=0.035)

    # d -- SNORD116 AS vs PWS
    ax = fig.add_subplot(gs[1, 1])
    d = mol[(mol["region"] == "SNORD116") & mol["group"].isin(["AS-DEL", "PWS-DEL"]) & (mol["end"] - mol["start"] >= MIN_SPAN_D)]
    for k, g in enumerate(["PWS-DEL", "AS-DEL"]):
        v = d.loc[d["group"] == g, "mean_methylation"].to_numpy()
        if len(v) > 1:
            parts = ax.violinplot(v, positions=[k], widths=0.7, showextrema=False)
            for body in parts["bodies"]:
                body.set_facecolor(cohort_lib.COLOR[g])
                body.set_alpha(0.3)
                body.set_edgecolor("none")
        per = d[d["group"] == g].groupby("sample_id")["mean_methylation"].median()
        ax.scatter(np.full(len(per), k) + rng.uniform(-0.08, 0.08, len(per)), per, s=22, marker=cohort_lib.MARKER[g],
                   color=cohort_lib.COLOR[g], edgecolor="white", lw=0.5, zorder=3)
        if len(v):
            ax.hlines(np.median(v), k - 0.25, k + 0.25, color=style.INK, lw=1.2, zorder=4)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"PWS-DEL\nretained maternal\n{int(sn['molecules_PWS'])} molecules",
                        f"AS-DEL\nretained paternal\n{int(sn['molecules_AS'])} molecules"])
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("molecule mean methylation")
    ax.text(0.5, 1.0, f"Δ median {sn['delta_median_molecules']:+.2f} [{style.fmt(sn['delta_ci_low'], 2, True)}, "
            f"{style.fmt(sn['delta_ci_high'], 2, True)}]; genomes p = {style.fmt_p(sn['genome_label_permutation_p'])}",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=style.BASE_FONT - 1.5, color=style.INK)
    style.panel_label(fig, ax, "d", "SNORD116: retained paternal vs maternal copy", dx=-0.07, dy=0.035)
    return style.save(fig, outdir, "Figure4")


def report(outdir, mol, ent, sn, assign, meta, notes, paths_out) -> Path:
    r = style.Report("Figure 4 report: single-molecule methylation at SNORD116")
    r.p("Generated by `scripts/figures/FIGURE_4.py` from the phased BAMs. Entropy = Shannon entropy (bits) of each "
        "molecule's methylated fraction (CpGs with probability >= 0.5).")
    r.h("Regions")
    r.table(pd.DataFrame([{"region": REGION_LABEL[k], "start": v[0], "end": v[1], "min CpGs per molecule": MIN_CPGS[k]}
                          for k, v in meta["regions"].items()]))
    r.h("Molecules")
    r.table(mol.groupby(["region", "group"]).agg(molecules=("read_id", "size"), genomes=("sample_id", "nunique"),
                                                 median_cpgs=("n_cpg", "median"),
                                                 median_entropy=("entropy_bits", "median"),
                                                 median_methylation=("mean_methylation", "median")).reset_index())
    r.h("Parental assignment of control molecules")
    r.table(assign)
    lab = mol[(mol["region"] == "SNORD116") & mol["group"].isin(["Control", "DiGeorge"])]
    r.table(lab.groupby(["assignment_source", "parent_label"]).size().rename("molecules").reset_index())
    r.h("Panel c: entropy (tests on genome medians; molecule tests descriptive)")
    r.table(ent)
    r.h("Panel d: SNORD116 methylation, AS-DEL vs PWS-DEL")
    r.table(pd.DataFrame([sn]).T.reset_index().rename(columns={"index": "quantity", 0: "value"}))
    if notes:
        r.h("Notes")
        for n in notes:
            r.p(f"- {n}")
    r.h("Caption draft")
    r.p(f"**Single-molecule methylation at SNORD116.** **a,b,** Paternal-like and maternal-like control HiFi molecules "
        f"across the SNORD116 cluster (heatmap, {DISPLAY_ROWS} molecules sampled evenly across the molecule-methylation "
        f"range for legibility; columns, CpG bins; red methylated, blue unmethylated, grey not observed); top, mean "
        f"methylation of all assigned molecules (the row subsample is visual only); bottom, SNORD116 copies. Haplotypes were oriented "
        f"by imprinting-centre methylation within the same phase block; molecules without HP/PS tags were classified "
        f"relative to the retained-copy AS-DEL and PWS-DEL methylation distributions and explicitly flagged. **c,** "
        f"Per-molecule methylation entropy at the "
        f"imprinting centre, SNORD116 and a non-imprinted control (OCA2); symbols, genome medians; q, Kruskal–Wallis "
        f"permutation test on genome medians, Benjamini–Hochberg across regions. **d,** SNORD116 molecule methylation "
        f"in AS-DEL (retained paternal) and PWS-DEL (retained maternal) genomes; Δ median {sn['delta_median_molecules']:+.2f} "
        f"(molecule bootstrap 95% CI {style.fmt(sn['delta_ci_low'], 2, True)} to {style.fmt(sn['delta_ci_high'], 2, True)}); "
        f"exact label permutation on genome medians p = {style.fmt_p(sn['genome_label_permutation_p'])} "
        f"({sn['genomes_AS']} vs {sn['genomes_PWS']} genomes).")
    r.p("")
    r.p("Files: " + ", ".join(p.name for p in paths_out))
    return r.write(outdir, "Figure4_report.md")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    paths_lib.add_common_arguments(ap)
    ap.add_argument("--include-22q", action="store_true", help="panels a,b: also 22q11.2DS molecules")
    ap.add_argument("--render-only", action="store_true", help="redraw from the tables of a previous run")
    a = ap.parse_args(argv)
    paths = paths_lib.resolve(a)
    outdir = paths.figure_dir(4, a.outdir)
    if a.render_only:
        mol = style.read_table(outdir, "Figure4_molecules.tsv.gz")
        calls = style.read_table(outdir, "Figure4_snord116_control_cpg_calls.tsv.gz")
        ent = style.read_table(outdir, "Figure4c_entropy_statistics.tsv")
        sn = style.read_table(outdir, "Figure4d_snord116_statistics.tsv").iloc[0].to_dict()
        assign = style.read_table(outdir, "Figure4_control_assignment.tsv")
        meta = json.loads((outdir / "tables" / "Figure4_meta.json").read_text())
        paths_out = render(outdir, mol, calls, ent, sn, meta)
        print(report(outdir, mol, ent, sn, assign, meta, meta.get("notes", []), paths_out))
        return
    cis = paths_lib.import_analysis_package()
    cohort = cohort_lib.load(paths.metadata)
    regs = regions(paths)
    mol, calls, ic_ps, notes = extract(paths, cohort, regs, a.include_22q, cis)
    if mol.empty:
        raise SystemExit("no molecules with methylation tags in the regions: check results/04_phasing")
    mol, assign = assign_parents(mol, ic_ps, paths, cohort, cis)
    groups = [g for g in cohort_lib.GROUPS if (mol["group"] == g).any()]
    ent = entropy_stats(mol, groups)
    sn = snord116_stats(mol)
    genes = ann.read_gtf_genes(paths.gtf, CHROM, *regs["SNORD116"])
    genes = genes[genes["gene"].str.match(r"^SNORD116")] if len(genes) else genes
    meta = {"regions": {k: list(v) for k, v in regs.items()}, "snord116_window": list(regs["SNORD116"]),
            "snord116_copies": genes[["start", "end"]].values.tolist() if len(genes) else [],
            "groups": groups, "notes": notes}
    style.write_table(mol, outdir, "Figure4_molecules.tsv.gz")
    style.write_table(calls, outdir, "Figure4_snord116_control_cpg_calls.tsv.gz")
    style.write_table(ent, outdir, "Figure4c_entropy_statistics.tsv")
    style.write_table(pd.DataFrame([sn]), outdir, "Figure4d_snord116_statistics.tsv")
    style.write_table(assign, outdir, "Figure4_control_assignment.tsv")
    (outdir / "tables" / "Figure4_meta.json").write_text(json.dumps(meta, indent=1))
    paths_out = render(outdir, mol, calls, ent, sn, meta)
    print(report(outdir, mol, ent, sn, assign, meta, notes, paths_out))


if __name__ == "__main__":
    main()
