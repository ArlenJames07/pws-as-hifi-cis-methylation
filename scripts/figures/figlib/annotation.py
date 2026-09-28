"""Genomic landmarks in T2T-CHM13v2.0 (0-based, half-open), and track drawing.

BP1-BP3 are the segmental-duplication clusters of the manuscript Methods. The
duplication blocks themselves come from scripts/duplicons (reference/window.sd.bed),
else from params 'segdup_bed', else the three clusters are drawn as blocks."""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

from . import style

CHROM = "chr15"
IC = (22_691_258, 22_693_494)
BP_CLUSTERS = {"BP1": (17_691_439, 20_454_275), "BP2": (20_753_698, 21_183_655),
               "BP3": (25_875_912, 26_632_507)}
FAMILY = re.compile(r"^(SNORD\d+|SNORA\d+)[-_.]?\d*$")
TRACK_GENES = ("MKRN3", "MAGEL2", "NDN", "SNRPN", "SNHG14", "SNORD116", "IPW", "SNORD115", "UBE3A",
               "ATP10A", "GABRB3", "GABRA5", "GABRG3", "OCA2", "HERC2", "CHRNA7", "TRPM1", "APBA2", "TJP1")


def read_gtf_genes(gtf: Path | None, chrom: str, start: int, end: int) -> pd.DataFrame:
    """gene features (or transcripts collapsed by name) overlapping chrom:[start, end)."""
    cols = ["gene", "start", "end", "strand", "type"]
    if gtf is None or not Path(gtf).is_file():
        return pd.DataFrame(columns=cols)
    opener = gzip.open if str(gtf).endswith(".gz") else open
    rows = []
    with opener(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[0] != chrom or f[2] not in ("gene", "transcript"):
                continue
            s, e = int(f[3]) - 1, int(f[4])
            if e < start or s > end:
                continue
            m = re.search(r'gene_name "([^"]+)"', f[8]) or re.search(r'gene "([^"]+)"', f[8]) \
                or re.search(r'gene_id "([^"]+)"', f[8])
            t = re.search(r'gene_(?:bio)?type "([^"]+)"', f[8])
            rows.append((f[2], m.group(1) if m else ".", s, e, f[6], t.group(1) if t else "."))
    g = pd.DataFrame(rows, columns=["feature"] + cols)
    if g.empty:
        return pd.DataFrame(columns=cols)
    if (g["feature"] == "gene").any():
        g = g[g["feature"] == "gene"]
    return (g.groupby("gene").agg(start=("start", "min"), end=("end", "max"), strand=("strand", "first"),
                                  type=("type", "first")).reset_index().sort_values("start")[cols])


def collapse_families(genes: pd.DataFrame) -> pd.DataFrame:
    """SNORD116-1 ... SNORD116-30 -> one 'SNORD116 cluster' span (same for SNORD115, ...)."""
    if genes.empty:
        return genes
    g = genes.copy()
    fam = g["gene"].map(lambda n: FAMILY.match(str(n)).group(1) if FAMILY.match(str(n)) else None)
    single = g[fam.isna()]
    clusters = (g[fam.notna()].assign(gene=fam[fam.notna()]).groupby("gene")
                .agg(start=("start", "min"), end=("end", "max"), strand=("strand", "first"),
                     type=("type", "first"), copies=("start", "size")).reset_index())
    clusters["gene"] = clusters["gene"] + np.where(clusters["copies"] > 1, " cluster", "")
    return pd.concat([single, clusters.drop(columns="copies")], ignore_index=True).sort_values("start")


def landmark_genes(paths, start: int, end: int, names=TRACK_GENES) -> pd.DataFrame:
    """Named genes for tracks: analysis/annotation_genes.tsv if present, else the GTF."""
    table = paths.analysis / "annotation_genes.tsv"
    if table.is_file():
        a = pd.read_csv(table, sep="\t").rename(columns={"label": "gene"})
        a = a[(a["end"] > start) & (a["start"] < end) & (a["gene"] != "IC")]
        a["gene"] = a["gene"].replace({"SNORD116": "SNORD116 cluster", "SNORD115": "SNORD115 cluster"})
        if len(a):
            return a[["gene", "start", "end"]].sort_values("start").reset_index(drop=True)
    g = collapse_families(read_gtf_genes(paths.gtf, CHROM, start, end))
    if g.empty:
        return pd.DataFrame(columns=["gene", "start", "end"])
    keep = g["gene"].str.replace(" cluster", "", regex=False).isin(names)
    return g[keep][["gene", "start", "end"]].reset_index(drop=True)


def read_bed(path: Path | None, chrom: str = CHROM) -> pd.DataFrame:
    if path is None or not Path(path).is_file():
        return pd.DataFrame(columns=["start", "end", "name"])
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[0] == chrom:
                rows.append((int(f[1]), int(f[2]), f[3] if len(f) > 3 else "."))
    return pd.DataFrame(rows, columns=["start", "end", "name"]).sort_values("start").reset_index(drop=True)


def merge_intervals(df: pd.DataFrame, gap: int = 0) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=["start", "end"])
    out = []
    for s, e in df.sort_values("start")[["start", "end"]].itertuples(index=False):
        if out and s <= out[-1][1] + gap:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return pd.DataFrame(out, columns=["start", "end"])


def sd_blocks(paths) -> tuple[pd.DataFrame, str]:
    """(merged duplication blocks, source)."""
    for path, source in ((paths.duplicons / "reference" / "window.sd.bed", "scripts/duplicons window.sd.bed"),
                         (paths.segdup_bed, "params segdup_bed")):
        blocks = read_bed(path)
        if len(blocks):
            return merge_intervals(blocks, gap=1_000), source
    return (pd.DataFrame([(s, e) for s, e in BP_CLUSTERS.values()], columns=["start", "end"]),
            "BP1-BP3 clusters (Methods); no duplication BED found")


def cluster_of(pos: float, slack: int = 50_000) -> str:
    for name, (s, e) in BP_CLUSTERS.items():
        if s - slack <= pos <= e + slack:
            return name
    return "."


def prespecified_regions(paths) -> pd.DataFrame:
    path = paths.analysis / "prespecified_regions.tsv"
    return pd.read_csv(path, sep="\t") if path.is_file() else pd.DataFrame()


# ------------------------------------------------------------------ drawing
def _repel(centres: np.ndarray, widths: np.ndarray, low: float, high: float) -> np.ndarray:
    pos = centres.astype(float).copy()
    for _ in range(400):
        moved = False
        for i in range(1, len(pos)):
            need = (widths[i] + widths[i - 1]) / 2
            if pos[i] - pos[i - 1] < need:
                shift = (need - (pos[i] - pos[i - 1])) / 2
                pos[i - 1] -= shift
                pos[i] += shift
                moved = True
        pos = np.clip(pos, low + widths / 2, high - widths / 2)
        if not moved:
            break
    return pos


def draw_gene_track(ax, genes: pd.DataFrame, xlim: tuple[float, float], fontsize: float | None = None,
                    color: str = "#5A5A5A", label_y: float = 0.05) -> None:
    """Genes as grey bars (Mb x axis) with non-overlapping labels below them."""
    fontsize = fontsize or style.BASE_FONT - 1
    low, high = xlim
    g = genes[(genes["end"] / 1e6 > low) & (genes["start"] / 1e6 < high)].sort_values("start")
    items = []
    for r in g.itertuples(index=False):
        s, e = max(r.start / 1e6, low), min(r.end / 1e6, high)
        ax.add_patch(matplotlib_rect((s, 0.62), max(e - s, (high - low) * 0.0015), 0.22, color))
        items.append(((s + e) / 2, str(r.gene).replace(" cluster", "")))
    if items:
        box = ax.get_position()
        inches = box.width * ax.figure.get_figwidth()
        char = fontsize * 0.6 / 72 * (high - low) / max(inches, 1e-3)
        centres = np.array([c for c, _ in items])
        widths = np.array([len(t) * char + 2 * char for _, t in items])
        pos = _repel(centres, widths, low, high)
        for (c, t), p in zip(items, pos):
            ax.text(p, label_y, t, ha="center", va="bottom", fontsize=fontsize, color=style.INK, style="italic")
            if abs(p - c) > char:
                ax.plot([c, c, p], [0.6, 0.48, 0.36], color=style.LIGHT, lw=0.4)
    ax.set_xlim(low, high)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)


def draw_sd_track(ax, blocks: pd.DataFrame, xlim: tuple[float, float], label_clusters: bool = True,
                  fontsize: float | None = None) -> None:
    fontsize = fontsize or style.BASE_FONT - 1
    low, high = xlim
    for r in blocks.itertuples(index=False):
        s, e = r.start / 1e6, r.end / 1e6
        if e < low or s > high:
            continue
        ax.add_patch(matplotlib_rect((max(s, low), 0.25), min(e, high) - max(s, low), 0.5, style.SD))
    if label_clusters:
        for name, (s, e) in BP_CLUSTERS.items():
            c = (s + e) / 2e6
            if low <= c <= high:
                ax.text(c, 0.5, name, ha="center", va="center", fontsize=fontsize, color=style.INK,
                        fontweight="bold")
    ax.set_xlim(low, high)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)


def matplotlib_rect(xy, w, h, color, **kw):
    from matplotlib.patches import Rectangle
    return Rectangle(xy, w, h, facecolor=color, edgecolor="none", lw=0, **kw)
