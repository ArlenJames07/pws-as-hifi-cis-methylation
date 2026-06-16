#!/usr/bin/env python3
"""
Self-contained Figure 5 generator for the hifi_multiomics_pipeline layout.

This script keeps the original Figure 5 outputs intact and creates a second,
reviewer-oriented version that separates:
1. chr15 deletion architecture,
2. non-chr15 genome-wide burden,
3. genome-wide SV burden,
4. breakpoint-coordinate-aligned methylation distance-decay,
5. compact methylation effect-size summaries.
"""

from __future__ import annotations

import argparse
import itertools
import math
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Optional, Sequence

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy import stats


DEFAULT_OUTDIR = Path("/home/rare/arlen/pws-as-hifi-cis-methylation/scripts/hifi_multiomics_pipeline/06_results")
DEFAULT_FASTA = Path("/home/rare/arlen/reference/chm13v22.fasta")
DEFAULT_GTF = Path("/home/rare/arlen/reference/chm13v22.sorted.gtf")

COHORT = [
    ("001P", "Prader-Willi syndrome", "PWS_DEL"),
    ("002P", "Prader-Willi syndrome", "PWS_DEL"),
    ("005P", "Prader-Willi syndrome", "PWS_DEL"),
    ("006P", "Prader-Willi syndrome", "PWS_DEL"),
    ("007P", "Prader-Willi syndrome", "PWS_DEL"),
    ("004P", "Prader-Willi syndrome", "PWS_mUPD"),
    ("013A", "Angelman syndrome", "AS_DEL"),
    ("014A", "Angelman syndrome", "AS_DEL"),
    ("016A", "Angelman syndrome", "AS_DEL"),
    ("017C", "Unaffected control", "CONTROL"),
    ("018C", "Unaffected control", "CONTROL"),
]

CHROM = "chr15"
PLOT_START = 17_000_000
PLOT_END = 33_000_000
PANEL_GROUP_ORDER = ["PWS_DEL", "AS_DEL", "PWS_mUPD", "CONTROL"]
GROUP_LABEL = {
    "PWS_DEL": "PWS-DEL",
    "AS_DEL": "AS-DEL",
    "PWS_mUPD": "PWS-UPD",
    "CONTROL": "Control",
}
METH_GROUP_LABEL = {
    "PWS_DEL": "PWS-DEL",
    "AS_DEL": "AS-DEL",
    "PWS_mUPD": "PWS-UPD BP-matched",
}
SAMPLE_DISPLAY_PREFIX = {
    "PWS_DEL": "PW",
    "AS_DEL": "AS",
    "PWS_mUPD": "UPD",
    "CONTROL": "CTRL",
}
GROUP_COLORS = {
    "PWS_DEL": "#c03a3a",
    "AS_DEL": "#2f6fb0",
    "PWS_mUPD": "#7a52b3",
    "CONTROL": "#4d4d4d",
}
GROUP_FILLS = {
    "PWS_DEL": "#f4c7c7",
    "AS_DEL": "#cfe0f4",
    "PWS_mUPD": "#e1d5f5",
    "CONTROL": "#d9d9d9",
}
SV_METRIC_COLORS = {
    "DEL": "#c03a3a",
    "INS": "#4c78a8",
    "DUP": "#e08b32",
    "INV": "#7a52b3",
    "BND": "#5a5a5a",
    "TOTAL_COUNT": "#222222",
    "TOTAL_SPAN_MB": "#8c564b",
}
BREAKPOINTS_MB = {
    "BP1": 20.94,
    "BP2": 21.07,
    "BP3": 26.05,
    "BP4": 26.46,
    "BP5": 31.84,
}
PANEL_A_BREAKPOINTS = {
    "BP1": 20_940_000,
    "BP2": 21_070_000,
    "BP3": 26_050_000,
    "BP4": 26_460_000,
    "BP5": 31_840_000,
}
DISTANCE_BINS = [
    ("0-10 kb", 0, 10_000),
    ("10-25 kb", 10_000, 25_000),
    ("25-50 kb", 25_000, 50_000),
    ("50-100 kb", 50_000, 100_000),
]
BP12_GENE_ORDER = ["NIPA1", "NIPA2", "CYFIP1", "TUBGCP5"]
PWS_CORE_GENE_ORDER = ["MKRN3", "MAGEL2", "NDN", "SNRPN", "SNORD116", "IPW", "SNORD115", "UBE3A", "GABRB3"]
EXTRA_007P_FIGURE_GENES = ["APBA2", "TJP1", "FAN1", "TRPM1", "OTUD7A", "CHRNA7", "FMN1", "RYR3", "AVEN"]
SNORD_CLUSTER_PREFIXES = {"SNORD116", "SNORD115"}
PANEL_A_COLORS = {
    "total": "#404040",
    "hap1": "#1B9E77",
    "hap2": "#D95F02",
    "unphased": "#BDBDBD",
    "deleted_pws": "#FDEDEC",
    "deleted_as": "#EAF2FB",
    "boundary_pws": "#C0392B",
    "boundary_as": "#2471A3",
    "canonical_bp": "#8B6914",
    "extended_bp": "#B08968",
    "core_gene": "#1F4E79",
    "extra_gene": "#8C3D1E",
    "core_gene_fill": "#DCEAF6",
    "extra_gene_fill": "#F7E1D6",
}


@dataclass(frozen=True)
class SampleInfo:
    sample_id: str
    syndrome: str
    group: str
    methylation_combined: Optional[Path]
    methylation_hap1: Optional[Path]
    methylation_hap2: Optional[Path]


@dataclass(frozen=True)
class GeneInterval:
    name: str
    start: int
    end: int
    strand: str


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def build_sample_display_labels() -> dict[str, str]:
    labels: dict[str, str] = {}
    counts: dict[str, int] = {}
    for sample_id, _syndrome, group in COHORT:
        counts[group] = counts.get(group, 0) + 1
        labels[sample_id] = f"{SAMPLE_DISPLAY_PREFIX[group]}-{counts[group]}"
    return labels


DISPLAY_SAMPLE_LABELS = build_sample_display_labels()
GROUP_N = {
    group: sum(1 for _sid, _syn, g in COHORT if g == group)
    for group in PANEL_GROUP_ORDER
}


def load_input_inventory(path: Path) -> dict[str, SampleInfo]:
    df = pd.read_csv(path, sep="\t")
    samples: dict[str, SampleInfo] = {}
    for _, row in df.iterrows():
        group = str(row["group"]).replace("-", "_").upper()
        group = {"PWS_UPD": "PWS_mUPD", "CONTROL": "CONTROL"}.get(group, group)
        group = row["group"]
        if group == "PWS-DEL":
            group = "PWS_DEL"
        elif group == "AS-DEL":
            group = "AS_DEL"
        elif group == "PWS-UPD":
            group = "PWS_mUPD"
        elif group == "Control":
            group = "CONTROL"
        samples[str(row["sample_id"])] = SampleInfo(
            sample_id=str(row["sample_id"]),
            syndrome=str(row["syndrome"]),
            group=group,
            methylation_combined=Path(row["methylation_combined"]) if str(row["methylation_combined"]).strip() else None,
            methylation_hap1=Path(row["methylation_hap1"]) if str(row["methylation_hap1"]).strip() else None,
            methylation_hap2=Path(row["methylation_hap2"]) if str(row["methylation_hap2"]).strip() else None,
        )
    return samples


def parse_group(group: str) -> str:
    if group in GROUP_LABEL:
        return group
    if group == "PWS-DEL":
        return "PWS_DEL"
    if group == "AS-DEL":
        return "AS_DEL"
    if group == "PWS-UPD":
        return "PWS_mUPD"
    if group == "Control":
        return "CONTROL"
    return group


def format_pvalue(value: float) -> str:
    if not np.isfinite(value):
        return "n/a"
    if value < 1e-3:
        return f"{value:.1e}"
    return f"{value:.3f}"


def format_effect(value: float) -> str:
    if not np.isfinite(value):
        return "n/a"
    return f"{value:+.3f}"


def weighted_average(values: pd.Series, weights: pd.Series) -> float:
    mask = values.notna() & weights.notna() & (weights > 0)
    if not mask.any():
        return np.nan
    return float(np.average(values[mask], weights=weights[mask]))


def bootstrap_ci(values: Sequence[float], n_boot: int = 5000, seed: int = 13) -> tuple[float, float]:
    arr = np.asarray([v for v in values if np.isfinite(v)], dtype=float)
    if len(arr) < 2:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    draws = rng.choice(arr, size=(n_boot, len(arr)), replace=True)
    means = draws.mean(axis=1)
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def bh_adjust(pvalues: Sequence[float]) -> np.ndarray:
    arr = np.asarray(pvalues, dtype=float)
    q = np.full(len(arr), np.nan, dtype=float)
    valid = np.isfinite(arr)
    if not valid.any():
        return q
    sub = arr[valid]
    order = np.argsort(sub)
    ranked = np.empty_like(sub)
    prev = 1.0
    m = len(sub)
    for j, idx in enumerate(order[::-1], start=1):
        raw = sub[idx] * m / (m - j + 1)
        prev = min(prev, raw)
        ranked[idx] = prev
    q[valid] = ranked
    return q


def kruskal_h(values: Sequence[float], labels: Sequence[str]) -> float:
    vals = np.asarray(values, dtype=float)
    labs = np.asarray(labels)
    ranks = stats.rankdata(vals, method="average")
    unique_labels = []
    groups = []
    for label in labs:
        if label not in unique_labels:
            unique_labels.append(label)
            groups.append(np.where(labs == label)[0])
    n = len(vals)
    h = (12.0 / (n * (n + 1.0))) * sum((ranks[idx].sum() ** 2) / len(idx) for idx in groups) - 3.0 * (n + 1.0)
    _, counts = np.unique(vals, return_counts=True)
    if n > 1:
        tie = 1.0 - ((counts**3 - counts).sum() / float(n**3 - n))
        if tie > 0:
            h /= tie
    return float(h)


def iter_group_partitions(indices: tuple[int, ...], group_sizes: Sequence[int]) -> Iterator[list[tuple[int, ...]]]:
    if len(group_sizes) == 1:
        yield [indices]
        return
    first = group_sizes[0]
    for combo in itertools.combinations(indices, first):
        remaining = tuple(idx for idx in indices if idx not in combo)
        for tail in iter_group_partitions(remaining, group_sizes[1:]):
            yield [tuple(combo)] + tail


def exact_permutation_kruskal(values: Sequence[float], labels: Sequence[str]) -> tuple[float, float]:
    vals = np.asarray(values, dtype=float)
    labs = np.asarray(labels)
    first_seen: list[str] = []
    group_sizes: list[int] = []
    for label in labs:
        if label not in first_seen:
            first_seen.append(label)
            group_sizes.append(int((labs == label).sum()))
    observed = kruskal_h(vals, labs)
    ge = 0
    total = 0
    index_order = tuple(range(len(vals)))
    for groups in iter_group_partitions(index_order, group_sizes):
        perm_labels = np.empty(len(vals), dtype=object)
        for label, group_idx in zip(first_seen, groups):
            for idx in group_idx:
                perm_labels[idx] = label
        stat = kruskal_h(vals, perm_labels)
        ge += int(stat >= observed - 1e-12)
        total += 1
    return observed, ge / float(total)


def epsilon_squared(h_stat: float, n: int, k: int) -> float:
    if n <= k:
        return np.nan
    return max(0.0, float((h_stat - k + 1.0) / (n - k)))


def exact_sign_flip_pvalue(differences: Sequence[float]) -> float:
    diffs = np.asarray([d for d in differences if np.isfinite(d)], dtype=float)
    n = len(diffs)
    if n < 2:
        return np.nan
    observed = abs(diffs.mean())
    ge = 0
    total = 0
    for bits in itertools.product([-1.0, 1.0], repeat=n):
        signed = diffs * np.asarray(bits, dtype=float)
        ge += int(abs(signed.mean()) >= observed - 1e-12)
        total += 1
    return ge / float(total)


def exact_rank_sum_pvalue(left: Sequence[float], right: Sequence[float]) -> float:
    left_vals = np.asarray([v for v in left if np.isfinite(v)], dtype=float)
    right_vals = np.asarray([v for v in right if np.isfinite(v)], dtype=float)
    if len(left_vals) == 0 or len(right_vals) == 0:
        return np.nan
    combined = np.concatenate([left_vals, right_vals])
    ranks = stats.rankdata(combined, method="average")
    n_left = len(left_vals)
    observed = abs(ranks[:n_left].mean() - ranks[n_left:].mean())
    ge = 0
    total = 0
    for left_idx in itertools.combinations(range(len(combined)), n_left):
        mask = np.zeros(len(combined), dtype=bool)
        mask[list(left_idx)] = True
        stat = abs(ranks[mask].mean() - ranks[~mask].mean())
        ge += int(stat >= observed - 1e-12)
        total += 1
    return ge / float(total)


def read_region_bed(
    path: Optional[Path],
    chrom: str,
    start: int,
    end: int,
    meth_col_1based: int = 9,
    cov_col_1based: int = 6,
) -> pd.DataFrame:
    if path is None or not path.exists():
        return pd.DataFrame(columns=["chrom", "start", "end", "mid", "meth", "coverage"])
    proc = subprocess.run(
        [
            "awk",
            "-v",
            f"chrom={chrom}",
            "-v",
            f"start={start}",
            "-v",
            f"end={end}",
            '$1==chrom && $2>=start && $2<end {print}',
            str(path),
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=240,
    )
    rows = []
    meth_idx = meth_col_1based - 1
    cov_idx = cov_col_1based - 1
    for line in proc.stdout.splitlines():
        fields = line.split()
        if len(fields) <= max(2, meth_idx, cov_idx):
            continue
        try:
            start0 = int(float(fields[1]))
            end0 = int(float(fields[2]))
            coverage = float(fields[cov_idx])
            meth = float(fields[meth_idx])
        except Exception:
            continue
        if meth > 1.5:
            meth /= 100.0
        rows.append(
            {
                "chrom": fields[0],
                "start": start0,
                "end": end0,
                "mid": int((start0 + end0) / 2),
                "meth": min(1.0, max(0.0, meth)),
                "coverage": coverage,
            }
        )
    return pd.DataFrame(rows)


def resolve_retained_track(
    sample_id: str,
    group: str,
    rq1_labels: pd.DataFrame,
    phase2_labels: pd.DataFrame,
) -> str:
    if group not in {"PWS_DEL", "AS_DEL"}:
        return "combined"
    if not rq1_labels.empty and {"sample", "layer", "allele_label"}.issubset(rq1_labels.columns):
        subset = rq1_labels[rq1_labels["sample"].astype(str) == sample_id]
        if not subset.empty:
            for _, row in subset.iterrows():
                label = str(row["allele_label"]).lower()
                if group == "PWS_DEL" and "retained_maternal" in label:
                    return str(row["layer"])
                if group == "AS_DEL" and "retained_paternal" in label:
                    return str(row["layer"])
    if not phase2_labels.empty and {"sample", "haplotype1_label", "haplotype2_label"}.issubset(phase2_labels.columns):
        subset = phase2_labels[phase2_labels["sample"].astype(str) == sample_id]
        if not subset.empty:
            row = subset.iloc[0]
            if group == "PWS_DEL" and "retained_maternal" in str(row["haplotype1_label"]).lower():
                return "hap1"
            if group == "PWS_DEL" and "retained_maternal" in str(row["haplotype2_label"]).lower():
                return "hap2"
            if group == "AS_DEL" and "retained_paternal" in str(row["haplotype1_label"]).lower():
                return "hap1"
            if group == "AS_DEL" and "retained_paternal" in str(row["haplotype2_label"]).lower():
                return "hap2"
    return "combined"


def load_chrom_sizes(fasta: Path) -> pd.DataFrame:
    fai = Path(str(fasta) + ".fai")
    rows = []
    with open(fai) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            if chrom.startswith("chr") and (chrom[3:].isdigit() or chrom in {"chrX", "chrY"}):
                rows.append({"chrom": chrom, "length": int(fields[1])})
    order = {f"chr{i}": i for i in range(1, 23)}
    order.update({"chrX": 23, "chrY": 24})
    df = pd.DataFrame(rows)
    df["order"] = df["chrom"].map(order)
    return df.dropna(subset=["order"]).sort_values("order").reset_index(drop=True)


def parse_gtf_attrs(attr: str) -> dict[str, str]:
    out = {}
    for part in attr.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if " " in part:
            key, val = part.split(" ", 1)
            out[key] = val.strip().strip('"')
        elif "=" in part:
            key, val = part.split("=", 1)
            out[key] = val.strip().strip('"')
    return out


def load_chr15_gene_catalog(gtf_path: Path, chrom: str = CHROM) -> pd.DataFrame:
    rows = []
    if not gtf_path.exists():
        return pd.DataFrame(columns=["gene_name", "start", "end", "strand"])
    seen_target = False
    with open(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            if fields[0] == chrom:
                seen_target = True
            elif seen_target:
                break
            else:
                continue
            if fields[2] != "gene":
                continue
            attrs = parse_gtf_attrs(fields[8])
            gene_name = attrs.get("gene") or attrs.get("gene_name") or attrs.get("gene_id")
            if not gene_name:
                continue
            rows.append(
                {
                    "gene_name": gene_name,
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                }
            )
    return pd.DataFrame(rows).sort_values(["start", "end", "gene_name"]).reset_index(drop=True)


def build_gene_interval(gene_catalog: pd.DataFrame, name: str) -> Optional[GeneInterval]:
    if gene_catalog.empty:
        return None
    if name in SNORD_CLUSTER_PREFIXES:
        sub = gene_catalog[gene_catalog["gene_name"].astype(str).str.startswith(name)].copy()
    else:
        sub = gene_catalog[gene_catalog["gene_name"].eq(name)].copy()
    if sub.empty:
        return None
    strand_mode = sub["strand"].mode()
    return GeneInterval(
        name=name,
        start=int(sub["start"].min()),
        end=int(sub["end"].max()),
        strand=str(strand_mode.iat[0]) if not strand_mode.empty else ".",
    )


def build_gene_track(gene_catalog: pd.DataFrame, ordered_names: Sequence[str]) -> list[GeneInterval]:
    genes = []
    for name in ordered_names:
        gene = build_gene_interval(gene_catalog, name)
        if gene is not None:
            genes.append(gene)
    return genes


def assign_gene_rows(genes: Sequence[GeneInterval], pad_bp: int) -> list[tuple[GeneInterval, int]]:
    row_ends: list[int] = []
    placements: list[tuple[GeneInterval, int]] = []
    for gene in sorted(genes, key=lambda item: (item.start, item.end, item.name)):
        for row_idx, last_end in enumerate(row_ends):
            if gene.start > last_end + pad_bp:
                row_ends[row_idx] = gene.end
                placements.append((gene, row_idx))
                break
        else:
            row_ends.append(gene.end)
            placements.append((gene, len(row_ends) - 1))
    return placements


def add_panel_a_breakpoint_guides(ax: plt.Axes, label_names: Sequence[str] = ()) -> None:
    for bp_name, bp_pos in PANEL_A_BREAKPOINTS.items():
        is_canonical = bp_name in {"BP1", "BP2", "BP3"}
        ax.axvline(
            bp_pos / 1e6,
            color=PANEL_A_COLORS["canonical_bp"] if is_canonical else PANEL_A_COLORS["extended_bp"],
            linestyle=":" if is_canonical else "--",
            linewidth=0.75 if is_canonical else 0.65,
            alpha=0.85 if is_canonical else 0.55,
            zorder=0,
        )
    y_offsets = {"BP1": 1.04, "BP2": 1.18, "BP3": 1.04, "BP4": 1.04, "BP5": 1.04}
    x_offsets = {"BP1": -7, "BP2": 7, "BP3": 0, "BP4": 0, "BP5": 0}
    x_align = {"BP1": "right", "BP2": "left", "BP3": "center", "BP4": "center", "BP5": "center"}
    for bp_name in label_names:
        if bp_name not in PANEL_A_BREAKPOINTS:
            continue
        is_canonical = bp_name in {"BP1", "BP2", "BP3"}
        ax.annotate(
            bp_name,
            xy=(PANEL_A_BREAKPOINTS[bp_name] / 1e6, y_offsets.get(bp_name, 1.04)),
            xycoords=ax.get_xaxis_transform(),
            xytext=(x_offsets.get(bp_name, 0), 0),
            textcoords="offset points",
            ha=x_align.get(bp_name, "center"),
            va="bottom",
            fontsize=6.5,
            color=PANEL_A_COLORS["canonical_bp"] if is_canonical else PANEL_A_COLORS["extended_bp"],
            fontweight="bold",
        )


def plot_panel_a_gene_axis(
    ax: plt.Axes,
    genes: Sequence[GeneInterval],
    line_color: str,
    fill_color: str,
    label_names: Sequence[str],
) -> None:
    placements = assign_gene_rows(genes, pad_bp=250_000)
    max_row = max((row_idx for _, row_idx in placements), default=0)
    label_fontsize = 4.0 if max_row >= 3 else 4.8
    for gene, row_idx in placements:
        y = float(max_row - row_idx)
        start_mb = gene.start / 1e6
        end_mb = gene.end / 1e6
        width_mb = max(end_mb - start_mb, 0.008)
        ax.add_patch(
            Rectangle(
                (start_mb, y - 0.18),
                width_mb,
                0.36,
                facecolor=fill_color,
                edgecolor=line_color,
                linewidth=0.85,
                alpha=0.9,
                zorder=2,
            )
        )
        if gene.name in SNORD_CLUSTER_PREFIXES or gene.strand not in {"+", "-"}:
            ax.plot([start_mb, end_mb], [y, y], color=line_color, linewidth=1.2, zorder=3)
        else:
            x0, x1 = (start_mb, end_mb) if gene.strand == "+" else (end_mb, start_mb)
            ax.annotate(
                "",
                xy=(x1, y),
                xytext=(x0, y),
                arrowprops=dict(arrowstyle="-|>", lw=1.1, color=line_color, shrinkA=0, shrinkB=0),
                zorder=4,
            )
        ax.text(
            start_mb + width_mb / 2.0,
            y + 0.18,
            gene.name,
            ha="center",
            va="bottom",
            fontsize=label_fontsize,
            color=line_color,
            fontweight="bold",
            zorder=5,
        )
    add_panel_a_breakpoint_guides(ax, label_names=label_names)
    ax.set_ylim(-0.45, max_row + 0.82)
    ax.set_xlim(PLOT_START / 1e6, PLOT_END / 1e6)
    ax.set_yticks([])
    ax.tick_params(axis="x", labelbottom=False, length=2)
    ax.spines[["left", "right", "top"]].set_visible(False)


def panel_a_deletion_label(row: pd.Series) -> str:
    deletion_type = str(row.get("deletion_type", ""))
    detail = str(row.get("classification_detail", ""))
    if deletion_type == "type I":
        label = "BP1-BP3 type I"
    elif deletion_type == "type II":
        label = "BP2-BP3 type II"
    elif deletion_type == "atypical":
        label = "atypical"
    else:
        label = deletion_type or "not called"
    if deletion_type in {"type I", "type II"} and ("extension" in detail or "BP3/BP4" in detail):
        label += " ext."
    return label


def add_genome_offsets(chrom_sizes: pd.DataFrame) -> tuple[dict[str, float], dict[str, float], dict[str, float]]:
    offsets: dict[str, float] = {}
    centers: dict[str, float] = {}
    ends: dict[str, float] = {}
    offset = 0.0
    for _, row in chrom_sizes.iterrows():
        chrom = str(row["chrom"])
        length = float(row["length"])
        offsets[chrom] = offset
        centers[chrom] = offset + length / 2.0
        offset += length
        ends[chrom] = offset
    return offsets, centers, ends


def prepare_deletion_panel_data(deletion_df: pd.DataFrame) -> pd.DataFrame:
    z = deletion_df.copy()
    z = z[z["patient_id"].astype(str).ne("")]
    z["group"] = z["group"].map(parse_group)
    z["display_label"] = z["patient_id"].map(DISPLAY_SAMPLE_LABELS)
    z["size_mb"] = pd.to_numeric(z["deletion_size"], errors="coerce") / 1e6
    type_label = {
        "type I": "BP1-BP3 type I",
        "type II": "BP2-BP3 type II",
        "atypical": "Atypical",
    }
    z["class_label"] = z["deletion_type"].map(type_label).fillna(z["deletion_type"])
    order = {sid: i for i, (sid, _syn, grp) in enumerate(COHORT) if grp in {"PWS_DEL", "AS_DEL"}}
    z["plot_order"] = z["patient_id"].map(order)
    return z.sort_values("plot_order").reset_index(drop=True)


def plot_panel_a_original(fig: plt.Figure, subplot_spec, deletion_df: pd.DataFrame, coverage_df: pd.DataFrame, gtf_path: Path) -> None:
    del_rows = deletion_df.copy()
    del_rows = del_rows[del_rows["patient_id"].isin([sid for sid, _syn, grp in COHORT if grp in {"PWS_DEL", "AS_DEL"}])]
    gene_catalog = load_chr15_gene_catalog(gtf_path)
    bp12_genes = build_gene_track(gene_catalog, BP12_GENE_ORDER)
    pws_core_genes = build_gene_track(gene_catalog, PWS_CORE_GENE_ORDER)
    extra_genes = build_gene_track(gene_catalog, EXTRA_007P_FIGURE_GENES)
    panel_a_ratios = [1.60, 1.12, 1.02, 0.98] + [1.28] * max(1, len(del_rows))
    gs_a = GridSpecFromSubplotSpec(
        len(panel_a_ratios),
        2,
        subplot_spec=subplot_spec,
        height_ratios=panel_a_ratios,
        width_ratios=[0.34, 1.0],
        hspace=0.16,
        wspace=0.07,
    )
    ax_a_head = fig.add_subplot(gs_a[0, :])
    ax_a_head.axis("off")
    ax_a_head.text(0.0, 0.94, "A. chr15 HiFi coverage and deletion classes", ha="left", va="center", fontsize=12.5, fontweight="bold")
    ax_a_head.text(
        0.0,
        0.66,
        "Dotted guides mark BP1-BP3; dashed guides mark BP4-BP5; colored dashed lines mark sample-specific deletion boundaries.",
        ha="left",
        va="center",
        fontsize=8.3,
        color="#222222",
    )
    panel_a_legend = [
        Line2D([0], [0], color=PANEL_A_COLORS["total"], lw=1.8, label="Total depth"),
        Line2D([0], [0], color=PANEL_A_COLORS["hap1"], lw=1.5, label="hap1 HP-tagged"),
        Line2D([0], [0], color=PANEL_A_COLORS["hap2"], lw=1.5, label="hap2 HP-tagged"),
        Line2D([0], [0], color=PANEL_A_COLORS["unphased"], lw=1.3, ls="--", label="Unphased"),
    ]
    ax_a_head.legend(handles=panel_a_legend, loc="lower left", bbox_to_anchor=(0.0, 0.02), ncol=4, fontsize=8.0, frameon=False, handlelength=2.2, columnspacing=1.2)

    gene_row_labels = {
        1: "BP1-BP2 genes\n(CHM13 GTF)",
        2: "Core PWS/AS\nregion genes",
        3: "Distal genes\nin 007P",
    }
    for row_idx, label in gene_row_labels.items():
        blank = fig.add_subplot(gs_a[row_idx, 0])
        blank.axis("off")
        blank.text(0.0, 0.55, label, ha="left", va="center", fontsize=8.0, fontweight="bold", linespacing=1.15)
    ax_bp12 = fig.add_subplot(gs_a[1, 1])
    ax_core = fig.add_subplot(gs_a[2, 1], sharex=ax_bp12)
    ax_extra = fig.add_subplot(gs_a[3, 1], sharex=ax_bp12)
    plot_panel_a_gene_axis(ax_bp12, bp12_genes, PANEL_A_COLORS["core_gene"], PANEL_A_COLORS["core_gene_fill"], ["BP1", "BP2", "BP3"])
    plot_panel_a_gene_axis(ax_core, pws_core_genes, PANEL_A_COLORS["core_gene"], PANEL_A_COLORS["core_gene_fill"], [])
    plot_panel_a_gene_axis(ax_extra, extra_genes, PANEL_A_COLORS["extra_gene"], PANEL_A_COLORS["extra_gene_fill"], ["BP4", "BP5"])
    ax_extra.set_xticks(np.arange(18, 33, 2))
    ax_extra.tick_params(axis="x", labelbottom=True, labelsize=7, pad=1)

    axes_a = []
    for idx, (_, row) in enumerate(del_rows.iterrows()):
        label_ax = fig.add_subplot(gs_a[idx + 4, 0])
        label_ax.axis("off")
        ax = fig.add_subplot(gs_a[idx + 4, 1], sharex=ax_bp12)
        axes_a.append(ax)
        sid = row["patient_id"]
        z = coverage_df[coverage_df["sample_id"].eq(sid)].copy()
        group = str(row.get("group", ""))
        boundary_color = PANEL_A_COLORS["boundary_as"] if group == "AS_DEL" else PANEL_A_COLORS["boundary_pws"]
        deleted_fill = PANEL_A_COLORS["deleted_as"] if group == "AS_DEL" else PANEL_A_COLORS["deleted_pws"]
        if pd.notna(row["breakpoint_5prime"]) and pd.notna(row["breakpoint_3prime"]):
            ax.axvspan(row["breakpoint_5prime"] / 1e6, row["breakpoint_3prime"] / 1e6, color=deleted_fill, alpha=0.80, zorder=0)
        add_panel_a_breakpoint_guides(ax)
        if not z.empty:
            x = z["bin_mid"] / 1e6
            ax.plot(x, z["total_depth"], color=PANEL_A_COLORS["total"], lw=1.0)
            ax.plot(x, z["hap1_depth"], color=PANEL_A_COLORS["hap1"], lw=0.9)
            ax.plot(x, z["hap2_depth"], color=PANEL_A_COLORS["hap2"], lw=0.9)
            if "unphased_depth" in z:
                ax.plot(x, z["unphased_depth"], color=PANEL_A_COLORS["unphased"], lw=0.8, ls="--")
            ymax = max(1.0, np.nanpercentile(z["total_depth"], 98.5) * 1.18)
            ax.set_ylim(0, ymax)
        if pd.notna(row["breakpoint_5prime"]) and pd.notna(row["breakpoint_3prime"]):
            ax.axvline(row["breakpoint_5prime"] / 1e6, color=boundary_color, ls="--", lw=0.9)
            ax.axvline(row["breakpoint_3prime"] / 1e6, color=boundary_color, ls="--", lw=0.9)
        deletion_mb = pd.to_numeric(row.get("deletion_size"), errors="coerce") / 1e6
        label_ax.text(0.00, 0.74, DISPLAY_SAMPLE_LABELS.get(sid, sid), ha="left", va="center", fontsize=9.2, color=boundary_color, fontweight="bold")
        label_ax.text(0.00, 0.45, panel_a_deletion_label(row), ha="left", va="center", fontsize=6.6, color=boundary_color, fontweight="bold")
        label_ax.text(0.00, 0.17, f"{deletion_mb:.3f} Mb loss", ha="left", va="center", fontsize=6.8, color="#333333")
        ax.set_xlim(PLOT_START / 1e6, PLOT_END / 1e6)
        ax.set_xticks(np.arange(18, 33, 2))
        ax.tick_params(axis="y", labelsize=6, length=2)
        ax.grid(axis="y", color="#eeeeee", lw=0.45)
        ax.spines[["top", "right"]].set_visible(False)
        if idx < len(del_rows) - 1:
            ax.tick_params(axis="x", labelbottom=False)
        else:
            ax.set_xlabel("chr15 position (Mb, CHM13)")
            ax.tick_params(axis="x", labelbottom=True, labelsize=8)
    if axes_a:
        axes_a[len(axes_a) // 2].set_ylabel("Depth per 50 kb bin")


def prepare_cnv_burden(cnv_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    z = cnv_df.copy()
    z["group"] = z["group"].map(parse_group)
    burden = (
        z[
            z["size_bp"].ge(2_000_000)
            & ~z["chrom"].isin(["chrX", "chrY"])
            & ~z["is_canonical_chr15_deletion"].astype(bool)
        ]
        .groupby(["sample_id", "group"], as_index=False)
        .agg(
            nonchr15_large_cnv_count=("chrom", "size"),
            nonchr15_large_cnv_total_mb=("size_bp", lambda s: s.sum() / 1e6),
            nonchr15_large_cnv_max_mb=("size_bp", lambda s: s.max() / 1e6),
        )
    )
    meta = pd.DataFrame(COHORT, columns=["sample_id", "syndrome", "group"])
    burden = meta.merge(burden, on=["sample_id", "group"], how="left").fillna(0)
    burden["display_label"] = burden["sample_id"].map(DISPLAY_SAMPLE_LABELS)
    rows = []
    for metric in ["nonchr15_large_cnv_count", "nonchr15_large_cnv_total_mb"]:
        observed_h, pvalue = exact_permutation_kruskal(burden[metric].to_numpy(), burden["group"].to_numpy())
        rows.append(
            {
                "metric": metric,
                "n_samples": int(len(burden)),
                "n_groups": 4,
                "exact_permutation_p": pvalue,
                "kruskal_h": observed_h,
                "epsilon_squared": epsilon_squared(observed_h, len(burden), 4),
            }
        )
    stats_df = pd.DataFrame(rows)
    stats_df["q_value"] = bh_adjust(stats_df["exact_permutation_p"])
    return burden, stats_df


def prepare_sv_burden(sv_burden_df: pd.DataFrame, sv_calls_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    burden = sv_burden_df.copy()
    burden["group"] = burden["group"].map(parse_group)
    span = sv_calls_df.copy()
    span["group"] = span["group"].map(parse_group)
    span = span.groupby(["sample_id", "group"], as_index=False)["size_bp"].sum().rename(columns={"size_bp": "TOTAL_SPAN_BP"})
    burden = burden.merge(span, on=["sample_id", "group"], how="left")
    burden["TOTAL_COUNT"] = burden[["DEL", "INS", "INV", "DUP", "BND"]].sum(axis=1)
    burden["TOTAL_SPAN_MB"] = burden["TOTAL_SPAN_BP"] / 1e6
    burden["display_label"] = burden["sample_id"].map(DISPLAY_SAMPLE_LABELS)
    metrics = ["DEL", "INS", "DUP", "INV", "BND", "TOTAL_COUNT", "TOTAL_SPAN_MB"]
    rows = []
    for metric in metrics:
        observed_h, pvalue = exact_permutation_kruskal(burden[metric].to_numpy(), burden["group"].to_numpy())
        rows.append(
            {
                "metric": metric,
                "n_samples": int(len(burden)),
                "n_groups": 4,
                "exact_permutation_p": pvalue,
                "kruskal_h": observed_h,
                "epsilon_squared": epsilon_squared(observed_h, len(burden), 4),
            }
        )
    stats_df = pd.DataFrame(rows)
    stats_df["q_value"] = bh_adjust(stats_df["exact_permutation_p"])
    return burden, stats_df


def summarize_region_by_distance(region_df: pd.DataFrame, breakpoint_position: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if region_df.empty:
        for bin_label, low, high in DISTANCE_BINS:
            rows.append(
                {
                    "distance_bin": bin_label,
                    "distance_start_bp": low,
                    "distance_end_bp": high,
                    "distance_mid_kb": (low + high) / 2000.0,
                    "n_cpg": 0,
                    "total_coverage": 0.0,
                    "methylation": np.nan,
                }
            )
        return rows
    tmp = region_df.copy()
    tmp["abs_distance_bp"] = (tmp["mid"] - breakpoint_position).abs()
    for bin_label, low, high in DISTANCE_BINS:
        z = tmp[(tmp["abs_distance_bp"] >= low) & (tmp["abs_distance_bp"] < high)]
        rows.append(
            {
                "distance_bin": bin_label,
                "distance_start_bp": low,
                "distance_end_bp": high,
                "distance_mid_kb": (low + high) / 2000.0,
                "n_cpg": int(len(z)),
                "total_coverage": float(z["coverage"].sum()) if not z.empty else 0.0,
                "methylation": weighted_average(z["meth"], z["coverage"]) if not z.empty else np.nan,
            }
        )
    return rows


def prepare_methylation_distance_decay(profile_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    z = profile_df.copy()
    z["group"] = z["group"].map(parse_group)
    z = z[z["group"].isin(["PWS_DEL", "AS_DEL", "PWS_mUPD"])].copy()
    z["abs_distance_bp"] = pd.to_numeric(z["relative_mid_bp"], errors="coerce").abs()
    z["coverage_weight"] = pd.to_numeric(z["retained_total_coverage"], errors="coerce").fillna(0.0)
    z["track_used"] = z["retained_track"].fillna("combined")
    z["distance_bin"] = pd.cut(
        z["abs_distance_bp"],
        bins=[0, 10_000, 25_000, 50_000, 100_000],
        labels=[label for label, _low, _high in DISTANCE_BINS],
        right=False,
        include_lowest=True,
    )
    z = z[z["distance_bin"].notna()].copy()
    bin_lookup = {label: (low, high) for label, low, high in DISTANCE_BINS}
    z["distance_start_bp"] = z["distance_bin"].map(lambda label: bin_lookup[str(label)][0])
    z["distance_end_bp"] = z["distance_bin"].map(lambda label: bin_lookup[str(label)][1])
    z["distance_mid_kb"] = z["distance_bin"].map(lambda label: (bin_lookup[str(label)][0] + bin_lookup[str(label)][1]) / 2000.0)

    agg_rows: list[dict[str, object]] = []
    for keys, unit_df in z.groupby(
        ["sample_id", "comparison_patient_id", "group", "breakpoint_side", "distance_bin", "distance_start_bp", "distance_end_bp", "distance_mid_kb", "track_used"],
        dropna=False,
        observed=False,
    ):
        sample_id, comparison_id, group, side, bin_label, low, high, mid_kb, track_used = keys
        agg_rows.append(
            {
                "sample_id": sample_id,
                "comparison_patient_id": comparison_id,
                "group": group,
                "group_label": METH_GROUP_LABEL.get(group, GROUP_LABEL[group]),
                "breakpoint_side": side,
                "distance_bin": str(bin_label),
                "distance_start_bp": int(low),
                "distance_end_bp": int(high),
                "distance_mid_kb": float(mid_kb),
                "track_used": track_used,
                "n_cpg": int(pd.to_numeric(unit_df["retained_n_cpg"], errors="coerce").fillna(0).sum()),
                "total_coverage": float(unit_df["coverage_weight"].sum()),
                "methylation": weighted_average(unit_df["retained_methylation"], unit_df["coverage_weight"]),
                "control_mean_methylation": weighted_average(unit_df["control_mean_methylation"], unit_df["coverage_weight"]),
                "control_sd_methylation": float(unit_df["control_sd_methylation"].mean()),
                "control_mean_cpg": float(pd.to_numeric(unit_df["control_mean_cpg"], errors="coerce").mean()),
                "delta_vs_control": weighted_average(unit_df["delta_vs_control"], unit_df["coverage_weight"]),
            }
        )
    merged = pd.DataFrame(agg_rows)
    merged["abs_delta_vs_control"] = merged["delta_vs_control"].abs()

    curve_rows: list[dict[str, object]] = []
    for (group, side, bin_label, distance_mid_kb), z in merged.groupby(
        ["group", "breakpoint_side", "distance_bin", "distance_mid_kb"],
        as_index=False,
    ):
        if group == "PWS_mUPD":
            sample_level = (
                z.groupby(["comparison_patient_id"], as_index=False)
                .agg(delta_vs_control=("delta_vs_control", "mean"))
                .dropna(subset=["delta_vs_control"])
            )
            values = sample_level["delta_vs_control"].to_numpy(dtype=float)
            mean_delta = float(np.nanmean(values)) if len(values) else np.nan
            ci_low, ci_high = (np.nan, np.nan)
            n_samples = 1 if len(values) else 0
            n_regions = int(len(values))
        else:
            sample_level = (
                z.groupby(["sample_id"], as_index=False)
                .agg(delta_vs_control=("delta_vs_control", "mean"))
                .dropna(subset=["delta_vs_control"])
            )
            values = sample_level["delta_vs_control"].to_numpy(dtype=float)
            mean_delta = float(np.nanmean(values)) if len(values) else np.nan
            ci_low, ci_high = bootstrap_ci(values) if len(values) >= 2 else (np.nan, np.nan)
            n_samples = int(len(values))
            n_regions = int(len(values))
        curve_rows.append(
            {
                "group": group,
                "group_label": METH_GROUP_LABEL.get(group, GROUP_LABEL[group]),
                "breakpoint_side": side,
                "distance_bin": bin_label,
                "distance_mid_kb": distance_mid_kb,
                "mean_delta_vs_control": mean_delta,
                "ci_low": ci_low,
                "ci_high": ci_high,
                "n_samples": n_samples,
                "n_regions": n_regions,
            }
        )

    effect_rows: list[dict[str, object]] = []
    for group in ["PWS_DEL", "AS_DEL", "PWS_mUPD"]:
        for side in ["5prime", "3prime"]:
            z = merged[(merged["group"] == group) & (merged["breakpoint_side"] == side)].copy()
            if group == "PWS_mUPD":
                per_unit_rows = []
                for comparison_id, unit_df in z.groupby("comparison_patient_id"):
                    near = unit_df[unit_df["distance_end_bp"] <= 25_000]
                    far = unit_df[unit_df["distance_start_bp"] >= 50_000]
                    if near.empty or far.empty:
                        continue
                    per_unit_rows.append(
                        {
                            "comparison_patient_id": comparison_id,
                            "near_abs_delta": float(near["abs_delta_vs_control"].mean()),
                            "far_abs_delta": float(far["abs_delta_vs_control"].mean()),
                            "near_signed_delta": float(near["delta_vs_control"].mean()),
                            "far_signed_delta": float(far["delta_vs_control"].mean()),
                        }
                    )
                per_unit = pd.DataFrame(per_unit_rows)
                near_abs_mean = float(per_unit["near_abs_delta"].mean()) if not per_unit.empty else np.nan
                far_abs_mean = float(per_unit["far_abs_delta"].mean()) if not per_unit.empty else np.nan
                near_signed_mean = float(per_unit["near_signed_delta"].mean()) if not per_unit.empty else np.nan
                far_signed_mean = float(per_unit["far_signed_delta"].mean()) if not per_unit.empty else np.nan
                effect_rows.append(
                    {
                        "group": group,
                        "group_label": METH_GROUP_LABEL[group],
                        "breakpoint_side": side,
                        "n_units": 1,
                        "n_reference_regions": int(len(per_unit)),
                        "near_abs_delta_mean": near_abs_mean,
                        "far_abs_delta_mean": far_abs_mean,
                        "near_minus_far_abs_delta": near_abs_mean - far_abs_mean if np.isfinite(near_abs_mean) and np.isfinite(far_abs_mean) else np.nan,
                        "ci_low": np.nan,
                        "ci_high": np.nan,
                        "exact_signflip_p": np.nan,
                        "near_signed_delta_mean": near_signed_mean,
                        "far_signed_delta_mean": far_signed_mean,
                        "note": "Descriptive only; BP-coordinate-matched regions in one PWS-UPD sample.",
                    }
                )
                continue

            per_unit_rows = []
            for sample_id, sample_df in z.groupby("sample_id"):
                near = sample_df[sample_df["distance_end_bp"] <= 25_000]
                far = sample_df[sample_df["distance_start_bp"] >= 50_000]
                if near.empty or far.empty:
                    continue
                per_unit_rows.append(
                    {
                        "sample_id": sample_id,
                        "near_abs_delta": float(near["abs_delta_vs_control"].mean()),
                        "far_abs_delta": float(far["abs_delta_vs_control"].mean()),
                        "near_signed_delta": float(near["delta_vs_control"].mean()),
                        "far_signed_delta": float(far["delta_vs_control"].mean()),
                    }
                )
            per_unit = pd.DataFrame(per_unit_rows)
            if per_unit.empty:
                effect_rows.append(
                    {
                        "group": group,
                        "group_label": METH_GROUP_LABEL[group],
                        "breakpoint_side": side,
                        "n_units": 0,
                        "n_reference_regions": 0,
                        "near_abs_delta_mean": np.nan,
                        "far_abs_delta_mean": np.nan,
                        "near_minus_far_abs_delta": np.nan,
                        "ci_low": np.nan,
                        "ci_high": np.nan,
                        "exact_signflip_p": np.nan,
                        "near_signed_delta_mean": np.nan,
                        "far_signed_delta_mean": np.nan,
                        "note": "Insufficient data",
                    }
                )
                continue
            diff = per_unit["near_abs_delta"] - per_unit["far_abs_delta"]
            ci_low, ci_high = bootstrap_ci(diff.to_numpy(dtype=float))
            effect_rows.append(
                {
                    "group": group,
                    "group_label": METH_GROUP_LABEL[group],
                    "breakpoint_side": side,
                    "n_units": int(len(per_unit)),
                    "n_reference_regions": int(len(per_unit)),
                    "near_abs_delta_mean": float(per_unit["near_abs_delta"].mean()),
                    "far_abs_delta_mean": float(per_unit["far_abs_delta"].mean()),
                    "near_minus_far_abs_delta": float(diff.mean()),
                    "ci_low": ci_low,
                    "ci_high": ci_high,
                    "exact_signflip_p": exact_sign_flip_pvalue(diff.to_numpy(dtype=float)),
                    "near_signed_delta_mean": float(per_unit["near_signed_delta"].mean()),
                    "far_signed_delta_mean": float(per_unit["far_signed_delta"].mean()),
                    "note": "Formal exact sign-flip test across samples.",
                }
            )

    effect_df = pd.DataFrame(effect_rows)
    qvals = bh_adjust(effect_df.loc[effect_df["group"].isin(["PWS_DEL", "AS_DEL"]), "exact_signflip_p"].to_numpy(dtype=float))
    effect_df["q_value"] = np.nan
    effect_df.loc[effect_df["group"].isin(["PWS_DEL", "AS_DEL"]), "q_value"] = qvals
    return merged, pd.DataFrame(curve_rows), effect_df


def panel_label(text: str, title: str) -> str:
    return f"{text}. {title}"


CLASS_MARKERS = {
    "BP1-BP3 type I": "o",
    "BP2-BP3 type II": "s",
    "Atypical": "D",
}


def plot_panel_a_schematic(ax: plt.Axes, deletion_df: pd.DataFrame) -> None:
    plot_df = deletion_df.copy()
    plot_df["y"] = np.arange(len(plot_df))[::-1]
    for bp_name, bp_pos in PANEL_A_BREAKPOINTS.items():
        is_core = bp_name in {"BP1", "BP2", "BP3"}
        ax.axvline(
            bp_pos / 1e6,
            color="#8a7a3a" if is_core else "#b79d73",
            lw=1.0 if is_core else 0.8,
            ls=":" if is_core else "--",
            zorder=0,
        )
        ax.text(bp_pos / 1e6, 1.01, bp_name, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=8, color="#7a6a2d", fontweight="bold")
    for _, row in plot_df.iterrows():
        y = float(row["y"])
        group = row["group"]
        start_mb = float(row["breakpoint_5prime"]) / 1e6
        end_mb = float(row["breakpoint_3prime"]) / 1e6
        ax.add_patch(
            Rectangle(
                (start_mb, y - 0.28),
                end_mb - start_mb,
                0.56,
                facecolor=GROUP_FILLS[group],
                edgecolor=GROUP_COLORS[group],
                linewidth=1.4,
                zorder=3,
            )
        )
        ax.text(16.92, y, f"{row['display_label']}  {GROUP_LABEL[group]}", ha="right", va="center", fontsize=9.3, color=GROUP_COLORS[group], fontweight="bold")
        ax.text(33.12, y, f"{row['class_label']} | {row['size_mb']:.2f} Mb", ha="left", va="center", fontsize=8.3, color="#333333")
    ax.set_title(panel_label("A", "Recurrent and atypical chr15q11–q13 deletion architectures"), loc="left", fontsize=13, fontweight="bold")
    ax.set_xlim(17.0, 33.0)
    ax.set_ylim(-0.8, len(plot_df) - 0.2)
    ax.set_yticks([])
    ax.set_xlabel("chr15 position (Mb, CHM13)")
    ax.grid(axis="x", color="#ededed", lw=0.6)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.text(
        0.0,
        -0.16,
        "Most deletion carriers fall into BP1/BP2-to-BP3/BP4-like classes; 007P is the only clear atypical extended deletion.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.8,
        color="#333333",
    )


def plot_deletion_size_panel(ax: plt.Axes, deletion_df: pd.DataFrame) -> None:
    plot_df = deletion_df[deletion_df["group"].isin(["PWS_DEL", "AS_DEL"])].copy()
    group_order = ["PWS_DEL", "AS_DEL"]
    positions = {group: idx for idx, group in enumerate(group_order)}
    rng = np.random.default_rng(21)
    for group in group_order:
        z = plot_df[plot_df["group"] == group]
        if z.empty:
            continue
        x = np.full(len(z), positions[group], dtype=float) + rng.uniform(-0.10, 0.10, size=len(z))
        for xi, (_, row) in zip(x, z.iterrows()):
            marker = CLASS_MARKERS.get(row["class_label"], "o")
            ax.scatter(
                xi,
                row["size_mb"],
                s=70 if marker != "D" else 95,
                marker=marker,
                color=GROUP_COLORS[group],
                edgecolor="white",
                linewidth=0.7,
                zorder=4,
            )
            ax.text(xi + 0.03, row["size_mb"] + 0.08, row["display_label"], fontsize=7.6, color="#333333", ha="left", va="bottom")
        ax.hlines(float(z["size_mb"].median()), positions[group] - 0.20, positions[group] + 0.20, color="#111111", lw=1.2, zorder=3)
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#666666", label="Type I", markersize=7),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#666666", label="Type II", markersize=7),
        Line2D([0], [0], marker="D", color="w", markerfacecolor="#666666", label="Atypical", markersize=7),
    ]
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=7.6, ncol=3, columnspacing=1.0, handletextpad=0.4)
    ax.set_title("B1. chr15 deletion size per sample", loc="left", fontsize=11.5, fontweight="bold")
    ax.set_ylabel("Deletion size (Mb)")
    ax.set_xticks([positions[g] for g in group_order])
    ax.set_xticklabels([f"{GROUP_LABEL[g]}\n(n={sum(plot_df['group']==g)})" for g in group_order], fontsize=8.2)
    ax.grid(axis="y", color="#ededed", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)


def plot_group_distribution(
    ax: plt.Axes,
    df: pd.DataFrame,
    value_col: str,
    title: str,
    y_label: str,
    group_order: Sequence[str] = PANEL_GROUP_ORDER,
    kw_pvalue: Optional[float] = None,
    direct_note: Optional[str] = None,
) -> None:
    positions = np.arange(len(group_order), dtype=float)
    rng = np.random.default_rng(31 + len(value_col))
    for pos, group in zip(positions, group_order):
        z = df[df["group"] == group][value_col].to_numpy(dtype=float)
        z = z[np.isfinite(z)]
        if len(z) >= 2:
            box = ax.boxplot(
                [z],
                positions=[pos],
                widths=0.44,
                patch_artist=True,
                showfliers=False,
                medianprops=dict(color="#111111", linewidth=1.1),
                whiskerprops=dict(color="#888888", linewidth=0.8),
                capprops=dict(color="#888888", linewidth=0.8),
                boxprops=dict(edgecolor="#888888", linewidth=0.8),
            )
            box["boxes"][0].set_facecolor(GROUP_FILLS[group])
            box["boxes"][0].set_alpha(0.65)
        x = np.full(len(z), pos, dtype=float) + rng.uniform(-0.11, 0.11, size=len(z))
        marker = "D" if len(z) == 1 else "o"
        ax.scatter(
            x,
            z,
            s=64 if len(z) == 1 else 52,
            marker=marker,
            color=GROUP_COLORS[group],
            edgecolor="white",
            linewidth=0.6,
            zorder=4,
        )
        if len(z) == 1:
            ax.text(pos, z[0] + max(np.nanmax(z) * 0.02, 0.05), "descriptive", ha="center", va="bottom", fontsize=6.8, color="#555555")
    ax.set_title(title, loc="left", fontsize=11.5, fontweight="bold")
    ax.set_ylabel(y_label)
    ax.set_xticks(positions)
    ax.set_xticklabels([f"{GROUP_LABEL[g]}\n(n={GROUP_N[g]})" for g in group_order], fontsize=8.0)
    ax.grid(axis="y", color="#ededed", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    if kw_pvalue is not None:
        ax.text(
            0.01,
            0.97,
            f"Kruskal-Wallis p = {format_pvalue(kw_pvalue)}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.2,
            bbox=dict(facecolor="white", edgecolor="#dddddd", linewidth=0.4, alpha=0.9),
        )
    if direct_note:
        ax.text(0.99, 0.03, direct_note, transform=ax.transAxes, ha="right", va="bottom", fontsize=8.0, color="#444444")


def plot_sv_stats_table(ax: plt.Axes, stats_df: pd.DataFrame) -> None:
    ax.axis("off")
    ax.text(0.0, 1.0, "Type-specific SV burden summary", ha="left", va="top", fontsize=9.8, fontweight="bold")
    ax.text(0.00, 0.77, "SV type", fontsize=8.6, fontweight="bold")
    ax.text(0.32, 0.77, "KW p", fontsize=8.6, fontweight="bold")
    ax.text(0.52, 0.77, "FDR q", fontsize=8.6, fontweight="bold")
    ax.text(0.74, 0.77, "Interpretation", fontsize=8.6, fontweight="bold")
    row_y = 0.69
    for metric in ["DEL", "INS", "DUP", "INV", "BND"]:
        row = stats_df[stats_df["metric"] == metric].iloc[0]
        interp = "borderline" if float(row["exact_permutation_p"]) < 0.10 else "NS"
        color = "#666666" if interp == "NS" else "#333333"
        ax.text(0.00, row_y, metric, fontsize=8.4)
        ax.text(0.32, row_y, format_pvalue(float(row["exact_permutation_p"])), fontsize=8.4)
        ax.text(0.52, row_y, format_pvalue(float(row["q_value"])), fontsize=8.4)
        ax.text(0.74, row_y, interp, fontsize=8.4, color=color)
        row_y -= 0.10
    ax.text(
        0.0,
        0.10,
        "Main figure emphasizes total SV count and total SV span.\nType-specific detail can be cited as supplementary support.\nPWS-UPD is shown as a single descriptive sample, not a distribution.",
        ha="left",
        va="bottom",
        fontsize=8.4,
        color="#333333",
        linespacing=1.35,
    )


def plot_methylation_decay_main(ax: plt.Axes, curve_df: pd.DataFrame, side: str) -> None:
    subset = curve_df[curve_df["breakpoint_side"] == side].copy()
    x_positions = np.arange(len(DISTANCE_BINS), dtype=float)
    x_labels = [label.replace(" kb", "") for label, _low, _high in DISTANCE_BINS]
    title = "5' breakpoint" if side == "5prime" else "3' breakpoint"
    ax.set_title(title, loc="left", fontsize=11.2, fontweight="bold")
    ax.axhline(0, color="#777777", lw=0.9, ls="--", zorder=0)
    for group in ["PWS_DEL", "AS_DEL", "PWS_mUPD"]:
        z = subset[subset["group"] == group].copy()
        if z.empty:
            continue
        z["bin_order"] = z["distance_bin"].map({label: i for i, (label, _low, _high) in enumerate(DISTANCE_BINS)})
        z = z.sort_values("bin_order")
        y = z["mean_delta_vs_control"].to_numpy(dtype=float)
        if group != "PWS_mUPD":
            ax.fill_between(x_positions, z["ci_low"].to_numpy(dtype=float), z["ci_high"].to_numpy(dtype=float), color=GROUP_COLORS[group], alpha=0.16, zorder=1)
            ax.plot(x_positions, y, color=GROUP_COLORS[group], lw=2.0, marker="o", ms=4.0, zorder=3)
        else:
            ax.plot(x_positions, y, color=GROUP_COLORS[group], lw=1.8, marker="o", ms=4.0, ls="--", zorder=3)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, fontsize=7.4)
    ax.set_ylabel("Δ methylation vs controls")
    ax.set_xlabel("Absolute distance bin from breakpoint (kb)")
    ax.grid(axis="y", color="#ededed", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(-0.08, 0.08)


def plot_effect_forest(ax_plot: plt.Axes, effect_df: pd.DataFrame) -> None:
    order = [
        ("PWS_DEL", "5prime"),
        ("PWS_DEL", "3prime"),
        ("AS_DEL", "5prime"),
        ("AS_DEL", "3prime"),
        ("PWS_mUPD", "5prime"),
        ("PWS_mUPD", "3prime"),
    ]
    rows = []
    for group, side in order:
        match = effect_df[(effect_df["group"] == group) & (effect_df["breakpoint_side"] == side)]
        if not match.empty:
            rows.append(match.iloc[0])
    plot_df = pd.DataFrame(rows).reset_index(drop=True)
    y_positions = np.arange(len(plot_df))[::-1]
    ax_plot.axvline(0, color="#777777", lw=0.9, ls="--", zorder=0)
    for y, (_, row) in zip(y_positions, plot_df.iterrows()):
        color = GROUP_COLORS[row["group"]]
        estimate = row["near_minus_far_abs_delta"]
        if np.isfinite(row["ci_low"]) and np.isfinite(row["ci_high"]):
            ax_plot.hlines(y, row["ci_low"], row["ci_high"], color=color, lw=2.0, zorder=2)
        ax_plot.scatter(estimate, y, s=60, color=color, edgecolor="white", linewidth=0.5, zorder=3, marker="D" if row["group"] == "PWS_mUPD" else "o")
        side_label = str(row["breakpoint_side"]).replace("prime", "'")
        ax_plot.text(-0.115, y, f"{row['group_label']} {side_label}", ha="left", va="center", fontsize=8.5, color="#222222")
    ax_plot.set_title(panel_label("E", "Near-versus-far breakpoint methylation effects are small"), loc="left", fontsize=13, fontweight="bold")
    ax_plot.set_xlabel("Near-minus-far Δ methylation")
    ax_plot.set_xlim(-0.12, 0.08)
    ax_plot.set_ylim(-0.8, len(plot_df) - 0.2)
    ax_plot.set_yticks([])
    ax_plot.grid(axis="x", color="#ededed", lw=0.6)
    ax_plot.spines[["top", "right", "left"]].set_visible(False)

def add_pairwise_brackets(
    ax: plt.Axes,
    values_by_group: dict[str, np.ndarray],
    positions: dict[str, float],
    comparisons: Sequence[tuple[str, str]],
    fontsize: float = 6.6,
) -> None:
    finite_vals = np.concatenate([vals[np.isfinite(vals)] for vals in values_by_group.values() if len(vals[np.isfinite(vals)]) > 0]) if values_by_group else np.array([])
    if len(finite_vals) == 0:
        return
    ymin = float(np.nanmin(finite_vals))
    ymax = float(np.nanmax(finite_vals))
    yrange = max(ymax - ymin, 1.0)
    bracket_height = yrange * 0.04
    step = yrange * 0.12
    start = ymax + yrange * 0.08
    for idx, (left_group, right_group) in enumerate(comparisons):
        left = values_by_group.get(left_group, np.array([], dtype=float))
        right = values_by_group.get(right_group, np.array([], dtype=float))
        pvalue = exact_rank_sum_pvalue(left, right)
        x1 = positions[left_group]
        x2 = positions[right_group]
        y = start + idx * step
        ax.plot([x1, x1, x2, x2], [y, y + bracket_height, y + bracket_height, y], color="#222222", lw=0.8, clip_on=False)
        ax.text((x1 + x2) / 2.0, y + bracket_height + yrange * 0.015, f"p={format_pvalue(pvalue)}", ha="center", va="bottom", fontsize=fontsize, color="#222222")
    ax.set_ylim(ymin - yrange * 0.05, start + len(comparisons) * step + yrange * 0.10)


def plot_cnv_panel(
    ax_manhattan: plt.Axes,
    cnv_df: pd.DataFrame,
    cnv_burden: pd.DataFrame,
    cnv_stats: pd.DataFrame,
    chrom_sizes: pd.DataFrame,
    deletion_df: pd.DataFrame,
) -> None:
    offsets, centers, ends = add_genome_offsets(chrom_sizes)
    autosomes = [f"chr{i}" for i in range(1, 23)]
    plot_df = cnv_df[
        cnv_df["size_bp"].ge(2_000_000)
        & ~cnv_df["chrom"].isin(["chrX", "chrY"])
    ].copy()
    plot_df["group"] = plot_df["group"].map(parse_group)
    label_lookup = (
        deletion_df[["patient_id", "display_label"]]
        .drop_duplicates()
        .rename(columns={"patient_id": "sample_id"})
    )
    plot_df = plot_df.merge(label_lookup, on="sample_id", how="left")
    if "genome_mid" not in plot_df.columns:
        plot_df["genome_mid"] = plot_df.apply(lambda r: offsets[str(r["chrom"])] + float(r["mid"]), axis=1)
    plot_df["size_mb"] = plot_df["size_bp"] / 1e6
    ax_manhattan.set_title(panel_label("B", "chr15 deletion dominates large CNV signal"), loc="left", fontsize=12.5, fontweight="bold")
    if "chr15" in offsets:
        ax_manhattan.axvspan(offsets["chr15"], ends["chr15"], color="#f4efe0", alpha=0.9, zorder=0)
        ax_manhattan.text(
            centers["chr15"],
            1.02,
            "chr15",
            transform=ax_manhattan.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=8.3,
            color="#7a6a2d",
            fontweight="bold",
        )
    for chrom in autosomes[:-1]:
        ax_manhattan.axvline(ends[chrom], color="#efefef", lw=0.6, zorder=0)
    background = plot_df[~plot_df["is_canonical_chr15_deletion"].astype(bool)].copy()
    if not background.empty:
        ax_manhattan.scatter(
            background["genome_mid"],
            background["size_mb"],
            s=np.clip(background["size_mb"] * 7.0, 16, 54),
            color="#b7b7b7",
            alpha=0.75,
            edgecolors="none",
            zorder=2,
        )
    canonical_rows: list[pd.Series] = []
    canonical = plot_df[plot_df["is_canonical_chr15_deletion"].astype(bool)].copy()
    for _, row in canonical.iterrows():
        group = str(row["group"])
        ax_manhattan.scatter(
            float(row["genome_mid"]),
            float(row["size_mb"]),
            s=float(np.clip(float(row["size_mb"]) * 16.0, 110, 260)),
            marker="o",
            color=GROUP_COLORS[group],
            edgecolor="#202020",
            linewidth=0.6,
            alpha=0.95,
            zorder=6,
        )
        canonical_rows.append(row)
    count_row = cnv_stats[cnv_stats["metric"] == "nonchr15_large_cnv_count"].iloc[0]
    span_row = cnv_stats[cnv_stats["metric"] == "nonchr15_large_cnv_total_mb"].iloc[0]
    ax_manhattan.text(
        0.01,
        0.96,
        "Non-chr15 burden\n"
        f"count p = {format_pvalue(float(count_row['exact_permutation_p']))}, q = {format_pvalue(float(count_row['q_value']))}\n"
        f"span p = {format_pvalue(float(span_row['exact_permutation_p']))}, q = {format_pvalue(float(span_row['q_value']))}",
        transform=ax_manhattan.transAxes,
        ha="left",
        va="top",
        fontsize=7.8,
        bbox=dict(facecolor="white", edgecolor="#dddddd", linewidth=0.4, alpha=0.92),
    )
    ax_manhattan.set_xticks([centers[c] for c in autosomes])
    ax_manhattan.set_xticklabels([c.replace("chr", "") for c in autosomes], fontsize=6.8)
    ax_manhattan.set_ylabel("CNV size (Mb)")
    ax_manhattan.set_xlabel("Chromosome")
    ax_manhattan.grid(axis="y", color="#ededed", lw=0.6)
    ax_manhattan.spines[["top", "right"]].set_visible(False)
    ax_manhattan.set_ylim(0, max(16, float(plot_df["size_mb"].max()) * 1.08 if not plot_df.empty else 16))
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=GROUP_COLORS["PWS_DEL"], markeredgecolor="#202020", label="PWS-DEL chr15 deletion", markersize=7),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=GROUP_COLORS["AS_DEL"], markeredgecolor="#202020", label="AS-DEL chr15 deletion", markersize=7),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#b7b7b7", label="Other autosomal CNVs >= 2 Mb", markersize=6),
    ]
    ax_manhattan.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.00), frameon=False, fontsize=7.3, ncol=1, columnspacing=0.9, handletextpad=0.5, borderaxespad=0.0)


def plot_sv_facet(ax: plt.Axes, burden_df: pd.DataFrame, stats_df: pd.DataFrame, metric: str, show_xticks: bool) -> None:
    positions = np.arange(len(PANEL_GROUP_ORDER), dtype=float)
    title_map = {
        "DEL": "DEL",
        "INS": "INS",
        "DUP": "DUP",
        "INV": "INV",
        "BND": "BND",
        "TOTAL_COUNT": "Total count",
        "TOTAL_SPAN_MB": "Total span",
    }
    rng = np.random.default_rng(7 + sum(ord(ch) for ch in metric))
    for idx, group in enumerate(PANEL_GROUP_ORDER):
        z = burden_df[burden_df["group"] == group][metric].to_numpy(dtype=float)
        z = z[np.isfinite(z)]
        if len(z) >= 2:
            box = ax.boxplot(
                [z],
                positions=[positions[idx]],
                widths=0.48,
                patch_artist=True,
                showfliers=False,
                medianprops=dict(color="#111111", linewidth=1.15),
                whiskerprops=dict(color="#888888", linewidth=0.8),
                capprops=dict(color="#888888", linewidth=0.8),
                boxprops=dict(edgecolor="#888888", linewidth=0.8),
            )
            box["boxes"][0].set_facecolor(GROUP_FILLS[group])
            box["boxes"][0].set_alpha(0.6)
        x = np.full(len(z), positions[idx]) + rng.uniform(-0.12, 0.12, size=len(z))
        ax.scatter(
            x,
            z,
            s=48 if len(z) == 1 else 36,
            marker="D" if len(z) == 1 else "o",
            color=GROUP_COLORS[group],
            alpha=0.9,
            edgecolor="white",
            linewidth=0.4,
            zorder=4,
        )
    ax.set_title(title_map.get(metric, metric.replace("_", " ")), fontsize=9.5, fontweight="bold", pad=6)
    ax.grid(axis="y", color="#ededed", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    if show_xticks:
        ax.set_xticks(positions)
        ax.set_xticklabels([f"{GROUP_LABEL[g]}\n(n={GROUP_N[g]})" for g in PANEL_GROUP_ORDER], fontsize=7.6)
    else:
        ax.set_xticks(positions)
        ax.set_xticklabels([])


def plot_sv_panel(fig: plt.Figure, subplot_spec, burden_df: pd.DataFrame, stats_df: pd.DataFrame) -> None:
    gs = GridSpecFromSubplotSpec(2, 2, subplot_spec=subplot_spec, height_ratios=[0.28, 1.0], hspace=0.05, wspace=0.30)
    metrics = ["TOTAL_COUNT", "TOTAL_SPAN_MB"]
    header_left = fig.add_subplot(gs[0, 0])
    header_right = fig.add_subplot(gs[0, 1])
    header_left.axis("off")
    header_right.axis("off")
    axes = []
    for idx, metric in enumerate(metrics):
        ax = fig.add_subplot(gs[1, idx])
        plot_sv_facet(ax, burden_df, stats_df, metric, show_xticks=True)
        if metric == "TOTAL_SPAN_MB":
            ax.set_ylabel("SV span (Mb)")
        elif metric == "TOTAL_COUNT":
            ax.set_ylabel("Count per sample")
        axes.append(ax)
    header_left.text(
        0.0,
        0.98,
        panel_label("C", "Global SV burden does not clearly separate diagnostic groups"),
        ha="left",
        va="top",
        fontsize=12.5,
        fontweight="bold",
    )
    header_left.text(
        0.0,
        0.12,
        f"Total count: p = {format_pvalue(float(stats_df[stats_df['metric'] == 'TOTAL_COUNT']['exact_permutation_p'].iloc[0]))}",
        ha="left",
        va="bottom",
        fontsize=7.3,
        color="#333333",
    )
    header_right.text(
        0.0,
        0.12,
        f"Total span: p = {format_pvalue(float(stats_df[stats_df['metric'] == 'TOTAL_SPAN_MB']['exact_permutation_p'].iloc[0]))}",
        ha="left",
        va="bottom",
        fontsize=7.3,
        color="#333333",
    )


def plot_methylation_distance_decay(
    ax: plt.Axes,
    curve_df: pd.DataFrame,
    side: str,
) -> None:
    subset = curve_df[curve_df["breakpoint_side"] == side].copy()
    x_positions = np.arange(len(DISTANCE_BINS), dtype=float)
    x_labels = [label.replace(" kb", "") for label, _low, _high in DISTANCE_BINS]
    title = "5' breakpoint-aligned methylation decay" if side == "5prime" else "3' breakpoint-aligned methylation decay"
    ax.set_title(title, loc="left", fontsize=10.8, fontweight="bold")
    ax.axhline(0, color="#666666", lw=0.9, ls="--", zorder=0)
    for group in ["PWS_DEL", "AS_DEL", "PWS_mUPD"]:
        z = subset[subset["group"] == group].copy()
        if z.empty:
            continue
        z["bin_order"] = z["distance_bin"].map({label: i for i, (label, _low, _high) in enumerate(DISTANCE_BINS)})
        z = z.sort_values("bin_order")
        y = z["mean_delta_vs_control"].to_numpy(dtype=float)
        color = GROUP_COLORS[group]
        if group != "PWS_mUPD":
            low = z["ci_low"].to_numpy(dtype=float)
            high = z["ci_high"].to_numpy(dtype=float)
            ax.fill_between(x_positions, low, high, color=color, alpha=0.15, zorder=1)
            ax.plot(x_positions, y, color=color, lw=2.0, marker="o", ms=4.0, zorder=3)
        else:
            ax.plot(x_positions, y, color=color, lw=1.9, marker="o", ms=4.0, ls="--", zorder=3)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, fontsize=7.2)
    ax.set_ylabel("Δ methylation vs controls")
    ax.set_xlabel("Absolute distance bin from breakpoint (kb)")
    ax.grid(axis="y", color="#ededed", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(-0.12, 0.12)


def plot_methylation_effect_summary(ax: plt.Axes, effect_df: pd.DataFrame) -> None:
    order = [
        ("PWS_DEL", "5prime"),
        ("PWS_DEL", "3prime"),
        ("AS_DEL", "5prime"),
        ("AS_DEL", "3prime"),
        ("PWS_mUPD", "5prime"),
        ("PWS_mUPD", "3prime"),
    ]
    rows = []
    for group, side in order:
        match = effect_df[(effect_df["group"] == group) & (effect_df["breakpoint_side"] == side)]
        if not match.empty:
            rows.append(match.iloc[0])
    plot_df = pd.DataFrame(rows).reset_index(drop=True)
    y_positions = np.arange(len(plot_df))[::-1]
    ax.axvline(0, color="#777777", lw=0.9, ls="--", zorder=0)
    for y, (_, row) in zip(y_positions, plot_df.iterrows()):
        color = GROUP_COLORS[row["group"]]
        estimate = row["near_minus_far_abs_delta"]
        if np.isfinite(row["ci_low"]) and np.isfinite(row["ci_high"]):
            ax.hlines(y, row["ci_low"], row["ci_high"], color=color, lw=2.0, zorder=2)
        ax.scatter(estimate, y, s=56, color=color, edgecolor="white", linewidth=0.5, zorder=3)
        side_label = str(row["breakpoint_side"]).replace("prime", "'")
        left_label = f"{row['group_label']} {side_label}"
        right_label = (
            f"n={int(row['n_units'])} | near={row['near_abs_delta_mean']:.3f} | far={row['far_abs_delta_mean']:.3f}\n"
            f"p={format_pvalue(row['exact_signflip_p'])} | q={format_pvalue(row['q_value'])}"
            if row["group"] != "PWS_mUPD"
            else f"n=1 sample | refs={int(row['n_reference_regions'])}\ndescriptive only"
        )
        ax.text(-0.102, y, left_label, ha="left", va="center", fontsize=8.1, color="#222222")
        ax.text(0.082, y, right_label, ha="left", va="center", fontsize=7.3, color="#333333")
    ax.set_title(panel_label("E", "Breakpoint-associated methylation effect summary"), loc="left", fontsize=12.5, fontweight="bold")
    ax.set_xlabel("Near-minus-far |Δ methylation| vs controls")
    ax.set_xlim(-0.11, 0.16)
    ax.set_ylim(-0.8, len(plot_df) - 0.2)
    ax.set_yticks([])
    ax.grid(axis="x", color="#ededed", lw=0.6)
    ax.spines[["top", "right", "left"]].set_visible(False)


def render_figure(
    out_png: Path,
    out_pdf: Path,
    deletion_df: pd.DataFrame,
    coverage_df: pd.DataFrame,
    cnv_df: pd.DataFrame,
    cnv_burden: pd.DataFrame,
    cnv_stats: pd.DataFrame,
    chrom_sizes: pd.DataFrame,
    sv_burden: pd.DataFrame,
    sv_stats: pd.DataFrame,
    methyl_curve: pd.DataFrame,
    methyl_effect: pd.DataFrame,
    gtf_path: Path,
) -> None:
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
        }
    )
    fig = plt.figure(figsize=(18.0, 22.4), constrained_layout=False)
    outer = GridSpec(5, 1, figure=fig, height_ratios=[2.80, 0.92, 0.86, 0.95, 0.70], hspace=0.42)

    plot_panel_a_original(fig, outer[0, 0], deletion_df, coverage_df, gtf_path)

    ax_b = fig.add_subplot(outer[1, 0])
    plot_cnv_panel(ax_b, cnv_df, cnv_burden, cnv_stats, chrom_sizes, deletion_df)

    plot_sv_panel(fig, outer[2, 0], sv_burden, sv_stats)

    gs_d = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[3, 0], width_ratios=[1.0, 1.0], wspace=0.25)
    ax_d1 = fig.add_subplot(gs_d[0, 0])
    ax_d2 = fig.add_subplot(gs_d[0, 1])
    plot_methylation_decay_main(ax_d1, methyl_curve, "5prime")
    plot_methylation_decay_main(ax_d2, methyl_curve, "3prime")
    ax_d1.text(-0.16, 1.10, "D.", transform=ax_d1.transAxes, fontsize=13, fontweight="bold", ha="left", va="bottom")
    ax_d1.text(0.0, 1.10, "Limited evidence for breakpoint-proximal methylation decay", transform=ax_d1.transAxes, fontsize=13, fontweight="bold", ha="left", va="bottom")
    handles = [
        Line2D([0], [0], color=GROUP_COLORS["PWS_DEL"], lw=2.0, marker="o", label="PWS-DEL"),
        Line2D([0], [0], color=GROUP_COLORS["AS_DEL"], lw=2.0, marker="o", label="AS-DEL"),
        Line2D([0], [0], color=GROUP_COLORS["PWS_mUPD"], lw=1.8, marker="o", ls="--", label="PWS-UPD BP-matched"),
    ]
    ax_d2.legend(handles=handles, loc="upper right", frameon=False, fontsize=8.0)

    ax_e = fig.add_subplot(outer[4, 0])
    plot_effect_forest(ax_e, methyl_effect)

    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def render_supplementary_coverage_figure(
    out_png: Path,
    out_pdf: Path,
    deletion_df: pd.DataFrame,
    coverage_df: pd.DataFrame,
    gtf_path: Path,
) -> None:
    fig = plt.figure(figsize=(18.0, 13.5), constrained_layout=False)
    outer = GridSpec(1, 1, figure=fig)
    plot_panel_a_original(fig, outer[0, 0], deletion_df, coverage_df, gtf_path)
    fig.suptitle("Supplementary Figure. Full chr15 coverage tracks for all deletion carriers", x=0.01, y=0.995, ha="left", fontsize=13, fontweight="bold")
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_report(
    out_path: Path,
    deletion_df: pd.DataFrame,
    cnv_burden: pd.DataFrame,
    cnv_stats: pd.DataFrame,
    sv_stats: pd.DataFrame,
    methyl_effect: pd.DataFrame,
) -> None:
    count_stats = cnv_stats[cnv_stats["metric"] == "nonchr15_large_cnv_count"].iloc[0]
    total_stats = cnv_stats[cnv_stats["metric"] == "nonchr15_large_cnv_total_mb"].iloc[0]
    sv_total = sv_stats[sv_stats["metric"] == "TOTAL_COUNT"].iloc[0]
    sv_span = sv_stats[sv_stats["metric"] == "TOTAL_SPAN_MB"].iloc[0]
    bnd_stat = sv_stats[sv_stats["metric"] == "BND"].iloc[0]
    atypical = deletion_df.loc[deletion_df["deletion_type"] == "atypical", "patient_id"].tolist()
    lines = []
    lines.append("# Figure 5 v7 report\n")
    lines.append("## Revised plotting strategy\n")
    lines.append("- Build the main figure around a single message: chr15q11-q13 deletion architecture is clear, whereas genome-wide CNV/SV burden and breakpoint-aligned methylation effects are modest.")
    lines.append("- Move dense full-sample chr15 coverage tracks out of the main figure and into a supplementary companion panel.")
    lines.append("- Keep only sample-level burden summaries that are interpretable at a glance: chr15 deletion size, non-chr15 large CNV burden, total SV count, and total SV span.")
    lines.append("- Treat PWS-UPD as a breakpoint-coordinate-matched descriptive reference rather than as a true breakpoint-flanking deletion analysis.\n")
    lines.append("## Figure layout\n")
    lines.append("- Panel A: chr15 HiFi coverage and deletion classes using the v1-style coverage track layout.")
    lines.append("- Panel B: single Manhattan-style CNV panel with chr15 highlighted and compact non-chr15 burden statistics.")
    lines.append("- Panel C: total SV count and total SV span per sample, with p-values displayed in a separate header band above the plots.")
    lines.append("- Panel D: breakpoint-aligned methylation difference versus controls across 0-10, 10-25, 25-50, and 50-100 kb bins.")
    lines.append("- Panel E: forest-style near-versus-far breakpoint methylation effect summary.")
    lines.append("- Supplementary coverage figure: full chr15 HiFi coverage tracks for all deletion carriers.\n")
    lines.append("## Panel-specific recommendations implemented\n")
    lines.append("- Panel A retains the coverage-track view but uses cleaner left-side labels and spacing.")
    lines.append("- Panel B now uses one condensed Manhattan-style panel with group-colored chr15 deletion points, separated sample callouts, and the chr15 label above the plotting area.")
    lines.append("- Panel C now keeps only total count and total span, and moves p-values into a dedicated header band rather than inside the plotting panels.")
    lines.append("- Panels D-E use conservative language and emphasize effect sizes, confidence intervals, and null-crossing intervals rather than implying strong distance-decay.")
    lines.append("- Multiple-testing correction is applied across SV metrics and across the four formal methylation near-versus-far tests.\n")
    lines.append("## Revised script structure\n")
    lines.append("1. Load existing Figure 5 input tables and methylation file inventory.")
    lines.append("2. Recompute sample-level non-chr15 CNV burden and SV burden statistics.")
    lines.append("3. Re-bin the existing breakpoint-coordinate-aligned methylation table into distance-decay intervals.")
    lines.append("4. Build distance-bin summaries, exact permutation/sign-flip statistics, and bootstrap confidence intervals.")
    lines.append("5. Render a simplified main figure plus a supplementary full-coverage figure and write analysis tables and narrative report.\n")
    lines.append("## Suggested panel titles\n")
    lines.append("- A. chr15 HiFi coverage and deletion classes")
    lines.append("- B. chr15 deletion dominates large CNV signal")
    lines.append("- C. Global SV burden does not clearly separate diagnostic groups")
    lines.append("- D. Limited evidence for breakpoint-proximal methylation decay")
    lines.append("- E. Near-versus-far breakpoint methylation effects are small\n")
    lines.append("## Key quantitative observations\n")
    lines.append(
        f"- Non-chr15 CNV count burden remains weakly separated across groups "
        f"(`exact permutation p={format_pvalue(count_stats['exact_permutation_p'])}`, "
        f"`q={format_pvalue(count_stats['q_value'])}`, `epsilon^2={count_stats['epsilon_squared']:.2f}`)."
    )
    lines.append(
        f"- Non-chr15 CNV total span is similarly non-dominant "
        f"(`exact permutation p={format_pvalue(total_stats['exact_permutation_p'])}`)."
    )
    lines.append(
        f"- Total SV count does not separate groups strongly "
        f"(`p={format_pvalue(sv_total['exact_permutation_p'])}`, `q={format_pvalue(sv_total['q_value'])}`, "
        f"`epsilon^2={sv_total['epsilon_squared']:.2f}`), and the same is true for total SV span "
        f"(`p={format_pvalue(sv_span['exact_permutation_p'])}`)."
    )
    lines.append(
        f"- `BND` burden is the only nominal SV signal "
        f"(`p={format_pvalue(bnd_stat['exact_permutation_p'])}`), but it does not remain significant after FDR "
        f"(`q={format_pvalue(bnd_stat['q_value'])}`)."
    )
    if atypical:
        lines.append(
            f"- The deletion architecture remains dominated by recurrent BP1/BP2-to-BP3/BP4-like events with `{', '.join(atypical)}` as the only clearly atypical extended deletion."
        )
    lines.append("- Methylation interpretation is intentionally conservative because parental-haplotype assignments are incomplete for most deletion carriers; the v7 analysis uses retained haplotype where labeled and combined methylation otherwise.\n")
    lines.append("## Suggested caption text\n")
    lines.append(
        "Figure 5. Structural deletion architecture, genome-wide structural burden, and breakpoint-associated methylation. "
        "(A) chr15 deletion intervals are shown schematically per sample with canonical BP1-BP5 guides, separating recurrent BP1/BP2-to-BP3/BP4-like classes from the single atypical extended deletion. "
        "(B) chr15 deletion size is shown per sample together with non-chr15 large autosomal CNV burden (`>=2 Mb`); canonical chr15 deletions are the dominant CNV events, whereas non-chr15 burden overlaps across groups. "
        "(C) Global SV burden is summarized as total SV count, total SV span, and a compact type-specific statistics table; no SV burden metric survives FDR correction. "
        "(D) Breakpoint-aligned methylation is summarized as signed delta methylation relative to matched controls across distance bins from `0-10 kb` to `50-100 kb`. "
        "(E) Near-versus-far breakpoint methylation effects are summarized per group and breakpoint side. PWS-UPD is shown as a breakpoint-coordinate-matched descriptive reference rather than a true breakpoint-flanking deletion analysis. "
        "A supplementary panel provides the full chr15 HiFi coverage tracks for all deletion carriers. Across panels, the figure supports a recurrent chr15 structural mechanism with limited evidence that global genome-wide CNV/SV burden or broad breakpoint-flanking methylation change is the primary discriminating signal."
    )
    lines.append("\n## Results-ready interpretation template\n")
    lines.append(
        "Deletion carriers showed a predominantly recurrent chr15 architecture, with most samples mapping to BP1/BP2-to-BP3/BP4-like classes and a single atypical extended deletion. "
        "Outside the canonical chr15 event, large autosomal CNV burden did not separate groups strongly (`exact permutation p="
        f"{format_pvalue(count_stats['exact_permutation_p'])}` for count burden), and genome-wide SV burden showed similarly shallow differences (`total SV count p={format_pvalue(sv_total['exact_permutation_p'])}`; all SV burden `q>=0.05`). "
        "Breakpoint-coordinate-aligned methylation differences relative to controls remained small overall and were most consistent with weak, local deviations rather than a broad or uniform epigenetic bleed effect. "
        "The UPD sample was analyzed only at canonical breakpoint-matched coordinates and is therefore interpreted descriptively rather than as evidence for true breakpoint-flanking methylation change."
    )
    lines.append("\n## Methylation effect summary\n")
    for _, row in methyl_effect.iterrows():
        side_label = row["breakpoint_side"].replace("prime", "'")
        lines.append(
            f"- {row['group_label']} {side_label}: near `|Δ|={row['near_abs_delta_mean']:.3f}`, "
            f"far `|Δ|={row['far_abs_delta_mean']:.3f}`, near-minus-far `{format_effect(row['near_minus_far_abs_delta'])}`, "
            f"p={format_pvalue(row['exact_signflip_p'])}, q={format_pvalue(row['q_value'])}; {row['note']}"
        )
    out_path.write_text("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    parser.add_argument("--gtf", type=Path, default=DEFAULT_GTF)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = ensure_dir(args.outdir)
    tables_dir = ensure_dir(outdir / "tables")
    figures_dir = ensure_dir(outdir / "figures")
    reports_dir = ensure_dir(outdir / "reports")

    deletion_df = prepare_deletion_panel_data(pd.read_csv(tables_dir / "Figure5A_deletion_breakpoint_characterization.tsv", sep="\t"))
    coverage_df = pd.read_csv(tables_dir / "Figure5A_haplotype_coverage_tracks.tsv.gz", sep="\t")
    cnv_df = pd.read_csv(tables_dir / "Figure5B_genomewide_cnv_calls.tsv.gz", sep="\t")
    sv_calls_df = pd.read_csv(tables_dir / "Figure5C_sv_calls.tsv.gz", sep="\t")
    sv_burden_df = pd.read_csv(tables_dir / "Figure5C_sv_burden_by_sample.tsv", sep="\t")
    methyl_profile_df = pd.read_csv(tables_dir / "Figure5D_breakpoint_flanking_methylation_profile.tsv.gz", sep="\t")

    cnv_burden, cnv_stats = prepare_cnv_burden(cnv_df)
    sv_burden, sv_stats = prepare_sv_burden(sv_burden_df, sv_calls_df)
    methyl_bins, methyl_curve, methyl_effect = prepare_methylation_distance_decay(methyl_profile_df)
    chrom_sizes = load_chrom_sizes(args.fasta)

    deletion_df.to_csv(tables_dir / "Figure5_v7_chr15_deletion_classes.tsv", sep="\t", index=False)
    cnv_burden.to_csv(tables_dir / "Figure5_v7_nonchr15_cnv_burden.tsv", sep="\t", index=False)
    cnv_stats.to_csv(tables_dir / "Figure5_v7_nonchr15_cnv_stats.tsv", sep="\t", index=False)
    sv_burden.to_csv(tables_dir / "Figure5_v7_sv_burden_sample_level.tsv", sep="\t", index=False)
    sv_stats.to_csv(tables_dir / "Figure5_v7_sv_stats.tsv", sep="\t", index=False)
    methyl_bins.to_csv(tables_dir / "Figure5_v7_methylation_distance_bins.tsv", sep="\t", index=False)
    methyl_curve.to_csv(tables_dir / "Figure5_v7_methylation_distance_decay_summary.tsv", sep="\t", index=False)
    methyl_effect.to_csv(tables_dir / "Figure5_v7_methylation_effect_summary.tsv", sep="\t", index=False)

    render_figure(
        figures_dir / "Figure5_v7.png",
        figures_dir / "Figure5_v7.pdf",
        deletion_df,
        coverage_df,
        cnv_df,
        cnv_burden,
        cnv_stats,
        chrom_sizes,
        sv_burden,
        sv_stats,
        methyl_curve,
        methyl_effect,
        args.gtf,
    )
    render_supplementary_coverage_figure(
        figures_dir / "Figure5_v7_supplementary_coverage.png",
        figures_dir / "Figure5_v7_supplementary_coverage.pdf",
        deletion_df,
        coverage_df,
        args.gtf,
    )
    write_report(
        reports_dir / "Figure5_v7_report.md",
        deletion_df,
        cnv_burden,
        cnv_stats,
        sv_stats,
        methyl_effect,
    )


if __name__ == "__main__":
    main()
