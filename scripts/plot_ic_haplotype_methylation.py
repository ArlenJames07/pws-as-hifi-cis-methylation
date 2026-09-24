#!/usr/bin/env python3
"""Plot and audit haplotype methylation across the PWS imprinting center.

The script reads pb-CpG-tools BED files from one directory per sample.  BED
column 4 is treated as percent methylation and column 6 as CpG coverage, which
matches the files produced by this project.  Maternal-like and paternal-like
reference centroids are learned from unaffected controls; HP1/HP2 labels are
never assumed to encode parental origin.

Deletion samples commonly have no hap1/hap2 calls across the deleted/retained
IC.  For those samples the two requested haplotype tracks remain visible as
missing, and the combined BED is plotted and assessed as the retained-copy
fallback.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "results" / "06_methylation"
DEFAULT_METADATA = PROJECT_ROOT / "assets" / "metadata.csv"
DEFAULT_OUTPUT = PROJECT_ROOT / "results" / "08_ic_haplotype_profiles"
DEFAULT_BED_TEMPLATE = "{sample}.cpg.{label}.bed"

DEFAULT_CHROM = "chr15"
DEFAULT_START = 22_691_258
DEFAULT_END = 22_693_494

MIN_STATE_CPGS = 3
MIN_STATE_MEAN_COVERAGE = 2.0
BOOTSTRAP_REPLICATES = 2_000
BOOTSTRAP_SEED = 15_226_912

COLORS = {
    "hap1": "#0072B2",
    "hap2": "#D55E00",
    "combined": "#7A3E9D",
    "maternal": "#B2182B",
    "paternal": "#2166AC",
    "boundary": "#666666",
}


@dataclass
class Track:
    label: str
    path: Path
    positions: np.ndarray
    beta: np.ndarray
    coverage: np.ndarray
    mean_beta: float | None
    mean_coverage: float | None
    ci_low: float | None
    ci_high: float | None

    @property
    def n_cpgs(self) -> int:
        return int(self.positions.size)

    @property
    def estimable(self) -> bool:
        return (
            self.mean_beta is not None
            and self.n_cpgs >= MIN_STATE_CPGS
            and self.mean_coverage is not None
            and self.mean_coverage >= MIN_STATE_MEAN_COVERAGE
        )


@dataclass(frozen=True)
class Reference:
    maternal: float
    paternal: float
    controls: tuple[str, ...]

    @property
    def boundary(self) -> float:
        return (self.maternal + self.paternal) / 2.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create per-sample hap1/hap2 methylation profiles over the PWS IC."
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--bed-template",
        default=DEFAULT_BED_TEMPLATE,
        help=("BED filename template relative to each sample directory; "
              "available fields: {sample}, {label} "
              "(default: %(default)s)"),
    )
    parser.add_argument("--chrom", default=DEFAULT_CHROM)
    parser.add_argument("--start", type=int, default=DEFAULT_START,
                        help="0-based BED start (default: %(default)s)")
    parser.add_argument("--end", type=int, default=DEFAULT_END,
                        help="0-based BED end, exclusive (default: %(default)s)")
    parser.add_argument("--smooth-cpgs", type=int, default=9,
                        help="Centered rolling window in CpGs (default: %(default)s)")
    parser.add_argument("--bootstrap-replicates", type=int,
                        default=BOOTSTRAP_REPLICATES)
    return parser.parse_args()
def load_metadata(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Metadata file not found: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required = {"sample_id", "clinical_diagnosis", "molecular_mechanism"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Metadata must contain {sorted(required)}: {path}")
    return {row["sample_id"]: row for row in rows}


def bed_path(sample_dir: Path, sample: str, label: str, template: str) -> Path:
    try:
        filename = template.format(sample=sample, label=label)
    except (KeyError, ValueError) as error:
        raise ValueError(
            "--bed-template must use only the {sample} and {label} fields"
        ) from error
    return sample_dir / filename


def discover_samples(input_dir: Path, bed_template: str) -> list[str]:
    if not input_dir.is_dir():
        raise NotADirectoryError(f"Input directory not found: {input_dir}")
    samples = []
    for directory in sorted(input_dir.iterdir()):
        if not directory.is_dir() or directory.name == "DSS":
            continue
        sample = directory.name
        if any(bed_path(directory, sample, label, bed_template).exists()
               for label in ("hap1", "hap2", "combined")):
            samples.append(sample)
    if not samples:
        raise RuntimeError(f"No sample BED files found under {input_dir}")
    return samples


def bootstrap_ci(beta: np.ndarray, coverage: np.ndarray, path: Path,
                 start: int, end: int, replicates: int) -> tuple[float | None, float | None]:
    if beta.size < 2 or replicates <= 0:
        return None, None
    digest = hashlib.sha256(f"{path}|{start}|{end}".encode()).hexdigest()
    rng = np.random.default_rng((BOOTSTRAP_SEED + int(digest[:8], 16)) % (2**32 - 1))
    values = np.empty(replicates, dtype=float)
    n = beta.size
    for i in range(replicates):
        index = rng.integers(0, n, size=n)
        b = beta[index]
        w = coverage[index]
        values[i] = float(np.average(b, weights=w)) if w.sum() > 0 else float(b.mean())
    low, high = np.quantile(values, [0.025, 0.975])
    return float(low), float(high)


def read_track(path: Path, label: str, chrom: str, start: int, end: int,
               bootstrap_replicates: int) -> Track:
    # Aggregate duplicate positions by coverage before calculating summaries.
    by_position: dict[int, list[float]] = {}
    if path.exists():
        chrom_prefix = f"{chrom}\t"
        seen_chrom = False
        with path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                # pb-CpG-tools BEDs are coordinate sorted and can be tens of
                # gigabytes across a cohort.  Reject other chromosomes before
                # splitting/parsing columns, then stop after the target block.
                if not line.startswith(chrom_prefix):
                    if seen_chrom:
                        break
                    continue
                seen_chrom = True
                fields = line.split()
                if len(fields) < 6:
                    continue
                try:
                    row_start = int(fields[1])
                    row_end = int(fields[2])
                    beta = float(fields[3]) / 100.0
                    coverage = float(fields[5])
                except ValueError:
                    continue
                if row_start >= end:
                    break
                if row_end <= start:
                    continue
                if not (math.isfinite(beta) and math.isfinite(coverage)) or coverage < 0:
                    continue
                if beta < 0 or beta > 1:
                    raise ValueError(
                        f"Methylation outside [0,1] after percent conversion at "
                        f"{path}:{line_number}"
                    )
                weighted_sum, weight = by_position.setdefault(row_start, [0.0, 0.0])
                by_position[row_start] = [weighted_sum + beta * coverage, weight + coverage]

    positions = np.asarray(sorted(by_position), dtype=int)
    beta_values = np.asarray([
        by_position[pos][0] / by_position[pos][1] if by_position[pos][1] > 0 else np.nan
        for pos in positions
    ], dtype=float)
    coverage_values = np.asarray([by_position[pos][1] for pos in positions], dtype=float)
    valid = np.isfinite(beta_values)
    positions = positions[valid]
    beta_values = beta_values[valid]
    coverage_values = coverage_values[valid]

    if beta_values.size:
        mean_beta = (float(np.average(beta_values, weights=coverage_values))
                     if coverage_values.sum() > 0 else float(beta_values.mean()))
        mean_coverage = float(coverage_values.mean())
        ci_low, ci_high = bootstrap_ci(
            beta_values, coverage_values, path, start, end, bootstrap_replicates
        )
    else:
        mean_beta = mean_coverage = ci_low = ci_high = None

    return Track(label, path, positions, beta_values, coverage_values,
                 mean_beta, mean_coverage, ci_low, ci_high)


def learn_reference(sample_tracks: dict[str, dict[str, Track]],
                    metadata: dict[str, dict[str, str]],
                    exclude: str | None = None) -> Reference:
    highs: list[float] = []
    lows: list[float] = []
    controls: list[str] = []
    for sample, row in metadata.items():
        if row["molecular_mechanism"] != "Control" or sample == exclude:
            continue
        tracks = sample_tracks.get(sample)
        if not tracks:
            continue
        h1, h2 = tracks["hap1"], tracks["hap2"]
        if not h1.estimable or not h2.estimable:
            continue
        low, high = sorted((float(h1.mean_beta), float(h2.mean_beta)))
        lows.append(low)
        highs.append(high)
        controls.append(sample)
    if not highs:
        raise RuntimeError(
            "No unaffected control has two estimable IC haplotypes; "
            "parental-like reference centroids cannot be learned."
        )
    maternal = float(np.median(highs))
    paternal = float(np.median(lows))
    if not maternal > paternal:
        raise RuntimeError("Control-derived maternal centroid must exceed paternal centroid")
    return Reference(maternal, paternal, tuple(controls))


def reference_for_sample(sample: str, mechanism: str,
                         tracks: dict[str, dict[str, Track]],
                         metadata: dict[str, dict[str, str]],
                         global_reference: Reference) -> Reference:
    if mechanism == "Control":
        try:
            return learn_reference(tracks, metadata, exclude=sample)
        except RuntimeError:
            pass
    return global_reference


def classify(track: Track, reference: Reference) -> str:
    if track.mean_beta is None:
        return "missing"
    if not track.estimable:
        return "low-support"
    if track.ci_low is not None and track.ci_high is not None:
        if track.ci_low > reference.boundary:
            return "maternal-like"
        if track.ci_high < reference.boundary:
            return "paternal-like"
        return "uncertain"
    return "maternal-like" if track.mean_beta > reference.boundary else "paternal-like"


def expected_pattern(mechanism: str) -> str:
    return {
        "PWS-DEL": "retained copy maternal-like",
        "AS-DEL": "retained copy paternal-like",
        "PWS-mUPD": "hap1 and hap2 maternal-like",
        "Control": "one maternal-like and one paternal-like (HP order arbitrary)",
        "DiGeorge": "one maternal-like and one paternal-like (HP order arbitrary)",
    }.get(mechanism, "not defined (metadata missing or unsupported mechanism)")


def assess(mechanism: str, calls: dict[str, str]) -> tuple[str, str, str]:
    """Return sample status, hap1 status and track used for sample assessment."""
    h1, h2, combined = calls["hap1"], calls["hap2"], calls["combined"]
    if mechanism in {"PWS-DEL", "AS-DEL"}:
        target = "maternal-like" if mechanism == "PWS-DEL" else "paternal-like"
        if combined in {"missing", "low-support", "uncertain"}:
            sample_status = "not_evaluable"
        else:
            sample_status = "match" if combined == target else "mismatch"
        if h1 in {"missing", "low-support"}:
            hap1_status = "not_evaluable_use_combined"
        else:
            hap1_status = "match" if h1 == target else "mismatch"
        return sample_status, hap1_status, "combined_retained_copy"

    if mechanism == "PWS-mUPD":
        if h1 not in {"maternal-like", "paternal-like"} or h2 not in {
            "maternal-like", "paternal-like"
        }:
            status = "not_evaluable"
        else:
            status = "match" if h1 == h2 == "maternal-like" else "mismatch"
        hap1_status = "match" if h1 == "maternal-like" else (
            "not_evaluable" if h1 in {"missing", "low-support", "uncertain"} else "mismatch"
        )
        return status, hap1_status, "hap1_and_hap2"

    if mechanism in {"Control", "DiGeorge"}:
        resolved = sorted((h1, h2))
        if any(call not in {"maternal-like", "paternal-like"} for call in (h1, h2)):
            status = "not_evaluable"
        else:
            status = "match" if resolved == ["maternal-like", "paternal-like"] else "mismatch"
        if h1 not in {"maternal-like", "paternal-like"}:
            hap1_status = "not_evaluable"
        else:
            hap1_status = "matches_biparental_pair" if status == "match" else "pair_mismatch"
        return status, hap1_status, "hap1_and_hap2"

    return "not_assessed", "not_assessed", "none"


def smooth_track(track: Track, window: int) -> np.ndarray:
    if track.n_cpgs == 0:
        return np.asarray([], dtype=float)
    window = max(1, int(window))
    half = window // 2
    out = np.empty(track.n_cpgs, dtype=float)
    for i in range(track.n_cpgs):
        left = max(0, i - half)
        right = min(track.n_cpgs, i + half + 1)
        b = track.beta[left:right]
        w = track.coverage[left:right]
        out[i] = float(np.average(b, weights=w)) if w.sum() > 0 else float(b.mean())
    return out


def fmt(value: float | None, digits: int = 3) -> str:
    return "" if value is None else f"{value:.{digits}f}"


def legend_label(track: Track, call: str) -> str:
    if track.mean_beta is None:
        return f"{track.label}: no IC CpGs"
    return (f"{track.label}: beta={track.mean_beta:.3f}, n={track.n_cpgs}, "
            f"{call}")


def add_track(ax: plt.Axes, track: Track, call: str, start: int,
              smooth_cpgs: int, color: str, linestyle: str = "-") -> None:
    label = legend_label(track, call)
    if track.n_cpgs == 0:
        ax.plot([], [], color=color, lw=1.8, linestyle=linestyle, label=label)
        return
    x = (track.positions - start) / 1000.0
    ax.scatter(x, track.beta, s=9, color=color, alpha=0.25, linewidths=0)
    ax.plot(x, smooth_track(track, smooth_cpgs), color=color, lw=1.8,
            linestyle=linestyle, label=label)


def plot_sample(sample: str, mechanism: str, tracks: dict[str, Track],
                calls: dict[str, str], reference: Reference, status: str,
                hap1_status: str, chrom: str, start: int, end: int,
                smooth_cpgs: int, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    add_track(ax, tracks["hap1"], calls["hap1"], start, smooth_cpgs, COLORS["hap1"])
    add_track(ax, tracks["hap2"], calls["hap2"], start, smooth_cpgs, COLORS["hap2"])
    if mechanism in {"PWS-DEL", "AS-DEL"} or (
        tracks["hap1"].n_cpgs == 0 and tracks["hap2"].n_cpgs == 0
    ):
        add_track(ax, tracks["combined"], calls["combined"], start, smooth_cpgs,
                  COLORS["combined"], "--")

    ax.axhline(reference.maternal, color=COLORS["maternal"], lw=0.9, ls=":",
               label=f"maternal reference {reference.maternal:.3f}")
    ax.axhline(reference.paternal, color=COLORS["paternal"], lw=0.9, ls=":",
               label=f"paternal reference {reference.paternal:.3f}")
    ax.axhline(reference.boundary, color=COLORS["boundary"], lw=0.8, ls="--",
               label=f"decision boundary {reference.boundary:.3f}")
    ax.set(xlim=(0, (end - start) / 1000.0), ylim=(-0.04, 1.04),
           xlabel=f"Position within PWS IC (kb from {chrom}:{start:,}; BED coordinates)",
           ylabel="CpG methylation beta")
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.grid(axis="y", color="#DDDDDD", lw=0.6)
    clinical = mechanism if mechanism != "Unknown" else "metadata unavailable"
    ax.set_title(f"{sample} | {clinical} | PWS IC {chrom}:{start:,}-{end:,}", loc="left")
    fig.text(0.10, 0.015,
             f"Expected: {expected_pattern(mechanism)} | sample assessment: {status} | "
             f"hap1 assessment: {hap1_status}", fontsize=8.5)
    ax.legend(loc="upper right", fontsize=7.5, frameon=False, ncol=2)
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_cohort(samples: list[str], metadata: dict[str, dict[str, str]],
                all_tracks: dict[str, dict[str, Track]],
                all_calls: dict[str, dict[str, str]], statuses: dict[str, str],
                references: dict[str, Reference], start: int, end: int,
                smooth_cpgs: int, output_stem: Path) -> None:
    ncols = 4
    nrows = math.ceil(len(samples) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 2.65 * nrows),
                             sharex=True, sharey=True)
    axes_flat = np.atleast_1d(axes).ravel()
    for ax, sample in zip(axes_flat, samples):
        mechanism = metadata.get(sample, {}).get("molecular_mechanism", "Unknown")
        tracks, calls, reference = all_tracks[sample], all_calls[sample], references[sample]
        for label in ("hap1", "hap2"):
            track = tracks[label]
            if track.n_cpgs:
                x = (track.positions - start) / 1000.0
                ax.plot(x, smooth_track(track, smooth_cpgs), lw=1.2,
                        color=COLORS[label], label=label)
        if mechanism in {"PWS-DEL", "AS-DEL"}:
            track = tracks["combined"]
            if track.n_cpgs:
                x = (track.positions - start) / 1000.0
                ax.plot(x, smooth_track(track, smooth_cpgs), lw=1.2, ls="--",
                        color=COLORS["combined"], label="combined")
        ax.axhline(reference.boundary, color="#888888", lw=0.6, ls=":")
        ax.set_title(f"{sample} | {mechanism} | {statuses[sample]}", fontsize=8, loc="left")
        ax.grid(axis="y", color="#E7E7E7", lw=0.4)
    for ax in axes_flat[len(samples):]:
        ax.set_visible(False)
    for ax in axes_flat:
        ax.set_xlim(0, (end - start) / 1000.0)
        ax.set_ylim(-0.04, 1.04)
    fig.supxlabel(f"Position within PWS IC (kb from chr15:{start:,})", fontsize=10)
    fig.supylabel("CpG methylation beta", fontsize=10)
    handles = [
        plt.Line2D([], [], color=COLORS["hap1"], label="hap1"),
        plt.Line2D([], [], color=COLORS["hap2"], label="hap2"),
        plt.Line2D([], [], color=COLORS["combined"], ls="--",
                   label="combined retained-copy fallback"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 1.005))
    fig.suptitle("Haplotype methylation across the PWS imprinting center",
                 y=1.025, fontsize=14)
    fig.tight_layout()
    # Path.with_suffix() would incorrectly discard the final ".hap1_hap2"
    # portion of this deliberately dotted output stem.
    fig.savefig(Path(f"{output_stem}.png"), dpi=180, bbox_inches="tight")
    fig.savefig(Path(f"{output_stem}.pdf"), bbox_inches="tight")
    plt.close(fig)


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    if args.end <= args.start:
        raise ValueError("--end must be greater than --start")
    metadata = load_metadata(args.metadata)
    samples = discover_samples(args.input_dir, args.bed_template)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    per_sample_dir = args.output_dir / "per_sample"
    per_sample_dir.mkdir(exist_ok=True)

    all_tracks: dict[str, dict[str, Track]] = {}
    for sample in samples:
        sample_dir = args.input_dir / sample
        all_tracks[sample] = {
            label: read_track(bed_path(sample_dir, sample, label, args.bed_template), label,
                              args.chrom, args.start, args.end,
                              args.bootstrap_replicates)
            for label in ("hap1", "hap2", "combined")
        }

    global_reference = learn_reference(all_tracks, metadata)
    all_calls: dict[str, dict[str, str]] = {}
    references: dict[str, Reference] = {}
    statuses: dict[str, str] = {}
    summary_rows: list[dict[str, object]] = []
    profile_rows: list[dict[str, object]] = []

    for sample in samples:
        row = metadata.get(sample, {})
        mechanism = row.get("molecular_mechanism", "Unknown")
        reference = reference_for_sample(
            sample, mechanism, all_tracks, metadata, global_reference
        )
        references[sample] = reference
        tracks = all_tracks[sample]
        calls = {label: classify(track, reference) for label, track in tracks.items()}
        all_calls[sample] = calls
        status, hap1_status, assessment_track = assess(mechanism, calls)
        statuses[sample] = status

        summary: dict[str, object] = {
            "sample_id": sample,
            "clinical_diagnosis": row.get("clinical_diagnosis", "metadata unavailable"),
            "molecular_mechanism": mechanism,
            "expected_pattern": expected_pattern(mechanism),
            "sample_pattern_status": status,
            "hap1_expected_pattern_status": hap1_status,
            "assessment_track": assessment_track,
            "maternal_reference_beta": fmt(reference.maternal),
            "paternal_reference_beta": fmt(reference.paternal),
            "decision_boundary_beta": fmt(reference.boundary),
            "reference_controls": ",".join(reference.controls),
        }
        for label, track in tracks.items():
            summary.update({
                f"{label}_n_cpgs": track.n_cpgs,
                f"{label}_mean_beta": fmt(track.mean_beta),
                f"{label}_ci_low": fmt(track.ci_low),
                f"{label}_ci_high": fmt(track.ci_high),
                f"{label}_mean_coverage": fmt(track.mean_coverage),
                f"{label}_pattern": calls[label],
            })
            for position, beta, coverage in zip(
                track.positions, track.beta, track.coverage, strict=True
            ):
                profile_rows.append({
                    "sample_id": sample,
                    "molecular_mechanism": mechanism,
                    "track": label,
                    "chrom": args.chrom,
                    "position_0based": int(position),
                    "position_1based": int(position) + 1,
                    "beta": f"{float(beta):.6f}",
                    "coverage": f"{float(coverage):.1f}",
                })
        summary_rows.append(summary)
        plot_sample(sample, mechanism, tracks, calls, reference, status,
                    hap1_status, args.chrom, args.start, args.end,
                    args.smooth_cpgs,
                    per_sample_dir / f"{sample}.PWS_IC.hap1_hap2.png")

    summary_fields = list(summary_rows[0])
    write_tsv(args.output_dir / "ic_haplotype_pattern_summary.tsv",
              summary_rows, summary_fields)
    profile_fields = ["sample_id", "molecular_mechanism", "track", "chrom",
                      "position_0based", "position_1based", "beta", "coverage"]
    profile_path = args.output_dir / "ic_cpg_profile_source_data.tsv.gz"
    with gzip.open(profile_path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=profile_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(profile_rows)

    plot_cohort(samples, metadata, all_tracks, all_calls, statuses, references,
                args.start, args.end, args.smooth_cpgs,
                args.output_dir / "all_samples.PWS_IC.hap1_hap2")

    counts = {name: sum(row["sample_pattern_status"] == name for row in summary_rows)
              for name in ("match", "mismatch", "not_evaluable", "not_assessed")}
    missing_metadata = [sample for sample in samples if sample not in metadata]
    readme = f"""# PWS IC haplotype methylation profiles

Generated by `scripts/plot_ic_haplotype_methylation.py` from
`{args.input_dir}`.

- Interval: `{args.chrom}:{args.start}-{args.end}` (0-based, half-open BED coordinates)
- BED filename template: `{args.bed_template}`
- Samples plotted: {len(samples)}
- Expected-pattern matches: {counts['match']}
- Expected-pattern mismatches: {counts['mismatch']}
- Not evaluable: {counts['not_evaluable']}
- Not assessed because the mechanism is unknown: {counts['not_assessed']}
- Maternal-like control centroid: {global_reference.maternal:.3f}
- Paternal-like control centroid: {global_reference.paternal:.3f}
- Global decision boundary: {global_reference.boundary:.3f}

HP1 and HP2 are technical phase-block labels, not parental labels. For Control
and DiGeorge samples, either HP order is expected as long as the pair contains
one maternal-like and one paternal-like state. PWS-mUPD expects two
maternal-like states. PWS-DEL and AS-DEL haplotype BEDs have no IC CpGs in this
dataset, so their retained-copy state is evaluated using the combined BED and
shown as a dashed purple fallback.

Missing metadata: {', '.join(missing_metadata) if missing_metadata else 'none'}.
These samples are plotted and classified relative to the control centroids, but
no expected biological pattern is assigned.
"""
    (args.output_dir / "README.md").write_text(readme, encoding="utf-8")

    print(f"Wrote {len(samples)} per-sample plots to {per_sample_dir}")
    print(f"Summary: {args.output_dir / 'ic_haplotype_pattern_summary.tsv'}")
    print("Status counts: " + ", ".join(f"{key}={value}" for key, value in counts.items()))


if __name__ == "__main__":
    main()
