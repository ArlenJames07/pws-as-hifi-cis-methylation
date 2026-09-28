#!/usr/bin/env python3
"""
07 -- Parent-of-origin methylation of repeat elements across 15q11-q13.

The question is epigenetic: a repeat's SEQUENCE is the same whichever parent
transmitted the chromosome, but its methylation can differ. The analysis
follows the evidence rules of scripts/analysis/README.md:

  direct (primary)  PWS-DEL retained (maternal) minus AS-DEL retained (paternal)
                    methylation per element, only for elements that lie fully
                    inside the common reciprocal CN=1 core. Participants are the
                    unit; participant-bootstrap intervals.
  phase-invariant   |H1 - H2| per element in Control and DiGeorge genomes
                    (unoriented allele-specific methylation, ASM). No maternal or
                    paternal label is assigned to H1/H2.
  mUPD check        |H1 - H2| in PWS-mUPD: both chromosomes are maternal, so
                    elements with a direct parental difference should show little
                    haplotype difference there.
  validation        do elements with consistent ASM in biparental genomes carry a
                    direct parental difference, more often than under a
                    circular-shift null that keeps the spatial clustering?

Optional sensitivity analyses (off by default; outputs are labelled "inferred",
never "observed"):
  IC-anchored       H1/H2 oriented by IC methylation, only inside the HiPhase
                    block that contains the IC.
  contig-anchored   maternal-like/paternal-like tracks from the assembly contig
                    carrying the IC (04_assembly_methylation.py).

Outputs: results/08_duplicons/repeat_methylation/
  element_by_participant.tsv.gz, element_summary.tsv, validation.tsv,
  class_summary.tsv, repeat_methylation.{png,pdf}, analysis_manifest.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
from duplicon_analysis import cis_analysis  # noqa: E402
from duplicon_analysis.assembly import read_bed_intervals  # noqa: E402
from duplicon_analysis.methylation import element_means as _element_means  # noqa: E402
from duplicon_analysis.methylation import find_track, read_track  # noqa: E402
from duplicon_analysis.repeats import read_repeat_annotation  # noqa: E402

METADATA_PATH = C.METADATA_PATH
METHYLATION_DIR = C.METHYLATION_DIR
PHASING_DIR = C.PHASING_DIR
ASSEMBLY_DIR = C.ASSEMBLY_DIR
STRUCTURAL_EVIDENCE_PATH = C.STRUCTURAL_EVIDENCE_PATH
# RepeatMasker on the CHM13 window from 01 (class/family included). Any RepeatMasker
# .out, UCSC rmsk.txt or BED with classes also works.
REPEATS_PATH = C.REFERENCE_DIR / "window.fa.out"
SD_BED_PATH = C.REFERENCE_DIR / "window.sd.bed"
OUTPUT_DIR = C.REPEAT_METHYLATION_DIR

CHROM = C.CHROM
DOMAIN_START, DOMAIN_END = C.DOMAIN_START, C.DOMAIN_END    # same domain as analysis 01-03
PWS_IC_START, PWS_IC_END = C.PWS_IC_START, C.PWS_IC_END
IC_EXCLUDE_BP = 5_000              # IC-proximal elements are excluded from tests
CN1_BREAKPOINT_BUFFER = C.CN1_BREAKPOINT_BUFFER
MIN_COVERAGE = C.METHYLATION_MIN_COVERAGE
MIN_CPGS = 3
EFFECT = 0.20                      # |delta beta| that defines a parent-of-origin element
ASM_CONSISTENCY = 0.80             # fraction of biparental genomes that must reach EFFECT
MIN_BIPARENTAL = 3
BOOTSTRAP_REPLICATES = 2_000
CIRCULAR_SHIFTS = 2_000
RANDOM_SEED = 20260925
RUN_IC_ANCHORED_SENSITIVITY = False
RUN_CONTIG_ANCHORED_SENSITIVITY = False

BIPARENTAL = ("Control", "DiGeorge")


# ----------------------------------------------------------------- helpers
def element_means(position, beta, starts, ends, lo, hi):
    return _element_means(position, beta, starts, ends, lo, hi, MIN_CPGS)


def load(sample, kind):
    path = find_track(METHYLATION_DIR, sample, kind)
    if path is None:
        return None
    track = read_track(path, CHROM, DOMAIN_START, DOMAIN_END)
    keep = track.coverage >= MIN_COVERAGE
    return track.position[keep], track.beta[keep]


def ic_block(sample):
    """HiPhase block containing the IC for this sample, or None."""
    path = PHASING_DIR / f"{sample}.blocks.tsv"
    if not path.is_file():
        return None
    table = pd.read_csv(path, sep="\t")
    table.columns = [c.lstrip("#") for c in table.columns]
    table = table[table["chrom"] == CHROM]
    hit = table[(table["start"] <= PWS_IC_START) & (table["end"] >= PWS_IC_END)]
    return (int(hit.iloc[0]["start"]), int(hit.iloc[0]["end"])) if len(hit) else None


def bootstrap_difference(mat, pat, rng):
    """Participant bootstrap of median(maternal-retained) - median(paternal-retained)."""
    if len(mat) == 0 or len(pat) == 0:
        return np.nan, np.nan
    m = np.median(mat[rng.integers(0, len(mat), (BOOTSTRAP_REPLICATES, len(mat)))], axis=1)
    p = np.median(pat[rng.integers(0, len(pat), (BOOTSTRAP_REPLICATES, len(pat)))], axis=1)
    return tuple(np.percentile(m - p, [2.5, 97.5]))


# --------------------------------------------------------------------- run
def run() -> dict[str, Path]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED)
    lib = cis_analysis()
    cohort = lib.load_cohort(METADATA_PATH)
    deletion_samples = cohort.samples_for("PWS-DEL") + cohort.samples_for("AS-DEL")
    deletion_map = lib.load_deletion_map(STRUCTURAL_EVIDENCE_PATH, deletion_samples, CN1_BREAKPOINT_BUFFER)
    core = (deletion_map.common_start, deletion_map.common_end)

    el = read_repeat_annotation(REPEATS_PATH, CHROM, DOMAIN_START, DOMAIN_END)
    el.insert(0, "element_id", [f"{CHROM}:{a}-{b}:{n}" for a, b, n in zip(el["start"], el["end"], el["name"])])
    starts, ends = el["start"].to_numpy(), el["end"].to_numpy()
    el["in_sd"] = False
    for s, e, _ in read_bed_intervals(SD_BED_PATH, CHROM):
        el.loc[(el["start"] < e) & (el["end"] > s), "in_sd"] = True
    el["near_ic"] = (el["end"] > PWS_IC_START - IC_EXCLUDE_BP) & (el["start"] < PWS_IC_END + IC_EXCLUDE_BP)
    el["in_common_core"] = (el["start"] >= core[0]) & (el["end"] <= core[1])
    el["testable"] = ~el["in_sd"] & ~el["near_ic"]

    rows, notes = [], []
    for sample in cohort.samples:
        mech = cohort.mechanism(sample)
        if mech in ("PWS-DEL", "AS-DEL"):
            track = load(sample, "combined")
            interval = deletion_map.interval(sample)
            if track is None or interval is None:
                notes.append(f"{sample}: missing combined track or CN=1 interval")
                continue
            v, n = element_means(*track, starts, ends, *interval)
            rows.append(pd.DataFrame({"element_id": el["element_id"], "sample_id": sample,
                                      "mechanism": mech, "measure": "retained_combined",
                                      "value": v, "n_cpg": n}))
        elif mech in BIPARENTAL or mech == "PWS-mUPD":
            h1, h2 = load(sample, "hap1"), load(sample, "hap2")
            if h1 is None or h2 is None:
                notes.append(f"{sample}: missing hap1/hap2 track")
                continue
            v1, n1 = element_means(*h1, starts, ends, DOMAIN_START, DOMAIN_END)
            v2, n2 = element_means(*h2, starts, ends, DOMAIN_START, DOMAIN_END)
            rows.append(pd.DataFrame({"element_id": el["element_id"], "sample_id": sample,
                                      "mechanism": mech, "measure": "abs_h1_minus_h2",
                                      "value": np.abs(v1 - v2), "n_cpg": np.minimum(n1, n2)}))
            comb = load(sample, "combined")
            if comb is not None and mech in BIPARENTAL:
                vc, nc = element_means(*comb, starts, ends, DOMAIN_START, DOMAIN_END)
                rows.append(pd.DataFrame({"element_id": el["element_id"], "sample_id": sample,
                                          "mechanism": mech, "measure": "scaffold_combined",
                                          "value": vc, "n_cpg": nc}))
    long = pd.concat(rows, ignore_index=True)
    long.to_csv(OUTPUT_DIR / "element_by_participant.tsv.gz", sep="\t", index=False,
                float_format="%.4f", na_rep="NA")

    # ---- per element
    def pivot(measure, mech):
        sub = long[(long["measure"] == measure) & long["mechanism"].isin(mech)]
        return sub.pivot_table(index="element_id", columns="sample_id", values="value", aggfunc="first")

    mat = pivot("retained_combined", ["PWS-DEL"]).reindex(el["element_id"])
    pat = pivot("retained_combined", ["AS-DEL"]).reindex(el["element_id"])
    asm = pivot("abs_h1_minus_h2", list(BIPARENTAL)).reindex(el["element_id"])
    mupd = pivot("abs_h1_minus_h2", ["PWS-mUPD"]).reindex(el["element_id"])
    scaf = pivot("scaffold_combined", list(BIPARENTAL)).reindex(el["element_id"])

    S = el.set_index("element_id").copy()
    S["n_pws"] = mat.notna().sum(axis=1)
    S["n_as"] = pat.notna().sum(axis=1)
    S["maternal_retained_median"] = mat.median(axis=1)
    S["paternal_retained_median"] = pat.median(axis=1)
    S["delta_direct_raw"] = S["maternal_retained_median"] - S["paternal_retained_median"]
    both = S["in_common_core"] & (S["n_pws"] >= 1) & (S["n_as"] >= 1)
    offset = float(S.loc[both & S["testable"], "delta_direct_raw"].median()) if (both & S["testable"]).any() else 0.0
    S["delta_direct"] = np.where(both, S["delta_direct_raw"] - offset, np.nan)
    S["complete_separation"] = both & ((mat.min(axis=1) > pat.max(axis=1)) | (mat.max(axis=1) < pat.min(axis=1)))
    boot_mask = (both & S["testable"]).to_numpy()
    mat_values, pat_values = mat.to_numpy(), pat.to_numpy()
    ci = [bootstrap_difference(mat_values[j][np.isfinite(mat_values[j])],
                               pat_values[j][np.isfinite(pat_values[j])], rng)
          if boot_mask[j] else (np.nan, np.nan) for j in range(len(S))]
    S["delta_direct_ci_low"] = [c[0] - offset if np.isfinite(c[0]) else np.nan for c in ci]
    S["delta_direct_ci_high"] = [c[1] - offset if np.isfinite(c[1]) else np.nan for c in ci]
    S["n_biparental"] = asm.notna().sum(axis=1)
    S["asm_mean"] = asm.mean(axis=1)
    S["asm_fraction_ge_effect"] = (asm >= EFFECT).sum(axis=1) / S["n_biparental"].replace(0, np.nan)
    S["scaffold_mean"] = scaf.mean(axis=1)
    S["mupd_abs_h1_minus_h2"] = mupd.mean(axis=1)
    S["direct_parent_of_origin"] = S["testable"] & S["complete_separation"] & (S["delta_direct"].abs() >= EFFECT)
    S["asm_candidate"] = (S["testable"] & (S["n_biparental"] >= MIN_BIPARENTAL)
                          & (S["asm_mean"] >= EFFECT) & (S["asm_fraction_ge_effect"] >= ASM_CONSISTENCY))
    S["direction"] = np.where(S["direct_parent_of_origin"],
                              np.where(S["delta_direct"] > 0, "maternal_higher", "paternal_higher"),
                              np.where(S["asm_candidate"] & ~S["in_common_core"],
                                       "asm_unoriented_outside_core", ""))
    S.reset_index().to_csv(OUTPUT_DIR / "element_summary.tsv", sep="\t", index=False,
                           float_format="%.4g", na_rep="NA")

    # ---- validation between designs, inside the common core
    V = S[S["testable"] & both & (S["n_biparental"] >= MIN_BIPARENTAL)].sort_values("start")
    cand = V["asm_candidate"].to_numpy()
    direct = (V["complete_separation"] & (V["delta_direct"].abs() >= EFFECT)).to_numpy()

    def agreement(direct_flags):
        return float(direct_flags[cand].mean()) if cand.any() else np.nan

    observed = agreement(direct)
    null = np.array([agreement(np.roll(direct, rng.integers(10, len(V) - 10)))
                     for _ in range(CIRCULAR_SHIFTS)]) if cand.any() and len(V) > 20 else np.array([])
    p_value = (1 + np.sum(null >= observed)) / (1 + null.size) if null.size else np.nan
    rho = V["asm_mean"].corr(V["delta_direct"].abs(), method="spearman") if len(V) > 5 else np.nan
    dpo = V[V["direct_parent_of_origin"]]
    validation = pd.DataFrame([
        ("parental offset removed (PWS-DEL minus AS-DEL, median over elements)", offset),
        ("elements testable in both designs (common core, outside SDs, >5 kb from IC)", len(V)),
        ("biparental ASM candidates", int(cand.sum())),
        ("... of which with a direct parent-of-origin difference", int((cand & direct).sum())),
        ("agreement (fraction of ASM candidates with direct difference)", observed),
        ("circular-shift null median agreement", float(np.median(null)) if null.size else np.nan),
        ("circular-shift p", p_value),
        ("direct parent-of-origin elements", len(dpo)),
        ("... of which ASM candidates in biparental genomes", int(dpo["asm_candidate"].sum())),
        ("median |H1-H2| in PWS-mUPD at direct elements (expected low)", float(dpo["mupd_abs_h1_minus_h2"].median()) if len(dpo) else np.nan),
        ("median |H1-H2| in biparental genomes at direct elements (expected high)", float(dpo["asm_mean"].median()) if len(dpo) else np.nan),
        ("Spearman rho(ASM, |direct delta|)", rho),
    ], columns=["quantity", "value"])
    validation.to_csv(OUTPUT_DIR / "validation.tsv", sep="\t", index=False, float_format="%.4g")

    T = S[S["testable"]].copy()
    cls = (T.groupby("class")
           .agg(elements=("start", "size"),
                in_common_core=("in_common_core", "sum"),
                direct_parent_of_origin=("direct_parent_of_origin", "sum"),
                asm_candidates=("asm_candidate", "sum"),
                median_asm=("asm_mean", "median"),
                median_abs_direct_delta=("delta_direct", lambda x: float(np.nanmedian(np.abs(x))) if x.notna().any() else np.nan))
           .reset_index())
    cls.to_csv(OUTPUT_DIR / "class_summary.tsv", sep="\t", index=False, float_format="%.4g")

    if RUN_IC_ANCHORED_SENSITIVITY:
        ic_anchored_sensitivity(S, starts, ends)
    if RUN_CONTIG_ANCHORED_SENSITIVITY:
        contig_anchored_sensitivity(S, starts, ends)

    plot(S, V, core, OUTPUT_DIR, observed, p_value)
    manifest = {k: v for k, v in globals().items() if k.isupper() and isinstance(v, (int, float, str, bool, tuple))}
    manifest.update({"common_core": list(core), "notes": notes})
    (OUTPUT_DIR / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2, default=str))
    for note in notes:
        print(f"[07] {note}", flush=True)
    print(f"[07] {int(S['direct_parent_of_origin'].sum())} direct parent-of-origin elements; "
          f"ASM agreement {observed:.2f} (circular-shift p = {p_value:.3g})", flush=True)
    return {"summary": OUTPUT_DIR / "element_summary.tsv", "validation": OUTPUT_DIR / "validation.tsv"}


# -------------------------------------------------- optional, inferred only
def _signed_table(rows, name):
    if not rows:
        return
    table = pd.concat(rows, ignore_index=True)
    table.to_csv(OUTPUT_DIR / f"{name}.tsv.gz", sep="\t", index=False, float_format="%.4f", na_rep="NA")


def ic_anchored_sensitivity(S, starts, ends):
    """Orient H1/H2 by IC methylation inside the HiPhase block containing the IC.
    Labels are inferred ('maternal_like'), valid only inside that block."""
    cohort = cis_analysis().load_cohort(METADATA_PATH)
    rows = []
    for sample in cohort.samples:
        if cohort.mechanism(sample) not in BIPARENTAL:
            continue
        block, h1, h2 = ic_block(sample), load(sample, "hap1"), load(sample, "hap2")
        if block is None or h1 is None or h2 is None:
            continue
        ic1 = h1[1][(h1[0] >= PWS_IC_START) & (h1[0] < PWS_IC_END)].mean()
        ic2 = h2[1][(h2[0] >= PWS_IC_START) & (h2[0] < PWS_IC_END)].mean()
        if not (np.isfinite(ic1) and np.isfinite(ic2)) or abs(ic1 - ic2) < 0.3:
            continue
        m, p = (h1, h2) if ic1 > ic2 else (h2, h1)
        vm, _ = element_means(*m, starts, ends, *block)
        vp, _ = element_means(*p, starts, ends, *block)
        rows.append(pd.DataFrame({"element_id": S.index, "sample_id": sample, "scope_start": block[0],
                                  "scope_end": block[1], "maternal_like_minus_paternal_like": vm - vp}))
    _signed_table(rows, "ic_anchored_sensitivity_inferred")


def contig_anchored_sensitivity(S, starts, ends):
    cohort = cis_analysis().load_cohort(METADATA_PATH)
    rows = []
    for sample in cohort.samples:
        if cohort.mechanism(sample) not in BIPARENTAL:
            continue
        paths = [ASSEMBLY_DIR / sample / f"{sample}.{lab}.chm13.bed" for lab in ("maternal_like", "paternal_like")]
        if not all(p.is_file() and p.stat().st_size for p in paths):
            continue
        tracks = [read_track(p, CHROM, DOMAIN_START, DOMAIN_END) for p in paths]
        vm, _ = element_means(tracks[0].position, tracks[0].beta, starts, ends, DOMAIN_START, DOMAIN_END)
        vp, _ = element_means(tracks[1].position, tracks[1].beta, starts, ends, DOMAIN_START, DOMAIN_END)
        rows.append(pd.DataFrame({"element_id": S.index, "sample_id": sample,
                                  "maternal_like_minus_paternal_like": vm - vp}))
    _signed_table(rows, "contig_anchored_sensitivity_inferred")


# -------------------------------------------------------------------- plot
def plot(S, V, core, outdir, observed, p_value):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ink, ink2, grid, base, hi, direct_c = "#0b0b0b", "#52514e", "#e4e3df", "#b9b8b3", "#eb6834", "#2a78d6"
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1, 1.6]})
    other = V[~V["direct_parent_of_origin"]]
    call = V[V["direct_parent_of_origin"]]
    a1.axhline(0, color=grid, lw=0.8)
    for y in (-EFFECT, EFFECT):
        a1.axhline(y, color=grid, lw=0.6, ls=(0, (2, 2)))
    a1.axvline(EFFECT, color=grid, lw=0.6, ls=(0, (2, 2)))
    a1.scatter(other["asm_mean"], other["delta_direct"], s=8, lw=0, color=base, label="other elements")
    a1.scatter(call["asm_mean"], call["delta_direct"], s=16, lw=0, color=hi, label="direct parent-of-origin")
    a1.set_xlim(0, 1); a1.set_ylim(-1, 1)
    a1.set_xlabel("biparental genomes: mean |H1 − H2| (unoriented)", fontsize=8, color=ink)
    a1.set_ylabel("PWS-DEL − AS-DEL retained (Δβ, centred)", fontsize=8, color=ink)
    title = (f"ASM candidates with a direct difference: {observed:.0%} (circular-shift p = {p_value:.3g})"
             if np.isfinite(observed) else "No ASM candidates in the common core")
    a1.set_title(title, fontsize=8, color=ink)
    a1.legend(fontsize=7, frameon=False, loc="lower right")

    Y = S[S["testable"]]
    a2.axvspan(core[0] / 1e6, core[1] / 1e6, color="#f3f2ee", lw=0, zorder=0)
    a2.axhline(0, color=grid, lw=0.8)
    a2.scatter(Y["start"] / 1e6, Y["asm_mean"], s=4, lw=0, color=base, label="|H1 − H2|, biparental (unoriented)")
    d = Y[Y["delta_direct"].notna()]
    a2.scatter(d["start"] / 1e6, d["delta_direct"], s=4, lw=0, color=direct_c, alpha=0.7,
               label="PWS-DEL − AS-DEL (common CN=1 core)")
    c = Y[Y["direct_parent_of_origin"]]
    a2.scatter(c["start"] / 1e6, c["delta_direct"], s=14, lw=0, color=hi, label="direct parent-of-origin")
    a2.axvline((PWS_IC_START + PWS_IC_END) / 2e6, color=ink2, lw=0.6, ls=(0, (2, 2)))
    a2.text((PWS_IC_START + PWS_IC_END) / 2e6, 1.01, "IC", ha="center", va="bottom", fontsize=7,
            color=ink, transform=a2.get_xaxis_transform())
    a2.set_xlim(DOMAIN_START / 1e6, DOMAIN_END / 1e6); a2.set_ylim(-1, 1)
    a2.set_xlabel("T2T-CHM13 chr15 (Mb); shaded = common reciprocal CN=1 core", fontsize=8, color=ink)
    a2.set_ylabel("Δβ", fontsize=8, color=ink)
    a2.legend(fontsize=7, frameon=False, loc="lower left", markerscale=2)
    for ax in (a1, a2):
        ax.tick_params(labelsize=7, colors=ink2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(grid)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(outdir / f"repeat_methylation.{ext}", dpi=200)
    plt.close(fig)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ic-anchored", action="store_true", help="also write the IC-anchored sensitivity table")
    ap.add_argument("--contig-anchored", action="store_true", help="also use the 04 contig labels (sensitivity)")
    a = ap.parse_args(argv)
    global RUN_IC_ANCHORED_SENSITIVITY, RUN_CONTIG_ANCHORED_SENSITIVITY
    RUN_IC_ANCHORED_SENSITIVITY = RUN_IC_ANCHORED_SENSITIVITY or a.ic_anchored
    RUN_CONTIG_ANCHORED_SENSITIVITY = RUN_CONTIG_ANCHORED_SENSITIVITY or a.contig_anchored
    run()


if __name__ == "__main__":
    main()
