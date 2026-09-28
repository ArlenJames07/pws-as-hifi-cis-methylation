#!/usr/bin/env python3
"""
08 -- Follow-up checks and one status table: what is confirmed, what is still missing.

Reads the outputs of 03, 05/05b and 07 (nothing is recomputed) and adds three checks:
  * split reads: reads with one alignment at a deletion edge and a supplementary alignment
    (SA tag) at the other edge -> the junction at bp level when one side is unique
    sequence. Counted in each carrier and, as a control, in every genome with two
    copies of chr15 at the same positions.
  * assembled fusions: hifiasm contigs that jump across the deletion (03), with the last
    and first confidently placed blocks on each side.
    Each fusion is tested against the directly oriented SD pairs: a switch between the two
    copies within 20 kb of homologous positions is NAHR confirmed by the assembly.
  * genes: genes of the CHM13 GTF (params.local.yml 'gtf') in every 250-kb bin of 07's
    regional table, and at the ends of the regions with a significant parental difference.
  * parental contrasts: per bin, PWS-DEL, AS-DEL and PWS-mUPD minus the biparental genomes
    (combined methylation). PWS-mUPD was not used to find the regions, so a shift of the
    same sign as PWS-DEL inside a region, and not outside it, validates the region
    (needs 07 from v6 on, which stores the PWS-mUPD combined track).

Output: results/08_duplicons/followup/
  breakpoint_status.tsv    one row per deletion carrier: edges, SD pair, each line of
                           evidence, verdict ('status') and what is missing ('missing')
  split_reads.tsv          carrier edge pair x genome: split reads and the junction
  assembled_fusions.tsv    contigs across a deletion
  methylation_regions.tsv  250-kb bins: raw and centred difference, permutation p,
                           the three contrasts against biparental genomes, genes
  methylation_domains.tsv  runs of significant bins: genes at their two ends, contrasts
                           inside vs outside
  followup_report.md       the same, in words

Usage: python3 scripts/duplicons/08_followup.py [--skip-split-reads] [--window 50000] [--gtf FILE]
  --skip-split-reads reuses followup/split_reads.tsv from an earlier run when it exists.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
from duplicon_analysis import bams as bam_lookup  # noqa: E402
from duplicon_analysis import cis_analysis, followup  # noqa: E402
from duplicon_analysis.assembly import parse_region  # noqa: E402
from duplicon_analysis.junction_reads import read_direct_pairs  # noqa: E402
from duplicon_analysis.tools import log  # noqa: E402

OUT = C.OUT / "followup"


def _read(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def _num(x):
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else x


def edges_by_sample(bp: pd.DataFrame) -> dict:
    out = {}
    for r in bp.itertuples(index=False):
        pos = _num(getattr(r, "sunk_pos", None))
        out.setdefault(r.sample_id, {})[r.edge] = {
            "pos": int(pos) if pos is not None else None,
            "hificnv": int(r.hificnv_pos),
            "cluster": getattr(r, "bp_cluster", "."),
            "note": getattr(r, "note", "") if isinstance(getattr(r, "note", ""), str) else ""}
    return out


def run_split_reads(edges: dict, cohort, window: int) -> pd.DataFrame:
    table = bam_lookup.read_bam_table(C.BAMS_PATH)
    bams, rows = {}, []
    for s in cohort.samples:
        try:
            found = bam_lookup.find_bams(s, C.ALIGNMENT_DIR, table)
        except bam_lookup.BamLookupError as error:
            log(f"split reads: {error}")
            continue
        found = [b for b in found if bam_lookup.reference_check(C.SAMTOOLS, b, C.REFERENCE, C.CHROM)[0]
                 and bam_lookup.bam_index(b)]
        if found:
            bams[s] = found
    panel = [s for s in bams if cohort.mechanism(s) in C.MECHANISMS_DIPLOID_CHR15]
    for carrier, e in edges.items():
        a = e.get("proximal", {})
        b = e.get("distal", {})
        pa = a.get("pos") or a.get("hificnv")
        pb = b.get("pos") or b.get("hificnv")
        if pa is None or pb is None:
            continue
        for genome in [carrier] + panel:
            if genome not in bams:
                continue
            df = followup.split_reads(C.SAMTOOLS, bams[genome], C.CHROM, int(pa), int(pb), window)
            rows.append({"carrier": carrier, "genome": genome, "genome_group": cohort.mechanism(genome),
                         "edge_proximal_used": int(pa), "edge_distal_used": int(pb),
                         **followup.summarise_split(df)})
        log(f"split reads: {carrier} done")
    return pd.DataFrame(rows)


def annotate_fusions(fusions: pd.DataFrame) -> pd.DataFrame:
    """Pair and distance to homologous positions for every assembled junction."""
    path = C.REFERENCE_DIR / "window.sd_pairs.bedpe"
    if fusions.empty or not path.is_file():
        return fusions
    pairs = read_direct_pairs(path)
    got = [followup.fusion_nahr(int(r.ref_left), int(r.ref_right), pairs) for r in fusions.itertuples(index=False)]
    fusions = fusions.copy()
    fusions["fusion_pair"] = [g[0] if g else "." for g in got]
    fusions["fusion_homolog_of_left_in_copy_b"] = [g[1] if g else np.nan for g in got]
    fusions["fusion_homolog_distance_bp"] = [g[2] if g else np.nan for g in got]
    return fusions


def methylation_tables(gtf: Path | None):
    reg = _read(C.REPEAT_METHYLATION_DIR / "direct_delta_by_region.tsv")
    val = _read(C.REPEAT_METHYLATION_DIR / "validation.tsv")
    if reg.empty:
        return pd.DataFrame(), None
    offset = 0.0
    if not val.empty:
        hit = val[val["quantity"].str.startswith("parental offset removed")]
        if len(hit):
            offset = float(hit.iloc[0]["value"])
    cols = ["start", "end", "size_mb", "direction", "bins", "median_raw_delta"]
    reg = reg.dropna(subset=["observed_median_delta"]).copy()
    if reg.empty:
        return reg, pd.DataFrame(columns=cols)
    reg["raw_median_delta"] = reg["observed_median_delta"] + offset
    bin_bp = int(reg["bin_end"].iloc[0] - reg["bin_start"].iloc[0])
    ebp_path = C.REPEAT_METHYLATION_DIR / "element_by_participant.tsv.gz"
    summ = _read(C.REPEAT_METHYLATION_DIR / "element_summary.tsv")
    contr = pd.DataFrame()
    if ebp_path.is_file() and not summ.empty:
        ebp = pd.read_csv(ebp_path, sep="\t")
        contr = followup.parental_contrasts(ebp, summ, bin_bp, C.MECHANISMS_BIPARENTAL)
        # outer: PWS-mUPD and AS-DEL/PWS-DEL contrasts also cover bins outside the common core
        reg = reg.merge(contr, on="bin_start", how="outer").sort_values("bin_start").reset_index(drop=True)
        reg["bin_end"] = reg["bin_start"] + bin_bp
        if contr["mupd_minus_biparental"].isna().all():
            log("no PWS-mUPD combined values in 07's output: rerun 07 (v6) for the mUPD check")
    genes = pd.DataFrame()
    if gtf and gtf.is_file():
        genes = followup.read_gtf_genes(gtf, C.CHROM, int(reg["bin_start"].min()) - 250_000,
                                        int(reg["bin_end"].max()) + 250_000)
        reg["genes"] = followup.genes_in_bins(reg, genes)
    else:
        log("no GTF ('gtf' in params.local.yml or --gtf): genes not annotated")
    doms = []
    for s, e, direction in followup.regional_domains(reg):
        row = {"start": s, "end": e, "size_mb": (e - s) / 1e6, "direction": direction,
               "bins": int(((reg["bin_start"] >= s) & (reg["bin_end"] <= e)).sum()),
               "median_raw_delta": float(reg.loc[(reg["bin_start"] >= s) & (reg["bin_end"] <= e),
                                                 "raw_median_delta"].median())}
        row.update(followup.contrasts_in_domain(contr, s, e))
        if len(genes):
            for side, pos in (("genes_at_start", s), ("genes_at_end", e)):
                near = genes[(genes["end"] > pos - 150_000) & (genes["start"] < pos + 150_000)]
                row[side] = followup.gene_names(list(near["gene"]), 20)
        doms.append(row)
    return reg, pd.DataFrame(doms, columns=None if doms else cols)


def _fmt(x, spec: str = ",") -> str:
    """Number for the report; 'NA' for None/NaN."""
    if x is None or (isinstance(x, float) and np.isnan(x)) or (isinstance(x, str) and x in ("", ".", "NA", "nan")):
        return "NA"
    if isinstance(x, str):
        return x
    if spec == ",":
        return f"{int(round(float(x))):,}"
    return format(float(x), spec)


def report(status: pd.DataFrame, fusions: pd.DataFrame, doms: pd.DataFrame, path: Path,
           split_ran: bool = True) -> None:
    lines = ["# Follow-up: breakpoints and repeat methylation", ""]
    lines.append("## Deletion carriers")
    lines.append("")
    for r in status.to_dict("records"):
        pair = r.get("sd_pair")
        pair_txt = f"; SD pair {pair}" if pair not in (None, ".", "") and pair == pair else ""
        split_txt = (f" Split reads: {_fmt(r.get('split_reads'))} in the carrier, "
                     f"max {_fmt(r.get('split_reads_panel_max'))} in the panel." if split_ran
                     else " Split reads: not counted.")
        lines.append(f"- **{r['sample']}** ({r['mechanism']}): **{r['status']}**. "
                     f"Edges {_fmt(r.get('proximal_edge'))} ({_fmt(r.get('proximal_cluster'))}) - "
                     f"{_fmt(r.get('distal_edge'))} ({_fmt(r.get('distal_cluster'))}){pair_txt}. "
                     f"Crossover: {_fmt(r.get('crossover_estimate'))}.{split_txt} "
                     f"Evidence: {str(r.get('evidence', '.')).rstrip('.') or 'none'}. "
                     f"Missing: {str(r.get('missing', '.')).rstrip('.') or 'nothing'}.")
    lines += ["", "## Assembled fusions", ""]
    if fusions.empty:
        lines.append("None with both ends in WINDOW and >= 50% of the HiFiCNV interval.")
    else:
        for r in fusions.to_dict("records"):
            nahr = ""
            if r.get("fusion_pair") not in (None, ".") and r.get("fusion_homolog_distance_bp") == r.get("fusion_homolog_distance_bp"):
                nahr = (f"; copies of {r['fusion_pair']}, {_fmt(r['fusion_homolog_distance_bp'])} bp "
                        "from homologous positions")
            lines.append(f"- {r['sample']} {r['hap']} {r['contig']}: {_fmt(r['ref_left'])} -> {_fmt(r['ref_right'])} "
                         f"(fraction of CNV deletion {_fmt(r['fraction_of_cnv_deletion'], '.2f')}{nahr})")
    lines += ["", "## Regions with a parental difference in repeat methylation", ""]
    if doms is None:
        lines.append("07's regional table (direct_delta_by_region.tsv) is missing: run 07.")
    elif doms.empty:
        lines.append("No run of 250-kb bins with label-permutation p <= 0.05 and |raw difference| >= 0.02.")
    else:
        for r in doms.to_dict("records"):
            gs, ge = (x if isinstance(x, str) and x else "none" for x in (r.get("genes_at_start"), r.get("genes_at_end")))
            genes = (f" Genes at start: {gs}; at end: {ge}."
                     if "genes_at_start" in r else "")
            contr = ""
            if "pws_del_inside" in r:
                contr = (" Minus biparental, inside / outside: "
                         f"PWS-DEL {_fmt(r.get('pws_del_inside'), '+.3f')} / {_fmt(r.get('pws_del_outside'), '+.3f')}, "
                         f"AS-DEL {_fmt(r.get('as_del_inside'), '+.3f')} / {_fmt(r.get('as_del_outside'), '+.3f')}, "
                         f"PWS-mUPD {_fmt(r.get('mupd_inside'), '+.3f')} / {_fmt(r.get('mupd_outside'), '+.3f')}.")
            lines.append(f"- chr15:{int(r['start']):,}-{int(r['end']):,} ({r['size_mb']:.2f} Mb, {r['bins']} bins): "
                         f"{r['direction']}, median raw difference {r['median_raw_delta']:+.3f}.{genes}{contr}")
        lines += ["", "Reading the contrasts: if the paternal chromosome is more methylated inside the region by d "
                  "(raw difference -d), PWS-DEL and PWS-mUPD should drop by about d/2 inside relative to outside "
                  "and AS-DEL should rise by about d/2 (columns *_inside_minus_outside in methylation_domains.tsv). "
                  "PWS-mUPD did not take part in finding the region, so it is the independent check."]
    path.write_text("\n".join(lines) + "\n")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--skip-split-reads", action="store_true")
    ap.add_argument("--window", type=int, default=50_000, help="bp around each edge for split reads")
    ap.add_argument("--gtf", help="CHM13 GTF (default: 'gtf' in params.local.yml)")
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    cohort = cis_analysis().load_cohort(C.METADATA_PATH)

    bp = _read(C.BREAKPOINT_DIR / "breakpoints_vs_hificnv.tsv")
    nahr = _read(C.BREAKPOINT_DIR / "nahr_test.tsv")
    enr = _read(C.BREAKPOINT_DIR / "junction_enrichment.tsv")
    if bp.empty or "bp_cluster" not in bp.columns:
        raise SystemExit("breakpoints_vs_hificnv.tsv is missing or from an old version: "
                         "run 05b_breakpoint_summary.py first")
    edges = edges_by_sample(bp)

    if a.skip_split_reads:
        split = _read(OUT / "split_reads.tsv")        # keep the counts of an earlier run
        log(f"split reads: {'reusing ' + str(len(split)) + ' rows of' if len(split) else 'none in'} "
            f"{OUT / 'split_reads.tsv'}")
    else:
        split = run_split_reads(edges, cohort, a.window)
        split.to_csv(OUT / "split_reads.tsv", sep="\t", index=False, na_rep="NA")

    w = parse_region(C.WINDOW)
    fusions = annotate_fusions(followup.assembled_fusions(C.ASSEMBLY_DIR, (w[1], w[2])))
    fusions.to_csv(OUT / "assembled_fusions.tsv", sep="\t", index=False, na_rep="NA")

    rows = []
    for s, e in edges.items():
        n = nahr[nahr["sample_id"] == s].iloc[0].to_dict() if len(nahr) and (nahr["sample_id"] == s).any() else {}
        n = {k: _num(v) for k, v in n.items()}
        sp = split[(split["carrier"] == s) & (split["genome"] == s)] if len(split) else pd.DataFrame()
        sp = sp.iloc[0].to_dict() if len(sp) else {}
        pan = split[(split["carrier"] == s) & (split["genome"] != s)] if len(split) else pd.DataFrame()
        pmax = int(pan.get("split_reads_deletion_type", pd.Series([0])).fillna(0).max()) if len(pan) else 0
        fu = fusions[fusions["sample"] == s] if len(fusions) else pd.DataFrame()
        fu = fu.sort_values("fraction_of_cnv_deletion", ascending=False).iloc[0].to_dict() if len(fu) else None
        ex = enr[(enr["sample"] == s) & (enr["excess"].astype(str) == "True")] if len(enr) else pd.DataFrame()
        ex = ex.sort_values(["pair_spans_carrier_edges", "poisson_p"], ascending=[False, True]).iloc[0].to_dict() \
            if len(ex) else None
        rows.append(followup.status_row(s, cohort.mechanism(s), e, n, {k: _num(v) for k, v in sp.items()},
                                        pmax, fu, ex))
    status = pd.DataFrame(rows)
    status.to_csv(OUT / "breakpoint_status.tsv", sep="\t", index=False, na_rep="NA")

    gtf = Path(a.gtf) if a.gtf else (Path(C.PARAMS["gtf"]) if C.PARAMS.get("gtf") else None)
    if gtf is not None and not gtf.is_absolute():
        gtf = C.PROJECT_ROOT / gtf
    reg, doms = methylation_tables(gtf)
    reg.to_csv(OUT / "methylation_regions.tsv", sep="\t", index=False, float_format="%.4g", na_rep="NA")
    (doms if doms is not None else pd.DataFrame()).to_csv(OUT / "methylation_domains.tsv", sep="\t", index=False, float_format="%.4g", na_rep="NA")
    report(status, fusions, doms, OUT / "followup_report.md", split_ran=len(split) > 0)
    with pd.option_context("display.width", 220, "display.max_columns", None, "display.max_colwidth", 70):
        print(status[["sample", "mechanism", "status", "crossover_estimate", "evidence"]].to_string(index=False))
    log(f"written: {OUT}")


if __name__ == "__main__":
    main()
