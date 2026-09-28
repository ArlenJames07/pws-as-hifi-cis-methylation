"""
contig_origin.py -- give hifiasm contigs a parental label from their own
imprinting-centre methylation, and project their CpG methylation onto CHM13.

Why: read-based phasing (HP/PS tags) breaks into blocks; orienting HP1/HP2 by the
IC is only valid inside the phase block that contains the IC. A hifiasm contig is
itself a phase unit and in unique sequence often spans Mb, so the contig that
carries the IC carries a parental label along its whole length (it still stops at
the BP blocks). This gives maternal-like / paternal-like tracks in biparental
genomes that do not depend on read-phase-block length.

Inputs
  --paf hap1=S.hap1.paf --paf hap2=S.hap2.paf   contigs -> CHM13 (minimap2 -c)
  --meth S.self.combined.bed                     per-CpG methylation on the diploid
                                                 assembly itself (pb-CpG-tools run on
                                                 reads mapped back to hap1+hap2; use
                                                 MAPQ >= 1 so reads in identical
                                                 hap1/hap2 stretches are excluded)
Outputs (--outdir)
  S.contig_origin.tsv              contig, IC methylation, label, CHM13 span
  S.maternal_like.chm13.bed        chrom start end meth cov (CHM13 coordinates)
  S.paternal_like.chm13.bed
Labels are "maternal_like"/"paternal_like": they come from IC methylation, so any
statement about IC methylation itself is circular; everything away from the IC is not.
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from .assembly import load_paf, parse_region  # noqa: E402
from .methylation import read_track  # noqa: E402


def read_contig_methylation(path: Path, contigs: dict[str, int], min_cov: int) -> pd.DataFrame:
    """pb-CpG-tools BED on the self-assembly, restricted to the given contigs. One pass
    keeps only those contigs (plus the header, so the percent scale is detected)."""
    keep = set(contigs)
    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as tmp:
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as handle:
            for line in handle:
                if line.startswith("#") or line.split("\t", 1)[0] in keep:
                    tmp.write(line)
        tmp_path = Path(tmp.name)
    frames = []
    try:
        for contig, length in contigs.items():
            track = read_track(tmp_path, contig, 0, length)
            ok = track.coverage >= min_cov
            frames.append(pd.DataFrame({"chrom": contig, "pos": track.position[ok],
                                        "meth": track.beta[ok], "cov": track.coverage[ok]}))
    finally:
        tmp_path.unlink(missing_ok=True)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
        columns=["chrom", "pos", "meth", "cov"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", action="append", required=True)
    ap.add_argument("--meth", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--group", default="")
    ap.add_argument("--ic", default="chr15:22691258-22693494")
    ap.add_argument("--window", default="chr15:17600000-28000000")
    ap.add_argument("--min-cov", type=int, default=4)
    ap.add_argument("--min-ic-cpgs", type=int, default=5)
    ap.add_argument("--min-ic-diff", type=float, default=0.30,
                    help="the two IC contigs must differ by this much to be labelled")
    ap.add_argument("--min-block", type=int, default=5000)
    ap.add_argument("--min-mapq", type=int, default=1)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    pre = os.path.join(a.outdir, a.sample)
    ic_c, ic_s, ic_e = parse_region(a.ic)
    w_c, w_s, w_e = parse_region(a.window)

    blocks = []
    for spec in a.paf:
        hap, path = spec.split("=", 1)
        blocks += [b for b in load_paf(path, hap, a.min_block) if b.t == w_c and b.mapq >= a.min_mapq]
    by_contig = {}
    for b in blocks:
        by_contig.setdefault((b.hap, b.q), []).append(b)
    in_win = {k: v for k, v in by_contig.items() if any(b.ts < w_e and b.te > w_s for b in v)}
    names = {q: bl[0].qlen for (_, q), bl in in_win.items()}
    meth = read_contig_methylation(Path(a.meth), names, a.min_cov)
    if meth.empty:
        sys.exit("no CpGs on the 15q contigs -- contig names in --meth must match the assembly")

    # project every CpG of every window contig to CHM13
    proj = []
    for (hap, q), bl in in_win.items():
        m = meth[meth["chrom"] == q]
        if m.empty:
            continue
        pos = m["pos"].to_numpy()
        for b in bl:
            sel = (pos >= b.qs) & (pos < b.qe)
            if not sel.any():
                continue
            t = b.q_to_t_array(pos[sel], cpg=True)
            ok = t >= 0
            proj.append(pd.DataFrame({"hap": hap, "contig": q, "tpos": t[ok],
                                      "meth": m["meth"].to_numpy()[sel][ok],
                                      "cov": m["cov"].to_numpy()[sel][ok]}))
    proj = pd.concat(proj, ignore_index=True) if proj else pd.DataFrame(
        columns=["hap", "contig", "tpos", "meth", "cov"])

    # IC methylation per contig
    rows = []
    for (hap, q), bl in in_win.items():
        pi = proj[(proj["contig"] == q) & (proj["tpos"] >= ic_s) & (proj["tpos"] < ic_e)]
        wb = [b for b in bl if b.ts < w_e and b.te > w_s]
        rows.append({"sample": a.sample, "group": a.group, "hap": hap, "contig": q,
                     "contig_len": bl[0].qlen,
                     "chm13_start": max(w_s, min(b.ts for b in wb)),
                     "chm13_end": min(w_e, max(b.te for b in wb)),
                     "spans_ic": any(b.ts <= ic_s and b.te >= ic_e for b in bl),
                     "ic_cpgs": len(pi), "ic_meth": pi["meth"].mean() if len(pi) else np.nan})
    tab = pd.DataFrame(rows)
    tab["origin"] = "unassigned"
    ic_tab = tab[tab["spans_ic"] & (tab["ic_cpgs"] >= a.min_ic_cpgs)].sort_values("ic_meth")
    note = ""
    if len(ic_tab) >= 2:
        lo, hi = ic_tab.iloc[0], ic_tab.iloc[-1]
        if lo["hap"] != hi["hap"] and hi["ic_meth"] - lo["ic_meth"] >= a.min_ic_diff:
            tab.loc[tab["contig"] == hi["contig"], "origin"] = "maternal_like"
            tab.loc[tab["contig"] == lo["contig"], "origin"] = "paternal_like"
        else:
            note = ("IC contigs not separable (same hap, or IC methylation difference "
                    f"{hi['ic_meth'] - lo['ic_meth']:.2f} < {a.min_ic_diff}); a homozygous IC or a "
                    "hemizygous locus duplicated into both haps looks like this")
    else:
        note = f"{len(ic_tab)} contig(s) span the IC with >= {a.min_ic_cpgs} CpGs; cannot orient"
    tab.to_csv(f"{pre}.contig_origin.tsv", sep="\t", index=False, float_format="%.4f")

    for lab in ("maternal_like", "paternal_like"):
        cs = set(tab.loc[tab["origin"] == lab, "contig"])
        d = proj[proj["contig"].isin(cs) & (proj["tpos"] >= w_s) & (proj["tpos"] < w_e)]
        d = d.sort_values("tpos")
        with open(f"{pre}.{lab}.chm13.bed", "w") as fh:
            for t, mm, cv in zip(d["tpos"], d["meth"], d["cov"]):
                fh.write(f"{w_c}\t{t}\t{t + 2}\t{mm:.4f}\t{cv:g}\n")
    lab = tab[tab["origin"] != "unassigned"]
    span = "; ".join(f"{r.origin}: {r.contig} {r.chm13_start:,}-{r.chm13_end:,}"
                     for r in lab.itertuples())
    print(f"[{a.sample}] {span or 'no contig oriented'}{(' -- ' + note) if note else ''}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
