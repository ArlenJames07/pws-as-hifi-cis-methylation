"""
repeat_content.py -- repeats predicted on the assembled 15q11-q13 haplotypes,
placed on CHM13, per BP block, and the repeat content of every insertion or
deletion of the haplotype relative to CHM13.

For each sample and haplotype it needs:
  * RepeatMasker .out of the extracted haplotype pieces (S.hapN.15q.fa), and
  * the PAF of those pieces against the CHM13 window (S.hapN.15q.vs_window.paf,
    minimap2 -c --eqx, written by 03_place_contigs.py).
Plus the RepeatMasker .out of the CHM13 window itself, the BP landmarks and SDs.
Optional: contig_origin.tsv to label contigs maternal_like / paternal_like.

Outputs (--outdir, prefix = sample)
  S.repeat_blocks.tsv  per haplotype x feature (window, each BP, SD union) x repeat
                       class: bp on the haplotype vs bp in CHM13 over the SAME
                       aligned reference interval, so assembly gaps are not read
                       as repeat loss
  S.sv_repeats.tsv     insertions/deletions >= --min-sv vs CHM13 with their repeat
                       annotation (e.g. an AluY-sized insertion, a deleted L1)

What differences between haplotypes mean: a haplotype's repeats are inherited DNA
sequence; they are not parent-of-origin specific. Maternal-like vs paternal-like
comparisons in controls are the null expectation for structural polymorphism
between two homologs. In deletion carriers the retained chromosome came from the
NON-transmitting parent, so it samples the population, not the predisposing allele.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from .assembly import load_paf  # noqa: E402

CLASSES = ["SINE", "LINE", "LTR", "DNA", "Simple_repeat", "Low_complexity", "Satellite", "other"]


def cls_of(c):
    c = c.split("/")[0]
    return c if c in CLASSES else "other"


def read_rm_out(path):
    """RepeatMasker .out -> DataFrame(seq, start0, end, name, class, family, strand, div)."""
    rows = []
    with open(path) as fh:
        for i, line in enumerate(fh):
            f = line.split()
            if i < 3 or len(f) < 11 or not f[0].isdigit():
                continue
            cf = f[10].split("/")
            rows.append((f[4], int(f[5]) - 1, int(f[6]), f[9], cls_of(f[10]),
                         cf[1] if len(cf) > 1 else cf[0], "+" if f[8] == "+" else "-", float(f[1])))
    return pd.DataFrame(rows, columns=["seq", "start", "end", "name", "class", "family", "strand", "div"])


def window_offset(name):
    m = re.match(r"^(.+):(\d+)-(\d+)$", name)
    return (m.group(1), int(m.group(2)) - 1) if m else (name, 0)


def read_bed(path, chrom):
    out = []
    if not path:
        return out
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3 and f[0] == chrom:
            out.append((int(f[1]), int(f[2]), f[3] if len(f) > 3 else f"{f[1]}-{f[2]}"))
    return sorted(out)


def merge(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def overlap_bp(merged, a, b):
    return sum(max(0, min(e, b) - max(s, a)) for s, e in merged)


def annotate(hits, a, b):
    """repeat bp by class and top names for hits overlapping [a, b)."""
    h = hits[(hits["start"] < b) & (hits["end"] > a)]
    if h.empty:
        return 0, ".", "."
    ov = np.minimum(h["end"], b) - np.maximum(h["start"], a)
    by = ov.groupby(h["class"]).sum().sort_values(ascending=False)
    names = (ov.groupby(h["name"]).sum().sort_values(ascending=False).index[:3])
    return int(ov.sum()), ";".join(f"{c}:{int(v)}" for c, v in by.items()), ",".join(names)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--hap", action="append", required=True,
                    help="HAP=pieces.fa.out,pieces.vs_window.paf (repeatable)")
    ap.add_argument("--ref-out", required=True, help="RepeatMasker .out of the CHM13 window")
    ap.add_argument("--landmarks-bed")
    ap.add_argument("--sd-bed")
    ap.add_argument("--origin", help="contig_origin.tsv for this sample (optional)")
    ap.add_argument("--min-sv", type=int, default=50)
    ap.add_argument("--min-block", type=int, default=2000)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    pre = os.path.join(a.outdir, a.sample)

    ref = read_rm_out(a.ref_out)
    if ref.empty:
        sys.exit("empty reference RepeatMasker .out")
    chrom, off = window_offset(ref["seq"].iloc[0])
    ref["start"] += off
    ref["end"] += off
    lms = read_bed(a.landmarks_bed, chrom)
    sds = read_bed(a.sd_bed, chrom)
    origin = {}
    if a.origin:
        o = pd.read_csv(a.origin, sep="\t")
        origin = dict(zip(o["contig"], o["origin"]))

    block_rows, sv_rows = [], []
    for spec in a.hap:
        hap, rest = spec.split("=", 1)
        out_path, paf_path = rest.split(",", 1)
        hits = read_rm_out(out_path)
        blocks = load_paf(paf_path, hap, a.min_block)
        if not blocks:
            print(f"[{a.sample}] {hap}: no alignments in {paf_path}", file=sys.stderr)
            continue
        # every piece is named contig:start-end[/rc]; the contig name carries the origin
        piece_origin = {b.q: origin.get(b.q.split(":")[0], "unassigned") for b in blocks}

        # ---- project repeat midpoints of each piece onto CHM13
        proj = []
        for b in blocks:
            h = hits[hits["seq"] == b.q]
            if h.empty:
                continue
            mid = ((h["start"] + h["end"]) // 2).to_numpy()
            sel = (mid >= b.qs) & (mid < b.qe)
            t = b.q_to_t_array(mid[sel])
            ok = t >= 0
            hh = h[sel][ok].copy()
            hh["tmid"] = t[ok] + off
            hh["origin"] = piece_origin[b.q]
            proj.append(hh)
        proj = pd.concat(proj, ignore_index=True) if proj else pd.DataFrame(
            columns=list(hits.columns) + ["tmid", "origin"])
        proj["len"] = proj["end"] - proj["start"]

        # ---- per origin (within this hap) x feature x class
        for org in sorted(set(piece_origin.values())):
            bl = [b for b in blocks if piece_origin[b.q] == org]
            cov = merge([(b.ts + off, b.te + off) for b in bl])
            feats = [("window", cov[0][0] if cov else 0, cov[-1][1] if cov else 0)]
            feats += [(n, s, e) for s, e, n in lms]
            if sds:
                feats.append(("SD_union", min(s for s, _, _ in sds), max(e for _, e, _ in sds)))
            for fname, fs, fe in feats:
                cov_f = [[max(s, fs), min(e, fe)] for s, e in cov if s < fe and e > fs]
                if fname == "SD_union":
                    sdm = merge([(s, e) for s, e, _ in sds])
                    cov_f = [[max(s, cs), min(e, ce)] for cs, ce in cov_f for s, e in sdm if s < ce and e > cs]
                cbp = sum(e - s for s, e in cov_f)
                if cbp == 0:
                    continue
                p = proj[(proj["origin"] == org)]
                p = p[np.any([(p["tmid"] >= s) & (p["tmid"] < e) for s, e in cov_f], axis=0)] if len(p) else p
                r = ref[(ref["end"] > fs) & (ref["start"] < fe)]
                for c in CLASSES:
                    rc = r[r["class"] == c]
                    rbp = sum(overlap_bp(cov_f, s, e) for s, e in zip(rc["start"], rc["end"]))
                    hbp = int(p.loc[p["class"] == c, "len"].sum()) if len(p) else 0
                    hn = int((p["class"] == c).sum()) if len(p) else 0
                    if rbp == 0 and hbp == 0:
                        continue
                    block_rows.append({"sample": a.sample, "hap": hap, "origin": org, "feature": fname,
                                       "aligned_ref_bp": cbp, "class": c, "hap_repeat_bp": hbp,
                                       "hap_n": hn, "ref_repeat_bp": int(rbp),
                                       "delta_bp": hbp - int(rbp)})

        # ---- insertions / deletions vs CHM13 with repeat annotation
        by_piece = defaultdict(list)
        for b in blocks:
            by_piece[b.q].append(b)
            for kind, t, qa, qb, n in b.indels(a.min_sv):
                if kind == "INS":
                    rbp, rcls, rnames = annotate(hits[hits["seq"] == b.q], qa, qb)
                else:
                    rbp, rcls, rnames = annotate(ref, t + off, t + off + n)
                sv_rows.append((a.sample, hap, piece_origin[b.q], b.q, kind, "cigar", t + off,
                                t + off + (n if kind == "DEL" else 0), qa, qb, n, rbp,
                                round(rbp / n, 3), rcls, rnames))
        # larger events split alignments into consecutive blocks
        for q, bl in by_piece.items():
            bl = sorted(bl, key=lambda x: x.qs)
            for u, v in zip(bl, bl[1:]):
                if u.strand != v.strand:
                    continue
                qgap = v.qs - u.qe
                tgap = (v.ts - u.te) if u.strand == "+" else (u.ts - v.te)
                d = qgap - tgap
                if abs(d) < a.min_sv or tgap < -a.min_sv * 10:
                    continue
                if d > 0:      # extra sequence on the haplotype
                    rbp, rcls, rnames = annotate(hits[hits["seq"] == q], u.qe, v.qs)
                    tpos = (u.te if u.strand == "+" else u.ts) + off
                    sv_rows.append((a.sample, hap, piece_origin[q], q, "INS", "block_gap", tpos, tpos,
                                    u.qe, v.qs, d, rbp, round(rbp / max(1, qgap), 3), rcls, rnames))
                else:          # reference sequence missing from the haplotype
                    t0 = (u.te if u.strand == "+" else v.te) + off
                    t1 = t0 + max(tgap, -d)
                    rbp, rcls, rnames = annotate(ref, t0, t1)
                    sv_rows.append((a.sample, hap, piece_origin[q], q, "DEL", "block_gap", t0, t1,
                                    u.qe, v.qs, -d, rbp, round(rbp / max(1, t1 - t0), 3), rcls, rnames))

    B = pd.DataFrame(block_rows)
    B.to_csv(f"{pre}.repeat_blocks.tsv", sep="\t", index=False)
    SV = pd.DataFrame(sv_rows, columns=["sample", "hap", "origin", "piece", "type", "evidence",
                                        "ref_start", "ref_end", "hap_start", "hap_end", "length",
                                        "repeat_bp", "repeat_fraction", "repeat_classes", "top_repeats"])
    SV["landmark"] = pd.Series(".", index=SV.index, dtype=object)
    for s, e, n in lms:
        SV.loc[(SV["ref_start"] < e) & (SV["ref_end"] >= s), "landmark"] = n
    SV["in_sd"] = pd.Series(False, index=SV.index, dtype=bool)
    for s, e, _ in sds:
        SV.loc[(SV["ref_start"] < e) & (SV["ref_end"] >= s), "in_sd"] = True
    SV.to_csv(f"{pre}.sv_repeats.tsv", sep="\t", index=False)
    n_te = int((SV["repeat_fraction"] >= 0.8).sum()) if len(SV) else 0
    print(f"[{a.sample}] {len(SV)} indels >= {a.min_sv} bp vs CHM13, {n_te} of them >=80% repeat; "
          f"outputs {pre}.repeat_blocks.tsv / .sv_repeats.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
