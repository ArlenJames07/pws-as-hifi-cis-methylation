"""
placement.py -- place hifiasm contigs on the 15q11-q13 BP1-BP5 window of
T2T-CHM13 and report what the ASSEMBLY can and cannot say about the deletion.

Input : minimap2 PAF (contigs -> CHM13, run with -c so the cg:Z: CIGAR is present),
        one PAF per haplotype (hap1 / hap2 / primary), plus optional SD and
        BP-landmark BED files in CHM13 coordinates.

Output (in --outdir):
  <prefix>.blocks.tsv        every alignment block touching the window
  <prefix>.contigs.tsv       one row per contig: span, strand flips, where it ends
  <prefix>.junctions.tsv     adjacent blocks of one contig that jump >= --min-jump
                             on the reference (deletion-junction or inversion candidates)
  <prefix>.window_cover.tsv  per haplotype: fraction of window / SD / each BP covered
  <prefix>.<hap>.extract.plus.txt   samtools faidx region lists per haplotype:
  <prefix>.<hap>.extract.minus.txt  forward, and reverse-complement (`samtools faidx -i`),
                                    so every extracted piece is in CHM13 orientation

Reading the output
  * A contig that ENDS inside an SD is an assembly break, not a breakpoint.
  * A junction row with a jump matching the HiFiCNV deletion and a small query gap
    is a candidate junction-spanning contig. Rare with HiFi-only data; if found,
    confirm with the junction reads of 05_duplicon_breakpoints.py.
  * Strand flips inside a BP block on an intact haplotype are inversion candidates
    (e.g. BP2-BP3); they are the predisposing configurations described in the
    literature, but only matter for "propensity" if they sit on the transmitting
    parent's chromosome.

Needs numpy only for coordinate projection.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

from .assembly import (  # noqa: E402
    covered_bp as covered_bp_iv,
    load_paf,
    merge_intervals as merge,
    names_at as hits,
    names_overlapping as overlap_names,
    parse_region,
)


def read_bed(path, chrom):
    from .assembly import read_bed_intervals
    return read_bed_intervals(path, chrom)


def covered_bp(merged, a, b):
    return covered_bp_iv(merged, a, b)


# ------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", action="append", required=True,
                    help="HAP=path.paf (repeatable), e.g. hap1=PW1.hap1.paf")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--window", required=True,
                    help="CHM13 window, e.g. chr15:19500000-31500000 (BP1..BP5 plus flanks)")
    ap.add_argument("--sd-bed", help="CHM13 segmental duplications (BED, chrom/start/end[/name])")
    ap.add_argument("--landmarks-bed", help="BP1..BP5 intervals in CHM13 (BED with names)")
    ap.add_argument("--deletion", help="this sample's HiFiCNV deletion, chr15:START-END (optional)")
    ap.add_argument("--min-block", type=int, default=5000,
                    help="ignore alignment blocks shorter than this (bp)")
    ap.add_argument("--min-mapq", type=int, default=1,
                    help="blocks below this MAPQ are reported but flagged low_mapq")
    ap.add_argument("--min-jump", type=int, default=500_000,
                    help="reference jump between adjacent blocks to call a junction (bp)")
    ap.add_argument("--max-qgap", type=int, default=100_000,
                    help="max contig gap/overlap between adjacent blocks for a junction (bp)")
    ap.add_argument("--end-slop", type=int, default=2000,
                    help="a contig 'ends' at a block if the block reaches within this of the contig end")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args(argv)

    chrom, W0, W1 = parse_region(a.window)
    sds = read_bed(a.sd_bed, chrom)
    lms = read_bed(a.landmarks_bed, chrom)
    dele = parse_region(a.deletion) if a.deletion else None
    os.makedirs(a.outdir, exist_ok=True)
    pre = os.path.join(a.outdir, a.sample)

    # ---- load
    by_contig = defaultdict(list)
    for spec in a.paf:
        if "=" not in spec:
            sys.exit(f"--paf expects HAP=path, got '{spec}'")
        hap, path = spec.split("=", 1)
        for b in load_paf(path, hap, a.min_block):
            if b.t == chrom:
                by_contig[(hap, b.q)].append(b)

    # keep contigs with >=1 block overlapping the window
    in_win = {k: v for k, v in by_contig.items()
              if any(b.ts < W1 and b.te > W0 for b in v)}
    if not in_win:
        sys.exit("no contig aligns to the window -- check chromosome naming (chr15 vs 15) "
                 "and that contigs, not scaffolds split at N-gaps, were aligned")

    # ---- blocks table
    with open(f"{pre}.blocks.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "hap", "contig", "contig_len", "q_start", "q_end", "strand",
                    "chrom", "t_start", "t_end", "aln_len", "identity", "mapq", "low_mapq",
                    "type", "sd_at_t_start", "sd_at_t_end", "landmarks_overlapped"])
        for (hap, q), bl in sorted(in_win.items()):
            for b in sorted(bl, key=lambda x: x.qs):
                w.writerow([a.sample, hap, q, b.qlen, b.qs, b.qe, b.strand, b.t, b.ts, b.te,
                            b.alen, f"{b.identity:.5f}", b.mapq, int(b.mapq < a.min_mapq),
                            b.tp, ";".join(hits(sds, b.ts)) or ".",
                            ";".join(hits(sds, b.te - 1)) or ".",
                            ";".join(overlap_names(lms, b.ts, b.te)) or "."])

    # ---- per contig summary, junctions, extraction lists
    contig_rows, junction_rows = [], []
    plus_regions, minus_regions = defaultdict(list), defaultdict(list)
    for (hap, q), bl in sorted(in_win.items()):
        bl = sorted(bl, key=lambda x: x.qs)
        qlen = bl[0].qlen
        win_bl = [b for b in bl if b.ts < W1 and b.te > W0]
        strands = [b.strand for b in win_bl]
        main_strand = max(set(strands), key=lambda s: sum(b.alen for b in win_bl if b.strand == s))
        n_flips = sum(1 for i in range(1, len(win_bl)) if win_bl[i].strand != win_bl[i - 1].strand)

        # where does the contig end on the reference? only meaningful if the terminal
        # block reaches the physical contig end (otherwise the end is unaligned sequence)
        first, last = bl[0], bl[-1]
        ends = []
        for side, blk, qpos in (("contig_start", first, first.qs), ("contig_end", last, qlen - last.qe)):
            reaches = qpos <= a.end_slop
            if side == "contig_start":
                tpos = blk.ts if blk.strand == "+" else blk.te - 1
            else:
                tpos = blk.te - 1 if blk.strand == "+" else blk.ts
            in_sd = hits(sds, tpos)
            ends.append((side, reaches, tpos, in_sd, hits(lms, tpos)))

        # clipped extraction interval on the contig (window only)
        qcoords = []
        for b in win_bl:
            t0, t1 = max(b.ts, W0), min(b.te, W1)
            qa, qb = b.t_to_q(t0), b.t_to_q(t1)
            qcoords += [qa, qb]
        qx0, qx1 = max(0, min(qcoords)), min(qlen, max(qcoords))
        tx0 = min(max(b.ts, W0) for b in win_bl)
        tx1 = max(min(b.te, W1) for b in win_bl)
        region = f"{q}:{qx0 + 1}-{qx1}"
        (plus_regions if main_strand == "+" else minus_regions)[hap].append(region)

        contig_rows.append([
            a.sample, hap, q, qlen, len(win_bl), main_strand, n_flips,
            tx0, tx1, qx0, qx1, region,
            *[x for side in ends for x in (int(side[1]), side[2],
                                             ";".join(side[3]) or ".",
                                             ";".join(side[4]) or ".")]])

        # junction scan over ALL blocks of the contig (they may leave the window)
        for i in range(1, len(bl)):
            u, v = bl[i - 1], bl[i]
            if u.t != v.t:
                continue
            qgap = v.qs - u.qe
            if abs(qgap) > a.max_qgap:
                continue
            if u.strand == v.strand == "+":
                jump, kind = v.ts - u.te, None
            elif u.strand == v.strand == "-":
                jump, kind = u.ts - v.te, None
            else:
                jump, kind = abs(v.ts - u.ts), "strand_switch"
            if kind is None:
                if jump >= a.min_jump:
                    kind = "deletion_like"
                elif jump <= -a.min_jump:
                    kind = "duplication_or_reorder"
                else:
                    continue
            # strand switches are always reported: inside a BP block they are the
            # inversion candidates, whatever their size
            left = u.te if u.strand == "+" else u.ts
            right = v.ts if v.strand == "+" else v.te
            lo, hi = min(left, right), max(left, right)
            match = ""
            if dele and kind == "deletion_like":
                dchr, ds, de_ = dele
                ov = max(0, min(hi, de_) - max(lo, ds))
                match = f"{ov / max(1, de_ - ds):.2f}"
            junction_rows.append([
                a.sample, hap, q, kind, u.strand, v.strand, left, right, hi - lo, qgap,
                f"{u.identity:.4f}", f"{v.identity:.4f}", u.mapq, v.mapq,
                ";".join(hits(sds, left)) or ".", ";".join(hits(sds, right)) or ".",
                ";".join(hits(lms, left)) or ".", ";".join(hits(lms, right)) or ".", match])

    with open(f"{pre}.contigs.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "hap", "contig", "contig_len", "n_blocks_in_window",
                    "main_strand", "n_strand_flips", "t_start", "t_end", "q_start", "q_end",
                    "faidx_region",
                    "start_is_contig_end", "start_ref_pos", "start_in_sd", "start_landmark",
                    "end_is_contig_end", "end_ref_pos", "end_in_sd", "end_landmark"])
        w.writerows(contig_rows)

    with open(f"{pre}.junctions.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "hap", "contig", "kind", "strand_left", "strand_right",
                    "ref_left", "ref_right", "ref_span", "contig_gap", "id_left", "id_right",
                    "mapq_left", "mapq_right", "sd_left", "sd_right",
                    "landmark_left", "landmark_right", "fraction_of_cnv_deletion"])
        w.writerows(junction_rows)

    # ---- how much of the window / SDs / each BP did each haplotype assemble?
    haps = sorted({h for h, _ in in_win})
    sd_m = merge([(max(s, W0), min(e, W1)) for s, e, _ in sds if s < W1 and e > W0])
    with open(f"{pre}.window_cover.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "hap", "feature", "feature_bp", "covered_bp", "fraction",
                    "covered_bp_mapq_ok", "fraction_mapq_ok", "copies_in_assembly"])
        for hap in haps:
            allb = [b for (h, _), bl in in_win.items() if h == hap for b in bl]
            cov = merge([(max(b.ts, W0), min(b.te, W1)) for b in allb if b.ts < W1 and b.te > W0])
            covq = merge([(max(b.ts, W0), min(b.te, W1)) for b in allb
                          if b.ts < W1 and b.te > W0 and b.mapq >= a.min_mapq])
            feats = [("window", [[W0, W1]])]
            if dele:
                # how much of the HiFiCNV CN=1 interval this assembly carries
                feats.append(("cn1_deletion_interval", [[dele[1], dele[2]]]))
            if sd_m:
                feats.append(("all_SD_in_window", sd_m))
            for s, e, n in lms:
                feats.append((n, [[s, e]]))
            blocks_iv = [(max(b.ts, W0), min(b.te, W1)) for b in allb if b.ts < W1 and b.te > W0]
            for name, iv in feats:
                tot = sum(e - s for s, e in iv)
                c = sum(covered_bp(cov, s, e) for s, e in iv)
                cq = sum(covered_bp(covq, s, e) for s, e in iv)
                # >1 means the same reference stretch is present more than once in
                # this assembly (duplicated or collapsed-then-split sequence)
                mult = sum(max(0, min(be, e) - max(bs, s)) for bs, be in blocks_iv for s, e in iv)
                w.writerow([a.sample, hap, name, tot, c, f"{c / max(1, tot):.3f}",
                            cq, f"{cq / max(1, tot):.3f}", f"{mult / max(1, c):.2f}"])

    for hap in haps:
        for tag, dd in (("plus", plus_regions), ("minus", minus_regions)):
            with open(f"{pre}.{hap}.extract.{tag}.txt", "w") as fh:
                fh.write("".join(r + "\n" for r in dd[hap]))

    n_del = sum(1 for r in junction_rows if r[3] == "deletion_like")
    n_inv = sum(1 for r in junction_rows if r[3] == "strand_switch")
    print(f"[{a.sample}] {len(in_win)} contigs on {a.window}; "
          f"{n_del} deletion-like junction(s), {n_inv} strand switch(es). "
          f"Outputs: {pre}.*", file=sys.stderr)


if __name__ == "__main__":
    main()
