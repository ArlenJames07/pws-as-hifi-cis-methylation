"""One table of deletion edges for Figures 5 and 7, with the best available position of
each edge and its source:

  1  assembled contig switching between the two copies of a direct SD pair within 20 kb of
     homologous positions (scripts/duplicons 08: assembled_fusions / breakpoint_status)
  2  assembled contig joining the edge to another locus (another chromosome, or chr15 > 1 Mb
     away) -- found here in the hifiasm alignments of scripts/duplicons 03
  3  paralog-specific copy number (SUNK edge, scripts/duplicons 05)
  4  HiFiCNV CN=1 interval (scripts/analysis 01)

Class from the BP clusters of the two edges: BP1-BP3 class I, BP2-BP3 class II, else atypical."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import annotation as ann

CLASS_NAME = {("BP1", "BP3"): "class I (BP1–BP3)", ("BP2", "BP3"): "class II (BP2–BP3)"}
JUNCTION_SEARCH_BP = 100_000     # contig junction within this of an edge refines it
OTHER_EDGE_BP = 1_000_000        # a chr15 partner this close to the deletion's other edge is the
                                 # deletion junction itself (scripts/duplicons 08), not another locus
MAX_QUERY_GAP = 100_000          # as scripts/duplicons placement --max-qgap
MICROHOMOLOGY_MAX = 50           # larger overlaps of the two alignments are reported as overlap
MIN_BLOCK_BP = 20_000
MIN_MAPQ = 20


def _read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.is_file() and path.stat().st_size else pd.DataFrame()


def read_paf(path: Path) -> pd.DataFrame:
    cols = ["qname", "qlen", "qs", "qe", "strand", "tname", "tlen", "ts", "te", "nmatch", "alen", "mapq"]
    rows = []
    with open(path) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            tp = next((x[5:] for x in f[12:] if x.startswith("tp:A:")), "P")
            if tp not in ("P", "I"):
                continue
            rows.append([f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[6]), int(f[7]), int(f[8]),
                         int(f[9]), int(f[10]), int(f[11])])
    return pd.DataFrame(rows, columns=cols)


def contig_junctions(paths, sample: str, edges: list[int]) -> pd.DataFrame:
    """Consecutive confident alignment blocks of one contig (MAPQ >= 20, >= 20 kb; shorter or
    ambiguous blocks between them are counted in filtered_blocks_between), one on chr15 within
    JUNCTION_SEARCH_BP of an edge and the next on another chromosome, or on chr15 > 1 Mb away
    and not at the deletion's other edge (a contig across the deletion itself is scripts/duplicons
    08's assembled fusion)."""
    rows = []
    edges = [float(e) for e in edges if e is not None and np.isfinite(e)]
    folder = paths.duplicons / "assembly" / sample
    for paf in sorted(folder.glob(f"{sample}.*.paf")) if folder.is_dir() else []:
        if paf.name.endswith(".vs_window.paf"):
            continue
        hap = paf.name[len(sample) + 1:-4]
        raw = read_paf(paf)
        b = raw[(raw["mapq"] >= MIN_MAPQ) & (raw["alen"] >= MIN_BLOCK_BP)]
        for contig, d in b.groupby("qname"):
            d = d.sort_values("qs").reset_index(drop=True)
            r = raw[raw["qname"] == contig]
            for i in range(len(d) - 1):
                x, y = d.iloc[i], d.iloc[i + 1]
                gap = int(y["qs"] - x["qe"])
                if abs(gap) > MAX_QUERY_GAP:
                    continue
                between = int(((r["qs"] >= x["qe"] - 1_000) & (r["qe"] <= y["qs"] + 1_000)
                               & ~((r["qs"] == x["qs"]) & (r["qe"] == x["qe"]))
                               & ~((r["qs"] == y["qs"]) & (r["qe"] == y["qe"]))).sum())
                for first, second, order in ((x, y, "chr15 first"), (y, x, "chr15 second")):
                    if first["tname"] != ann.CHROM:
                        continue
                    if second["tname"] == ann.CHROM and abs(second["ts"] - first["ts"]) < 1_000_000:
                        continue
                    # position of the junction in each block (the side facing the other block)
                    if order == "chr15 first":
                        pos = first["ts"] if first["strand"] == "-" else first["te"]
                        partner = second["ts"] if second["strand"] == "+" else second["te"]
                    else:
                        pos = first["te"] if first["strand"] == "-" else first["ts"]
                        partner = second["te"] if second["strand"] == "+" else second["ts"]
                    near = [e for e in edges if abs(pos - e) <= JUNCTION_SEARCH_BP]
                    if not near:
                        continue
                    others = [e for e in edges if abs(pos - e) > JUNCTION_SEARCH_BP]
                    if second["tname"] == ann.CHROM and any(abs(partner - e) <= OTHER_EDGE_BP for e in others):
                        continue
                    overlap = max(0, -gap)
                    rows.append({"sample_id": sample, "assembly": hap, "contig": contig, "chr15_pos": int(pos),
                                 "partner_chrom": second["tname"], "partner_pos": int(partner),
                                 "chr15_block": f"{first['ts']}-{first['te']}({first['strand']})",
                                 "partner_block": f"{second['ts']}-{second['te']}({second['strand']})",
                                 "partner_identity": second["nmatch"] / max(1, second["alen"]),
                                 "microhomology_bp": overlap if overlap <= MICROHOMOLOGY_MAX else np.nan,
                                 "alignment_overlap_bp": overlap, "inserted_bp": max(0, gap),
                                 "filtered_blocks_between": between,
                                 "edge_refined": int(min(near, key=lambda e: abs(pos - e)))})
    return pd.DataFrame(rows)


def deletion_table(paths, cohort) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(one row per deletion carrier, contig junctions found)."""
    ev = _read(paths.analysis / "01_evidence_matrix" / "chr15_structural_evidence.tsv")
    status = _read(paths.duplicons / "followup" / "breakpoint_status.tsv")
    bvh = _read(paths.duplicons / "breakpoints" / "breakpoints_vs_hificnv.tsv")
    rows, junctions = [], []
    for s in cohort.of("PWS-DEL", "AS-DEL"):
        e = ev[ev["sample_id"] == s].iloc[0] if len(ev) and (ev["sample_id"] == s).any() else None
        h0 = float(e["cn_event_start"]) if e is not None and pd.notna(e["cn_event_start"]) else np.nan
        h1 = float(e["cn_event_end"]) if e is not None and pd.notna(e["cn_event_end"]) else np.nan
        row = {"sample_id": s, "group": cohort.group(s), "label": cohort.label(s), "hificnv_start": h0, "hificnv_end": h1,
               "start": h0, "end": h1, "start_source": "HiFiCNV", "end_source": "HiFiCNV",
               "status": ".", "crossover_estimate": ".", "sd_pair": ".",
               "analysis01_type": e["deletion_type"] if e is not None and "deletion_type" in e else "."}
        st = status[status["sample"] == s].iloc[0] if len(status) and (status["sample"] == s).any() else None
        if st is not None:
            row["status"] = st.get("status", ".")
            row["crossover_estimate"] = st.get("crossover_estimate", ".")
            row["sd_pair"] = st.get("sd_pair", ".")
            for side, col in (("start", "proximal_edge"), ("end", "distal_edge")):
                v = pd.to_numeric(st.get(col), errors="coerce")
                if pd.notna(v):
                    row[side], row[f"{side}_source"] = float(v), "SUNK copy number"
            fd = pd.to_numeric(st.get("fusion_homolog_distance_bp"), errors="coerce")
            if pd.notna(fd) and fd <= 20_000:
                # reference order: on a minus-strand contig ref_left is the distal edge
                a_, b_ = float(st["fusion_ref_left"]), float(st["fusion_ref_right"])
                row["start"], row["end"] = min(a_, b_), max(a_, b_)
                row["start_source"] = row["end_source"] = f"assembled contig ({st.get('fusion_contig', '.')})"
        elif len(bvh):
            for side, edge in (("start", "proximal"), ("end", "distal")):
                r = bvh[(bvh["sample_id"] == s) & (bvh["edge"] == edge)]
                if len(r) and pd.notna(r["sunk_pos"].iloc[0]):
                    row[side], row[f"{side}_source"] = float(r["sunk_pos"].iloc[0]), "SUNK copy number"
        j = contig_junctions(paths, s, [row["start"], row["end"], h0, h1])
        if len(j):
            junctions.append(j)
            for side in ("start", "end"):
                ref = row[side]
                near = j[(j["chr15_pos"] - ref).abs() <= JUNCTION_SEARCH_BP] if np.isfinite(ref) else j.iloc[0:0]
                if len(near) and not str(row[f"{side}_source"]).startswith("assembled"):
                    best = near.iloc[(near["chr15_pos"] - ref).abs().argmin()]
                    row[side] = float(best["chr15_pos"])
                    row[f"{side}_source"] = (f"contig junction to {best['partner_chrom']}:{best['partner_pos']:,} "
                                             f"({best['assembly']} {best['contig']})")
        c0 = ann.cluster_of(row["start"]) if np.isfinite(row["start"]) else "."
        c1 = ann.cluster_of(row["end"]) if np.isfinite(row["end"]) else "."
        row.update({"proximal_cluster": c0, "distal_cluster": c1,
                    "deletion_class": CLASS_NAME.get((c0, c1), "atypical"),
                    "size_mb": (row["end"] - row["start"]) / 1e6})
        rows.append(row)
    return pd.DataFrame(rows), (pd.concat(junctions, ignore_index=True) if junctions else pd.DataFrame())


def split_read_support(bam_path: Path, chrom: str, pos: int, partner_chrom: str, partner_pos: int,
                       window: int = 5_000, partner_window: int = 50_000, min_mapq: int = 1) -> int:
    """Reads near chrom:pos whose supplementary alignment lies near partner_chrom:partner_pos,
    plus reads near the partner whose supplementary alignment lies near chrom:pos (each read
    counted once). Both pieces MAPQ >= min_mapq; duplicates, QC failures and secondary
    alignments skipped (as scripts/duplicons 08)."""
    from .modbam import pysam_or_exit
    pysam = pysam_or_exit()
    names = set()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for here, there_chrom, there_pos, w_here, w_there in (
                ((chrom, pos), partner_chrom, partner_pos, window, partner_window),
                ((partner_chrom, partner_pos), chrom, pos, partner_window, window)):
            if here[0] not in bam.references:
                continue
            for read in bam.fetch(here[0], max(0, here[1] - w_here), here[1] + w_here):
                if (read.is_unmapped or read.is_secondary or read.is_duplicate or read.is_qcfail
                        or read.mapping_quality < min_mapq or not read.has_tag("SA")):
                    continue
                for sa in str(read.get_tag("SA")).rstrip(";").split(";"):
                    f = sa.split(",")          # rname,pos(1-based),strand,CIGAR,mapQ,NM
                    if (len(f) >= 5 and f[0] == there_chrom and int(f[4]) >= min_mapq
                            and abs(int(f[1]) - 1 - there_pos) <= w_there):
                        names.add(read.query_name)
                        break
    return len(names)
