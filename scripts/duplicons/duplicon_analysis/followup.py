"""
Follow-up checks on the breakpoint and repeat-methylation results (used by 08_followup.py).

  split_reads()          reads whose primary alignment lies at one deletion edge and whose
                         supplementary alignment (SA tag) lies at the other: the junction at
                         base-pair level when at least one side is unique sequence
  assembled_fusions()    hifiasm contigs that jump across a deletion (03's junctions.tsv),
                         with the nearest confidently placed blocks on each side
  fusion_nahr()          is an assembled junction at homologous positions of a direct SD pair?
  parental_contrasts()   PWS-DEL, AS-DEL and PWS-mUPD minus biparental genomes, per bin
  genes_in_bins()        genes overlapping each 250-kb bin of 07's regional table
  status_row()           per-carrier verdict and what is still missing
"""
from __future__ import annotations

import gzip
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

CIGAR_REF = re.compile(r"(\d+)([MIDNSHP=X])")


def ref_length(cigar: str) -> int:
    """Reference bases spanned by a CIGAR string."""
    return sum(int(n) for n, op in CIGAR_REF.findall(cigar) if op in "MDN=X")


def parse_sa(tag_value: str) -> list[tuple[str, int, str, str, int]]:
    """SA:Z value -> [(chrom, pos 1-based, strand, cigar, mapq)]."""
    out = []
    for item in tag_value.strip(";").split(";"):
        f = item.split(",")
        if len(f) >= 5:
            out.append((f[0], int(f[1]), f[2], f[3], int(f[4])))
    return out


def junction_from_pieces(a: tuple[int, int, str], b: tuple[int, int, str]):
    """Two alignment pieces (start0, end0, strand) of one read -> (kind, proximal bp, distal bp).
    Same strand = deletion-type junction: proximal piece ends at the proximal breakpoint and
    the distal piece starts at the distal breakpoint."""
    (s1, e1, t1), (s2, e2, t2) = sorted([a, b])
    kind = "deletion" if t1 == t2 else "inversion"
    return kind, e1, s2


def split_reads(samtools: str, bams: list[Path], chrom: str, edge_a: int, edge_b: int,
                window: int = 50_000, min_mapq: int = 1) -> pd.DataFrame:
    """Reads with one alignment piece within +/- window of edge_a and another within
    +/- window of edge_b (either piece may be the primary). One row per read."""
    rows, seen = [], set()
    for here, there in ((edge_a, edge_b), (edge_b, edge_a)):
        region = f"{chrom}:{max(1, here - window)}-{here + window}"
        for bam in bams:
            proc = subprocess.run([samtools, "view", "-F", "0x904", str(bam), region],
                                  capture_output=True, text=True)
            if proc.returncode != 0:
                continue
            for line in proc.stdout.splitlines():
                f = line.split("\t")
                if len(f) < 11 or f[0] in seen or int(f[4]) < min_mapq:
                    continue
                sa = next((x[5:] for x in f[11:] if x.startswith("SA:Z:")), None)
                if not sa:
                    continue
                start0 = int(f[3]) - 1
                primary = (start0, start0 + ref_length(f[5]), "-" if int(f[1]) & 16 else "+")
                for c, pos, strand, cigar, mq in parse_sa(sa):
                    s0 = pos - 1
                    if c != chrom or mq < min_mapq or abs(s0 - there) > window + 50_000:
                        continue
                    piece = (s0, s0 + ref_length(cigar), strand)
                    kind, bp1, bp2 = junction_from_pieces(primary, piece)
                    rows.append({"read": f[0], "kind": kind, "proximal_bp": bp1, "distal_bp": bp2,
                                 "mapq_primary": int(f[4]), "mapq_supplementary": mq})
                    seen.add(f[0])
                    break
    return pd.DataFrame(rows, columns=["read", "kind", "proximal_bp", "distal_bp",
                                       "mapq_primary", "mapq_supplementary"])


def summarise_split(df: pd.DataFrame) -> dict:
    if df.empty:
        return {"split_reads": 0}
    d = df[df["kind"] == "deletion"]
    out = {"split_reads": int(len(df)), "split_reads_deletion_type": int(len(d))}
    if len(d):
        out.update({"split_junction_proximal": int(d["proximal_bp"].median()),
                    "split_junction_distal": int(d["distal_bp"].median()),
                    "split_junction_spread_bp": int(max(d["proximal_bp"].max() - d["proximal_bp"].min(),
                                                        d["distal_bp"].max() - d["distal_bp"].min()))})
    return out


def assembled_fusions(assembly_dir: Path, window: tuple[int, int], min_fraction: float = 0.5,
                      anchor_identity: float = 0.995, anchor_mapq: int = 20) -> pd.DataFrame:
    """deletion_like junctions covering >= min_fraction of the HiFiCNV interval, both ends in
    WINDOW, with the last/first confidently placed block (identity, MAPQ) on each side."""
    rows = []
    for jpath in sorted(assembly_dir.glob("*/*.junctions.tsv")):
        j = pd.read_csv(jpath, sep="\t")
        if j.empty:
            continue
        j = j[(j["kind"] == "deletion_like")
              & (pd.to_numeric(j["fraction_of_cnv_deletion"], errors="coerce") >= min_fraction)
              & j["ref_left"].between(*window) & j["ref_right"].between(*window)]
        if j.empty:
            continue
        bpath = jpath.with_name(jpath.name.replace(".junctions.tsv", ".blocks.tsv"))
        blocks = pd.read_csv(bpath, sep="\t") if bpath.is_file() else pd.DataFrame()
        for r in j.itertuples(index=False):
            row = r._asdict()
            if not blocks.empty:
                b = blocks[(blocks["hap"] == r.hap) & (blocks["contig"] == r.contig)]
                good = b[(b["identity"] >= anchor_identity) & (b["mapq"] >= anchor_mapq)]
                left = good[good["t_end"] <= r.ref_left + 1_000]
                right = good[good["t_start"] >= r.ref_right - 1_000]
                row["last_confident_proximal_pos"] = int(left["t_end"].max()) if len(left) else np.nan
                row["first_confident_distal_pos"] = int(right["t_start"].min()) if len(right) else np.nan
                row["blocks_between_anchors"] = int(len(b) - len(left) - len(right)) if len(b) else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def fusion_nahr(ref_left: int, ref_right: int, pairs) -> tuple[str, int, int] | None:
    """Directly oriented pair with one copy at each side of an assembled junction ->
    (pair, homologous position of ref_left in copy B, distance to ref_right)."""
    best = None
    for pr in pairs:
        if pr.a0 - 5_000 <= ref_left <= pr.a1 + 5_000 and pr.b0 - 5_000 <= ref_right <= pr.b1 + 5_000:
            frac = (min(max(ref_left, pr.a0), pr.a1) - pr.a0) / max(1, pr.a1 - pr.a0)
            hom = int(round(pr.b0 + frac * (pr.b1 - pr.b0)))
            if best is None or abs(hom - ref_right) < best[2]:
                best = (pr.name, hom, abs(hom - ref_right))
    return best


CONTRASTS = ("pws_del_minus_biparental", "as_del_minus_biparental", "mupd_minus_biparental")


def parental_contrasts(ebp: pd.DataFrame, summary: pd.DataFrame, bin_bp: int = 250_000,
                       biparental=("Control", "DiGeorge")) -> pd.DataFrame:
    """Per bin, medians over testable elements of (group - biparental) combined
    methylation, for PWS-DEL (maternal copy only), AS-DEL (paternal copy only) and
    PWS-mUPD (two maternal copies). If the paternal chromosome is more methylated by d,
    PWS-DEL and PWS-mUPD sit d/2 below the biparental genomes and AS-DEL d/2 above.
    PWS-mUPD was not used to find the region, so it is an independent check.
    Columns are NaN when 07 did not store that group (mUPD needs 07 from v6 on)."""
    cols = ["bin_start", *CONTRASTS, "elements_pws", "elements_as", "elements_mupd"]
    if ebp.empty or summary.empty:
        return pd.DataFrame(columns=cols)
    s = summary[summary["testable"].astype(str).isin(["True", "true", "1"])].set_index("element_id")
    e = ebp[ebp["element_id"].isin(s.index)]

    def per_element(measure, mechs):
        sub = e[(e["measure"] == measure) & e["mechanism"].isin(mechs)]
        return sub.groupby("element_id")["value"].median()

    table = pd.DataFrame(index=s.index)
    table["bin_start"] = (s["start"].astype(int) // bin_bp) * bin_bp
    bip = per_element("scaffold_combined", list(biparental)).reindex(table.index)
    for name, measure, mechs in (("pws_del_minus_biparental", "retained_combined", ["PWS-DEL"]),
                                 ("as_del_minus_biparental", "retained_combined", ["AS-DEL"]),
                                 ("mupd_minus_biparental", "scaffold_combined", ["PWS-mUPD"])):
        table[name] = per_element(measure, mechs).reindex(table.index) - bip
    g = table.groupby("bin_start")
    out = g[list(CONTRASTS)].median()
    for name, short in zip(CONTRASTS, ("pws", "as", "mupd")):
        out[f"elements_{short}"] = g[name].count()
    return out.reset_index()[cols]


def contrasts_in_domain(contr: pd.DataFrame, start: int, end: int) -> dict:
    """Median of the per-bin contrasts inside [start, end) and in all other bins, and the
    difference (inside - outside). The difference removes offsets shared by the whole
    window (e.g. between deletion carriers and biparental genomes)."""
    out = {}
    if contr.empty:
        return out
    inside = (contr["bin_start"] >= start) & (contr["bin_start"] < end)
    for c in CONTRASTS:
        a = contr.loc[inside, c].median()
        b = contr.loc[~inside, c].median()
        short = c.replace("_minus_biparental", "")
        out[f"{short}_inside"] = a
        out[f"{short}_outside"] = b
        out[f"{short}_inside_minus_outside"] = a - b if a == a and b == b else np.nan
    return out


def read_gtf_genes(gtf: Path, chrom: str, start: int, end: int) -> pd.DataFrame:
    """Genes (feature 'gene', or transcripts collapsed by name when no gene lines) in a region."""
    opener = gzip.open if str(gtf).endswith(".gz") else open
    rows = []
    with opener(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != chrom or f[2] not in ("gene", "transcript"):
                continue
            s, e = int(f[3]) - 1, int(f[4])
            if e < start or s > end:
                continue
            m = re.search(r'gene_name "([^"]+)"', f[8]) or re.search(r'gene_id "([^"]+)"', f[8])
            t = re.search(r'gene_(?:bio)?type "([^"]+)"', f[8])
            rows.append({"feature": f[2], "gene": m.group(1) if m else ".", "start": s, "end": e,
                         "strand": f[6], "type": t.group(1) if t else "."})
    g = pd.DataFrame(rows, columns=["feature", "gene", "start", "end", "strand", "type"])
    if g.empty:
        return g
    if (g["feature"] == "gene").any():
        g = g[g["feature"] == "gene"]
    return (g.groupby("gene").agg(start=("start", "min"), end=("end", "max"), strand=("strand", "first"),
                                  type=("type", "first")).reset_index().sort_values("start"))


FAMILY = re.compile(r"^(.+?)[-_.](\d+)$")


def gene_names(names, max_names: int = 25) -> str:
    """Comma-separated names in order, numbered copies of one family collapsed
    (SNORD116-1 ... SNORD116-30 -> 'SNORD116 x30'), so that a snoRNA cluster does not
    push the protein-coding genes out of the list."""
    order, count = [], {}
    for n in names:
        m = FAMILY.match(n)
        key = m.group(1) if m and m.group(1).upper().startswith(("SNORD", "SNORA", "MIR")) else n
        if key not in count:
            order.append(key)
        count[key] = count.get(key, 0) + 1
    labels = [f"{k} x{count[k]}" if count[k] > 1 else k for k in order]
    return ",".join(labels[:max_names]) + (f",+{len(labels) - max_names}" if len(labels) > max_names else "")


def genes_in_bins(bins: pd.DataFrame, genes: pd.DataFrame, max_names: int = 25) -> list[str]:
    out = []
    for b in bins.itertuples(index=False):
        hit = genes[(genes["end"] > b.bin_start) & (genes["start"] < b.bin_end)] if len(genes) else genes
        out.append(gene_names(list(hit["gene"]) if len(hit) else [], max_names))
    return out


def regional_domains(reg: pd.DataFrame, alpha: float = 0.05, min_abs_raw: float = 0.02) -> list[tuple[int, int, str]]:
    """Runs of consecutive bins with permutation p <= alpha, |raw difference| >= min_abs_raw
    and the same sign of the RAW difference. The centred value has the global offset
    removed, so a region can look shifted only because another region is; shifts smaller
    than min_abs_raw are of the size of that offset. Returns (start, end, direction)."""
    runs, cur = [], None
    for r in reg.sort_values("bin_start").itertuples(index=False):
        ok = (r.permutation_p_two_sided == r.permutation_p_two_sided and r.permutation_p_two_sided <= alpha
              and abs(r.raw_median_delta) >= min_abs_raw)
        sign = "paternal_higher" if r.raw_median_delta < 0 else "maternal_higher"
        if ok and cur and cur[2] == sign and cur[1] == r.bin_start:
            cur = (cur[0], r.bin_end, sign)
        elif ok:
            if cur:
                runs.append(cur)
            cur = (r.bin_start, r.bin_end, sign)
        else:
            if cur:
                runs.append(cur)
            cur = None
    if cur:
        runs.append(cur)
    return runs


def status_row(sample: str, mechanism: str, edges: dict, nahr: dict, split: dict,
               panel_split_max: int, fusion: dict | None, enrichment_at_edges: dict | None) -> dict:
    """Combine the evidence of one carrier into a verdict and the next step."""
    prox, dist = edges.get("proximal", {}), edges.get("distal", {})
    fusion_d = fusion.get("fusion_homolog_distance_bp") if fusion else None
    fusion_d = None if fusion_d is None or fusion_d != fusion_d else fusion_d
    fusion_nahr_ok = fusion_d is not None and fusion_d <= 20_000
    row = {"sample": sample, "mechanism": mechanism,
           "proximal_edge": prox.get("pos"), "proximal_cluster": prox.get("cluster", "."),
           "distal_edge": dist.get("pos"), "distal_cluster": dist.get("cluster", "."),
           "sd_pair": nahr.get("pair", "."), "nahr_homolog_distance_bp": nahr.get("distance_bp"),
           "nahr_consistent": nahr.get("consistent_with_NAHR"),
           "junction_reads_at_crossover": nahr.get("junction_reads_at_crossover"),
           "junction_reads_expected": nahr.get("expected_from_panel_at_crossover"),
           "junction_p_at_crossover": nahr.get("poisson_p_at_crossover"),
           "assembled_fusion": bool(fusion), "split_reads": split.get("split_reads") or 0,
           "split_reads_panel_max": panel_split_max}
    if fusion:
        row.update({"fusion_contig": f"{fusion['hap']}:{fusion['contig']}",
                    "fusion_ref_left": fusion["ref_left"], "fusion_ref_right": fusion["ref_right"],
                    "fusion_pair": fusion.get("fusion_pair", "."), "fusion_homolog_distance_bp": fusion_d})
    if (split.get("split_reads_deletion_type") or 0) > 0:
        row.update({"split_junction": f"{split['split_junction_proximal']}-{split['split_junction_distal']}"})
    if enrichment_at_edges:
        row.update({"enriched_pair_at_edges": enrichment_at_edges["pair"],
                    "enriched_pair_p": enrichment_at_edges["poisson_p"]})

    evidence, missing = [], []
    nahr_ok = bool(nahr.get("consistent_with_NAHR"))
    jr_ok = (nahr.get("poisson_p_at_crossover") is not None
             and nahr.get("poisson_p_at_crossover") == nahr.get("poisson_p_at_crossover")
             and nahr.get("poisson_p_at_crossover") < 1e-3)
    split_ok = (split.get("split_reads_deletion_type") or 0) >= 2 and panel_split_max == 0
    if nahr_ok:
        evidence.append("edges at homologous positions of one SD pair")
    if jr_ok:
        evidence.append("junction-read excess at its own crossover")
    if fusion_nahr_ok:
        evidence.append(f"assembled contig switches between the copies of {fusion.get('fusion_pair')} "
                        f"{int(fusion_d)} bp from homologous positions")
    elif fusion:
        evidence.append("assembled contig across the deletion (its switch is not at homologous positions "
                        "of a direct pair)")
    if split_ok:
        evidence.append("split reads give the junction (bp)")
    if prox.get("pos") is None or dist.get("pos") is None:
        missing.append("an edge lies outside WINDOW: analyse a supplementary window (README) "
                       "or rely on split reads")
    pair = nahr.get("pair") if nahr.get("pair") not in (None, ".") else None
    edge_missing = prox.get("pos") is None or dist.get("pos") is None
    if split_ok:
        status = "junction resolved at bp level"
    elif (nahr_ok and jr_ok) or fusion_nahr_ok:
        status = "confirmed NAHR"
    elif fusion:
        status = "confirmed by assembly"
        missing.append("the assembled junction is not at homologous positions of a direct SD pair "
                       "(fusion_pair/fusion_homolog_distance_bp): check the contig blocks and inverted pairs")
    elif nahr_ok:
        status = "compatible with NAHR"
        missing.append("independent confirmation: assembled fusion or junction-read excess at the crossover")
    elif edge_missing:
        status = "edge outside WINDOW"
    elif pair:
        status = "edges near one SD pair, not homologous"
        d = nahr.get("distance_bp")
        missing.append(f"edges are {d / 1000:.1f} kb from homologous positions in {pair} (limit 20 kb): "
                       "look at both zoom plots; confirm with an assembled fusion or split reads")
    else:
        status = "unresolved"
        if enrichment_at_edges:
            missing.append(f"junction-read excess in {enrichment_at_edges['pair']} but its copies are not at the "
                           "SUNK edges: look at both zoom plots and the HiFiCNV edges")
        else:
            missing.append("no directly oriented SD pair links the two edges: look at the zoom plots, the CN of "
                           "the flanks, and inverted pairs")
    if prox.get("note") or dist.get("note"):
        missing.append("flank note: " + "; ".join(x for x in (prox.get("note"), dist.get("note")) if x))
    # best position of the crossover. An assembled switch that is not at homologous positions
    # is where the aligner left copy A inside a duplicon (it may place the whole hybrid copy on
    # copy B), not the crossover itself; when the SUNK edges are homologous they are better.
    sunk_ok = prox.get("pos") is not None and dist.get("pos") is not None
    if split_ok:
        row["crossover_estimate"] = f"{row.get('split_junction')} (split reads, bp)"
    elif fusion_nahr_ok:
        row["crossover_estimate"] = f"{int(fusion['ref_left'])}-{int(fusion['ref_right'])} (assembly)"
    elif nahr_ok and sunk_ok:
        row["crossover_estimate"] = f"{int(prox['pos'])}-{int(dist['pos'])} (SUNK, several kb)"
    elif fusion:
        row["crossover_estimate"] = (f"{int(fusion['ref_left'])}-{int(fusion['ref_right'])} "
                                     "(assembly, not homologous)")
    elif sunk_ok:
        row["crossover_estimate"] = f"{int(prox['pos'])}-{int(dist['pos'])} (SUNK, several kb)"
    else:
        row["crossover_estimate"] = "."
    row.update({"status": status, "evidence": "; ".join(evidence) or ".",
                "missing": "; ".join(missing) or "."})
    return row
