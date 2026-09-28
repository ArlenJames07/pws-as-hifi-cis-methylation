#!/usr/bin/env python3
"""
03 -- Place each participant's hifiasm contigs on T2T-CHM13 and extract 15q11-q13.

Input: assemblies.local.csv at the repository root (see config.ASSEMBLIES_PATH):
    sample,hap1,hap2,primary
    017C,/data/asm/017C.bp.hap1.p_ctg.gfa,/data/asm/017C.bp.hap2.p_ctg.gfa,
    001P,/data/asm/001P.bp.hap1.p_ctg.gfa,/data/asm/001P.bp.hap2.p_ctg.gfa,/data/asm/001P.bp.p_ctg.gfa
Any non-empty subset of the three columns is accepted. Biparental genomes: hap1 +
hap2. PWS-DEL / AS-DEL: add the primary contigs (*.bp.p_ctg.gfa) -- inside the
hemizygous interval the one existing chromosome is always in p_ctg, whereas hifiasm
may drop, fragment or duplicate it in hap1/hap2. Check which assembly carries it in
<sample>.window_cover.tsv (feature cn1_deletion_interval). Outside that interval
p_ctg is a mosaic of both parental chromosomes.
Use the hifiasm CONTIGS (GFA or FASTA, optionally gzipped), not RagTag scaffolds:
scaffolds are ordered and oriented by the reference and cannot show a breakpoint
or an inversion.

Per participant (results/08_duplicons/assembly/<sample>/):
  <sample>.<asm>.paf                 contigs vs CHM13 (minimap2 -cx asm5 --eqx); <asm> = hap1/hap2/primary
  <sample>.{contigs,blocks,junctions,window_cover}.tsv
                                     where contigs end (inside SDs = assembly break),
                                     strand switches (inversion candidates), contigs
                                     that jump across the deletion, BP-block coverage
  <sample>.hapN.15q.fa               the window, per haplotype, in CHM13 orientation
  <sample>.hapN.15q.vs_window.paf    those pieces vs the CHM13 window
  <sample>.hapN.15q.fa.out           RepeatMasker on the pieces (if installed)
The whole-genome contig FASTA is deleted after extraction (04 rebuilds it).
"""
from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
from duplicon_analysis import placement  # noqa: E402
from duplicon_analysis.tools import log, parse_samples_arg, require, run, to_fasta, up_to_date  # noqa: E402

ASSEMBLY_COLUMNS = ("hap1", "hap2", "primary")


def read_assemblies(path: Path) -> dict[str, dict[str, Path]]:
    """sample -> {'hap1': path, 'hap2': path, 'primary': path} (only the columns filled in)."""
    if not path.is_file():
        raise SystemExit(f"{path} not found: create it with columns sample,hap1,hap2[,primary]")
    rows = {}
    with open(path) as handle:
        for row in csv.DictReader(handle):
            entry = {}
            for column in ASSEMBLY_COLUMNS:
                value = (row.get(column) or "").strip()
                if not value:
                    continue
                p = Path(value)
                entry[column] = p if p.is_absolute() else path.parent / p
                if not entry[column].is_file():
                    raise SystemExit(f"{row['sample']} {column}: {entry[column]} not found")
            if not entry:
                raise SystemExit(f"{row['sample']}: no assembly given (hap1, hap2 or primary)")
            rows[row["sample"].strip()] = entry
    return rows


def deletion_interval(sample: str) -> str | None:
    """HiFiCNV CN=1 event of a PWS-DEL/AS-DEL participant, from analysis 01."""
    if not C.STRUCTURAL_EVIDENCE_PATH.is_file():
        return None
    table = pd.read_csv(C.STRUCTURAL_EVIDENCE_PATH, sep="\t")
    row = table[table["sample_id"] == sample]
    if row.empty or pd.isna(row.iloc[0].get("cn_event_start")):
        return None
    return f"{C.CHROM}:{int(row.iloc[0]['cn_event_start'])}-{int(row.iloc[0]['cn_event_end'])}"


def place_sample(sample: str, haps: dict[str, Path], threads: int, force: bool) -> None:
    d = C.ASSEMBLY_DIR / sample
    d.mkdir(parents=True, exist_ok=True)
    window = C.REFERENCE_DIR / "window.fa"
    present = [h for h in ASSEMBLY_COLUMNS if h in haps]
    final = [d / f"{sample}.{h}.15q.vs_window.paf" for h in present]
    if up_to_date(final, force):
        log(f"{sample}: placement up to date")
    else:
        for hap in present:
            fasta = d / f"{sample}.{hap}.contigs.fa"
            to_fasta(haps[hap], fasta)
            run([C.SAMTOOLS, "faidx", fasta])
            run([C.MINIMAP2, "-t", threads, "-cx", "asm5", "--eqx", "--secondary=no", C.REFERENCE, fasta],
                stdout=d / f"{sample}.{hap}.paf")
        args = [f"--paf={h}={d / f'{sample}.{h}.paf'}" for h in present]
        args += ["--sample", sample, "--window", C.WINDOW, "--outdir", str(d),
                 "--sd-bed", str(C.REFERENCE_DIR / "window.sd.bed")]
        deletion = deletion_interval(sample)
        if deletion:
            args += ["--deletion", deletion]
        placement.main(args)
        for hap in present:
            fasta = d / f"{sample}.{hap}.contigs.fa"
            pieces = d / f"{sample}.{hap}.15q.fa"
            with open(pieces, "w"):
                pass
            for strand, extra in (("plus", []), ("minus", ["-i"])):
                regions = d / f"{sample}.{hap}.extract.{strand}.txt"
                if regions.is_file() and regions.stat().st_size:
                    part = d / f"{sample}.{hap}.{strand}.fa"
                    run([C.SAMTOOLS, "faidx", *extra, fasta, "-r", regions], stdout=part)
                    with open(pieces, "a") as out, open(part) as src:
                        shutil.copyfileobj(src, out)
                    part.unlink()
            run([C.MINIMAP2, "-t", threads, "-cx", "asm5", "--eqx", window, pieces],
                stdout=d / f"{sample}.{hap}.15q.vs_window.paf")
            fasta.unlink(missing_ok=True)
            Path(f"{fasta}.fai").unlink(missing_ok=True)

    if shutil.which(C.REPEATMASKER) is None:
        log(f"{sample}: {C.REPEATMASKER} not found; RepeatMasker on the pieces skipped (needed by 06)")
        return
    for hap in present:
        pieces = d / f"{sample}.{hap}.15q.fa"
        out = d / f"{sample}.{hap}.15q.fa.out"
        if pieces.stat().st_size and not up_to_date([out], force):
            run([C.REPEATMASKER, "-species", "human", "-xsmall", "-pa", threads, "-dir", d, pieces])


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples")
    ap.add_argument("--threads", type=int, default=C.THREADS)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args(argv)
    require(C.SAMTOOLS, C.MINIMAP2)
    if not (C.REFERENCE_DIR / "window.fa").is_file():
        raise SystemExit("run 01_prepare_reference.py first")
    wanted = parse_samples_arg(a.samples)
    for sample, haps in read_assemblies(C.ASSEMBLIES_PATH).items():
        if wanted is None or sample in wanted:
            place_sample(sample, haps, a.threads, a.force)


if __name__ == "__main__":
    main()
