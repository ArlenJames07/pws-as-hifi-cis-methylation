#!/usr/bin/env python3
"""
04 -- (optional) Parental label for the contig that carries the imprinting centre.

Sensitivity analysis only. Reads (with their MM/ML methylation tags) are mapped to
each participant's own diploid assembly, pb-CpG-tools is run on that alignment,
and the contig spanning the IC in each haplotype is labelled maternal_like
(methylated IC) or paternal_like. A hifiasm contig is a phase unit, so the label
holds along the whole contig, not only inside the read-phase block of the IC.
The labels are inferred from IC methylation: never report them as observed
parental states, and any statement about IC methylation itself is circular.

By default only Control and DiGeorge genomes are processed.

Output: results/08_duplicons/assembly/<sample>/
  <sample>.contig_origin.tsv
  <sample>.{maternal_like,paternal_like}.chm13.bed   CpGs projected to CHM13
"""
from __future__ import annotations

import argparse
import importlib
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
from duplicon_analysis import bams as bam_lookup  # noqa: E402
from duplicon_analysis import cis_analysis, origin  # noqa: E402
from duplicon_analysis.tools import log, parse_samples_arg, require, run, run_pipeline, to_fasta, up_to_date  # noqa: E402

place = importlib.import_module("03_place_contigs")


def label_sample(sample: str, mechanism: str, haps: dict[str, Path], bams: list[Path],
                 threads: int, force: bool) -> None:
    d = C.ASSEMBLY_DIR / sample
    outputs = [d / f"{sample}.contig_origin.tsv", d / f"{sample}.maternal_like.chm13.bed"]
    if up_to_date(outputs, force):
        log(f"{sample}: up to date")
        return
    pafs = [d / f"{sample}.{h}.paf" for h in ("hap1", "hap2")]
    if not all(p.is_file() for p in pafs):
        log(f"{sample}: run 03_place_contigs.py first; skipped")
        return
    dip = d / f"{sample}.dip.fa"
    with open(dip, "w") as out:
        for hap in ("hap1", "hap2"):
            part = d / f"{sample}.{hap}.tmp.fa"
            to_fasta(haps[hap], part)
            with open(part) as src:
                shutil.copyfileobj(src, out)
            part.unlink()
    run([C.SAMTOOLS, "faidx", dip])
    self_bam = d / f"{sample}.self.bam"
    # MAPQ >= 1 later drops reads in stretches where hap1 and hap2 are identical
    # reads are re-mapped here, so the reference of the original BAM does not matter
    reads = bam_lookup.fastq_command(C.SAMTOOLS, bams, threads, tags="MM,ML")
    run_pipeline(f"{reads} "
                 f"| {C.MINIMAP2} -t {threads} -y -ax map-hifi {dip} - "
                 f"| {C.SAMTOOLS} sort -@ {threads} -o {self_bam}")
    run([C.SAMTOOLS, "index", self_bam])
    prefix = d / f"{sample}.self.cpg"
    run([C.PBCPGTOOLS, "--bam", self_bam, "--output-prefix", prefix, "--model", C.CPG_MODEL,
         "--pileup-mode", "model", "--modsites-mode", "denovo",
         "--min-coverage", C.METHYLATION_MIN_COVERAGE, "--min-mapq", "1", "--threads", threads])
    combined = next(p for p in sorted(d.glob(f"{sample}.self.cpg.combined.bed*")) if not p.name.endswith(".tbi"))
    origin.main([f"--paf=hap1={pafs[0]}", f"--paf=hap2={pafs[1]}", "--meth", str(combined),
                 "--sample", sample, "--group", mechanism, "--window", C.WINDOW,
                 "--ic", f"{C.CHROM}:{C.PWS_IC_START}-{C.PWS_IC_END}", "--outdir", str(d)])
    for tmp in (self_bam, Path(f"{self_bam}.bai"), dip, Path(f"{dip}.fai")):
        tmp.unlink(missing_ok=True)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", help="default: all Control and DiGeorge participants")
    ap.add_argument("--threads", type=int, default=C.THREADS)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args(argv)
    require(C.SAMTOOLS, C.MINIMAP2, C.PBCPGTOOLS)
    cohort = cis_analysis().load_cohort(C.METADATA_PATH)
    wanted = parse_samples_arg(a.samples)
    assemblies = place.read_assemblies(C.ASSEMBLIES_PATH)
    table = bam_lookup.read_bam_table(C.BAMS_PATH)
    for sample, haps in assemblies.items():
        mechanism = cohort.mechanism(sample)
        chosen = sample in wanted if wanted else mechanism in C.MECHANISMS_BIPARENTAL
        if not chosen:
            continue
        if not {"hap1", "hap2"} <= set(haps):
            log(f"{sample}: needs hap1 and hap2 assemblies; skipped")
            continue
        try:
            bams = bam_lookup.find_bams(sample, C.ALIGNMENT_DIR, table)
        except bam_lookup.BamLookupError as error:
            log(f"SKIPPED {error}")
            continue
        if not bams:
            log(f"{sample}: no BAM found; skipped")
            continue
        if bam_lookup.has_methylation_tags(C.SAMTOOLS, bams[0], f"{C.CHROM}:{C.DOMAIN_START}-{C.DOMAIN_END}") is False:
            log(f"{sample}: reads in {bams[0].name} carry no MM/ML tags; skipped")
            continue
        label_sample(sample, mechanism, haps, bams, a.threads, a.force)


if __name__ == "__main__":
    main()
