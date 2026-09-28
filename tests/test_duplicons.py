from __future__ import annotations

import importlib.util
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DUPLICONS = ROOT / "scripts" / "duplicons"
sys.path.insert(0, str(DUPLICONS))

from duplicon_analysis import junction_reads  # noqa: E402
from duplicon_analysis.assembly import Block  # noqa: E402

COMP = str.maketrans("ACGT", "TGCA")


def load_script(filename: str, module_name: str):
    path = DUPLICONS / filename
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def random_seq(rng, n):
    return "".join(rng.choice(list("ACGT"), n))


def mutate(rng, seq, rate):
    chars = np.array(list(seq))
    for i in np.nonzero(rng.random(len(chars)) < rate)[0]:
        chars[i] = rng.choice([b for b in "ACGT" if b != chars[i]])
    return "".join(chars)


def canonical_kmers(seq, k):
    out = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        rc = kmer.translate(COMP)[::-1]
        out.append(min(kmer, rc))
    return out


class PafProjectionTests(unittest.TestCase):
    def test_projection_on_both_strands_and_indels(self) -> None:
        # '+' block with a 5-bp insertion after 100 matched bases
        plus = Block("hap1", "c1\t1000\t0\t305\t+\tchr15\t5000\t1000\t1300\t300\t305\t60\ttp:A:P\tcg:Z:100=5I200=".split("\t"))
        self.assertEqual(plus.q_to_t_array([50, 102, 150]).tolist(), [1050, -1, 1145])
        self.assertEqual(list(plus.indels(min_len=5)), [("INS", 1100, 100, 105, 5)])
        # '-' block: contig base qe-1 sits on the first reference base
        minus = Block("hap2", "c2\t500\t0\t300\t-\tchr15\t5000\t2000\t2300\t300\t300\t60\ttp:A:P\tcg:Z:300=".split("\t"))
        self.assertEqual(minus.q_to_t_array([299, 0]).tolist(), [2000, 2299])
        # a CpG starting at contig q on a '-' block is reported at the reference C (t of q+1)
        self.assertEqual(minus.q_to_t_array([198], cpg=True).tolist(), [2100])


class SunkBreakpointTests(unittest.TestCase):
    def test_nahr_switches_recovered_at_homologous_positions(self) -> None:
        rng = np.random.default_rng(11)
        k = 31
        sd_a = random_seq(rng, 20_000)
        sd_b = mutate(rng, sd_a, 0.005)
        ref = random_seq(rng, 40_000) + sd_a + random_seq(rng, 80_000) + sd_b + random_seq(rng, 40_000)
        a0, b0 = 40_000, 140_000
        control = random_seq(rng, 60_000)
        cross = 7_000
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            offset = 1_000_000
            (root / "window.fa").write_text(f">chr15:{offset + 1}-{offset + len(ref)}\n{ref}\n")
            (root / "control.fa").write_text(f">chr1:1-{len(control)}\n{control}\n")
            (root / "sd.bed").write_text(
                f"chr15\t{offset + a0}\t{offset + a0 + 20000}\tSDA\nchr15\t{offset + b0}\t{offset + b0 + 20000}\tSDB\n")
            (root / "pairs.bedpe").write_text(
                f"chr15\t{offset + a0}\t{offset + a0 + 20000}\tchr15\t{offset + b0}\t{offset + b0 + 20000}\tp\t0.995\t+\t+\n")
            ref_counts = Counter(canonical_kmers(ref, k) + canonical_kmers(control, k))
            sunks = {kmer for kmer, n in ref_counts.items() if n == 1}
            (root / "sunks.txt").write_text("".join(f"{s}\t1\n" for s in sunks))
            rows = ["sample\tgroup\tcounts"]

            def write_sample(name, group, haplotypes):
                counts = Counter()
                for hap in list(haplotypes) + [control, control]:
                    counts.update(canonical_kmers(hap, k))
                with (root / f"{name}.txt").open("w") as handle:
                    for kmer in sunks:
                        value = rng.poisson(10 * counts.get(kmer, 0))
                        if value:
                            handle.write(f"{kmer}\t{value}\n")
                rows.append(f"{name}\t{group}\t{root / f'{name}.txt'}")

            for i in range(4):
                write_sample(f"C{i}", "Control", [mutate(rng, ref, 0.001), mutate(rng, ref, 0.001)])
            deleted = ref[:a0 + cross] + ref[b0 + cross:]
            write_sample("P1", "PWS-DEL", [mutate(rng, ref, 0.001), deleted])
            (root / "samples.tsv").write_text("\n".join(rows) + "\n")
            subprocess.run(
                [sys.executable, "-m", "duplicon_analysis.sunk",
                 "--window-fasta", str(root / "window.fa"), "--control-fasta", str(root / "control.fa"),
                 "--sunks", str(root / "sunks.txt"), "--samples", str(root / "samples.tsv"),
                 "--panel-groups", "Control", "--sd-bed", str(root / "sd.bed"),
                 "--sd-pairs", str(root / "pairs.bedpe"), "--bin-bp", "2000",
                 "--outdir", str(root / "out")],
                check=True, capture_output=True, cwd=DUPLICONS)
            transitions = pd.read_csv(root / "out" / "sunk15q.transitions.tsv", sep="\t")
            p1 = transitions[transitions["sample"] == "P1"].sort_values("refined_pos")
            self.assertEqual(p1[["from_cn", "to_cn"]].values.tolist(), [[2, 1], [1, 2]])
            self.assertLess(abs(p1.iloc[0]["refined_pos"] - (offset + a0 + cross)), 1_000)
            self.assertLess(abs(p1.iloc[1]["refined_pos"] - (offset + b0 + cross)), 1_000)
            self.assertTrue(transitions[transitions["sample"] != "P1"].empty)
            nahr = pd.read_csv(root / "out" / "sunk15q.nahr_pairs.tsv", sep="\t")
            self.assertTrue(bool(nahr.iloc[0]["consistent_with_NAHR"]))


class JunctionReadTests(unittest.TestCase):
    def test_fusion_read_gives_crossover_and_normal_reads_do_not(self) -> None:
        rng = np.random.default_rng(3)
        k = 31
        sd_a = random_seq(rng, 20_000)
        sd_b = mutate(rng, sd_a, 0.01)
        ref = random_seq(rng, 30_000) + sd_a + random_seq(rng, 50_000) + sd_b + random_seq(rng, 30_000)
        a0, b0 = 30_000, 100_000
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "window.fa").write_text(f">chr15:1-{len(ref)}\n{ref}\n")
            counts = Counter(canonical_kmers(ref, k))
            (root / "sunks.txt").write_text("".join(f"{s}\t1\n" for s, n in counts.items() if n == 1))
            (root / "pairs.bedpe").write_text(
                f"chr15\t{a0}\t{a0 + 20000}\tchr15\t{b0}\t{b0 + 20000}\tp1\t0.99\t+\t+\n")
            pairs, index = junction_reads.prepare(root / "window.fa", root / "sunks.txt", root / "pairs.bedpe", k)
            hybrid = ref[:a0 + 8_000] + ref[b0 + 8_000:]
            fusion_read = hybrid[a0:a0 + 15_000]                       # crosses the crossover at offset 8,000
            rc = fusion_read[::-1].translate(COMP)
            normal_read = ref[a0 + 2_000:a0 + 17_000]
            for read in (fusion_read, rc):
                call = junction_reads.classify_read(read, index, k)[0]
                self.assertTrue(call["ordered_A_then_B"] and call["offsets_consistent"])
                self.assertLessEqual(call["crossover_offset_low"], 8_000)
                self.assertGreaterEqual(call["crossover_offset_high"], 8_000 - k)
                self.assertLess(call["crossover_offset_high"] - call["crossover_offset_low"], 1_000)
            self.assertEqual(junction_reads.classify_read(normal_read, index, k), [])


class RepeatMethylationTests(unittest.TestCase):
    def test_direct_contrast_and_unoriented_asm(self) -> None:
        rng = np.random.default_rng(5)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cohort = ([(f"P{i}", "PWS-DEL") for i in range(3)] + [(f"A{i}", "AS-DEL") for i in range(2)]
                      + [(f"C{i}", "Control") for i in range(2)] + [(f"D{i}", "DiGeorge") for i in range(2)]
                      + [("U0", "PWS-mUPD")])
            pd.DataFrame(cohort, columns=["sample_id", "molecular_mechanism"]).to_csv(root / "metadata.csv", index=False)
            pd.DataFrame([{"sample_id": s, "ic_deletion_status": "confirmed", "cn_event_start": 1_000,
                           "cn_event_end": 99_000, "cn_event_mean_cn": 1.0, "evidence_basis": "CN"}
                          for s, m in cohort if m in ("PWS-DEL", "AS-DEL")]).to_csv(
                root / "structural.tsv", sep="\t", index=False)
            starts = np.arange(2_000, 90_000, 2_000)
            (root / "repeats.bed").write_text("".join(
                f"chr15\t{s}\t{s + 1000}\trep{j}\tSINE\tAlu\n" for j, s in enumerate(starts)))
            (root / "sd.bed").write_text("")
            parental = {int(starts[5]): (0.9, 0.1), int(starts[20]): (0.1, 0.9)}
            cpg = np.arange(1_500, 95_000, 50)

            def values(maternal: bool):
                v = np.full(cpg.size, 0.7)
                for s, (m, p) in parental.items():
                    inside = (cpg >= s) & (cpg < s + 1000)
                    v[inside] = m if maternal else p
                return np.clip(v + rng.normal(0, 0.03, cpg.size), 0, 1)

            def write(sample, kind, beta):
                d = root / "meth" / sample
                d.mkdir(parents=True, exist_ok=True)
                lines = ["##pb-cpg-tools-version=3.0.0",
                         "#chrom\tbegin\tend\tmod_score\ttype\tcov\test_mod_count\test_unmod_count"]
                lines += [f"chr15\t{p}\t{p + 2}\t{100 * b:.1f}\t.\t10\t{10 * b:.1f}\t{10 - 10 * b:.1f}"
                          for p, b in zip(cpg, beta)]
                (d / f"{sample}.cpg.{kind}.bed").write_text("\n".join(lines) + "\n")

            for sample, mech in cohort:
                if mech == "PWS-DEL":
                    write(sample, "combined", values(True))
                elif mech == "AS-DEL":
                    write(sample, "combined", values(False))
                elif mech == "PWS-mUPD":
                    write(sample, "hap1", values(True)); write(sample, "hap2", values(True))
                else:
                    m, p = values(True), values(False)
                    h1, h2 = (m, p) if rng.random() < 0.5 else (p, m)   # unoriented
                    write(sample, "hap1", h1); write(sample, "hap2", h2); write(sample, "combined", (m + p) / 2)

            script = load_script("07_repeat_methylation.py", "repeat_methylation_test")
            script.METADATA_PATH = root / "metadata.csv"
            script.METHYLATION_DIR = root / "meth"
            script.STRUCTURAL_EVIDENCE_PATH = root / "structural.tsv"
            script.REPEATS_PATH = root / "repeats.bed"
            script.SD_BED_PATH = root / "sd.bed"
            script.OUTPUT_DIR = root / "out"
            script.PHASING_DIR = root / "phasing"
            script.ASSEMBLY_DIR = root / "assembly"
            script.DOMAIN_START, script.DOMAIN_END = 1_000, 99_000
            script.PWS_IC_START, script.PWS_IC_END = 97_000, 97_500
            script.CN1_BREAKPOINT_BUFFER = 0
            script.MIN_BIPARENTAL = 3
            script.BOOTSTRAP_REPLICATES = 200
            script.CIRCULAR_SHIFTS = 200
            script.run()
            summary = pd.read_csv(root / "out" / "element_summary.tsv", sep="\t").set_index("name")
            self.assertEqual(summary.loc["rep5", "direction"], "maternal_higher")
            self.assertEqual(summary.loc["rep20", "direction"], "paternal_higher")
            self.assertTrue(summary.loc[["rep5", "rep20"], "asm_candidate"].all())
            self.assertEqual(int(summary["direct_parent_of_origin"].sum()), 2)
            self.assertLess(summary.loc["rep5", "mupd_abs_h1_minus_h2"], 0.1)
            self.assertGreater(summary.loc["rep5", "asm_mean"], 0.7)


class InputLookupTests(unittest.TestCase):
    """BAM and methylation-track lookup for workflow names and SMRT Link names."""

    def test_bam_names(self) -> None:
        from duplicon_analysis import bams
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp) / "01_alignment"
            d.mkdir()
            for name in ("08_1_A01_bc2043_001P.bam", "08_1_A01_bc2044_002P.bam",
                         "08_1_A01_bc2044v2_002P.bam", "x_1005P.bam", "017C.aligned.bam",
                         "08_1_C01_bc2059_017C.bam"):
                (d / name).write_text("")
                (d / f"{name}.bai").write_text("")
            (d / "08_1_B01_bc2048_006P.bam").symlink_to(d / "does_not_exist.bam")
            self.assertEqual([p.name for p in bams.find_bams("001P", d, {})], ["08_1_A01_bc2043_001P.bam"])
            # the workflow name wins over a SMRT Link name
            self.assertEqual([p.name for p in bams.find_bams("017C", d, {})], ["017C.aligned.bam"])
            # the ID must be a separate token: 005P is not x_1005P
            self.assertEqual(bams.find_bams("005P", d, {}), [])
            with self.assertRaisesRegex(bams.BamLookupError, "2 BAMs match"):
                bams.find_bams("002P", d, {})
            with self.assertRaisesRegex(bams.BamLookupError, "broken link"):
                bams.find_bams("006P", d, {})
            table_path = Path(tmp) / "bams.local.csv"
            table_path.write_text("sample,bam\n002P,01_alignment/08_1_A01_bc2044v2_002P.bam\n"
                                  f"013A,{d}/08_1_A01_bc2043_001P.bam;{d}/017C.aligned.bam\n")
            table = bams.read_bam_table(table_path)
            self.assertEqual([p.name for p in bams.find_bams("002P", d, table)], ["08_1_A01_bc2044v2_002P.bam"])
            self.assertEqual(len(bams.find_bams("013A", d, table)), 2)
            self.assertIsNotNone(bams.bam_index(d / "08_1_A01_bc2043_001P.bam"))
            command = bams.fastq_command("samtools", bams.find_bams("013A", d, table), 4, tags="MM,ML")
            self.assertTrue(command.startswith("{ samtools fastq -@ 4 -T MM,ML "))

    def test_methylation_track_names(self) -> None:
        from duplicon_analysis.methylation import find_track
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "001P").mkdir()
            (root / "001P" / "001P.cpg.combined.bed.gz").write_text("")
            (root / "smrt").mkdir()
            for name in ("08_1_A01_bc2055_013A.combined.bed.gz", "08_1_A01_bc2055_013A.combined.bed.gz.tbi",
                         "08_1_A01_bc2055_013A.hap1.bed.gz", "08_1_A01_bc2044_002P.hap1.bed.gz",
                         "08_1_A01_bc2044v2_002P.hap1.bed.gz"):
                (root / "smrt" / name).write_text("")
            self.assertEqual(find_track(root, "001P", "combined").name, "001P.cpg.combined.bed.gz")
            self.assertEqual(find_track(root, "013A", "combined").name, "08_1_A01_bc2055_013A.combined.bed.gz")
            self.assertIsNone(find_track(root, "013A", "hap2"))
            with self.assertRaises(ValueError):
                find_track(root, "002P", "hap1")


if __name__ == "__main__":
    unittest.main()
