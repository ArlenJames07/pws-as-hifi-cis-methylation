from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "analysis"))

from cis_analysis import classify_evidence, load_deletion_map, read_track


def load_analysis_script(filename: str, module_name: str):
    path = ROOT / "scripts" / "analysis" / filename
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class AnalysisCoreTests(unittest.TestCase):
    def test_region_catalog_without_pandas_runtime(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gtf = root / "genes.gtf"
            names = ("MAGEL2", "NDN", "SNRPN", "SNHG14", "SNORD116-1", "UBE3A", "GABRB3", "GABRA5", "GABRG3", "OCA2")
            lines = [
                f'chr15\ttest\tgene\t{1000 + i * 1000}\t{1500 + i * 1000}\t.\t+\t.\tgene_id "G{i}"; gene_name "{name}";'
                for i, name in enumerate(names)
            ]
            gtf.write_text("\n".join(lines) + "\n")
            region_script = load_analysis_script(
                "00_build_prespecified_regions.py", "analysis_regions_test"
            )
            region_script.GTF_PATH = gtf
            region_script.OUTPUT_PATH = root / "prespecified_regions.tsv"
            rows = region_script.build_regions()
            self.assertEqual(len(rows), 7)
            self.assertEqual(rows[-1]["region_id"], "PWS/AS imprinting centre")
            region_script.main()
            self.assertTrue(region_script.OUTPUT_PATH.is_file())

    def test_pbcpg_bed_columns(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.combined.bed"
            path.write_text(
                "##pb-cpg-tools-version=3.0.0\n"
                "#chrom\tbegin\tend\tmod_score\ttype\tcov\test_mod_count\test_unmod_count\n"
                "chr15\t100\t102\t80\t.\t4\t3.2\t0.8\n"
                "chr15\t110\t112\t20\t.\t8\t1.6\t6.4\n"
                "chr16\t100\t102\t90\t.\t20\t18\t2\n"
            )
            track = read_track(path, "chr15", 90, 120)
            np.testing.assert_array_equal(track.position, [100, 110])
            np.testing.assert_allclose(track.beta, [0.8, 0.2])
            np.testing.assert_allclose(track.coverage, [4, 8])

    def test_parental_direction_is_cn1_limited(self) -> None:
        inside = classify_evidence(
            "PWS-DEL", "combined", 200, 300, (100, 500), (150, 450)
        )
        outside = classify_evidence(
            "PWS-DEL", "combined", 500, 600, (100, 500), (150, 450)
        )
        control = classify_evidence(
            "Control", "hap1", 200, 300, None, (150, 450)
        )
        self.assertTrue(inside["parental_direction_observed"])
        self.assertEqual(inside["parental_state"], "maternal_retained")
        self.assertFalse(outside["parental_direction_observed"])
        self.assertIsNone(outside["parental_state"])
        self.assertEqual(control["evidence_class"], "diploid_unoriented")
        self.assertIsNone(control["parental_state"])

    def test_common_cn1_intersection(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "provenance.tsv"
            pd.DataFrame(
                [
                    {
                        "sample_id": "P1",
                        "ic_deletion_status": "confirmed",
                        "cn_event_start": 100,
                        "cn_event_end": 1000,
                        "cn_event_mean_cn": 1.0,
                        "evidence_basis": "CN",
                    },
                    {
                        "sample_id": "A1",
                        "ic_deletion_status": "confirmed",
                        "cn_event_start": 200,
                        "cn_event_end": 900,
                        "cn_event_mean_cn": 1.0,
                        "evidence_basis": "CN",
                    },
                ]
            ).to_csv(path, sep="\t", index=False)
            deletion_map = load_deletion_map(path, ("P1", "A1"), 50)
            self.assertEqual((deletion_map.common_start, deletion_map.common_end), (250, 850))

    def test_end_to_end_analysis_contract(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            methylation = root / "methylation"
            samples = {
                "P1": "PWS-DEL",
                "P2": "PWS-DEL",
                "P3": "PWS-DEL",
                "A1": "AS-DEL",
                "A2": "AS-DEL",
                "C1": "Control",
                "C2": "Control",
                "D1": "DiGeorge",
            }
            metadata = root / "metadata.csv"
            pd.DataFrame(
                {
                    "sample_id": list(samples),
                    "molecular_mechanism": list(samples.values()),
                }
            ).to_csv(metadata, index=False)
            provenance = root / "provenance.tsv"
            pd.DataFrame(
                [
                    {
                        "sample_id": sample,
                        "ic_deletion_status": "confirmed",
                        "cn_event_start": 100,
                        "cn_event_end": 900,
                        "cn_event_mean_cn": 1.0,
                        "evidence_basis": "synthetic CN=1",
                    }
                    for sample, mechanism in samples.items()
                    if mechanism in {"PWS-DEL", "AS-DEL"}
                ]
            ).to_csv(provenance, sep="\t", index=False)
            for sample, mechanism in samples.items():
                kinds = ["combined"]
                if mechanism in {"Control", "DiGeorge"}:
                    kinds.extend(["hap1", "hap2"])
                sample_dir = methylation / sample
                sample_dir.mkdir(parents=True)
                for kind in kinds:
                    if mechanism == "PWS-DEL":
                        score = 80
                    elif mechanism == "AS-DEL":
                        score = 20
                    elif kind == "hap1":
                        score = 80
                    elif kind == "hap2":
                        score = 20
                    else:
                        score = 50
                    lines = [
                        "##pb-cpg-tools-version=3.0.0",
                        "##pileup-mode=model",
                        "##modsites-mode=denovo",
                        "##min-coverage=4",
                        "##min-mapq=1",
                    ]
                    for window_start in range(100, 900, 100):
                        for offset in (10, 30, 50):
                            position = window_start + offset
                            lines.append(
                                f"chr15\t{position}\t{position + 2}\t{score}\t.\t8\t0\t0\t{score}"
                            )
                    (sample_dir / f"{sample}.cpg.{kind}.bed").write_text("\n".join(lines) + "\n")
            out = root / "analysis"
            first = out / "01_evidence_matrix"
            second = out / "02_depth_sensitivity"
            third = out / "03_cis_architecture"
            matrix_path = first / "chr15_window_evidence_matrix.tsv.gz"
            matrix_script = load_analysis_script(
                "01_build_chr15_evidence_matrix.py", "analysis_matrix_test"
            )
            matrix_script.METADATA_PATH = metadata
            matrix_script.METHYLATION_DIR = methylation
            matrix_script.DELETION_PROVENANCE_PATH = provenance
            matrix_script.OUTPUT_DIR = first
            matrix_script.CHROM = "chr15"
            matrix_script.DOMAIN_START = 100
            matrix_script.DOMAIN_END = 900
            matrix_script.WINDOW_SIZE = 100
            matrix_script.CN1_BREAKPOINT_BUFFER = 0
            matrix_script.MIN_CPGS = 3
            matrix_script.DEPTH_CAPS = (5.0, 10.0, 15.0, 20.0)
            matrix_script.build_matrix()
            depth_script = load_analysis_script(
                "02_depth_missingness_sensitivity.py", "analysis_depth_test"
            )
            depth_script.EVIDENCE_MATRIX_PATH = matrix_path
            depth_script.OUTPUT_DIR = second
            depth_script.MIN_CPG_GRID = (3, 5, 10)
            depth_script.PRIMARY_MIN_CPGS = 3
            depth_script.MIN_PWS_PARTICIPANTS = 3
            depth_script.MIN_AS_PARTICIPANTS = 2
            depth_script.run()
            regions_path = out / "prespecified_regions.tsv"
            pd.DataFrame(
                [
                    {
                        "region_id": "synthetic_domain",
                        "chrom": "chr15",
                        "start": 100,
                        "end": 900,
                    }
                ]
            ).to_csv(regions_path, sep="\t", index=False)
            architecture_script = load_analysis_script(
                "03_reciprocal_cis_architecture.py", "analysis_architecture_test"
            )
            architecture_script.EVIDENCE_MATRIX_PATH = matrix_path
            architecture_script.REGIONS_PATH = regions_path
            architecture_script.OUTPUT_DIR = third
            architecture_script.MIN_CPGS = 3
            architecture_script.MIN_PWS_PARTICIPANTS = 3
            architecture_script.MIN_AS_PARTICIPANTS = 2
            architecture_script.BOOTSTRAP_REPLICATES = 100
            architecture_script.RANDOM_SEED = 20260924
            architecture_script.FOCAL_EFFECT_THRESHOLD = 0.10
            architecture_script.EQUIVALENCE_MARGIN = 0.10
            architecture_script.run()
            architecture = pd.read_csv(
                out / "03_cis_architecture" / "cis_architecture_windows.tsv.gz", sep="\t"
            )
            np.testing.assert_allclose(architecture["delta_maternal_minus_paternal"], 0.6)
            np.testing.assert_allclose(architecture["scaffold_beta"], 0.5)
            np.testing.assert_allclose(architecture["absolute_asm_mean"], 0.6)


if __name__ == "__main__":
    unittest.main()
