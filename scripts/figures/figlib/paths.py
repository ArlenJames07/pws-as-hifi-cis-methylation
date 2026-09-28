"""Where the pipeline wrote its outputs, and where each figure is written.

Everything is derived from params.local.yml (or $PARAMS_FILE) and the results
directory, so the figure scripts need no hard-coded server paths:

  results/01_alignment            pbmm2 BAMs
  results/03_structural_variants  pbsv calls
  results/04_phasing              HiPhase BAMs (HP/PS tags) and phased VCFs
  results/05_cnv/<sample>/        HiFiCNV
  results/06_methylation/<sample> pb-CpG-tools tracks
  results/analysis                scripts/analysis (00-03)
  results/08_duplicons            scripts/duplicons (01-08)
  results/07_figures/figure_N     this package's outputs
"""
from __future__ import annotations

import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]
ANALYSIS_SCRIPTS = PROJECT_ROOT / "scripts" / "analysis"
DUPLICON_SCRIPTS = PROJECT_ROOT / "scripts" / "duplicons"


def read_params(path: str | Path | None = None) -> dict[str, str]:
    """Flat `key: value` pairs of the Nextflow params YAML (comments and nesting ignored)."""
    path = Path(path or os.environ.get("PARAMS_FILE", PROJECT_ROOT / "params.local.yml"))
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    out: dict[str, str] = {}
    if not path.is_file():
        return out
    for line in path.read_text().splitlines():
        m = re.match(r"^([A-Za-z0-9_]+)\s*:\s*(.*?)\s*$", line)
        if not m:
            continue
        value = re.sub(r"\s+#.*$", "", m.group(2)).strip().strip("'\"")
        if value and value not in ("|", ">", "null", "~"):
            out[m.group(1)] = value
    return out


def _as_path(value: str | Path | None, default: Path | None = None) -> Path | None:
    if value in (None, "", "null", "None"):
        return default
    p = Path(value).expanduser()
    return p if p.is_absolute() else PROJECT_ROOT / p


@dataclass
class Paths:
    results: Path
    params: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        r = self.results
        self.alignment = r / "01_alignment"
        self.sv_calls = r / "03_structural_variants" / "calls"
        self.phasing = r / "04_phasing"
        self.cnv = r / "05_cnv"
        self.methylation = r / "06_methylation"
        self.figures = r / "07_figures"
        self.analysis = r / "analysis"
        self.duplicons = _as_path(os.environ.get("DUPLICON_OUT"), r / "08_duplicons")
        p = self.params
        self.reference = _as_path(p.get("reference"))
        self.gtf = _as_path(p.get("gtf"))
        self.metadata = _as_path(p.get("metadata"), PROJECT_ROOT / "assets" / "metadata.csv")
        self.icr_bed = _as_path(p.get("icr_bed"))
        self.segdup_bed = _as_path(p.get("segdup_bed"))

    def figure_dir(self, number: int | str, outdir: str | Path | None = None) -> Path:
        out = _as_path(outdir) if outdir else self.figures / f"figure_{number}"
        for sub in ("figures", "tables", "reports"):
            (out / sub).mkdir(parents=True, exist_ok=True)
        return out

    def fai(self) -> Path | None:
        if self.reference is None:
            return None
        fai = Path(f"{self.reference}.fai")
        return fai if fai.is_file() else None


def add_common_arguments(parser) -> None:
    parser.add_argument("--results", help="pipeline results directory (default: params 'outdir' or results/)")
    parser.add_argument("--params", help="params YAML (default: $PARAMS_FILE or params.local.yml)")
    parser.add_argument("--outdir", help="output directory (default: results/07_figures/figure_N)")
    parser.add_argument("--metadata", help="cohort metadata CSV (default: params 'metadata')")
    parser.add_argument("--gtf", help="CHM13 GTF (default: params 'gtf')")


def resolve(args) -> Paths:
    params = read_params(getattr(args, "params", None))
    results = _as_path(getattr(args, "results", None) or os.environ.get("FIGURES_RESULTS") or params.get("outdir"),
                       PROJECT_ROOT / "results")
    paths = Paths(results=results, params=params)
    if getattr(args, "metadata", None):
        paths.metadata = _as_path(args.metadata)
    if getattr(args, "gtf", None):
        paths.gtf = _as_path(args.gtf)
    return paths


def import_analysis_package():
    """scripts/analysis/cis_analysis (cohort, deletion map, pb-CpG-tools reader)."""
    if str(ANALYSIS_SCRIPTS) not in sys.path:
        sys.path.insert(0, str(ANALYSIS_SCRIPTS))
    import cis_analysis  # noqa: E402

    return cis_analysis
