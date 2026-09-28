"""
Helpers for the 15q11-q13 duplicon, breakpoint and repeat analyses
(scripts/duplicons/01-07). Cohort metadata and CN=1 intervals are read with the
repository's own library (scripts/analysis/cis_analysis) so the evidence rules
stay in one place.
"""
from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[2] / "analysis"


def cis_analysis():
    """Import scripts/analysis/cis_analysis (load_cohort, load_deletion_map)."""
    if str(ANALYSIS_DIR) not in sys.path:
        sys.path.insert(0, str(ANALYSIS_DIR))
    import cis_analysis as module
    return module
