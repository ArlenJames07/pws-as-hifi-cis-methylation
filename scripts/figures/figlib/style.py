"""Matplotlib settings, colours, panel labels, saving, tables and markdown reports.

Maternal and paternal copies use the Figure 1 state colours (vermillion / blue);
groups use the Figure 1 group colours and markers (cohort.py). Neutral grey marks
genomes with both parental copies. Figures are built at 7.09 in (double column)."""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

MATERNAL = "#D55E00"
PATERNAL = "#0072B2"
BIPARENTAL = "#7F7F7F"
INK = "#222222"
MUTED = "#6B6B6B"
LIGHT = "#BDBDBD"
GRID = "#E6E6E6"
TRACE = "#6A2C91"          # smoothed differences
CORE = "#C49A00"           # shared boundary core / highlighted interval
SD = "#CDBB95"             # segmental duplications
DOMAIN = "#F3E9C6"         # shading of a highlighted genomic domain
BASE_FONT = 7.0
WIDTH = 7.09               # inches, double column
FORMATS = ("png", "pdf", "svg")
DPI = 600
FONT_SANS = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]


def setup() -> None:
    plt.rcParams.update({
        # Figure 1 uses the same Arial-compatible stack.  Keeping one ordered
        # fallback list also makes the actual face deterministic on machines
        # without a licensed Arial installation (Liberation Sans is preferred
        # before Matplotlib's DejaVu fallback).
        "font.family": "sans-serif", "font.sans-serif": FONT_SANS,
        "font.size": BASE_FONT, "axes.titlesize": BASE_FONT + 0.5,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Liberation Sans",
        "mathtext.it": "Liberation Sans:italic",
        "mathtext.bf": "Liberation Sans:bold",
        "axes.labelsize": BASE_FONT, "xtick.labelsize": BASE_FONT - 0.5, "ytick.labelsize": BASE_FONT - 0.5,
        "legend.fontsize": BASE_FONT - 0.5, "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
        "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK, "legend.frameon": False,
        "pdf.fonttype": 42, "svg.fonttype": "none", "figure.facecolor": "white", "axes.facecolor": "white",
    })


def mb(x):
    return np.asarray(x, dtype=float) / 1e6 if np.ndim(x) else float(x) / 1e6


def panel_label(fig, ax, letter: str, title: str = "", dx: float = -0.06, dy: float = 0.012) -> None:
    """Bold letter and a short title above the top-left corner of `ax`, in figure coordinates."""
    box = ax.get_position()
    x = max(0.004, box.x0 + dx)
    fig.text(x, box.y1 + dy, letter, fontsize=BASE_FONT + 3, fontweight="bold", va="bottom", ha="left")
    if title:
        fig.text(x + 0.025, box.y1 + dy + 0.002, title, fontsize=BASE_FONT + 0.5, va="bottom", ha="left")


def save(fig, outdir: Path, stem: str, formats=FORMATS, dpi: int = DPI) -> list[Path]:
    fig.canvas.draw()
    width = fig.get_tightbbox(fig.canvas.get_renderer()).width
    if width > WIDTH + 0.3:
        warnings.warn(f"{stem}: exported width {width:.2f} in exceeds {WIDTH} in")
    paths = []
    for fmt in formats:
        path = Path(outdir) / "figures" / f"{stem}.{fmt}"
        fig.savefig(path, dpi=dpi if fmt == "png" else None, bbox_inches="tight", pad_inches=0.04)
        paths.append(path)
    plt.close(fig)
    return paths


def write_table(df: pd.DataFrame, outdir: Path, name: str) -> Path:
    """TSV (gzip if the name ends in .gz). Float columns holding whole numbers (coordinates
    that became float through NaN) are written as integers, other floats with 8 digits."""
    path = Path(outdir) / "tables" / name
    df = df.copy()
    for col in df.columns:
        if pd.api.types.is_float_dtype(df[col]):
            v = df[col].to_numpy(dtype=float)
            fin = np.isfinite(v)
            if fin.any() and np.all(np.mod(v[fin], 1) == 0) and np.abs(v[fin]).max() < 2 ** 53:
                df[col] = df[col].astype("Int64")
    df.to_csv(path, sep="\t", index=False, float_format="%.8g", na_rep="NA")
    return path


def read_table(outdir: Path, name: str) -> pd.DataFrame:
    path = Path(outdir) / "tables" / name
    if not path.is_file():
        raise SystemExit(f"--render-only: {path} is missing; run without --render-only first")
    return pd.read_csv(path, sep="\t", low_memory=False)


def fmt(value, digits: int = 3, signed: bool = False) -> str:
    try:
        v = float(value)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(v):
        return "NA"
    return f"{v:+.{digits}f}" if signed else f"{v:.{digits}f}"


def fmt_p(p) -> str:
    try:
        p = float(p)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(p):
        return "NA"
    return f"{p:.1e}" if p < 1e-3 else f"{p:.3f}"


def _cell(v) -> str:
    if isinstance(v, (bool, np.bool_)):
        return "yes" if v else "no"
    if isinstance(v, (int, np.integer)):
        return f"{int(v):,}" if abs(int(v)) >= 10_000 else str(int(v))
    if isinstance(v, (float, np.floating)):
        if not np.isfinite(v):
            return "NA"
        if float(v).is_integer():
            return f"{int(v):,}" if abs(v) >= 10_000 else str(int(v))
        return f"{v:.4g}"
    return str(v).replace("|", "\\|")


def md_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_(empty)_"
    cols = list(df.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for row in df.itertuples(index=False):
        lines.append("| " + " | ".join(_cell(v) for v in row) + " |")
    return "\n".join(lines)


class Report:
    """Markdown report written next to the figure: inputs, numbers behind every claim,
    a caption draft. Nothing in it is hard-coded; every value comes from the tables."""

    def __init__(self, title: str) -> None:
        self.lines = [f"# {title}", ""]

    def h(self, text: str) -> None:
        self.lines += ["", f"## {text}", ""]

    def p(self, text: str = "") -> None:
        self.lines.append(text)

    def table(self, df: pd.DataFrame) -> None:
        self.lines += [md_table(df), ""]

    def write(self, outdir: Path, name: str) -> Path:
        path = Path(outdir) / "reports" / name
        path.write_text("\n".join(self.lines) + "\n")
        return path


def group_legend_handles(groups, counts: dict | None = None, filled=True):
    from matplotlib.lines import Line2D

    from .cohort import COLOR, DISPLAY, MARKER
    handles = []
    for g in groups:
        n = f" (n = {counts[g]})" if counts and g in counts else ""
        handles.append(Line2D([], [], marker=MARKER[g], ls="", ms=4.5, mec=COLOR[g],
                              mfc=COLOR[g] if filled else "white", mew=0.9, label=f"{DISPLAY[g]}{n}"))
    return handles
