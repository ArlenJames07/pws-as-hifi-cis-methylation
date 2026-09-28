"""
Small helpers for calling the command-line tools from Python: logging, checks
that a program exists, running commands and pipelines, and resumable steps.
"""
from __future__ import annotations

import gzip
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Iterable, Sequence


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", file=sys.stderr, flush=True)


def require(*programs: str) -> None:
    missing = [p for p in programs if shutil.which(p) is None and not Path(p).is_file()]
    if missing:
        raise SystemExit(
            "Not found on PATH: " + ", ".join(missing)
            + ". Install them or set the *_bin entries in params.local.yml (see config.py)."
        )


def run(command: Sequence[str | Path], stdout: Path | None = None) -> None:
    """Run one program; raise on failure. stdout is written to a file if given."""
    command = [str(c) for c in command]
    log("$ " + " ".join(command) + (f" > {stdout}" if stdout else ""))
    if stdout is None:
        subprocess.run(command, check=True)
    else:
        tmp = Path(f"{stdout}.part")
        with open(tmp, "w") as handle:
            subprocess.run(command, check=True, stdout=handle)
        tmp.replace(stdout)


def run_pipeline(pipeline: str) -> None:
    """Run a shell pipeline under bash with pipefail, so a failing step is not hidden."""
    log("$ " + pipeline)
    subprocess.run(["bash", "-o", "pipefail", "-c", pipeline], check=True)


def up_to_date(outputs: Iterable[Path], force: bool = False) -> bool:
    """True when every output exists and is non-empty (the step can be skipped)."""
    if force:
        return False
    outputs = list(outputs)
    return bool(outputs) and all(Path(o).exists() and (Path(o).is_dir() or Path(o).stat().st_size > 0)
                                 for o in outputs)


def to_fasta(source: Path, target: Path) -> None:
    """hifiasm GFA (S lines) or FASTA, optionally gzipped -> plain FASTA."""
    opener = gzip.open if str(source).endswith(".gz") else open
    name = str(source).removesuffix(".gz")
    tmp = Path(f"{target}.part")
    with opener(source, "rt") as src, open(tmp, "w") as dst:
        if name.endswith(".gfa"):
            for line in src:
                if line.startswith("S\t"):
                    fields = line.rstrip("\n").split("\t")
                    dst.write(f">{fields[1]}\n{fields[2]}\n")
        else:
            shutil.copyfileobj(src, dst)
    tmp.replace(target)


def parse_samples_arg(value: str | None) -> set[str] | None:
    return {s.strip() for s in value.split(",") if s.strip()} if value else None
