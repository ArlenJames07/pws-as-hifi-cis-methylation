from __future__ import annotations

import numpy as np
import pandas as pd


def fixed_windows(chrom: str, start: int, end: int, width: int) -> pd.DataFrame:
    if start < 0 or end <= start or width <= 0:
        raise ValueError("Invalid genomic interval or window width")
    starts = np.arange(start, end, width, dtype=np.int64)
    ends = np.minimum(starts + width, end)
    return pd.DataFrame(
        {
            "chrom": chrom,
            "window_id": [f"{chrom}:{left}-{right}" for left, right in zip(starts, ends)],
            "start": starts,
            "end": ends,
            "mid": (starts + ends) / 2.0,
        }
    )
