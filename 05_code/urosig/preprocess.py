"""Within-cohort expression normalisation."""
from __future__ import annotations

import numpy as np
import pandas as pd


def quantile_normalise(expr: pd.DataFrame) -> pd.DataFrame:
    """Quantile normalisation across samples (columns); NaNs stay NaN and are excluded."""
    x = expr.to_numpy(float)
    n = x.shape[0]
    sorted_cols = []
    for j in range(x.shape[1]):
        v = np.sort(x[~np.isnan(x[:, j]), j])
        sorted_cols.append(np.interp(np.linspace(0, 1, n), np.linspace(0, 1, len(v)), v))
    ref = np.mean(sorted_cols, axis=0)
    out = np.full_like(x, np.nan)
    for j in range(x.shape[1]):
        ok = ~np.isnan(x[:, j])
        r = pd.Series(x[ok, j]).rank(method="average").to_numpy()
        pos = (r - 1) / max(ok.sum() - 1, 1)
        out[ok, j] = np.interp(pos, np.linspace(0, 1, n), ref)
    return pd.DataFrame(out, index=expr.index, columns=expr.columns)


def fpkm_to_log2tpm(fpkm: pd.DataFrame) -> pd.DataFrame:
    tpm = fpkm.div(fpkm.sum(axis=0), axis=1) * 1e6
    return np.log2(tpm + 1)
