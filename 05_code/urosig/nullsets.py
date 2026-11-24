"""Expression-matched random gene sets: the 'random signature' baseline."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .cox import cox_fit


def matched_random_sets(sig_genes, mean_expr: pd.Series, n_sets: int, seed: int = 0, n_bins: int = 10):
    """Random gene sets of the signature's size, each gene drawn from the same
    mean-expression decile as the signature gene it replaces. Signature genes excluded."""
    sig = [g for g in sig_genes if g in mean_expr.index]
    bins = pd.qcut(mean_expr.rank(method="first"), n_bins, labels=False)
    pools = {b: np.array([g for g in bins.index[bins == b] if g not in set(sig)]) for b in range(n_bins)}
    rng = np.random.default_rng(seed)
    need = pd.Series(bins[sig].values).value_counts().to_dict()
    out = []
    for _ in range(n_sets):
        s = []
        for b, k in sorted(need.items()):
            s.extend(rng.choice(pools[b], size=k, replace=False).tolist())
        out.append(s)
    return out


def orient_genes(expr: pd.DataFrame, T, E, genes) -> pd.Series:
    """Sign of each gene's univariate Cox coefficient in the development cohort (+1 = risk)."""
    signs = {}
    for g in genes:
        x = expr.loc[g].values.astype(float)
        x = (x - x.mean()) / (x.std(ddof=1) or 1.0)
        beta, _ = cox_fit(T, E, x[:, None])
        signs[g] = 1 if beta[0] >= 0 else -1
    return pd.Series(signs)
