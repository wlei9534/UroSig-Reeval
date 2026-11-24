"""Random-effects pooling (REML tau², optional Hartung–Knapp CI)."""
from __future__ import annotations

import numpy as np
from scipy import stats


def logit_c(c: float, se_c: float):
    """Logit-transformed C and its delta-method variance."""
    return float(np.log(c / (1 - c))), float(se_c ** 2 / (c * (1 - c)) ** 2)


def inv_logit(x):
    return 1 / (1 + np.exp(-np.asarray(x, float)))


def _reml_tau2(yi, vi, max_iter=200, tol=1e-12):
    t2 = max(0.0, np.var(yi, ddof=1) - vi.mean())
    for _ in range(max_iter):
        w = 1 / (vi + t2)
        mu = (w * yi).sum() / w.sum()
        # Fisher scoring step on the restricted likelihood
        p_resid = w * (yi - mu)
        num = (p_resid ** 2).sum() - (w.sum() - (w ** 2).sum() / w.sum())
        den = (w ** 2).sum() - 2 * (w ** 3).sum() / w.sum() + ((w ** 2).sum() / w.sum()) ** 2
        new = max(0.0, t2 + num / den)
        if abs(new - t2) < tol:
            t2 = new
            break
        t2 = new
    return t2


def random_effects(yi, vi, hk: bool = True, level: float = 0.95) -> dict:
    yi = np.asarray(yi, float); vi = np.asarray(vi, float)
    k = len(yi)
    t2 = _reml_tau2(yi, vi) if k > 1 else 0.0
    w = 1 / (vi + t2)
    mu = (w * yi).sum() / w.sum()
    wf = 1 / vi
    muf = (wf * yi).sum() / wf.sum()
    Q = (wf * (yi - muf) ** 2).sum()
    i2 = max(0.0, (Q - (k - 1)) / Q) if Q > 0 and k > 1 else 0.0
    a = 1 - (1 - level) / 2
    if hk and k > 1:
        q = (w * (yi - mu) ** 2).sum() / (k - 1)
        se = np.sqrt(q / w.sum())
        crit = stats.t.ppf(a, k - 1)
    else:
        se = np.sqrt(1 / w.sum())
        crit = stats.norm.ppf(a)
    pi_crit = stats.t.ppf(a, k - 2) if k > 2 else np.nan
    pi_se = np.sqrt(t2 + se ** 2)
    return {"est": float(mu), "se": float(se), "lo": float(mu - crit * se), "hi": float(mu + crit * se),
            "tau2": float(t2), "i2": float(i2), "Q": float(Q), "k": k,
            "pi_lo": float(mu - pi_crit * pi_se), "pi_hi": float(mu + pi_crit * pi_se)}
