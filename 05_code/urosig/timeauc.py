"""IPCW discrimination measures in numpy (no scikit-survival on the analysis machine).

Uno's C (truncated at tau) and the cumulative/dynamic time-dependent AUC, both with inverse-probability-of-
censoring weights from the Kaplan-Meier estimate of the censoring distribution in the same cohort.
Written to match scikit-survival 0.28 (concordance_index_ipcw, cumulative_dynamic_auc); the test in
05_code/tests/test_timeauc.py checks the agreement on simulated data.
Convention: higher score = higher risk.
"""
from __future__ import annotations

import numpy as np


def censoring_km(time, event):
    """Kaplan-Meier of the censoring distribution: returns (unique times, G(t) right after each time)."""
    t = np.asarray(time, float)
    c = 1 - np.asarray(event, int)  # censoring indicator
    ut = np.unique(t)
    ev = 1 - c
    # as scikit-survival's reverse KM: at a tied time, events leave the risk set before censorings
    at_risk = np.array([(t >= u).sum() - ev[t == u].sum() for u in ut], float)
    d = np.array([c[t == u].sum() for u in ut], float)
    ratio = np.divide(d, at_risk, out=np.zeros_like(d), where=d != 0)
    g = np.cumprod(1 - ratio)
    return ut, g


def _g_at(ut, g, x):
    """G(x): KM survival of censoring evaluated at x (step function, right-continuous)."""
    idx = np.searchsorted(ut, x, side="right") - 1
    out = np.where(idx >= 0, g[np.clip(idx, 0, None)], 1.0)
    return out


def uno_c(time, event, score, tau: float | None = None) -> float:
    t = np.asarray(time, float); e = np.asarray(event, int); s = np.asarray(score, float)
    ut, g = censoring_km(t, e)
    if tau is None:
        tau = t.max()
    gi = _g_at(ut, g, t)
    w = np.where((e == 1) & (t < tau) & (gi > 0), 1.0 / np.square(np.where(gi > 0, gi, 1)), 0.0)
    num = den = 0.0
    for i in np.flatnonzero(w > 0):
        comp = (t > t[i]) | ((t == t[i]) & (e == 0))  # a censored case at the same time is comparable
        if not comp.any():
            continue
        si, sj = s[i], s[comp]
        conc = (si > sj).sum() + 0.5 * (si == sj).sum()
        num += w[i] * conc
        den += w[i] * comp.sum()
    return float(num / den) if den > 0 else float("nan")


def td_auc(time, event, score, at: float) -> float:
    """Cumulative/dynamic AUC at time `at`: cases T<=at with event, controls T>at, cases weighted 1/G(T)."""
    t = np.asarray(time, float); e = np.asarray(event, int); s = np.asarray(score, float)
    ut, g = censoring_km(t, e)
    case = (t <= at) & (e == 1)
    ctrl = t > at
    if case.sum() == 0 or ctrl.sum() == 0:
        return float("nan")
    gi = _g_at(ut, g, t[case])
    w = np.where(gi > 0, 1.0 / np.where(gi > 0, gi, 1), 0.0)
    sc, sk = s[case], s[ctrl]
    conc = (sc[:, None] > sk[None, :]).sum(1) + 0.5 * (sc[:, None] == sk[None, :]).sum(1)
    return float((w * conc).sum() / (w.sum() * ctrl.sum()))
