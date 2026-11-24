"""Discrimination and incremental value. Convention: higher score = higher risk."""
from __future__ import annotations

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from scipy import stats

from .cox import cox_lp


def harrell_c(time, event, score) -> float:
    """Harrell's C where a higher score means shorter survival."""
    return float(concordance_index(np.asarray(time), -np.asarray(score, float), np.asarray(event)))


def _fit_lp(df: pd.DataFrame, cols: list[str]):
    lp, ll = cox_lp(df["T"].values, df["E"].values, df[cols].values)
    return ll, lp


def _frame(time, event, score, clin):
    df = pd.DataFrame(np.asarray(clin, float) if not isinstance(clin, pd.DataFrame) else clin.values,
                      columns=list(clin.columns) if isinstance(clin, pd.DataFrame) else None)
    clin_cols = [f"clin_{c}" for c in df.columns]
    df.columns = clin_cols
    df["sig"] = np.asarray(score, float)
    df["T"] = np.asarray(time, float)
    df["E"] = np.asarray(event, int)
    return df.dropna().reset_index(drop=True), clin_cols


def incremental_value(time, event, score, clin: pd.DataFrame, n_boot: int = 1000, seed: int = 0) -> dict:
    """ΔC of (clinical + score) over clinical, with percentile bootstrap CI and LR test.

    The score enters as ONE covariate; its internal coefficients stay fixed.
    Apparent ΔC in the external cohort (both models refit there); bootstrap resamples patients.
    """
    df, cc = _frame(time, event, score, clin)
    m0, lp0 = _fit_lp(df, cc)
    m1, lp1 = _fit_lp(df, cc + ["sig"])
    c0 = harrell_c(df["T"], df["E"], lp0)
    c1 = harrell_c(df["T"], df["E"], lp1)
    lr = 2 * (m1 - m0)
    p = float(stats.chi2.sf(lr, 1))
    rng = np.random.default_rng(seed)
    deltas = []
    n = len(df)
    for _ in range(n_boot):
        b = df.iloc[rng.integers(0, n, n)].reset_index(drop=True)
        if b["E"].sum() < 2:
            continue
        try:
            _, b0 = _fit_lp(b, cc)
            _, b1 = _fit_lp(b, cc + ["sig"])
        except Exception:
            continue
        deltas.append(harrell_c(b["T"], b["E"], b1) - harrell_c(b["T"], b["E"], b0))
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    return {"c_clin": c0, "c_full": c1, "delta_c": c1 - c0, "delta_c_lo": float(lo),
            "delta_c_hi": float(hi), "lr_chi2": float(lr), "lr_p": p, "n": n,
            "events": int(df["E"].sum()), "n_boot_ok": len(deltas)}


def hr_per_sd(time, event, score) -> dict:
    s = np.asarray(score, float)
    df = pd.DataFrame({"z": (s - s.mean()) / s.std(ddof=1), "T": time, "E": event}).dropna()
    cph = CoxPHFitter().fit(df, duration_col="T", event_col="E")
    b, se = cph.params_["z"], cph.standard_errors_["z"]
    return {"hr": float(np.exp(b)), "lo": float(np.exp(b - 1.96 * se)), "hi": float(np.exp(b + 1.96 * se)),
            "p": float(cph.summary.loc["z", "p"])}
