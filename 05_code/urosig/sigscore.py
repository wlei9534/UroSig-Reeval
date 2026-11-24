"""Risk scores for locked signatures in one cohort (SAP v2 section 5). No outcome data is touched here.

Signature kinds
  linear      sum(beta_g * x_g)                       roles main / sens_figure / sens_hr / sens_lasso
  pair        sum(beta_p * I(rule(A, B)))             gene-pair signatures (HGNC symbol 'A|B')
  equal       mean(direction_g * z_g)                 roles equal_weight / equal_signed
Scaling of x_g (linear kind)
  'z'         within-cohort z-score (primary)
  'raw'       log2 expression as processed (sensitivity)
  'tcga'      within-cohort z-score x SD of the gene in TCGA-PRAD (sensitivity; TCGA-trained signatures only)
Missing genes contribute 0 and the remaining weights are not rescaled. Higher score = higher risk.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

# pair rule per paper: 'gt' = 1 if A > B; 'lt' = 1 if A < B (locked extraction notes)
PAIR_RULE = {"42553741": "gt", "35143414": "gt", "32593277": "lt",
             "31961962": "gt", "38191471": "gt", "34876840": "gt"}  # v2: 74-GPS (Table 2 note), k-TSP (inferred), lncRNA pairs (Methods)


@dataclass
class Score:
    sig_id: str
    cohort: str
    mode: str
    values: pd.Series | None
    n_terms: int
    n_used: int
    gene_cov: float
    weight_cov: float | None
    kind: str
    note: str = ""
    missing: list = field(default_factory=list)


def _z(x: pd.DataFrame) -> pd.DataFrame:
    sd = x.std(axis=1, ddof=1).replace(0, np.nan)
    return x.sub(x.mean(axis=1), axis=0).div(sd, axis=0).fillna(0.0)


def score_signature(sig_id: str, cohort: str, terms: pd.DataFrame, expr: pd.DataFrame, mode: str = "z",
                    directions: dict | None = None, tcga_sd: pd.Series | None = None) -> Score:
    """terms: rows of signature_terms_locked for this signature; expr: genes x samples, HGNC symbols."""
    role = terms["分析角色"].iloc[0]
    syms = terms["HGNC 现行符号"].tolist()
    idx = expr.index
    if any("|" in s for s in syms):
        kind = "pair"
    elif role in ("equal_weight", "equal_signed"):
        kind = "equal"
    else:
        kind = "linear"

    if kind == "pair":
        rule = PAIR_RULE[sig_id.split("_")[0]]
        coefs = terms["系数"].astype(float).to_numpy()
        ok = np.array([all(p in idx for p in s.split("|")) for s in syms])
        total = np.abs(coefs).sum()
        val = pd.Series(0.0, index=expr.columns)
        for s, b, k in zip(syms, coefs, ok):
            if not k:
                continue
            a, c = s.split("|")
            ind = (expr.loc[a] > expr.loc[c]) if rule == "gt" else (expr.loc[a] < expr.loc[c])
            val += b * ind.astype(float)
        miss = [s for s, k in zip(syms, ok) if not k]
        return Score(sig_id, cohort, mode, val, len(syms), int(ok.sum()), float(ok.mean()),
                     float(np.abs(coefs[ok]).sum() / total), kind, f"pair rule {rule}", miss)

    if kind == "equal":
        if role == "equal_signed":
            dirs = terms["系数"].astype(float).to_numpy()
        else:
            d = directions or {}
            dirs = np.array([d.get(g, 0) for g in terms["基因（原文写法）"]], float)
        use = np.array([s in idx and dv != 0 for s, dv in zip(syms, dirs)])
        found = np.array([s in idx for s in syms])
        if use.sum() == 0:
            return Score(sig_id, cohort, mode, None, len(syms), 0, float(found.mean()), None, kind, "no directional gene found")
        z = _z(expr.loc[[s for s, u in zip(syms, use) if u]].astype(float))
        val = (z.T * dirs[use]).mean(axis=1)
        return Score(sig_id, cohort, mode, val, len(syms), int(use.sum()), float(found.mean()), None, kind,
                     f"{int((dirs == 0).sum())} genes without direction dropped", [s for s, f in zip(syms, found) if not f])

    coefs = terms["系数"].astype(float).to_numpy()
    found = np.array([s in idx for s in syms])
    total = np.abs(coefs).sum()
    kept = [s for s, f in zip(syms, found) if f]
    x = expr.loc[kept].astype(float)
    note = ""
    if mode == "z":
        x = _z(x)
    elif mode == "tcga":
        if tcga_sd is None:
            raise ValueError("tcga mode needs tcga_sd")
        sd = tcga_sd.reindex(kept)
        if sd.isna().any():
            note = f"no TCGA SD for {sd.isna().sum()} genes (set to cohort-z scale)"
        x = _z(x).mul(sd.fillna(1.0), axis=0)
    elif mode == "raw":
        # 23 missing cells in 7 genes (GSE70769/GSE70768): gene's cohort mean, as z = 0 does in the primary analysis
        x = x.T.fillna(x.mean(axis=1)).T
    else:
        raise ValueError(mode)
    val = (x.T * coefs[found]).sum(axis=1)
    return Score(sig_id, cohort, mode, val, len(syms), int(found.sum()), float(found.mean()),
                 float(np.abs(coefs[found]).sum() / total) if total else None, kind, note,
                 [s for s, f in zip(syms, found) if not f])
