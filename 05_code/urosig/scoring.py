"""Signature scores computed exactly as published: fixed genes and coefficients, no refitting."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd


@dataclass
class ScoreResult:
    evaluable: bool
    score: pd.Series | None = None
    missing: list[str] = field(default_factory=list)
    reason: str | None = None
    frac_missing: float = 0.0
    weight_missing: float = 0.0


def score_signature(
    expr: pd.DataFrame,
    genes: list[str],
    coefs: list[float],
    zscore: bool = True,
    max_missing_frac: float = 0.2,
    max_missing_weight: float = 0.2,
) -> ScoreResult:
    """Linear signature score.

    expr: genes x samples (log scale), index = HGNC symbols.
    Missing genes are dropped and their weight is NOT redistributed. The cohort is
    not evaluable if > max_missing_frac of genes, or > max_missing_weight of total
    |coef|, are missing.
    """
    if len(genes) != len(coefs):
        raise ValueError("genes and coefs differ in length")
    if len(set(genes)) != len(genes):
        raise ValueError("duplicate gene symbols in signature")
    coefs = np.asarray(coefs, dtype=float)
    present = np.array([g in expr.index for g in genes])
    missing = [g for g, p in zip(genes, present) if not p]
    frac = 1 - present.mean()
    wfrac = np.abs(coefs[~present]).sum() / np.abs(coefs).sum() if np.abs(coefs).sum() > 0 else 0.0
    if frac > max_missing_frac + 1e-12:
        return ScoreResult(False, None, missing, "missing_genes", frac, wfrac)
    if wfrac > max_missing_weight + 1e-12:
        return ScoreResult(False, None, missing, "missing_weight", frac, wfrac)
    kept = [g for g, p in zip(genes, present) if p]
    x = expr.loc[kept].astype(float)
    if zscore:
        sd = x.std(axis=1, ddof=1).replace(0, np.nan)
        x = x.sub(x.mean(axis=1), axis=0).div(sd, axis=0).fillna(0.0)
    score = (x.T * coefs[present]).sum(axis=1)
    score.name = "score"
    return ScoreResult(True, score, missing, None, frac, wfrac)
