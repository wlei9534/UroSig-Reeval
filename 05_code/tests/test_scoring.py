import numpy as np
import pandas as pd
import pytest
from urosig.scoring import score_signature

@pytest.fixture
def expr():
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(10)]
    return pd.DataFrame(rng.normal(5, 2, size=(10, 30)), index=genes,
                        columns=[f"S{j}" for j in range(30)])

def test_zscore_weighted_sum_matches_manual(expr):
    genes, coefs = ["G1", "G2", "G3"], [0.5, -1.0, 2.0]
    res = score_signature(expr, genes, coefs, zscore=True)
    z = expr.loc[genes].sub(expr.loc[genes].mean(1), axis=0).div(expr.loc[genes].std(1, ddof=1), axis=0)
    manual = (z.T * np.array(coefs)).sum(1)
    assert res.evaluable
    np.testing.assert_allclose(res.score.values, manual.values)

def test_raw_scale_when_zscore_false(expr):
    res = score_signature(expr, ["G1", "G2"], [1.0, 1.0], zscore=False)
    np.testing.assert_allclose(res.score.values, (expr.loc["G1"] + expr.loc["G2"]).values)

def test_missing_genes_dropped_without_reweighting(expr):
    genes, coefs = ["G1", "G2", "G3", "G4", "G5", "MISSING"], [1, 1, 1, 1, 1, 0.1]
    res = score_signature(expr, genes, coefs, zscore=False)
    assert res.evaluable
    assert res.missing == ["MISSING"]
    np.testing.assert_allclose(res.score.values, expr.loc[["G1","G2","G3","G4","G5"]].sum(0).values)

def test_not_evaluable_when_too_many_genes_missing(expr):
    res = score_signature(expr, ["G1", "G2", "G3", "X1"], [1, 1, 1, 1], max_missing_frac=0.2)
    assert not res.evaluable and res.reason == "missing_genes"

def test_not_evaluable_when_missing_weight_too_large(expr):
    genes = ["G1","G2","G3","G4","G5","G6","G7","G8","G9","X1"]
    coefs = [0.1]*9 + [5.0]
    res = score_signature(expr, genes, coefs, max_missing_frac=0.2, max_missing_weight=0.2)
    assert not res.evaluable and res.reason == "missing_weight"

def test_duplicate_gene_symbols_rejected(expr):
    with pytest.raises(ValueError):
        score_signature(expr, ["G1", "G1"], [1, 2])
