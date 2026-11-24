import numpy as np
import pandas as pd
from urosig.nullsets import matched_random_sets, orient_genes

def test_random_sets_same_size_exclude_signature_and_match_deciles():
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(2000)]
    mean_expr = pd.Series(rng.uniform(0, 15, 2000), index=genes)
    sig = ["G1", "G5", "G10", "G500"]
    sets = matched_random_sets(sig, mean_expr, n_sets=50, seed=1)
    bins = pd.qcut(mean_expr, 10, labels=False)
    target = sorted(bins[sig])
    assert len(sets) == 50
    for s in sets:
        assert len(s) == 4 and len(set(s)) == 4
        assert not set(s) & set(sig)
        assert sorted(bins[s]) == target

def test_random_sets_reproducible_with_seed():
    mean_expr = pd.Series(np.linspace(0, 10, 300), index=[f"G{i}" for i in range(300)])
    a = matched_random_sets(["G1", "G2"], mean_expr, 5, seed=7)
    b = matched_random_sets(["G1", "G2"], mean_expr, 5, seed=7)
    assert a == b

def test_orient_genes_recovers_direction():
    rng = np.random.default_rng(2)
    n = 600
    expr = pd.DataFrame(rng.normal(size=(3, n)), index=["risk", "protect", "null"])
    lam = np.exp(0.8 * expr.loc["risk"] - 0.8 * expr.loc["protect"])
    t = rng.exponential(1 / lam.values); c = rng.exponential(2, n)
    T = np.minimum(t, c); E = (t <= c).astype(int)
    s = orient_genes(expr, T, E, ["risk", "protect"])
    assert s["risk"] == 1 and s["protect"] == -1
