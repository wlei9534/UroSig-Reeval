import numpy as np
import pandas as pd
from urosig.preprocess import quantile_normalise, fpkm_to_log2tpm

def test_quantile_normalise_gives_identical_column_distributions():
    rng = np.random.default_rng(0)
    e = pd.DataFrame(rng.normal(size=(200, 4)) * [1, 2, 3, 4] + [0, 1, 2, 3])
    q = quantile_normalise(e)
    s = np.sort(q.values, axis=0)
    np.testing.assert_allclose(s[:, 0], s[:, 3])
    assert (q.rank().values == e.rank().values).all()

def test_quantile_normalise_ignores_missing_values():
    e = pd.DataFrame({"a": [1.0, 2.0, np.nan, 4.0], "b": [10.0, 20.0, 30.0, 40.0]})
    q = quantile_normalise(e)
    assert np.isnan(q.loc[2, "a"]) and q["a"].notna().sum() == 3

def test_fpkm_to_log2tpm_columns_sum_to_million_before_log():
    f = pd.DataFrame({"s1": [1.0, 3.0, 0.0], "s2": [2.0, 2.0, 4.0]})
    out = fpkm_to_log2tpm(f)
    np.testing.assert_allclose((2 ** out - 1).sum(0).values, [1e6, 1e6])
