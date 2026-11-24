import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from urosig.cox import cox_fit

def test_cox_fit_matches_lifelines_breslow():
    rng = np.random.default_rng(4)
    n = 300
    X = rng.normal(size=(n, 3))
    t = rng.exponential(1 / np.exp(X @ [0.5, -0.3, 0.0]))
    t = np.round(t, 1) + 0.1          # induce ties
    c = rng.exponential(2, n)
    T = np.minimum(t, c); E = (t <= c).astype(int)
    beta, loglik = cox_fit(T, E, X)
    df = pd.DataFrame(X, columns=list("abc")); df["T"] = T; df["E"] = E
    ref = CoxPHFitter().fit(df, "T", "E")
    np.testing.assert_allclose(beta, ref.params_.values, atol=1e-4)
    np.testing.assert_allclose(loglik, ref.log_likelihood_, atol=1e-4)
