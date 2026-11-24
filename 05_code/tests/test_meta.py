import numpy as np
import pytest
from scipy import optimize, stats
from urosig.meta import random_effects, logit_c

def test_homogeneous_studies_give_zero_tau2_and_iv_mean():
    yi = np.array([0.3, 0.3, 0.3, 0.3]); vi = np.array([0.01, 0.02, 0.03, 0.04])
    r = random_effects(yi, vi)
    assert r["tau2"] == pytest.approx(0, abs=1e-8)
    assert r["est"] == pytest.approx(0.3)

def test_reml_tau2_maximizes_restricted_likelihood():
    yi = np.array([0.10, 0.45, -0.05, 0.60, 0.30, 0.20]); vi = np.array([0.02, 0.03, 0.015, 0.04, 0.02, 0.05])
    def nrll(t2):
        w = 1 / (vi + t2); mu = (w * yi).sum() / w.sum()
        return 0.5 * (np.log(vi + t2).sum() + np.log(w.sum()) + (w * (yi - mu) ** 2).sum())
    ref = optimize.minimize_scalar(nrll, bounds=(0, 5), method="bounded", options={"xatol": 1e-10}).x
    r = random_effects(yi, vi)
    assert r["tau2"] == pytest.approx(ref, abs=1e-6)

def test_hartung_knapp_uses_t_with_k_minus_1_df():
    yi = np.array([0.10, 0.45, -0.05, 0.60, 0.30, 0.20]); vi = np.array([0.02, 0.03, 0.015, 0.04, 0.02, 0.05])
    r = random_effects(yi, vi, hk=True)
    w = 1 / (vi + r["tau2"]); mu = (w * yi).sum() / w.sum()
    q = (w * (yi - mu) ** 2).sum() / (len(yi) - 1)
    se = np.sqrt(q / w.sum())
    crit = stats.t.ppf(0.975, len(yi) - 1)
    assert r["lo"] == pytest.approx(mu - crit * se) and r["hi"] == pytest.approx(mu + crit * se)

def test_i2_and_prediction_interval_present():
    yi = np.array([0.10, 0.45, -0.05, 0.60]); vi = np.array([0.02, 0.03, 0.015, 0.04])
    r = random_effects(yi, vi)
    assert 0 <= r["i2"] <= 1 and r["pi_lo"] < r["lo"] and r["pi_hi"] > r["hi"]

def test_logit_c_delta_method():
    y, v = logit_c(0.7, 0.02)
    assert y == pytest.approx(np.log(0.7 / 0.3))
    assert v == pytest.approx(0.02 ** 2 / (0.7 * 0.3) ** 2)
