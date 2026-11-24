import numpy as np
import pandas as pd
import pytest
from lifelines.utils import concordance_index
from urosig.metrics import harrell_c, incremental_value, hr_per_sd

def sim(n=400, beta_clin=0.7, beta_sig=0.7, seed=1):
    rng = np.random.default_rng(seed)
    clin = rng.normal(size=n); sig = rng.normal(size=n)
    lam = np.exp(beta_clin * clin + beta_sig * sig)
    t = rng.exponential(1 / lam); c = rng.exponential(2, size=n)
    time = np.minimum(t, c); event = (t <= c).astype(int)
    return time, event, pd.DataFrame({"stage": clin}), pd.Series(sig)

def test_harrell_c_perfect_and_reversed():
    time = np.array([1, 2, 3, 4, 5.0]); event = np.ones(5, int)
    assert harrell_c(time, event, -time) == pytest.approx(1.0)
    assert harrell_c(time, event, time) == pytest.approx(0.0)

def test_harrell_c_matches_lifelines_high_score_high_risk():
    time, event, _, sig = sim()
    assert harrell_c(time, event, sig.values) == pytest.approx(concordance_index(time, -sig.values, event), abs=1e-9)

def test_incremental_value_positive_when_signature_informative():
    time, event, clin, sig = sim(beta_sig=0.8)
    r = incremental_value(time, event, sig, clin, n_boot=200, seed=0)
    assert r["delta_c"] > 0.03
    assert r["delta_c_lo"] > 0
    assert r["lr_p"] < 0.001

def test_incremental_value_null_when_signature_is_noise():
    time, event, clin, _ = sim(beta_sig=0.0, seed=3)
    noise = pd.Series(np.random.default_rng(9).normal(size=len(time)))
    r = incremental_value(time, event, noise, clin, n_boot=200, seed=0)
    assert r["delta_c_lo"] < 0 < r["delta_c_hi"] or abs(r["delta_c"]) < 0.01

def test_hr_per_sd_recovers_effect_direction():
    time, event, _, sig = sim(beta_clin=0.0, beta_sig=0.6, n=1500)
    r = hr_per_sd(time, event, sig)
    assert 1.5 < r["hr"] < 2.2 and r["lo"] > 1
