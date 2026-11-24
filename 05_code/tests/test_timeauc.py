"""Agreement of urosig.timeauc with scikit-survival 0.28 on simulated cohorts with tied times.

Run where scikit-survival is installed (the cloud workspace): python tests/test_timeauc.py
Last run 2026-11-01: max absolute difference 1e-15 over 30 cohorts (Uno C at tau=60; td-AUC at 24/36/60).
"""
import sys
from pathlib import Path

import numpy as np
from sksurv.metrics import concordance_index_ipcw, cumulative_dynamic_auc
from sksurv.util import Surv

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from urosig.timeauc import td_auc, uno_c  # noqa: E402

rng = np.random.default_rng(1)
mx = 0.0
for _ in range(30):
    n = int(rng.integers(60, 300))
    x = rng.normal(size=n)
    T = rng.exponential(np.exp(-0.7 * x)) * 60
    C = rng.uniform(5, 120, n)
    t = np.round(np.minimum(T, C), 0)
    e = (T <= C).astype(int)
    s = np.round(x + rng.normal(scale=0.5, size=n), 1)
    y = Surv.from_arrays(e.astype(bool), t)
    mx = max(mx, abs(concordance_index_ipcw(y, y, s, tau=60)[0] - uno_c(t, e, s, tau=60)))
    for at in (24, 36, 60):
        try:
            ref = cumulative_dynamic_auc(y, y, s, [at])[0][0]
        except ValueError:
            continue
        mx = max(mx, abs(ref - td_auc(t, e, s, at)))
print("max abs diff", mx)
assert mx < 1e-10
