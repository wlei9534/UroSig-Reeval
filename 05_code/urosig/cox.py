"""Small, fast Cox PH fitter (Efron ties, Newton-Raphson) for bootstrap loops."""
from __future__ import annotations

import numpy as np


def _prep(T, E, X):
    T = np.asarray(T, float); E = np.asarray(E, int); X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    o = np.argsort(T, kind="mergesort")
    return T[o], E[o], X[o]


def _stats(T, E, X, beta):
    n, p = X.shape
    eta = X @ beta
    eta = eta - eta.max()
    w = np.exp(eta)
    wx = w[:, None] * X
    wxx = wx[:, :, None] * X[:, None, :]
    # risk-set sums: suffix cumsums (rows sorted by time ascending)
    S0c = np.cumsum(w[::-1])[::-1]
    S1c = np.cumsum(wx[::-1], axis=0)[::-1]
    S2c = np.cumsum(wxx[::-1], axis=0)[::-1]
    ev_times, first_idx = np.unique(T, return_index=True)
    ev = E == 1
    times_with_ev = np.unique(T[ev])
    pos = np.searchsorted(ev_times, times_with_ev)
    fi = first_idx[pos]
    S0, S1, S2 = S0c[fi], S1c[fi], S2c[fi]
    grp = np.searchsorted(times_with_ev, T[ev])
    k = len(times_with_ev)
    d = np.bincount(grp, minlength=k).astype(float)
    T0 = np.bincount(grp, weights=w[ev], minlength=k)
    T1 = np.zeros((k, p)); np.add.at(T1, grp, wx[ev])
    T2 = np.zeros((k, p, p)); np.add.at(T2, grp, wxx[ev])
    sx = X[ev].sum(0)
    ll = eta[ev].sum(); grad = sx.copy(); hess = np.zeros((p, p))
    for l in range(int(d.max())):
        m = d > l
        f = l / d[m]
        den = S0[m] - f * T0[m]
        num1 = S1[m] - f[:, None] * T1[m]
        num2 = S2[m] - f[:, None, None] * T2[m]
        ll -= np.log(den).sum()
        grad -= (num1 / den[:, None]).sum(0)
        hess -= (num2 / den[:, None, None]).sum(0) - np.einsum("ki,kj->ij", num1 / den[:, None], num1 / den[:, None])
    # undo the eta shift in the log-likelihood
    return ll, grad, hess


def cox_fit(T, E, X, max_iter: int = 50, tol: float = 1e-9, ridge: float = 1e-8):
    """Return (beta, log partial likelihood). Ties handled with Efron's method."""
    T, E, X = _prep(T, E, X)
    p = X.shape[1]
    beta = np.zeros(p)
    ll_old = None
    for _ in range(max_iter):
        ll, g, h = _stats(T, E, X, beta)
        step = np.linalg.solve(-h + ridge * np.eye(p), g)
        # step halving for stability
        for _h in range(20):
            nb = beta + step
            ll_new = _stats(T, E, X, nb)[0]
            if ll_new >= ll - 1e-12:
                break
            step /= 2
        beta = nb
        if ll_old is not None and abs(ll_new - ll_old) < tol:
            break
        ll_old = ll_new
    ll_final = _loglik_exact(T, E, X, beta)
    return beta, ll_final


def _loglik_exact(T, E, X, beta):
    """Efron log partial likelihood without the max-shift (for LR tests vs lifelines)."""
    eta = X @ beta
    w = np.exp(eta)
    S0c = np.cumsum(w[::-1])[::-1]
    ev = E == 1
    uniq, first_idx = np.unique(T, return_index=True)
    times_with_ev = np.unique(T[ev])
    fi = first_idx[np.searchsorted(uniq, times_with_ev)]
    grp = np.searchsorted(times_with_ev, T[ev])
    k = len(times_with_ev)
    d = np.bincount(grp, minlength=k).astype(float)
    T0 = np.bincount(grp, weights=w[ev], minlength=k)
    S0 = S0c[fi]
    ll = eta[ev].sum()
    for l in range(int(d.max())):
        m = d > l
        ll -= np.log(S0[m] - (l / d[m]) * T0[m]).sum()
    return float(ll)


def cox_lp(T, E, X):
    """Fit and return (linear predictor in input order, loglik)."""
    Xa = np.asarray(X, float)
    if Xa.ndim == 1:
        Xa = Xa[:, None]
    beta, ll = cox_fit(T, E, Xa)
    return Xa @ beta, ll
