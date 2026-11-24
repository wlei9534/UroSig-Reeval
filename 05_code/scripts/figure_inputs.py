"""Aggregate, non-identifying inputs for figures that need patient-level data (run where the cohorts live).

Usage: python figure_inputs.py <project_root>
Writes <ev>/figure_inputs/:
  random_hist.csv   per cohort: histogram (bin width 0.01) of the C-index of all random gene sets (1,000 per pair)
  km_curves.csv     Kaplan-Meier step coordinates for BCR-free survival, pooled over the 4 primary cohorts:
                    ISUP grade group (1, 2, 3, 4-5) and within-cohort tertiles of two signatures chosen by a fixed
                    rule among signatures evaluable in all 4 cohorts: highest pooled C ("top") and median pooled C
                    ("median"). No selection on how the curves look.
  km_meta.csv       per panel: n, events, log-rank P (stratified by cohort), numbers at risk at 0/24/48/72/96 months
  rho_grade.csv     per signature x cohort: Spearman correlation of the score with ISUP grade group
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from scipy import stats

ROOT = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(ROOT / "05_code"))
from urosig.paths import p as P  # noqa: E402
from urosig.sigscore import score_signature  # noqa: E402

EV = P(ROOT, "ev")
OUT = EV / "figure_inputs"
OUT.mkdir(parents=True, exist_ok=True)
MAIN4 = ["GSE54460", "GSE70769", "GSE70768", "MSKCC_GSE21032"]
terms = pd.read_csv(P(ROOT, "terms"), dtype=str).fillna("")
dirs_df = pd.read_csv(P(ROOT, "dirs"), dtype={"sig_id": str})
pairs = pd.read_csv(EV / "run_primary" / "pairs.csv", dtype={"sig_id": str})
pairs = pairs[pairs["status"] == "ok"]
sig = pd.read_csv(EV / "primary" / "signatures.csv", dtype={"sig_id": str})

# ---- random gene-set distributions
z = np.load(EV / "run_primary" / "random_c.npz", allow_pickle=True)
bins = np.round(np.arange(0.25, 0.951, 0.01), 2)
rows = []
for c in MAIN4:
    v = np.concatenate([z[k] for k in z.keys() if k.endswith("|" + c)])
    h, _ = np.histogram(v, bins=bins)
    for b, n in zip(bins[:-1], h):
        rows.append({"cohort": c, "bin_left": b, "count": int(n), "n_total": len(v), "q025": np.quantile(v, .025),
                     "q50": np.quantile(v, .5), "q95": np.quantile(v, .95), "q975": np.quantile(v, .975)})
pd.DataFrame(rows).to_csv(OUT / "random_hist.csv", index=False)

# ---- scores for every pair (z mode, as in the primary analysis)
full = sig[sig["k_cohorts"] == 4].sort_values("pooled_c").reset_index(drop=True)
pick = {"top": full["sig_id"].iloc[-1], "median": full["sig_id"].iloc[len(full) // 2]}
pooled, rho = [], []
for c in MAIN4:
    with gzip.open(ROOT / "04_cohorts/processed" / c / "expr_log2.csv.gz", "rt", encoding="utf-8") as f:
        expr = pd.read_csv(f, index_col=0)
    clin = pd.read_csv(ROOT / "04_cohorts/processed" / c / "clinical.csv", index_col=0)
    clin = clin[clin["time_bcr"].notna() & clin["bcr"].notna()]
    clin = clin.loc[clin.index.intersection(expr.columns)]
    expr = expr[clin.index]
    d0 = clin[["time_bcr", "bcr", "isup", "isup_cat"]].copy()
    d0["cohort"] = c
    for sid in pairs.loc[pairs["cohort"] == c, "sig_id"]:
        tt = terms[terms["sig_id"] == sid]
        d = dict(zip(dirs_df.loc[dirs_df["sig_id"] == sid, "基因（原文写法）"], dirs_df.loc[dirs_df["sig_id"] == sid, "方向"]))
        sc = score_signature(sid, c, tt, expr, "z", d, None)
        if sc.values is None:
            continue
        s = sc.values.reindex(clin.index)
        ok = s.notna() & clin["isup"].notna()
        rho.append({"sig_id": sid, "cohort": c, "rho_grade": stats.spearmanr(s[ok], clin.loc[ok, "isup"]).statistic, "n": int(ok.sum())})
        for name, pid in pick.items():
            if sid == pid:
                d0["tert_" + name] = pd.qcut(s, 3, labels=["Low", "Middle", "High"]).astype(str)
    pooled.append(d0)
pd.DataFrame(rho).to_csv(OUT / "rho_grade.csv", index=False)
D = pd.concat(pooled, ignore_index=True)

curves, meta = [], []
for panel, col, order in (("grade", "isup_cat", ["1", "2", "3", "4-5"]), ("top", "tert_top", ["Low", "Middle", "High"]),
                          ("median", "tert_median", ["Low", "Middle", "High"])):
    d = D[D[col].isin(order)].copy()
    d["g"] = d[col].map({k: i for i, k in enumerate(order)})
    cph = CoxPHFitter().fit(d[["time_bcr", "bcr", "g", "cohort"]], "time_bcr", "bcr", strata=["cohort"], formula="C(g)")
    p = float(cph.log_likelihood_ratio_test().p_value)
    for g in order:
        x = d[d[col] == g]
        km = KaplanMeierFitter().fit(x["time_bcr"], x["bcr"])
        sf = km.survival_function_.reset_index()
        sf.columns = ["time", "surv"]
        sf = sf[sf["time"] <= 120]
        curves.append(sf.assign(panel=panel, group=g))
        meta.append({"panel": panel, "group": g, "n": len(x), "events": int(x["bcr"].sum()), "p_stratified_lr": p,
                     "sig_id": pick.get(panel, ""), **{f"at_risk_{t}": int((x["time_bcr"] >= t).sum()) for t in (0, 24, 48, 72, 96)}})
pd.concat(curves).to_csv(OUT / "km_curves.csv", index=False)
pd.DataFrame(meta).to_csv(OUT / "km_meta.csv", index=False)
print(pick, len(D), pd.DataFrame(meta)[["panel", "group", "n", "events", "p_stratified_lr"]].to_string())
