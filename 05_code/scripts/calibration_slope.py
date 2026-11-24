"""V5 (post hoc, added 2026-11-18 for TRIPOD+AI item on calibration): calibration slope of the published linear predictors.

Usage: python calibration_slope.py <project_root>
For every primary signature x cohort pair of kind 'linear' (published coefficients; not gene pairs), the published
linear predictor is computed on unstandardised log2 expression (the S2a 'raw' scale, the only scale on which the
published coefficients keep their units) and entered as the single covariate of a Cox model for BCR. The fitted
coefficient is the calibration slope: 1 = published weights are correctly scaled in the new cohort, <1 = predictions
too extreme (overfitting), <=0 = no or reversed association. Wald 95% CI from the observed information.
Per signature: inverse-variance (fixed-effect) pooled slope over its cohorts, for description only.
Limits: no paper gives a baseline hazard, so calibration-in-the-large / absolute risk cannot be assessed; platforms
and normalisation differ from the development data, so the slope mixes overfitting with scale differences.
Outputs: <ev>/posthoc/V5_calibration_slope_pairs.csv, V5_calibration_slope_signatures.csv, summary_V5.json
"""
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(ROOT / "05_code"))
from urosig.cox import _prep, _stats, cox_fit  # noqa: E402
from urosig.paths import p as P  # noqa: E402
from urosig.sigscore import score_signature  # noqa: E402

EV = P(ROOT, "ev")
OUT = EV / "posthoc"
OUT.mkdir(parents=True, exist_ok=True)
terms = pd.read_csv(P(ROOT, "terms"), dtype=str).fillna("")
dirs_df = pd.read_csv(P(ROOT, "dirs"), dtype={"sig_id": str})
pairs = pd.read_csv(EV / "run_raw" / "pairs.csv", dtype={"sig_id": str})
pairs = pairs[(pairs["status"] == "ok") & (pairs["kind"] == "linear")]

rows = []
for c, sub in pairs.groupby("cohort"):
    with gzip.open(ROOT / "04_cohorts/processed" / c / "expr_log2.csv.gz", "rt", encoding="utf-8") as f:
        expr = pd.read_csv(f, index_col=0)
    clin = pd.read_csv(ROOT / "04_cohorts/processed" / c / "clinical.csv", index_col=0)
    clin = clin[clin["time_bcr"].notna() & clin["bcr"].notna()]
    clin = clin.loc[clin.index.intersection(expr.columns)]
    expr = expr[clin.index]
    T, E = clin["time_bcr"].to_numpy(float), clin["bcr"].to_numpy(int)
    for sid in sub["sig_id"]:
        tt = terms[terms["sig_id"] == sid]
        d = dict(zip(dirs_df.loc[dirs_df["sig_id"] == sid, "基因（原文写法）"], dirs_df.loc[dirs_df["sig_id"] == sid, "方向"]))
        sc = score_signature(sid, c, tt, expr, "raw", d, None)
        if sc.values is None:
            continue
        s = sc.values.reindex(clin.index).to_numpy(float)
        ok = ~np.isnan(s)
        t, e, s = T[ok], E[ok], s[ok]
        sd = s.std(ddof=1)
        if not np.isfinite(sd) or sd == 0:
            continue
        z = (s - s.mean()) / sd  # fit on the standardised score for numerical stability, then rescale
        b, _ = cox_fit(t, e, z[:, None])
        _, _, h = _stats(*_prep(t, e, z[:, None]), b)
        se = float(np.sqrt(-1 / h[0, 0]))
        rows.append({"sig_id": sid, "cohort": c, "n": len(t), "events": int(e.sum()), "lp_sd": sd,
                     "slope": b[0] / sd, "slope_se": se / sd, "slope_lo": (b[0] - 1.96 * se) / sd, "slope_hi": (b[0] + 1.96 * se) / sd})
r = pd.DataFrame(rows)
r.to_csv(OUT / "V5_calibration_slope_pairs.csv", index=False, encoding="utf-8-sig")


def pool(g):
    w = 1 / g["slope_se"] ** 2
    m = float((w * g["slope"]).sum() / w.sum())
    se = float(np.sqrt(1 / w.sum()))
    return pd.Series({"k": len(g), "slope": m, "slope_lo": m - 1.96 * se, "slope_hi": m + 1.96 * se, "median_lp_sd": g["lp_sd"].median()})


s = r.groupby("sig_id").apply(pool).reset_index()
s.to_csv(OUT / "V5_calibration_slope_signatures.csv", index=False, encoding="utf-8-sig")


def block(d):
    return {"n": int(len(d)), "median_slope": float(d["slope"].median()), "iqr": [float(d["slope"].quantile(.25)), float(d["slope"].quantile(.75))],
            "share_slope_below_1": float((d["slope"] < 1).mean()), "share_ci_below_1": float((d["slope_hi"] < 1).mean()),
            "share_ci_includes_1": float(((d["slope_lo"] <= 1) & (d["slope_hi"] >= 1)).mean()),
            "share_ci_above_1": float((d["slope_lo"] > 1).mean()), "share_slope_le_0": float((d["slope"] <= 0).mean()),
            "share_ci_includes_0": float(((d["slope_lo"] <= 0) & (d["slope_hi"] >= 0)).mean())}


summ = {"pairs": block(r), "signatures": block(s), "by_cohort": {c: block(g) for c, g in r.groupby("cohort")},
        "median_lp_sd_pairs": float(r["lp_sd"].median())}
(OUT / "summary_V5.json").write_text(json.dumps(summ, ensure_ascii=False, indent=1), encoding="utf-8")
print(json.dumps({k: summ[k] for k in ("pairs", "signatures")}, ensure_ascii=False, indent=1))
