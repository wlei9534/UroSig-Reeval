"""Pool signature x cohort results per signature and summarise the signature population (SAP v2 section 8).

Usage: python pool_results.py <project_root> --tag primary
Reads 06_results/evaluation/<tag>/pairs.csv + random_c.npz; writes signatures.csv and summary.json there.
Per signature: random-effects meta-analysis of logit(Harrell C) over its evaluable cohorts (REML, Hartung-Knapp;
k = 1 -> that cohort's C). Grade benchmark pooled with the same weights. Random baseline: for each random draw r
the random C of every cohort is pooled with the signature's weights, giving a pooled random distribution and an
empirical p; Benjamini-Hochberg over signatures. Paper-level clustering: 32264916 and 32264916_OS count once in
proportions.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from urosig.paths import p as P, im_globs, _V as PATHS  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("root")
ap.add_argument("--tag", default="primary", help="output tag (and input, unless --src is given)")
ap.add_argument("--src", default="", help="comma list of run tags whose pairs are pooled together")
ap.add_argument("--cohorts", default="", help="keep only these cohorts")
ap.add_argument("--cov", default="", choices=["", "0.8", "1.0"], help="keep only pairs eligible at this coverage")
ap.add_argument("--endpoint", default="", help="keep only signatures of this original-endpoint class (S6: BCR/RFS)")
ap.add_argument("--roles", default="", help="keep only these analysis roles")
args = ap.parse_args()
ROOT = Path(args.root)
sys.path.insert(0, str(ROOT / "05_code"))
from scipy import stats  # noqa: E402

from urosig.meta import inv_logit, random_effects  # noqa: E402

EV = P(ROOT, "ev")
D = EV / args.tag
D.mkdir(parents=True, exist_ok=True)
src = args.src.split(",") if args.src else [args.tag]
p = pd.concat([pd.read_csv(EV / t / "pairs.csv", dtype={"sig_id": str}) for t in src], ignore_index=True)
p = p[p["status"] == "ok"].copy()
p = p.drop_duplicates(["sig_id", "cohort"], keep="last")  # a later --src run replaces the same pair (S10)
rand = {}
for t in src:
    rand.update(dict(np.load(EV / t / "random_c.npz")))
elig = pd.read_csv(P(ROOT, "elig"), dtype={"sig_id": str})
p = p.merge(elig[["sig_id", "cohort", "分析角色", "可评估_1.0"]], on=["sig_id", "cohort"], how="left")
if args.cohorts:
    p = p[p["cohort"].isin(args.cohorts.split(","))]
if args.cov == "1.0":
    p = p[p["可评估_1.0"] == True]  # noqa: E712
if args.roles:
    p = p[p["分析角色"].isin(args.roles.split(","))]
if args.endpoint:
    ep = pd.read_csv(P(ROOT, "endpoint"), dtype={"sig_id": str})
    p = p[p["sig_id"].isin(ep.loc[ep["终点类别"] == args.endpoint, "sig_id"])]
filters = {k: getattr(args, k) for k in ("src", "cohorts", "cov", "endpoint", "roles") if getattr(args, k)}
lg = lambda c: np.log(c / (1 - c))  # noqa: E731


def bh(pv):
    pv = np.asarray(pv, float)
    q = np.full_like(pv, np.nan)
    ok = ~np.isnan(pv)
    x = pv[ok]
    if len(x):
        o = np.argsort(x)
        r = x[o] * len(x) / np.arange(1, len(x) + 1)
        r = np.minimum.accumulate(r[::-1])[::-1]
        q_ok = np.empty_like(x); q_ok[o] = np.minimum(r, 1)
        q[ok] = q_ok
    return q


rows = []
for sid, g in p.groupby("sig_id"):
    se = g["logit_c_se"].to_numpy(float)
    y = g["logit_c"].to_numpy(float)
    k = len(g)
    if k > 1:
        re = random_effects(y, se ** 2)
        est, lo, hi, w = re["est"], re["lo"], re["hi"], 1 / (se ** 2 + re["tau2"])
        i2, tau2 = re["i2"], re["tau2"]
    else:
        est, lo, hi, w, i2, tau2 = y[0], y[0] - 1.96 * se[0], y[0] + 1.96 * se[0], np.array([1.0]), np.nan, np.nan
    w = w / w.sum()
    grade = float(inv_logit((w * lg(g["grade_c"].to_numpy(float))).sum()))
    sig_on_grade = float(inv_logit((w * lg(g["sig_c_on_grade_pts"].to_numpy(float))).sum()))
    r = {"sig_id": sid, "PMID": sid.split("_")[0], "k_cohorts": k, "cohorts": ";".join(g["cohort"]),
         "pooled_c": float(inv_logit(est)), "pooled_c_lo": float(inv_logit(lo)), "pooled_c_hi": float(inv_logit(hi)),
         "i2": i2, "tau2": tau2, "pooled_grade_c": grade, "pooled_sig_c_same_pts": sig_on_grade,
         "beats_grade": sig_on_grade > grade,
         "mean_delta_c_corrected": float(np.average(g["delta_c_corrected"], weights=w)) if g["delta_c_corrected"].notna().all() else np.nan,
         "min_lr_p": float(g["lr_p"].min()) if "lr_p" in g else np.nan,
         # per-signature incremental-value test over its cohorts: independent patients, so the LR chi2 add up (df = k)
         "lr_chi2_sum": float(g["lr_chi2"].sum()) if g["lr_chi2"].notna().all() else np.nan,
         "lr_p_combined": float(stats.chi2.sf(g["lr_chi2"].sum(), k)) if g["lr_chi2"].notna().all() else np.nan,
         "median_uno_c": float(g["uno_c"].median()), "median_auc36": float(g["auc36"].median()),
         "median_hr_sd": float(g["hr_sd"].median()), "ci_excludes_0.5": bool(inv_logit(lo) > 0.5)}
    keys = [f"{sid}|{c}" for c in g["cohort"]]
    if all(kk in rand for kk in keys):
        R = np.vstack([lg(np.clip(rand[kk], 1e-6, 1 - 1e-6)) for kk in keys])  # cohorts x draws
        pooled_rand = (w[:, None] * R).sum(0)
        obs = (w * y).sum()
        r["rand_p_pooled"] = float((1 + (pooled_rand >= obs).sum()) / (1 + R.shape[1]))
        r["beats_random_95"] = bool(obs > np.percentile(pooled_rand, 95))
    rows.append(r)
s = pd.DataFrame(rows)
if "rand_p_pooled" in s:
    s["rand_q_bh"] = bh(s["rand_p_pooled"])
s["lr_q_bh"] = bh(s["min_lr_p"])  # as first coded: smallest p over cohorts (anti-conservative; kept for transparency)
s["lr_q_bh_combined"] = bh(s["lr_p_combined"])
s.to_csv(D / "signatures.csv", index=False, encoding="utf-8-sig")

# paper-level proportions: one row per PMID (32264916 and its OS model count once, using the mean of the indicators)
pl = s.groupby("PMID").agg(pooled_c=("pooled_c", "mean"), beats_grade=("beats_grade", "mean"),
                           beats_random=("beats_random_95", "mean") if "beats_random_95" in s else ("beats_grade", "size"))
summ = {"tag": args.tag, "filters": filters, "n_signatures": int(len(s)), "n_papers": int(len(pl)), "n_pairs": int(len(p)),
        "pooled_c_median": float(s["pooled_c"].median()),
        "pooled_c_iqr": [float(s["pooled_c"].quantile(0.25)), float(s["pooled_c"].quantile(0.75))],
        "share_beats_grade_papers": float(pl["beats_grade"].mean()),
        "share_beats_random_papers": float(pl["beats_random"].mean()) if "beats_random_95" in s else None,
        "n_random_q_lt_0.05": int((s.get("rand_q_bh", pd.Series(dtype=float)) < 0.05).sum()),
        "n_lr_q_lt_0.05": int((s["lr_q_bh"] < 0.05).sum()),
        "n_lr_combined_q_lt_0.05": int((s["lr_q_bh_combined"] < 0.05).sum()),
        "n_signatures_k2plus": int((s["k_cohorts"] >= 2).sum()),
        "n_signatures_with_random": int(s["rand_p_pooled"].notna().sum()) if "rand_p_pooled" in s else 0,
        "n_ci_excludes_0.5": int(s["ci_excludes_0.5"].sum()),
        "median_delta_c_corrected": float(s["mean_delta_c_corrected"].median())}
(D / "summary.json").write_text(json.dumps(summ, ensure_ascii=False, indent=1), encoding="utf-8")
print(json.dumps(summ, ensure_ascii=False, indent=1))
