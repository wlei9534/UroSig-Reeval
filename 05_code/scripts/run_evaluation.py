"""External evaluation of locked signatures (SAP v2 sections 5-7). Run only after the SAP is final.

Usage: python run_evaluation.py <project_root> --mode z --cov 0.8 --cohorts GSE54460,GSE70769,GSE70768,MSKCC_GSE21032
                                [--roles main] [--nboot 1000] [--nrand 1000] [--seed 20261101] [--tag primary]
Per eligible signature x cohort (eligibility table, threshold --cov), with the score computed in --mode:
  Harrell C, Uno C (tau 60 mo), td-AUC 36/60 mo, each with a 1,000x patient bootstrap percentile CI and the
  bootstrap SD of logit(C) for pooling; HR per SD (Cox);
  benchmark: grade C (ISUP; Gleason-sum groups where ISUP is incomplete) on the same patients;
  incremental value: extended clinical Cox vs + score, apparent and optimism-corrected delta C, LR test;
  random baseline: --nrand expression-matched random gene sets carrying the signature's coefficients.
Writes 06_results/evaluation/<tag>/pairs.csv (+ random_c.npz with the random C per pair for pooling).
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
import time as _time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from urosig.paths import p as P, im_globs, _V as PATHS  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("root")
ap.add_argument("--mode", default="z", choices=["z", "raw", "tcga"])
ap.add_argument("--cov", default="0.8", choices=["0.8", "1.0"])
ap.add_argument("--cohorts", default="GSE54460,GSE70769,GSE70768,MSKCC_GSE21032")
ap.add_argument("--roles", default="main")
ap.add_argument("--nboot", type=int, default=1000)
ap.add_argument("--nrand", type=int, default=1000)
ap.add_argument("--seed", type=int, default=20261101)
ap.add_argument("--tag", default="primary")
ap.add_argument("--limit", type=int, default=0, help="debug: first N pairs only")
ap.add_argument("--full-formula", action="store_true", help="S9: add the published clinical terms")
ap.add_argument("--endpoint", default="bcr", choices=["bcr", "met"], help="E3: metastasis (GSE116918 only)")
ap.add_argument("--drop-repeats", action="store_true", help="S10: GSE54460 without the 6 patients with repeat libraries")
ap.add_argument("--only", default="", help="comma list of sig_ids to run (S9)")
ap.add_argument("--shard", default="", help="'i/n': run every n-th pair starting at i (seeds are per pair, so results do not depend on n)")
args = ap.parse_args()

ROOT = Path(args.root)
sys.path.insert(0, str(ROOT / "05_code"))
from lifelines.utils import concordance_index  # noqa: E402
from scipy import stats  # noqa: E402

from urosig.cox import _prep, _stats, cox_fit  # noqa: E402
from urosig.nullsets import matched_random_sets  # noqa: E402
from urosig.sigscore import score_signature  # noqa: E402
from urosig.timeauc import td_auc, uno_c  # noqa: E402

OUT = P(ROOT, "ev") / args.tag
OUT.mkdir(parents=True, exist_ok=True)
TAU, T_AUC = 60.0, (36.0, 60.0)
terms = pd.read_csv(P(ROOT, "terms"), dtype=str).fillna("")
elig = pd.read_csv(P(ROOT, "elig"), dtype={"sig_id": str})
dirs_df = pd.read_csv(P(ROOT, "dirs"), dtype={"sig_id": str})
leak = pd.read_csv(P(ROOT, "leak"), dtype={"sig_id": str}).set_index("sig_id")
roles = args.roles.split(",")
cohorts = args.cohorts.split(",")
pairs = elig[elig[f"可评估_{args.cov}"] & elig["cohort"].isin(cohorts) & elig["分析角色"].isin(roles)]
pairs = pairs[["sig_id", "cohort"]].reset_index(drop=True)
if args.only:
    pairs = pairs[pairs["sig_id"].isin(args.only.split(","))].reset_index(drop=True)
if args.limit:
    pairs = pairs.head(args.limit)
SHARD = ""
if args.shard:
    si, sn = map(int, args.shard.split("/"))
    pairs = pairs.iloc[si::sn]
    SHARD = f"_shard{si:02d}of{sn:02d}"

# published clinical terms (S9 only): sig_id -> {covariate: coef}
CLIN_TERMS = {"34002865": {"gleason_sum": 0.5, "stage": 0.32}, "39832058": {"gleason_sum": 0.63}}


def load_cohort(c):
    with gzip.open(ROOT / "04_cohorts/processed" / c / "expr_log2.csv.gz", "rt", encoding="utf-8") as f:
        expr = pd.read_csv(f, index_col=0)
    clin = pd.read_csv(ROOT / "04_cohorts/processed" / c / "clinical.csv", index_col=0)
    if args.endpoint == "met":  # E3: metastasis-free interval in place of BCR (same code path)
        clin = clin.drop(columns=["time_bcr", "bcr"]).rename(columns={"time_met": "time_bcr", "met": "bcr"})
    clin = clin[clin["time_bcr"].notna() & clin["bcr"].notna()]
    clin = clin.loc[clin.index.intersection(expr.columns)]
    if c == "GSE54460" and args.drop_repeats:  # S10: drop the 6 patients whose repeat library had other follow-up
        clin = clin[clin["flag_inconsistent_repeat"] != 1]
    return expr[clin.index], clin


def harrell(t, e, s):
    return float(concordance_index(t, -np.asarray(s, float), e))


def grade_var(clin):
    if clin["isup"].notna().mean() >= 0.9:
        return clin["isup"].astype(float), "ISUP grade group"
    g = clin["gleason_sum"].astype(float)
    return pd.Series(np.select([g <= 6, g == 7], [1.0, 2.0], 3.0), index=clin.index).where(g.notna()), "Gleason sum 6/7/8-10"


def ext_clin(clin):
    cols = {}
    gv, _ = grade_var(clin)
    cols["grade"] = gv
    if clin["pt_cat"].notna().mean() >= 0.5:
        cols["pT3plus"] = clin["pt_cat"].map(lambda v: np.nan if pd.isna(v) else float(str(v) != "T2"))
    if clin["log_psa"].notna().mean() >= 0.5:
        cols["log_psa"] = clin["log_psa"].astype(float)
    if clin["margin"].notna().mean() >= 0.5:
        cols["margin"] = clin["margin"].astype(float)
    return pd.DataFrame(cols, index=clin.index)


def cox_c(t, e, X):
    beta, ll = cox_fit(t, e, X)
    return X @ beta, ll


def incremental(t, e, X0, s, rng, nboot):
    # the score enters the Cox model standardised: C and the LR test are unchanged by an affine rescaling, and
    # very large raw scores (e.g. coefficients read from figures, |b| ~ 10) no longer overflow exp() in the fit
    s = (s - s.mean()) / s.std(ddof=1)
    X1 = np.column_stack([X0, s])
    lp0, ll0 = cox_c(t, e, X0)
    lp1, ll1 = cox_c(t, e, X1)
    c0, c1 = harrell(t, e, lp0), harrell(t, e, lp1)
    lr = 2 * (ll1 - ll0)
    opt = []
    n = len(t)
    for _ in range(nboot):
        i = rng.integers(0, n, n)
        if e[i].sum() < 3:
            continue
        try:
            b0, _ = cox_fit(t[i], e[i], X0[i]); b1, _ = cox_fit(t[i], e[i], X1[i])
            d_boot = harrell(t[i], e[i], X1[i] @ b1) - harrell(t[i], e[i], X0[i] @ b0)
            d_orig = harrell(t, e, X1 @ b1) - harrell(t, e, X0 @ b0)
        except Exception:  # noqa: BLE001  (non-converged bootstrap fit: skip this replicate)
            continue
        opt.append(d_boot - d_orig)
    d = c1 - c0
    return {"c_clin": c0, "c_clin_sig": c1, "delta_c": d, "delta_c_corrected": d - (np.mean(opt) if opt else np.nan),
            "lr_chi2": lr, "lr_p": float(stats.chi2.sf(lr, 1)), "n_inc": n, "ev_inc": int(e.sum()), "clin_vars": ""}


def boot_metrics(t, e, s, rng, nboot):
    n = len(t)
    out = {"harrell": [], "uno": [], "auc36": [], "auc60": []}
    for _ in range(nboot):
        i = rng.integers(0, n, n)
        if e[i].sum() < 2:
            continue
        out["harrell"].append(harrell(t[i], e[i], s[i]))
        out["uno"].append(uno_c(t[i], e[i], s[i], tau=TAU))
        out["auc36"].append(td_auc(t[i], e[i], s[i], 36.0))
        out["auc60"].append(td_auc(t[i], e[i], s[i], 60.0))
    return {k: np.array(v, float) for k, v in out.items()}


def auc_ok(t, e, at):
    return ((t > at).sum() >= 10) and (((t <= at) & (e == 1)).sum() >= 10)


tcga_sd = None
if args.mode == "tcga":
    with gzip.open(ROOT / "04_cohorts/processed/TCGA-PRAD/expr_log2.csv.gz", "rt", encoding="utf-8") as f:
        tcga_sd = pd.read_csv(f, index_col=0).std(axis=1, ddof=1)

rows, rand_store = [], {}
t0 = _time.time()
for c in cohorts:
    sub = pairs[pairs["cohort"] == c]
    if sub.empty:
        continue
    expr, clin = load_cohort(c)
    T = clin["time_bcr"].to_numpy(float); E = clin["bcr"].to_numpy(int)
    mean_expr = expr.mean(axis=1)
    gvar, gname = grade_var(clin)
    xc = ext_clin(clin)
    for _, p in sub.iterrows():
        sid = p["sig_id"]
        tt = terms[terms["sig_id"] == sid]
        d = dict(zip(dirs_df.loc[dirs_df["sig_id"] == sid, "基因（原文写法）"], dirs_df.loc[dirs_df["sig_id"] == sid, "方向"]))
        if args.mode == "tcga" and leak.loc[sid, "TCGA"] != "用过":
            continue  # TCGA-SD rescaling applies to signatures whose paper used TCGA
        sc = score_signature(sid, c, tt, expr, args.mode, d, tcga_sd)
        if sc.values is None:
            rows.append({"sig_id": sid, "cohort": c, "status": "not scorable", "note": sc.note})
            continue
        s = sc.values.reindex(clin.index).to_numpy(float)
        if args.full_formula and sid in CLIN_TERMS:
            for k, b in CLIN_TERMS[sid].items():
                v = clin["gleason_sum"].astype(float) if k == "gleason_sum" else clin["pt_cat"].map(lambda x: np.nan if pd.isna(x) else float(str(x) != "T2"))
                s = s + b * v.to_numpy(float)
        ok = ~np.isnan(s)
        t, e, s = T[ok], E[ok], s[ok]
        seed = zlib.crc32(f"{args.seed}|{sid}|{c}".encode())  # stable across runs and machines
        rng = np.random.default_rng(seed)
        r = {"sig_id": sid, "cohort": c, "status": "ok", "mode": args.mode, "kind": sc.kind, "n": len(t), "events": int(e.sum()),
             "n_terms": sc.n_terms, "n_used": sc.n_used, "gene_cov": sc.gene_cov, "weight_cov": sc.weight_cov, "note": sc.note,
             "harrell_c": harrell(t, e, s), "uno_c": uno_c(t, e, s, tau=TAU),
             "auc36": td_auc(t, e, s, 36.0) if auc_ok(t, e, 36.0) else np.nan,
             "auc60": td_auc(t, e, s, 60.0) if auc_ok(t, e, 60.0) else np.nan, "seed": seed}
        bm = boot_metrics(t, e, s, rng, args.nboot)
        for k, key in (("harrell", "harrell_c"), ("uno", "uno_c"), ("auc36", "auc36"), ("auc60", "auc60")):
            v = bm[k][~np.isnan(bm[k])]
            r[key + "_lo"], r[key + "_hi"] = (np.percentile(v, [2.5, 97.5]) if len(v) and not np.isnan(r[key]) else (np.nan, np.nan))
        hb = bm["harrell"]
        hb = hb[(hb > 0) & (hb < 1)]
        r["logit_c"] = float(np.log(r["harrell_c"] / (1 - r["harrell_c"])))
        r["logit_c_se"] = float(np.std(np.log(hb / (1 - hb)), ddof=1))
        z = (s - s.mean()) / s.std(ddof=1)
        beta, _ = cox_fit(t, e, z[:, None])
        # Wald SE from observed information
        _, _, h = _stats(*_prep(t, e, z[:, None]), beta)
        se = float(np.sqrt(-1 / h[0, 0]))
        r.update(hr_sd=float(np.exp(beta[0])), hr_sd_lo=float(np.exp(beta[0] - 1.96 * se)), hr_sd_hi=float(np.exp(beta[0] + 1.96 * se)),
                 hr_sd_p=float(2 * stats.norm.sf(abs(beta[0] / se))))
        # grade benchmark on the same patients
        g = gvar.to_numpy(float)[ok]
        gm = ~np.isnan(g)
        r["grade_var"] = gname
        r["grade_c"] = harrell(t[gm], e[gm], g[gm])
        r["sig_c_on_grade_pts"] = harrell(t[gm], e[gm], s[gm])
        # incremental value over the extended clinical model (complete cases)
        X0 = xc.to_numpy(float)[ok]
        cc = ~np.isnan(X0).any(axis=1)
        if e[cc].sum() >= 5:
            inc = incremental(t[cc], e[cc], X0[cc], s[cc], rng, args.nboot)
            inc["clin_vars"] = "+".join(xc.columns)
            r.update(inc)
        # random baseline: expression-matched gene sets with this signature's coefficients (linear kinds)
        if sc.kind in ("linear", "equal") and args.nrand:
            used = tt[tt["HGNC 现行符号"].isin(expr.index)]
            if sc.kind == "equal":
                w = used["系数"].astype(float) if used["分析角色"].iloc[0] == "equal_signed" else used["基因（原文写法）"].map(lambda g: d.get(g, 0)).astype(float)
                used = used[w.to_numpy() != 0]
                coefs = w[w != 0].to_numpy()
            else:
                coefs = used["系数"].astype(float).to_numpy()
            genes = used["HGNC 现行符号"].tolist()
            sets = matched_random_sets(genes, mean_expr, args.nrand, seed=seed)
            rc = np.empty(len(sets))
            for k, gs in enumerate(sets):
                x = expr.loc[gs].astype(float)
                if args.mode in ("z", "tcga") or sc.kind == "equal":
                    x = x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1, ddof=1).replace(0, np.nan), axis=0).fillna(0)
                    if args.mode == "tcga":
                        x = x.mul(tcga_sd.reindex(gs).fillna(1.0).to_numpy(), axis=0)
                else:  # raw: missing cells -> gene's cohort mean (as in sigscore)
                    x = x.T.fillna(x.mean(axis=1)).T
                rs = (x.T.to_numpy() @ coefs)
                rc[k] = harrell(T, E, rs)
            rand_store[f"{sid}|{c}"] = rc
            r["rand_pct"] = float((rc < r["harrell_c"]).mean() * 100)
            r["rand_p"] = float((1 + (rc >= r["harrell_c"]).sum()) / (1 + len(rc)))
            r["beats_random_95"] = bool(r["harrell_c"] > np.percentile(rc, 95))
        rows.append(r)
        print(f"{len(rows)}/{len(pairs)} {sid} {c} C={r['harrell_c']:.3f} ({_time.time() - t0:.0f}s)", flush=True)

res = pd.DataFrame(rows)
res.to_csv(OUT / f"pairs{SHARD}.csv", index=False, encoding="utf-8-sig")
np.savez_compressed(OUT / f"random_c{SHARD}.npz", **rand_store)
import lifelines, scipy  # noqa: E402,E401
meta = {"args": vars(args), "n_pairs": len(res), "seconds": round(_time.time() - t0),
        "versions": {m.__name__: getattr(m, "__version__", "") for m in (np, pd, scipy, lifelines)},
        "python": sys.version}
(OUT / f"run_info{SHARD}.json").write_text(json.dumps(meta, ensure_ascii=False, indent=1), encoding="utf-8")
print("saved", OUT)
