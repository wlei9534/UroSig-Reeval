"""Post-hoc analyses V1-V4 (specified in decision_log 2026-11-12 before running; V1 confirmatory of the implementation,
V2-V4 exploratory). Usage: python posthoc_v1_v4.py <project_root>

V1  implementation check: score every main signature in the cohorts its paper used (TCGA-PRAD and any GEO/MSKCC cohort
    excluded by the leakage rule), coverage >= 0.8; Harrell C and 36-month td-AUC for BCR; compare with the value the
    paper reported for the same cohort (C-index vs Harrell C; 3-year AUC vs td-AUC 36), whole-cohort roles only
    (training / combined / external; random-split internal subsets excluded). Three-step contrast per signature:
    reported training value -> our recomputation in that training cohort -> independent external pooled C.
V2  precision: minimum detectable signature-minus-grade C difference (80% power, two-sided 0.05: 2.8 x SE) from the
    E1 paired-bootstrap SEs, per cohort and for pooled estimates by number of cohorts.
V3  correlates of external performance: random-effects meta-regression of logit pooled C (weights 1/(v_i + tau2),
    v_i = fixed-effect variance over the signature's cohorts, tau2 by method of moments) on: TCGA-only training,
    external validation reported in the paper, log2 gene count, modelling method, median rho with meta-PCNA.
V4  redundancy: pairwise Spearman correlation of all main-signature scores in TCGA-PRAD (reference, no outcome) and the
    four primary cohorts; average-linkage clustering on 1 - rho; clusters at rho >= 0.5; cluster correlation with meta-PCNA.
Outputs: 06_results/evaluation/posthoc/
"""
import glob
import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from urosig.paths import p as P, im_globs, _V as PATHS  # noqa: E402
from lifelines.utils import concordance_index
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

ROOT = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(ROOT / "05_code"))
from urosig.sigscore import score_signature  # noqa: E402
from urosig.timeauc import td_auc  # noqa: E402

OUT = P(ROOT, "ev") / "posthoc"
OUT.mkdir(parents=True, exist_ok=True)
EV = P(ROOT, "ev")
COHORTS = ["TCGA-PRAD", "GSE54460", "GSE70769", "GSE70768", "MSKCC_GSE21032", "GSE116918"]
MAIN4 = ["GSE54460", "GSE70769", "GSE70768", "MSKCC_GSE21032"]
NAME_RX = {"TCGA-PRAD": r"TCGA", "GSE54460": r"GSE54460", "GSE70769": r"GSE70769|Stockholm",
           "GSE70768": r"GSE70768|Cambridge", "MSKCC_GSE21032": r"GSE21032|GSE21034|MSKCC|Taylor|MSK",
           "GSE116918": r"GSE116918|Belfast"}

terms = pd.read_csv(P(ROOT, "terms"), dtype=str).fillna("")
sig = pd.read_csv(P(ROOT, "sig"), dtype=str).fillna("")
main_ids = sig.loc[sig["分析角色"] == "main", "sig_id"].tolist()
dirs = pd.read_csv(P(ROOT, "dirs"), dtype={"sig_id": str})
elig = pd.read_csv(P(ROOT, "elig"), dtype={"sig_id": str})
prim = pd.read_csv(EV / "primary/signatures.csv", dtype={"sig_id": str}).set_index("sig_id")


def read_cohort(c):
    with gzip.open(ROOT / "04_cohorts/processed" / c / "expr_log2.csv.gz", "rt", encoding="utf-8") as f:
        expr = pd.read_csv(f, index_col=0)
    clin = pd.read_csv(ROOT / "04_cohorts/processed" / c / "clinical.csv", index_col=0)
    clin = clin.loc[clin.index.intersection(expr.columns)]
    return expr[clin.index], clin


def score(sid, c, expr):
    t = terms[terms["sig_id"] == sid]
    d = dict(zip(dirs.loc[dirs["sig_id"] == sid, "基因（原文写法）"], dirs.loc[dirs["sig_id"] == sid, "方向"]))
    return score_signature(sid, c, t, expr, "z", d, None)


COHORTS = [c for c in COHORTS if (ROOT / "04_cohorts/processed" / c / "expr_log2.csv.gz").exists()]
data = {c: read_cohort(c) for c in COHORTS}
cov_ok = elig.assign(wc=pd.to_numeric(elig["权重覆盖率"], errors="coerce"))
cov_ok = cov_ok[(cov_ok["基因覆盖率"] >= 0.8) & (cov_ok["wc"].isna() | (cov_ok["wc"] >= 0.8))]
cov_ok_set = set(zip(cov_ok["sig_id"], cov_ok["cohort"]))
used = set(zip(elig.loc[elig["原文用过"] == True, "sig_id"], elig.loc[elig["原文用过"] == True, "cohort"]))  # noqa: E712

# ---------------- scores for all main signatures in all cohorts (used by V1 and V4)
scores = {c: {} for c in COHORTS}
for c in COHORTS:
    expr, clin = data[c]
    for sid in main_ids:
        if (sid, c) not in cov_ok_set:
            continue
        sc = score(sid, c, expr)
        if sc.values is not None:
            scores[c][sid] = sc.values.reindex(clin.index)

# ---------------- V1
rows = []
for c in COHORTS:
    expr, clin = data[c]
    ok_bcr = clin["time_bcr"].notna() & clin["bcr"].notna()
    T, E = clin.loc[ok_bcr, "time_bcr"].to_numpy(float), clin.loc[ok_bcr, "bcr"].to_numpy(int)
    for sid, s in scores[c].items():
        if (sid, c) not in used:
            continue
        x = s[ok_bcr].to_numpy(float)
        m = ~np.isnan(x)
        rows.append({"sig_id": sid, "cohort": c, "n": int(m.sum()), "events": int(E[m].sum()),
                     "harrell_c": float(concordance_index(T[m], -x[m], E[m])),
                     "td_auc36": td_auc(T[m], E[m], x[m], 36.0)})
v1 = pd.DataFrame(rows)

rep = []
for f in sorted([x for g in im_globs(ROOT) for x in glob.glob(g)]):
    for line in open(f, encoding="utf-8"):
        r = json.loads(line)
        for mt in r.get("metrics", []):
            if mt.get("role") not in ("training", "combined", "external") or mt.get("metric") not in ("C-index", "AUC"):
                continue
            coh = next((k for k, rx in NAME_RX.items() if re.search(rx, str(mt.get("cohort", "")), re.I)), None)
            try:
                val = float(re.search(r"\d+(\.\d+)?", str(mt.get("value"))).group(0))
                val = val / 100 if val > 1 else val
                ty = float(mt.get("time_years")) if mt.get("time_years") not in (None, "", "NA") else np.nan
            except (AttributeError, ValueError, TypeError):
                continue
            if coh and (mt["metric"] == "C-index" or ty == 3.0):
                rep.append({"sig_id": r["sig_id"], "cohort": coh, "role": mt["role"], "metric": mt["metric"],
                            "reported": val, "reported_endpoint": mt.get("endpoint", ""), "from_figure": mt.get("from_figure")})
rep = pd.DataFrame(rep)
rep = rep.groupby(["sig_id", "cohort", "metric"], as_index=False).agg(
    reported=("reported", "median"), role=("role", "first"), reported_endpoint=("reported_endpoint", "first"))
cmp_ = rep.merge(v1, on=["sig_id", "cohort"], how="inner")
cmp_["ours"] = np.where(cmp_["metric"] == "C-index", cmp_["harrell_c"], cmp_["td_auc36"])
cmp_["diff_reported_minus_ours"] = cmp_["reported"] - cmp_["ours"]
cmp_["endpoint_bcr_like"] = cmp_["reported_endpoint"].str.contains(r"BCR|RFS|bRFS|BFS|biochem|recur|relapse|DFS|复发", case=False, na=False)
v1.to_csv(OUT / "V1_used_cohort_recomputation.csv", index=False, encoding="utf-8-sig")
cmp_.to_csv(OUT / "V1_reported_vs_recomputed.csv", index=False, encoding="utf-8-sig")

# three-step contrast: reported training value (E4 rule) -> recomputation in a training cohort -> external pooled C
e4 = pd.read_csv(EV / "exploratory/internal_reported.csv", dtype={"sig_id": str})
tr = cmp_[cmp_["role"] == "training"].copy()
tr["pri"] = tr["metric"].map({"C-index": 0, "AUC": 1})
tr = tr.sort_values(["sig_id", "pri"]).groupby("sig_id").first().reset_index()
three = e4[["sig_id", "reported_training", "training_metric"]].merge(
    tr[["sig_id", "cohort", "metric", "reported", "ours", "endpoint_bcr_like"]], on="sig_id", how="inner")
three = three.merge(prim[["pooled_c"]], left_on="sig_id", right_index=True, how="inner")
three.to_csv(OUT / "V1_three_step.csv", index=False, encoding="utf-8-sig")

# leakage contrast within signature: C in used vs unused GEO/MSKCC cohorts
pp = pd.read_csv(EV / "run_primary/pairs.csv", dtype={"sig_id": str})[["sig_id", "cohort", "harrell_c"]]
lk = v1[v1["cohort"].isin(MAIN4 + ["GSE116918"])].groupby("sig_id")["harrell_c"].median().rename("c_used_geo")
un = pp.groupby("sig_id")["harrell_c"].median().rename("c_unused")
lkc = pd.concat([lk, un], axis=1).dropna()


def desc(x):
    x = pd.Series(x).dropna()
    return {"n": int(len(x)), "median": float(x.median()), "iqr": [float(x.quantile(.25)), float(x.quantile(.75))]} if len(x) else {"n": 0}


s1 = {"n_used_pairs_scored": int(len(v1)), "by_cohort": v1.groupby("cohort").size().to_dict(),
      "recomputed_harrell_c_in_used_cohorts": desc(v1["harrell_c"]),
      "matched_comparisons": int(len(cmp_)),
      "diff_reported_minus_ours_all": desc(cmp_["diff_reported_minus_ours"]),
      "diff_C_index": desc(cmp_.loc[cmp_.metric == "C-index", "diff_reported_minus_ours"]),
      "diff_AUC3y": desc(cmp_.loc[cmp_.metric == "AUC", "diff_reported_minus_ours"]),
      "diff_bcr_like_endpoint": desc(cmp_.loc[cmp_.endpoint_bcr_like, "diff_reported_minus_ours"]),
      "diff_other_endpoint": desc(cmp_.loc[~cmp_.endpoint_bcr_like, "diff_reported_minus_ours"]),
      "spearman_reported_vs_ours": float(stats.spearmanr(cmp_["reported"], cmp_["ours"]).statistic) if len(cmp_) > 3 else None,
      "share_within_0.05": float((cmp_["diff_reported_minus_ours"].abs() <= 0.05).mean()) if len(cmp_) else None,
      "three_step": {"n": int(len(three)), "reported_training": desc(three["reported"]), "ours_in_training_cohort": desc(three["ours"]),
                     "external_pooled_c": desc(three["pooled_c"])},
      "leakage_within_signature": {"n": int(len(lkc)), "c_used_cohorts": desc(lkc["c_used_geo"]), "c_unused_cohorts": desc(lkc["c_unused"]),
                                   "median_paired_diff": float((lkc["c_used_geo"] - lkc["c_unused"]).median()) if len(lkc) else None,
                                   "wilcoxon_p": float(stats.wilcoxon(lkc["c_used_geo"] - lkc["c_unused"]).pvalue) if len(lkc) > 5 else None}}

# ---------------- V2
pe = pd.read_csv(EV / "exploratory/pairs_E1E2.csv", dtype={"sig_id": str})
se_e = pd.read_csv(EV / "exploratory/signatures_E1E2.csv", dtype={"sig_id": str})
v2c = pe.groupby("cohort").agg(n=("n", "max"), events=("events", "max"), median_se=("d_grade_se", "median"))
v2c["MDD_80pct"] = 2.8 * v2c["median_se"]
se_e["half"] = (se_e["d_grade_hi"] - se_e["d_grade_lo"]) / 2
v2k = se_e.groupby("k")["half"].median().rename("median_CI_halfwidth").to_frame()
v2c.to_csv(OUT / "V2_mdd_by_cohort.csv", encoding="utf-8-sig")
v2k.to_csv(OUT / "V2_pooled_ci_halfwidth_by_k.csv", encoding="utf-8-sig")
s2 = {"by_cohort": v2c.round(4).reset_index().to_dict("records"), "pooled_ci_halfwidth_by_k": v2k.round(4).reset_index().to_dict("records")}

# ---------------- V3
rev = {}
for f in sorted(glob.glob(str(ROOT / "03_extraction/review/out_[AM]*.jsonl"))) + [str(ROOT / "03_extraction/review/out_FIG.jsonl")]:
    for line in open(f, encoding="utf-8"):
        r = json.loads(line)
        rev[r["pmid"]] = r
# v2 signatures: training cohort and coefficient type from the reconciled dual extraction (keyed by sig_id)
_rv2 = ROOT / "03_extraction/v2/reconciled_v2.csv"
if _rv2.exists() and "检索记录 PMID" in sig.columns:
    _r2 = pd.read_csv(_rv2, dtype=str).fillna("").set_index("PMID")
    for _sid, _pm in zip(sig["sig_id"], sig["检索记录 PMID"]):
        if _pm in _r2.index and _sid not in rev:
            rev[_sid] = {"training": _r2.at[_pm, "training"], "coef_type": _r2.at[_pm, "coef_type"]}
OTHER_COH = r"GSE\d+|MSKCC|Taylor|DKFZ|CPC|CancerMap|CIT|Cambridge|Stockholm|Belfast|E-MTAB|ICGC|PRJ|in-house|医院|本院|cohort of \d+"
feat = []
for sid in prim.index:
    p = sid.split("_")[0]
    r = rev.get(p, {})
    trn = r.get("training", "")
    tcga_only = bool(re.search("TCGA", trn)) and not re.search(OTHER_COH, trn, re.I)
    ct = r.get("coef_type", "")
    method = "LASSO 类" if re.search(r"LASSO|lasso|Lasso|弹性|elastic|Ridge|ridge|CoxBoost|随机生存|RSF|机器学习", ct) else \
        ("多因素 Cox" if re.search(r"多因素|multivariate|multivariable|逐步|stepwise", ct, re.I) else "单因素或其他")
    ng = int(sig.loc[sig.sig_id == sid, "基因数"].astype(float).iloc[0]) if len(sig.loc[sig.sig_id == sid]) else np.nan
    feat.append({"sig_id": sid, "tcga_only_training": tcga_only, "method": method, "log2_genes": np.log2(ng) if ng else np.nan})
feat = pd.DataFrame(feat)
e4x = e4[["sig_id", "reported_external"]].assign(ext_reported=lambda d: d["reported_external"].notna())
rho = se_e[["sig_id", "median_rho_metapcna"]]
pairs_full = pd.read_csv(EV / "run_primary/pairs.csv", dtype={"sig_id": str})
v_fe = pairs_full.groupby("sig_id")["logit_c_se"].apply(lambda s: 1 / np.sum(1 / s.to_numpy(float) ** 2)).rename("v")
X = feat.merge(e4x[["sig_id", "ext_reported"]], on="sig_id").merge(rho, on="sig_id").merge(v_fe, left_on="sig_id", right_index=True)
X["y"] = np.log(prim.loc[X.sig_id, "pooled_c"].to_numpy() / (1 - prim.loc[X.sig_id, "pooled_c"].to_numpy()))
X.to_csv(OUT / "V3_features.csv", index=False, encoding="utf-8-sig")


def metareg(df, cols):
    d = df.dropna(subset=cols + ["y", "v"])
    Xm = pd.get_dummies(d[cols], drop_first=True, dtype=float)
    Xm.insert(0, "const", 1.0)
    A, y, v = Xm.to_numpy(float), d["y"].to_numpy(float), d["v"].to_numpy(float)
    w = 1 / v
    W = np.diag(w)
    XtWX_inv = np.linalg.pinv(A.T @ W @ A)
    b = XtWX_inv @ A.T @ W @ y
    Q = float(((y - A @ b) ** 2 * w).sum())
    kdf = len(y) - A.shape[1]
    P = W - W @ A @ XtWX_inv @ A.T @ W
    tau2 = max(0.0, (Q - kdf) / np.trace(P))
    w2 = 1 / (v + tau2)
    W2 = np.diag(w2)
    cov = np.linalg.pinv(A.T @ W2 @ A)
    b2 = cov @ A.T @ W2 @ y
    se = np.sqrt(np.diag(cov))
    z = b2 / se
    res = pd.DataFrame({"term": Xm.columns, "coef_logitC": b2, "se": se, "lo": b2 - 1.96 * se, "hi": b2 + 1.96 * se,
                        "p": 2 * stats.norm.sf(np.abs(z))})
    res["n"] = len(y)
    res["tau2"] = tau2
    return res


uni = pd.concat([metareg(X, [c]).assign(model=f"单变量:{c}") for c in
                 ["tcga_only_training", "ext_reported", "log2_genes", "method", "median_rho_metapcna"]], ignore_index=True)
multi = metareg(X, ["tcga_only_training", "ext_reported", "log2_genes", "method", "median_rho_metapcna"]).assign(model="多变量")
v3 = pd.concat([uni, multi], ignore_index=True)
v3 = v3[v3["term"] != "const"]
v3.to_csv(OUT / "V3_metaregression.csv", index=False, encoding="utf-8-sig")
s3 = {"n": int(len(X)), "tcga_only_training": int(X["tcga_only_training"].sum()), "ext_reported": int(X["ext_reported"].sum()),
      "method_counts": X["method"].value_counts().to_dict(),
      "multivariable": multi[multi.term != "const"][["term", "coef_logitC", "lo", "hi", "p"]].round(4).to_dict("records")}

# ---------------- V4
mp_genes = pd.read_csv(EV / "exploratory/metapcna_genes.csv")["gene"].tolist()
s4 = {}
for c in ["TCGA-PRAD"] + MAIN4:
    expr, clin = data[c]
    S = pd.DataFrame(scores[c]).dropna(axis=1, how="all")
    S = S.loc[:, S.std() > 0]
    R = S.rank().corr()  # Spearman
    iu = np.triu_indices(len(R), 1)
    vals = R.to_numpy()[iu]
    x = expr.loc[[g for g in mp_genes if g in expr.index]].astype(float)
    mp = x.sub(x.mean(1), axis=0).div(x.std(1), axis=0).mean().reindex(S.index)
    rho_mp = S.apply(lambda col: stats.spearmanr(col, mp, nan_policy="omit").statistic)
    entry = {"n_signatures": int(S.shape[1]), "n_patients": int(S.shape[0]), "median_abs_rho": float(np.median(np.abs(vals))),
             "share_abs_rho_ge_0.5": float((np.abs(vals) >= 0.5).mean()), "share_rho_ge_0.7": float((vals >= 0.7).mean()),
             "median_rho_with_metapcna": float(rho_mp.median()), "share_rho_metapcna_ge_0.5": float((rho_mp >= 0.5).mean())}
    if c == "TCGA-PRAD":
        D = 1 - R.to_numpy()
        np.fill_diagonal(D, 0)
        Z = linkage(squareform(np.clip(D, 0, 2), checks=False), method="average")
        cl = fcluster(Z, t=0.5, criterion="distance")
        cs = pd.Series(cl, index=R.index)
        sizes = cs.value_counts()
        big = sizes.index[0]
        entry.update({"n_clusters_at_rho_0.5": int(sizes.size), "largest_cluster_size": int(sizes.iloc[0]),
                      "clusters_size_ge_5": int((sizes >= 5).sum()),
                      "largest_cluster_median_rho_metapcna": float(rho_mp[cs[cs == big].index].median())})
        pd.DataFrame({"sig_id": R.index, "cluster": cl, "rho_metapcna": rho_mp.reindex(R.index).to_numpy()}).to_csv(
            OUT / "V4_clusters_TCGA.csv", index=False, encoding="utf-8-sig")
        R.to_csv(OUT / "V4_spearman_TCGA.csv", encoding="utf-8-sig")
    s4[c] = entry

summ = {"V1": s1, "V2": s2, "V3": s3, "V4": s4}
(OUT / "summary_posthoc.json").write_text(json.dumps(summ, ensure_ascii=False, indent=1, default=float), encoding="utf-8")
print(json.dumps(summ, ensure_ascii=False, indent=1, default=float))
