"""E4 (exploratory, pre-specified 2026-11-12): discrimination the papers reported for their own model vs this study.

Usage: python compare_internal.py <project_root>
Input: 03_extraction/review/out_IM*.jsonl (reported C-index / time-dependent AUC per signature, blind extraction),
       06_results/evaluation/primary/signatures.csv (pooled external C) and run_primary/pairs.csv (td-AUC 36 mo).
Rule (pre-specified): "reported training discrimination" = training-set C-index; if none, training 3-year AUC, then
5-year, then 1-year (median if several values at the chosen level). "Reported external" = same priority applied to the
paper's external cohorts, median over cohorts. Nomogram (signature + clinical) values are not used.
Secondary (added, descriptive): reported 3-year AUC vs this study's median 36-month td-AUC.
Writes 06_results/evaluation/exploratory/{internal_reported.csv, summary_E4.json}.
"""
import glob
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from urosig.paths import p as P, im_globs, _V as PATHS  # noqa: E402
from scipy.stats import spearmanr, wilcoxon

ROOT = Path(sys.argv[1])
OUT = P(ROOT, "ev") / "exploratory"
OUT.mkdir(parents=True, exist_ok=True)


def num(v):
    s = str(v).strip().replace("%", "")
    m = re.search(r"\d+(\.\d+)?", s)
    if not m:
        return np.nan
    x = float(m.group(0))
    return x / 100 if x > 1 else x


def yrs(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


rec, status = [], {}
for f in sorted([x for g in im_globs(ROOT) for x in glob.glob(g)]):
    for line in open(f, encoding="utf-8"):
        r = json.loads(line)
        status[r["sig_id"]] = r["status"]
        for m in r.get("metrics", []):
            rec.append({"sig_id": r["sig_id"], "role": m.get("role"), "metric": m.get("metric"), "t": yrs(m.get("time_years")),
                        "value": num(m.get("value")), "cohort": m.get("cohort"), "from_figure": m.get("from_figure")})
R = pd.DataFrame(rec)


def pick(g):
    for metric, t in (("C-index", None), ("AUC", 3.0), ("AUC", 5.0), ("AUC", 1.0)):
        h = g[(g["metric"] == metric) & (g["t"].isna() if t is None else g["t"] == t) & g["value"].between(0, 1)]
        if len(h):
            return float(h["value"].median()), metric if t is None else f"AUC{int(t)}y", len(h)
    return np.nan, "", 0


rows = []
for sid in status:
    g = R[R["sig_id"] == sid] if len(R) else R
    tr = pick(g[g["role"] == "training"]) if len(g) else (np.nan, "", 0)
    ex = pick(g[g["role"] == "external"]) if len(g) else (np.nan, "", 0)
    a3 = g[(g["role"] == "training") & (g["metric"] == "AUC") & (g["t"] == 3.0)]["value"].median() if len(g) else np.nan
    rows.append({"sig_id": sid, "提取状态": status[sid], "reported_training": tr[0], "training_metric": tr[1],
                 "reported_external": ex[0], "external_metric": ex[1], "reported_training_auc3y": a3})
d = pd.DataFrame(rows)
s = pd.read_csv(P(ROOT, "ev") / "primary/signatures.csv", dtype={"sig_id": str})
pairs = pd.read_csv(P(ROOT, "ev") / "run_primary/pairs.csv", dtype={"sig_id": str})
d = d.merge(s[["sig_id", "pooled_c", "pooled_c_lo", "pooled_c_hi", "k_cohorts"]], on="sig_id", how="left")
d = d.merge(pairs.groupby("sig_id")["auc36"].median().rename("our_median_auc36").reset_index(), on="sig_id", how="left")
d["gap_training"] = d["reported_training"] - d["pooled_c"]
d["gap_external"] = d["reported_external"] - d["pooled_c"]
d["reported_above_our_ci"] = d["reported_training"] > d["pooled_c_hi"]
d.to_csv(OUT / "internal_reported.csv", index=False, encoding="utf-8-sig")


def block(x, gap, ref="pooled_c"):
    x = x.dropna(subset=[gap])
    if len(x) < 3:
        return {"n": int(len(x))}
    rho = spearmanr(x[gap.replace("gap_", "reported_")] if gap.startswith("gap_") else x[gap], x[ref]).statistic
    return {"n": int(len(x)), "median_reported": float(x[gap.replace("gap_", "reported_")].median()),
            "median_ours": float(x[ref].median()), "median_gap": float(x[gap].median()),
            "gap_iqr": [float(x[gap].quantile(.25)), float(x[gap].quantile(.75))],
            "share_gap_positive": float((x[gap] > 0).mean()), "wilcoxon_p": float(wilcoxon(x[gap]).pvalue),
            "spearman_reported_vs_ours": float(rho)}


summ = {"n_signatures": int(len(d)), "status_counts": d["提取状态"].value_counts().to_dict(),
        "training_metric_used": d["training_metric"].value_counts().to_dict(),
        "training_all": block(d, "gap_training"),
        "training_C_only": block(d[d["training_metric"] == "C-index"], "gap_training"),
        "training_AUC_only": block(d[d["training_metric"].str.startswith("AUC")], "gap_training"),
        "external_all": block(d, "gap_external"),
        "share_reported_training_above_our_ci": float(d.loc[d["reported_training"].notna(), "reported_above_our_ci"].mean())}
x = d.dropna(subset=["reported_training_auc3y", "our_median_auc36"])
summ["secondary_auc3y_vs_our_auc36"] = {"n": int(len(x)), "median_reported": float(x["reported_training_auc3y"].median()),
                                         "median_ours": float(x["our_median_auc36"].median()),
                                         "median_gap": float((x["reported_training_auc3y"] - x["our_median_auc36"]).median())} if len(x) else {}
(OUT / "summary_E4.json").write_text(json.dumps(summ, ensure_ascii=False, indent=1), encoding="utf-8")
print(json.dumps(summ, ensure_ascii=False, indent=1))
