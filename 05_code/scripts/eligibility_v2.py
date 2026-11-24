"""Signature x cohort eligibility for the active lock (no outcome data read). Same rule as v1 build_eligibility.py.

Usage: python eligibility_v2.py <project_root>
Inputs : coverage (gene_coverage.py), leakage matrix (leakage_v2.py), equal-weight directions (lock_extraction_v2.py)
Eligible (external cohorts; TCGA-PRAD is reference only): not used by the paper AND gene coverage >= threshold
(and weight coverage >= threshold for coefficient signatures) AND, for equal_weight, >= 1 gene with a direction.
Also writes the original-endpoint class (S6) with the v1 first-term rule, using the extractors' endpoint field.
"""
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from urosig.paths import p as P  # noqa: E402

ROOT = Path(sys.argv[1])
COH = ["GSE54460", "GSE70769", "GSE70768", "MSKCC", "GSE116918", "TCGA"]
CMAP = {"MSKCC_GSE21032": "MSKCC", "TCGA-PRAD": "TCGA"}
lk = pd.read_csv(P(ROOT, "leak"), dtype=str).fillna("")
dr = pd.read_csv(P(ROOT, "dirs"), dtype={"sig_id": str})
n_dir = dr[dr["方向"].astype(int) != 0].groupby("sig_id").size()
cov = pd.read_csv(P(ROOT, "cov"), dtype={"sig_id": str, "PMID": str})
cov["coh"] = cov["cohort"].map(lambda c: CMAP.get(c, c))
cov = cov.merge(lk[["sig_id"] + COH], on="sig_id", how="left")
assert cov[COH].notna().all().all(), "signature without leakage row"
cov["原文用过"] = cov.apply(lambda r: r[r["coh"]] == "用过", axis=1)
cov["有可评分基因"] = cov.apply(lambda r: r["分析角色"] != "equal_weight" or n_dir.get(r["sig_id"], 0) > 0, axis=1)
for th in (0.8, 1.0):
    wc = pd.to_numeric(cov["权重覆盖率"], errors="coerce")
    cov[f"可评估_{th}"] = (cov["coh"] != "TCGA") & ~cov["原文用过"] & (cov["基因覆盖率"] >= th) & \
        (wc.isna() | (wc >= th)) & cov["有可评分基因"]
keep = ["sig_id", "PMID", "分析角色", "cohort", "原文用过", "基因项数", "找到", "基因覆盖率", "权重覆盖率", "有可评分基因",
        "可评估_0.8", "可评估_1.0", "缺失基因"]
cov[keep].to_csv(P(ROOT, "elig"), index=False, encoding="utf-8-sig")
summ = []
for th in (0.8, 1.0):
    e = cov[cov[f"可评估_{th}"]]
    for role, g in cov.groupby("分析角色"):
        ids = set(g["sig_id"])
        ok = e[e["sig_id"].isin(ids)].groupby("sig_id").size()
        summ.append({"门槛": th, "分析角色": role, "signature 数": len(ids), "≥1 个可评估外部队列": len(ok),
                     "≥2 个": int((ok >= 2).sum()), "无可评估队列": len(ids) - len(ok), "可评估 signature×队列": int(ok.sum())})
summ = pd.DataFrame(summ)
summ.to_csv(P(ROOT, "elig_summary"), index=False, encoding="utf-8-sig")
print(summ.to_string(index=False))

# original endpoint class (S6): v1 rows as locked, v2 rows from the extractors' endpoint field, first-term rule
BCR = r"BCR|BRFS|bRFS|\bBFS\b|BFFS|biochemical|PSA-free|\bRFS\b|\bDFS\b|recurrence|relapse|复发|生化"
OTHER = r"\bOS\b|overall|\bPFS\b|\bPFI\b|进展|\bMFS\b|metasta|转移|总生存|未明确"
v1e = pd.read_csv(ROOT / "06_results" / "eligibility" / "original_endpoint_class.csv", dtype=str).fillna("")
sig = pd.read_csv(P(ROOT, "sig"), dtype=str).fillna("")
v2 = sig[sig["来源版本"].str.startswith("v2") & (sig["分析角色"] != "excluded")]
rows = []
for _, s in v2.iterrows():
    txt = s["终点"]
    mb, mo = re.search(BCR, txt, re.I), re.search(OTHER, txt, re.I)
    cls = "BCR/RFS" if mb and (not mo or mb.start() < mo.start()) else ("other" if mo else "unclear")
    rows.append({"sig_id": s["sig_id"], "PMID": s["PMID"], "分析角色": s["分析角色"], "原文终点（复核员记录）": txt, "终点类别": cls})
ep = pd.concat([v1e, pd.DataFrame(rows)], ignore_index=True)
ep.to_csv(P(ROOT, "endpoint"), index=False, encoding="utf-8-sig")
print(pd.DataFrame(rows).groupby("终点类别").size().to_string())
