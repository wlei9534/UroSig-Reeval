"""Reconcile the two v2 extractions (X, Y) into one proposed record per paper, and build the review workbook.

Usage: python reconcile_v2.py <v2_dir>
Inputs : compare_v2.csv (from compare_v2.py), adjudication_v2.csv (Claude's source checks and proposals).
Outputs: reconciled_v2.csv (one row per PMID, proposed analysis role) and 裁定表_v2.xlsx (sheets for Wei/expert).

Role mapping follows the v1 lock (signatures_locked_2026-11-01.csv, SAP v2):
  text/table/supplement coefficients -> main            (primary)
  coefficients read from a figure     -> sens_figure     (S7)
  HR only                             -> sens_hr         (S7, ln HR)
  author-defined equal weights (+-1)  -> equal_signed    (S8)
  genes only, directions known        -> equal_weight    (S8)
  gene pairs with weights             -> main; equal-vote pairs -> equal_signed (S8)
  non-linear / excluded / unclear     -> excluded (reason kept)
Rule flags mark records whose role depends on a rule Wei and the expert still have to confirm.
"""
import re
import sys
from pathlib import Path

import pandas as pd

D = Path(sys.argv[1])
c = pd.read_csv(D / "compare_v2.csv", dtype=str).fillna("")
adj = pd.read_csv(D / "adjudication_v2.csv", dtype=str).fillna("")

ROLE = {"main": "main", "figure": "sens_figure", "hr": "sens_hr", "genes": "equal_weight",
        "nonlinear": "excluded", "exclude": "excluded", "unclear": "excluded"}

# rule questions (R-codes as in the workbook); regexes run on both extractors' status + note
RULES = {
    "R4 非生存目标/机制研究附带的 signature": r"机制|表型|谱系|干性|缺氧|LVI|药效|分化状态|治疗选择|非生存|附带|side analysis",
    "R5 PCA 评分（载荷未报告）": r"PCA|主成分|PC1",
    "R6 预定义基因集打分（无基因筛选）": r"预定义|MSigDB|GeneCards|现成|既有基因集|predefined",
    "R9 系数作用于二分/区间变量": r"二分后|二分变量|0/1|区间指示|离散化|step function|按中位数二分的",
    "R10 模型含临床变量": r"临床变量|加 Gleason|\+ ?Gleason",
    "R13 非肿瘤组织（癌旁/血液）": r"癌旁|TAN|外周血",
    "R14 作者把 HR 当权重写进公式": r"HR 当|把 HR|当系数写进|当权重",
    "R15 原文系数/HR/方向自相矛盾": r"对不上|自相矛盾|前后矛盾|矛盾",
    "缺补充材料": r"未获取|拿不到|没有补充|不在本批|本批次没有|本地没有|supp_dirs 为空|未下载",
}


def directions(note):
    m = re.search(r"directions[：:]\s*([^\n。]+)", note)
    if not m:
        return ""
    return m.group(1).strip()


rows = []
for _, r in c.iterrows():
    p = r["PMID"]
    agree = r["一致性"] == "完全一致"
    cls = r["X_类别"] if agree else ""
    pick = "X"
    role = ROLE.get(cls, "") if agree else "待裁定"
    is_pair = "|" in r["X_genes"] and "|" in r["Y_genes"]
    if agree and cls == "main" and is_pair and not r["X_coefs"].strip():
        role = "equal_signed"
    text = " ".join(r[k] for k in ("X_status", "Y_status", "X_note", "Y_note"))
    flags = [k for k, rx in RULES.items() if re.search(rx, text)] if (role not in ("excluded",) or not agree) else []
    rec = {"PMID": p, "批次": r["批次"], "一致性": r["一致性"], "X_status": r["X_status"], "Y_status": r["Y_status"],
           "拟定角色": role, "采用": pick if agree else "", "规则标记": "；".join(flags),
           "genes": r["X_genes"], "coefs": r["X_coefs"], "coef_type": r["X_coef_type"], "source": r["X_source"],
           "endpoint": r["X_endpoint"], "training": r["X_training"], "validation": r["X_validation"],
           "directions": directions(r["X_note"]) or directions(r["Y_note"]),
           "X_note": r["X_note"], "Y_note": r["Y_note"], "裁定": "", "裁定依据": "", "需人工": ""}
    rows.append(rec)
R = pd.DataFrame(rows).set_index("PMID")

# apply adjudication rows (later rows override earlier ones)
for _, a in adj.iterrows():
    p = a["PMID"]
    if p not in R.index:
        continue
    R.loc[p, "裁定"] = a["裁定"]
    R.loc[p, "裁定依据"] = a["依据"]
    R.loc[p, "需人工"] = a["需人工"]
    if a.get("角色", ""):
        R.loc[p, "拟定角色"] = a["角色"]
    if a.get("采用", "") in ("X", "Y"):
        k = a["采用"]
        R.loc[p, "采用"] = k
        row = c.set_index("PMID").loc[p]
        for f in ("genes", "coefs", "coef_type", "source", "endpoint", "training", "validation"):
            R.loc[p, f] = row[f"{k}_{f}"]
    if a.get("系数覆盖", ""):
        R.loc[p, "coefs"] = a["系数覆盖"]

R = R.reset_index()
R.to_csv(D / "reconciled_v2.csv", index=False, encoding="utf-8-sig")
print(R["拟定角色"].value_counts().to_string())
print("rule-flagged (non-excluded):", (R["规则标记"] != "").sum())
