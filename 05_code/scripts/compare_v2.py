"""Compare the two independent v2 extractions (X, Y) paper by paper.

Usage: python compare_v2.py <v2_dir>
Status is mapped to a class: main (coefficients from text/table/supplement), figure (read from figure), hr (HR only),
genes (genes only), nonlinear, exclude (all non-* statuses), unclear (无法判断).
Agreement levels: 完全一致 (same class, same gene set, coefficients equal to 1e-6 after normalising signs),
同类别不同内容 (same class, genes or coefficients differ), 类别不同.
Writes compare_v2.csv and 比对表_v2.xlsx (one row per PMID, both extractions side by side).
"""
import glob
import json
import re
import sys
from pathlib import Path

import pandas as pd

D = Path(sys.argv[1])


def cls(s):
    s = s or ""
    if s.startswith("已提取（从图片"):
        return "figure"
    if s.startswith("已提取"):
        return "main"
    if s.startswith("仅 HR"):
        return "hr"
    if s.startswith("仅基因"):
        return "genes"
    if s.startswith("非线性"):
        return "nonlinear"
    if s.startswith("无法判断"):
        return "unclear"
    return "exclude"


def num(x):
    x = str(x).strip().replace("−", "-").replace("–", "-").replace("—", "-").replace("‐", "-").replace("‑", "-").replace(" ", "").strip("()")
    try:
        return float(x)
    except ValueError:
        return None


def gset(g):
    return [re.sub(r"\s+", "", x).upper() for x in str(g or "").split(";") if x.strip()]


# preprint copies left in earlier batches; the published version (V15) is the one assessed (preprints_excluded.txt)
SUPERSEDED = {("37034687", "V05"), ("40458374", "V14")}


def load(code):
    rows = {}
    for f in sorted(glob.glob(str(D / "out" / f"out_V*_{code}.jsonl"))):
        b = re.search(r"out_(V\d+)_", f).group(1)
        for line in open(f, encoding="utf-8"):
            if line.strip():
                r = json.loads(line)
                r["batch"] = b
                if (str(r["pmid"]), b) in SUPERSEDED:
                    continue
                rows[str(r["pmid"])] = r
    return rows


X, Y = load("X"), load("Y")
out = []
for p in sorted(set(X) | set(Y)):
    x, y = X.get(p), Y.get(p)
    row = {"PMID": p, "批次": (x or y)["batch"]}
    for k, r in (("X", x), ("Y", y)):
        r = r or {}
        row.update({f"{k}_status": r.get("status", "（缺）"), f"{k}_类别": cls(r.get("status")) if r else "（缺）",
                    f"{k}_genes": r.get("genes", ""), f"{k}_coefs": r.get("coefs", ""), f"{k}_coef_type": r.get("coef_type", ""),
                    f"{k}_source": r.get("source", ""), f"{k}_endpoint": r.get("endpoint", ""), f"{k}_training": r.get("training", ""),
                    f"{k}_validation": r.get("validation", ""), f"{k}_note": r.get("note", "")})
    if not x or not y:
        row["一致性"] = "单方缺失"
    elif row["X_类别"] != row["Y_类别"]:
        row["一致性"] = "类别不同"
    elif row["X_类别"] in ("exclude", "nonlinear", "unclear"):
        row["一致性"] = "完全一致"
    else:
        gx, gy = gset(x.get("genes")), gset(y.get("genes"))
        same_genes = sorted(gx) == sorted(gy)
        same_coef = True
        bothblank = not str(x.get("coefs") or "").strip() and not str(y.get("coefs") or "").strip()
        if row["X_类别"] in ("main", "figure") and same_genes and not bothblank:
            cx = dict(zip(gx, [num(c) for c in str(x.get("coefs", "")).split(";")]))
            cy = dict(zip(gy, [num(c) for c in str(y.get("coefs", "")).split(";")]))
            same_coef = all(cx.get(g) is not None and cy.get(g) is not None and abs(cx[g] - cy[g]) < 1e-6 for g in gx)
        row["基因相同"] = same_genes
        row["系数相同"] = same_coef
        row["一致性"] = "完全一致" if same_genes and same_coef else "同类别不同内容"
    out.append(row)
c = pd.DataFrame(out)
c.to_csv(D / "compare_v2.csv", index=False, encoding="utf-8-sig")
with pd.ExcelWriter(D / "比对表_v2.xlsx", engine="openpyxl") as w:
    c.to_excel(w, sheet_name="逐篇比对", index=False)
    pd.crosstab(c["X_类别"], c["Y_类别"]).to_excel(w, sheet_name="类别交叉表")
print(pd.crosstab(c["X_类别"], c["Y_类别"]).to_string())
print(c["一致性"].value_counts().to_string())
