"""Merge the two blind title/abstract screeners (A, B) for the 719 records not screened before.

Usage: python merge_ta.py <ta_dual_dir>
Rule (decision_log 2026-11-12): both exclude -> excluded; both include/uncertain -> full text;
one exclude + one include/uncertain -> adjudication by Wei/expert.
Writes merged_TA.csv and 待裁定_TA.xlsx (conflicts with title/abstract and both reasons).
"""
import glob
import json
import sys
from pathlib import Path

import pandas as pd

D = Path(sys.argv[1])
recs = pd.concat([pd.read_csv(f, dtype=str) for f in sorted(D.glob("batch_TA*.csv"))], ignore_index=True)


def load(r):
    rows = {}
    for f in sorted(glob.glob(str(D / f"out_{r}_TA*.jsonl"))):
        for line in open(f, encoding="utf-8"):
            line = line.strip()
            if line:
                x = json.loads(line)
                rows[str(x["PMID"])] = x
    return pd.DataFrame(rows.values()).rename(columns={"decision": f"{r}_决定", "code": f"{r}_代码", "reason": f"{r}_理由"})


A, B = load("A"), load("B")
m = recs.merge(A[["PMID", "A_决定", "A_代码", "A_理由"]], on="PMID", how="left").merge(
    B[["PMID", "B_决定", "B_代码", "B_理由"]], on="PMID", how="left")
assert m["A_决定"].notna().all() and m["B_决定"].notna().all(), "screening incomplete"
inc = lambda d: d in ("纳入", "不确定")  # noqa: E731
m["结果"] = [("进入全文" if inc(a) and inc(b) else "排除" if (a == "排除" and b == "排除") else "待裁定")
            for a, b in zip(m["A_决定"], m["B_决定"])]
m.to_csv(D / "merged_TA.csv", index=False, encoding="utf-8-sig")
conf = m[m["结果"] == "待裁定"][["编号", "PMID", "年份", "期刊", "标题", "摘要", "A_决定", "A_代码", "A_理由", "B_决定", "B_代码", "B_理由"]].copy()
conf["裁定（纳入/排除）"] = ""
conf["裁定理由"] = ""
with pd.ExcelWriter(D / "待裁定_TA.xlsx", engine="openpyxl") as w:
    conf.to_excel(w, sheet_name="待裁定", index=False)
    ws = w.book["待裁定"]
    from openpyxl.styles import Alignment, PatternFill
    widths = {"A": 8, "B": 11, "C": 7, "D": 18, "E": 50, "F": 90}
    for col, wd in widths.items():
        ws.column_dimensions[col].width = wd
    for col in "GHIJKLMN":
        ws.column_dimensions[col].width = 16
    for row in ws.iter_rows(min_row=2):
        for c in row:
            c.alignment = Alignment(wrap_text=True, vertical="top")
    for c in ws["M"][1:] + ws["N"][1:]:
        c.fill = PatternFill("solid", fgColor="FFF2CC")
kappa_tab = pd.crosstab(m["A_决定"], m["B_决定"])
po = (m["A_决定"].map(inc) == m["B_决定"].map(inc)).mean()
pa, pb = m["A_决定"].map(inc).mean(), m["B_决定"].map(inc).mean()
pe = pa * pb + (1 - pa) * (1 - pb)
print(kappa_tab.to_string())
print(m["结果"].value_counts().to_string())
print(f"agreement (include/uncertain vs exclude) = {po:.3f}, Cohen kappa = {(po - pe) / (1 - pe):.3f}")
