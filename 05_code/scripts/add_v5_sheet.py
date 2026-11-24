"""Add the V5 (calibration slope) sheet to the v2 results workbook. Usage: python add_v5_sheet.py <project_root>"""
import sys
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
from openpyxl.styles import Font

ROOT = Path(sys.argv[1]); EV = ROOT / "06_results" / "evaluation_v2" / "posthoc"
f = ROOT / "06_results" / "UroSig_结果表_v2_2026-11-17.xlsx"
v = pd.read_csv(EV / "V5_calibration_slope_pairs.csv", dtype={"sig_id": str})


def pool(d):
    w = 1 / d["slope_se"] ** 2; b = (w * d["slope"]).sum() / w.sum(); se = np.sqrt(1 / w.sum())
    return pd.Series({"slope": b, "lo": b - 1.96 * se, "hi": b + 1.96 * se})


def summ(name, slope, hi, n_label, lp=None):
    return [name, n_label, round(slope.median(), 3), round(slope.quantile(.25), 3), round(slope.quantile(.75), 3), round(100 * (hi < 1).mean(), 1), round(100 * (slope <= 0).mean(), 1), "" if lp is None else round(lp, 3)]


rows = []
for c in ["GSE54460", "GSE70769", "GSE70768", "MSKCC_GSE21032"]:
    d = v[v.cohort == c]; rows.append(summ(c, d["slope"], d["slope_hi"], len(d), d["lp_sd"].median()))
allp = v.groupby("sig_id")[["slope", "slope_se"]].apply(pool); no54 = v[v.cohort != "GSE54460"].groupby("sig_id")[["slope", "slope_se"]].apply(pool)
rows.append(summ("全部队列合并（逆方差）", allp["slope"], allp["hi"], len(allp)))
rows.append(summ("去掉 GSE54460 后合并", no54["slope"], no54["hi"], len(no54)))
wb = openpyxl.load_workbook(f)
if "V5_校准斜率" in wb.sheetnames:
    del wb["V5_校准斜率"]
ws = wb.create_sheet("V5_校准斜率")
notes = ["V5 校准斜率（事后分析，2026-11-18 新增；2026-11-19 更正表述）",
         "定义：在未标准化 log2 表达上计算原文线性预测值，作为唯一协变量拟合 Cox，系数即校准斜率（理想值 1）。",
         "注意：合并值主要由 GSE54460 决定（逆方差权重中位数约 82%）；该队列线性预测值的 SD 约为其余队列的 2.4 倍，提示平台尺度差异。应逐队列报告，不以合并中位数 0.28 作为单一结论。",
         "来源：06_results/evaluation_v2/posthoc/V5_calibration_slope_pairs.csv"]
for t in notes:
    ws.append([t])
ws.append([])
ws.append(["队列", "signature 数", "斜率中位数", "Q1", "Q3", "95%CI 整体低于 1（%）", "斜率 ≤0（%）", "线性预测值 SD 中位数"])
for c in ws[ws.max_row]:
    c.font = Font(bold=True)
for r in rows:
    ws.append(r)
ws["A1"].font = Font(bold=True); ws.column_dimensions["A"].width = 28
for col in "BCDEFGH":
    ws.column_dimensions[col].width = 20
wb.save(f)
print("sheet added:", [r[:3] + r[5:6] for r in rows])
