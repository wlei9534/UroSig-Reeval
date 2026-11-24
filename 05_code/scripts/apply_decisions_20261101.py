"""Apply Wei's 2026-11-01 decisions after the blind review (runs after apply_review_to_workbook.py).

Usage: python apply_decisions_20261101.py <input_workbook> <output_workbook>
"""
import sys

import pandas as pd
from openpyxl import load_workbook

INP, OUT = sys.argv[1:3]
STAMP = "当前状态（2026-10-31）"
DEC = {  # PMID -> (status override or None, 建议用途 override or None, decision text)
    "42190028": (None, "排除（理由见状态）", "Wei：12 基因探索性评分不算，维持排除"),
    "32264916": (None, None, "Wei：DFS 22 基因与 OS 5 基因两个模型都纳入；OS 模型见 03_extraction/extra_signatures.jsonl（32264916_OS），文献层面汇总按同一篇聚类"),
    "40463716": ("非线性模型（作者主推 RSF；4 基因 Cox 未用于算分，Wei 决定排除）", "排除：非线性模型无法按加权和复现", "Wei：排除"),
    "33336591": (None, "主分析候选（系数来自正文/表格/补充材料）", "Wei：Figure 1B 数值在 PDF 文本层可检索，算正文数值，进主分析"),
}
sheets = pd.read_excel(INP, sheet_name=None, dtype=str)
per = sheets["逐篇"].fillna("")
per["Wei 决定（2026-11-01）"] = per["PMID"].map({k: v[2] for k, v in DEC.items()}).fillna("")
for p, (st, use, _) in DEC.items():
    k = per["PMID"] == p
    assert k.sum() == 1, p
    if st:
        per.loc[k, STAMP] = st
        per.loc[k, "复核后状态"] = st
    if use:
        per.loc[k, "建议用途"] = use
sheets["逐篇"] = per
sheets["说明"] = pd.concat([sheets["说明"], pd.DataFrame({"说明": [
    "", "2026-11-01 Wei 决定：AI 预提取 17 篇按子代理数据；精读 15 篇按更正后数据；其余不一致按各文主推模型取值；"
    "42190028 维持排除；32264916 两个模型都纳入（OS 模型见 extra_signatures.jsonl）；40463716 排除；33336591 算正文数值进主分析；"
    "其余条目专家意见与复核一致。"]})], ignore_index=True)
summ = per[STAMP].value_counts().rename_axis("当前状态").reset_index(name="篇数")
use = per["建议用途"].replace("", "（未提取）").value_counts().rename_axis("建议用途").reset_index(name="篇数")
old = sheets["汇总"]
with pd.ExcelWriter(OUT, engine="openpyxl") as w:
    for k, v in sheets.items():
        if k != "汇总":
            v.to_excel(w, sheet_name=k, index=False)
    summ.to_excel(w, sheet_name="汇总", index=False)
    use.to_excel(w, sheet_name="汇总", index=False, startrow=len(summ) + 3)
src, dst = load_workbook(INP), load_workbook(OUT)
for ws in dst.worksheets:  # carry over column widths / header styles from the input
    if ws.title in src.sheetnames:
        s = src[ws.title]
        for col, dim in s.column_dimensions.items():
            ws.column_dimensions[col].width = dim.width
        ws.freeze_panes = s.freeze_panes
        for a, b in zip(s[1], ws[1]):
            b.font, b.fill, b.alignment = a.font.copy(), a.fill.copy(), a.alignment.copy()
    ws.column_dimensions["A"].width = max(ws.column_dimensions["A"].width or 10, 12)
dst.save(OUT)
print(use.to_string(index=False))
