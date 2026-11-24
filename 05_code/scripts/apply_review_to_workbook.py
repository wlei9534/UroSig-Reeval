"""Bring the blind-review outcome into the extraction workbook (last step of the rebuild chain).

Usage: python apply_review_to_workbook.py <input_workbook> <comparison_xlsx> <output_workbook>
Adds to 逐篇 the reconciled status/genes/coefficients and the adjudication for the 287 reviewed papers,
updates 当前状态 for AI pre-extracted papers the review showed must be excluded, and copies the comparison
as a sheet. Rebuild chain: verify_manual -> merge_manual_extraction -> verify_pdf53 -> merge_pdf53
-> apply_decisions_20261031 -> compare_review -> this script.
"""
import sys

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill

INP, CMP, OUT = sys.argv[1:4]
STAMP = "当前状态（2026-10-31）"
sheets = pd.read_excel(INP, sheet_name=None, dtype=str)
cmp = pd.read_excel(CMP, sheet_name="对比", dtype=str).fillna("")
c = cmp.set_index("PMID")
per = sheets["逐篇"].fillna("")

per["盲法复核（2026-11-01）"] = per["PMID"].map(c["两者一致性"]).fillna("")
per["复核裁定（谁对）"] = per["PMID"].map(c["Claude 对照原文的裁定（谁对）"]).fillna("")
per["复核裁定依据"] = per["PMID"].map(c["裁定依据"]).fillna("")
per["复核后状态"] = per["PMID"].map(c["复核后状态"]).fillna("")
per["复核后基因"] = per["PMID"].map(c["复核后基因"]).fillna("")
per["复核后系数"] = per["PMID"].map(c["复核后系数"]).fillna("")

# AI pre-extracted papers: the reconciled status replaces "已预提取，待人工核对"
ai = per[STAMP].str.startswith("已预提取") & per["复核后状态"].ne("")
per.loc[ai, STAMP] = per.loc[ai, "复核后状态"] + "（AI 预提取经盲法复核）"
excl = ai & ~per["复核后状态"].str.startswith("已提取")
per.loc[ai & ~excl, "建议用途"] = "主分析候选（AI 预提取经盲法复核；待专家确认）"
per.loc[excl, "建议用途"] = "排除（理由见状态）"

summ = per[STAMP].value_counts().rename_axis("当前状态").reset_index(name="篇数")
use_summ = per["建议用途"].replace("", "（未提取）").value_counts().rename_axis("建议用途").reset_index(name="篇数")
rev_summ = cmp.groupby(["集合", "两者一致性"]).size().unstack(fill_value=0).reset_index()
adj_summ = cmp[cmp["两者一致性"] == "不一致"]["Claude 对照原文的裁定（谁对）"].value_counts().rename_axis("不一致条目的裁定").reset_index(name="篇数")
notes = sheets["说明"]["说明"].tolist() + [
    "",
    "2026-11-01 盲法复核：287 篇（AI 预提取 81 + 精读 206）由子代理只看原文与补充材料独立重提；系数只印在图上的 23 篇另给原始图片补读。",
    "逐篇页新增“盲法复核”“复核裁定”“复核后状态/基因/系数”列；“复核对比”页为四方对比表（专家 A/B 列待填）。",
    "复核后值 = Claude 对照原文裁定后的值：AI 预提取 81 篇中 17 篇按子代理更正；精读 206 篇中 15 篇更正（见 manual_extraction*.jsonl 备注，原记录存 *.before_review.jsonl）。",
    "“待专家判断”的条目复核后值暂按提取者1原值保留，最终以专家意见为准。",
]

with pd.ExcelWriter(OUT, engine="openpyxl") as w:
    pd.DataFrame({"说明": notes}).to_excel(w, sheet_name="说明", index=False)
    per.to_excel(w, sheet_name="逐篇", index=False)
    for k in ("人工精读提取", "人工精读长表", "待人工核对", "基因-系数长表（AI预提取）"):
        sheets[k].to_excel(w, sheet_name=k, index=False)
    cmp.to_excel(w, sheet_name="复核对比", index=False)
    summ.to_excel(w, sheet_name="汇总", index=False, startrow=0)
    use_summ.to_excel(w, sheet_name="汇总", index=False, startrow=len(summ) + 3)
    r0 = len(summ) + len(use_summ) + 6
    rev_summ.to_excel(w, sheet_name="汇总", index=False, startrow=r0)
    adj_summ.to_excel(w, sheet_name="汇总", index=False, startrow=r0 + len(rev_summ) + 3)

wb = load_workbook(OUT)
WIDE = {"标题": 50, "基因": 40, "系数": 40, "人工精读基因": 40, "人工精读系数": 40, "来源（原文位置）": 50,
        "人工精读来源": 50, "备注": 60, "人工精读备注": 60, "核对什么": 40, "说明": 120, "复核后基因": 40,
        "复核后系数": 40, "复核裁定依据": 50, "差异说明": 50, "裁定依据": 50}
for ws in wb.worksheets:
    ws.freeze_panes = "A2" if ws.title != "说明" else None
    for cell in ws[1]:
        cell.font = Font(bold=True)
        cell.fill = PatternFill("solid", fgColor="FFF2CC" if str(cell.value).startswith(("专家", "最终")) else "DDEBF7")
        cell.alignment = Alignment(wrap_text=True, vertical="top")
        ws.column_dimensions[cell.column_letter].width = WIDE.get(str(cell.value), 14)
    for row in ws.iter_rows(min_row=2):
        for cell in row:
            if isinstance(cell.value, str) and len(cell.value) > 30:
                cell.alignment = Alignment(wrap_text=True, vertical="top")
    if ws.max_row > 1 and ws.title != "汇总":
        ws.auto_filter.ref = ws.dimensions
wb.save(OUT)
print("saved", OUT)
print(summ.to_string(index=False))
