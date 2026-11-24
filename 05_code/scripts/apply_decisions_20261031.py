"""Apply Wei's 2026-10-31 decisions on top of the merged workbook (last step of the rebuild chain).

Rebuild chain:  verify_manual.py -> merge_manual_extraction.py (from archive/ai_prefill_PRAD_2026-10-26.xlsx)
                -> verify_pdf53.py -> merge_pdf53.py -> this script
Usage: python apply_decisions_20261031.py <input_workbook> <output_workbook>
"""
import sys

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill

INP, OUT = sys.argv[1], sys.argv[2]
STAMP = "当前状态（2026-10-31）"

STATUS = {  # PMID -> (new status, suggested use)
    "39399052": ("排除：39518911 的预印本（同一 7 基因模型，Wei 决定删除，只保留 39518911）", "排除（理由见状态）"),
    "38924839": ("非前列腺癌（题录阶段误纳入，Wei 确认排除）", "排除（理由见状态）"),
    "35983076": ("非前列腺癌（题录阶段误纳入，Wei 确认排除）", "排除（理由见状态）"),
    "39511540": ("非前列腺癌（题录阶段误纳入，Wei 确认排除）", "排除（理由见状态）"),
}
CHECKED = {  # 待人工核对 rows resolved by Wei's reviewers
    "33338276": "Wei 安排专人核对：补充表表头标反成立，CliSigGlu-22 / GleSigGlu-26 基因归属按基因个数判定无误",
    "36647875": "Wei 安排专人核对：PAQR6 系数原文印作 396，确认为 0.396 的笔误（使用时按 0.396）",
    "39760809": "Wei 安排专人核对：“MATN * 30.449” 确认为 MATN3 0.449",
    "39399052": "Wei 决定：39399052 为 39518911 的预印本，删除 39399052，只保留 39518911",
    "39518911": "Wei 决定：保留本篇（正式发表版），删除其预印本 39399052",
}
NOTES = [
    "",
    "2026-10-31 Wei 的决定：",
    "系数来源分级：正文/表格/补充材料（含 35414328 由作者数据重拟合）→ 主分析；图片读数、HR 取对数换算 → 仅敏感性分析；"
    "只有基因名、没有系数 → 等权评分（方向取原文 HR 或上下调方向）。",
    "35414328 进入主分析。39399052（39518911 的预印本）删除。38924839、35983076、39511540 题录阶段误纳入，确认排除。",
    "36647875 的 PAQR6 系数在本表中按原文记为 396，分析时按 0.396 使用。",
]

sheets = pd.read_excel(INP, sheet_name=None, dtype=str)
per = sheets["逐篇"]
for p, (st, use) in STATUS.items():
    k = per["PMID"] == p
    assert k.sum() == 1, p
    per.loc[k, STAMP] = st
    per.loc[k, "建议用途"] = use
pr = sheets["人工精读提取"]
for p, (st, use) in STATUS.items():
    k = pr["PMID"] == p
    pr.loc[k, "人工精读状态"] = st
    pr.loc[k, "建议用途"] = use
todo = sheets["待人工核对"]
todo["核对结果（2026-10-31）"] = todo["PMID"].map(CHECKED).fillna("")
summ = per[STAMP].value_counts().rename_axis("当前状态").reset_index(name="篇数")
use_summ = per["建议用途"].fillna("（AI 预提取以外、无人工精读记录）").value_counts().rename_axis("建议用途").reset_index(name="篇数")
note = pd.DataFrame({"说明": sheets["说明"]["说明"].tolist() + NOTES})

with pd.ExcelWriter(OUT, engine="openpyxl") as w:
    note.to_excel(w, sheet_name="说明", index=False)
    per.to_excel(w, sheet_name="逐篇", index=False)
    pr.to_excel(w, sheet_name="人工精读提取", index=False)
    sheets["人工精读长表"].to_excel(w, sheet_name="人工精读长表", index=False)
    todo.to_excel(w, sheet_name="待人工核对", index=False)
    sheets["基因-系数长表（AI预提取）"].to_excel(w, sheet_name="基因-系数长表（AI预提取）", index=False)
    summ.to_excel(w, sheet_name="汇总", index=False, startrow=0)
    use_summ.to_excel(w, sheet_name="汇总", index=False, startrow=len(summ) + 3)

wb = load_workbook(OUT)
WIDE = {"标题": 50, "基因": 40, "系数": 40, "人工精读基因": 40, "人工精读系数": 40, "来源（原文位置）": 50,
        "人工精读来源": 50, "备注": 60, "人工精读备注": 60, "核对什么": 40, "说明": 120, "公式原文（截取）": 40,
        "核对结果（2026-10-31）": 40}
for ws in wb.worksheets:
    ws.freeze_panes = "A2" if ws.title != "说明" else None
    for c in ws[1]:
        c.font = Font(bold=True)
        c.fill = PatternFill("solid", fgColor="DDEBF7")
        c.alignment = Alignment(wrap_text=True, vertical="top")
        ws.column_dimensions[c.column_letter].width = WIDE.get(str(c.value), 14)
    for row in ws.iter_rows(min_row=2):
        for c in row:
            if isinstance(c.value, str) and len(c.value) > 30:
                c.alignment = Alignment(wrap_text=True, vertical="top")
    if ws.max_row > 1:
        ws.auto_filter.ref = ws.dimensions
wb.save(OUT)
print("saved", OUT)
print(use_summ.to_string(index=False))
