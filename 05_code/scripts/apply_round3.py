"""Apply the round-3 adjudication (main-analysis candidates with any extractor disagreement) to the workbook.

Usage: python apply_round3.py <project_root> <input_workbook> <round3_form_xlsx> <output_workbook> <output_form_xlsx>
Inputs : 03_extraction/review/out_R3a.jsonl, out_R3b.jsonl (sub-agent round-3 rulings, each final coefficient
         checked verbatim against the source text by Claude except 35414328 [refit] and 39760809 [typo fix]).
Effect : 逐篇 复核后基因/系数 take the round-3 final values (higher precision where the paper prints it; BAK -> BAK1
         for 35669428); new columns 第三轮核定 / 第三轮依据; extra_signatures.jsonl updated for 32264916_OS;
         the round-3 form gets the sub-agent columns filled in.
"""
import json
import sys
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook

ROOT, INP, FORM, OUT, OUTFORM = Path(sys.argv[1]), *sys.argv[2:6]
X = ROOT / "03_extraction"
R = {}
for f in ("out_R3a.jsonl", "out_R3b.jsonl"):
    for l in open(X / "review" / f, encoding="utf-8"):
        r = json.loads(l)
        R[r["pmid"]] = r
assert len(R) == 50, len(R)

UNSURE = {"35222640": "NCF1 系数前原文缺正负号（公式印作“LRMP ∗0 .5567643∗ NCF1”），暂按 +0.5567643；分析时做正、负两种符号的敏感性分析"}

sheets = pd.read_excel(INP, sheet_name=None, dtype=str)
per = sheets["逐篇"].fillna("")
per["第三轮核定（2026-11-01）"] = ""
per["第三轮依据"] = ""
for p, r in R.items():
    if p.endswith("_OS"):
        continue
    k = per["PMID"] == p
    assert k.sum() == 1, p
    per.loc[k, "第三轮核定（2026-11-01）"] = r["verdict"]
    per.loc[k, "第三轮依据"] = (UNSURE.get(p, "") + "；" if p in UNSURE else "") + r["points"]
    if r["final_genes"]:
        g, c = r["final_genes"].split(";"), r["final_coefs"].split(";") if r["final_coefs"] else []
        assert not c or len(g) == len(c), p
        per.loc[k, "复核后基因"] = r["final_genes"]
        per.loc[k, "复核后系数"] = r["final_coefs"]
sheets["逐篇"] = per
sheets["说明"] = pd.concat([sheets["说明"], pd.DataFrame({"说明": [
    "", "2026-11-01 第三轮复核（主分析候选中两版提取有差异的 49 篇 + 32264916_OS）：采用值正确 40、需修改 9（8 篇改用原文最高精度、"
    "35669428 基因名 BAK→BAK1）、无法确定 1（35222640 NCF1 系数符号，按 + 记录并做正负敏感性分析）。复核后系数已按原文最高精度更新，"
    "每个系数均在原文逐字找到（35414328 为重拟合、39760809 为已确认的笔误纠正除外）。"]})], ignore_index=True)

with pd.ExcelWriter(OUT, engine="openpyxl") as w:
    for k, v in sheets.items():
        v.to_excel(w, sheet_name=k, index=False)
src, dst = load_workbook(INP), load_workbook(OUT)
for ws in dst.worksheets:
    if ws.title in src.sheetnames:
        s = src[ws.title]
        for col, dim in s.column_dimensions.items():
            ws.column_dimensions[col].width = dim.width
        ws.freeze_panes = s.freeze_panes
        for a, b in zip(s[1], ws[1]):
            b.font, b.fill, b.alignment = a.font.copy(), a.fill.copy(), a.alignment.copy()
dst.save(OUT)

# extra signature 32264916_OS
es = [json.loads(l) for l in open(X / "extra_signatures.jsonl", encoding="utf-8")]
for e in es:
    r = R.get(e["sig_id"])
    if r:
        e["genes"], e["coefs"] = r["final_genes"], r["final_coefs"]
        e["note"] += f"；第三轮核定：{r['verdict']}"
open(X / "extra_signatures.jsonl", "w", encoding="utf-8").write("".join(json.dumps(e, ensure_ascii=False) + "\n" for e in es))

# round-3 form for expert C
fm = pd.read_excel(FORM, sheet_name=None, dtype=str)
t = fm["第三轮复核"].fillna("")
t["子代理第三轮核定"] = t["PMID"].map({p: r["verdict"] for p, r in R.items()})
t["子代理第三轮最终基因"] = t["PMID"].map({p: r["final_genes"] for p, r in R.items()})
t["子代理第三轮最终系数"] = t["PMID"].map({p: r["final_coefs"] for p, r in R.items()})
t["子代理第三轮依据"] = t["PMID"].map({p: (UNSURE.get(p, "") + "；" if p in UNSURE else "") + r["points"] for p, r in R.items()})
fm["第三轮复核"] = t
with pd.ExcelWriter(OUTFORM, engine="openpyxl") as w:
    for k, v in fm.items():
        v.to_excel(w, sheet_name=k, index=False)
src, dst = load_workbook(FORM), load_workbook(OUTFORM)
for ws in dst.worksheets:
    s = src[ws.title]
    for col, dim in s.column_dimensions.items():
        ws.column_dimensions[col].width = dim.width
    ws.freeze_panes = s.freeze_panes
    for a, b in zip(s[1], ws[1]):
        b.font, b.fill, b.alignment = a.font.copy(), a.fill.copy(), a.alignment.copy()
    from openpyxl.styles import Alignment
    for row in ws.iter_rows(min_row=2):
        for c in row:
            if isinstance(c.value, str) and len(c.value) > 30:
                c.alignment = Alignment(wrap_text=True, vertical="top")
dst.save(OUTFORM)
print(pd.Series([r["verdict"] for r in R.values()]).value_counts().to_string())
