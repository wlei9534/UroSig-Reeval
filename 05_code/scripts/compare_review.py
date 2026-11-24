"""Side-by-side comparison of extractor 1 (Claude, main extraction) and extractor 2 (blind sub-agent review).

Usage: python compare_review.py <project_root> <workbook> <output_xlsx>
Inputs : <workbook> 逐篇 sheet (AI pre-extraction genes/coefs for the 81 "已预提取" papers),
         03_extraction/manual_extraction.jsonl + manual_extraction_pdf53.jsonl (manual extraction),
         03_extraction/review/out_*.jsonl (blind review), review/all_papers.csv (paper list)
Output : one row per paper, with empty columns for expert A and expert B and the final ruling.
"""
from __future__ import annotations

import glob
import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT, WB, OUT = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
X = ROOT / "03_extraction"
papers = pd.read_csv(X / "review" / "all_papers.csv", dtype=str).fillna("")
per = pd.read_excel(WB, sheet_name="逐篇", dtype=str).fillna("").set_index("PMID")

mine = {}
for _, r in papers[papers["set"] == "AI预提取81"].iterrows():
    w = per.loc[r.PMID]
    mine[r.PMID] = {"status": "已提取（AI 预提取）", "genes": w["基因"], "coefs": w["系数"], "source": w["来源"], "note": w["提示"]}
final = {}  # reconciled values after Claude's adjudication against the source text
for f in ("manual_extraction", "manual_extraction_pdf53"):
    before = X / f"{f}.before_review.jsonl"  # extractor 1 as originally recorded (pre-review)
    for l in open(before if before.exists() else X / f"{f}.jsonl", encoding="utf-8"):
        r = json.loads(l)
        mine[r["pmid"]] = {k: r.get(k, "") for k in ("status", "genes", "coefs", "source", "note")}
    for l in open(X / f"{f}.jsonl", encoding="utf-8"):
        r = json.loads(l)
        final[r["pmid"]] = {k: r.get(k, "") for k in ("status", "genes", "coefs")}
ADJ = {}
for f in ("adjudication_AI.json", "adjudication_M.json"):
    ADJ.update(json.load(open(X / "review" / f, encoding="utf-8")))
# Wei's later decisions override the status shown for extractor 1
for p in ("38924839", "35983076", "39511540"):
    mine[p]["status"] = "非前列腺癌（Wei 确认排除）"
if "39399052" in mine:
    mine["39399052"]["status"] = "排除：39518911 的预印本（Wei 决定）"

rev = {}
files = sorted(glob.glob(str(X / "review" / "out_*.jsonl")))
files = [f for f in files if not f.endswith("out_FIG.jsonl")] + [f for f in files if f.endswith("out_FIG.jsonl")]
FIGSET = set()
for f in files:  # the figure-enabled re-review (out_FIG) is read last and supersedes the no-figure pass
    for l in open(f, encoding="utf-8"):
        if l.strip():
            r = json.loads(l)
            rev[str(r["pmid"])] = r
            if f.endswith("out_FIG.jsonl"):
                FIGSET.add(str(r["pmid"]))


def cat(status: str) -> str:
    s = status or ""
    if s.startswith("已提取（基因对"):
        return "有系数-基因对"
    if s.startswith("已提取（从图片") or s.startswith("已提取（由 HR") or s.startswith("仅 HR"):
        return "有系数-图片/HR"
    if s.startswith("已提取"):
        return "有系数"
    if s.startswith("仅基因") or s.startswith("系数在补充材料") or s.startswith("系数仅在"):
        return "仅基因/无系数"
    if s.startswith("非线性") or s.startswith("非锁定"):
        return "非线性"
    if s.startswith("无法判断"):
        return "无法判断"
    if s.startswith("非") or s.startswith("排除"):
        return "排除"
    return "其他"


def norm(g: str) -> str:
    return re.sub(r"[\s_‐‑–—-]", "", g).upper()


def split(s: str):
    return [x.strip() for x in (s or "").split(";") if x.strip()]


def num(c: str):
    try:
        return float(str(c).replace("−", "-").replace("–", "-"))
    except ValueError:
        return None


def same_num(a: float, b: float) -> bool:
    """Equal, or equal after rounding to the shorter of the two printed precisions (e.g. 0.4025 vs 0.402463)."""
    if abs(a - b) <= 1e-9:
        return True
    return False


def compare(a: dict, b: dict) -> tuple[str, str]:
    ca, cb = cat(a["status"]), cat(b.get("status", ""))
    if not b:
        return "子代理未完成", ""
    diffs = []
    if ca != cb and not (ca == "有系数-图片/HR" and b.get("status", "").startswith("仅 HR") and a["status"].startswith("已提取（由 HR")):
        diffs.append(f"判定不同：{ca} vs {cb}")
    ga, gb = split(a["genes"]), split(b.get("genes", ""))
    na, nb = {norm(g): g for g in ga}, {norm(g): g for g in gb}
    only_a, only_b = [na[k] for k in na if k not in nb], [nb[k] for k in nb if k not in na]
    if ga or gb:
        if only_a:
            diffs.append("仅提取者1有：" + ",".join(only_a))
        if only_b:
            diffs.append("仅子代理有：" + ",".join(only_b))
    ka = dict(zip([norm(g) for g in ga], split(a["coefs"]))) if a["coefs"] else {}
    kb = dict(zip([norm(g) for g in gb], split(b.get("coefs", "")))) if b.get("coefs") else {}
    if b.get("status", "").startswith("仅 HR") and a["status"].startswith("已提取（由 HR"):
        import math
        hr_bad = []
        for k in set(ka) & set(kb):
            x, y = num(ka[k]), num(kb[k])
            if x is not None and y and y > 0 and abs(x - math.log(y)) > 0.01:
                hr_bad.append(f"{na.get(k, k)}(ln HR {math.log(y):.4f} vs {ka[k]})")
        if hr_bad:
            diffs.append("HR 换算后不同：" + "；".join(hr_bad))
        ka, kb = {}, {}
        if ca == "有系数-图片/HR":
            cb = ca
            diffs = [d for d in diffs if not d.startswith("判定不同")]
    bad, rounding = [], []
    for k in set(ka) & set(kb):
        x, y = num(ka[k]), num(kb[k])
        if x is None or y is None:
            continue
        if same_num(x, y):
            continue
        dx = len(ka[k].split(".")[-1]) if "." in ka[k] else 0
        dy = len(kb[k].split(".")[-1]) if "." in kb[k] else 0
        d = min(dx, dy)
        if round(x, d) == round(y, d):
            rounding.append(na.get(k, k))
        elif abs(x + y) <= 1e-9:
            bad.append(f"{na.get(k, k)} 符号相反({ka[k]} vs {kb[k]})")
        else:
            bad.append(f"{na.get(k, k)}({ka[k]} vs {kb[k]})")
    if bad:
        diffs.append("系数不同：" + "；".join(bad))
    if rounding:
        diffs.append("仅小数位数不同：" + ",".join(rounding))
    if bool(ka) != bool(kb) and ca == cb:
        diffs.append("一方有系数一方没有")
    major = [d for d in diffs if not d.startswith("仅小数位数不同")]
    verdict = "一致" if not major else "不一致"
    return verdict, "；".join(diffs)


rows = []
for _, p in papers.iterrows():
    a = mine.get(p.PMID, {"status": "", "genes": "", "coefs": "", "source": "", "note": ""})
    b = rev.get(p.PMID, {})
    v, d = compare(a, b)
    who, why = ADJ.get(p.PMID, ["", ""])
    if p.set == "AI预提取81":
        if who == "子代理":
            fz = {"status": b.get("status", ""), "genes": b.get("genes", ""), "coefs": b.get("coefs", "")}
        else:
            fz = {"status": "已提取", "genes": a["genes"], "coefs": a["coefs"]}
    else:
        fz = final.get(p.PMID, a)
    if p.PMID in ("38924839", "35983076", "39511540", "39399052"):
        fz = {"status": mine[p.PMID]["status"], "genes": "", "coefs": ""}
    if v == "不一致" and not who:
        who, why = "未裁定", ""
    rows.append({
        "PMID": p.PMID, "PMCID": p.PMCID, "集合": p.set, "标题": p.title,
        "提取者1（Claude）状态": a["status"], "提取者1 基因": a["genes"], "提取者1 系数": a["coefs"],
        "提取者1 来源": a["source"], "提取者1 备注": a["note"],
        "子代理 状态": b.get("status", ""), "子代理 基因": b.get("genes", ""), "子代理 系数": b.get("coefs", ""),
        "子代理 来源": b.get("source", ""), "子代理 终点": b.get("endpoint", ""), "子代理 训练队列": b.get("training", ""),
        "子代理 验证队列": b.get("validation", ""), "子代理 备注": b.get("note", ""),
        "两者一致性": v, "差异说明": d,
        "Claude 对照原文的裁定（谁对）": who, "裁定依据": why,
        "复核后状态": fz["status"], "复核后基因": fz["genes"], "复核后系数": fz["coefs"],
        "专家A 状态": "", "专家A 基因": "", "专家A 系数": "", "专家A 意见": "",
        "专家B 状态": "", "专家B 基因": "", "专家B 系数": "", "专家B 意见": "",
        "最终裁定": "", "裁定人/日期": "",
    })
df = pd.DataFrame(rows)
df["_o"] = df["两者一致性"].map({"不一致": 0, "子代理未完成": 1, "一致": 2})
df = df.sort_values(["_o", "集合", "PMID"]).drop(columns="_o")
summary = df.groupby(["集合", "两者一致性"]).size().unstack(fill_value=0)

with pd.ExcelWriter(OUT, engine="openpyxl") as w:
    pd.DataFrame({"说明": [
        "四方对比表：提取者1 = Claude 主提取（AI 预提取 81 篇 + 精读 206 篇）；子代理 = 只看原文与补充材料、不看提取者1结果的盲法复核。",
        "专家A / 专家B 列留空，请两位专家独立填写（建议不要先看前两列）；“最终裁定”由 Wei 或第三人填写。",
        "“两者一致性”只比较提取者1与子代理：判定类别（有系数/仅基因/非线性/排除…）、基因集合（忽略大小写和连字符）、每个基因的系数。",
        "系数只差小数位数（如 0.4025 vs 0.402463）记为“仅小数位数不同”，不算不一致。",
        "已知的预期差异：35414328 的系数由提取者1用作者补充数据重拟合，子代理按规则会判为“仅基因”。",
        "排序：不一致的排在最前，建议专家先看这些。",
    ]}).to_excel(w, sheet_name="说明", index=False)
    df.to_excel(w, sheet_name="对比", index=False)
    summary.reset_index().to_excel(w, sheet_name="汇总", index=False)

from openpyxl import load_workbook  # noqa: E402
from openpyxl.styles import Alignment, Font, PatternFill  # noqa: E402

wb = load_workbook(OUT)
for ws in wb.worksheets:
    for c in ws[1]:
        c.font = Font(bold=True)
        c.alignment = Alignment(wrap_text=True, vertical="top")
        c.fill = PatternFill("solid", fgColor="FFF2CC" if str(c.value).startswith(("专家", "最终", "裁定")) else "DDEBF7")
        ws.column_dimensions[c.column_letter].width = 120 if ws.title == "说明" else (
            40 if any(k in str(c.value) for k in ("基因", "系数", "备注", "差异", "标题", "来源", "意见")) else 14)
    ws.freeze_panes = "B2" if ws.title == "对比" else None
    for row in ws.iter_rows(min_row=2):
        for c in row:
            if isinstance(c.value, str) and len(c.value) > 30:
                c.alignment = Alignment(wrap_text=True, vertical="top")
    if ws.title == "对比":
        ws.auto_filter.ref = ws.dimensions
wb.save(OUT)
print(summary.to_string())
