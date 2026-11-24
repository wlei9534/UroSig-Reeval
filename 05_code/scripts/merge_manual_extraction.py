"""Merge the manual full-text extraction into the AI pre-extraction workbook.

Usage: python merge_manual_extraction.py <project_root> <input_workbook> <output_workbook>
Inputs : <input_workbook> (ai_prefill_PRAD.xlsx as last saved on the device),
         03_extraction/manual_extraction.jsonl, 03_extraction/manual_extraction_verify.csv,
         03_extraction/fulltext_xml/<PMCID>.xml (for the cohort-mention scan)
Output : workbook with the original sheets kept and these added/updated:
         逐篇            + 当前状态、人工精读字段、建议用途、原文提及的再评估队列
         人工精读提取    one row per paper (151)
         人工精读长表    one row per gene/coefficient with the literal-check result
         待人工核对      every item a human must check before the values are used
         汇总            recounted
"""
from __future__ import annotations

import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

ROOT, INP, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
X = ROOT / "03_extraction"

sheets = pd.read_excel(INP, sheet_name=None, dtype=str)
per = sheets["逐篇"].copy()
recs = pd.DataFrame([json.loads(l) for l in open(X / "manual_extraction.jsonl", encoding="utf-8")]).fillna("")
ver = pd.read_csv(X / "manual_extraction_verify.csv", dtype=str).fillna("")

USE = [  # (status prefix, suggested use)
    ("已提取（基因对 signature）", "主分析候选（基因对，系数来自正文/补充表）"),
    ("已提取（仅方向", "敏感性分析（仅方向，等权 ±1）"),
    ("已提取（从图片读取）", "敏感性分析（系数从图片读取，需对 PDF 原图核对）"),
    ("已提取（由 HR 换算）", "敏感性分析（系数由 HR 取对数换算）"),
    ("已提取", "主分析候选（系数来自正文/表格/补充材料）"),
    ("仅基因", "无系数：可做等权评分敏感性分析，或向作者索取系数"),
    ("非线性模型（附 LASSO", "敏感性分析（LASSO 阶段系数，非原模型）"),
    ("非线性模型", "排除：非线性模型无法按加权和复现"),
    ("系数仅在低分辨率图中", "待定：需看出版社高清 PDF"),
    ("非", "排除（理由见状态）"),
]


def use_of(status: str) -> str:
    for k, v in USE:
        if status.startswith(k):
            return v
    return ""


# Cohorts this project re-evaluates in; a mention is not proof of use, so this flags papers to check for leakage.
COH = {
    "GSE54460": r"GSE54460",
    "GSE70769": r"GSE70769|GSE70770",
    "GSE70768": r"GSE70768|GSE70770",
    "MSKCC": r"GSE21032|GSE21034|GSE21035|MSKCC|Memorial Sloan|prad_mskcc|Taylor et al",
    "GSE116918": r"GSE116918",
    "TCGA": r"TCGA",
}


def cohorts(pmcid: str) -> str:
    f = X / "fulltext_xml" / f"{pmcid}.xml"
    if not pmcid or not f.exists():
        return ""
    root = ET.fromstring(f.read_text(encoding="utf-8"))
    for parent in list(root.iter()):  # drop the reference list; author affiliations live in <front>, which is skipped
        for ch in list(parent):
            if ch.tag == "ref-list":
                parent.remove(ch)
    t = " ".join(" ".join(el.itertext()) for tag in ("body", "back") for el in root.iter(tag))
    return ";".join(k for k, rx in COH.items() if re.search(rx, t))


EP = {"BCR": r"biochemical(?:al)? recurrence|biochemical relapse|\\bBCR", "RFS": r"recurrence[- ]free survival|\\bRFS\\b",
      "DFS": r"disease[- ]free survival|\\bDFS\\b", "PFS/PFI": r"progression[- ]free (?:survival|interval)|\\bPFS\\b|\\bPFI\\b",
      "OS": r"overall survival|\\bOS\\b", "MFS": r"metastasis[- ]free survival|\\bMFS\\b"}


def endpoint_counts(pmcid: str) -> str:
    """Keyword counts in the body text: a pointer to the endpoint actually modelled, not a classification."""
    f = X / "fulltext_xml" / f"{pmcid}.xml"
    if not pmcid or not f.exists():
        return ""
    body = ET.fromstring(f.read_text(encoding="utf-8")).find(".//body")
    t = " ".join(body.itertext()) if body is not None else ""
    n = {k: len(re.findall(rx, t, flags=re.I if k != "OS" else 0)) for k, rx in EP.items()}
    return "；".join(f"{k} {v}" for k, v in sorted(n.items(), key=lambda kv: -kv[1]) if v)


recs["建议用途"] = recs["status"].map(use_of)
recs["终点关键词次数（正文检索）"] = recs["pmcid"].map(endpoint_counts)
recs["基因数"] = recs["genes"].map(lambda g: len(g.split(";")) if g else 0)
recs["原文提及的再评估队列（正文+附录检索，不含参考文献）"] = recs["pmcid"].map(cohorts)

# per-paper literal-check summary from the verify table
def lit(pmid: str) -> str:
    v = ver[ver["PMID"] == pmid]
    if v.empty:
        return ""
    if v["数值来源"].str.startswith("图片").all():
        return "数值来自图片/HR 换算，原文无此数；基因名在原文 %d/%d" % ((v["基因在原文"] == "True").sum(), len(v))
    nc = v[v["系数在原文"] != ""]
    return "系数在原文 %d/%d；基因名在原文 %d/%d；基因与系数相邻 %d/%d" % (
        (nc["系数在原文"] == "True").sum(), len(nc), (v["基因在原文"] == "True").sum(), len(v),
        (nc["基因与系数相邻(≤200字符)"] == "True").sum(), len(nc))


recs["逐字核对（人工精读）"] = recs["pmid"].map(lit)

# explained exceptions from the literal check (all checked by hand on 2026-10-31)
EXPLAINED = {
    ("34136527", ""): "蛋白名→基因符号为提取者映射（alpha-Catenin→CTNNA1、DJ1→PARK7、YB1→YBX1），数值已在原文找到",
    ("38240702", "SDC1"): "原文写作 SDC-1，公式逐字核对一致",
    ("41346575", "RRM2"): "公式为 MathML，紧接下一段落，导致脚本相邻判定失败；人工核对一致",
    ("40355680", "CYP4F8"): "原文误写 CYP4FB，提取时改为 CYP4F8",
    ("41495164", ""): "基因名只出现在 Figure 5B 图片中",
    ("42138194", ""): "基因与系数只出现在 Figure 3C 图片表格中",
    ("34900992", "LCMT1-AS2"): "正文误写为 LCMT1-SA2；Figure 1C 高清图为 LCMT1-AS2",
    ("37008945", "ACPP"): "公式里笔误写作 APCC，同段基因列表为 ACPP（现行符号 ACP3）",
}


def explain(pmid: str, gene: str) -> str:
    return EXPLAINED.get((pmid, gene)) or EXPLAINED.get((pmid, ""), "")


ver["异常说明"] = [explain(p, g) if (gf != "True" or (cf == "True" and nr != "True")) else ""
                 for p, g, gf, cf, nr in zip(ver["PMID"], ver["基因"], ver["基因在原文"], ver["系数在原文"],
                                              ver["基因与系数相邻(≤200字符)"])]

# ---- 逐篇: add merged status and manual fields
m = recs.set_index("pmid")
per["当前状态（2026-10-31）"] = per["PMID"].map(m["status"]).fillna("")
per.loc[per["当前状态（2026-10-31）"] == "", "当前状态（2026-10-31）"] = per["状态"].map(
    lambda s: "需人工获取全文（见 无全文文献清单_60篇.txt）" if s == "需人工获取全文" else s)
for col, src in [("人工精读基因数", "基因数"), ("人工精读基因", "genes"), ("人工精读系数", "coefs"),
                 ("人工精读来源", "source"), ("人工精读备注", "note"), ("逐字核对（人工精读）", "逐字核对（人工精读）"),
                 ("建议用途", "建议用途")]:
    per[col] = per["PMID"].map(m[src]).fillna("")
per.loc[per["当前状态（2026-10-31）"].str.startswith("已预提取"), "建议用途"] = "主分析候选（AI 预提取，待人工双核）"
per.loc[per["状态"] == "需人工获取全文", "建议用途"] = "需先获取全文（见 无全文文献清单_60篇.txt）"
per["原文提及的再评估队列（正文+附录检索，不含参考文献）"] = per["PMCID"].fillna("").map(cohorts)
per["终点关键词次数（正文检索）"] = per["PMCID"].fillna("").map(endpoint_counts)
cols = list(per.columns)
cols.insert(cols.index("状态") + 1, cols.pop(cols.index("当前状态（2026-10-31）")))
per = per[cols]

# ---- 人工精读提取
meta = sheets["逐篇"].set_index("PMID")[["年份", "期刊", "标题"]]
pr = recs.join(meta, on="pmid")
pr = pr.rename(columns={"pmid": "PMID", "pmcid": "PMCID", "status": "人工精读状态", "genes": "基因", "coefs": "系数",
                        "source": "来源（原文位置）", "note": "备注"})
pr = pr[["PMID", "PMCID", "年份", "期刊", "标题", "人工精读状态", "建议用途", "基因数", "基因", "系数", "来源（原文位置）",
         "备注", "终点关键词次数（正文检索）", "原文提及的再评估队列（正文+附录检索，不含参考文献）", "逐字核对（人工精读）"]]

# ---- 待人工核对
FLAG = re.compile(r"需核对|需对 ?PDF|不确定|末位|疑|映射|旧符号|现名|有误|矛盾|不一致|歧义|误写|目测|重复发表|拼作|不稳|极大|EPV|事件数仅|多一个")
todo = []
for _, r in recs.iterrows():
    why = []
    if "从图片读取" in r["status"]:
        why.append("系数从图片读取：对照出版社 PDF 原图核对每个数字")
    if "由 HR 换算" in r["status"]:
        why.append("系数由 HR 取对数换算：只用于敏感性分析，核对 HR 读数")
    if "系数仅在低分辨率图中" in r["status"]:
        why.append("PMC 图片太模糊：需看出版社 PDF 高清图")
    if FLAG.search(r["note"]):
        why.append("备注中的待核事项（见右栏）")
    if why:
        todo.append({"PMID": r["pmid"], "PMCID": r["pmcid"], "人工精读状态": r["status"], "核对什么": "；".join(why),
                     "备注": r["note"], "来源（原文位置）": r["source"]})
todo = pd.DataFrame(todo)

# ---- 人工精读长表 (+ current HGNC symbol, if an HGNC snapshot path is given as the 4th argument)
lt = ver.rename(columns={"状态": "人工精读状态"})
if len(sys.argv) > 4:
    from collections import defaultdict
    sys.path.insert(0, str(ROOT / "05_code"))
    from urosig.genes import build_symbol_map
    h = pd.read_csv(sys.argv[4], sep="\t", dtype=str, low_memory=False)
    smap = build_symbol_map(h)
    cand = defaultdict(set)
    for col in ("prev_symbol", "alias_symbol"):
        for sym, fld in zip(h["symbol"], h[col].fillna("")):
            for o in filter(None, (x.strip() for x in fld.split("|"))):
                cand[o].add(sym)

    def hgnc(sym: str) -> str:
        out = []
        for g in sym.split("|"):
            t = smap.get(g)
            if t is None:
                out.append("未映射")
                continue
            others = sorted(cand.get(g, set()) - {t}) if t != g else []
            out.append(t + (f"（亦为 {'/'.join(others)} 的旧名或别名）" if others else ""))
        return "|".join(out)

    lt.insert(lt.columns.get_loc("基因") + 1, "HGNC 现行符号（2026-10-26 快照）", lt["基因"].map(hgnc))

# ---- 汇总
summ = per["当前状态（2026-10-31）"].value_counts().rename_axis("当前状态").reset_index(name="篇数")
use_summ = per["建议用途"].replace("", "（AI 预提取以外、无人工精读记录）").value_counts().rename_axis("建议用途").reset_index(name="篇数")

note_lines = sheets["说明"]["说明"].tolist() + [
    "",
    "2026-10-31 人工精读补充（Claude 逐篇阅读全文、补充材料和图片完成，仍需第二人独立核对）：",
    "对象：AI 预提取中“未找到系数，需人工”142 篇 +“系数可能在补充材料”9 篇，共 151 篇，全部完成，见“人工精读提取”页。",
    "状态含义：已提取 = 系数逐字取自正文/表格/补充材料；已提取（从图片读取）= 系数印在图上、按图读出；"
    "已提取（由 HR 换算）= 原文只给多因素 HR，系数 = ln(HR)；仅基因（无系数）= 原文未给系数；"
    "非线性模型 = RSF、XGBoost、SVM 等无法按加权和复现；其余“非…”为不符合纳入标准的理由。",
    "逐字核对：manual_extraction_verify.csv 及“人工精读长表”页逐项检查基因名与系数是否在原文/补充材料中出现；"
    "所有取自正文/表格/补充材料的系数均在原文找到，例外见“异常说明”栏（已人工核对）。",
    "“待人工核对”页列出图片读数、HR 换算、符号映射、疑似重复发表、原文矛盾等，用于双人核对时逐项确认。",
    "“终点关键词次数”为正文中 BCR/RFS/DFS/PFS/OS/MFS 等词的出现次数，只提示原文建模终点，需人工确认。",
    "“原文提及的再评估队列”为正文和附录（不含作者单位、参考文献）检索结果（提到不等于用于建模），用于判断泄漏：原文建模/验证用过的队列，"
    "不能再用来评估该 signature。",
    "“建议用途”只是提取者的建议，纳入与否按方案由双人提取和 Wei 决定。",
    "“人工精读长表”的 HGNC 列按 2026-10-26 快照映射；“未映射”多为 GENCODE 克隆名 lncRNA（如 AC139749.1）或 ENSG 编号，需按 Ensembl 版本另行对应。",
    "无开放获取全文的 60 篇见 03_extraction/无全文文献清单_60篇.txt；找到全文后放入 03_extraction/fulltext_manual/<PMID>/。",
]

with pd.ExcelWriter(OUT, engine="openpyxl") as w:
    pd.DataFrame({"说明": note_lines}).to_excel(w, sheet_name="说明", index=False)
    per.to_excel(w, sheet_name="逐篇", index=False)
    pr.to_excel(w, sheet_name="人工精读提取", index=False)
    lt.to_excel(w, sheet_name="人工精读长表", index=False)
    todo.to_excel(w, sheet_name="待人工核对", index=False)
    sheets["基因-系数长表"].to_excel(w, sheet_name="基因-系数长表（AI预提取）", index=False)
    summ.to_excel(w, sheet_name="汇总", index=False, startrow=0)
    use_summ.to_excel(w, sheet_name="汇总", index=False, startrow=len(summ) + 3)

# formatting: header style, freeze panes, column widths, wrap long text
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill

wb = load_workbook(OUT)
WIDE = {"标题": 50, "基因": 40, "系数": 40, "人工精读基因": 40, "人工精读系数": 40, "来源（原文位置）": 50,
        "人工精读来源": 50, "备注": 60, "人工精读备注": 60, "核对什么": 40, "说明": 120, "公式原文（截取）": 40}
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
print(summ.to_string(index=False))
print(use_summ.to_string(index=False))
print("待人工核对", len(todo))
