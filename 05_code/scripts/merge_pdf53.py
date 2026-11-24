"""Merge the extraction from Wei's manually downloaded PDFs into the workbook.

Usage: python merge_pdf53.py <project_root> <input_workbook> <output_workbook> <pdf_text_dir> <hgnc_txt>
Inputs : <input_workbook>   ai_prefill_PRAD.xlsx as produced by merge_manual_extraction.py
         03_extraction/manual_extraction_pdf53.jsonl, 03_extraction/manual_extraction_pdf53_verify.csv
         <pdf_text_dir>/<PMID>.raw.txt (pdftotext of each supplied PDF, used for the cohort/endpoint scans)
Output : same sheets; the 60 rows that were "需人工获取全文" now carry the PDF extraction (or the reason
         there is none), and the per-paper / long / to-check sheets gain the new rows.
"""
from __future__ import annotations

import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT, INP, OUT, TXT, HGNC = (Path(a) for a in sys.argv[1:6])
X = ROOT / "03_extraction"
STAMP = "当前状态（2026-10-31）"

sheets = pd.read_excel(INP, sheet_name=None, dtype=str)
sheets = {k: v.fillna("") for k, v in sheets.items()}
per = sheets["逐篇"]
recs = pd.DataFrame([json.loads(l) for l in open(X / "manual_extraction_pdf53.jsonl", encoding="utf-8")]).fillna("")
ver = pd.read_csv(X / "manual_extraction_pdf53_verify.csv", dtype=str).fillna("")

# papers in the 60-list with no extraction, and why (Wei, 2026-10-31)
NOT_SUPPLIED = {
    "41568495": "未获取全文：Wei 未能找到全文，按其决定放弃",
    "41111380": "未获取全文：Wei 未能找到全文，按其决定放弃",
    "40916526": "未获取全文：中文文献，按 Wei 决定排除",
    "37846114": "未获取全文：中文文献，按 Wei 决定排除",
    "33504416": "未获取全文：中文文献，按 Wei 决定排除",
}

USE = [
    ("已提取（基因对 signature）", "主分析候选（基因对，系数来自正文/补充表）"),
    ("已提取（从图片读取）", "敏感性分析（系数从图片读取，需对 PDF 原图核对）"),
    ("已提取（由作者补充数据重拟合）", "主分析候选（系数由作者公布的训练数据重拟合，阈值与原文一致；待双人确认）"),
    ("已提取", "主分析候选（系数来自正文/表格/补充材料）"),
    ("仅基因", "无系数：可做等权评分敏感性分析，或向作者索取系数"),
    ("系数在补充材料", "待定：系数在出版社补充材料，需 Wei 下载"),
    ("非线性模型", "排除：非线性模型无法按加权和复现"),
    ("非", "排除（理由见状态）"),
    ("未获取全文", "无全文，未提取"),
]


def use_of(status: str) -> str:
    return next((v for k, v in USE if status.startswith(k)), "")


COH = {
    "GSE54460": r"GSE54460",
    "GSE70769": r"GSE70769|GSE70770",
    "GSE70768": r"GSE70768|GSE70770",
    "MSKCC": r"GSE21032|GSE21034|GSE21035|MSKCC|Memorial Sloan|prad_mskcc|Taylor(?! ?& ?Francis)",
    "GSE116918": r"GSE116918",
    "TCGA": r"TCGA",
}
EP = {"BCR": r"biochemical(?:al)? recurrence|biochemical relapse|\bBCR", "RFS": r"recurrence[- ]free survival|\bRFS\b",
      "DFS": r"disease[- ]free survival|\bDFS\b", "PFS/PFI": r"progression[- ]free (?:survival|interval)|\bPFS\b|\bPFI\b",
      "OS": r"overall survival|\bOS\b", "MFS": r"metastasis[- ]free survival|\bMFS\b"}


def body_text(pmid: str) -> str:
    """PDF text with the reference list cut off (from the last 'References' heading onward)."""
    f = TXT / f"{pmid}.raw.txt"
    if not f.exists():
        return ""
    t = f.read_text(encoding="utf-8", errors="replace")
    t = re.sub(r"-\n(?=[a-z])", "", t)
    m = list(re.finditer(r"\n\s*(?:References|REFERENCES|Reference|Bibliography)\s*\n", t))
    return t[: m[-1].start()] if m else t


def cohorts(pmid: str) -> str:
    t = body_text(pmid)
    return ";".join(k for k, rx in COH.items() if re.search(rx, t)) if t else ""


def endpoint_counts(pmid: str) -> str:
    t = body_text(pmid)
    n = {k: len(re.findall(rx, t, flags=re.I if k != "OS" else 0)) for k, rx in EP.items()}
    return "；".join(f"{k} {v}" for k, v in sorted(n.items(), key=lambda kv: -kv[1]) if v)


def lit(pmid: str) -> str:
    v = ver[ver["PMID"] == pmid]
    if v.empty:
        return ""
    if v["数值来源"].str.startswith("图片").all():
        return "数值来自图片或重拟合，原文无此数；基因名在原文 %d/%d" % ((v["基因在原文"] == "True").sum(), len(v))
    nc = v[v["系数在原文"] != ""]
    if nc.empty:
        return "基因名在原文 %d/%d（无系数）" % ((v["基因在原文"] == "True").sum(), len(v))
    return "系数在原文 %d/%d；基因名在原文 %d/%d；基因与系数相邻 %d/%d" % (
        (nc["系数在原文"] == "True").sum(), len(nc), (v["基因在原文"] == "True").sum(), len(v),
        (nc["基因与系数相邻(≤200字符)"] == "True").sum(), len(nc))


EXPLAINED = {  # literal-check exceptions, each checked by hand on 2026-10-31
    ("39760809", "MATN3"): "原文排版为 “MATN * 30.449”，判读为 MATN3 * 0.449（判读，需双人确认）",
    ("40295064", ""): "原文为科学计数法（‑8.4×10–1 等），脚本无法按数值匹配；人工核对一致",
    ("37193176", "BCL2L11"): "公式跨页（p2463→下一栏），中间夹页眉与图注，导致相邻判定失败；人工核对一致",
}


def explain(pmid: str, gene: str) -> str:
    return EXPLAINED.get((pmid, gene)) or EXPLAINED.get((pmid, ""), "")


ver["异常说明"] = [explain(p, g) if (gf != "True" or cf == "False" or (cf == "True" and nr != "True")) else ""
                 for p, g, gf, cf, nr in zip(ver["PMID"], ver["基因"], ver["基因在原文"], ver["系数在原文"],
                                              ver["基因与系数相邻(≤200字符)"])]

recs["建议用途"] = recs["status"].map(use_of)
recs["基因数"] = recs["genes"].map(lambda g: len(g.split(";")) if g else 0)
recs["终点关键词次数（正文检索）"] = recs["pmid"].map(endpoint_counts)
recs["原文提及的再评估队列（正文+附录检索，不含参考文献）"] = recs["pmid"].map(cohorts)
recs["逐字核对（人工精读）"] = recs["pmid"].map(lit)

# ---- 逐篇
m = recs.set_index("pmid")
target = per["状态"] == "需人工获取全文"
assert target.sum() == 60, target.sum()
assert set(m.index) | set(NOT_SUPPLIED) == set(per.loc[target, "PMID"]), "PMID sets differ"
for i in per.index[target]:
    p = per.at[i, "PMID"]
    if p in m.index:
        r = m.loc[p]
        per.at[i, STAMP] = r["status"]
        for col, src in [("人工精读基因数", "基因数"), ("人工精读基因", "genes"), ("人工精读系数", "coefs"),
                         ("人工精读来源", "source"), ("人工精读备注", "note"), ("逐字核对（人工精读）", "逐字核对（人工精读）"),
                         ("建议用途", "建议用途"), ("终点关键词次数（正文检索）", "终点关键词次数（正文检索）"),
                         ("原文提及的再评估队列（正文+附录检索，不含参考文献）", "原文提及的再评估队列（正文+附录检索，不含参考文献）")]:
            per.at[i, col] = str(r[src])
        per.at[i, "全文"] = "Wei 手动下载 PDF（fulltext/shouldbe60files）"
    else:
        per.at[i, STAMP] = NOT_SUPPLIED[p]
        per.at[i, "建议用途"] = use_of(NOT_SUPPLIED[p])

# ---- 人工精读提取
pr = sheets["人工精读提取"]
if "全文来源" not in pr.columns:
    pr.insert(pr.columns.get_loc("PMCID") + 1, "全文来源", "PMC 开放获取包")
meta = per.set_index("PMID")[["PMCID", "年份", "期刊", "标题"]]
new = recs.join(meta, on="pmid", rsuffix="_m").rename(columns={
    "pmid": "PMID", "status": "人工精读状态", "genes": "基因", "coefs": "系数", "source": "来源（原文位置）", "note": "备注"})
new["PMCID"] = new["PMCID"]
new["全文来源"] = new["fulltext"]
new["基因数"] = new["基因数"].astype(str)
pr = pd.concat([pr, new[pr.columns]], ignore_index=True)

# ---- 人工精读长表 (+ HGNC current symbol)
sys.path.insert(0, str(ROOT / "05_code"))
from urosig.genes import build_symbol_map  # noqa: E402

h = pd.read_csv(HGNC, sep="\t", dtype=str, low_memory=False)
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


lt = sheets["人工精读长表"]
nv = ver.rename(columns={"状态": "人工精读状态"})
nv["PMCID"] = nv["PMID"].map(meta["PMCID"]).fillna("")
nv["HGNC 现行符号（2026-10-26 快照）"] = nv["基因"].map(hgnc)
# papers whose Ensembl IDs are known (listed in the note): map by ID, which beats symbol lookup for clone-named lncRNAs
ens2sym = h.set_index("ensembl_gene_id")["symbol"].to_dict()
for _, r in recs.iterrows():
    mm = re.search(r"Ensembl ID（与基因同序）：([A-Z0-9;]+)", r["note"])
    if not mm:
        continue
    for g, e in zip(r["genes"].split(";"), mm.group(1).split(";")):
        k = (nv["PMID"] == r["pmid"]) & (nv["基因"] == g)
        nv.loc[k, "HGNC 现行符号（2026-10-26 快照）"] = f"{ens2sym[e]}（按 {e}）" if e in ens2sym else f"未映射（{e} 不在 HGNC）"
        nv.loc[k, "异常说明"] = f"按 Ensembl ID {e} 匹配表达（作者补充数据给出）"
# lncRNA symbols that HGNC now lists as previous symbols of protein-coding loci (checked against the snapshot):
# the formal mapping holds, but the transcript the paper measured may not be the coding one.
LNC_TO_PC = {"UBXN10-AS1": "PLA2G2C", "LINC01125": "C2orf92"}
for g, t in LNC_TO_PC.items():
    k = (nv["基因"] == g) & ~nv["异常说明"].str.startswith("按 Ensembl ID")
    nv.loc[k, "异常说明"] = (nv.loc[k, "异常说明"] + "；" if nv.loc[k, "异常说明"].ne("").any() else "") + \
        f"HGNC 快照中 {g} 为 {t}（蛋白编码基因）的旧符号；原文测的是 lncRNA，表达匹配需按 Ensembl ID 核对"
nv["数值来源"] = nv["数值来源"].str.replace("图片/重拟合（原文无此数）", "图片/HR换算/重拟合（原文无此数）", regex=False)
if "全文来源" not in lt.columns:
    lt.insert(lt.columns.get_loc("PMCID") + 1, "全文来源", "PMC 开放获取包")
nv["全文来源"] = "Wei 手动下载 PDF"
for c in lt.columns:
    if c not in nv.columns:
        nv[c] = ""
lt = pd.concat([lt, nv[lt.columns].astype(str)], ignore_index=True)

# ---- 待人工核对
FLAG = re.compile(r"需双人确认|需核对|需对 ?PDF|不确定|疑|映射|旧符号|有误|矛盾|不一致|存疑|歧义|误写|判读|排版|异常|谨慎|表述不一|未核验")
todo_new = []
for _, r in recs.iterrows():
    why = []
    if "从图片读取" in r["status"]:
        why.append("系数从图片读取：对照 PDF 原图核对每个数字")
    if "重拟合" in r["status"]:
        why.append("系数由作者补充数据重拟合：复跑 05_code/scripts/refit_35414328.py 核对")
    if "系数在补充材料" in r["status"] or "待下载" in r["status"]:
        why.append("系数/基因列表在出版社补充材料：需 Wei 下载后补提")
    if FLAG.search(r["note"]):
        why.append("备注中的待核事项（见右栏）")
    if why:
        todo_new.append({"PMID": r["pmid"], "PMCID": meta["PMCID"].get(r["pmid"], ""), "人工精读状态": r["status"],
                         "核对什么": "；".join(why), "备注": r["note"], "来源（原文位置）": r["source"]})
for p, why in NOT_SUPPLIED.items():
    if "待 Wei 确认" in why or "PDF 损坏" in why:
        todo_new.append({"PMID": p, "PMCID": meta["PMCID"].get(p, ""), "人工精读状态": why, "核对什么": why,
                         "备注": "", "来源（原文位置）": ""})
todo = pd.concat([sheets["待人工核对"], pd.DataFrame(todo_new)], ignore_index=True).fillna("")

# ---- 汇总
summ = per[STAMP].value_counts().rename_axis("当前状态").reset_index(name="篇数")
use_summ = per["建议用途"].replace("", "（AI 预提取以外、无人工精读记录）").value_counts().rename_axis("建议用途").reset_index(name="篇数")

notes = sheets["说明"]["说明"].tolist()
notes = [n for n in notes if not n.startswith("无开放获取全文的 60 篇见")] + [
    "",
    "2026-10-31 无开放获取全文 60 篇的补充提取（Wei 手动下载 PDF，Claude 逐篇阅读）：",
    "Wei 提供 55 个 PDF（03_extraction/fulltext/shouldbe60files），全部读完；未提取 5 篇：41568495、41111380（找不到全文）、"
    "40916526、37846114、33504416（中文文献），均按 Wei 决定放弃。39252915（medRxiv 预印本）按 Wei 决定排除。",
    "这些行在“人工精读提取”“人工精读长表”中的“全文来源”为“Wei 手动下载 PDF”；逐字核对以 pdftotext 抽取的 PDF 文本及已下载的补充材料为语料"
    "（05_code/scripts/verify_pdf53.py → 03_extraction/manual_extraction_pdf53_verify.csv）。",
    "补充材料：33338276 Table S3（Wiley）、35414328 Supplementary file 1–6、38070051 Table S3（Springer）已下载并提取，在 03_extraction/fulltext_manual/supp_<PMID>/。",
    "35414328 原文未给系数：用作者公布的训练集表达数据重拟合多因素 Cox（05_code/scripts/refit_35414328.py），风险评分中位数与原文阈值 −2.5775 一致；表达匹配须按备注中的 Ensembl ID。",
]

with pd.ExcelWriter(OUT, engine="openpyxl") as w:
    pd.DataFrame({"说明": notes}).to_excel(w, sheet_name="说明", index=False)
    per.to_excel(w, sheet_name="逐篇", index=False)
    pr.to_excel(w, sheet_name="人工精读提取", index=False)
    lt.to_excel(w, sheet_name="人工精读长表", index=False)
    todo.to_excel(w, sheet_name="待人工核对", index=False)
    sheets["基因-系数长表（AI预提取）"].to_excel(w, sheet_name="基因-系数长表（AI预提取）", index=False)
    summ.to_excel(w, sheet_name="汇总", index=False, startrow=0)
    use_summ.to_excel(w, sheet_name="汇总", index=False, startrow=len(summ) + 3)

from openpyxl import load_workbook  # noqa: E402
from openpyxl.styles import Alignment, Font, PatternFill  # noqa: E402

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
print("人工精读提取", len(pr), "| 长表", len(lt), "| 待人工核对", len(todo))
print("HGNC 未映射:", nv.loc[nv["HGNC 现行符号（2026-10-26 快照）"].str.contains("未映射"), ["PMID", "基因"]].values.tolist())
print("HGNC 改名:", nv.loc[(nv["HGNC 现行符号（2026-10-26 快照）"] != nv["基因"]) & ~nv["HGNC 现行符号（2026-10-26 快照）"].str.contains("未映射"),
                          ["PMID", "基因", "HGNC 现行符号（2026-10-26 快照）"]].values.tolist())
