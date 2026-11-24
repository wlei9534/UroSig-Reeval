"""Fetch PMC full text for candidate papers and pre-extract genes/coefficients (third-reader aid).

Usage: python prefill_signatures.py <project_root> <candidate_csv>
Outputs: 03_extraction/fulltext_xml/<PMCID>.xml and 03_extraction/ai_prefill_PRAD.xlsx
Every extracted value carries its source (paragraph text or table label) for checking against the paper.
"""
from __future__ import annotations

import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

ROOT = Path(sys.argv[1]); sys.path.insert(0, str(ROOT / "05_code"))
from urosig.fulltext import _SCORE_CTX, cited_formula, coef_tables, formula_terms, supplement_hint  # noqa: E402

CAND = pd.read_csv(sys.argv[2], dtype=str)
XD = ROOT / "03_extraction" / "fulltext_xml"; XD.mkdir(parents=True, exist_ok=True)
EF = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fetch(pmcids):
    """curl is used because urllib is throttled on this machine."""
    import subprocess
    url = f"{EF}?db=pmc&retmode=xml&id=" + ",".join(p.replace("PMC", "") for p in pmcids)
    r = subprocess.run(["curl.exe", "-s", "--retry", "3", "--max-time", "120", url], capture_output=True)
    return r.stdout.decode("utf-8", errors="replace")


todo = [p for p in CAND["pmcid"].dropna() if not (XD / f"{p}.xml").exists()]
for i in range(0, len(todo), 10):
    batch = todo[i:i + 10]
    xml = fetch(batch)
    for art in re.findall(r"<article[\s>].*?</article>", xml, flags=re.S):
        m = re.search(r'<article-id pub-id-type="pmc(?:id)?">(?:PMC)?(\d+)</article-id>', art)
        if m:
            (XD / f"PMC{m.group(1)}.xml").write_text(art, encoding="utf-8")
    print(f"fetched {min(i + 10, len(todo))}/{len(todo)}"); time.sleep(0.4)

rows, terms_long = [], []
for _, r in CAND.iterrows():
    rec = {"PMID": r["pmid"], "PMCID": r.get("pmcid") or "", "年份": r["year"], "期刊": r["journal"], "标题": r["title"]}
    f = XD / f"{r['pmcid']}.xml" if isinstance(r.get("pmcid"), str) else None
    if not f or not f.exists():
        rec.update({"全文": "无 PMC 全文", "状态": "需人工获取全文"}); rows.append(rec); continue
    jats = f.read_text(encoding="utf-8")
    try:
        root = ET.fromstring(jats)
    except ET.ParseError:
        rec.update({"全文": "XML 解析失败", "状态": "需人工"}); rows.append(rec); continue
    xml_pmid = [e.text for e in root.iter("article-id") if e.get("pub-id-type") == "pmid"]
    if xml_pmid and xml_pmid[0] != str(r["pmid"]):
        rec.update({"全文": f"PMC 全文 PMID 不符（{xml_pmid[0]}）", "状态": "需人工"}); rows.append(rec); continue
    body = root.find(".//body")
    if body is None or len("".join(body.itertext())) < 2000:
        rec.update({"全文": "PMC 仅元数据（非开放获取正文）", "状态": "需人工获取全文"}); rows.append(rec); continue
    paras = ["".join(p.itertext()) for p in body.iter("p")]
    best, src, n_cited = [], "", 0
    for p in paras:
        t = formula_terms(p)
        if not t:
            continue
        if cited_formula(p):
            n_cited += 1
            continue
        if len(t) > len(best):
            m = _SCORE_CTX.search(p)
            best, src = t, p[max(0, m.start() - 40): m.start() + 400] if m else p[:400]
    tabs = coef_tables(jats)
    tab = max(tabs, key=lambda x: len(x["terms"])) if tabs else None
    supp = supplement_hint(" ".join(paras))
    supp_files = [s.get("{http://www.w3.org/1999/xlink}href", "") for s in root.iter("supplementary-material")]
    supp_files += [m.get("{http://www.w3.org/1999/xlink}href", "") for s in root.iter("supplementary-material") for m in s.iter("media")]
    choice, terms = None, []
    if best and tab:
        same = {g for g, _ in best} == {g for g, _ in tab["terms"]}
        choice, terms = ("正文公式（与表格一致）" if same else "正文公式（与表格不一致，需核对）"), best
    elif best:
        choice, terms = "正文公式", best
    elif tab:
        choice, terms = f"{tab['label']}：{tab['caption'][:60]}", tab["terms"]
    notes = []
    if len({g for g, _ in terms}) < len(terms):
        notes.append("基因重复出现，可能是基因对（gene pair）signature，需按原文重建")
    if n_cited:
        notes.append(f"另有 {n_cited} 处公式引自其他研究，已跳过")
    if not best and tab and tab.get("mixed_uni_multi"):
        notes.append("系数来自含单因素结果的表格，需确认取的是多因素模型系数")
    note = "；".join(notes)
    status = "已预提取，待人工核对" if terms else ("系数可能在补充材料" if supp else "未找到系数，需人工")
    rec.update({"全文": "PMC 开放获取正文", "来源": choice or "", "基因数": len(terms) or "",
                "基因": ";".join(g for g, _ in terms), "系数": ";".join(f"{c:g}" for _, c in terms),
                "公式原文（截取）": src if best else "", "补充材料提示": "是" if supp else "",
                "补充材料文件": ";".join(dict.fromkeys(x for x in supp_files if x))[:300], "状态": status, "提示": note})
    rows.append(rec)
    for k, (g, c) in enumerate(terms, 1):
        terms_long.append({"PMID": r["pmid"], "PMCID": r["pmcid"], "序号": k, "基因": g, "系数": c, "来源": choice})

df = pd.DataFrame(rows)
cols = ["PMID", "PMCID", "年份", "期刊", "标题", "全文", "状态", "提示", "来源", "基因数", "基因", "系数", "公式原文（截取）", "补充材料提示", "补充材料文件"]
df = df.reindex(columns=cols)
out = ROOT / "03_extraction" / "ai_prefill_PRAD.xlsx"
with pd.ExcelWriter(out, engine="openpyxl") as w:
    pd.DataFrame({"说明": [
        "AI 预提取（第三方核对用，不替代双人人工提取）",
        "来源：PMC 开放获取全文 JATS XML，2026-10-26 抓取；候选 = 题录筛选表中 AI 预筛为“可能纳入”的 292 篇",
        "规则：先找正文中的风险评分公式，再找表头含 coefficient/β 的表格；两者都有时比对基因集合是否一致",
        "每一行的基因和系数都必须对照原文（公式原文 / 表格编号）核对后才能进入提取表；自动提取可能漏项、错位或误读符号",
        "未覆盖：补充材料里的系数表、图片形式的公式、非开放获取文章——这些标为“需人工”",
    ]}).to_excel(w, sheet_name="说明", index=False)
    df.to_excel(w, sheet_name="逐篇", index=False)
    pd.DataFrame(terms_long).to_excel(w, sheet_name="基因-系数长表", index=False)
    pd.DataFrame(df["状态"].value_counts()).reset_index().to_excel(w, sheet_name="汇总", index=False)
print(df["状态"].value_counts().to_string().encode("ascii", "backslashreplace").decode()); print("saved", out)
