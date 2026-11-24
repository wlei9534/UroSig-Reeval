"""Build per-paper reading digests (Markdown) for manual extraction of signature genes/coefficients.

Usage: python build_digests.py <project_root>
Reads 03_extraction/ai_prefill_PRAD.xlsx (status column) and 03_extraction/fulltext/<PMCID>/ (XML + supplements).
Writes 03_extraction/digests/<PMID>.md
"""
from __future__ import annotations

import csv
import io
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

ROOT = Path(sys.argv[1])
FT = ROOT / "03_extraction" / "fulltext"
OUT = ROOT / "03_extraction" / "digests"; OUT.mkdir(parents=True, exist_ok=True)
WANT = {"未找到系数，需人工", "系数可能在补充材料"}
KEY = re.compile(r"coefficient|LASSO|least absolute|risk ?score|riskscore|multivariate cox|multivariable cox|stepwise|signature|prognostic model|risk model|β|beta|formula", re.I)
GENEISH = re.compile(r"^[A-Z][A-Z0-9]*(?:[-.][A-Z0-9]+)*$")
NUMISH = re.compile(r"^[+\-−–]?\s*\d+(?:\.\d+)?(?:[eE][-+]?\d+)?$")


def txt(el):
    return re.sub(r"\s+", " ", "".join(el.itertext())).strip() if el is not None else ""


def table_rows(tw):
    tab = tw.find(".//table")
    if tab is None:
        return []
    return [[txt(c) for c in tr if c.tag in ("td", "th")] for tr in tab.iter("tr")]


def geneish_rows(rows):
    return sum(1 for r in rows if r and any(GENEISH.match(c) for c in r[:2]))


def md_table(rows, max_rows=60):
    rows = [r for r in rows if any(r)][: max_rows + 1]
    if not rows:
        return ""
    w = max(len(r) for r in rows)
    rows = [r + [""] * (w - len(r)) for r in rows]
    s = ["| " + " | ".join(c.replace("|", "/")[:40] for c in rows[0]) + " |", "|" + "---|" * w]
    s += ["| " + " | ".join(c.replace("|", "/")[:40] for c in r) + " |" for r in rows[1:]]
    return "\n".join(s)


def supp_tables(path: Path):
    """Yield (label, rows) for tabular supplements that look like gene x number tables."""
    try:
        if path.suffix.lower() in (".xlsx", ".xlsm"):
            for name, df in pd.read_excel(path, sheet_name=None, header=None, dtype=str).items():
                yield f"{path.name} / {name}", df.fillna("").values.tolist()
        elif path.suffix.lower() == ".xls":
            for name, df in pd.read_excel(path, sheet_name=None, header=None, dtype=str).items():
                yield f"{path.name} / {name}", df.fillna("").values.tolist()
        elif path.suffix.lower() in (".csv", ".tsv", ".txt") and path.stat().st_size < 5e6 and not path.stem.startswith("PMC"):
            sep = "\t" if path.suffix.lower() == ".tsv" else ","
            yield path.name, list(csv.reader(io.StringIO(path.read_text(errors="replace")), delimiter=sep))
        elif path.suffix.lower() == ".docx":
            with zipfile.ZipFile(path) as z:
                x = ET.fromstring(z.read("word/document.xml"))
            ns = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
            for k, t in enumerate(x.iter(ns + "tbl"), 1):
                rows = [["".join(n.text or "" for n in tc.iter(ns + "t")).strip() for tc in tr.iter(ns + "tc")] for tr in t.iter(ns + "tr")]
                yield f"{path.name} / table {k}", rows
    except Exception as e:  # noqa: BLE001
        yield f"{path.name} (unreadable: {type(e).__name__})", []


def interesting(rows):
    """Gene-like first columns with a numeric column, or a header mentioning coefficient."""
    if not rows:
        return False
    hdr = " ".join(str(c) for r in rows[:3] for c in r)
    g = sum(1 for r in rows if r and GENEISH.match(str(r[0]).strip()))
    num = sum(1 for r in rows if any(NUMISH.match(str(c).strip()) for c in r[1:]))
    return (re.search(r"coef|β|beta|lasso", hdr, re.I) and g >= 2) or (2 <= g <= 60 and num >= g * 0.6 and len(rows) <= 80)


x = pd.read_excel(ROOT / "03_extraction" / "ai_prefill_PRAD.xlsx", sheet_name="逐篇", dtype=str)
todo = x[x["状态"].isin(WANT)]
made = 0
for _, r in todo.iterrows():
    d = FT / str(r["PMCID"])
    xmlf = next(d.glob("PMC*.xml"), None) if d.exists() else None
    if xmlf is None:
        xmlf = ROOT / "03_extraction" / "fulltext_xml" / f"{r['PMCID']}.xml"
    root = ET.fromstring(xmlf.read_text(encoding="utf-8"))
    body = root.find(".//body")
    lines = [f"# PMID {r['PMID']} · {r['PMCID']}", f"**{r['标题']}** — {r['期刊']} ({r['年份']})", f"预提取状态：{r['状态']}", ""]
    paras = [p for p in body.iter("p")]
    hits = [txt(p) for p in paras if KEY.search(txt(p))]
    lines.append("## 相关段落")
    total = 0
    for h in hits:
        if total > 9000:
            lines.append("…（其余段落省略）"); break
        lines.append("- " + h[:1500]); total += min(len(h), 1500)
    lines.append("\n## 正文表格（含基因样行的）")
    for tw in root.iter("table-wrap"):
        rows = table_rows(tw)
        if geneish_rows(rows) >= 2 or re.search(r"coef|β|lasso|signature|risk", txt(tw.find('caption')), re.I):
            lines.append(f"**{txt(tw.find('label'))}** {txt(tw.find('caption'))[:200]}\n\n{md_table(rows)}\n")
    lines.append("\n## 图注（涉及系数/LASSO/森林图/风险评分）")
    for fg in root.iter("fig"):
        cap = txt(fg.find("caption"))
        if re.search(r"coef|lasso|forest|risk ?score|multivariate|signature|model", cap, re.I):
            g = fg.find(".//graphic")
            href = g.get("{http://www.w3.org/1999/xlink}href") if g is not None else ""
            lines.append(f"- **{txt(fg.find('label'))}** [{href}] {cap[:600]}")
    lines.append("\n## 补充材料")
    if d.exists():
        sup = [p for p in d.iterdir() if p.is_file() and not p.name.startswith("PMC") and p.suffix.lower() not in (".jpg", ".jpeg", ".png", ".gif")]
        lines.append("文件：" + ", ".join(p.name for p in sup) if sup else "（无可读补充文件）")
        for p in sup:
            for label, rows in supp_tables(p):
                if interesting(rows):
                    lines.append(f"\n**{label}**\n\n{md_table([[str(c) for c in rr] for rr in rows])}\n")
    else:
        lines.append("（未下载 OA 包）")
    (OUT / f"{r['PMID']}.md").write_text("\n".join(lines), encoding="utf-8")
    made += 1
print("digests", made)
