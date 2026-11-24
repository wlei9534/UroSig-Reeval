"""Literal check of manually extracted genes/coefficients against the source full text and supplements.

Usage: python verify_manual.py <project_root> [extra_supp_dir ...]
Reads 03_extraction/manual_extraction.jsonl; corpus = PMC JATS XML text + every readable supplement
(xlsx/xls/csv/tsv/txt/docx/pdf, including files inside zip archives) under 03_extraction/fulltext/<PMCID>/.
Writes 03_extraction/manual_extraction_verify.csv with one row per gene/coefficient term.

Checks per term:
  gene_found   gene symbol (or both genes of a pair) occurs as a whole word in the corpus
  coef_found   |coefficient| occurs as a number in the corpus (numeric match, so 0.0960 == 0.096)
  near         gene and |coefficient| occur within 200 characters of each other somewhere in the corpus
Records whose values were read from figures or converted from HR are flagged, not failed: their numbers
are not expected to be in the text.
"""
from __future__ import annotations

import csv
import io
import json
import re
import subprocess
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(sys.argv[1])
EXTRA = [Path(p) for p in sys.argv[2:]]
FT = ROOT / "03_extraction" / "fulltext"
XD = ROOT / "03_extraction" / "fulltext_xml"
NORM = str.maketrans({"−": "-", "–": "-", "—": "-", " ": " ", " ": " ", "\xa0": " ", "‐": "-", "‑": "-"})
NUM = re.compile(r"(?<![\w.])-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


def docx_text(data: bytes) -> str:
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        x = z.read("word/document.xml").decode("utf-8", "replace")
    x = re.sub(r"</w:p>|</w:tc>", " \n ", x)
    return re.sub(r"<[^>]+>", "", x)


def xlsx_text(data: bytes, name: str) -> str:
    import pandas as pd
    try:
        sheets = pd.read_excel(io.BytesIO(data), sheet_name=None, header=None, dtype=str)
    except Exception:  # noqa: BLE001  (some .xls are really tab text)
        return data.decode("utf-8", "replace")
    return "\n".join(" | ".join(str(c) for c in row) for df in sheets.values() for row in df.fillna("").values.tolist())


def pdf_text(data: bytes) -> str:
    r = subprocess.run(["pdftotext", "-layout", "-", "-"], input=data, capture_output=True)
    return r.stdout.decode("utf-8", "replace")


def file_text(name: str, data: bytes) -> str:
    low = name.lower()
    try:
        if low.endswith(".docx"):
            return docx_text(data)
        if low.endswith((".xlsx", ".xlsm", ".xls")):
            return xlsx_text(data, name)
        if low.endswith(".pdf"):
            return pdf_text(data)
        if low.endswith((".csv", ".tsv", ".txt")):
            return data.decode("utf-8", "replace")
        if low.endswith(".zip"):
            out = []
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                for m in z.namelist():
                    if not m.endswith("/") and z.getinfo(m).file_size < 60e6:
                        out.append(file_text(m, z.read(m)))
            return "\n".join(out)
    except Exception as e:  # noqa: BLE001
        return f"[unreadable {name}: {type(e).__name__}]"
    return ""


def corpus(pmcid: str) -> str:
    parts = []
    xf = XD / f"{pmcid}.xml"
    if xf.exists():
        root = ET.fromstring(xf.read_text(encoding="utf-8"))
        parts.append(" ".join(root.itertext()))
    for d in [FT / pmcid] + [e / pmcid for e in EXTRA]:
        if d.is_dir():
            for p in sorted(d.rglob("*")):
                if p.is_file() and not p.name.startswith(pmcid) and p.stat().st_size < 60e6:
                    parts.append(file_text(p.name, p.read_bytes()))
    return re.sub(r"[ \t]+", " ", "\n".join(parts).translate(NORM))


def numbers(text: str):
    """(position, |value|) for every number in the text."""
    out = []
    for m in NUM.finditer(text):
        try:
            out.append((m.start(), abs(float(m.group()))))
        except ValueError:
            pass
    return out


def same(a: float, b: float) -> bool:
    return abs(a - b) <= 1e-12 * max(1.0, abs(a))


FIG = ("从图片读取", "由 HR 换算")
recs = [json.loads(l) for l in open(ROOT / "03_extraction" / "manual_extraction.jsonl", encoding="utf-8")]
rows = []
for r in recs:
    if not r["genes"]:
        continue
    pmcid = r.get("pmcid", "")
    text = corpus(pmcid) if pmcid else ""
    nums = numbers(text)
    squashed = re.sub(r"\s+", "", text)  # MathML and table markup can split "ADH5" into "ADH 5"
    sq_nums = numbers(squashed)
    genes = r["genes"].split(";")
    coefs = r["coefs"].split(";") if r["coefs"] else [""] * len(genes)
    from_fig = any(k in r["status"] for k in FIG)
    for g, c in zip(genes, coefs):
        parts = g.split("|")
        gpos = []
        ok_gene = True
        for gg in parts:
            body = "[-_ ]?".join(re.escape(seg) for seg in gg.split("-"))
            pat = re.compile(rf"(?:(?<![A-Za-z0-9])|(?<=exp)|(?<=Exp)){body}(?![A-Za-z0-9])", re.I)
            hits = [m.start() for m in pat.finditer(text)]
            sq_hits = [m.start() for m in pat.finditer(squashed)]
            ok_gene &= bool(hits or sq_hits)
            gpos.append((hits, sq_hits))
        coef_found = near = ""
        if c and "仅方向" not in r["status"]:
            v = abs(float(c))
            cpos = [p for p, x in nums if same(x, v)]
            sq_cpos = [p for p, x in sq_nums if same(x, v)]
            coef_found = bool(cpos or sq_cpos)
            near = any(abs(p - q) <= 200 for p in cpos for q in gpos[-1][0]) or \
                any(abs(p - q) <= 120 for p in sq_cpos for q in gpos[-1][1])
        rows.append({"PMID": r["pmid"], "PMCID": pmcid, "状态": r["status"], "基因": g, "系数": c,
                     "基因在原文": ok_gene, "系数在原文": coef_found, "基因与系数相邻(≤200字符)": near,
                     "数值来源": "图片/HR换算（原文无此数，需对图核对）" if from_fig else "正文/表格/补充材料"})

out = ROOT / "03_extraction" / "manual_extraction_verify.csv"
with open(out, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
txt_rows = [x for x in rows if x["数值来源"].startswith("正文")]
print("terms", len(rows), "| text-sourced", len(txt_rows))
print("gene not found:", [(x["PMID"], x["基因"]) for x in rows if not x["基因在原文"]])
print("coef not found (text-sourced):", [(x["PMID"], x["基因"], x["系数"]) for x in txt_rows if x["系数"] != "" and x["系数在原文"] is False])
print("not adjacent (text-sourced):", [(x["PMID"], x["基因"], x["系数"]) for x in txt_rows if x["系数在原文"] is True and not x["基因与系数相邻(≤200字符)"]])
print("saved", out)
