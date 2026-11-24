"""Lock the signature extraction (end of the extraction phase, 2026-11-01).

Usage: python lock_extraction.py <project_root> <workbook> <hgnc_txt> <out_dir>
Writes to <out_dir>:
  signatures_locked_2026-11-01.csv       one row per signature (167 analysable + excluded papers for the flow diagram)
  signature_terms_locked_2026-11-01.csv  one row per gene term of every analysable signature
  MANIFEST.txt                           md5 of both files
Analysis role (from Wei's decisions):
  main            coefficients printed in the text/tables/supplements (incl. 35414328 refit, gene pairs)
  sens_figure     coefficients read from a figure
  sens_hr         coefficients = ln(HR)
  sens_lasso      LASSO-stage coefficients of a non-linear final model
  equal_signed    direction-only, equal weights +/-1
  equal_weight    genes only; equal-weight score (direction still to be taken from the paper)
  excluded        not analysable (reason in status)
"""
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT, WB, HGNC, OUT = Path(sys.argv[1]), sys.argv[2], sys.argv[3], Path(sys.argv[4])
OUT.mkdir(parents=True, exist_ok=True)
w = pd.read_excel(WB, sheet_name="逐篇", dtype=str).fillna("")

# final value fixes decided by Wei / experts
FIX_COEF = {("36647875", "PAQR6"): ("0.396", "原文印作 396，专家确认为 0.396 的笔误")}
FIX_NOTE = {"35222640": "NCF1 系数前原文缺符号，专家C 确认为正号；公式笔误 COLA1=COL5A1、PLXNDI=PLXND1（专家C 确认）"}


def role(use: str, status: str) -> str:
    if use.startswith("主分析候选"):
        return "main"
    if "从图片读取" in use:
        return "sens_figure"
    if "HR 取对数" in use:
        return "sens_hr"
    if "LASSO 阶段" in use:
        return "sens_lasso"
    if "仅方向" in use:
        return "equal_signed"
    if use.startswith("无系数"):
        return "equal_weight"
    return "excluded"


rows = []
for _, r in w.iterrows():
    reviewed = r["复核后状态"] != ""
    status = r["复核后状态"] if reviewed else r["当前状态（2026-10-31）"]
    genes = r["复核后基因"] if reviewed else (r["人工精读基因"] or r["基因"])
    coefs = r["复核后系数"] if reviewed else (r["人工精读系数"] or r["系数"])
    rl = role(r["建议用途"], status)
    if rl == "excluded":
        genes = coefs = ""
    note = FIX_NOTE.get(r["PMID"], "")
    if genes and coefs:
        g, c = genes.split(";"), coefs.split(";")
        assert len(g) == len(c), r["PMID"]
        for i, gg in enumerate(g):
            if (r["PMID"], gg) in FIX_COEF:
                c[i], why = FIX_COEF[(r["PMID"], gg)]
                note = (note + "；" if note else "") + f"{gg} {why}"
        coefs = ";".join(c)
    rows.append({"sig_id": r["PMID"], "PMID": r["PMID"], "PMCID": r["PMCID"], "年份": r["年份"], "期刊": r["期刊"],
                 "标题": r["标题"], "分析角色": rl, "最终状态": status, "建议用途": r["建议用途"],
                 "基因数": len(genes.split(";")) if genes else 0, "基因（原文写法）": genes, "系数": coefs,
                 "锁定备注": note, "Wei 决定": r.get("Wei 决定（2026-11-01）", ""),
                 "第三轮核定": r.get("第三轮核定（2026-11-01）", "")})
for l in open(ROOT / "03_extraction" / "extra_signatures.jsonl", encoding="utf-8"):
    e = json.loads(l)
    base = w.set_index("PMID").loc[e["pmid"]]
    rows.append({"sig_id": e["sig_id"], "PMID": e["pmid"], "PMCID": e["pmcid"], "年份": base["年份"], "期刊": base["期刊"],
                 "标题": base["标题"], "分析角色": "main", "最终状态": e["status"], "建议用途": "主分析候选（同文第二模型，Wei 决定纳入）",
                 "基因数": len(e["genes"].split(";")), "基因（原文写法）": e["genes"], "系数": e["coefs"],
                 "锁定备注": e["note"], "Wei 决定": "两个模型都纳入，文献层面按同一篇聚类", "第三轮核定": ""})
sig = pd.DataFrame(rows)

# HGNC mapping per term (approved > previous > alias; ambiguous previous/alias dropped)
sys.path.insert(0, str(ROOT / "05_code"))
from urosig.genes import build_symbol_map  # noqa: E402

h = pd.read_csv(HGNC, sep="\t", dtype=str, low_memory=False)
smap = build_symbol_map(h)
amb = defaultdict(set)
for col in ("prev_symbol", "alias_symbol"):
    for s, f in zip(h["symbol"], h[col].fillna("")):
        for o in filter(None, (x.strip() for x in f.split("|"))):
            amb[o].add(s)
ens = h.set_index("symbol")["ensembl_gene_id"].fillna("").to_dict()
ens2sym = {v: k for k, v in ens.items() if v}
ENSG_OVERRIDE = {}  # 35414328: Ensembl IDs from the authors' own data
for l in open(ROOT / "03_extraction" / "manual_extraction_pdf53.jsonl", encoding="utf-8"):
    r = json.loads(l)
    if r["pmid"] == "35414328":
        import re
        ids = re.search(r"Ensembl ID（与基因同序）：([A-Z0-9;]+)", r["note"]).group(1).split(";")
        ENSG_OVERRIDE = {("35414328", g): e for g, e in zip(r["genes"].split(";"), ids)}

terms = []
for _, s in sig[sig["分析角色"] != "excluded"].iterrows():
    if not s["基因（原文写法）"]:
        continue
    g = s["基因（原文写法）"].split(";")
    c = s["系数"].split(";") if s["系数"] else [""] * len(g)
    for i, (gg, cc) in enumerate(zip(g, c), 1):
        parts = gg.split("|")
        mapped = [smap.get(p, "") for p in parts]
        e = ENSG_OVERRIDE.get((s["PMID"], gg), "")
        if e:  # authors' Ensembl ID decides the symbol (e.g. their UBXN10-AS1 is not HGNC's PLA2G2C)
            mapped = [ens2sym.get(e, "")]
        terms.append({"sig_id": s["sig_id"], "PMID": s["PMID"], "分析角色": s["分析角色"], "序号": i,
                      "基因（原文写法）": gg, "HGNC 现行符号": "|".join(m or "未映射" for m in mapped),
                      "Ensembl ID": e or "|".join(ens.get(m, "") for m in mapped if m),
                      "Ensembl 来源": "作者补充数据" if e else ("HGNC 快照" if any(mapped) else ""),
                      "系数": cc,
                      "符号映射说明": ("按作者给出的 Ensembl ID 对应" if e else "；".join(f"{p}→{m}（旧名/别名）" for p, m in zip(parts, mapped) if m and m != p))})
terms = pd.DataFrame(terms)

f1 = OUT / "signatures_locked_2026-11-01.csv"
f2 = OUT / "signature_terms_locked_2026-11-01.csv"
sig.to_csv(f1, index=False, encoding="utf-8-sig")
terms.to_csv(f2, index=False, encoding="utf-8-sig")
man = "\n".join(f"{hashlib.md5(p.read_bytes()).hexdigest()}  {p.name}" for p in (f1, f2))
(OUT / "MANIFEST.txt").write_text("提取锁定 2026-11-01（Wei 定稿）\n" + man + "\n", encoding="utf-8")
print(sig["分析角色"].value_counts().to_string())
print("signatures analysable:", (sig["分析角色"] != "excluded").sum(), "| terms:", len(terms),
      "| unmapped terms:", (terms["HGNC 现行符号"].str.contains("未映射")).sum())
print(man)
