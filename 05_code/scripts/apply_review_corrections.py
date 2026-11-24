"""Patch extractor-1 records with values the blind review showed were missed, each re-checked against the source.

Usage: python apply_review_corrections.py <project_root>
Reads 03_extraction/review/out_*.jsonl (blind review) and rewrites, in place, the matching records of
03_extraction/manual_extraction.jsonl and manual_extraction_pdf53.jsonl. The previous versions are kept as
*.before_review.jsonl. Every corrected record gets a note saying what changed and why.
"""
import glob
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(sys.argv[1])
X = ROOT / "03_extraction"
STAMP = "2026-11-01 盲法复核更正"

# PMID -> (new status, reason). genes/coefs are taken from the blind review record (verified in the source text).
TAKE = {
    "37565774": ("已提取", "提取者1误判为无系数；PDF Results 3.2 有公式，pdftotext 丢了负号，渲染页面核对"),
    "32566639": ("已提取", "提取者1误判为无系数；摘要与 Results 3.3 有 CMU5 公式"),
    "35669428": ("已提取", "提取者1误判为无系数；Results 有 PRRS 公式（原文写 BAK，即 BAK1）"),
    "36248980": ("已提取", "提取者1误判为无系数；Results 有风险评分公式"),
    "36275733": ("已提取", "提取者1误判为无系数；Results 有 CRPGS 公式"),
    "37008945": ("已提取", "提取者1误判为无系数；Results 有公式（原文笔误 APCC=ACPP，摘要 TMEM388=TMEM38B）"),
    "37033267": ("已提取", "提取者1误判为无系数；Results 有 MSRPI 公式（XML 中基因名被拆成单个字母）"),
    "41261613": ("已提取", "提取者1误判为无系数；Results 3.4 有 CRGs-RS 公式"),
    "41758298": ("已提取", "提取者1误判为无系数；Results 有 M2GS 公式（另有常数项 −0.009，不影响排序）"),
    "39780215": ("已提取", "提取者1误取同文 LGG（胶质瘤）模型；前列腺癌模型为 GRB2、RAD51"),
    "34769200": ("已提取（仅方向，等权 ±1）", "等权 z 分数求和，DNTT 方向取反；改为方向记法"),
    "41281377": ("仅基因（无系数）", "补全 9 个基因名（见 Discussion）；Lasso+plsRcox 未报告权重"),
    # from the figure-enabled re-review (out_FIG)
    "32782607": ("已提取", "提取者1误用 Figure 4A 的 HR 换算；Results 正文本有风险评分公式（2–3 位小数），改用正文系数，进主分析"),
    "36312753": ("已提取", "提取者1误取 Figure 2M 的 3 个关键基因；主推模型为正文公式的 18 基因 LASSO 模型（注意交叉验证图纵轴为 Binomial Deviance，可能按 logistic 拟合）"),
    "35559023": ("已提取（由 HR 换算）", "提取者1误取 Figure 3C 的 m6A 中间步骤；主推为 6 基因 panel，TCGA 训练系数未报告，唯一数值为补充图 S1 合并队列（654 例、17 个事件）多因素 Cox HR，系数 = ln(HR)，仅敏感性分析"),
}
FROM_FIG = {"32782607", "36312753", "35559023"}
NORM = str.maketrans({"−": "-", "–": "-", "‐": "-", "‑": "-"})

rev = {}
for f in sorted(glob.glob(str(X / "review" / "out_*.jsonl"))):
    if f.endswith("out_FIG.jsonl"):
        continue
    for l in open(f, encoding="utf-8"):
        if l.strip():
            r = json.loads(l)
            rev[str(r["pmid"])] = r
for l in open(X / "review" / "out_FIG.jsonl", encoding="utf-8"):
    if l.strip():
        r = json.loads(l)
        if str(r["pmid"]) in FROM_FIG:
            rev[str(r["pmid"])] = r

done = []
for name in ("manual_extraction.jsonl", "manual_extraction_pdf53.jsonl"):
    p = X / name
    bak = p.with_suffix(".before_review.jsonl")
    if not bak.exists():
        shutil.copy2(p, bak)
    shutil.copy2(bak, p)  # always start from the pre-review version, so re-running is safe
    recs = [json.loads(l) for l in open(p, encoding="utf-8")]
    for r in recs:
        if r["pmid"] not in TAKE:
            continue
        st, why = TAKE[r["pmid"]]
        b = rev[r["pmid"]]
        old = f"原记录：{r['status']}｜{r['genes']}｜{r['coefs']}"
        r["status"] = st
        r["genes"] = b["genes"].translate(NORM).replace(" ", "")
        r["coefs"] = b["coefs"].translate(NORM).replace(" ", "") if b.get("coefs") else ""
        if b.get("status", "").startswith("仅 HR") and r["coefs"]:  # reviewer gave HR; store ln(HR)
            import math
            r["coefs"] = ";".join(f"{math.log(float(h)):.4f}" for h in r["coefs"].split(";"))
        if r["genes"] and r["coefs"]:
            assert len(r["genes"].split(";")) == len(r["coefs"].split(";")), r["pmid"]
        r["source"] = b["source"]
        r["note"] = f"【{STAMP}】{why}。{old}。" + (r.get("note") or "")
        done.append(r["pmid"])
    open(p, "w", encoding="utf-8").write("".join(json.dumps(r, ensure_ascii=False) + "\n" for r in recs))

missing = set(TAKE) - set(done)
print("corrected", len(done), sorted(done))
assert not missing, missing
