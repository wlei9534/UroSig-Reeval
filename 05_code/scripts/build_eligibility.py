"""Leakage matrix, gene directions for equal-weight signatures, and signature x cohort eligibility (no outcomes read).

Usage: python build_eligibility.py <project_root> <workbook> <coverage_csv>
Leakage evidence per signature x cohort:
  - body-text scan of the paper (reference list removed) for cohort identifiers / names, and
  - the blind reviewers' recorded training/validation cohorts;
  where the two agree, that is the call; where they disagree (21 signatures), the source-text verification in
  review/out_LEAK.jsonl decides (训练/验证/其他 = used; 未用 = not used).
Eligibility (external cohorts; TCGA-PRAD is reference only): not used by the paper AND gene coverage >= threshold
  (and, for coefficient signatures, weight coverage >= threshold). Computed at 0.8 and 1.0.
Outputs to 06_results/eligibility/.
"""
import glob
import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT, WB, COV = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
X = ROOT / "03_extraction"
OUT = ROOT / "06_results" / "eligibility"
OUT.mkdir(parents=True, exist_ok=True)
COH = {"GSE54460": r"GSE54460", "GSE70769": r"GSE70769|GSE70770|Stockholm", "GSE70768": r"GSE70768|GSE70770|Cambridge",
       "MSKCC": r"GSE21032|GSE21034|GSE21035|MSKCC|Memorial Sloan|prad_mskcc|Taylor(?! ?& ?Francis)",
       "GSE116918": r"GSE116918", "TCGA": r"TCGA"}
COH_OLD = dict(COH, GSE70769=r"GSE70769|GSE70770", GSE70768=r"GSE70768|GSE70770")  # as used to pick the 21 for review
CMAP = {"MSKCC_GSE21032": "MSKCC", "TCGA-PRAD": "TCGA"}

sig = pd.read_csv(X / "locked" / "signatures_locked_2026-11-01.csv", dtype=str).fillna("")
sig = sig[sig["分析角色"] != "excluded"]
w = pd.read_excel(WB, sheet_name="逐篇", dtype=str).fillna("").set_index("PMID")
rev = {}
for f in sorted(glob.glob(str(X / "review" / "out_[AM]*.jsonl"))) + [str(X / "review" / "out_FIG.jsonl")]:
    for l in open(f, encoding="utf-8"):
        r = json.loads(l)
        rev[r["pmid"]] = r
leak = {json.loads(l)["sig_id"]: json.loads(l) for l in open(X / "review" / "out_LEAK.jsonl", encoding="utf-8")}


def body_pdf(p):
    t = open(X / "pdf60" / f"{p}.raw.txt", encoding="utf-8", errors="replace").read()
    m = list(re.finditer(r"\n\s*(?:References|REFERENCES|Reference|Bibliography)\s*\n", t))
    return t[: m[-1].start()] if m else t


def body_xml(pmcid):
    import xml.etree.ElementTree as ET
    root = ET.fromstring((X / "fulltext_xml" / f"{pmcid}.xml").read_text(encoding="utf-8"))
    for parent in list(root.iter()):
        for ch in list(parent):
            if ch.tag == "ref-list":
                parent.remove(ch)
    return " ".join(" ".join(el.itertext()) for tag in ("body", "back") for el in root.iter(tag))


rows = []
for _, s in sig.iterrows():
    p = s["PMID"]
    pdf = w.loc[p, "全文"].startswith("Wei")
    text = body_pdf(p) if pdf else body_xml(s["PMCID"])
    r = rev.get(p, {})
    rtxt = r.get("training", "") + " | " + r.get("validation", "")
    scan_old = {k for k, rx in COH_OLD.items() if re.search(rx, text)}
    revc = {k for k, rx in COH_OLD.items() if re.search(rx, rtxt)}
    row = {"sig_id": s["sig_id"], "PMID": p, "分析角色": s["分析角色"]}
    if s["sig_id"] in leak:
        L = leak[s["sig_id"]]
        for k in COH:
            row[k] = "用过" if L[k] in ("训练", "验证", "其他") else "未用"
            row[k + "_用途"] = L[k]
        row["判定来源"] = "回原文核定（out_LEAK）"
        row["依据"] = L.get("evidence", "") + ("；" + L["note"] if L.get("note") else "")
    else:
        assert scan_old == revc or s["sig_id"] in leak, s["sig_id"]
        for k in COH:
            row[k] = "用过" if k in scan_old else "未用"
            row[k + "_用途"] = ""
        row["判定来源"] = "正文扫描与复核员记录一致"
        row["依据"] = f"复核员记录 训练：{r.get('training', '')}｜验证：{r.get('validation', '')}"
    rows.append(row)
lk = pd.DataFrame(rows)
lk.to_csv(OUT / "leakage_matrix.csv", index=False, encoding="utf-8-sig")

# directions for equal-weight signatures
dirs = {json.loads(l)["sig_id"]: json.loads(l) for l in open(X / "review" / "out_DIR.jsonl", encoding="utf-8")}
dr = []
for sid, d in dirs.items():
    for g, v in zip(d["genes"].split(";"), d["directions"].split(";")):
        dr.append({"sig_id": sid, "基因（原文写法）": g, "方向": int(v)})
dr = pd.DataFrame(dr)
dr.to_csv(OUT / "equal_weight_directions.csv", index=False, encoding="utf-8-sig")
n_dir = dr[dr["方向"] != 0].groupby("sig_id").size()

# eligibility
cov = pd.read_csv(COV, dtype={"sig_id": str, "PMID": str})
cov["coh"] = cov["cohort"].map(lambda c: CMAP.get(c, c))
cov = cov.merge(lk[["sig_id"] + list(COH)], on="sig_id", how="left")
cov["原文用过"] = cov.apply(lambda r: r[r["coh"]] == "用过", axis=1)
cov["有可评分基因"] = cov.apply(lambda r: r["分析角色"] != "equal_weight" or n_dir.get(r["sig_id"], 0) > 0, axis=1)
for th in (0.8, 1.0):
    wc = pd.to_numeric(cov["权重覆盖率"], errors="coerce")
    cov[f"可评估_{th}"] = (cov["coh"] != "TCGA") & ~cov["原文用过"] & (cov["基因覆盖率"] >= th) & \
        (wc.isna() | (wc >= th)) & cov["有可评分基因"]
keep = ["sig_id", "PMID", "分析角色", "cohort", "原文用过", "基因项数", "找到", "基因覆盖率", "权重覆盖率", "有可评分基因",
        "可评估_0.8", "可评估_1.0", "缺失基因"]
cov[keep].to_csv(OUT / "eligibility_signature_by_cohort.csv", index=False, encoding="utf-8-sig")

summ = []
for th in (0.8, 1.0):
    e = cov[cov[f"可评估_{th}"]]
    for role, g in cov.groupby("分析角色"):
        ids = set(g["sig_id"])
        ok = e[e["sig_id"].isin(ids)].groupby("sig_id").size()
        summ.append({"门槛": th, "分析角色": role, "signature 数": len(ids), "≥1 个可评估外部队列": len(ok),
                     "≥2 个": int((ok >= 2).sum()), "无可评估队列": len(ids) - len(ok),
                     "可评估 signature×队列": int(ok.sum())})
summ = pd.DataFrame(summ)
summ.to_csv(OUT / "eligibility_summary.csv", index=False, encoding="utf-8-sig")
used = lk[list(COH)].eq("用过").sum()
print("原文用过各队列的 signature 数：", used.to_dict())
print(summ.to_string(index=False))
print("equal_weight directions: zero-direction signatures:", sorted(set(dirs) - set(n_dir.index)))
