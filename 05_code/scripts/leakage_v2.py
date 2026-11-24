"""Leakage evidence for v2 signatures (which of our cohorts each paper used), same logic as v1 build_eligibility.py.

Usage: python leakage_v2.py <project_root>
Per v2 analysable signature x cohort: body-text scan of the paper (reference list removed) vs the two extractors'
recorded training/validation cohorts (X or Y). Where they agree, that is the call; disagreements are written to
03_extraction/v2/leak_check.csv for source verification (out_LEAK_v2.jsonl decides, as in v1).
Output: 06_results/eligibility/leakage_matrix_v2.csv (v1 rows + v2 rows) once all disagreements are verified.
"""
import glob
import json
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

ROOT = Path(sys.argv[1])
V2 = ROOT / "03_extraction" / "v2"
COH = {"GSE54460": r"GSE54460", "GSE70769": r"GSE70769|GSE70770|Stockholm", "GSE70768": r"GSE70768|GSE70770|Cambridge",
       "MSKCC": r"GSE21032|GSE21034|GSE21035|MSKCC|Memorial Sloan|prad_mskcc|Taylor(?! ?& ?Francis)",
       "GSE116918": r"GSE116918", "TCGA": r"TCGA"}
sig = pd.read_csv(ROOT / "03_extraction" / "locked" / "signatures_locked_v2_2026-11-14.csv", dtype=str).fillna("")
sig = sig[(sig["来源版本"].str.startswith("v2")) & (sig["分析角色"] != "excluded")]
b = pd.concat([pd.read_csv(f, dtype=str).fillna("").assign(batch=Path(f).stem[-3:]) for f in sorted(glob.glob(str(V2 / "batches" / "batch_V*.csv")))])
b = b[b["batch"] != "V16"].drop_duplicates("pmid", keep="last").set_index("pmid")  # V15 holds the published versions


def body(path):
    if path.endswith(".xml"):
        root = ET.fromstring(Path(path).read_text(encoding="utf-8"))
        for parent in list(root.iter()):
            for ch in list(parent):
                if ch.tag == "ref-list":
                    parent.remove(ch)
        return " ".join(" ".join(el.itertext()) for tag in ("body", "back") for el in root.iter(tag))
    t = subprocess.run(["pdftotext", "-layout", path, "-"], capture_output=True, text=True).stdout
    m = list(re.finditer(r"\n\s*(?:References|REFERENCES|Reference|Bibliography)\s*\n", t))
    return t[: m[-1].start()] if m else t


rows, check = [], []
for _, s in sig.iterrows():
    p = s["检索记录 PMID"]
    rows.append((s["sig_id"], p, s["分析角色"], body(b.loc[p, "main_text"]), ""))

cmp_ = pd.read_csv(V2 / "compare_v2.csv", dtype=str).fillna("").set_index("PMID")
out = []
for sid, p, role, text, _ in rows:
    c = cmp_.loc[p]
    rtxt = " | ".join([c["X_training"], c["X_validation"], c["Y_training"], c["Y_validation"]])
    scan = {k for k, rx in COH.items() if re.search(rx, text)}
    rev = {k for k, rx in COH.items() if re.search(rx, rtxt)}
    row = {"sig_id": sid, "PMID": sid, "分析角色": role}
    for k in COH:
        row[k] = "用过" if k in scan else "未用"
        row[k + "_用途"] = ""
    row["判定来源"] = "正文扫描与提取员记录一致" if scan == rev else "待回原文核定"
    row["依据"] = f"提取员记录 X 训练：{c['X_training']}｜X 验证：{c['X_validation']}｜Y 训练：{c['Y_training']}｜Y 验证：{c['Y_validation']}"
    out.append(row)
    if scan != rev:
        check.append({"sig_id": sid, "检索记录 PMID": p, "main_text": b.loc[p, "main_text"], "supp_dirs": b.loc[p, "supp_dirs"],
                      "仅正文扫描命中": ";".join(sorted(scan - rev)), "仅提取员记录": ";".join(sorted(rev - scan)),
                      "X 训练": c["X_training"], "X 验证": c["X_validation"], "Y 训练": c["Y_training"], "Y 验证": c["Y_validation"]})
lk = pd.DataFrame(out)
lk.to_csv(V2 / "leakage_v2_scan.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(check).to_csv(V2 / "leak_check.csv", index=False, encoding="utf-8-sig")
print("v2 analysable:", len(lk), "| agree:", (lk["判定来源"] != "待回原文核定").sum(), "| to verify:", len(check))
print(lk[list(COH)].eq("用过").sum().to_dict())

# apply source verification (out_LEAK_v2.jsonl) and write the combined v1 + v2 leakage matrix
vf = V2 / "out_LEAK_v2.jsonl"
if vf.exists():
    ver = {json.loads(l)["sig_id"]: json.loads(l) for l in open(vf, encoding="utf-8")}
    missing = set(pd.DataFrame(check)["sig_id"]) - set(ver) if check else set()
    assert not missing, f"unverified: {sorted(missing)}"
    for i, row in lk.iterrows():
        v = ver.get(row["sig_id"])
        if not v:
            continue
        for k, use in v["cohorts"].items():
            lk.at[i, k] = "用过" if use in ("训练", "验证", "其他") else "未用"
            lk.at[i, k + "_用途"] = use
        lk.at[i, "判定来源"] = "回原文核定（out_LEAK_v2）"
        lk.at[i, "依据"] = v["evidence"] + "｜" + row["依据"]
    v1 = pd.read_csv(ROOT / "06_results" / "eligibility" / "leakage_matrix.csv", dtype=str).fillna("")
    allk = pd.concat([v1, lk], ignore_index=True)
    allk.to_csv(ROOT / "06_results" / "eligibility" / "leakage_matrix_v2.csv", index=False, encoding="utf-8-sig")
    print("leakage_matrix_v2:", len(allk), "rows;", lk[list(COH)].eq("用过").sum().to_dict())
