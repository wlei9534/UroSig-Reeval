"""Lock extraction v2 (2026-11-14): v1 lock (292 keyword-prioritised papers) + 300 papers from the dual T/A screening.

Usage: python lock_extraction_v2.py <project_root> <hgnc_txt>
Inputs : 03_extraction/locked/signatures_locked_2026-11-01.csv, signature_terms_locked_2026-11-01.csv (v1, unchanged)
         03_extraction/v2/reconciled_v2.csv (one proposed record per v2 paper, after adjudication)
         02_screening/ta_dual/merged_TA.csv, 03_extraction/v2/v2_fulltext_set.csv (metadata)
Outputs: 03_extraction/locked/signatures_locked_v2_2026-11-14.csv
         03_extraction/locked/signature_terms_locked_v2_2026-11-14.csv
         06_results/eligibility/equal_weight_directions_v2.csv (v1 rows + v2 rows)
         03_extraction/locked/MANIFEST_v2.txt
Roles as in v1. sens_hr coefficients are ln(HR) of the printed multivariable HRs (HR_PRINTED below, both extractors agree).
equal_signed coefficients are the author-defined +-1 weights (1 for each equal-vote gene pair).
Directions for equal_weight: the extractors' 'directions:' lists; opposite signs -> 0 (gene dropped); one 0 and one +-1 -> the non-zero one.
"""
import hashlib
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT, HGNC = Path(sys.argv[1]), sys.argv[2]
L = ROOT / "03_extraction" / "locked"
v1s = pd.read_csv(L / "signatures_locked_2026-11-01.csv", dtype=str).fillna("")
v1t = pd.read_csv(L / "signature_terms_locked_2026-11-01.csv", dtype=str).fillna("")
v1d = pd.read_csv(ROOT / "06_results" / "eligibility" / "equal_weight_directions.csv", dtype=str).fillna("")
R = pd.read_csv(ROOT / "03_extraction" / "v2" / "reconciled_v2.csv", dtype=str).fillna("")
adj = pd.read_csv(ROOT / "03_extraction" / "v2" / "adjudication_v2.csv", dtype=str).fillna("")
ta = pd.read_csv(ROOT / "02_screening" / "ta_dual" / "merged_TA.csv", dtype=str).fillna("").set_index("PMID")
ft = pd.read_csv(ROOT / "03_extraction" / "v2" / "v2_fulltext_set.csv", dtype=str).fillna("").set_index("PMID")

# multivariable HRs printed in the papers (X and Y identical), gene order as in reconciled genes
HR_PRINTED = {
    "33189527": "2.12;1.92",
    "34876840": "0.648;0.819;0.901;1.670;0.514;0.506;0.493;0.943;0.278;0.386",
    "35592542": "0.85;1.20;1.26;1.18;0.80",
    "37559353": "0.77;9.05;4.63;2.83;0.44;0.69",
    "37686409": "4.6;2.4",
    "39068564": "1.286;1.528;1.713;2.064",
    "40900420": "0.355;0.575;1.995",
    "41668068": "1.7;2.7;1.7;2.4;1.8",
    "41878699": "1.58;1.86;1.56;1.43;0.54;1.71",
}
# screening record was a preprint; the published version is the one assessed
PUBLISHED = {"37034687": ("38191471", "PMC10774315", "2024", "Nature communications", "检索记录为预印本 PMID 37034687，按正式版评估"),
             "40654646": ("42143470", "", "2026", "Translational oncology", "检索记录为预印本 PMID 40654646，按正式版评估"),
             "40458374": ("40458374", "", "2026", "Molecular cancer", "检索记录为预印本；正式版 Mol Cancer 2026（Article in Press，DOI 10.1186/s12943-026-02795-3）尚无 PMID，沿用预印本 PMID 作标识")}


# obvious misprints noted by both extractors (original spelling kept in 基因（原文写法）)
MAP_FIX = {("30860081", "CRND"): ("CRNDE", "公式笔误，正文其余处均为 CRNDE"),
           ("31741711", "CACNAID"): ("CACNA1D", "摘要笔误，表格为 CACNA1D"),
           ("31961962", "OR,2T2"): ("OR2T2", "原文多一逗号"), ("31961962", "OR,2T11"): ("OR2T11", "原文多一逗号"),
           ("33197890", "DNRB"): ("EDNRB", "Table 1 笔误，正文为 EDNRB"),
           ("32735556", "HOXB8+"): ("HOXB8", "图 4B 标注符号"),
           ("35805202", "NF-κB"): ("RELA", "原文 NF-κB(p65)/RelA")}


def clean_gene(g):
    g = g.strip().replace("–", "-").replace("—", "-").replace("‐", "-")
    return re.sub(r"\s*\(.*?\)\s*", "", g).strip()


def num(x):
    return float(str(x).strip().replace("−", "-").replace("–", "-").replace("‐", "-").replace("—", "-"))


def parse_dirs(note):
    m = re.search(r"directions[：:]\s*([^\n（(。]+)", note)
    out = {}
    if not m:
        return out
    for part in m.group(1).split(";"):
        if ":" in part:
            g, v = part.rsplit(":", 1)
            v = v.strip().replace("−", "-")
            if v in ("+1", "1", "-1", "0"):
                out[clean_gene(g)] = int(v)
    return out


final_status = {}
for _, a in adj.iterrows():
    final_status[a["PMID"]] = a["裁定"]

rows, dir_rows, flags = [], [], []
for _, r in R.iterrows():
    p, role = r["PMID"], r["拟定角色"]
    pid, pmc, yr, jr, pnote = PUBLISHED.get(p, (p, ft["PMCID"].get(p, ""), ta["年份"].get(p, ""), ta["期刊"].get(p, ""), ""))
    genes = [clean_gene(g) for g in r["genes"].split(";") if g.strip()] if role != "excluded" else []
    coefs = []
    if role in ("main", "sens_figure"):
        coefs = [num(c) for c in r["coefs"].split(";") if c.strip()]
    elif role == "sens_hr":
        coefs = [round(math.log(num(h)), 6) for h in HR_PRINTED[p].split(";")]
    elif role == "equal_signed":
        if r["coefs"].strip():
            coefs = [num(c) for c in r["coefs"].split(";") if c.strip()]
        elif all("|" in g for g in genes):
            coefs = [1.0] * len(genes)
        else:
            d = parse_dirs(r["X_note"]) or parse_dirs(r["Y_note"])
            coefs = [float(d.get(g, 0)) for g in genes]
    if coefs:
        assert len(coefs) == len(genes), (p, len(coefs), len(genes))
    if role == "equal_weight":
        dx, dy = parse_dirs(r["X_note"]), parse_dirs(r["Y_note"])
        for g in genes:
            a, b = dx.get(g), dy.get(g)
            if a is not None and b is not None and a != 0 and b != 0 and a != b:
                v, src = 0, f"X={a}/Y={b} 相反，记 0"
                flags.append((pid, g, src))
            elif a is not None and b is not None and (a == 0) != (b == 0):
                v, src = (a or b), f"X={a}/Y={b}：一方未定方向，取另一方"
            else:
                v = a if a is not None else (b if b is not None else 0)
                src = "X、Y 一致" if (a is not None and b is not None) else ("仅一方给出" if (a is not None or b is not None) else "未给出")
            dir_rows.append({"sig_id": pid, "基因（原文写法）": g, "方向": v, "来源": src})
    status = final_status.get(p, "") or (r["X_status"] if r["一致性"] == "完全一致" else "")
    rows.append({"sig_id": pid, "PMID": pid, "PMCID": pmc, "年份": yr, "期刊": jr, "标题": ta["标题"].get(p, ""),
                 "分析角色": role, "最终状态": status, "建议用途": "", "基因数": len(genes),
                 "基因（原文写法）": ";".join(genes), "系数": ";".join(f"{c:g}" if role != "sens_hr" else f"{c}" for c in coefs),
                 "锁定备注": "；".join(x for x in [pnote, ("ln(HR)，HR=" + HR_PRINTED[p]) if role == "sens_hr" else ""] if x),
                 "Wei 决定": "", "第三轮核定": "", "检索记录 PMID": p, "来源版本": "v2（双人筛选补充）",
                 "终点": r["endpoint"], "训练队列": r["training"], "验证队列": r["validation"]})
v2s = pd.DataFrame(rows)
v1s = v1s.assign(**{"检索记录 PMID": v1s["PMID"], "来源版本": "v1（关键词优先 292 篇）", "终点": "", "训练队列": "", "验证队列": ""})
sig = pd.concat([v1s, v2s], ignore_index=True)
assert not sig["sig_id"].duplicated().any(), sig[sig["sig_id"].duplicated()]["sig_id"].tolist()

# HGNC mapping for v2 terms (same rules as v1)
sys.path.insert(0, str(ROOT / "05_code"))
from urosig.genes import build_symbol_map  # noqa: E402

h = pd.read_csv(HGNC, sep="\t", dtype=str, low_memory=False)
smap = build_symbol_map(h)
ens = h.set_index("symbol")["ensembl_gene_id"].fillna("").to_dict()
terms = []
for _, s in v2s[v2s["分析角色"] != "excluded"].iterrows():
    g = s["基因（原文写法）"].split(";")
    c = s["系数"].split(";") if s["系数"] else [""] * len(g)
    for i, (gg, cc) in enumerate(zip(g, c), 1):
        parts = gg.split("|")
        fixes = [MAP_FIX.get((s["sig_id"], x)) for x in parts]
        mapped = [smap.get(f[0], f[0]) if f else (smap.get(x, "") or smap.get(x.upper(), "")) for x, f in zip(parts, fixes)]
        fix = next(((x, f) for x, f in zip(parts, fixes) if f), None)
        fix = (fix[1][0], f"{fix[1][1]}") if fix else None
        terms.append({"sig_id": s["sig_id"], "PMID": s["PMID"], "分析角色": s["分析角色"], "序号": i, "基因（原文写法）": gg,
                      "HGNC 现行符号": "|".join(m or "未映射" for m in mapped),
                      "Ensembl ID": "|".join(ens.get(m, "") for m in mapped if m),
                      "Ensembl 来源": "HGNC 快照" if any(mapped) else "", "系数": cc,
                      "符号映射说明": ("；".join(f"{x}→{f[0]}（{f[1]}）" for x, f in zip(parts, fixes) if f) if fix else "；".join(f"{x}→{m}（旧名/别名）" for x, m in zip(parts, mapped) if m and m != x))})
terms = pd.concat([v1t, pd.DataFrame(terms)], ignore_index=True)
dirs = pd.concat([v1d.assign(来源="v1"), pd.DataFrame(dir_rows)], ignore_index=True)

f1 = L / "signatures_locked_v2_2026-11-14.csv"
f2 = L / "signature_terms_locked_v2_2026-11-14.csv"
f3 = ROOT / "06_results" / "eligibility" / "equal_weight_directions_v2.csv"
sig.to_csv(f1, index=False, encoding="utf-8-sig")
terms.to_csv(f2, index=False, encoding="utf-8-sig")
dirs.to_csv(f3, index=False, encoding="utf-8-sig")
man = "\n".join(f"{hashlib.md5(x.read_bytes()).hexdigest()}  {x.name}" for x in (f1, f2, f3))
(L / "MANIFEST_v2.txt").write_text("提取锁定 v2 2026-11-14（v1 锁定 + 双人筛选补充 300 篇；规则由 Wei 授权 Claude 确定）\n" + man + "\n", encoding="utf-8")
print(pd.crosstab(sig["分析角色"], sig["来源版本"]).to_string())
v2t = terms[terms["sig_id"].isin(v2s["sig_id"])]
print("v2 terms:", len(v2t), "| unmapped:", v2t["HGNC 现行符号"].str.contains("未映射").sum())
print("direction conflicts:", len(flags), flags[:10])
d2 = pd.DataFrame(dir_rows)
print("v2 equal_weight sigs:", d2["sig_id"].nunique(), "| sigs with no nonzero direction:",
      (d2.groupby("sig_id")["方向"].apply(lambda s: (s != 0).sum()) == 0).sum())
print(man)
