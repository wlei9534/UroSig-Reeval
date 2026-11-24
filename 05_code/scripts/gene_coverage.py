"""Gene coverage of every locked signature in every processed cohort (no outcome data is read).

Usage: python gene_coverage.py <project_root>
Inputs : 03_extraction/locked/signature_terms_locked_2026-11-01.csv
         04_cohorts/processed/<cohort>/expr_log2.csv.gz  (only the gene index column is read)
Output : 06_results/coverage/gene_coverage_by_signature.csv  one row per signature x cohort
         06_results/coverage/missing_terms.csv               every term not found, per cohort
A term is found when its current HGNC symbol is in the cohort's gene index; a gene pair needs both genes.
Weight coverage = sum |coef| of found terms / sum |coef| of all terms (coefficient signatures only).
"""
import gzip
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from urosig.paths import p as P  # noqa: E402

ROOT = Path(sys.argv[1])
COHORTS = ["GSE54460", "GSE70769", "GSE70768", "MSKCC_GSE21032", "GSE116918", "TCGA-PRAD"]
terms = pd.read_csv(P(ROOT, "terms"), dtype=str).fillna("")
out = ROOT / "06_results" / "coverage"
out.mkdir(parents=True, exist_ok=True)

index = {}
for c in COHORTS:
    with gzip.open(ROOT / "04_cohorts" / "processed" / c / "expr_log2.csv.gz", "rt", encoding="utf-8") as f:
        next(f)
        index[c] = {line.split(",", 1)[0].strip('"') for line in f}
    print(c, len(index[c]))


def found(sym: str, idx: set) -> bool:
    parts = sym.split("|")
    return all(p and p != "未映射" and p in idx for p in parts)


rows, miss = [], []
for sig, t in terms.groupby("sig_id", sort=False):
    role = t["分析角色"].iloc[0]
    coef = pd.to_numeric(t["系数"], errors="coerce").abs()
    for c in COHORTS:
        ok = t["HGNC 现行符号"].map(lambda s: found(s, index[c]))
        w = coef.sum()
        rows.append({"sig_id": sig, "PMID": t["PMID"].iloc[0], "分析角色": role, "cohort": c, "基因项数": len(t),
                     "找到": int(ok.sum()), "基因覆盖率": round(ok.mean(), 4),
                     "权重覆盖率": round(coef[ok].sum() / w, 4) if w and w == w else "",
                     "缺失基因": ";".join(t.loc[~ok, "基因（原文写法）"])})
        for _, r in t[~ok].iterrows():
            miss.append({"sig_id": sig, "cohort": c, "基因（原文写法）": r["基因（原文写法）"],
                         "HGNC 现行符号": r["HGNC 现行符号"], "Ensembl ID": r["Ensembl ID"]})
cov = pd.DataFrame(rows)
cov.to_csv(P(ROOT, "cov"), index=False, encoding="utf-8-sig")
pd.DataFrame(miss).to_csv(P(ROOT, "missing"), index=False, encoding="utf-8-sig")
piv = cov.pivot(index="sig_id", columns="cohort", values="基因覆盖率")
print("signatures:", len(piv))
for th in (1.0, 0.9, 0.8):
    print(f"gene coverage >= {th}:", {c: int((piv[c] >= th).sum()) for c in COHORTS})
