"""Collect the primary analysis and S1-S10 into one results workbook (SAP v2 sections 8-10).

Usage: python build_report_tables.py <project_root> <out.xlsx>
Sheets: 汇总对比 (primary + each sensitivity, with the paired difference to primary on common signatures),
主分析_逐个signature, 主分析_组合明细, 队列基准, 无法独立验证, S9_全公式, 探索性.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from urosig.paths import p as P, im_globs, _V as PATHS  # noqa: E402
from scipy import stats

ROOT, OUTX = Path(sys.argv[1]), sys.argv[2]
EV = P(ROOT, "ev")
sig = pd.read_csv(P(ROOT, "sig"), dtype=str).fillna("")
ep = pd.read_csv(P(ROOT, "endpoint"), dtype={"sig_id": str})
elig = pd.read_csv(P(ROOT, "elig"), dtype={"sig_id": str})
MAIN4 = ["GSE54460", "GSE70769", "GSE70768", "MSKCC_GSE21032"]

TAGS = [("primary", "主分析：4 个根治术队列，z 分数，覆盖 ≥80%"),
        ("S1_cov100", "S1 基因覆盖 100%"),
        ("S2a_raw_log2", "S2a 未标准化 log2"),
        ("S2b_tcga_sd", "S2b TCGA 标准差还原尺度（原文用过 TCGA 者）"),
        ("S3_plus_GSE116918", "S3 加入 GSE116918（放疗）"),
        ("S4_without_MSKCC", "S4 去掉 MSKCC"),
        ("S5_GSE54460_only", "S5 只用 RNA-seq 队列 GSE54460"),
        ("S6_endpoint_BCR_RFS", "S6 原文终点为 BCR/RFS 类"),
        ("S7_plus_low_precision", "S7 加入图片读数/HR 换算/LASSO 阶段系数"),
        ("S8_equal_weight", "S8 等权评分 signature 单独汇总"),
        ("S10_GSE54460_no_repeats", "S10 GSE54460 去掉 6 例重复文库病人")]


def fmt_med_iqr(x):
    return f"{np.median(x):.3f} ({np.percentile(x, 25):.3f}–{np.percentile(x, 75):.3f})"


prim = pd.read_csv(EV / "primary" / "signatures.csv", dtype={"sig_id": str}).set_index("sig_id")
rows = []
for tag, label in TAGS:
    f = EV / tag / "summary.json"
    if not (f.exists() and (EV / tag / "signatures.csv").exists()):
        rows.append({"分析": label, "状态": "未完成"})
        continue
    sm = json.loads(f.read_text(encoding="utf-8"))
    s = pd.read_csv(EV / tag / "signatures.csv", dtype={"sig_id": str}).set_index("sig_id")
    common = s.index.intersection(prim.index)
    d = (s.loc[common, "pooled_c"] - prim.loc[common, "pooled_c"]) if tag != "primary" else pd.Series(dtype=float)
    rows.append({
        "分析": label, "状态": "完成", "signature": sm["n_signatures"], "文章": sm["n_papers"], "signature×队列": sm["n_pairs"],
        "≥2 个队列": sm.get("n_signatures_k2plus"),
        "汇总 C 中位数（IQR）": fmt_med_iqr(s["pooled_c"]),
        "95%CI 不含 0.5": sm.get("n_ci_excludes_0.5"),
        "超过 Gleason/ISUP 分级的文章比例": round(sm["share_beats_grade_papers"], 3),
        "优于随机（>P95）的文章比例": None if sm["share_beats_random_papers"] is None else round(sm["share_beats_random_papers"], 3),
        "对随机 FDR<5% 的 signature": sm["n_random_q_lt_0.05"],
        "增量价值 合并LR FDR<5%": sm.get("n_lr_combined_q_lt_0.05"),
        "增量价值 最小p FDR<5%（原代码，偏宽松）": sm["n_lr_q_lt_0.05"],
        "校正 ΔC 中位数": None if sm.get("median_delta_c_corrected") is None else round(sm["median_delta_c_corrected"], 4),
        "与主分析共同 signature 数": len(common) if tag != "primary" else None,
        "汇总 C 与主分析之差 中位数（IQR）": fmt_med_iqr(d) if len(d) else None})
summary = pd.DataFrame(rows)

meta = sig.set_index("sig_id")[["PMID", "年份", "期刊", "标题", "分析角色", "基因数"]]
per = prim.join(meta.drop(columns="PMID"), how="left").join(ep.set_index("sig_id")[["终点类别"]], how="left")
per = per.sort_values("pooled_c", ascending=False).reset_index()
keep = ["sig_id", "PMID", "年份", "期刊", "标题", "基因数", "终点类别", "k_cohorts", "cohorts", "pooled_c", "pooled_c_lo", "pooled_c_hi",
        "i2", "pooled_grade_c", "pooled_sig_c_same_pts", "beats_grade", "median_uno_c", "median_auc36", "median_hr_sd",
        "mean_delta_c_corrected", "lr_p_combined", "lr_q_bh_combined", "rand_p_pooled", "rand_q_bh", "beats_random_95"]
per = per[[c for c in keep if c in per.columns]]

pairs = pd.read_csv(EV / "run_primary" / "pairs.csv", dtype={"sig_id": str})
pcols = ["sig_id", "cohort", "kind", "n", "events", "n_terms", "n_used", "gene_cov", "weight_cov", "harrell_c", "harrell_c_lo",
         "harrell_c_hi", "uno_c", "uno_c_lo", "uno_c_hi", "auc36", "auc60", "hr_sd", "hr_sd_lo", "hr_sd_hi", "hr_sd_p", "grade_var",
         "grade_c", "sig_c_on_grade_pts", "clin_vars", "c_clin", "c_clin_sig", "delta_c", "delta_c_corrected", "lr_chi2", "lr_p",
         "rand_pct", "rand_p", "beats_random_95"]
pairs = pairs[pcols]

coh = pairs.groupby("cohort").agg(病人数=("n", "first"), BCR事件=("events", "first"), 组合数=("sig_id", "size"),
                                  signature_C中位数=("harrell_c", "median"), 分级变量=("grade_var", "first"),
                                  分级C=("grade_c", "median"), 扩展临床模型变量=("clin_vars", "first"),
                                  扩展临床模型C=("c_clin", "median"), 报告60月AUC的组合比例=("auc60", lambda s: s.notna().mean())).reset_index()

# signatures with no evaluable main cohort, and why
m = elig[(elig["分析角色"] == "main") & elig["cohort"].isin(MAIN4)]
none = []
for sid, g in m.groupby("sig_id"):
    if g["可评估_0.8"].any():
        continue
    used = g["原文用过"].sum()
    lowcov = (~g["原文用过"] & ~g["可评估_0.8"]).sum()
    reason = "原文用尽全部 4 个队列" if used == 4 else ("原文用过的队列外，其余队列基因覆盖不足" if used else "4 个队列基因覆盖均不足")
    none.append({"sig_id": sid, "PMID": sid.split("_")[0], "原文用过的主分析队列数": int(used), "覆盖不足的队列数": int(lowcov),
                 "原因": reason, "可在 GSE116918 评估": bool(elig[(elig.sig_id == sid) & (elig.cohort == "GSE116918")]["可评估_0.8"].any())})
none = pd.DataFrame(none)

s9 = pd.DataFrame()
if (EV / "run_full_formula" / "pairs.csv").exists():
    a = pd.read_csv(EV / "run_full_formula" / "pairs.csv", dtype={"sig_id": str})
    b = pairs[pairs["sig_id"].isin(a["sig_id"])]
    s9 = a[["sig_id", "cohort", "harrell_c", "harrell_c_lo", "harrell_c_hi", "grade_c"]].merge(
        b[["sig_id", "cohort", "harrell_c"]], on=["sig_id", "cohort"], suffixes=("_基因+临床全公式", "_只算基因"))

# exploratory (SAP 8, not for conclusions): pooled C vs number of genes, year, original endpoint
ex = per.dropna(subset=["pooled_c"]).copy()
ex["基因数"] = pd.to_numeric(ex["基因数"], errors="coerce")
ex["年份"] = pd.to_numeric(ex["年份"], errors="coerce")
exr = []
for v in ("基因数", "年份"):
    ok = ex[v].notna()
    rho, p = stats.spearmanr(ex.loc[ok, v], ex.loc[ok, "pooled_c"])
    exr.append({"关系": f"汇总 C vs {v}", "统计量": f"Spearman ρ = {rho:.3f}", "p": p, "n": int(ok.sum())})
g1, g2 = ex.loc[ex["终点类别"] == "BCR/RFS", "pooled_c"], ex.loc[ex["终点类别"] != "BCR/RFS", "pooled_c"]
u = stats.mannwhitneyu(g1, g2)
exr.append({"关系": "汇总 C：原文终点 BCR/RFS 类 vs 其他", "统计量": f"中位数 {g1.median():.3f} vs {g2.median():.3f}（Mann-Whitney）",
            "p": u.pvalue, "n": len(ex)})
exr = pd.DataFrame(exr)

# exploratory E1-E4 (pre-specified 2026-11-12; not for conclusions)
XD = EV / "exploratory"
ex_sheets = {}
if (XD / "summary_E1E2.json").exists():
    e12 = json.loads((XD / "summary_E1E2.json").read_text(encoding="utf-8"))
    s12 = pd.read_csv(XD / "signatures_E1E2.csv", dtype={"sig_id": str})
    rows_e = [{"项目": "E1 signature C − 分级 C（同一批病人、配对 bootstrap，按 signature 汇总）",
               "signature 数": e12["n_signatures"], "汇总 ΔC 中位数": round(e12["E1_grade"]["median_pooled_dC"], 3),
               "CI 整体 >0（优于分级）": e12["E1_grade"]["n_ci_above_0"], "CI 整体 <0（不如分级）": e12["E1_grade"]["n_ci_below_0"],
               "CI 含 0": e12["E1_grade"]["n_ci_includes_0"]},
              {"项目": f"E2 signature C − 增殖评分 meta-PCNA C（{e12['metapcna_n_genes']} 个基因）",
               "signature 数": e12["n_signatures"], "汇总 ΔC 中位数": round(e12["E2_metapcna"]["median_pooled_dC_sig_minus_metapcna"], 3),
               "CI 整体 >0（优于分级）": e12["E2_metapcna"]["n_ci_above_0"], "CI 整体 <0（不如分级）": e12["E2_metapcna"]["n_ci_below_0"],
               "CI 含 0": e12["n_signatures"] - e12["E2_metapcna"]["n_ci_above_0"] - e12["E2_metapcna"]["n_ci_below_0"]}]
    e = pd.DataFrame(rows_e).rename(columns={"CI 整体 >0（优于分级）": "CI 整体 >0（signature 更好）", "CI 整体 <0（不如分级）": "CI 整体 <0（signature 更差）"})
    mpc = pd.read_csv(XD / "metapcna_cohorts.csv").rename(columns={"c_metapcna": "meta-PCNA C", "c_grade": "分级 C", "c_sig_median": "signature C 中位数"})
    ex_sheets["E1E2_对分级与增殖"] = (e, mpc, s12, e12)
if (EV / "E3_met_GSE116918" / "summary.json").exists():
    e3 = pd.read_csv(EV / "E3_met_GSE116918" / "signatures.csv", dtype={"sig_id": str}).sort_values("pooled_c", ascending=False)
    ex_sheets["E3_转移_GSE116918"] = e3
if (XD / "summary_E4.json").exists():
    ex_sheets["E4_原文自报vs外部"] = pd.read_csv(XD / "internal_reported.csv", dtype={"sig_id": str})

with pd.ExcelWriter(OUTX, engine="openpyxl") as w:
    summary.to_excel(w, sheet_name="汇总对比", index=False)
    per.to_excel(w, sheet_name="主分析_逐个signature", index=False)
    pairs.to_excel(w, sheet_name="主分析_组合明细", index=False)
    coh.to_excel(w, sheet_name="队列基准", index=False)
    none.to_excel(w, sheet_name="无法独立验证", index=False)
    if len(s9):
        s9.to_excel(w, sheet_name="S9_全公式", index=False)
    exr.to_excel(w, sheet_name="探索性（不作结论）", index=False)
    if "E1E2_对分级与增殖" in ex_sheets:
        e, mpc, s12, e12 = ex_sheets["E1E2_对分级与增殖"]
        e.to_excel(w, sheet_name="E1E2_对分级与增殖", index=False)
        mpc.to_excel(w, sheet_name="E1E2_对分级与增殖", index=False, startrow=len(e) + 3)
        s12.to_excel(w, sheet_name="E1E2_逐个signature", index=False)
    if "E3_转移_GSE116918" in ex_sheets:
        ex_sheets["E3_转移_GSE116918"].to_excel(w, sheet_name="E3_转移_GSE116918", index=False)
    if "E4_原文自报vs外部" in ex_sheets:
        ex_sheets["E4_原文自报vs外部"].to_excel(w, sheet_name="E4_原文自报vs外部", index=False)
    for ws in w.book.worksheets:
        ws.freeze_panes = "B2"
        for col in ws.columns:
            width = max(len(str(c.value or "")) for c in col[:60])
            ws.column_dimensions[col[0].column_letter].width = min(max(10, width * 1.1), 60)
print(summary.to_string())
print(coh.to_string())
print(none["原因"].value_counts().to_string() if len(none) else "none")
print(exr.to_string())
