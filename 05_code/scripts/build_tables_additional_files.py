"""Build Table 1, Table 2 and the spreadsheet additional files (3, 4 and 6) of the manuscript.

Numbering follows the order of first citation in the article: 3 = records, signatures and exclusions; 4 = use of the cohorts; 6 = results.

Usage: python build_tables_additional_files.py <project_root>

Everything written here is in English. Free-text fields that the extractors recorded in Chinese
(original endpoint, training and validation data, exclusion notes, basis of the cohort-use judgement)
are replaced by coded categories; the original records stay in 03_extraction/ and 06_results/eligibility/.
The script stops if any Chinese character is left in a workbook.
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(sys.argv[1]); OUT = ROOT / "07_submission" / "additional_files"; TAB = ROOT / "07_submission" / "tables"
OUT.mkdir(parents=True, exist_ok=True); TAB.mkdir(parents=True, exist_ok=True)
EV = ROOT / "06_results" / "evaluation_v2"; EL = ROOT / "06_results" / "eligibility"; L = ROOT / "03_extraction" / "locked"
ZH = re.compile("[\\u3000-\\u303f\\u4e00-\\u9fff\\uff00-\\uffef]"); HAN = re.compile("[\\u4e00-\\u9fff]")


def readme(rows):
    return pd.DataFrame(rows, columns=["sheet", "column", "definition"])


def write_book(path, sheets):
    """sheets: list of (name, DataFrame). Refuses to write Chinese text."""
    for name, d in sheets:
        # titles are kept verbatim from PubMed (one carries a full-width comma), so only ideographs count there
        bad = [c for c in d.columns if ZH.search(str(c)) or d[c].map(lambda x, rx=(HAN if c in ("title", "journal") else ZH): isinstance(x, str) and bool(rx.search(x))).any()]
        if bad:
            raise SystemExit("Chinese text left in %s / %s: %s" % (path.name, name, bad))
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        for name, d in sheets:
            d.to_excel(xw, sheet_name=name[:31], index=False)
            ws = xw.sheets[name[:31]]; ws.freeze_panes = "A2"
            for i, c in enumerate(d.columns):
                width = max([len(str(c))] + [len(str(v)) for v in d[c].head(200).tolist()])
                ws.column_dimensions[ws.cell(row=1, column=i + 1).column_letter].width = min(max(10, width + 2), 60)
    print(path.name, [(n, len(d)) for n, d in sheets])


# ---------------- Table 1: cohorts
ORDER = ["GSE54460", "GSE70769", "GSE70768", "MSKCC_GSE21032", "GSE116918", "TCGA-PRAD"]
if (ROOT / "04_cohorts" / "processed" / ORDER[0] / "clinical.csv").exists():
    rows = []
    for c in ORDER:
        d = pd.read_csv(ROOT / "04_cohorts" / "processed" / c / "clinical.csv"); q = json.load(open(ROOT / "04_cohorts" / "processed" / c / "qc.json", encoding="utf-8"))
        u = d[d["time_bcr"].notna() & d["bcr"].notna()]
        isup = pd.to_numeric(u["isup"], errors="coerce"); g = {k: int((isup == k).sum()) for k in (1, 2, 3)}; g45 = int((isup >= 4).sum())
        pt = u["pt_cat"].astype(str); pt_known = u["pt_cat"].notna().sum(); pt3 = int(pt.str.contains("T3|T4").sum()) if pt_known else None
        psa = pd.to_numeric(u["psa"], errors="coerce")
        rows.append({"cohort": c, "patients_with_bcr_data": len(u), "bcr_events": int(u["bcr"].sum()), "median_followup_censored_months": q["median_followup_censored"],
                     "grade_group_1": g[1], "grade_group_2": g[2], "grade_group_3": g[3], "grade_group_4_5": g45, "grade_group_missing": int(isup.isna().sum()),
                     "pT3_or_higher": pt3, "pT_known": int(pt_known), "psa_median": None if psa.notna().sum() == 0 else round(float(psa.median()), 1), "psa_known": int(psa.notna().sum()),
                     "gleason_sum_counts": u["gleason_sum"].value_counts().sort_index().to_dict() if c == "GSE116918" else "", "expression_note": q["notes"][0]})
    pd.DataFrame(rows).to_csv(TAB / "Table1_cohorts.csv", index=False); print("Table1_cohorts.csv written")

# ---------------- shared inputs
sig_all = pd.read_csv(L / "signatures_locked_v2_2026-11-14.csv", dtype=str).fillna("")
terms = pd.read_csv(L / "signature_terms_locked_v2_2026-11-14.csv", dtype=str).fillna("")
sig = sig_all[sig_all["分析角色"] != "excluded"].copy(); terms = terms[terms["sig_id"].isin(sig["sig_id"])].copy()
prim = pd.read_csv(EV / "primary" / "signatures.csv", dtype={"sig_id": str, "PMID": str})
feat = pd.read_csv(EV / "posthoc" / "V3_features.csv", dtype={"sig_id": str})
endp = pd.read_csv(EL / "original_endpoint_class_v2.csv", dtype=str).fillna("")
METHOD = {"多因素 Cox": "multivariable Cox", "LASSO 类": "LASSO-type", "单因素或其他": "univariable or other"}
def TR(x):
    x = METHOD.get(x, x)
    for k, v in (("多变量", "multivariable"), ("单变量:", "univariable: "), ("method_多因素 Cox", "method: multivariable Cox"), ("method_单因素或其他", "method: univariable or other"), ("method_LASSO 类", "method: LASSO-type")):
        x = x.replace(k, v)
    return x


STAGE = lambda x: "first stage" if x.startswith("v1") else "second stage"  # noqa: E731

# ---------------- Table 2: characteristics of the signatures in the primary analysis
p = prim.merge(feat, on="sig_id").merge(sig[["sig_id", "年份", "基因数"]], on="sig_id").merge(endp[["sig_id", "终点类别"]], on="sig_id")
assert len(p) == len(prim) == 211, len(p)
yr = p["年份"].astype(int); ng = p["基因数"].astype(int)
t2 = [("Signatures", "", len(p)), ("Articles", "", p["PMID"].nunique())]
t2 += [("Year of publication", lab, int(((yr >= a) & (yr <= b)).sum())) for lab, a, b in (("2018 or earlier", 0, 2018), ("2019-2020", 2019, 2020), ("2021-2022", 2021, 2022), ("2023-2024", 2023, 2024), ("2025-2026", 2025, 2026))]
t2 += [("Number of genes", lab, int(((ng >= a) & (ng <= b)).sum())) for lab, a, b in (("2-5", 2, 5), ("6-10", 6, 10), ("11-20", 11, 20), ("more than 20", 21, 10 ** 6))]
t2 += [("Number of genes", "median (interquartile range)", "%g (%g-%g)" % (ng.median(), ng.quantile(.25), ng.quantile(.75)))]
t2 += [("Training data", "TCGA-PRAD only", int((p["tcga_only_training"] == 1).sum())), ("Training data", "other", int((p["tcga_only_training"] != 1).sum()))]
t2 += [("Final model", METHOD[k], int((p["method"] == k).sum())) for k in METHOD]
t2 += [("Endpoint of the original model", "biochemical recurrence or recurrence-free survival", int((p["终点类别"] == "BCR/RFS").sum())), ("Endpoint of the original model", "other or unclear", int((p["终点类别"] != "BCR/RFS").sum()))]
t2 += [("External validation reported by the authors", "yes", int((p["ext_reported"] == 1).sum())), ("External validation reported by the authors", "no", int((p["ext_reported"] != 1).sum()))]
t2 += [("Independent cohorts available for re-evaluation", str(k), int((p["k_cohorts"] == k).sum())) for k in (4, 3, 2, 1)]
t2 = pd.DataFrame(t2, columns=["characteristic", "level", "n"]); t2["percent"] = [round(100 * v / len(p), 1) if isinstance(v, int) and lev else "" for v, lev in zip(t2["n"], t2["level"])]
t2.to_csv(TAB / "Table2_signature_characteristics.csv", index=False); print(t2.to_string(index=False))

# ---------------- Additional file 3: records, signatures, gene terms, excluded articles
s3 = sig[["sig_id", "PMID", "PMCID", "年份", "期刊", "标题", "分析角色", "来源版本", "基因数", "基因（原文写法）", "系数"]].rename(columns={
    "年份": "year", "期刊": "journal", "标题": "title", "分析角色": "analysis_role", "来源版本": "selection_stage", "基因数": "n_genes", "基因（原文写法）": "genes_as_written", "系数": "coefficients"})
s3["selection_stage"] = s3["selection_stage"].map(STAGE)
s3 = s3.merge(endp[["sig_id", "终点类别"]].rename(columns={"终点类别": "original_endpoint_class"}), on="sig_id", how="left")
s3["original_endpoint_class"] = s3["original_endpoint_class"].replace({"BCR/RFS": "BCR or RFS"})
f = feat[["sig_id", "tcga_only_training", "method", "ext_reported"]].copy(); f["method"] = f["method"].map(METHOD)
f = f.rename(columns={"tcga_only_training": "trained_on_TCGA_only", "method": "final_model", "ext_reported": "external_validation_reported"})
s3 = s3.merge(f, on="sig_id", how="left"); s3.insert(7, "in_primary_analysis", s3["sig_id"].isin(prim["sig_id"]).map({True: "yes", False: "no"}))
for c in ("trained_on_TCGA_only", "external_validation_reported"):
    s3[c] = s3[c].map({1: "yes", 0: "no", 1.0: "yes", 0.0: "no"}).fillna("")
s3["final_model"] = s3["final_model"].fillna("")


def note_en(x):
    if not x:
        return ""
    m = re.match(r"^(.*?)（(.*)）$", x); head, tail = (m.group(1), m.group(2)) if m else (x, "")
    if x.startswith("按作者给出的 Ensembl ID"):
        return "mapped by the Ensembl ID given by the authors"
    kind = {"旧名/别名": "previous symbol or alias", "原文多一逗号": "stray comma in the article", "图 4B 标注符号": "label symbol in the figure removed",
            "原文 NF-κB(p65)/RelA": "the article specifies NF-kB p65 (RelA)"}.get(tail)
    if kind is None:
        assert "笔误" in tail, tail
        kind = "misprint in the article; the symbol is given correctly elsewhere in the same article"
    return head.replace("→", " -> ") + " (" + kind + ")"


t3 = terms.rename(columns={"分析角色": "analysis_role", "序号": "term_number", "基因（原文写法）": "gene_as_written", "HGNC 现行符号": "hgnc_symbol", "Ensembl ID": "ensembl_id",
                           "Ensembl 来源": "ensembl_source", "系数": "coefficient", "符号映射说明": "mapping_note"})
t3["hgnc_symbol"] = t3["hgnc_symbol"].str.replace("未映射", "not mapped")
t3["ensembl_source"] = t3["ensembl_source"].replace({"HGNC 快照": "HGNC snapshot", "作者补充数据": "authors' supplementary data"})
t3["mapping_note"] = t3["mapping_note"].map(note_en)

RULES = [("未获取全文：中文", "Chinese-language article, not assessed"), ("未获取全文", "Full text not retrievable"), ("预印本", "Preprint without a peer-reviewed version, or duplicate of a published article"),
         ("非前列腺癌", "Not prostate cancer"), ("非线性模型", "Non-linear model or threshold rule; no linear score to recompute"), ("非锁定 signature", "Non-linear model or threshold rule; no linear score to recompute"),
         ("非新 signature", "Not a new signature (validation or reuse of an existing signature or gene set)"), ("非预后 signature", "Not a prognostic signature"),
         ("非肿瘤组织", "Not derived from tumour tissue"), ("非转录组", "Not transcriptomic"), ("非 mRNA 级", "Not at the mRNA gene level (for example circRNA, splicing events, tRNA fragments)"),
         ("非基因级线性", "Not a gene-level linear score (pathway or gene-set level)"), ("非多基因", "Not a multigene signature"),
         ("无法判断", "Not reproducible: genes or weights not available"), ("无法复现", "Not reproducible: genes or weights not available"), ("原文图表", "Not reproducible: the article contradicts itself")]


def reason(s):
    s = s.replace("排除：", "")
    for k, v in RULES:
        if s.startswith(k) or (k == "预印本" and k in s):
            return v
    raise SystemExit("unmapped exclusion note: " + s)


ex = sig_all[sig_all["分析角色"] == "excluded"].copy(); ex["reason"] = ex["最终状态"].map(reason); ex["stage"] = ex["来源版本"].map(STAGE)
e3 = ex[["PMID", "年份", "期刊", "标题", "stage", "reason"]].rename(columns={"年份": "year", "期刊": "journal", "标题": "title", "stage": "selection_stage", "reason": "reason_for_exclusion"})
print(pd.crosstab(e3["reason_for_exclusion"], e3["selection_stage"]).to_string())

sheets3 = []
scr_path = ROOT / "08_public_release" / "UroSig-Reeval" / "02_screening" / "screening_decisions.csv"
if scr_path.exists():
    scr = pd.read_csv(scr_path, dtype=str).fillna("")
    by_rec = {}
    for _, r in sig_all.iterrows():
        by_rec.setdefault(r["检索记录 PMID"] or r["PMID"], []).append(r)
    TA = {"排除": "excluded by both screeners", "排除（裁定）": "excluded after adjudication of a disagreement", "进入全文": "to full text", "": "to full text by the keyword rule (first stage)"}
    rec_rows = []
    for _, r in scr.iterrows():
        hits = by_rec.get(r["PMID"], []); ta = TA[r["结果"]]
        if hits:
            ok = [h for h in hits if h["分析角色"] != "excluded"]
            out = "analysable signature (%d)" % len(ok) if ok else "excluded at full text: " + reason(hits[0]["最终状态"])
            assessed = hits[0]["PMID"] if hits[0]["PMID"] != r["PMID"] else ""
        elif r["结果"] == "进入全文":
            out = "Chinese-language article, not assessed" if r["期刊"].startswith("Zhonghua") else ("duplicate preprint of an article already in the search results" if "bioRxiv" in r["期刊"] or "medRxiv" in r["期刊"] else "full text not retrievable")
            assessed = ""
        else:
            out, assessed = "", ""
        rec_rows.append({"record": r["编号"], "PMID": r["PMID"], "DOI": r["DOI"], "year": r["年份"], "journal": r["期刊"], "title": r["标题"],
                         "selection_stage": "first stage" if r["筛选阶段"].startswith("第一") else "second stage",
                         "screener_A": {"排除": "exclude", "纳入": "include", "不确定": "unsure", "": ""}[r["A_决定"]], "screener_B": {"排除": "exclude", "纳入": "include", "不确定": "unsure", "": ""}[r["B_决定"]],
                         "title_abstract_result": ta, "full_text_outcome": out, "assessed_as_PMID": assessed})
    rec = pd.DataFrame(rec_rows)
    print(rec["title_abstract_result"].value_counts().to_dict()); print(rec["full_text_outcome"].str.replace(r"\(\d\)", "", regex=True).str.split(":").str[0].value_counts().to_dict())
    assert len(rec) == 1011 and rec["full_text_outcome"].str.startswith("analysable").sum() == 374
    sheets3.append(("Records", rec))
else:
    print("screening_decisions.csv not found: the Records sheet is left out")
rd3 = readme([
    ("Records", "(sheet)", "One row per PubMed record (1,011). The two screeners of the second stage were two independent runs of an AI model; see Additional file 1."),
    ("Records", "title_abstract_result", "Result of title and abstract screening. First-stage records went to full text by a keyword rule without duplicate screening."),
    ("Records", "full_text_outcome", "Outcome of full-text assessment; empty for records excluded on title and abstract."),
    ("Records", "assessed_as_PMID", "PMID of the peer-reviewed version when a preprint record was assessed as its published article."),
    ("Signatures", "(sheet)", "One row per analysable signature (375 signatures from 374 articles)."),
    ("Signatures", "sig_id", "Signature identifier: PMID of the article, with a suffix when an article contributes more than one signature."),
    ("Signatures", "analysis_role", "main: coefficients published; sens_hr: weights are log hazard ratios; sens_figure: coefficients read from a figure; sens_lasso: LASSO coefficients before the final model; equal_weight: genes only, direction from the article; equal_signed: genes with signs only."),
    ("Signatures", "in_primary_analysis", "yes for the 211 signatures with role main and at least one independent cohort with 80% gene coverage or more."),
    ("Signatures", "genes_as_written / coefficients", "Gene symbols as printed in the article and the matching coefficients, separated by semicolons."),
    ("Signatures", "original_endpoint_class", "Endpoint of the original model: BCR or RFS, other, or unclear."),
    ("Signatures", "trained_on_TCGA_only / final_model / external_validation_reported", "Coded for the 211 signatures of the primary analysis (meta-regression covariates); empty otherwise."),
    ("Gene terms", "(sheet)", "One row per gene term of each analysable signature, with the current HGNC symbol and Ensembl gene ID used for scoring."),
    ("Gene terms", "mapping_note", "Why the symbol used differs from the symbol printed in the article."),
    ("Excluded at full text", "(sheet)", "The 218 articles excluded at full-text assessment (70 in the first stage, 148 in the second), with the reason category. The verbatim notes of the extractors are in the public repository."),
])
sheets3 += [("Signatures", s3), ("Gene terms", t3), ("Excluded at full text", e3), ("README", rd3)]
write_book(OUT / "Additional_file_3_signatures.xlsx", sheets3)

# ---------------- Additional file 6: results
SIGDEF = [("sig_id / PMID", "Signature identifier and PMID of the article."), ("k_cohorts / cohorts", "Number and names of the cohorts in which the signature was evaluated."),
          ("pooled_c, pooled_c_lo, pooled_c_hi", "Random-effects pooled Harrell C-index (pooled on the logit scale) and 95% confidence interval; for one cohort, the cohort estimate."),
          ("i2 / tau2", "Between-cohort heterogeneity of the logit C-index."), ("pooled_grade_c", "Pooled C-index of ISUP grade group in the same cohorts, same weights."),
          ("pooled_sig_c_same_pts", "Pooled C-index of the signature in the patients with grade available."), ("beats_grade", "pooled_sig_c_same_pts greater than pooled_grade_c."),
          ("mean_delta_c_corrected", "Weighted mean over cohorts of the optimism-corrected change in C-index when the signature is added to the clinical model."),
          ("min_lr_p / lr_chi2_sum / lr_p_combined", "Likelihood-ratio test of the signature added to the clinical model: smallest P over cohorts, sum of chi-square statistics and its P value (degrees of freedom = number of cohorts)."),
          ("median_uno_c / median_auc36 / median_hr_sd", "Median over cohorts of Uno's C, the time-dependent AUC at 36 months and the hazard ratio per standard deviation of the score."),
          ("ci_excludes_0.5", "Lower limit of the pooled confidence interval above 0.5."),
          ("rand_p_pooled / beats_random_95 / rand_q_bh", "Comparison with expression-matched random gene sets carrying the same coefficients: empirical P of the pooled C-index, whether it exceeds the 95th percentile of the random sets, and the Benjamini-Hochberg q value."),
          ("lr_q_bh / lr_q_bh_combined", "Benjamini-Hochberg q values of min_lr_p (as first coded; anti-conservative, kept for transparency) and of lr_p_combined (reported).")]
PAIRDEF = [("status / mode / kind", "Evaluation status; scoring mode (z: genes standardised within cohort); kind of score (linear: weighted sum; pair: gene-pair rule)."),
           ("n / events", "Patients and biochemical recurrences in the evaluation."), ("n_terms / n_used / gene_cov / weight_cov", "Gene terms of the signature, terms found in the cohort, share of terms found and share of the sum of absolute coefficients found."),
           ("harrell_c, uno_c, auc36, auc60 (with _lo, _hi)", "Harrell's C, Uno's C and time-dependent AUC at 36 and 60 months with bootstrap 95% confidence intervals."),
           ("logit_c / logit_c_se", "Logit of Harrell's C and its bootstrap standard error, as pooled."), ("hr_sd, hr_sd_lo, hr_sd_hi, hr_sd_p", "Hazard ratio per standard deviation of the score, Wald 95% confidence interval and P value."),
           ("grade_var / grade_c / sig_c_on_grade_pts", "Grade variable, its C-index, and the C-index of the signature in the patients with grade available."),
           ("clin_vars / c_clin / c_clin_sig / delta_c / delta_c_corrected", "Variables of the clinical model, its C-index, the C-index after adding the signature, the difference and the optimism-corrected difference (bootstrap)."),
           ("lr_chi2 / lr_p / n_inc / ev_inc", "Likelihood-ratio test of adding the signature to the clinical model, and the patients and events with complete clinical variables."),
           ("rand_pct / rand_p / beats_random_95", "Percentile of the signature among expression-matched random gene sets, empirical P, and whether it exceeds their 95th percentile."), ("seed", "Seed of the random-number generator for this evaluation.")]
sheets4 = [("Primary signatures", prim), ("Primary evaluations", pd.read_csv(EV / "run_primary" / "pairs.csv", dtype={"sig_id": str}))]
rd4 = [("Primary signatures", c, d) for c, d in SIGDEF] + [("Primary evaluations", c, d) for c, d in PAIRDEF]
LABEL = {"S1_cov100": "Sensitivity: signatures with 100% gene coverage only", "S2a_raw_log2": "Sensitivity: scores from unstandardised log2 expression", "S2b_tcga_sd": "Sensitivity: genes rescaled to their TCGA-PRAD standard deviation",
         "S3_plus_GSE116918": "Sensitivity: radiotherapy cohort GSE116918 added", "S4_without_MSKCC": "Sensitivity: without the MSKCC cohort", "S5_GSE54460_only": "Sensitivity: GSE54460 only",
         "S6_endpoint_BCR_RFS": "Sensitivity: signatures whose original endpoint was BCR or RFS", "S7_plus_low_precision": "Sensitivity: signatures with weights from hazard ratios, figures or LASSO added",
         "S8_equal_weight": "Sensitivity: signatures published without coefficients, scored with equal weights", "S10_GSE54460_no_repeats": "Sensitivity: GSE54460 without the six patients who had repeat libraries with discordant follow-up",
         "E3_met_GSE116918": "Exploratory: metastasis as the outcome in GSE116918"}
for d in sorted(q for q in EV.iterdir() if q.is_dir() and (q / "signatures.csv").exists() and q.name[0] in "SE" and q.name != "primary"):
    sheets4.append((d.name, pd.read_csv(d / "signatures.csv", dtype={"sig_id": str, "PMID": str}))); rd4.append((d.name, "(sheet)", LABEL.get(d.name, d.name) + ". Columns as in Primary signatures."))
e12 = pd.read_csv(EV / "exploratory" / "signatures_E1E2.csv", dtype={"sig_id": str}); sheets4.append(("E1 E2 grade and proliferation", e12))
rd4 += [("E1 E2 grade and proliferation", "d_grade, d_grade_lo, d_grade_hi", "Pooled difference in C-index, signature minus ISUP grade group, with 95% confidence interval."),
        ("E1 E2 grade and proliferation", "d_mp, d_mp_lo, d_mp_hi", "Pooled difference in C-index, signature minus the meta-PCNA proliferation index, with 95% confidence interval."),
        ("E1 E2 grade and proliferation", "median_rho_metapcna", "Median over cohorts of the Spearman correlation between the signature score and the meta-PCNA index.")]
for fn, name, defs in (("V1_reported_vs_recomputed.csv", "V1 reported vs recomputed", "Post hoc: performance reported by the authors against our recomputation."), ("V1_three_step.csv", "V1 three step", "Post hoc: reported training performance, our recomputation in the training cohort, and external performance."),
                       ("V1_used_cohort_recomputation.csv", "V1 used vs unused cohorts", "Post hoc: C-index in cohorts the article had used against cohorts it had not used."),
                       ("V2_mdd_by_cohort.csv", "V2 minimum detectable C", "Post hoc: smallest C-index each cohort could distinguish from 0.5."), ("V2_pooled_ci_halfwidth_by_k.csv", "V2 CI half-width by cohorts", "Post hoc: half-width of the pooled confidence interval by number of cohorts."),
                       ("V3_features.csv", "V3 covariates", "Post hoc: covariates of the meta-regression (y: logit pooled C-index; v: its variance)."), ("V3_metaregression.csv", "V3 meta-regression", "Post hoc: random-effects meta-regression of the logit pooled C-index."),
                       ("V4_clusters_TCGA.csv", "V4 clusters in TCGA", "Post hoc: clusters of correlated signature scores in TCGA-PRAD."), ("V5_calibration_slope_pairs.csv", "V5 calibration by cohort", "Post hoc: calibration slope of each signature in each cohort (lp_sd: standard deviation of the linear predictor)."),
                       ("V5_calibration_slope_signatures.csv", "V5 calibration pooled", "Post hoc: calibration slope pooled over cohorts.")):
    if (EV / "posthoc" / fn).exists():
        d = pd.read_csv(EV / "posthoc" / fn, dtype={"sig_id": str, "PMID": str})
        d = d.rename(columns=lambda c: {"分析角色": "analysis_role", "提取状态": "extraction_status"}.get(c, c))
        if "reported_endpoint" in d:  # endpoint as named by the article; the extractors' Chinese remarks in brackets are dropped
            d["reported_endpoint"] = d["reported_endpoint"].fillna("").map(lambda x: re.sub(r"（[^）]*）", "", x).replace("未说明", "not stated").replace("原文写 ", "").replace("BCR或远处转移", "BCR or distant metastasis").strip())
        for c in d.columns:  # value labels; written without relying on the column dtype (object in pandas 2, str in pandas 3)
            d[c] = d[c].map(lambda x: TR(x) if isinstance(x, str) else x)
        sheets4.append((name, d)); rd4.append((name, "(sheet)", defs))
sheets4.append(("README", readme(rd4)))
write_book(OUT / "Additional_file_6_results.xlsx", sheets4)

# ---------------- Additional file 4: use of the cohorts by each article, gene coverage and evaluability
lk = pd.read_csv(EL / "leakage_matrix_v2.csv", dtype=str).fillna("")
SRC = {"正文扫描与复核员记录一致": "full-text scan and extraction record agree", "正文扫描与提取员记录一致": "full-text scan and extraction record agree",
       "回原文核定（out_LEAK_v2）": "disagreement resolved by re-reading the article", "回原文核定（out_LEAK）": "disagreement resolved by re-reading the article"}
COH = ["GSE54460", "GSE70769", "GSE70768", "MSKCC", "GSE116918", "TCGA"]
w6 = lk[["sig_id", "PMID", "分析角色"]].rename(columns={"分析角色": "analysis_role"})
for c in COH:
    w6["TCGA-PRAD" if c == "TCGA" else c] = lk[c].map({"用过": "used", "未用": "not used"})
w6["n_cohorts_used"] = (w6[[("TCGA-PRAD" if c == "TCGA" else c) for c in COH]] == "used").sum(axis=1); w6["source_of_judgement"] = lk["判定来源"].map(SRC)
assert w6.notna().all().all()
el = pd.read_csv(EL / "eligibility_signature_by_cohort_v2.csv", dtype=str).fillna("")
l6 = el.rename(columns={"分析角色": "analysis_role", "原文用过": "used_by_article", "基因项数": "n_terms", "找到": "n_found", "基因覆盖率": "gene_coverage", "权重覆盖率": "weight_coverage",
                        "有可评分基因": "has_scorable_genes", "可评估_0.8": "evaluable_at_80pct", "可评估_1.0": "evaluable_at_100pct", "缺失基因": "missing_genes"})
rd6 = readme([("Use of cohorts by article", "(sheet)", "One row per analysable signature: whether the article had used each of the six cohorts in any way (training, validation or other). A cohort used by the article is not treated as independent for that signature."),
              ("Use of cohorts by article", "source_of_judgement", "Agreement between an automated scan of the full text for cohort identifiers and the record of the extractors; disagreements were resolved by re-reading the article."),
              ("Signature by cohort", "(sheet)", "One row per signature and cohort: gene coverage and whether the combination was evaluable."),
              ("Signature by cohort", "gene_coverage / weight_coverage", "Share of gene terms found in the cohort, and share of the sum of absolute coefficients found."),
              ("Signature by cohort", "evaluable_at_80pct / evaluable_at_100pct", "Cohort not used by the article and coverage at or above the threshold (80% in the primary analysis).")])
write_book(OUT / "Additional_file_4_cohort_use.xlsx", [("Use of cohorts by article", w6), ("Signature by cohort", l6), ("README", rd6)])
print("used counts", {c: int((w6[c] == "used").sum()) for c in w6.columns[3:9]})
