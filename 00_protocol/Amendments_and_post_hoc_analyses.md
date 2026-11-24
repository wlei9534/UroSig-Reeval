# UroSig-Reeval: amendments to the statistical analysis plan and post hoc analyses

This document lists everything that was decided or changed after the statistical analysis plan (SAP v2) was finalised on 11 November 2026, and states for each item whether it happened before or after results had been seen. It is a summary in English of the dated project log (`logs/decision_log.md`, in Chinese), which is released with the code.

## Version of this document

- Version 1, 20 November 2026 (https://doi.org/10.5281/zenodo.24385588).
- Version 2, 20 November 2026. Two changes: the first version called this record an "OSF project", although the documents are deposited in Zenodo, and that wording is corrected in section 1; and the correction of 20 November to post hoc analysis V2 is added in sections 3 and 4. The original plan (PDF) is unchanged and has the same checksum.
- Version 3, November 2026. Additions only: a note on the software used for the Cox models and the meta-analysis (section 2), and the clarifications of reporting made when the manuscript was checked before submission (section 7). No analysis and no result was changed. The original plan (PDF) and its translation are unchanged.

## 1. Status of the plan

- The SAP was drafted on 1 November 2026 and finalised on 11 November 2026, before any outcome or score result of any cohort had been read. The finalised PDF (SHA-256 `dc9965d5662bfbb408ffbd68f5ef90c57991e4a9765f0a326f7f54745694f110`) is posted unchanged.
- The plan was **not registered prospectively**. Its first public record is this Zenodo record (https://doi.org/10.5281/zenodo.24385587), first published on 20 November 2026 after the analyses were complete. The Zenodo timestamp therefore shows only that the documents existed on that date.
- The first extraction lock (1 November 2026: 223 analysable signatures) preceded the SAP. The primary analysis of those signatures was run on 11 November 2026.

## 2. Implementation notes made after finalisation and before outcomes were read (11 November 2026)

- Uno's C and the time-dependent AUC were computed with our own NumPy implementation instead of scikit-survival (not installed on the analysis computer). On simulated data the two agree to within 1e-15 (scikit-survival 0.28; `05_code/tests/test_timeauc.py`).
- The evaluation was split into parallel shards with fixed seeds.
- Software. The SAP named pandas and lifelines. Harrell's C was computed with lifelines. The Cox models (hazard ratio per standard deviation, clinical models and, later, the calibration slope) were fitted with our own NumPy implementation (Newton–Raphson, Efron's method for ties; `05_code/urosig/cox.py`), and the random-effects meta-analysis (restricted maximum likelihood, Hartung–Knapp) was implemented in NumPy (`05_code/urosig/meta.py`). Both files date from 26 October 2026, before the SAP was drafted, so this is not a change made after results were seen; it is recorded here because the SAP did not mention it. The released tests compare the Cox fit with lifelines (coefficients and log-likelihood within 1e-4 on simulated data) and the REML estimate with direct maximisation of the restricted likelihood. Versions used: Python 3.12.0, NumPy 2.3.5, pandas 3.0.1, SciPy 1.16.1, lifelines 0.30.0. (Note added in version 3.)

## 3. Changes made after the primary results had been seen

| Date (2026) | Change | Reason |
|---|---|---|
| 11 Nov | In the unstandardised-expression sensitivity analysis, 23 missing expression values (7 genes, GSE70768 and GSE70769) are replaced by the gene's cohort mean. | The run stopped on missing values. Equivalent to the value 0 used for z-scores in the primary analysis, which is unaffected. |
| 11 Nov | The reported test of incremental value is a combined likelihood-ratio test (sum of the per-cohort χ² statistics, degrees of freedom equal to the number of cohorts). The originally coded version (smallest per-cohort p value, then Benjamini–Hochberg) is kept in the supplementary tables. | Noticed after seeing the summary figures that the smallest p value is too lenient when a signature has several cohorts. |
| 12 Nov | Two signatures with a pooled C below 0.5 (PMIDs 35143414 and 42418409) were re-read against the source articles by two independent AI runs blinded to the result, by W.L. and by a second urologist. Extraction and calculation were confirmed; the locked data were not changed. | Check requested after seeing the result. |
| 14–17 Nov | In the incremental-value bootstrap the score is standardised within each resample and non-converged resamples are skipped. | Numerical overflow for one signature with very large coefficients read from a figure. |
| 18 Nov | Two pooling corrections in the rerun: the endpoint classification table of the second lock is used for sensitivity analysis S6, and the meta-regression covariates are read for all signatures. | Stale file paths found when checking the output. |
| 20 Nov | In the figure of post hoc analysis V2, the smallest C that a cohort can distinguish from 0.5 is computed from the bootstrap standard error of logit C (0.62 to 0.71 across the four cohorts). It had been drawn as 0.5 plus the minimum detectable difference between signature and grade (0.64 to 0.73), which is a different quantity. The minimum detectable difference from grade (0.14 to 0.23) is unchanged and is reported as such. | Mislabelled quantity found when the additional files of the article were prepared. No other result is affected. |

## 4. Analyses added after the primary results had been seen

All were written into the log before they were run. None is used for the three primary population-level results.

- **Exploratory, specified 12 November.** E1: paired bootstrap difference between the C of the signature and the C of the Gleason grade group in the same resample, pooled per signature. E2: a proliferation benchmark (meta-PCNA; the 168 genes most correlated with PCNA in TCGA-PRAD, the reference cohort, without use of outcome; after Venet et al. 2011, PMID 22028643), scored with equal weights.
- **Exploratory, listed in the SAP and run on 12 November.** E3: metastasis in GSE116918. E4: relation between the discrimination reported by each article and the external result. The reported values were extracted for this purpose after the primary analysis.
- **Post hoc, specified 12 November.** V1: check of the implementation by recomputing each signature in the cohorts its article had used, and comparison of reported training performance, recomputed training performance and independent performance. V2: precision available in each cohort (minimum detectable difference in C between signature and grade, from the paired bootstrap of E1; and, from 20 November, the smallest C distinguishable from 0.5, see section 3). V3: meta-regression of the pooled C on characteristics of the signatures. V4: correlation between signature scores in TCGA-PRAD.
- **Post hoc, added 18 November.** V5: calibration slope, the Cox coefficient of the published linear predictor computed on unstandardised log2 expression. The SAP stated that calibration would not be assessed. Absolute risk still cannot be assessed, and the slope depends on the expression scale of each platform, so it is reported by cohort and interpreted with that limitation.
- Kaplan–Meier curves in the figures are descriptive, as stated in the SAP.

## 5. Completion of the literature screening and second extraction lock (12 to 14 November 2026)

This was decided on 12 November 2026, **after the primary results for the first 223 signatures were known**.

- **What was found.** While the flow diagram was being checked, we found that the 292 articles taken to full text in the first stage had been selected by a keyword rule and that the planned duplicate screening of titles and abstracts had not been carried out. The other 719 of the 1,011 records had not been screened.
- **Screening.** The 719 records were screened independently by two runs of an AI model (Claude, Anthropic), blinded to each other and to the keyword label, against the original criteria. Agreement was 98.6% (Cohen's κ = 0.97). Both excluded 386 records, 323 went to full text, and 10 disagreements were settled by W.L. and a second urologist (all excluded).
- **Full text.** Of the 323 records, 19 could not be retrieved, 2 were preprints duplicating included articles and 2 were not in English; 300 were assessed and 148 excluded, leaving 152 analysable signatures. A preprint with a peer-reviewed version was assessed using the published version and counted once.
- **Extraction.** Genes and coefficients were extracted independently in duplicate by two runs of the AI model, compared item by item, and differences were settled against the source article.
- **Verification by the authors.** Screening decisions and extracted data were then checked item by item against the source articles by three authors (J.Z., S.X. and X.W.).
- **Additional eligibility rules.** The second stage met situations not covered by the first. Rules for them were proposed by the AI model and adopted under the authority of W.L. on 14 November 2026, before any second-stage signature was scored; they are recorded in the log. In brief: scores from mechanistic studies are included if the combined score was tested against survival; principal-component scores only if loadings are reported; direct scoring of predefined gene sets, gene-set or cell-score models, cluster or centroid classifiers, and models applied to dichotomised or binned variables are excluded; for models with clinical variables the gene part is used only if it has separate coefficients; coefficients read from figures or given only as hazard ratios enter sensitivity analysis S7; author-defined equal weights enter S8.
- **Second lock, 14 November 2026.** 375 analysable signatures (247 with published coefficients, 85 equal-weight, 18 from hazard ratios, 16 read from figures, 8 with author-defined equal weights, 1 LASSO-stage). The data of the first lock were not changed.
- **Rerun.** All analyses of the SAP were rerun on the second lock: 211 signatures and 569 signature × cohort evaluations in the primary analysis (first lock: 142 and 384). The median pooled C was 0.575 on the first lock and 0.573 on the second. Both can be reproduced from the released code (`UROSIG_LOCK=v1` or `v2`).

## 6. Use of artificial intelligence

An AI model (Claude, Anthropic) was used for the keyword pre-screening, the duplicate screening of titles and abstracts, the extraction and re-extraction of genes and coefficients, writing and running the analysis code, and drafting project documents including this one. Decisions on eligibility rules, adjudication of disagreements and verification against the source articles are recorded in the log with the person responsible. The authors take responsibility for the content.

## 7. Clarifications of reporting made before submission (23 November 2026)

Before submission every number in the manuscript was recomputed from the result files and every statement of method was checked against the code. No analysis was rerun and no result changed. The following points of reporting were clarified in the manuscript and its additional files.

- **Random gene sets.** The comparison was made for signatures scored as a weighted sum (208 signatures, 564 signature × cohort evaluations). It was not made for the three gene-pair signatures (PMIDs 32593277, 35143414 and 42553741; 5 evaluations), because random sets were generated only for weighted-sum scores. This was so from the first run but had not been stated. The proportion "better than random" is therefore 59.7% of 207 articles (124 of 208 signatures; the two models of PMID 32264916 count as half an article each), and 190 of 564 evaluations (34%) exceeded the 95th percentile. An earlier figure label gave 33%, from a slightly different count (percentile above 95 instead of the prespecified criterion). With Benjamini–Hochberg control at 5%, as prespecified, 106 of 208 signatures differed from random gene sets.
- **Reported against external discrimination (E4, V1).** Of 158 signatures with a reported training value, 129 reported a 3-year AUC, 5 an AUC at another time point and 24 a C-index. The comparison with the pooled external C-index is kept, and the like-for-like comparisons that were computed at the same time are now also reported: reported C-index against external C-index (24 signatures, medians 0.73 and 0.59) and reported 3-year AUC against our AUC at 36 months in independent cohorts (146 signatures, medians 0.76 and 0.60). In V1 the recomputation in the training cohort used the measure the article had reported and biochemical recurrence as recorded in our data. The comparison between cohorts used and not used by an article (151 signatures) leaves TCGA-PRAD out and uses a Wilcoxon signed-rank test.
- **Clinical model.** Its C-index is one value per cohort (0.77 in GSE54460, 0.76 in GSE70769, 0.80 in GSE70768 and 0.83 in MSKCC) and is now reported as that range instead of a median over evaluations. Grade group enters as a linear term; a variable is used in a cohort if it is recorded for at least half of the patients.
- **Prior use of cohorts.** The percentages of prior use of each cohort are computed over the 211 signatures, not over articles; the wording of the figure and the text was corrected.
- **Calibration slope (V5).** Computed for the 564 evaluations of the 208 weighted-sum signatures.
- **Meta-regression (V3).** The model also contains the modelling method, which was not associated with performance; the between-signature variance is estimated by the method of moments.
