External validation of 211 published transcriptomic signatures for biochemical recurrence of prostate cancer: comparison with Gleason grade and random gene sets.

This file gives the search strategy, the eligibility criteria, the flow of records with reasons for exclusion, details of the cohorts and of the scoring, the sensitivity, exploratory and post hoc analyses, and every amendment to the statistical analysis plan. The plan itself (original and English translation) is archived at <https://doi.org/10.5281/zenodo.24385587>. Folder names such as `03_extraction/locked` refer to the public repository.

# 1. Search strategy

PubMed was searched on 26 October 2026. The search returned 1011 records. No other database was searched.

```
("prostate cancer"[tiab] OR "prostate adenocarcinoma"[tiab] OR PRAD[tiab]
   OR "Prostatic Neoplasms"[mh])
AND (signature*[tiab] OR "gene panel"[tiab] OR "risk score*"[tiab]
   OR "risk model*"[tiab] OR "prognostic model*"[tiab] OR riskscore[tiab]
   OR classifier*[tiab])
AND (prognos*[tiab] OR survival[tiab] OR recurrence[tiab]
   OR "Kaplan-Meier"[tiab] OR "Cox"[tiab])
AND (gene*[tiab] OR transcriptom*[tiab] OR mRNA*[tiab] OR lncRNA*[tiab]
   OR "RNA-seq"[tiab] OR microarray[tiab] OR TCGA[tiab]
   OR "Gene Expression Omnibus"[tiab] OR GEO[tiab])
AND ("biochemical recurrence"[tiab] OR BCR[tiab] OR "recurrence-free"[tiab]
   OR "metastasis-free"[tiab] OR survival[tiab])
NOT (review[pt] OR "systematic review"[pt] OR meta-analysis[pt]
   OR editorial[pt] OR comment[pt] OR letter[pt] OR "published erratum"[pt]
   OR "retracted publication"[pt] OR "case reports"[pt])
```

# 2. Eligibility criteria

**Inclusion.** Original research that proposed a new prognostic signature, risk score or risk model for prostate cancer, built from the expression of two or more genes (messenger RNA or long non-coding RNA) in tumour tissue, with any prognostic endpoint (biochemical recurrence, metastasis, progression-free, disease-free, cancer-specific or overall survival). Pan-cancer studies were eligible only if they reported a separate model for prostate cancer. Gene-pair models were eligible.

**Exclusion at title and abstract** (codes used by the screeners): E1, not original research; E2, not prostate cancer, or pan-cancer without a separate prostate model; E3, not a tissue transcriptomic signature (microRNA, methylation, protein or immunohistochemistry, radiomics, pathomics, blood, urine or other body-fluid markers, germline variants); E4, single-gene marker; E5, diagnosis, subtyping or staging only, with no prognostic endpoint; E6, not a human study, or cell lines or animals only; E7, other. Screening at this stage was inclusive: a record was marked include, unsure or exclude, and only records excluded by both screeners were dropped without adjudication.

**Exclusion at full text.** Not a new signature (validation or reuse of an existing or commercial classifier, or direct scoring of a predefined gene set); not a prognostic signature; final model not linear in gene expression (for example random survival forest, gradient boosting, support vector machine) or applied to dichotomised or binned expression, so that no weighted sum could be recomputed; not at the gene level (pathway, gene-set or cell-fraction scores); not transcriptomic, or not at the mRNA gene level (circular RNA, splicing events, transfer RNA fragments); single gene; not prostate cancer; not derived from tumour tissue; not reproducible because the genes or the weights could not be obtained from the article and its supplementary material; preprint without a peer-reviewed version; full text not retrievable; not in English.

**Rules added in the second stage.** The second stage met situations that the first had not. The following rules were proposed by the AI model and adopted under the authority of W.L. on 14 November 2026, before any second-stage signature was scored. (i) A multigene score from a mainly mechanistic study was included if the combined score had been tested against survival. (ii) Principal-component scores were included only if the loadings were reported. (iii) Direct scoring of a predefined gene set was not a new signature. (iv) Models whose coefficients applied to dichotomised or binned variables, threshold rules and scores that were not weighted sums were excluded as non-linear, and no further analysis was added for them. (v) For models with clinical variables, the gene part was used if it had separate coefficients; otherwise the model was excluded. (vi) Models built on tissue adjacent to the tumour were excluded; genes selected in blood but scored in tumour tissue were eligible. (vii) Weights read from figures or given only as hazard ratios entered sensitivity analysis S7, and equal weights defined by the authors entered S8.

# 3. Selection process

The selection was done in two stages, and the second was not planned.

**First stage.** A keyword rule applied by the AI model to titles and abstracts labelled 292 of the 1011 records as likely to report a signature, a modelling method and transcriptomic data. These 292 went to full-text assessment. The duplicate screening of titles and abstracts that the protocol had planned was not carried out at this stage. Genes and coefficients were extracted by the AI model and re-extracted by an independent run blinded to the first; differences were settled against the source article. The first lock (1 November 2026) held 223 analysable signatures from 222 articles.

**Second stage.** On 12 November 2026, while the flow diagram was being checked and after the primary results for the first 223 signatures were known, we found that the other 719 records had not been screened. They were then screened independently by two runs of the AI model (Claude, Anthropic), blinded to each other and to the keyword label. The two runs agreed on 98.6% of records (Cohen's kappa 0.97). Both excluded 386 records; 10 disagreements were settled by W.L. and a second urologist, and all 10 were excluded; 323 records went to full text. Of these, 19 could not be retrieved, 2 were preprints that duplicated articles already among the records and 2 were not in English. The remaining 300 were assessed, with genes and coefficients extracted independently in duplicate by two runs of the AI model and compared item by item. A preprint with a peer-reviewed version was assessed in its published form and counted once. The second lock (14 November 2026) added 152 analysable signatures from 152 articles; the data of the first lock were not changed.

**Verification by the authors.** All screening decisions and all extracted genes and coefficients of both stages were checked item by item against the source articles by three authors (J.Z., S.X. and X.W.).

**Table S1. Reasons for exclusion at full-text assessment**

| Reason | First stage | Second stage |
|:-------------------------------------------------------|------:|------:|
| Not a new signature | 8 | 58 |
| Non-linear model or threshold rule | 20 | 29 |
| Not a prognostic signature | 4 | 21 |
| Not reproducible: genes or weights not available | 0 | 14 |
| Not a multigene signature | 13 | 1 |
| Not transcriptomic | 4 | 7 |
| Not at the mRNA gene level | 4 | 7 |
| Not a gene-level linear score | 5 | 3 |
| Not prostate cancer | 5 | 3 |
| Preprint | 2 | 3 |
| Not in English, not assessed | 3 | 0 |
| Full text not retrievable | 2 | 0 |
| Not reproducible: the article contradicts itself | 0 | 1 |
| Not derived from tumour tissue | 0 | 1 |
| **Total** | **70** | **148** |

In the second stage, records that were not retrievable (19), duplicate preprints (2) and non-English records (2) were removed before full-text assessment and are not in the 148. Every record, with its decision at each step, is listed in Additional file 3.

**Table S2. Analysable signatures by source of the weights**

| Analysis role | Source of the weights | Signatures | Use |
|:-----------------|:------------------------------------------------|-----------:|:---------------------------|
| `main` | Coefficients printed in the text, tables or supplementary material | 247 | Primary analysis |
| `sens_hr` | Hazard ratios only; the weight is ln(HR) | 18 | Sensitivity analysis S7 |
| `sens_figure` | Coefficients read from a figure | 16 | Sensitivity analysis S7 |
| `sens_lasso` | LASSO-step coefficients of a non-linear final model | 1 | Sensitivity analysis S7 |
| `equal_weight` | Gene names only; direction taken from the article | 85 | Sensitivity analysis S8 |
| `equal_signed` | Signs only, as defined by the authors | 8 | Sensitivity analysis S8 |
| Total | | 375 | |

Of the 247 signatures with published coefficients, 211 had at least one independent cohort with enough gene coverage and entered the primary analysis; 36 had none (Additional file 5: Fig. S1). In general only the final model promoted by the authors was taken from each article. One article (PMID 32264916) contributed two models, which count as one article in article-level proportions.

# 4. Cohorts

All cohorts are public. The processed expression matrices are on a log2 scale; microarray cohorts were quantile-normalised again within the cohort. Probes or transcripts were mapped to current HGNC symbols (snapshot of 26 October 2026).

- **GSE54460.** RNA sequencing of formalin-fixed, paraffin-embedded prostatectomy tissue from three centres. FPKM values were converted to TPM and to log2(TPM + 1). Six patients had two libraries; the original library was kept, and because the repeat entries carried a different follow-up time these patients were removed in sensitivity analysis S10. This cohort was the development cohort of one signature (PMID 24713434), which was therefore not evaluated in it.
- **GSE70768 and GSE70769.** Illumina HumanHT-12 v4 microarrays; values as normalised by the submitters. Benign samples were removed from GSE70768. An article that used the series GSE70770 was treated as having used both.
- **MSKCC (GSE21032).** Affymetrix Human Exon 1.0 ST arrays, gene-level series matrix. The outcome is the disease-free survival field distributed through cBioPortal, which is not identical to biochemical recurrence; sensitivity analysis S4 repeats the analysis without this cohort. PSA and margin status were not available, so the clinical model in this cohort contains grade group and pathological stage only.
- **GSE116918.** Almac Prostate Disease Specific Arrays of diagnostic biopsies from men treated with radical radiotherapy. Only the Gleason score is available, so the grade benchmark uses score groups 6, 7 and 8 to 10, and there is no pathological stage or margin. This cohort entered sensitivity analysis S3 and the exploratory analysis of metastasis (22 events).
- **TCGA-PRAD.** RNA sequencing, STAR TPM values from the Genomic Data Commons, primary tumours, one aliquot per patient. Used as a reference only (for the scaling in sensitivity analysis S2, the proliferation index and the correlation between signatures) and never to evaluate a signature. Its PSA field is the most recent value, not the preoperative one, and was not used.

# 5. Scoring rules

For signature *s* in cohort *c*, the score of patient *i* is the sum over the genes found in the cohort of the published coefficient multiplied by the z-score of the gene within the cohort. The direction is that of the original article; no sign was changed and nothing was estimated from the outcome of the evaluation cohort.

- A signature was evaluated in a cohort only if at least 80% of its gene terms were found and the terms found carried at least 80% of the sum of absolute coefficients. Missing genes had a coefficient of zero and the other coefficients were not reweighted. For a gene pair both genes had to be present.
- Gene-pair models: each pair was scored 0 or 1 by the rule of the article, multiplied by its coefficient and summed.
- Equal-weight scores: the mean of the z-scores of the genes with a stated direction, each multiplied by +1 or −1. Signatures with no direction for any gene were not scored.
- Models with clinical terms (PMIDs 34002865 and 39832058): the gene part only in the primary analysis; the full formula in sensitivity analysis S9.
- Constant terms were omitted because they do not affect ranking.
- Misprints settled against the article before the lock are listed in Additional file 3 (sheet Gene terms).

# 6. Performance measures, benchmarks and pooling

- Harrell's C-index used all follow-up. Uno's C-index was truncated at 60 months. The time-dependent AUC used the cumulative/dynamic definition with inverse probability of censoring weights, at 36 and 60 months, and was not reported at a time point with fewer than 10 patients at risk or fewer than 10 events. Confidence intervals: 1000 bootstrap resamples of patients, percentile method. The hazard ratio per standard deviation came from a univariable Cox model with a Wald interval.
- Clinical model: ISUP grade group (as a linear term), pathological stage T3 or higher, log PSA and surgical margin status in GSE54460, GSE70768 and GSE70769 (465 evaluations); grade group and stage in MSKCC (104 evaluations). A variable entered the model of a cohort if it was recorded for at least half of the patients. Patients with missing clinical variables were left out (complete cases: 92 of 100 patients with 45 events in GSE54460, 86 of 92 with 43 in GSE70769, 110 of 111 with 19 in GSE70768 and 130 of 131 with 27 in MSKCC). The C-index of the clinical model alone was 0.77, 0.76, 0.80 and 0.83 in these four cohorts. The gain in C-index when the score was added was corrected for optimism by refitting both models in each of 1000 bootstrap resamples.
- Incremental value was tested by a likelihood-ratio test in each cohort; the statistics were summed over the cohorts of a signature (degrees of freedom equal to the number of cohorts).
- Random gene sets: 1000 per signature and cohort, each gene drawn from the genes outside the signature in the same decile of mean expression in that cohort, scored with the signature's own coefficients. For the pooled comparison, the random C-indices of draw *r* in each cohort were pooled with the signature's weights. The comparison was made for signatures scored as a weighted sum (208 signatures, 564 evaluations). It was not made for the three gene-pair signatures (PMIDs 32593277, 35143414 and 42553741; 5 evaluations), which are therefore not in the denominators: 124 of 208 signatures, 59.7% of 207 articles (the two models of PMID 32264916 counting as half each) and 190 of 564 evaluations (34%) exceeded the 95th percentile.
- Pooling: random-effects meta-analysis of the logit C-index, restricted maximum likelihood, Hartung–Knapp interval; with one cohort, the estimate of that cohort.
- Multiplicity: Benjamini–Hochberg at 5% across signatures, for the incremental-value test (61 of 211 signatures with *q* < 0.05) and for the empirical *P* value against random gene sets (106 of 208).
- Software: Python 3.12 with NumPy 2.3, pandas 3.0, SciPy 1.16 and lifelines 0.30. Harrell's C-index was computed with lifelines. The Cox models (hazard ratio per standard deviation, clinical models, calibration slope) were fitted with our own NumPy implementation (Newton–Raphson, Efron's method for ties), and the random-effects meta-analysis was also implemented in NumPy. Both were written on 26 October 2026, before the plan was drafted. The tests released with the code compare the Cox fit with lifelines (coefficients and log-likelihood within 1e-4 on simulated data) and the restricted maximum likelihood estimate with direct maximisation of the restricted likelihood.

\newpage

# 7. Sensitivity analyses

**Table S3. Sensitivity analyses and the three population-level results**

| | Analysis | n | Median pooled C (IQR) | Above grade, % | Above random, % |
|:---|:-------------------------------------------|-----:|:--------------------|--------:|--------:|
| | Primary analysis | 211 | 0.573 (0.535–0.618) | 1.4 | 59.7 |
| S1 | Complete gene coverage only | 203 | 0.573 (0.529–0.617) | 2.5 | 58.5 |
| S2a | Unstandardised log2 expression | 211 | 0.584 (0.539–0.624) | 1.4 | 62.6 |
| S2b | Genes scaled to TCGA-PRAD standard deviations | 202 | 0.590 (0.540–0.626) | 2.0 | 62.9 |
| S3 | Radiotherapy cohort GSE116918 added | 213 | 0.571 (0.536–0.612) | 1.4 | 59.6 |
| S4 | Without MSKCC | 208 | 0.568 (0.525–0.613) | 2.9 | 54.2 |
| S5 | GSE54460 (RNA sequencing) only | 190 | 0.549 (0.494–0.592) | 6.3 | 36.6 |
| S6 | Original endpoint BCR or recurrence-free survival | 149 | 0.581 (0.544–0.621) | 1.3 | 64.4 |
| S7 | Weights from figures, hazard ratios or the LASSO step added | 236 | 0.573 (0.528–0.616) | 1.3 | 57.5 |
| S8 | Signatures without coefficients, equal weights | 65 | 0.578 (0.526–0.622) | 0.0 | 65.6 |
| S10 | GSE54460 without six patients with repeat libraries | 211 | 0.577 (0.536–0.622) | 1.0 | 59.7 |

n: number of signatures. Above grade: articles whose signature had a higher pooled C-index than the Gleason grade group. Above random: articles whose signature exceeded the 95th percentile of expression-matched random gene sets, among the articles whose signature could be compared with random gene sets (gene-pair signatures excluded; 207 articles in the primary analysis). IQR: interquartile range.

**S9, full formula for models with clinical terms.** Two signatures were published with clinical terms. Adding those terms to the gene part changed the C-index from 0.49 to 0.58 (GSE54460), 0.61 to 0.70 (GSE70768) and 0.62 to 0.73 (GSE70769) for PMID 34002865, and from 0.47 to 0.50, 0.68 to 0.68 and 0.57 to 0.65 for PMID 39832058.

# 8. Exploratory and post hoc analyses

None of these is used for the three prespecified population-level results. All were written into the project log before they were run.

- **E1** (specified 12 November 2026). Paired bootstrap difference between the C-index of the signature and that of the grade group in the same resample, pooled for each signature.
- **E2** (specified 12 November 2026). Proliferation benchmark: the 168 genes most correlated with *PCNA* in TCGA-PRAD, chosen without use of outcome and scored with equal weights (after Venet et al., reference 12 of the article).
- **E3** (listed in the plan). Metastasis in GSE116918.
- **E4** (listed in the plan). Relation between the discrimination reported in each article and the external result. The reported values were extracted for this purpose after the primary analysis. Of the 211 signatures, 158 had a reported measure of discrimination in the training set (a 3-year AUC for 129, an AUC at another time point for 5 and a C-index for 24). The reported value was compared with the pooled external C-index. Like was also compared with like: the 24 reported C-indices with the external C-index (medians 0.73 and 0.59), and 146 reported 3-year AUCs with our AUC at 36 months in independent cohorts (medians 0.76 and 0.60).
- **V1** (post hoc, 12 November 2026). Check of the implementation: each signature was recomputed in the cohorts its article had used, and reported training performance, recomputed performance in the training cohort and performance in independent cohorts were compared. The recomputation used the measure reported by the article (the time-dependent AUC at 36 months for a 3-year AUC, 123 signatures; Harrell's C-index for a C-index, 21 signatures) and biochemical recurrence as recorded in our data, whatever the endpoint of the article (108 of the 144 articles had modelled BCR or recurrence-free survival). The comparison between cohorts the article had used and cohorts it had not (151 signatures) left TCGA-PRAD out and used a Wilcoxon signed-rank test.
- **V2** (post hoc, 12 November 2026). Precision available in each cohort. Two quantities are reported, each as 2.8 times a median standard error (80% power, two-sided 5% level): the smallest difference in C-index between a signature and the grade group that a cohort could detect, from the paired bootstrap of E1 (0.14 in GSE54460, 0.15 in GSE70769, 0.23 in GSE70768 and 0.18 in MSKCC); and the smallest C-index that a cohort could distinguish from 0.5, from the bootstrap standard error of the logit C-index (0.62, 0.63, 0.71 and 0.67). The second quantity was first drawn from the standard error of the difference with grade, which overstated it (0.64 to 0.73); this was corrected on 20 November 2026 when the additional files were prepared.
- **V3** (post hoc, 12 November 2026). Random-effects meta-regression of the logit pooled C-index on training on TCGA-PRAD only, external validation reported, number of genes (log2), final modelling method and correlation with the proliferation index (all 211 signatures; between-signature variance by the method of moments). The modelling method was not associated with performance: against LASSO-type models, the difference in logit C-index was −0.022 (95% CI −0.091 to 0.047) for multivariable Cox models and +0.142 (−0.047 to 0.331) for univariable or other models.
- **V4** (post hoc, 12 November 2026). Spearman correlation between signature scores in TCGA-PRAD, with hierarchical clustering.
- **V5** (post hoc, 18 November 2026). Calibration slope: the Cox coefficient of the published linear predictor computed on unstandardised log2 expression. The plan had stated that calibration would not be assessed. Absolute risk still cannot be assessed, and the slope depends on the expression scale of each platform; it is therefore reported by cohort. It was computed for the 564 evaluations of the 208 signatures scored as a weighted sum.
- Kaplan–Meier curves are descriptive, as stated in the plan.

# 9. Amendments to the statistical analysis plan

The plan was drafted on 1 November 2026 and finalised on 11 November 2026, before any outcome or score of any cohort had been read. It was not registered prospectively; its first public record is the Zenodo deposit made on 20 November 2026, after the analyses were complete. The first extraction lock preceded the plan, and the primary analysis of the first 223 signatures was run on 11 November 2026.

**Before outcomes were read (11 November 2026).** Uno's C-index and the time-dependent AUC were computed with our own NumPy implementation instead of scikit-survival, which was not installed on the analysis computer; the two agree to within 1e-15 on simulated data. The evaluation was split into parallel runs with fixed seeds. The plan named lifelines for the Cox models; they were fitted with our own implementation (section 6), which had been written and tested against lifelines before the plan was drafted.

**Table S4. Changes made after the primary results had been seen**

| Date (2026) | Change | Reason |
|:----------|:-------------------------------------------|:-----------------------------------|
| 11 Nov | In the unstandardised-expression analysis (S2a), 23 missing expression values (7 genes, GSE70768 and GSE70769) were replaced by the cohort mean of the gene. | The run stopped on missing values. The primary analysis is unaffected. |
| 11 Nov | The reported test of incremental value is the combined likelihood-ratio test. The version first coded (smallest *P* value over cohorts, then Benjamini–Hochberg) is kept in Additional file 6. | The smallest *P* value is too lenient when a signature has several cohorts. |
| 12 Nov | Two signatures with a pooled C-index below 0.5 (PMIDs 35143414 and 42418409) were read again against the source articles. Extraction and calculation were confirmed; the locked data were not changed. | Check requested after seeing the result. |
| 12–14 Nov | Screening of the remaining 719 records and second extraction lock (section 3). | The planned duplicate screening had not been carried out in the first stage. |
| 14–17 Nov | In the bootstrap for incremental value, the score is standardised within each resample and resamples that do not converge are skipped. | Numerical overflow for one signature with very large coefficients read from a figure. |
| 18 Nov | Two corrections in the rerun: the endpoint classification of the second lock is used for S6, and the meta-regression covariates are read for all signatures. | Outdated file paths found when the output was checked. |
| 20 Nov | The smallest C-index distinguishable from 0.5 (Fig. 6d of the article) is computed from the standard error of the C-index, not from that of its difference with grade. | Error in a post hoc descriptive quantity found when the additional files were prepared. No other result is affected. |

All analyses of the plan were rerun on the second lock. The median pooled C-index was 0.575 on the first lock (142 signatures, 384 evaluations) and 0.573 on the second (211 signatures, 569 evaluations). Both can be reproduced from the released code.

# 10. Use of artificial intelligence

An AI model (Claude, Anthropic) was used for the keyword pre-screening, the duplicate screening of titles and abstracts, the extraction and re-extraction of genes and coefficients, writing and running the analysis code, and drafting project documents, this file and the manuscript. Eligibility rules, the adjudication of disagreements and the verification against the source articles were the responsibility of the authors and are recorded with the name of the person responsible in the project log. The authors take responsibility for the content.
