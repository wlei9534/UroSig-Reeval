# UroSig-Reeval: Statistical Analysis Plan, version 2 (English translation)

**This is an English translation of the original Chinese document** `SAP_v2_定稿_2026-11-11.pdf` (6 pages), which is the authoritative version and is posted alongside this file unchanged.

- Drafted 1 November 2026; finalised 11 November 2026.
- SHA-256 of the original PDF: `dc9965d5662bfbb408ffbd68f5ef90c57991e4a9765f0a326f7f54745694f110`
- MD5 of the original PDF: `e475832f1907454c96a520c25abe3755`
- The plan was not registered prospectively. It was first made public in this Zenodo record (https://doi.org/10.5281/zenodo.24385587) on 20 November 2026. Everything decided after 11 November 2026 is listed in `Amendments_and_post_hoc_analyses.md`.
- Version 2 of this file (20 November 2026): the first version called this record an "OSF project"; only that wording in this header was corrected. The translation of the plan is unchanged.
- Paths such as `03_extraction/locked/` refer to folders of the project repository. "Wei" is the investigator who took each decision (Wei Lei). "Sub-agent" means an independent run of the AI model used for screening and extraction.

---

## 1. Aim and design

This study re-evaluates published transcriptomic prognostic signatures for prostate cancer in independent public cohorts using one uniform method. It asks three questions: how well do the signatures discriminate in external data, do they outperform routine clinical variables, and do they outperform random gene sets of the same size.

- Design: external validation of published models (Original Article). No signature is retrained.
- Main population: localised prostate cancer after radical prostatectomy.
- Primary endpoint: biochemical recurrence (BCR).
- This plan was finalised before any outcome or score result of any cohort was read. All earlier work used only the literature, gene names and the gene index of each cohort.
- This plan replaces the analysis section of protocol v1 (26 October 2026). All changes made after protocol v1 are recorded in `logs/decision_log.md`.

## 2. Source and locking of signatures

The extraction was locked on 1 November 2026: 223 analysable signatures and 1,624 gene terms, in `03_extraction/locked/` (with MD5 checksums), not to be modified afterwards.

- Search: PubMed, 26 October 2026, 1,011 records for prostate cancer; 292 articles went to full-text extraction.
- Extraction: AI pre-extraction for 81 articles; close reading by Claude for 206 articles (151 from PMC full text, 55 from PDFs downloaded manually by Wei); full text not obtained for 5 articles.
- Verification: round 1, blinded re-extraction of 287 articles by a sub-agent; round 2, checking by experts A and B; round 3, the 50 primary-analysis candidates on which the two versions differed were settled by sub-agents and expert C; all coefficients were checked verbatim against the original articles.
- In general only the final model promoted by the authors is taken from each article. The exception is PMID 32264916 (both the 22-gene DFS model and the 5-gene OS model are included, and are clustered as one article in article-level summaries).

| Analysis role | Source of weights | Signatures | Use |
|---|---|---|---|
| `main` | Coefficients printed in the text, tables or supplementary material (including gene-pair models and PMID 35414328, refitted from the authors' data) | 167 | Primary analysis |
| `equal_weight` | Gene names only; direction taken from the article | 35 (27 with a usable direction) | Equal-weight score, reported separately |
| `equal_signed` | The article gives only the direction, ±1 | 2 | Equal-weight score, reported separately |
| `sens_figure` | Coefficients read from a figure | 9 | Sensitivity analysis |
| `sens_hr` | Coefficient = ln(HR) | 9 | Sensitivity analysis |
| `sens_lasso` | LASSO-stage coefficients of a non-linear final model | 1 | Sensitivity analysis |

## 3. External cohorts and endpoints

The primary analysis uses four radical prostatectomy cohorts (140 BCR events in total). The radiotherapy cohort GSE116918 enters sensitivity analyses only. TCGA-PRAD is a reference only and is not counted in any external evaluation.

| Cohort | Population | Samples usable for BCR | BCR events | Median follow-up of censored patients (months) | Role |
|---|---|---|---|---|---|
| GSE54460 | Radical prostatectomy, FFPE, 3 centres | 100 | 49 | 73.8 | Primary analysis |
| GSE70769 | Radical prostatectomy | 92 | 45 | 79.4 | Primary analysis |
| GSE70768 | Radical prostatectomy | 111 | 19 | 33.9 | Primary analysis |
| MSKCC (GSE21032) | Radical prostatectomy | 131 | 27 | 50.4 | Primary analysis |
| GSE116918 | Radical radiotherapy (diagnostic biopsy) | 248 | 56 | 86.0 | Sensitivity analysis |
| TCGA-PRAD | Radical prostatectomy | 420 | 52 | 31.6 | Reference only (used by the original articles of 217 signatures) |

- Primary endpoint, BCR: defined by the recurrence field and follow-up time (months) provided in each cohort's public data, harmonised in `04_cohorts/processed`. The MSKCC outcome is the DFS field from cBioPortal, which is not identical to BCR; a sensitivity analysis repeats the analysis without this cohort.
- Secondary endpoint, metastasis: only GSE116918 has metastasis follow-up (248 patients, 22 events, radiotherapy population). The other cohorts have sample type only. With few events in a single cohort, the metastasis analysis is exploratory and descriptive and is not used for conclusions.
- In GSE54460, only the original library is kept for 6 patients with repeat libraries. The PSA field in TCGA contains post-operative information and must not be used as a baseline covariate. Both points are recorded in the decision log.

## 4. Evaluability rules

Under the three pre-specified rules below, 142 of the 167 primary-analysis signatures can be evaluated in at least one primary cohort (121 in two or more), giving 384 signature × cohort combinations; with GSE116918, 143 signatures and 461 combinations.

1. Leakage criterion (Wei, 1 November 2026): if the original article used a cohort in any way (training, validation, differential expression, candidate gene screening, single-gene analysis and so on), that cohort is not used to evaluate the signature. A mention only in the discussion or reference list does not count. GSE70770 is treated as use of both GSE70768 and GSE70769. Judgements are in `06_results/eligibility/leakage_matrix.csv`.
2. Gene coverage (Wei, 1 November 2026): the primary analysis requires that at least 80% of the gene terms are found in the cohort and that the found terms carry at least 80% of the sum of absolute coefficients. The sensitivity analysis requires 100%. Genes are matched by current HGNC symbol (snapshot of 26 October 2026); for a gene pair both genes must be present.
3. Direction for equal-weight scores (Wei, 1 November 2026): only genes whose direction is given in the article are used. Six articles with no direction for any gene (PMIDs 28361034, 37280589, 38129557, 38577608, 39331250, 39391236) and two articles without gene names (41003824, 39789383) are not evaluated.

Signatures with no usable cohort are not discarded. They are listed in the results with the reason (all cohorts used by the article, or insufficient gene coverage) and reported as a finding: "cannot be validated independently".

## 5. Calculation of the risk score

For each signature in each evaluable cohort one continuous risk score is calculated. Its direction is fixed by the original article (higher score, higher risk); no sign is flipped and nothing is re-estimated from the cohort's outcome.

S_i = Σ_{g ∈ G_found} β_g · z_ig

- Expression data: the processed log2 expression matrices in `04_cohorts/processed` (microarray cohorts quantile-normalised within cohort).
- Standardisation (Wei, 1 November 2026): in the primary analysis each gene is converted to a z-score within the cohort and then multiplied by its coefficient, to remove differences in gene scale between microarray and sequencing platforms. Two sensitivity analyses: unstandardised log2 expression; and, for signatures trained on TCGA, rescaling by the gene's standard deviation in TCGA-PRAD (β × SD_TCGA × z). The C-indices from the three approaches are reported side by side, and signatures whose result changes substantially with the approach are flagged.
- Missing genes: coefficient set to 0; the remaining coefficients are not re-weighted.
- Gene-pair models: each pair is scored 0/1 by the article's rule (in most, 1 if A is higher than B; in PMID 32593277 the reverse, 1 if A is lower than B), then multiplied by the coefficient and summed.
- Equal-weight score: the mean of the z-scores of the genes with a direction, each multiplied by +1 or −1.
- Models containing clinical terms (Gleason and stage terms in PMID 34002865; Gleason term in PMID 39832058): the primary analysis uses the gene part only.
- Constant terms (for example −0.009 in PMID 41758298) do not affect ranking and are omitted.
- Settled special cases: PAQR6 in PMID 36647875 is taken as 0.396; NCF1 in PMID 35222640 is positive; PMID 35414328 is matched by the Ensembl IDs given by the authors.

## 6. Performance measures

The primary measure is Harrell's C-index for each signature × cohort (BCR, continuous score). The score is never dichotomised.

| Measure | Definition | Confidence interval |
|---|---|---|
| Harrell C (primary) | Concordance between score and time to BCR, all follow-up | Patient-level bootstrap, 1,000 resamples, percentile method |
| Uno C | IPCW-weighted, truncated at 60 months | As above |
| Time-dependent AUC | At 36 and 60 months, cumulative/dynamic definition, IPCW | As above |
| HR per SD | Univariable Cox model, score scaled by its within-cohort standard deviation | Wald 95% CI |

- If at a time point a cohort has fewer than 10 patients at risk or fewer than 10 events, the AUC at that time point is not reported. The median follow-up of censored patients in GSE70768 is only 33.9 months, so its 60-month AUC may not be reported.
- Kaplan–Meier tests after a median split are not used as conclusions; if shown, they are descriptive only.
- Calibration is not assessed: most original articles do not report a baseline hazard, and the external cohorts are on a different expression scale, so absolute risk cannot be recovered.

## 7. Benchmarks

Each signature is compared with three benchmarks: clinical grade, the increment when the score is added to a clinical model, and random gene sets of the same size.

1. Clinical benchmark: the C-index of the Gleason grade group (ISUP 1–5), compared with the signature in the same patients. GSE116918 does not separate Gleason 3+4 from 4+3, so Gleason sum groups (6 / 7 / 8–10) are used there. Extended clinical model: grade group + pT stage + pre-operative PSA + surgical margin, using only the variables available in the cohort (MSKCC has no PSA or margin; GSE116918 has no pT stage or margin).
2. Incremental value: within-cohort Cox models, "extended clinical" versus "extended clinical + score"; ΔC is reported (bootstrap, 1,000 resamples, refitted in each resample with optimism correction) together with a likelihood-ratio test. Missing clinical variables are handled by complete-case analysis and the sample size is reported.
3. Random gene sets: for each signature × cohort, 1,000 random gene sets of the same size are generated (each gene drawn from the non-signature genes in the same decile of mean expression in that cohort). They use the signature's coefficient vector and the C-index is calculated in the same way. The percentile of the observed C in the random distribution is reported; above the 95th percentile counts as "better than random". Random seeds are fixed and recorded.

The weights of commercial tests (Decipher, Prolaris, Oncotype DX GPS) are not public; this study does not reproduce them and does not use them as benchmarks.

## 8. Pooling across cohorts and multiplicity

Conclusions rest on the distribution across the population of signatures, not on whether any single signature is "significant".

- Single signature: with two or more evaluable cohorts, a random-effects meta-analysis of logit(C) (REML, confidence interval with the Hartung–Knapp adjustment), reporting the pooled C, τ² and I²; with one cohort, only that cohort's C.
- Three primary population-level results:
  1. the distribution of the pooled C (median, interquartile range);
  2. the proportion of signatures whose pooled C exceeds the C of the Gleason grade group in the same cohorts;
  3. the proportion of signatures "better than random" (the random C values of draw r in each cohort are pooled with the signature's weights to give a pooled random distribution; with 1,000 draws the smallest empirical p is 0.001).
- Multiplicity: per-signature tests (likelihood-ratio test of incremental value; empirical p against random sets) are controlled at a false discovery rate of 5% by the Benjamini–Hochberg method.
- Clustering by article: the two models of PMID 32264916 count as one article in proportion-type results.
- Exploratory analyses (not for conclusions): relation of the pooled C with the number of genes, year of publication, type of original endpoint, and the internal C/AUC reported by the article.

## 9. Sensitivity analyses

Each item below recomputes the three population-level results of section 8, reported alongside the primary analysis.

| # | Change | Purpose |
|---|---|---|
| S1 | Gene coverage required to be 100% | Whether missing genes lower performance |
| S2 | Score computed from unstandardised log2 expression; and, for signatures trained on TCGA, rescaled by TCGA standard deviations | Effect of the standardisation method |
| S3 | Add GSE116918 (radiotherapy population) | Generalisation across population and treatment |
| S4 | Remove MSKCC | Its outcome is DFS, not strictly BCR |
| S5 | Use only GSE54460, the RNA-seq cohort | Microarray versus sequencing platform |
| S6 | Keep only signatures whose original endpoint was BCR/RFS | Performance when the original modelling endpoint matches this study |
| S7 | Add the `sens_figure`, `sens_hr` and `sens_lasso` signatures | Signatures with lower-precision coefficients |
| S8 | `equal_weight` and `equal_signed` signatures pooled separately | Performance of signatures without coefficients |
| S9 | Models with clinical terms computed with the full formula (genes + clinical) | Comparison with the gene part alone |
| S10 | GSE54460 without the 6 patients with repeat libraries | Effect of inconsistent follow-up time |

## 10. Reporting, software and reproducibility

The plan becomes final once confirmed by Wei: a PDF is exported, its checksum and date are entered in the decision log and archived by e-mail; any later change is recorded as a deviation.

- Reporting standard: external validation of prediction models reported according to TRIPOD; the literature source is given as a PRISMA-style flow diagram (1,011 records → 292 full-text extractions → 223 analysable → 142 evaluable in the primary analysis), with reasons and numbers for exclusion at each step.
- Software: Python (pandas, lifelines; time-dependent AUC and Uno C with scikit-survival), versions written into the result files; code in `05_code`, all run on Wei's own computer.
- Randomness: seeds for the bootstrap and the random gene sets are fixed and recorded.
- Deviations: all changes are written to `logs/decision_log.md`, stating whether they were made after results had been seen.
- Data and code sharing: the locked signature library, the evaluability tables and the analysis code are released with the article; all cohort data are public.

Open items: none. To be finalised once Wei has read and confirmed it.
