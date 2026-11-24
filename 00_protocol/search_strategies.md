# 检索式（PubMed，2026-10-26 实测）

共用模块
- SIG = (signature*[tiab] OR "gene panel"[tiab] OR "risk score*"[tiab] OR "risk model*"[tiab] OR "prognostic model*"[tiab] OR riskscore[tiab] OR classifier*[tiab])
- PROG = (prognos*[tiab] OR survival[tiab] OR recurrence[tiab] OR "Kaplan-Meier"[tiab] OR "Cox"[tiab])
- OMIC = (gene*[tiab] OR transcriptom*[tiab] OR mRNA*[tiab] OR lncRNA*[tiab] OR "RNA-seq"[tiab] OR microarray[tiab] OR TCGA[tiab] OR "Gene Expression Omnibus"[tiab] OR GEO[tiab])
- NOTPT = NOT (review[pt] OR "systematic review"[pt] OR meta-analysis[pt] OR editorial[pt] OR comment[pt] OR letter[pt] OR "published erratum"[pt] OR "retracted publication"[pt] OR "case reports"[pt])

| 癌种 | 疾病模块 | 额外限定 | 2026-10-26 命中 |
|---|---|---|---|
| ccRCC | ("clear cell renal"[tiab] OR ccRCC[tiab] OR KIRC[tiab] OR "kidney renal clear"[tiab] OR "Carcinoma, Renal Cell"[mh]) | — | 1,651 |
| 膀胱癌 | ("bladder cancer"[tiab] OR "bladder carcinoma"[tiab] OR "bladder urothelial"[tiab] OR BLCA[tiab] OR "Urinary Bladder Neoplasms"[mh]) | — | 1,483 |
| 前列腺癌 | ("prostate cancer"[tiab] OR "prostate adenocarcinoma"[tiab] OR PRAD[tiab] OR "Prostatic Neoplasms"[mh]) | AND ("biochemical recurrence"[tiab] OR BCR[tiab] OR "recurrence-free"[tiab] OR "metastasis-free"[tiab] OR survival[tiab]) | 1,011 |

完整检索式 = 疾病模块 AND SIG AND PROG AND OMIC [额外限定] NOTPT

待办：Embase（Emtree 翻译）、Web of Science 检索式；检索截止日在 G1 时写定。
