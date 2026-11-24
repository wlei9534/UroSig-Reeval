# Data dictionary

Column headers and category labels that appear in Chinese in the tables.

## Column headers

| Header | English | First seen in |
|---|---|---|
| AI 理由 | reason given by the keyword pre-screen | `02_screening/screening_decisions.csv` |
| AI 预筛 | label of the keyword pre-screen | `02_screening/screening_decisions.csv` |
| A_代码 | reason code of screener A | `02_screening/screening_decisions.csv` |
| A_决定 | decision of screener A (AI run) | `02_screening/screening_decisions.csv` |
| A_理由 | reason of screener A | `02_screening/screening_decisions.csv` |
| B_代码 | reason code of screener B | `02_screening/screening_decisions.csv` |
| B_决定 | decision of screener B (AI run) | `02_screening/screening_decisions.csv` |
| B_理由 | reason of screener B | `02_screening/screening_decisions.csv` |
| Ensembl 来源 | source of the Ensembl ID | `03_extraction/locked/signature_terms_locked_2026-11-01.csv` |
| GSE116918_用途 | use of GSE116918 by the article | `03_extraction/v2/leakage_v2_scan.csv` |
| GSE54460_用途 | use of GSE54460 by the article | `03_extraction/v2/leakage_v2_scan.csv` |
| GSE70768_用途 | use of GSE70768 by the article | `03_extraction/v2/leakage_v2_scan.csv` |
| GSE70769_用途 | use of GSE70769 by the article | `03_extraction/v2/leakage_v2_scan.csv` |
| HGNC 现行符号 | current HGNC symbol | `03_extraction/locked/signature_terms_locked_2026-11-01.csv` |
| MSKCC_用途 | use of MSKCC by the article | `03_extraction/v2/leakage_v2_scan.csv` |
| TCGA_用途 | use of TCGA-PRAD by the article | `03_extraction/v2/leakage_v2_scan.csv` |
| Wei 决定 | decision by W.L. | `03_extraction/locked/signatures_locked_2026-11-01.csv` |
| X 训练 | training cohort (extractor X) | `03_extraction/v2/leak_check.csv` |
| X 验证 | validation cohort (extractor X) | `03_extraction/v2/leak_check.csv` |
| X_类别 | category (extractor X) | `03_extraction/v2/compare_v2.csv` |
| Y 训练 | training cohort (extractor Y) | `03_extraction/v2/leak_check.csv` |
| Y 验证 | validation cohort (extractor Y) | `03_extraction/v2/leak_check.csv` |
| Y_类别 | category (extractor Y) | `03_extraction/v2/compare_v2.csv` |
| signature 数 | number of signatures | `06_results/eligibility/eligibility_summary.csv` |
| ≥1 个可评估外部队列 | at least one evaluable external cohort | `06_results/eligibility/eligibility_summary.csv` |
| ≥2 个 | at least two evaluable external cohorts | `06_results/eligibility/eligibility_summary.csv` |
| 一致性 | agreement | `03_extraction/v2/compare_v2.csv` |
| 仅提取员记录 | recorded by the extractors only | `03_extraction/v2/leak_check.csv` |
| 仅正文扫描命中 | found by the full-text scan only | `03_extraction/v2/leak_check.csv` |
| 依据 | basis | `03_extraction/v2/adjudication_v2.csv` |
| 分析角色 | analysis role | `03_extraction/locked/signatures_locked_2026-11-01.csv` |
| 分歧类型 | type of disagreement | `03_extraction/v2/adjudication_v2.csv` |
| 判定来源 | source of the judgement | `03_extraction/v2/leakage_v2_scan.csv` |
| 原文用过 | cohort used by the article | `06_results/eligibility/eligibility_signature_by_cohort.csv` |
| 原文终点（复核员记录） | endpoint of the original model (recorded by the checker) | `06_results/eligibility/original_endpoint_class.csv` |
| 可评估 signature×队列 | evaluable signature x cohort combinations | `06_results/eligibility/eligibility_summary.csv` |
| 可评估_0.8 | evaluable at 80% coverage | `06_results/eligibility/eligibility_signature_by_cohort.csv` |
| 可评估_1.0 | evaluable at 100% coverage | `06_results/eligibility/eligibility_signature_by_cohort.csv` |
| 基因 | gene | `03_extraction/manual_extraction_pdf53_verify.csv` |
| 基因与系数相邻(≤200字符) | gene and coefficient within 200 characters in the text | `03_extraction/manual_extraction_pdf53_verify.csv` |
| 基因在原文 | gene found in the article text | `03_extraction/manual_extraction_pdf53_verify.csv` |
| 基因数 | number of genes | `03_extraction/locked/signatures_locked_2026-11-01.csv` |
| 基因相同 | genes identical | `03_extraction/v2/compare_v2.csv` |
| 基因覆盖率 | gene coverage | `06_results/coverage/gene_coverage_by_signature.csv` |
| 基因项数 | number of gene terms | `06_results/coverage/gene_coverage_by_signature.csv` |
| 基因（原文写法） | genes as written in the article | `03_extraction/locked/signatures_locked_2026-11-01.csv` |
| 备注 | note | `02_screening/screening_decisions.csv` |
| 年份 | year of publication | `03_extraction/locked/signatures_locked_2026-11-01.csv` |
| 序号 | serial number | `03_extraction/locked/signature_terms_locked_2026-11-01.csv` |
| 建议用途 | suggested use | `03_extraction/locked/signatures_locked_2026-11-01.csv` |
| 批次 | batch | `03_extraction/v2/adjudication_v2.csv` |
| 找到 | terms found in the cohort | `06_results/coverage/gene_coverage_by_signature.csv` |
| 拟定角色 | proposed analysis role | `03_extraction/v2/reconciled_v2.csv` |
| 提取状态 | extraction status | `06_results/evaluation/exploratory/internal_reported.csv` |
| 数值来源 | source of the values | `03_extraction/manual_extraction_pdf53_verify.csv` |
| 方向 | direction | `06_results/eligibility/equal_weight_directions.csv` |
| 无可评估队列 | no evaluable cohort | `06_results/eligibility/eligibility_summary.csv` |
| 最终状态 | final extraction status | `03_extraction/locked/signatures_locked_2026-11-01.csv` |
| 有可评分基因 | has scorable genes | `06_results/eligibility/eligibility_signature_by_cohort.csv` |
| 期刊 | journal | `03_extraction/locked/signatures_locked_2026-11-01.csv` |
| 权重覆盖率 | coverage of the sum of absolute coefficients | `06_results/coverage/gene_coverage_by_signature.csv` |
| 来源 | source | `06_results/eligibility/equal_weight_directions_v2.csv` |
| 来源版本 | extraction stage | `03_extraction/locked/signatures_locked_v2_2026-11-14.csv` |
| 标题 | title | `03_extraction/locked/signatures_locked_2026-11-01.csv` |
| 检索记录 PMID | PMID of the search record | `03_extraction/locked/signatures_locked_v2_2026-11-14.csv` |
| 状态 | status | `03_extraction/manual_extraction_pdf53_verify.csv` |
| 符号映射说明 | note on symbol mapping | `03_extraction/locked/signature_terms_locked_2026-11-01.csv` |
| 第三轮核定 | third-round check | `03_extraction/locked/signatures_locked_2026-11-01.csv` |
| 筛选阶段 | screening stage | `02_screening/screening_decisions.csv` |
| 系数 | coefficients | `03_extraction/locked/signatures_locked_2026-11-01.csv` |
| 系数在原文 | coefficient found in the article text | `03_extraction/manual_extraction_pdf53_verify.csv` |
| 系数相同 | coefficients identical | `03_extraction/v2/compare_v2.csv` |
| 系数覆盖 | coefficient coverage | `03_extraction/v2/adjudication_v2.csv` |
| 终点 | endpoint of the original model | `03_extraction/locked/signatures_locked_v2_2026-11-14.csv` |
| 终点类别 | endpoint class | `06_results/eligibility/original_endpoint_class.csv` |
| 终裁 | final decision | `02_screening/screening_decisions.csv` |
| 结果 | result | `03_extraction/v2/v2_fulltext_set.csv` |
| 编号 | record number | `02_screening/screening_decisions.csv` |
| 缺失基因 | missing genes | `06_results/coverage/gene_coverage_by_signature.csv` |
| 裁定 | adjudication | `03_extraction/v2/adjudication_v2.csv` |
| 裁定人 | adjudicator | `03_extraction/v2/adjudication_v2.csv` |
| 裁定依据 | basis of the adjudication | `03_extraction/v2/reconciled_v2.csv` |
| 规则标记 | rule flag | `03_extraction/v2/reconciled_v2.csv` |
| 角色 | role | `03_extraction/v2/adjudication_v2.csv` |
| 训练队列 | training cohort | `03_extraction/locked/signatures_locked_v2_2026-11-14.csv` |
| 采用 | adopted value | `03_extraction/v2/adjudication_v2.csv` |
| 锁定备注 | note at locking | `03_extraction/locked/signatures_locked_2026-11-01.csv` |
| 门槛 | threshold | `06_results/eligibility/eligibility_summary.csv` |
| 需人工 | needs manual review | `03_extraction/v2/adjudication_v2.csv` |
| 验证队列 | validation cohort | `03_extraction/locked/signatures_locked_v2_2026-11-14.csv` |

## Category labels

| Label | English |
|---|---|
| 用过 | used by the article |
| 未用 | not used |
| 多因素 Cox | multivariable Cox |
| LASSO 类 | LASSO-type |
| 单因素或其他 | univariable or other |
| 未映射 | not mapped |
| 已提取 | extracted |
| 第一阶段（关键词规则） | stage 1 (keyword rule) |
| 第二阶段（双份独立筛选） | stage 2 (independent duplicate screening) |
| v1（关键词优先 292 篇） | stage 1 (292 keyword-prioritised records) |
| v2（双人筛选补充） | stage 2 (duplicate screening of the remaining records) |
