# 决策日志

> 公开版本说明：本文件收录与研究设计、数据、分析和结果报告有关的决定与偏离。与分析无关的投稿事务（题目与格式、投稿信、审稿人建议等）的记录未收入本文件。


| 日期 | 决定 | 原因 | 决定人 |
|---|---|---|---|
| 2026-10-26 | 选题定为泌尿肿瘤已发表预后 signature 系统性再评估 | 2026 年 PubMed 全景中 signature 文章饱和、跨研究外部再评估稀缺 | Wei |
| 2026-10-26 | 暂按 ccRCC 主分析、膀胱癌/前列腺癌重复臂推进 | 方案默认；在线文档中留有评论待确认 | 默认，待 Wei 确认 |
| 2026-10-26 | 分析代码用 Python 实现 | 当前计算环境无 R；Python 生存分析栈可覆盖 C 指数、Cox、bootstrap、随机效应汇总 | Claude 执行 |
| 2026-10-26 | PubMed 检索执行：ccRCC 1,651、膀胱癌 1,483、前列腺癌 1,011，合计 4,145 条，按 PMID 去重 4,102 条 | 原始导出 01_search/raw/pubmed_2026-10-26.json.gz（SHA-256 b99e8a32…cf0b） | Claude 执行 |
| 2026-10-26 | 筛选表加入关键词规则"AI 预筛"列（可能纳入 1,600 / 不确定 2,181 / 可能排除 321），只作排序提示 | 加快双人筛选；不替代人工判断 | Claude 执行 |
| 2026-10-26 | 外部队列核验：ccRCC 公开队列 OS 事件数很少（GSE29609 17、GSE40912 16；E-MTAB-1980 生存数据不在库中；CPTAC-3 RCC 33 例死亡），膀胱癌与前列腺癌队列明显更充足 | 见 04_cohorts/cohort_registry.csv | 需 Wei 决定主分析癌种 |
| 2026-10-26 | 研究范围改为只做前列腺癌；ccRCC 与膀胱癌不纳入 | ccRCC 公开外部队列事件数不足；Wei 选择前列腺癌单癌种 | Wei |
| 2026-10-26 | 主分析人群 = 根治术后局限期前列腺癌，主要终点 BCR；放疗队列（GSE116918）作敏感性分析，转移为次要终点 | 方案更新 | Wei 决定范围，Claude 执行 |
| 2026-10-26 | 文章类型按 Original Article 投稿 | 核心产出是在独立个体数据上的新分析，非文献结果汇总 | Wei 认可 |
| 2026-10-26 | 筛选改用前列腺癌专用表（1,011 条）；三癌种总表移入 02_screening/archive_3cancers | 范围变更 | Claude 执行 |
| 2026-10-26 | 数据获取改为在 Wei 本机直接下载（Desktop Commander，curl），预处理在本机 Python 运行 | 云端环境不能直连 NCBI，浏览器中转单次上限约 260 KB | Claude 执行 |
| 2026-10-26 | 6 个队列完成统一预处理（04_cohorts/processed）；主分析 = GSE54460、GSE70769、MSKCC、GSE70768，合计 140 例事件 | MSKCC 经 cBioPortal 补到结局（DFS，131 例 27 事件） | Claude 执行 |
| 2026-10-26 | 偏离方案：Illumina / Almac / exon 芯片使用 GEO series matrix 的提交者处理值，再在队列内做分位数标准化，而非从原始 CEL/idat 重做 | 原始文件不全或体量大（GSE21032 原始数据 24 GB）；作为已知偏离记录，必要时以敏感性分析补做原始数据处理 | Claude 执行，待 Wei 确认 |
| 2026-10-26 | GSE54460 的 6 个重复文库在 GEO 中随访时间不一致：保留原始文库（标题无 -2），加 flag 做剔除敏感性分析 | 数据源不一致 | Claude 执行 |
| 2026-10-26 | TCGA-PRAD 的 PSA 字段为最近一次 PSA（与 BCR 的 C=0.88），禁止作为基线协变量 | 含术后信息，属泄漏 | Claude 执行 |
| 2026-10-26 | HGNC 符号快照固定为 2026-10-26 下载版本（04_cohorts/raw/reference） | 可重复性 | Claude 执行 |
| 2026-10-26 | 发现并修正检索导出缺陷：PubMed XML 解析时 DOI/PMCID 误取了参考文献列表中的标识符（4,102 条中 DOI 错 2,520 条、PMCID 错 3,264 条；PMID、标题、摘要、期刊不受影响）。改用 esummary 重取每篇文章自身的 DOI/PMCID，重建两份筛选表（无人开始填写前替换） | 修正文件：01_search/pubmed_dedup_2026-10-26_idfixed.csv、01_search/raw/pubmed_2026-10-26_idfixed.json.gz；原始导出保留作溯源 | Claude 执行 |
| 2026-10-26 | 错误 PMCID 抓到的 239 篇全文移入 03_extraction/_to_delete/（未删除） | 与目标文献不符 | Claude 执行，待 Wei 删除 |
| 2026-10-26 | 生成 AI 预提取表 03_extraction/ai_prefill_PRAD.xlsx：292 篇候选中 81 篇自动提取到基因与系数（79 篇逐字核对通过），142 篇未找到、60 篇无开放获取全文、9 篇系数在补充材料 | 第三方核对用，不替代双人提取；已排除引自他人研究的公式、识别基因对与单/多因素混合表 | Claude 执行 |
| 2026-10-31 | 本机全文下载完成并逐个核对大小：PMC 开放获取包内 3,679 个文件（XML/PDF/补充材料/图片，含 145 个 TIFF 补充图）全部到位；另有 11 篇在 PMC 只有题录、无开放获取包 | 03_extraction/fulltext/<PMCID>/；清单 fulltext_inventory.csv、tiff_expected.csv | Claude 执行 |
| 2026-10-31 | 无开放获取全文的 60 篇（49 篇 PMC 无全文 + 11 篇 PMC 仅题录）导出为 03_extraction/无全文文献清单_60篇.txt（PMID、标题、DOI；4 篇无 DOI 给 PubMed 链接），DOI 与 esummary 记录逐条比对一致 | 供 Wei 自行获取全文，放入 03_extraction/fulltext_manual/<PMID>/ 后再提取 | Claude 执行 |
| 2026-10-31 | 151 篇（142 篇“未找到系数”+ 9 篇“系数可能在补充材料”）逐篇精读全文、补充材料与图片完成：已提取 46、基因对 2、仅方向 1、图片读数 9、HR 换算 10、仅基因 40、非线性模型 11（另 1 篇附 LASSO 系数）、图片过模糊 2、其余 29 篇按理由排除 | 03_extraction/manual_extraction.jsonl；合并进 ai_prefill_PRAD.xlsx 的“人工精读提取”等页（旧版存 03_extraction/archive/） | Claude 执行，待第二人独立核对 |
| 2026-10-31 | 逐字核对：取自正文/表格/补充材料的 695 个基因-系数项，系数全部在原文找到；例外均已人工确认（蛋白名→基因符号映射、MathML 排版、原文笔误 CYP4FB） | 05_code/scripts/verify_manual.py → 03_extraction/manual_extraction_verify.csv | Claude 执行 |
| 2026-10-31 | 系数来源分级：正文/表格/补充材料 → 主分析候选；图片读数与 HR 取对数换算 → 仅作敏感性分析；无系数者可做等权评分敏感性分析或联系作者 | 图片读数与 HR 两位小数换算精度不足 | Claude 建议，待 Wei 决定 |
| 2026-10-31 | 42553741 的 36 个基因对系数在出版社网站补充材料（Dove Press），已下载至 03_extraction/fulltext_manual/supp_42553741/ | PMC 开放获取包未收录该补充文件 | Claude 执行 |
| 2026-10-31 | 疑似重复发表：39518911 与 39399052 为同一 7 基因模型，纳入时只计一次（取系数精度更高的 39399052） | 避免同一 signature 重复计数 | 待 Wei 确认 |
| 2026-10-31 | 题录阶段误纳入、全文确认排除：38924839（食管鳞癌）、35983076（肾癌/泛癌）、39511540（泛癌单基因） | 非前列腺癌 | Claude 执行，待双人确认 |
| 2026-10-31 | 无全文 60 篇补充提取：Wei 提供 53 个 PDF，52 篇逐篇读完（已提取 24、图片读数 1、基因对 1、仅基因 6、系数在补充材料待下载 2、非线性 7、其余 11 篇按理由排除）；38070051 的系数从 Springer 补充表 Table S3 取得 | 03_extraction/manual_extraction_pdf53.jsonl；合并进 ai_prefill_PRAD.xlsx（旧版存 03_extraction/archive/ai_prefill_PRAD_2026-10-31_before_pdf53.xlsx） | Claude 执行，待第二人独立核对 |
| 2026-10-31 | 60 篇中未提取 8 篇：41568495、41111380（找不到全文）、40916526、37846114、33504416（中文）按 Wei 决定放弃；42361668、32882325 不在提供文件中（非中文）；32898860 的 PDF 被截断无法打开 | 逐篇“当前状态”已标注 | Wei 决定 / 后 3 篇待 Wei 确认或重新下载 |
| 2026-10-31 | 逐字核对（PDF 文本 + 已下载补充材料）：262 项，取自正文/表格/补充材料的系数除 4 项外均在原文找到、另 1 项不相邻；例外已人工确认（科学计数法 40295064、排版错位 39760809 MATN3、图片中的基因名 35414328 LINC02562、跨页公式 37193176） | 05_code/scripts/verify_pdf53.py → 03_extraction/manual_extraction_pdf53_verify.csv | Claude 执行 |
| 2026-10-31 | 原文疑似笔误，按原文记录并标出待判：36647875 PAQR6 系数印作 396（疑为 0.396）；39760809 “MATN * 30.449” 判读为 MATN3 0.449 | 不擅自改数；见“待人工核对”页 | 待双人确认 |
| 2026-10-31 | 补齐 60 篇：Wei 重新下载 32898860 并补交 42361668、32882325（两篇非中文，原先漏放）；三篇读完：32882325、42361668 已提取，32898860 为基因组选择方法学论文、GENE160 在每个队列内重拟合 BLUP 权重，无固定模型，排除 | 03_extraction/manual_extraction_pdf53.jsonl（现 55 条） | Claude 执行，待第二人独立核对 |
| 2026-10-31 | 39252915（medRxiv 预印本 Meta-Score）排除 | 未经同行评议，45 基因列表取不到 | Wei |
| 2026-10-31 | 33338276：Wei 下载 Wiley 补充表 Table S3，提取 CliSigGlu-22（主记录）与 GleSigGlu-26（写在备注）；补充表两列表头疑似标反，按基因个数（22/26）和正文“重叠 17 个”判定归属 | fulltext_manual/supp_33338276；见“待人工核对” | Claude 执行，待双人确认 |
| 2026-10-31 | 35414328：补充材料只有训练集表达矩阵、无系数；用作者数据重拟合多因素 Cox 得到 17 个系数，ENSG→名称经单因素系数逐一匹配（差 <5e-5），风险评分中位数 −2.5775 与原文 Fig 3c 阈值一致；表达匹配必须用 Ensembl ID（原文的 UBXN10-AS1 = ENSG00000225986，不是 HGNC 旧名所指的 PLA2G2C） | 05_code/scripts/refit_35414328.py | Claude 执行，建议作主分析候选，待 Wei 决定 |
| 2026-10-31 | 队列检索正则修正：“Taylor & Francis”（出版社名）不再被误判为 MSKCC/Taylor 队列 | merge_pdf53.py | Claude 执行 |
| 2026-10-31 | 三处疑点经专人核对确认：33338276 补充表表头标反（基因归属按个数判定成立）；36647875 PAQR6 系数 396 为 0.396 的笔误，使用时按 0.396；39760809 “MATN * 30.449” 即 MATN3 0.449 | 记入 ai_prefill_PRAD.xlsx“待人工核对”页“核对结果”列 | Wei 安排核对 |
| 2026-10-31 | 系数来源分级定稿：正文/表格/补充材料（含 35414328 重拟合系数）进主分析；图片读数、HR 换算仅作敏感性分析；只有基因名的做等权评分 | 按 Claude 建议 | Wei |
| 2026-10-31 | 35414328 进入主分析；39399052 为 39518911 的预印本，删除 39399052、保留 39518911；38924839、35983076、39511540 确认排除 | apply_decisions_20261031.py 写入工作簿 | Wei |
| 2026-10-31 | 用出版社 PDF/高清图复核两篇“图太模糊”的文章：36550819 正文本就有 CPRS 公式（先前精读漏看，属提取错误，已更正为“已提取”，6 基因系数逐字见于原文）；34900992 由高清 Figure 1C 读出 11 个 lncRNA 的多因素 HR，系数 = ln(HR)，仅作敏感性分析（OS，10 个事件） | 文件移至 03_extraction/fulltext_manual/supp_36550819、supp_34900992 | Claude 执行 |
| 2026-10-31 | 工作簿改为可重建流程：verify_manual → merge_manual_extraction（以 archive/ai_prefill_PRAD_2026-10-26.xlsx 为输入）→ verify_pdf53 → merge_pdf53 → apply_decisions_20261031 | 05_code/scripts | Claude 执行 |
| 2026-11-01 | 盲法复核：287 篇（AI 预提取 81 + 精读 206）由 9 个子代理只看原文与补充材料独立重提（不看已有结果）；系数只印在图上的 23 篇另给 PMC 原始图片补读 | 03_extraction/review/（BRIEF、批次、out_*.jsonl） | Claude 执行（Wei 要求） |
| 2026-11-01 | 复核结果：AI 预提取 81 篇中 64 篇一致、17 篇不一致，17 篇对照原文均为 AI 预提取错（漏基因 8、丢负号 1、基因名截断或粘上标 4、照抄笔误 1、本应排除 3） | adjudication_AI.json | Claude 裁定，待专家确认 |
| 2026-11-01 | 复核结果：精读 206 篇中不一致 55 篇，对照原文后：提取者1（Claude）错 15 篇（9 篇正文有公式却判为“无系数”、1 篇误取 LGG 模型、3 篇取错模型或误用 HR、1 篇方向记法、1 篇补基因名），已更正；其余为无实质差异、预期差异、分级规则差异或待专家判断（42190028、32264916、40463716、33336591） | adjudication_M.json；manual_extraction*.jsonl（原记录 *.before_review.jsonl）；apply_review_corrections.py | Claude 裁定，待专家确认 |
| 2026-11-01 | 更正后逐字核对：所有取自正文的更正系数均在原文找到 | verify_manual.py / verify_pdf53.py | Claude 执行 |
| 2026-11-01 | 四方对比表 03_extraction/review/四方对比表_2026-11-01.xlsx：提取者1原值、子代理、Claude 裁定与复核后值，专家 A/B 列待填；工作簿同步加入“复核对比”页与复核后列 | 复核后：主分析候选 166、仅基因（等权评分）35、敏感性分析 22 | Claude 执行 |
| 2026-11-01 | 采纳盲法复核：AI 预提取 17 篇按子代理数据；精读 15 篇按更正后数据；其余不一致按各文主推模型取值；专家意见与复核一致 | ai_prefill_PRAD.xlsx “复核后”列 | Wei |
| 2026-11-01 | 42190028 维持排除（12 基因探索性评分不算）；40463716 排除；33336591 的 Figure 1B 数值算正文数值，进主分析 | apply_decisions_20261101.py | Wei |
| 2026-11-01 | 32264916 两个模型都纳入：DFS 22 基因为主记录，OS 5 基因为附加 signature（32264916_OS）；文献层面汇总时按同一篇聚类，避免重复计数 | 03_extraction/extra_signatures.jsonl | Wei 决定，Claude 定处理方式 |
| 2026-11-01 | 主分析候选定为 166 篇（167 个 signature）；对其中两版提取有差异的 49 篇及新增 OS 模型共 50 条做第三轮复核：2 个子代理对照原文核定当前采用值，Wei 另请专家C 复核 | review/BRIEF_R3.md、batch_R3a/b.csv；第三轮复核表_主分析候选_2026-11-01.xlsx | Wei 要求，Claude 执行 |
| 2026-11-01 | 第三轮复核（2 个子代理，50 条）：采用值正确 40、需修改 9（8 篇改用原文最高精度系数；35669428 基因名 BAK→BAK1）、无法确定 1（35222640 NCF1 系数前原文缺正负号，暂按 + 记录，分析时做正负两种符号的敏感性分析） | review/out_R3a.jsonl、out_R3b.jsonl；apply_round3.py | Claude 执行，待专家C 与 Wei 定稿 |
| 2026-11-01 | 第三轮最终系数逐个回原文核对，均逐字找到（35414328 为作者数据重拟合，子代理独立复算一致、中位数 −2.5775 与原文阈值吻合；39760809 为已确认笔误） | Claude 核对 | Claude 执行 |
| 2026-11-01 | 第三轮复核定稿：专家C 与第三轮结果一致；35222640 的 NCF1 系数确认为正号，公式笔误 COLA1=COL5A1、PLXNDI=PLXND1 确认 | 专家C 核对 | Wei |
| 2026-11-01 | 提取锁定：223 个可分析 signature（主分析 167〔166 篇 + 32264916_OS〕、等权评分 35、方向等权 2、图片读数 9、HR 换算 9、LASSO 阶段 1），1,624 个基因项；36647875 PAQR6 按 0.396；35414328 按作者 Ensembl ID 对应基因 | 03_extraction/locked/（signatures_locked_2026-11-01.csv、signature_terms_locked_2026-11-01.csv、MANIFEST.txt）；lock_extraction.py | Wei 定稿，Claude 执行 |
| 2026-11-01 | 基因覆盖：按锁定的 1,624 个基因项（HGNC 现行符号）在 6 个处理后队列的基因索引中查找，只读基因名、不读结局 | 06_results/coverage/；gene_coverage.py（在 Wei 本机运行） | Claude 执行 |
| 2026-11-01 | 泄漏矩阵：正文扫描（去参考文献，含 Stockholm/Cambridge/Taylor 别名）与盲法复核员记录一致者直接采用（202 个）；不一致的 21 个由子代理回原文核定（训练/验证/其他=用过，仅讨论或引文提及=未用） | 06_results/eligibility/leakage_matrix.csv；review/out_LEAK.jsonl | Claude 执行 |
| 2026-11-01 | 等权评分基因方向：33 个仅基因 signature 由子代理从原文提取每个基因的风险方向（+1/−1/0）；6 篇全部基因方向不明（28361034、37280589、38129557、38577608、39331250、39391236） | 06_results/eligibility/equal_weight_directions.csv；review/out_DIR.jsonl | Claude 执行 |
| 2026-11-01 | 可评估表（草案，待 Wei 定门槛）：覆盖≥80% 时主分析 167 个中 143 个至少有 1 个可评估外部队列（461 个 signature×队列），24 个无；覆盖 100% 时 141 个（402 个） | 06_results/eligibility/eligibility_summary.csv | 待 Wei 决定 |
| 2026-11-01 | 泄漏标准（预先规定）：原文以任何方式用过某外部队列（训练、验证，或差异表达、候选基因筛选、单基因分析等其他分析）即不用该队列评估该 signature；仅在讨论或参考文献中提及不算 | 06_results/eligibility/leakage_matrix.csv | Wei |
| 2026-11-01 | 基因覆盖门槛（预先规定）：主分析要求基因覆盖 ≥80% 且系数绝对值权重覆盖 ≥80%；敏感性分析要求 100% 覆盖。按此主分析 167 个 signature 中 143 个可评估（461 个 signature×队列） | eligibility_signature_by_cohort.csv（可评估_0.8 / 可评估_1.0） | Wei |
| 2026-11-01 | 等权评分：只用原文给出方向的基因（+1/−1），无方向基因剔除；全部基因无方向的 6 篇（28361034、37280589、38129557、38577608、39331250、39391236）不评估，流程图注明 | equal_weight_directions.csv | Wei |
| 2026-11-01 | 统计分析计划 SAP v2 起草（Claude Docs：UroSig-Reeval 统计分析计划），取代方案 v1 的分析部分；定稿方式：导出 PDF、校验码与日期记入本日志、邮件存档（不做 OSF 注册） | https://claude.ai/artifact/NjvjhZA5xBYUpWUH3b75XC | Wei |
| 2026-11-01 | 评分标准化（预先规定）：主分析每个基因在队列内转 z 分数再乘系数；敏感性分析：未标准化 log2 表达；原文用 TCGA 训练的 signature 按 TCGA-PRAD 标准差还原尺度（β×SD_TCGA×z）；三种做法并列报告 | SAP 第 5、9 节 | Wei |
| 2026-11-01 | 主分析队列只含 4 个根治术队列：142 个主分析 signature 可评估（121 个 ≥2 队列），384 个 signature×队列；此前报告的 143 个 / 461 个含 GSE116918 | eligibility_signature_by_cohort.csv | Claude 更正 |
| 2026-11-01 | 转移终点只在 GSE116918 有随访（248 例、22 个事件）；其余队列只有样本类型字段。转移分析只作探索性描述 | 04_cohorts/processed/*/clinical.csv | Claude 核对 |
| 2026-11-11 | 统计分析计划 SAP v2 定稿（Wei 确认）：导出 PDF 00_protocol/SAP_v2_定稿_2026-11-11.pdf，SHA-256 dc9965d5662bfbb408ffbd68f5ef90c57991e4a9765f0a326f7f54745694f110（MD5 e475832f1907454c96a520c25abe3755），文档版本 rev 20；此后改动均作为偏离记录 | 00_protocol/；待 Wei 邮件存档 | Wei |
| 2026-11-11 | 实现说明（定稿后、读结局前）：① Uno C 与时间依赖 AUC 用自写 numpy 实现（05_code/urosig/timeauc.py），本机无 scikit-survival，已在模拟数据上与 scikit-survival 0.28 比对、最大差 1e-15（05_code/tests/test_timeauc.py）——属 SAP 第 10 节“用 scikit-survival”的实现偏离，结果等价；② 为在本机 16 核并行，run_evaluation.py 增加按组合分片（--shard），随机种子按 signature×队列固定，分片与否结果一致（模拟数据核对最大差 4e-16）；③ S6 的“原文终点为 BCR/RFS 类”按复核员记录的终点文字、取最先出现的终点词判定（规则见 classify_endpoint.py），主分析 167 个中 117 个归 BCR/RFS 类、50 个为 OS/PFS/PFI/MFS/未明确 | 05_code/scripts/classify_endpoint.py、merge_shards.py；06_results/eligibility/original_endpoint_class.csv | Claude 执行 |
| 2026-11-11 | 主分析完成（384 个组合全部成功，用时 33 分钟）。随后 S2a（未标准化 log2）在随机基因集上因 GSE70769/GSE70768 共 23 个缺失表达值（7 个基因）报错中断；修正：raw 模式下缺失值用该基因队列均值填补（与主分析 z 分数下缺失记 0 等价），S2–S10 重新运行。主分析不受影响 | sigscore.py、run_evaluation.py | Claude 执行（偏离：实现修正） |
| 2026-11-11 | 汇总程序修正（看到主分析汇总数字后发现，透明记录）：原代码把每个 signature 各队列似然比检验的最小 p 值拿去做 BH，这在多队列时偏宽松；改为各队列 LR χ² 相加、自由度 = 队列数的合并检验（各队列病人互不重叠），作为“增量价值检验”的报告值；原最小 p 版本一并保留在结果表中。另补充描述性汇总（Uno C、AUC36、HR/SD 中位数，CI 不含 0.5 的个数） | pool_results.py | Claude 执行（偏离：统计实现修正，待 Wei 知悉） |
| 2026-11-11 | 主分析与敏感性分析 S1–S10 全部完成（Wei 本机，14 进程并行，约 2.5 小时；所有分片无报错）；全部汇总用同一版 pool_results.py 统一重跑；结果表 06_results/UroSig_结果表_主分析与敏感性分析_2026-11-11.xlsx（build_report_tables.py） | 06_results/evaluation/ | Claude 执行 |
| 2026-11-11 | 结果核对中注意到（未改动任何锁定数据，仅记录，待 Wei 决定是否回原文复查）：35143414（基因对，规则 A>B 记 1）与 42418409 汇总 C 的 95%CI 整体低于 0.5，即与原文方向相反；只有 2 个队列的 signature 用 Hartung-Knapp 时置信区间极宽（中位宽 0.73），属预先规定方法的已知性质 | primary/signatures.csv | Claude 记录 |
| 2026-11-11 | 未做（SAP 列出、本轮未覆盖）：转移终点（仅 GSE116918）的探索性描述；汇总 C 与原文自报内部 C/AUC 的关系（原文内部 C/AUC 尚未提取） | — | 待办 |
| 2026-11-12 | 方向复核（Wei 要求，看到结果后进行）：35143414 与 42418409 由 2 个子代理只看原文独立复核（事先不告知外部结果），两人一致：方向明确，锁定的公式/系数/配对规则（IRG1>IRG2 记 1）与原文一致；35143414 为膀胱+前列腺+肾三癌联合模型。外部 C<0.5 判定为真实结果，非提取或计算错误。专家D 与 Wei 待填 | review/BRIEF_DIR2.md、out_DIR2a/b.jsonl、方向复核表_35143414_42418409_2026-11-12.xlsx | 子代理×2；待专家D、Wei |
| 2026-11-12 | 探索性分析预先规定（运行前写入，均不作结论）：E1 配对 bootstrap：同一次重抽样中算 signature C − 分级 C（ΔC），1,000 次，按 signature 对 ΔC 做随机效应汇总（REML+HK），报告 ΔC>0 / <0 / 含 0 的个数；E2 增殖基准 meta-PCNA（Venet 2011，PMID 22028643）：在 TCGA-PRAD（参照队列，不用结局）中取与 PCNA 正相关最强的前 1% 基因（Pearson，不含 PCNA 本身），各队列内 z 分数取均值为增殖评分；报告其 C、signature C − meta-PCNA C 的配对 bootstrap ΔC 汇总、signature 评分与 meta-PCNA 的 Spearman 相关；E3 转移终点：GSE116918（放疗队列）用 time_met/met，按与 BCR 相同的可评估规则与方法计算，只做描述；E4 原文自报区分度：子代理提取原文对同一模型报告的 C 指数与时间依赖 AUC；“原文训练集区分度”取训练集 C 指数，无则训练集 3 年 AUC，再无则 5 年、1 年；外部验证取原文外部队列同类指标的中位数；与本研究汇总 C 比较（差值中位数、Spearman 相关），C 与 AUC 分开报告 | 05_code/scripts/exploratory.py（待写） | Claude 预先规定 |
| 2026-11-12 | 新颖性检索（PubMed，9 条检索式，545 条记录）：未发现在前列腺癌中用原文系数、控制数据泄漏、以临床分级/增量价值/随机基因集为基准、预先锁定分析计划的大规模重评；最接近者为把已发表 signature 当作新模型对照的研究（如 PMID 38395379、38623253）与只做叙述汇总的综述（PMID 38112859）；他癌种有先例（肺癌 PMID 28200038、卵巢癌 PMID 24700801、乳腺癌 PMID 22028643）| 01_search/novelty/ | Claude 执行 |
| 2026-11-12 | 探索性分析完成（按预先规定执行，均不作结论）：E1 配对 bootstrap（signature C − 分级 C）：142 个中 73 个 CI 整体 <0、69 个含 0、0 个 >0，汇总 ΔC 中位数 −0.143；E2 meta-PCNA（TCGA 中与 PCNA 相关最强的 168 个基因）各队列 C 0.47–0.66，signature − meta-PCNA ΔC 中位数 +0.022，16 个 CI>0、1 个 <0；E3 转移（GSE116918，248 例 22 事件）：77 个 signature C 中位数 0.571，Gleason 总分组 C 0.599，11 个 CI>0.5；E4 原文自报区分度（5 个子代理盲法提取，142 个中 102 个有训练集指标，87 个为 3 年 AUC）：原文训练集中位 0.767 vs 本研究外部汇总 C 0.583，差值中位 0.20（99% 为正），两者无相关（Spearman −0.16）；原文外部验证（62 个）中位 0.709，差 0.12 | 06_results/evaluation/exploratory/、E3_met_GSE116918/；结果表 UroSig_结果表_主分析与敏感性分析_2026-11-12.xlsx | Claude 执行 |
| 2026-11-12 | 实现说明：meta-PCNA 候选基因限定为 4 个主分析队列都测到的基因（16,775 个，取前 1% = 168 个）；E4 另加描述性比较“原文 3 年 AUC vs 本研究 36 个月 td-AUC 中位数”（97 个，0.769 vs 0.601） | exploratory.py、compare_internal.py | Claude 执行 |
| 2026-11-12 | 方向复核定稿：35143414、42418409 由 2 个子代理、Wei 与另一位泌尿外科专家（专家D）独立复核，四方一致：提取与计算正确，外部 C<0.5 为真实结果，锁定数据不改 | 方向复核表_35143414_42418409_2026-11-12.xlsx | Wei |
| 2026-11-12 | Yin 2024（PMID 38395379）全文已读：102 个已发表 signature = Li 2022（PMID 35358302）整理的 30 个 + PubMed 近 4 年 72 个（含 Decipher、Prolaris、GPS）；在 7 个长随访队列（含 TCGA-PRAD，即多数 signature 的训练集）中“计算评分、算 C 指数”，未说明是否用原文系数，未处理数据泄漏，未与临床分级或随机基因集比较，只用于给自己的 9 基因模型排名；新模型的候选基因又是在这 7 个队列中筛出来的 | 01_search/novelty/fulltext/38395379_Yin2024_CancerLett.pdf | Claude 阅读 |
| 2026-11-12 | 新颖性检索补漏：经 Yin 2024 参考文献追溯发现 Li R 等 2022 Cancer Res（PMID 35358302）——30 个已发表前列腺癌 signature、15 种机器学习算法、10 个数据集 1,558 例的“系统评估”，此前 PubMed 检索式未命中（题目用 “gene expression signatures… prognosis”，无 published/re-evaluat 等词）。摘要显示其做法是用各 signature 的基因列表重新训练模型，评估的是“基因集+算法”，而非原文模型本身；全文未开放，待 Wei 下载核实 | 01_search/novelty/新颖性检索报告 | Claude 执行 |
| 2026-11-12 | Li 2022 Cancer Res（PMID 35358302）全文已读（Wei 下载，存 01_search/novelty/fulltext/35358302_Li2022_CancerRes.pdf）：30 个 signature 只用基因列表、以 15 种算法重新训练（数据集内 10 折 CV / 数据集间训练-测试），原文明言“不复现原文预测模型”；10 个队列；无临床变量比较；未处理训练队列重叠；10 个随机基因集。与本研究为互补关系（其为 2020 年前经典 signature 的基因集潜力，本研究为 2020 年后 signature 的原文模型可迁移性；本研究 98% 发表于 2020 年后） | 新颖性检索报告（更新） | Claude 阅读 |
| 2026-11-12 | 补检（替代 WoS/Scopus，Wei 无账号）：Europe PMC 7 条检索式 + 5 篇种子文献双向引用追踪（332 条）+ Semantic Scholar 施引文献（57 条）；未发现新的同类前列腺癌研究；补充他癌种先例 2 篇预印本（三阴乳腺癌 62 个 signature、头颈癌 5 个）| 01_search/novelty/openalex/ | Claude 执行 |
| 2026-11-12 | 目标期刊定为 BMC Medicine；增量价值检验正文主报“合并 LR χ²”版本，“各队列最小 p”版本放附表（按 Claude 建议，Wei 未提异议） | — | Wei |
| 2026-11-12 | 事后分析预先写明（看到主结果后新增，运行前记录；V1 为验证性，V2–V4 为探索性，均不作主要结论）：V1 实现正确性：在原文用过的队列（TCGA-PRAD 及被泄漏规则排除的 GEO/MSKCC 队列）中按同一实现计算 Harrell C 与 36 个月 td-AUC（基因覆盖 ≥80%），与原文在同一队列报告的同类指标逐一比对（原文 C-index 对 Harrell C，原文 3 年 AUC 对 36 个月 td-AUC），报告差值中位数与一致性；并给出“原文自报 → 本研究在原文队列重算 → 独立外部队列”三段对比；V2 精度：用 E1 配对 bootstrap 的标准误，按 80% 把握度、双侧 α=0.05 计算各队列及汇总层面可检出的 signature−分级 C 差（MDD = 2.8×SE）；V3 外部表现的相关特征：以 logit 汇总 C 为因变量、1/方差加权，自变量为训练数据是否仅 TCGA（复核员记录的训练队列文字判定）、原文是否报告外部验证指标（E4）、基因数（log）、建模方法（复核员记录 coef_type 归为 LASSO 类/多因素 Cox/单因素或其他）、与 meta-PCNA 的相关（E2 中位 ρ），先单变量后多变量；V4 冗余：在 TCGA-PRAD（参照队列，不涉结局）及 4 个主分析队列中计算全部可评分主分析 signature 评分的两两 Spearman 相关，报告 |ρ| 分布、层次聚类（平均连接，距离 1−ρ）及各簇与 meta-PCNA 的相关 | 05_code/scripts/posthoc_v1_v4.py（待写） | Claude 预先规定 |
| 2026-11-12 | 事后分析结果：V1 实现正确性——在原文用过的队列中重算 338 个组合，Harrell C 中位数 0.642；与原文同队列同类指标逐一比对 165 处，原文减本研究差值中位 0.052（IQR 0.017–0.122），48% 相差 ≤0.05，BCR 类终点差 0.048、其他终点 0.080；三段对比（95 个）：原文训练集自报 0.761 → 本研究在原文训练队列重算 0.694 → 独立外部汇总 0.584；同一 signature 在原文用过与未用过的 GEO/MSKCC 队列中 C 中位 0.626 vs 0.597（配对差 +0.024，Wilcoxon p<0.001）。V2 精度——单队列可检出的 signature−分级 C 差（80% 把握度）0.14–0.23；汇总 CI 半宽 k=3/4 约 0.13。V3（探索性）——仅用 TCGA 训练（121/142）汇总 C 偏低（logit 系数 −0.115，95%CI −0.235 至 0.005），原文报告过外部验证者偏高（+0.093，0.009–0.178），基因数、建模方法、与增殖相关性无明显关联。V4（探索性）——TCGA 中 154 个评分两两 |ρ| 中位 0.30，18% 的配对 |ρ|≥0.5；按 ρ≥0.5 聚为 51 簇，最大簇 54 个 signature，其与 meta-PCNA 的 ρ 中位 0.49 | 06_results/evaluation/posthoc/ | Claude 执行 |
| 2026-11-12 | 流程图核对中发现重大缺口（需 Wei 决定）：题录阶段的双人人工筛选从未执行（筛选表“筛选者1/2 决定、终裁”全空），进入全文的 292 篇实为关键词规则“AI 预筛=可能纳入”的结果；其余 719 条（不确定 632、可能排除 87）未经任何人工筛选。按标题/摘要关键词粗查，其中约 79 条可能符合纳入标准（如 PMID 42729135、41246829、40860343）。当前样本不能称为检索结果的完整纳入集 | 02_screening/screening_title_abstract_PRAD.xlsx | Claude 发现，待 Wei 决定 |
| 2026-11-12 | 补完整题录筛选（Wei 决定）：未经人工筛选的 719 条（规则“不确定”632、“可能排除”87）由 2 个子代理独立盲筛（不看 AI 预筛标签、互不可见），按筛选表原纳入/排除标准；两人一致排除者排除，一致纳入或不确定者进入全文，分歧交 Wei/专家裁定；新进入全文者按原流程提取、盲法复核、专家核对后锁定 v2，再按 SAP v2 全部重跑（含敏感性与探索性分析）。原 292 篇已全部经全文评估，流程图中如实分列 | 02_screening/ta_dual/ | Wei 决定，Claude 执行 |
| 2026-11-12 | 补筛完成：719 条由 2 个子代理独立盲筛（3 批×2）。一致性：纳入/不确定 vs 排除 98.6%，Cohen κ = 0.972。结果：两人均排除 386、两人均纳入或不确定 323（进入全文）、分歧 10（待 Wei/专家裁定）。核查：6 个输出文件 PMID 与批次逐行一致、无重复；筛选过程中子代理共用临时目录导致脚本互相覆盖，各自发现后作废受影响判定并重做，最终文件核对无串批；A/B 理由文字完全相同 64 条，分散分布、多为短标签，判断为同模型措辞趋同而非互相抄录 | 02_screening/ta_dual/（out_*.jsonl、merged_TA.csv、待裁定_TA.xlsx） | 子代理×6，Claude 核查 |
| 2026-11-12 | 10 条分歧由 Wei 与另一位专家裁定：全部排除（理由与其中一位子代理的排除理由基本一致）。题录阶段最终：719 条中进入全文 323、排除 396 | merged_TA.csv | Wei + 专家 |
| 2026-11-12 | v2 全文：2 篇中文文献排除（37846626、37847083，语言，与第一版 3 篇中文文献处理一致）；Wei 手动下载 41 个 PDF，经 DOI/标题逐个核对对应 40 篇（1 个重复：42092422 两份完全相同，副本移至 03_extraction/v2/_to_delete/），已按 PMID 命名移入 03_extraction/v2/fulltext_manual/；Wei 报告 15 篇找不到全文；另有 5 篇未下载也未报告（38581330、32253914、38179769、32687727、41975070） | move_manual_v2.ps1 | Wei 下载，Claude 核对 |
| 2026-11-12 | v2 全文获取：PMC 开放获取（S3）240 篇 + Europe PMC 补到 1 篇 + Wei 手动下载 48 篇（DOI 逐个核对），共 289/321 篇有全文进入提取；32 篇无法获取全文（PMC 非开放 20 篇、Wei 找不到 12 篇；含 Wei 认为不必要的 38581330 预警期刊），流程图列为“全文无法获取”。按摘要判断其中可能合格的有 37559989、29511883、35968320、33934205、35731271、40654646、15931389、15067324、26195723、41246829 等，已在局限中说明 | 03_extraction/v2/v2_availability.csv；v2_需手动下载全文_第2批.xlsx | Wei 下载，Claude 核对 |
| 2026-11-12 | Wei 补下第 2 批 10 篇（DOI/标题核对 9 篇一致；40654646 下到的是其期刊正式版 Transl Oncol 2026，PubMed 中正式版为 42143470，不在 2026-10-26 检索结果内）。全文阶段排除非同行评议预印本（沿用第一版 39252915 的处理）：40654646、40458374、42922617、38559135、38076845、37205576、37034687、36798177（bioRxiv/medRxiv/Research Square），共 8 篇 | 03_extraction/v2/preprints_excluded.txt | Claude 按既有规则执行 |
| 2026-11-12 | 预印本规则修正（Wei）：检索命中的预印本若已有同行评议正式版，用正式版全文评估，以检索命中的记录计数（不另算一篇）。40654646（bioRxiv）用正式版 Transl Oncol 2026（PMID 42143470，DOI 10.1016/j.tranon.2026.102819）评估，文件 fulltext_manual/40654646_published_42143470.pdf（重复副本移入 _to_delete）。其余 7 篇预印本待逐篇查是否有正式版。38581330（Arginase 2，Alternative Therapies in Health and Medicine）Wei 已下载，纳入评估 | 03_extraction/v2/ | Wei |
| 2026-11-12 | 其余 7 篇预印本查正式版（PubMed UpdateIn 链接 + 作者/标题检索）：37034687 → 38191471（Nat Commun 2024，不在检索结果内，按 Wei 规则用正式版评估，PMC10774315）；38559135 → 40553565、38076845 → 39102549、37205576 → 41037000（正式版本身已在检索结果中并单独评估，预印本作重复排除）；40458374、42922617、36798177 未找到正式版，作“未经同行评议的预印本”排除 | 03_extraction/v2/preprints_excluded.txt | Claude 执行 |
| 2026-11-12 | 40458374 找到正式版（Mol Cancer 2026，Article in Press，DOI 10.1186/s12943-026-02795-3，Wei 下载，存 fulltext_manual/40458374_published_MolCancer2026.pdf），改为用正式版评估；42922617、36798177 未找到正式版，排除（Wei 确认） | preprints_excluded.txt | Wei |
## 2026-11-14 v2 双人提取（待同步到 D:\article\UroSig-Reeval\logs\decision_log.md）
- 15 批 302 行（300 个 PMID）X/Y 两份提取全部完成。
- 37034687、40458374 的预印本副本被同时放进 V05/V14，正式版在 V15。按 Wei 规则（有正式版用正式版）以 V15 为准，V05/V14 中这两条作废；compare_v2.py 加 SUPERSEDED 排除。
- compare_v2.py 修正：基因对 signature 双方均无系数时不比系数（31961962 误报）；负号统一（−、–、‐、- 归一）后比较。
- 逐条裁定见 03_extraction/v2/adjudication_v2.csv。
- 裁定表 03_extraction/v2/裁定表_v2.xlsx 已生成：R1/R2/R3/R7/R8/R15 沿用 v1 锁定做法；R4/R5/R6/R9/R10/R13/R14/R17 为 v1 未遇到的情形，待 Wei 与专家决定。不一致 24 篇中 18 条已由 Claude 回原文或按 v1 规则裁定（见“已裁定”表）。
- 2026-11-14 Wei：待裁定 8 篇（28724614、30127947、31936761、35405009、36607288、37108678、39418984、42415249）全部排除。
- 2026-11-14 Wei 授权 Claude 按最优做法确定 v1 未遇到的规则，Claude 定为：
  - R4 机制研究附带/非生存目标的多基因评分：组合评分做过生存检验即纳入（无系数进 S8）；未做生存检验的排除（37729607）。
  - R5 PCA 评分：报告载荷的按载荷作线性权重（32951005，原文高分=预后好，取负号）；未报告载荷的排除为“无法复现”（29556354、31107158、35495154、36203594、41828396）。
  - R6 预定义基因集直接打分：排除为非新 signature。
  - R9 系数作用于二分/分箱变量、阈值规则分组、或非加权和评分：排除为非线性模型（30637711、35204174、38576998、30024105、37349324、38302875、35531101、37176106、30108134、30540935），不另做事后敏感性，避免分析扩张。
  - R10 含临床变量：基因部分有单独系数的只评基因部分（24713434），拆不开的排除。
  - R13 癌旁组织建模的 34631510 排除；从血液选基因、在肿瘤组织打分的照常纳入。
  - R14 34804963 用 Table 2 的 B（公式印的是 exp(B)）；35592542 只有 HR，ln(HR) 进 S7。
  - R17 CCP 首篇（21310658）算新；29511883（既有基因集重训）不算新；38302875 按 R9 排除。
- 2026-11-14 Wei 提供 6 份补充材料（40216135 Table S2、28408174 Supp Table 13、34876840 Figure S1、39413958 Fig S1E、34131070 Table S5、21289624 Supp Tables 4–5），作为 V16 批由 X/Y 双人盲法重新提取。
- 2026-11-14 V16 补充材料双人提取：28408174、34131070、21289624 两人系数逐一相同 → main；34876840 只有 HR → sens_hr；39413958、40216135 系数只在柱/棒棒糖图上，两人独立估读一致（差 ≤0.01 / ≤0.2）→ sens_figure（按 R2）。
- 2026-11-14 37205576、42922617、36798177 的排除理由按预印本规则统一记为“预印本”。
- 2026-11-14 提取锁定 v2（03_extraction/locked/*_v2_2026-11-14.csv，MANIFEST_v2.txt）：v1 223 + v2 152 = 375 个可分析 signature（main 247）。sens_hr 系数 = 原文多因素 HR 取 ln；作者定义等权/等票 → ±1；equal_weight 方向取两人记录，相反记 0、一方未定取另一方；明显笔误按两人一致意见映射（CRND→CRNDE、CACNAID→CACNA1D、OR,2T2→OR2T2、OR,2T11→OR2T11、DNRB→EDNRB、HOXB8+→HOXB8、NF-κB(p65)→RELA）。基因对规则新增 31961962 gt、38191471 gt、34876840 gt。
- 2026-11-14 泄漏矩阵 v2：正文扫描与两位提取员记录一致 125 个；不一致 27 个由 Claude 回原文核定（out_LEAK_v2.jsonl，多数为 Abcam/Broad 等机构所在地 “Cambridge” 误命中；“其他”用途按 v1 规则计为用过）。
- 2026-11-14 代码按锁定版本切换路径（05_code/urosig/paths.py，UROSIG_LOCK=v1 可复现 11-01 结果；v2 结果写入 06_results/evaluation_v2）。覆盖率与可评估性 v2 已算出（0.8 门槛：main 213 个 signature 有 ≥1 个可评估外部队列）。
- 2026-11-14 14:08 按 SAP v2 在 Wei 电脑上启动 v2 全部重跑（run_all：主分析 + S1–S10，14 并行）。
- 2026-11-14 v2 新增 80 个 main signature 的原文自报区分度提取（IM2_1–4，供 E4/V1–V3）已派出。
- 2026-11-17 v2 重跑：11-14 与 11-17 16:28 两次由 Claude 远程启动的运行均在桌面端组件重启时被终止（无 Python 报错）；18:08 的运行完成主分析评估（569 对，全部成功），但 run_all 仍按 v1 写死的 384 对核对而停在合并步骤。修正：run_all 从 v2 可评估性表计算应有对数，并可从已完成分片续跑（--force 才重算）；19:26 续跑 S2 起。v2 主分析：211 个 signature / 569 对，合并 C 中位数 0.573（IQR 0.535–0.618；v1 0.575），优于分级的论文 1.4%（v1 1.4%），LR 合并检验 FDR<0.05 共 61 个（v1 39）。
- 2026-11-17 数值修正：S7 中一个从图读系数的 signature（系数绝对值约 10）在增量 Cox 拟合时 exp 溢出导致 shard 失败。修正为在 incremental() 内先把评分标准化（C 指数与 LR 检验对仿射变换不变，不改变定义），并跳过不收敛的 bootstrap 重抽样。已完成的 run 理论上不受影响。
- 2026-11-18 v2 全部重跑完成（主分析、S1–S10、E1–E4、V1–V4）。两处汇总修正：S6 首次汇总误读 v1 终点分类表（pool_results 路径未同步），已用 v2 表重算（149 个 signature）；V3 的训练队列/建模方法协变量原只读 v1 提取记录，v2 的 69 个 signature 被归入缺省类，已改为同时读 reconciled_v2.csv 后重算。结果表：06_results/UroSig_结果表_v2_2026-11-17.xlsx。

- 2026-11-18 Wei 决定（C 阶段）：补校准斜率；SAP 上传 OSF、代码与派生数据上 GitHub+Zenodo；医院不要求出具伦理豁免证明；资助/利益冲突/作者贡献最后再写；图用 Python。
- 2026-11-18 事后分析 V5（校准斜率，calibration_slope.py）：对主分析中有系数的线性 signature，在未标准化 log2 表达上计算原文线性预测值，作为唯一协变量拟合 Cox，系数即校准斜率。564 对 / 208 个 signature：逐 signature 合并斜率中位数 0.28（IQR 0.08–0.78），65% 的 95% CI 完全低于 1，4% 高于 1，17% 斜率 ≤0。局限：原文无基线风险，无法评估绝对风险校准；斜率混合了过拟合与平台尺度差异。
- 2026-11-18 主图 1–5 初稿（make_figures.py）。Wei 反馈"简陋、AI 味重、图太少"，该版作废，未同步到本机。
- 2026-11-19 主图重做为 6 张（make_figures_v2.py，输出 07_submission\figures，含 PDF/SVG/TIFF 600 dpi/PNG 与 source_data 27 个 CSV）。按 BMC Medicine 图片要求：宽 170 mm、高 ≤225 mm、图例键在图内、标题与图注不进图、字体嵌入（Arial）、线宽 ≥0.5 pt、单文件 ≤10 MB。全套固定配色语义：蓝=已发表 signature，朱红=Gleason 分级分组，灰=随机基因集，赭=增殖指数（meta-PCNA），炭黑=临床模型，紫=原文报告值。
  图1 研究设计与证据库；图2 外部区分度全景（排序图+逐队列热图+注释轨+KM 曲线）；图3 基准比较（随机基因集、临床模型增量、分级、增殖指数）；图4 报告值与外部表现、校准斜率；图5 冗余性、基因复现、Meta 回归、发表年份；图6 稳健性（敏感性分析、指标一致性、转移终点、检验效能）。
- 2026-11-19 新增 figure_inputs.py（在本机运行，输出 06_results\evaluation_v2\figure_inputs）：随机基因集 C 的逐队列直方图；4 队列合并 KM 曲线及在险人数（ISUP 分级分组；排名第一 39524226 与排名居中 35891787 两个 signature 的队列内三分位，仅取 4 个队列均可评估者；分层 log-rank P：分级 6.4e-19，第一名 6.8e-06，居中 0.50）；signature 评分与 ISUP 的 Spearman 相关（中位数 0.15）。属描述性作图输入，不改变任何已报告的统计结果。
- 2026-11-19 图5b 取"被 ≥5 个 signature 使用的基因"（21 个），并列按基因名字母序，保证不同环境下结果一致。
- 2026-11-19 Wei 认可主图新样式。附图 S1–S10 完成（make_supp_figures.py，输出 07_submission\supplementary_figures）。make_figures_v2.py 末尾加 `if __name__ == "__main__"` 以便附图脚本复用其样式与数据表。
- 2026-11-19 校准斜率表述修正（未改 V5 计算）：逐队列斜率中位数 GSE54460 0.13、GSE70769 0.69、GSE70768 0.93、MSKCC 0.81；合并值 0.28 主要由 GSE54460 决定（逆方差权重中位数 82%），去掉该队列后合并中位数 0.81、CI 低于 1 占 31%。GSE54460 的线性预测值 SD 约为其余队列的 2.4 倍，提示平台尺度差异。主图 4e 改为逐队列展示；正文应报告逐队列数值并在局限性中说明。此为 Claude 按"最优策略自主执行"所做的呈现方式调整，待 Wei 过目。
- 2026-11-19 作者与单位（Wei 提供）：赵瑾（共同第一）、薛松（共同第一）、苏瑞平、雷炜、袁图斌、万象（通讯）。赵、薛、苏、万：西安市人民医院（西安市第四医院）；袁：贵州茅台医院；雷：深圳市坪山区中心医院。赵、万为肾脏内科，其余为泌尿外科。英文名：深圳市政府英文网站写作 Pingshan District Central Hospital of Shenzhen（已核验）；西安与茅台两家官网未见英文名，暂用 Xi'an People's Hospital (Xi'an Fourth Hospital)、Kweichow Moutai Hospital，未核验，待作者本人确认。
- 2026-11-19 Wei 说明：筛选与提取结果已由赵瑾、薛松、万象对照原文逐条核对。公开材料按此表述；核对记录请作者留存。
- 2026-11-19 公开方案（Wei 决定）：SAP 现在公开，代码与派生数据投稿时公开。OSF 自 2026-11-16 起停止新建项目（帮助文档 help.osf.io/article/727-osf-projects-transition），改为在 Zenodo 单独存一份 SAP（原件 PDF、英文译本、偏离与事后分析说明，位于 08_public_release\osf）。
- 2026-11-19 GitHub 私有仓库已建：https://github.com/wlei9534/UroSig-Reeval（空仓库，待整理公开包后推送）。
- 2026-11-19 主图 1 流程图文字由 "by two independent reviewers" 改为 "independently in duplicate"（补筛由 AI 双份独立完成，避免被理解为两位人工筛选者）。
- 2026-11-19 单位英文名按 PubMed 近期论文的署名核定（Wei 要求；E-utilities 实时检索，2026-11-20 UTC）：① 西安：2025–2026 年署该院的 396 篇中取最新 200 篇，643 条署名里 "Xi'an People's Hospital (Xi'an Fourth Hospital)" 449 条，为主流写法；薛松 2026 年论文（PMID 42253846）、苏瑞平 2026 年论文（PMID 42558341）均用此写法并加 Department of Urology；肾脏内科见 "Department of Nephrology, Xi'an People's Hospital (Xi'an Fourth Hospital), Affiliated People's Hospital of Northwest University"。是否在医院名后加挂西北大学（各篇写法不一）留待西安作者自定。② 茅台：署 "Moutai Hospital" 的 120 篇，主流写法 "Kweichow Moutai Hospital, Renhuai, Guizhou, China"；未检到袁图斌本人的 PubMed 记录。③ 深圳：61 篇，"Pingshan District Central Hospital of Shenzhen" 60 条，与深圳市政府英文网站一致。此前记为"未核验"的两家改为已按 PubMed 核定。
- 2026-11-19 Zenodo 草稿已建（未发布）：https://zenodo.org/uploads/24385588，类型 Publication / Other，CC BY 4.0，六位作者及单位已填。文件上传工具无法读取本机文件，三份文件待 Wei 手动拖入；发布须 Wei 确认。
- 2026-11-20 Zenodo 记录已发布（Wei 上传三份文件并确认无误后发布）：https://zenodo.org/records/24385588 ，DOI 10.5281/zenodo.24385588（所有版本 10.5281/zenodo.24385587），CC BY 4.0。内容：SAP v2 原件 PDF、英文译本、偏离与事后分析说明。非前瞻性注册。
- 2026-11-20 公开包已整理并推送到 GitHub 私有仓库 wlei9534/UroSig-Reeval（首个提交 874 个文件、约 23 MB；08_public_release\build_release.py 生成，本机副本在 08_public_release\UroSig-Reeval）。不含文献全文、队列数据、检索导出与题录筛选表。公开前待办：GEO/TCGA 再利用条款核对；题录筛选表是否公开。
- 2026-11-20 Wei 决定：西安作者的单位不加挂西北大学；GEO 与 TCGA 的再利用条款无问题（仓库本就不分发队列数据）；README 中队列构建脚本的说明以后再核。
- 2026-11-20 题录筛选表是否公开（Wei 交 Claude 决定）：公开筛选决定，不公开摘要。02_screening 的原表含摘要原文，仓库只放导出的 screening_decisions.csv（1,011 条：PMID、DOI、年份、期刊、标题、关键词预筛标签与理由、两次 AI 筛选的决定与理由、结果、终裁）及给筛选者的说明 BRIEF_TA.md。
- 2026-11-20 更正：此前说"结果摘要和 Excel 结果表里仍写着校准斜率 0.28"有误，两者原本都没有 V5 内容。现已补入：结果表新增 V5_校准斜率 页（逐队列与合并），《结果摘要_v2锁定》补 V5 一节。
- 2026-11-20 资助与声明（Wei 提供）：深圳市坪山区卫生健康系统科研项目，编号 2025201；无利益冲突；作者贡献：赵瑾、薛松、万象负责筛选与提取结果的核对，苏瑞平、雷炜、袁图斌参与数据清洗与撰写。Declarations 草稿见 07_submission\Declarations_draft_2026-11-20.md。资助项目英文名未核验（Europe PMC 接口超时、网页搜索无结果），待 Wei 按立项文件或本院已发表论文核对。
- 2026-11-20 Wei 确认：资助项目英文写法 Pingshan District Health System Research Project, Shenzhen (grant number 2025201)；资助方未参与研究设计、分析与撰稿。作者贡献（Wei 说明）：万象提出并设计研究、指导监督；赵瑾、薛松、万象核对筛选与提取结果；雷炜、袁图斌、苏瑞平完成分析并参与数据清洗与撰稿。Declarations 草稿已按此更新；仓库 README 补资助信息。
- 2026-11-20 备注：GitHub 私有仓库的提交历史中留有一版早期 Declarations 草稿（提交 48afb7e，已在 b8e96db 移出工作区）。公开仓库前需由 Wei 决定：重写历史，或另建干净的公开仓库。
- 2026-11-20 作者贡献补一句 "WL contributed to the study design"（Wei 同意），与 Zenodo 已公开材料中 W.L. 决定纳排规则、裁定分歧的记录一致。
- 2026-11-20 开始写作（Wei 指示）。正文初稿 v1 写成 Claude 文档：https://claude.ai/artifact/K6h3mhBr2o5pR8QAHutd5w （题名页、摘要、Background、Methods、Results、Discussion、Conclusions、缩写、Declarations、参考文献 34 条、图题图注、表 1、附加文件清单）。参考文献逐条经 PubMed E-utilities 实时核对（第 32、33 条不在 PubMed，经 Crossref 核对）；第 25–29 条只核对了第一作者，作者列表待排版时由 PMID 生成。待办：定标题；表 1 补平台与分级分期、表 2 生成；附加文件 1、3–6 汇编；语言润色；TRIPOD+AI 清单；投稿前核对期刊对 AI 披露的要求。

## 2026-11-20（续）表格与附加文件；Fig. 6d 更正

- **表与附加文件**：`05_code/scripts/build_tables_additional_files.py` 生成 `07_submission/tables/Table1_cohorts.csv`、`Table2_signature_characteristics.csv`，以及 `07_submission/additional_files/` 下的 Additional file 3（1,011 条记录的筛选与全文结局、375 个可分析 signature、4,218 个基因项、218 篇全文排除及理由）、Additional file 4（主分析与全部敏感性/探索性/事后分析结果）、Additional file 6（各文章对 6 个队列的使用情况、signature×队列的基因覆盖与可评估性）。三个工作簿全部为英文；提取员用中文写的自由文本（原文终点、训练/验证队列、排除备注、泄漏判定依据）不进入附加文件，改用编码类别，原始记录仍在仓库的 `03_extraction/`、`06_results/eligibility/`。脚本末尾有中文字符检查，残留即报错。
- **全文排除理由的英文归类**：按 `signatures_locked_v2_2026-11-14.csv` 的“最终状态”前缀归入 14 类；第二阶段各类计数与 `06_results/flow_v2_2026-11-14.md` 完全一致（非新 58、非线性 29、非预后 21、无法复现 14、非转录组 7、非 mRNA 级 7、非基因级线性 3、非前列腺癌 3、预印本 3、原文矛盾 1、非多基因 1、非肿瘤组织 1，合计 148）；第一阶段 70 篇：非线性 20、非多基因 13、非新 8、非前列腺癌 5、非基因级线性 5、非预后 4、非 mRNA 级 4、非转录组 4、中文文献 3、预印本 2、全文未获取 2。
- **Additional file 1、2**：由 Claude 在云端用 pandoc + XeLaTeX 生成 PDF（源文件为同名 .md）。AF1 = 检索式、纳排标准、筛选过程、排除理由表、队列处理、打分规则、敏感性分析结果表、探索性与事后分析、SAP 修订；AF2 = 附图 S1–S10 及图注（矢量图取自 `07_submission/supplementary_figures/`）。
- **更正（结果已知之后，事后描述性指标）**：Fig. 6d 原来把“0.5 + 2.8 × signature 与分级之差的自助法标准误”标成“80% 把握度下可检出的最小 C-index”（0.64–0.73，达到者 2%–15%），两者不是同一个量，阈值被高估。现改为按 C-index 自身的标准误计算：expit(2.8 × 各队列 logit C 自助法标准误的中位数)，GSE54460 0.62、GSE70769 0.63、GSE70768 0.71、MSKCC 0.67，达到者分别为 18%、23%、7%、16%。`V2_mdd_by_cohort.csv`（与分级之差的最小可检出值 0.14–0.23）不变，正文按其本义另行报告。`make_figures_v2.py` 的 Fig. 6d 段、图注、正文 Robustness 段、AF1 第 8、9 节同步修改；其余结果不受影响。发现人：Claude（准备附加文件时核对代码定义）。
- **正文同步的事实性修改**：Methods 写明 MSKCC 结局为 cBioPortal 的 DFS 字段、与 BCR 不完全等同；Results 的 meta 回归补报另外三个协变量（基因数每翻倍 +0.044，95% CI 0.008–0.079；与增殖指数相关性 −0.129，95% CI −0.256 至 −0.002；建模方法无关联），避免只报显著项中的两项；正文 Table 1、Table 2 填入完整内容。
- **待 Wei 决定**：Zenodo 上已发布的两份英文文件（SAP 英文译本、Amendments）正文里写的是“this OSF project”，实际存放在 Zenodo；文件已发布不可改，如需更正须在同一 concept DOI 下发新版本。

## 2026-11-20（续）Zenodo 第二版

- **Wei 决定**：在同一 concept DOI 下发布第二版，更正两份英文文件中“this OSF project”的说法（文件实际存放在 Zenodo）。
- **第二版**：https://zenodo.org/records/24390581 ，版本 DOI 10.5281/zenodo.24390581；concept DOI 10.5281/zenodo.24385587 不变（正文引用的即此 DOI，始终指向最新版）。三个文件的 MD5：`Amendments_and_post_hoc_analyses.md` 635f47ed58acb632bbadbafa4546c10f；`SAP_v2_English_translation.md` e81a532a3c5eb5de50fa02a7d72a48a8；`SAP_v2_original_Chinese_2026-11-11.pdf` e475832f1907454c96a520c25abe3755（与第一版相同，未改动）。
- **与第一版的差异**：两份英文文件的“OSF project”改为“Zenodo record”，并各加一段版本说明；Amendments 第 3、4 节加入 Fig. 6d（V2）更正。元数据（标题、作者、单位、许可）沿用第一版。
- **发布前核对**：首次上传的 Amendments 与准备稿不一致（电脑上的文件在 00:02 被改动，缺第 6 节“Use of artificial intelligence”），未发布；经 Wei 确认恢复完整版，重新上传并核对 MD5 后发布。
- 本地文件：第二版在 `08_public_release/zenodo_v2/`，第一版保留在 `08_public_release/osf/`；`build_release.py` 改为从 `zenodo_v2/` 复制到仓库 `00_protocol/`。

## 2026-11-20（续）投稿格式核对；附加文件重新编号；TRIPOD+AI 清单

- **依据**：BMC Medicine 投稿须知（link.springer.com/journal/12916/submission-guidelines 及其 research-articles 页，2026-11-20 读取）：结构式摘要不超过 350 词；关键词 3–10 个；图、表标题不超过 15 词，图注不超过 300 词；表格不得用颜色或底纹；附加文件须按在正文中首次引用的顺序编号，并逐个列出文件名、格式、标题、说明；正文双倍行距并带行号、页码；LLM 的使用须写入 Methods。
- **附加文件按正文首次引用顺序重新编号**（原编号 → 新编号）：补充方法 1 → 1；TRIPOD+AI 清单 5 → 2；记录/signature/排除 3 → 3；队列使用与基因覆盖 6 → 4；附图 S1–S10 2 → 5；结果表 4 → 6。正文、AF1 内的交叉引用、文件名、生成脚本同步修改；旧文件名的文件移入 `_superseded/`。正文 Methods 中依次加入各附加文件的首次引用。
- **Additional file 2（TRIPOD+AI 清单）**：按 TRIPOD+AI 声明（PMID 38626948）的条目编号逐条注明报告位置（按章节与段落，不用页码）。条目主题为缩写，未照抄官方原文；仅适用于模型开发的条目标注“Development only”。未报告或未涉及的条目如实注明：3c、14（健康不平等、公平性，未涉及）、5b（入组与随访日期，见各队列原始文献）、20c（与开发数据的正式比较未做）。
- **因清单核对而补入正文的两句**：Methods 增加“患者与公众未参与本研究的设计、实施与报告”（条目 19）；Discussion 局限性增加“未按种族等社会人口学分组考察表现，这些变量不在提取的队列变量中”（条目 3c、14、23a）。
- **参考文献**：第 25–29 条作者经 PubMed E-utilities 核对后补全（Cerami E, Gao J, Dogrusoz U, et al.；Gao J, Aksoy BA, Dogrusoz U, et al.；Harrell FE Jr, Lee KL, Mark DB；Uno H, Cai T, Pencina MJ, et al.；Blanche P, Dartigues JF, Jacqmin-Gadda H）；第 32、33 条 DOI 改为 https://doi.org/ 形式；删去参考文献节开头的工作说明。
- **Table 1 平台名称**经 GEO（E-utilities，db=gds）核对：GPL10558 = Illumina HumanHT-12 V4.0 expression beadchip；GPL10264 = Affymetrix Human Exon 1.0 ST Array；GPL25318 = Almac Diagnostics Prostate Disease Specific Array (DSA)。
- **缩略语表**补入 FFPE、ISUP、LASSO、PSA、RFS。
- **Word 稿**：`05_code/scripts/format_for_submission.py` 对从在线文档导出的 Word 文件设置 A4、双倍行距、连续行号、页码，去除表格底纹，并报告摘要词数与残留的方括号工作文字。
- **仍待 Wei 决定或提供**：题目三选一；通讯作者邮箱用机构邮箱还是 Gmail；赵瑾、万象、袁图斌的科室英文写法；GitHub 仓库公开方式（历史中早期 Declarations 草稿的处理）及随后的代码归档 DOI；附加表格是否需要另存为 .xls（须知写的可接受表格格式为 .xls 或 .csv，现为 .xlsx）。

## 2026-11-23 投稿前总核对（Wei 要求：做不可撤销的步骤之前，全部内容先核对到无需再改）

**范围**：在线正文稿（rev 83）、表 1–2、主图 6 张及图注、Additional file 1–6、Zenodo v2 已发布的两份文件、代码与方法表述、参考文献 15–17 的描述。数字全部用 `06_results/evaluation_v2` 的结果文件重算；方法表述逐条对代码。

**核对无误（不改）**：211 个 signature / 210 篇 / 569 个评估；434 例 140 个事件；合并 C 中位数 0.573（0.535–0.618）；CI 高于 0.5 的 51 个；分级 0.72；优于分级 1.4%；E1 −0.144（110/101/0）；增量 0.0005、FDR 显著 61 个；E2 +0.025（22/5/184）；自报 vs 外部 0.755/0.581、99%、ρ −0.10；三步 0.75/0.69/0.58；+0.023（P 0.0002）；外部验证自报 0.71 vs 0.59（96%）；校准斜率 0.13/0.69/0.93/0.81；1360 个基因、78%、NOX4 9 个；元回归各系数；敏感性分析表；流程图各数；表 1、表 2 各格；图 2、4、5、6 上的全部数字。

**发现并已更正（均为表述或标注，分析结果未变）**

1. **随机基因集比较的分母**。3 个基因对 signature（PMID 32593277、35143414、42553741，共 5 个评估）没有做随机基因集比较（代码只对加权和类型的评分生成随机集），此前正文未说明。现写明：59.7% 的分母是 207 篇文章（124/208 个 signature；PMID 32264916 的两个模型各计半篇，即 123.5/207）；逐评估的比例应为 190/564 = 34%（原写 33%，图 3b 标注用的是 `rand_pct > 95` 的 187 个，与判定标准 `beats_random_95` 的 190 个差 3 个，现统一用后者）。另补报按 SAP 预先规定的 BH 校正结果：106/208 个（51%）q < 0.05。
2. **图 1d 与正文的“articles”**。图中百分比是按 211 个 signature 算的（95.7、36.5、33.2、30.3、24.6、10.0），坐标轴、图注和正文却写成“文章”。数字不动，坐标轴、图注、正文统一改为按 signature 表述。
3. **自报与外部表现的比较混用了指标**。158 个有训练集自报值的 signature 中，129 个报的是 3 年 AUC、5 个是其他时间点的 AUC、24 个是 C 指数，而外部值是合并 C 指数。E4 当时已算好同口径结果（`summary_E4.json`），现补入正文和 AF1：自报 C 指数（24 个）0.73 vs 外部 C 指数 0.59；自报 3 年 AUC（146 个）0.76 vs 本研究外部 36 个月 AUC 0.60。三步对比（144 个）的第二步本来就是按原文所报指标在训练队列重算（123 个 AUC、21 个 C 指数），现在正文、图注、AF1 中写明；并写明重算用的是本研究数据中的 BCR 结局（144 个中 108 个原文终点为 BCR/RFS）。Discussion 中“可能反映所报指标不同”一句相应改为“终点定义与数据处理不同，或选择性报告”。
4. **临床模型的 C 指数**。每个队列内是常数（GSE54460 0.767、GSE70769 0.757、GSE70768 0.803、MSKCC 0.831），原文“中位数 0.77”是对评估取的中位数，改为“0.76 至 0.83”；图 3c 标注同步修改。AF1 补充：分级按线性项进入；变量在该队列至少半数病人有记录才纳入；完整病例数（92/100、86/92、110/111、130/131）。
5. **软件表述有误**。原文写了 statsmodels，实际代码从未调用。实际情况：Harrell C 用 lifelines；Cox 模型（每 SD 的 HR、临床模型、校准斜率）和随机效应 meta 分析是自己用 NumPy 写的（`urosig/cox.py`、`urosig/meta.py`，文件日期 2026-10-26，早于 SAP 起草），测试中 Cox 与 lifelines 的系数、对数似然相差不超过 1e-4，REML 与直接最大化限制性似然一致。SAP 第 10 节软件一项写的是“pandas、lifelines”，未提自写的 Cox 实现；这不是看到结果后的改动，但属于与计划文字不一致，已写入正文 Software 和 AF1 第 6、9 节。运行环境版本取自 `run_primary/run_info.json`：Python 3.12.0、NumPy 2.3.5、pandas 3.0.1、SciPy 1.16.1、lifelines 0.30.0。
6. **探索性分析的时间表述**。正文 Methods 原写四项探索性分析都是“主结果之后规定”，而 AF1 和 Zenodo 文件写的是 E3（转移）、E4（自报 vs 外部）列在 SAP 中、E1 和 E2 为主结果之后规定。以 SAP 为准，改正文。
7. **“用过的队列表现更好”一句**。原句读起来像 151 个 signature 全都更高；实际是 151 个中 63%（95 个）更高，中位差 +0.023，Wilcoxon 符号秩检验 P = 0.0002；“用过的队列”不含 TCGA-PRAD。正文、图注、AF1 已写明。
8. **元回归**。模型里有建模方法一项（多因素 Cox 对 LASSO 类 −0.022，95% CI −0.091 至 0.047；单因素或其他 +0.142，−0.047 至 0.331），图 5c 未画。图注写明“模型还含建模方法，未显示”，系数写入 AF1。
9. **校准斜率**的评估数是 564 个（208 个 signature；不含基因对模型），图注和 AF1 写明。
10. **补一篇应引而未引的前期工作**：Manjang 等 2021（Front Genet 12:649429，PMID 34367234，PMC8340016）。该文用 32 个前列腺癌 signature 的基因集（基因集第一主成分，不用原文模型）在 TCGA-PRAD 单一队列中与 1,000 个随机基因集比较。它在 11 月 12 日的新颖性检索报告中已列出，但成文时没有引用；本研究的结论之一正是“与随机基因集难以区分”，不引不妥。现作为第 18 条参考文献加入 Background 和 Discussion，其后各条顺延（原 18–34 变为 19–35），正文、表 1、AF2 中的编号已同步，并逐条程序核对。同时把 Li 2022 的描述补全为“用基因列表重新训练，并与 10 个随机基因集比较”（依据：11 月 12 日读过的作者接受稿全文；今天又在其 bioRxiv 预印本上核对了一次）。Discussion 中“以随机基因集为基准”相应改为“与每个 signature 匹配的随机基因集”，不再暗示前人没做过随机基准。

**核对了但不需要改的**：参考文献 15–17 的描述（Li 2022 重新训练、不与分级比较；Yin 2024 把 102 个已发表 signature 当对照；An 2023 汇总 71 个 signature 及其自报结果）与 11 月 12 日全文阅读记录及今天的复查一致。Li 2022 的卷期页码（82(9):1832–1843）与 Europe PMC 记录一致。

**改动的文件**
- 在线正文稿：rev 83 → rev 125（23 处内容修改 + 参考文献插入与重新编号）。摘要未动（319 词）。
- `05_code/scripts/make_figures_v2.py`：图 1d 坐标轴文字、图 3b 百分比的算法、图 3c 临床模型 C 的标注。主图 1、3 已在电脑上重出；其余主图和附图未动。
- `07_submission/additional_files/`：AF1（6 页 → 7 页）与 AF2（参考文献编号）的源文件和 PDF。AF3–AF6 未动。
- `07_submission/Manuscript_BMC_Medicine_v3.docx`（28 页）：本次 Chrome 扩展未连接、内置浏览器未登录，无法走“在线文档导出”的老办法；改用 `05_code/scripts/update_export.py`，以 11 月 20 日的导出稿为底，把在线文档当前文本逐段写入（样式、编号、表格沿用导出稿），再跑 `format_for_submission.py`。已程序核对：Word 稿 152 段和两张表的文字与在线文档 rev 125 完全一致。v2 及更早的 Word 稿已过期。

**未做、留待最终一次性处理的（不可撤销）**：Zenodo 已发布的 Amendments 文件（v2）没有第 5 条的软件说明；内容与 AF1 不矛盾，但不如 AF1 完整。按 Wei 的要求，等全部内容定稿后再决定是否发 v3，现在不动。
