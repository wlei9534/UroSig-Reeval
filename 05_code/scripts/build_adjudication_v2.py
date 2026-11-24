"""Build 裁定表_v2.xlsx for Wei and the expert from reconciled_v2.csv and the rule list below.

Usage: python build_adjudication_v2.py <v2_dir>
"""
import re
import sys
from pathlib import Path

import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill

D = Path(sys.argv[1])
R = pd.read_csv(D / "reconciled_v2.csv", dtype=str).fillna("")
adj = pd.read_csv(D / "adjudication_v2.csv", dtype=str).fillna("")

# ---- rules: (code, rule, v1 precedent, Claude's proposal, affected PMIDs, needs Wei) ----
RULES = [
    ("R1", "作者明确定义的等权评分（全部 +1 或 ±1，含等票基因对投票）", "有：v1 记 equal_signed，进 S8",
     "按 v1 先例记 equal_signed（S8）；基因对等票投票按投票数计分", "21672909;39402324;31961962;37034687;42351650;21310658", "否（沿用）"),
    ("R2", "系数只印在图上、从图读出", "有：v1 记 sens_figure，进 S7",
     "按 v1 先例；两人读数差 ≤0.01 时取 X", "42729135;41988531;41584032;32687727", "否（沿用）"),
    ("R3", "只报告多因素 HR", "有：v1 取 ln(HR) 记 sens_hr，进 S7", "按 v1 先例",
     "33189527;37559353;37686409;39068564;40900420;41668068;41878699;35592542", "否（沿用）"),
    ("R7", "模型变量是基因集/通路/免疫细胞评分而不是基因", "有：v1 排除为“非基因级线性 signature”", "按 v1 先例排除",
     "23819917;36918939;42635880;34128342", "否（沿用）"),
    ("R8", "模板相关匹配、聚类、质心分类器等", "有：v1 排除为“非线性模型”", "按 v1 先例排除", "15067324;15931389;17183660", "否（沿用）"),
    ("R15", "原文系数与 HR/方向/正文描述自相矛盾", "有：v1 方向复核时保持原文所印数值",
     "照原文公式/表格所印数值，不自行改正；在逐篇表中标注供专家知悉",
     "35732871;38179769;41591692;35741849;33197890;32735556;37746266;41988531", "否（沿用）"),
    ("R4", "signature 不是对着生存终点建的（表型/谱系/缺氧/分化等），或机制研究里附带的多基因评分，只做了预后关联检验",
     "无（v1 的 292 篇以典型预后模型为主）",
     "纳入：只要文中把这个多基因组合评分作为预后指标、并对生存/复发做了检验。没有系数的进 S8（不影响主分析）。"
     "排除：只有单基因的预后检验、或组合评分没做生存分析",
     "见“逐篇裁定”表 R4 列", "是"),
    ("R5", "评分是 PCA 主成分", "无",
     "原文报告了载荷 → 按载荷作线性系数（方向按原文“高分=高/低风险”取正负）；未报告载荷 → 排除为“无法复现（PCA 载荷未报告）”",
     "32951005（有载荷，高分=预后好，取负）;29556354;31107158;35495154;36203594;41828396;36607288;36814131;38988627", "是"),
    ("R6", "对预定义基因集（MSigDB/GeneCards 等）直接打分，没有基因筛选或建模", "无",
     "排除为“非新 signature（预定义基因集）”", "42415249;37108678;28460430;41090793;37547288", "是"),
    ("R9", "系数作用于二分/分箱后的指示变量，或评分不是加权和（如加权表达的中位数）", "无",
     "不进主分析（SAP §5 的线性评分只定义了连续表达的加权和）；记“非线性（阈值/分箱规则）”。"
     "可选：对阈值按队列内定义的（中位数、SD、min-max 分箱）另做事后敏感性，需标明 post hoc",
     "30637711;35204174;38576998;30024105;37349324;38302875;35531101;37176106", "是"),
    ("R10", "最终模型含临床变量", "部分：v1 对无法拆分的混合模型排除",
     "基因部分有单独印出的系数 → 只评基因部分（24713434 的 24 基因 RNA score）；无法拆分 → 排除为“非转录组 signature”",
     "24713434;36110544", "是"),
    ("R13", "模型建在非肿瘤组织上", "无",
     "34631510（癌旁组织建模并验证）排除为“非肿瘤组织 signature”；从血液挑基因、但在肿瘤组织打分的（36291815、37509358）照常纳入",
     "34631510;36291815;37509358", "是"),
    ("R14", "作者把 HR（exp(B)）当权重写进公式", "无",
     "原文同时印了 B → 主分析用 B，公式权重写入备注（34804963）；只有 HR → ln(HR) 进 S7（35592542）",
     "34804963;35592542", "是"),
    ("R17", "作者对自己既往 signature 的改良版/原始推导", "部分：v1 排除商业分类器的验证研究",
     "原始推导论文算新 signature（21310658 CCP 首次提出）；本文重新筛基因的改良版算新（38302875）；只是换算法重训同一基因集的不算新（29511883）",
     "21310658;38302875;29511883", "是"),
]
rules = pd.DataFrame(RULES, columns=["编号", "规则", "v1 先例", "Claude 建议", "涉及 PMID", "需 Wei/专家决定"])
rules["Wei 决定"] = ["沿用 v1" if n.startswith("否") else "按 Claude 建议执行（Wei 2026-11-14 授权 Claude 定规则）" for n in rules["需 Wei/专家决定"]]
rules["专家意见"] = ""

# ---- per-paper sheet: pending + R4/R5/R6/R9/R10/R13/R14/R17 papers ----
rule_of = {}
for code, _, _, _, ids, _ in RULES:
    for p in re.findall(r"\d{8}", ids):
        rule_of.setdefault(p, []).append(code)
r4_rx = r"机制|表型|谱系|干性|缺氧|LVI|药效|分化状态|治疗选择|非生存|附带"
for _, r in R.iterrows():
    if r["拟定角色"] in ("equal_weight", "equal_signed", "main", "sens_hr", "sens_figure", "待裁定") and \
            re.search(r4_rx, r["X_note"] + r["Y_note"] + r["X_status"] + r["Y_status"]):
        rule_of.setdefault(r["PMID"], [])
        if "R4" not in rule_of[r["PMID"]]:
            rule_of[r["PMID"]].append("R4")
need = adj[adj["需人工"].str.startswith("是")]["PMID"].tolist()
ids = sorted(set([p for p, v in rule_of.items() if set(v) - {"R1", "R2", "R3", "R7", "R8", "R15"}] + need))
P = R[R["PMID"].isin(ids)].copy()
P["涉及规则"] = P["PMID"].map(lambda p: "、".join(rule_of.get(p, [])) or "原文图表矛盾")
P = P.merge(adj[["PMID", "裁定", "依据"]].rename(columns={"裁定": "Claude 建议", "依据": "Claude 依据"}), on="PMID", how="left").fillna("")
P["Wei/专家裁定（纳入角色或排除理由）"] = ""
P["备注"] = ""
pcols = ["PMID", "批次", "涉及规则", "拟定角色", "X_status", "Y_status", "Claude 建议", "Claude 依据", "genes", "coefs",
         "endpoint", "X_note", "Y_note", "Wei/专家裁定（纳入角色或排除理由）", "备注"]
P = P[pcols]

# ---- Claude-resolved disagreements ----
done = adj[~adj["需人工"].str.startswith("是")].copy()

# ---- missing supplements ----
SUPP = [
    ("34131070", "Supplementary Table S5", "217 基因 + elastic-net Cox 系数", "高：拿到即可进主分析"),
    ("28408174", "Supplementary Table 13", "70 个转录本 + PLS 权重", "高：拿到即可进主分析"),
    ("21289624", "Supplementary Tables 4–5", "4 基因多因素 Cox 系数", "高：拿到即可进主分析"),
    ("40216135", "Supplementary Table S2", "系数（CoxBoost+plsRcox）", "高"),
    ("39413958", "Supplementary Fig. S1E", "9 基因 LASSO 系数", "高"),
    ("34876840", "Supplementary Figure S1", "10 个 lncRNA 对的系数", "高"),
    ("37105783", "补充材料（fixed weights）", "4 基因权重", "中：不确定补充材料里有"),
    ("37758909", "Supplementary Table S6", "可能含系数", "中：不确定"),
    ("42092422", "Script S1 / Table S3", "可能含 LASSO 系数", "中：不确定"),
    ("40312179", "补充材料 / GitHub 代码", "RIDGE Cox 系数（OS 模型；主推为 GBM）", "低：主推模型为非线性"),
    ("37559353", "勘误（erratum）", "核对 HR/系数是否更正", "低"),
    ("36607288", "Table S4", "14 基因列表（取决于 R4）", "取决于 R4"),
    ("35405009", "Figure 5C 原图", "21 基因列表（取决于 R4）", "取决于 R4"),
    ("33334054", "Table S4/S5（Excel）", "C1 基质 50 基因（取决于 R4）", "取决于 R4"),
]
S = pd.DataFrame(SUPP, columns=["PMID", "需要的文件", "里面有什么", "优先级"])

intro = pd.DataFrame({"说明": [
    "v2 双人独立提取已完成：15 批 300 篇（两篇预印本以正式版计）。",
    f"两人完全一致 {(R['一致性'] == '完全一致').sum()} 篇；不一致 {(R['一致性'] != '完全一致').sum()} 篇。",
    "不一致中，能回原文核实或可直接套 v1 规则的已由 Claude 裁定（见“已裁定”表），其余进“逐篇裁定”。",
    "“规则裁定”表：R1/R2/R3/R7/R8/R15 沿用 v1 已锁定的做法，只需知悉；R4/R5/R6/R9/R10/R13/R14/R17 是 v1 没遇到过的情形，需要 Wei 和专家定。",
    "规则定下后，“逐篇裁定”表里大部分行可以按规则自动落定；表里只需写例外。",
    "“缺补充材料”表：前 6 篇拿到后可以直接进主分析，请优先下载。",
    "“全部300篇”表是按当前建议算出的拟定角色（main=主分析；sens_figure/sens_hr=S7；equal_weight/equal_signed=S8；excluded=排除）。",
]})

with pd.ExcelWriter(D / "裁定表_v2.xlsx", engine="openpyxl") as w:
    intro.to_excel(w, sheet_name="说明", index=False)
    rules.to_excel(w, sheet_name="规则裁定", index=False)
    P.to_excel(w, sheet_name="逐篇裁定", index=False)
    done.to_excel(w, sheet_name="已裁定（Claude）", index=False)
    S.to_excel(w, sheet_name="缺补充材料", index=False)
    R.to_excel(w, sheet_name="全部300篇_拟定角色", index=False)
    yellow = PatternFill("solid", fgColor="FFF2CC")
    for ws in w.book.worksheets:
        for col in ws.columns:
            ws.column_dimensions[col[0].column_letter].width = 18
        for row in ws.iter_rows():
            for c in row:
                c.alignment = Alignment(wrap_text=True, vertical="top")
        for c in ws[1]:
            c.font = Font(bold=True)
    w.book["说明"].column_dimensions["A"].width = 120
    ws = w.book["规则裁定"]
    for k, wd in zip("ABCDEFGH", (6, 40, 26, 60, 40, 12, 24, 24)):
        ws.column_dimensions[k].width = wd
    for c in list(ws["G"][1:]) + list(ws["H"][1:]):
        c.fill = yellow
    ws = w.book["逐篇裁定"]
    for k, wd in zip("ABCDEFGHIJKLMNO", (10, 6, 12, 12, 22, 22, 40, 40, 30, 20, 14, 60, 60, 30, 20)):
        ws.column_dimensions[k].width = wd
    for c in list(ws["N"][1:]) + list(ws["O"][1:]):
        c.fill = yellow
print("rules", len(rules), "per-paper", len(P), "resolved", len(done), "supp", len(S))
print(P["涉及规则"].value_counts().head(15).to_string())
