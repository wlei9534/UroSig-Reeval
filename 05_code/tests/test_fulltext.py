from urosig.fulltext import formula_terms, coef_tables, supplement_hint

def test_formula_terms_coef_times_gene_with_unicode_minus():
    t = "Risk score = (0.1234 × FOXM1) + (−0.2056 × PTEN) + (0.0871 × AC009032.1)."
    assert formula_terms(t) == [("FOXM1", 0.1234), ("PTEN", -0.2056), ("AC009032.1", 0.0871)]

def test_formula_terms_gene_times_coef_and_star():
    t = "riskScore = TOP2A*0.452 + CDK1 * -0.13 - KLK3*0.5"
    assert formula_terms(t) == [("TOP2A", 0.452), ("CDK1", -0.13), ("KLK3", -0.5)]

def test_formula_terms_expression_wording():
    t = "risk score = 0.311 × expression of BIRC5 − 0.274 × expression level of NKX3-1"
    assert formula_terms(t) == [("BIRC5", 0.311), ("NKX3-1", -0.274)]

def test_formula_terms_ignores_text_without_formula():
    assert formula_terms("Patients were divided at the median risk score of 1.2 (P < 0.05).") == []

JATS = """<article><body><table-wrap id="t2"><label>Table 2</label><caption><p>Genes and coefficients of the signature</p></caption>
<table><thead><tr><th>Gene</th><th>Coefficient</th><th>HR</th></tr></thead>
<tbody><tr><td>MKI67</td><td>0.321</td><td>1.4</td></tr><tr><td>LINC00261</td><td>−0.112</td><td>0.8</td></tr><tr><td>Total</td><td></td><td></td></tr></tbody></table></table-wrap>
<table-wrap id="t1"><label>Table 1</label><caption><p>Baseline characteristics</p></caption><table><thead><tr><th>Variable</th><th>n</th></tr></thead><tbody><tr><td>Age</td><td>62</td></tr></tbody></table></table-wrap></body></article>"""

def test_coef_tables_reads_gene_and_coefficient_columns_only_from_coef_tables():
    out = coef_tables(JATS)
    assert len(out) == 1
    assert out[0]["label"] == "Table 2"
    assert out[0]["terms"] == [("MKI67", 0.321), ("LINC00261", -0.112)]

def test_supplement_hint():
    assert supplement_hint("The coefficients of each gene are listed in Supplementary Table S3.")
    assert not supplement_hint("Supplementary Table S1 lists patient characteristics.")

def test_formula_terms_tolerates_thin_and_nbsp_spaces():
    t = "Risk score = (+ 0.25 × AR) + (− 0.5 × FOXA1)"
    assert formula_terms(t) == [("AR", 0.25), ("FOXA1", -0.5)]

def test_formula_terms_never_raises_on_odd_numbers():
    formula_terms("risk score = (+  2 × AB1) + (3e × CD2) + (0.1 × EF3)")

def test_formula_terms_exp_prefix_glued_to_gene():
    t = "Risk score = 0.27761 × ExpNELL2 − 0.39196 × ExpUNC5A + 0.23357 × ExpTNC."
    assert formula_terms(t) == [("NELL2", 0.27761), ("UNC5A", -0.39196), ("TNC", 0.23357)]

def test_formula_terms_rejects_placeholder_symbols():
    assert formula_terms("risk score = 1 × X1 + 2 × X2 + 3 × X3") == []
    assert formula_terms("Risk score = Coef gene1 × expression gene1 + Coef gene2 × expression gene2") == []

def test_gene_symbol_not_glued_to_negative_coefficient():
    t = "Risk score = FGD2-0.3572187 + VAV1-0.3522158 + NKX3-1*0.2 + RP11-497G19.3*0.5"
    genes = [g for g, _ in formula_terms(t)]
    assert "FGD2-0" not in genes and "NKX3-1" in genes and "RP11-497G19.3" in genes

def test_consecutive_integer_coefficients_are_rejected():
    assert formula_terms("risk score = 1 × BUB1 + 2 × PTGS2 + 3 × RAC3 + 4 × IRX1") == []

def test_formula_from_cited_study_is_flagged():
    from urosig.fulltext import cited_formula
    assert cited_formula("Liu et al. reported a 13-gene signature, with the formula: Risk score = 0.2 × A1 + 0.3 × B2")
    assert not cited_formula("The risk score was calculated as follows: Risk score = 0.2 × A1 + 0.3 × B2")

def test_univariate_table_is_flagged():
    from urosig.fulltext import coef_tables
    j = JATS.replace("Genes and coefficients of the signature", "Univariable and multivariate Cox regression analysis")
    assert coef_tables(j)[0]["mixed_uni_multi"] is True

def test_gene_regex_rejects_numeric_suffixes():
    import re
    from urosig.fulltext import GENE
    assert not re.fullmatch(GENE, "FGD2-0.3572187")
    for ok in ["NKX3-1", "AC009032.1", "RP11-497G19.3", "HLA-DRA", "MIR100HG", "ZEB1-AS1"]:
        assert re.fullmatch(GENE, ok), ok
