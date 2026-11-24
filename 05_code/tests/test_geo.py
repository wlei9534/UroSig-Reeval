import gzip
import pandas as pd
from urosig.geo import read_series_matrix

TXT = '''!Series_title\t"Toy"
!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"
!Sample_title\t"t1"\t"t2"\t"t3"
!Sample_characteristics_ch1\t"tissue: tumor"\t"tissue: tumor"\t"tissue: normal"
!Sample_characteristics_ch1\t"os time (months): 12.5"\t"os time (months): 30"\t""
!Sample_characteristics_ch1\t"os status: dead"\t"os status: alive"\t""
!series_matrix_table_begin
"ID_REF"\t"GSM1"\t"GSM2"\t"GSM3"
"p1"\t1.0\t2.0\t3.0
"p2"\t4.0\t5.0\tnull
!series_matrix_table_end
'''

def test_read_series_matrix_parses_expression_and_characteristics(tmp_path):
    f = tmp_path / "toy_series_matrix.txt.gz"
    with gzip.open(f, "wt") as fh:
        fh.write(TXT)
    expr, pheno = read_series_matrix(f)
    assert list(expr.columns) == ["GSM1", "GSM2", "GSM3"]
    assert expr.loc["p1", "GSM2"] == 2.0 and pd.isna(expr.loc["p2", "GSM3"])
    assert pheno.loc["GSM1", "tissue"] == "tumor"
    assert pheno.loc["GSM2", "os time (months)"] == "30"
    assert pheno.loc["GSM1", "title"] == "t1"

def test_collapse_keeps_highest_mean_probe():
    from urosig.geo import collapse_to_genes
    expr = pd.DataFrame({"S1": [1.0, 5.0, 2.0], "S2": [1.0, 7.0, 2.0]}, index=["a", "b", "c"])
    out = collapse_to_genes(expr, pd.Series({"a": "TP53", "b": "TP53", "c": "VHL"}))
    assert out.loc["TP53", "S2"] == 7.0 and set(out.index) == {"TP53", "VHL"}
