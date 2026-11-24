import pandas as pd
from urosig.genes import build_symbol_map, harmonise_index

HGNC = pd.DataFrame({"symbol": ["CCDC175", "TP53", "AMER1", "MARCHF1"],
                     "prev_symbol": ["C14orf38", "", "FAM123B|WTX", "MARCH1"],
                     "alias_symbol": ["", "LFS1", "", ""]})

def test_symbol_map_resolves_previous_and_alias():
    m = build_symbol_map(HGNC)
    assert m["C14orf38"] == "CCDC175" and m["WTX"] == "AMER1" and m["LFS1"] == "TP53" and m["TP53"] == "TP53"

def test_symbol_map_drops_ambiguous_alias():
    h = pd.concat([HGNC, pd.DataFrame({"symbol": ["GENEX"], "prev_symbol": [""], "alias_symbol": ["LFS1"]})])
    assert "LFS1" not in build_symbol_map(h)

def test_harmonise_index_renames_and_keeps_highest_mean_on_collision():
    m = build_symbol_map(HGNC)
    e = pd.DataFrame({"S1": [1.0, 9.0, 3.0], "S2": [1.0, 9.0, 3.0]}, index=["C14orf38", "CCDC175", "UNKNOWNGENE"])
    out, report = harmonise_index(e, m)
    assert list(out.index) == ["CCDC175"] and out.loc["CCDC175", "S1"] == 9.0
    assert report["unmapped"] == 1 and report["collisions"] == 1
