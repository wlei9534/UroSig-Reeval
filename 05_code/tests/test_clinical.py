import numpy as np
import pandas as pd
import pytest
from urosig.clinical import parse_gleason, isup_group, parse_pt, bcr_time, to_float

@pytest.mark.parametrize("raw,exp", [("7=4+3", (4, 3)), ("6=3+3", (3, 3)), ("10=5+5", (5, 5)),
                                     ("347", (3, 4)), ("437", (4, 3)), ("549", (5, 4)), ("4+3", (4, 3))])
def test_parse_gleason_pattern_forms(raw, exp):
    assert parse_gleason(raw) == exp

@pytest.mark.parametrize("raw", ["unknown", "N/A", "", None, "7"])
def test_parse_gleason_returns_none_without_patterns(raw):
    assert parse_gleason(raw) is None

@pytest.mark.parametrize("p,s,isup", [(3, 3, 1), (3, 2, 1), (2, 4, 1), (3, 4, 2), (4, 3, 3), (4, 4, 4), (3, 5, 4), (5, 3, 4), (4, 5, 5), (5, 5, 5)])
def test_isup_from_patterns(p, s, isup):
    assert isup_group(p, s) == isup

@pytest.mark.parametrize("raw,exp", [("pT3a N0Mx", "T3a"), ("T2bMxNx", "T2"), ("pT2C", "T2"), ("pT1C", "T2"),
                                     ("T3bM1bNx", "T3b+"), ("pT4 N0Mx", "T3b+"), ("pT3", "T3"), ("TxM0N0", None),
                                     ("UNKNOWN", None), ("T3", "T3")])
def test_parse_pt_categories(raw, exp):
    assert parse_pt(raw) == exp

def test_bcr_time_uses_event_time_for_events_and_followup_otherwise():
    ev = pd.Series(["Y", "N", "N/A", "Y"])
    t_bcr = pd.Series(["12", "N/A", "N/A", "NA"])
    t_fu = pd.Series(["60", "65", "70", "50"])
    T, E = bcr_time(ev, t_bcr, t_fu, event_values=("Y",), censor_values=("N",))
    assert list(E[:2]) == [1, 0] and np.isnan(E[2]) and np.isnan(E[3])
    assert T[0] == 12 and T[1] == 65 and np.isnan(T[2])

def test_to_float_handles_placeholders():
    out = to_float(pd.Series(["7.2", "NA", "UNKNOWN", "", "10"]))
    assert out[0] == 7.2 and out[4] == 10 and out[1:4].isna().all()
