"""Harmonise prostate-cancer clinical fields across cohorts."""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

_MISSING = {"", "na", "n/a", "nan", "unknown", "none", "-", "null"}


def to_float(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s.astype(str).str.strip().where(~s.astype(str).str.strip().str.lower().isin(_MISSING)), errors="coerce")


def parse_gleason(raw):
    """Return (primary, secondary) Gleason patterns, or None if patterns are not recoverable.

    Accepts '7=4+3', '4+3' and the three-digit GSE54460 code '347' (primary, secondary, sum).
    """
    if raw is None:
        return None
    s = str(raw).strip()
    m = re.search(r"(\d)\s*\+\s*(\d)", s)
    if m:
        return int(m.group(1)), int(m.group(2))
    m = re.fullmatch(r"(\d)(\d)(\d{1,2})", s)
    if m and int(m.group(1)) + int(m.group(2)) == int(m.group(3)):
        return int(m.group(1)), int(m.group(2))
    return None


def isup_group(primary: int, secondary: int) -> int:
    total = primary + secondary
    if total <= 6:
        return 1
    if total == 7:
        return 2 if primary == 3 else 3
    if total == 8:
        return 4
    return 5


def parse_pt(raw):
    """Pathological T category: 'T2' (incl. T1/T2), 'T3a', 'T3' (unspecified), 'T3b+' (T3b/T4); None if unknown."""
    if raw is None:
        return None
    s = str(raw).strip().lower()
    m = re.search(r"p?t\s*([1-4x])\s*([abc])?", s)
    if not m or m.group(1) == "x":
        return None
    n, sub = m.group(1), m.group(2)
    if n in ("1", "2"):
        return "T2"
    if n == "4":
        return "T3b+"
    if sub == "a":
        return "T3a"
    if sub in ("b", "c"):
        return "T3b+"
    return "T3"


def bcr_time(event, t_event, t_followup, event_values=("1",), censor_values=("0",)):
    """Assemble (time, event): time = time to BCR for events, total follow-up for censored."""
    ev = event.astype(str).str.strip()
    E = pd.Series(np.where(ev.isin(event_values), 1.0, np.where(ev.isin(censor_values), 0.0, np.nan)), index=event.index)
    te, tf = to_float(t_event), to_float(t_followup)
    T = pd.Series(np.where(E == 1, te, np.where(E == 0, tf, np.nan)), index=event.index)
    bad = T.isna()
    E[bad] = np.nan
    return T, E
