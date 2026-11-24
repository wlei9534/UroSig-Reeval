"""Gene symbol harmonisation to current HGNC approved symbols."""
from __future__ import annotations

from collections import defaultdict

import pandas as pd


def build_symbol_map(hgnc: pd.DataFrame) -> dict:
    """Map approved, previous and alias symbols to the approved symbol.

    Priority: approved > previous > alias. A previous/alias symbol pointing to more
    than one approved symbol is ambiguous and dropped.
    """
    approved = set(hgnc["symbol"].astype(str))
    m = {s: s for s in approved}
    for col in ("prev_symbol", "alias_symbol"):
        cand = defaultdict(set)
        for sym, field in zip(hgnc["symbol"].astype(str), hgnc[col].fillna("").astype(str)):
            for old in filter(None, (x.strip() for x in field.split("|"))):
                cand[old].add(sym)
        for old, targets in cand.items():
            if old in m:
                continue
            if len(targets) == 1:
                m[old] = next(iter(targets))
    return m


def harmonise_index(expr: pd.DataFrame, symbol_map: dict):
    """Rename gene index to approved symbols; on collisions keep the row with the highest mean."""
    new = expr.index.map(lambda s: symbol_map.get(str(s)))
    keep = pd.notna(new)
    e = expr[keep].copy()
    e.index = new[keep]
    n_before = len(e)
    e = e.assign(_m=e.mean(axis=1)).sort_values("_m", ascending=False)
    e = e[~e.index.duplicated(keep="first")].drop(columns="_m")
    report = {"input": len(expr), "unmapped": int((~keep).sum()), "collisions": n_before - len(e), "output": len(e)}
    return e.sort_index(), report
