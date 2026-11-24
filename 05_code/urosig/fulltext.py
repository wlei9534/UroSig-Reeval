"""Pre-extraction of signature genes and coefficients from JATS full text (third-reader aid, not a substitute for manual extraction)."""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET

GENE = r"[A-Z][A-Z0-9]*(?:-(?!\d+\.\d)[A-Z0-9]+|\.\d{1,2}(?!\d))*"
NUM = r"\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"
_MINUS = str.maketrans({"−": "-", "–": "-", "—": "-", "×": "*", "∗": "*", "⋅": "*", "·": "*", " ": " "})
_FILLER = r"(?:\s*(?:the\s+)?(?:relative\s+)?(?:expression(?:\s+(?:level|value))?(?:\s+of)?\s|exp\s|Exp(?=[A-Z0-9]))\s*)?"
_STOP = {"HR", "CI", "OS", "BCR", "DFS", "RFS", "PFS", "TCGA", "GEO", "PSA", "RNA", "DNA", "AUC", "ROC", "COX", "LASSO", "SD", "P", "N", "T", "M"}

_SCORE_CTX = re.compile(r"(risk\s*score|riskscore|risk\s+value|prognostic\s+index|score)\s*[=:]", re.I)
_TERM_CG = re.compile(rf"([+-]?)\s*\(?\s*([+-]?\s*{NUM})\s*\*\s*{_FILLER}\(?\s*({GENE})\s*\)?\s*\)?")
_TERM_GC = re.compile(rf"([+-]?)\s*\(?\s*{_FILLER}({GENE})\s*\)?\s*\*\s*\(?\s*([+-]?\s*{NUM})\s*\)?")


def _num(sign: str, s: str):
    try:
        v = float(re.sub(r"\s", "", s))
    except ValueError:
        return None
    return -v if sign == "-" else v


def formula_terms(text: str):
    """(gene, coefficient) pairs from a written risk-score formula; [] if no formula is present."""
    t = text.translate(_MINUS)
    m = _SCORE_CTX.search(t)
    if not m:
        return []
    seg = t[m.end():m.end() + 4000]
    seg = re.split(r"(?<![\d.])\.\s+(?=[A-Z][a-z])|;\s|\n", seg)[0]
    terms = []
    for rx, order in ((_TERM_CG, "cg"), (_TERM_GC, "gc")):
        found = []
        for mm in rx.finditer(seg):
            if order == "cg":
                sign, c, g = mm.group(1), mm.group(2), mm.group(3)
            else:
                sign, g, c = mm.group(1), mm.group(2), mm.group(3)
            if g in _STOP or len(g) < 2 or re.fullmatch(r"X\d+|GENE\d*|MRNA\d*", g, re.I):
                continue
            v = _num(sign, c)
            if v is not None:
                found.append((g, v))
        if len(found) > len(terms):
            terms = found
    seen, out = set(), []
    for g, c in terms:
        if g not in seen:
            seen.add(g)
            out.append((g, c))
    if len(out) >= 2 and [c for _, c in out] == [float(i) for i in range(1, len(out) + 1)]:
        return []  # numbering, not coefficients
    return out if len(out) >= 2 else []


def cited_formula(text: str) -> bool:
    """True if the formula text appears to quote another study's signature."""
    head = text.translate(_MINUS)
    m = _SCORE_CTX.search(head)
    ctx = head[max(0, (m.start() if m else 0) - 300):(m.start() if m else 300)]
    return bool(re.search(r"et al\.?|previously (?:reported|published|developed)|reported by|proposed by|published signature", ctx, re.I))


def _text(el) -> str:
    return "".join(el.itertext()).strip() if el is not None else ""


_COEF_HDR = re.compile(r"coef|β|beta|\bb\b|weight", re.I)


def coef_tables(jats: str):
    """Tables whose header names a coefficient column; returns [{'label','caption','terms'}]."""
    try:
        root = ET.fromstring(jats)
    except ET.ParseError:
        return []
    out = []
    for tw in root.iter("table-wrap"):
        tab = tw.find(".//table")
        if tab is None:
            continue
        rows = [[_text(c).translate(_MINUS) for c in tr if c.tag in ("td", "th")] for tr in tab.iter("tr")]
        rows = [r for r in rows if r]
        if len(rows) < 3:
            continue
        hdr = rows[0]
        ci = next((i for i, h in enumerate(hdr) if _COEF_HDR.search(h)), None)
        if ci is None:
            continue
        gi = 0 if ci != 0 else 1
        terms = []
        for r in rows[1:]:
            if len(r) <= max(ci, gi):
                continue
            g, c = r[gi].strip(), r[ci].strip()
            if re.fullmatch(GENE, g) and g not in _STOP and re.fullmatch(rf"[+-]?\s*{NUM}", c):
                v = _num("", c)
                if v is not None:
                    terms.append((g, v))
        if len(terms) >= 2:
            cap = _text(tw.find("caption"))
            mixed = bool(re.search(r"univariat|univariabl", cap + " " + " ".join(hdr), re.I))
            out.append({"label": _text(tw.find("label")), "caption": cap[:200], "terms": terms, "mixed_uni_multi": mixed})
    return out


def supplement_hint(text: str) -> bool:
    """True if the text says coefficients are given in supplementary material."""
    for m in re.finditer(r"supplement\w*\s+(?:table|file|data|material)\w*\s*\w*", text, re.I):
        window = text[max(0, m.start() - 200): m.end() + 100]
        if re.search(r"coefficient|β|beta|weight", window, re.I):
            return True
    return False
