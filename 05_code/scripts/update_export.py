"""Bring an earlier Word export of the online manuscript up to date with the current text of the online document.
Usage: python update_export.py <old_export.docx> <doc.xml> <new_export.docx>
Styles, numbering and tables of the export are kept; only paragraphs whose text changed are rewritten."""
import copy, sys
import xml.etree.ElementTree as ET
from docx import Document
from docx.oxml.ns import qn
from docx.text.paragraph import Paragraph
from docx.table import Table

old, xmlf, new = sys.argv[1:4]
doc = Document(old); root = ET.parse(xmlf).getroot()
body = [c for c in doc.element.body if c.tag in (qn("w:p"), qn("w:tbl"))]
pos = 0; changed = []; added = []

def runs_of(par):
    out = []
    for t in par:
        if t.tag == "text":
            b = t.find("bold")
            out.append((b.text if b is not None else (t.text or ""), b is not None))
        elif t.tag == "date": out.append((None, False))      # keep the export's own rendering
        elif t.tag == "mention": out.append((None, False))
    return out

def set_par(p_el, runs):
    par = Paragraph(p_el, None)
    for r in list(p_el.findall(qn("w:r"))): p_el.remove(r)
    for txt, bold in runs:
        r = par.add_run(txt); r.bold = True if bold else None

def do_par(par_xml, list_item=False):
    global pos
    runs = runs_of(par_xml)
    if any(t is None for t, _ in runs):  # date / mention line: leave as exported
        pos += 1; return
    text = "".join(t for t, _ in runs)
    p_el = body[pos]
    assert p_el.tag == qn("w:p"), (pos, text[:40])
    cur = Paragraph(p_el, None).text
    if list_item and cur != text and text.startswith("Manjang K"):
        new_el = copy.deepcopy(body[pos - 1]); body[pos - 1].addnext(new_el); body.insert(pos, new_el)
        set_par(new_el, runs); added.append(text[:50]); pos += 1; return
    if cur != text:
        set_par(p_el, runs); changed.append((cur[:45], text[:45]))
    pos += 1

for blk in root:
    if blk.tag == "paragraph": do_par(blk)
    elif blk.tag == "list":
        for li in blk:
            for par in li: do_par(par, True)
    elif blk.tag == "table":
        t_el = body[pos]; assert t_el.tag == qn("w:tbl"); tbl = Table(t_el, doc); pos += 1
        rows = list(blk)
        assert len(rows) == len(tbl.rows), (len(rows), len(tbl.rows))
        for rx, rw in zip(rows, tbl.rows):
            cx = list(rx); assert len(cx) == len(rw.cells)
            for c_xml, c_w in zip(cx, rw.cells):
                pars = list(c_xml); assert len(pars) == len(c_w.paragraphs)
                for px, pw in zip(pars, c_w.paragraphs):
                    runs = runs_of(px); text = "".join(t for t, _ in runs)
                    if pw.text != text:
                        hb = any(r.bold for r in pw.runs)
                        set_par(pw._p, [(t, b or hb) for t, b in runs]); changed.append((pw.text[:30], text[:30]))
assert pos == len(body), (pos, len(body))
doc.save(new)
print("paragraphs rewritten:", len(changed), "| added:", added)
for a, b in changed: print("  -", a, "=>", b)
