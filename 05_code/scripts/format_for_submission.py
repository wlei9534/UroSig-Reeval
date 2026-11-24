"""Format the Word export of the manuscript for submission to BMC Medicine.

Usage: python format_for_submission.py <exported.docx> <output.docx>

Applies what the journal asks of the main text file: A4, double line spacing, continuous line numbers,
page numbers, no colour or shading in tables (tables are set in 8 or 9 pt, single spaced). Removes the working header of the online document
(its title line and the date line) and reports any leftover working text in square brackets.
The wording of the manuscript is not touched.
"""
import re
import sys

from docx import Document
from docx.enum.text import WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt

src, dst = sys.argv[1], sys.argv[2]
doc = Document(src)

# working header of the online document: first heading (document name) and the date / owner line under it
body = doc.element.body
paras = doc.paragraphs
removed = []
if paras and paras[0].text.strip().startswith("UroSig-Reeval manuscript"):
    removed.append(paras[0].text.strip()); body.remove(paras[0]._p)
    if re.match(r"^\s*(\d{4}-\d{2}-\d{2}|[A-Z][a-z]+ \d{1,2}, \d{4}|\d{1,2} [A-Z][a-z]+ \d{4})", paras[1].text.strip()):
        removed.append(paras[1].text.strip()); body.remove(paras[1]._p)

# page set-up, line numbers, page numbers
for s in doc.sections:
    s.page_width, s.page_height = Cm(21.0), Cm(29.7)
    s.left_margin = s.right_margin = s.top_margin = s.bottom_margin = Cm(2.5)
    sp = s._sectPr
    for old in sp.findall(qn("w:lnNumType")):
        sp.remove(old)
    ln = OxmlElement("w:lnNumType"); ln.set(qn("w:countBy"), "1"); ln.set(qn("w:restart"), "continuous")
    # schema order: lnNumType comes after pgMar and before pgNumType / cols
    anchor = sp.find(qn("w:pgMar"))
    (anchor.addnext(ln) if anchor is not None else sp.append(ln))
    fp = s.footer.paragraphs[0] if s.footer.paragraphs else s.footer.add_paragraph()
    for r in list(fp.runs):
        r._r.getparent().remove(r._r)
    fp.alignment = 1
    run = fp.add_run()
    for tag, text in (("begin", None), (None, " PAGE "), ("end", None)):
        if tag:
            el = OxmlElement("w:fldChar"); el.set(qn("w:fldCharType"), tag)
        else:
            el = OxmlElement("w:instrText"); el.set(qn("xml:space"), "preserve"); el.text = text
        run._r.append(el)

# font and double spacing everywhere, including inside tables
normal = doc.styles["Normal"]
normal.font.name = "Times New Roman"; normal.font.size = Pt(12)
normal.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")


def double(p):
    p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE


for p in doc.paragraphs:
    double(p)
n_shade = 0
TOTAL = int((21.0 - 5.0) / 2.54 * 1440)  # text width in twips
for t in doc.tables:
    # tables: 9 pt, single spacing, column widths in proportion to their content, so that words are not broken
    ncol = len(t.columns); need = [0] * ncol; length = [1] * ncol
    size = 8 if ncol > 7 else 9  # points
    for row in t.rows:
        for i, c in enumerate(row.cells[:ncol]):
            txt = c.text.strip()
            need[i] = max(need[i], max([len(w) for w in txt.split()] + [1]) * 11 * size + 140)  # longest word, plus cell padding
            length[i] = max(length[i], len(txt))
    spare = TOTAL - sum(need)
    if spare >= 0:  # share the spare width by the amount of text in each column
        widths = [n + int(spare * l / sum(length)) for n, l in zip(need, length)]
    else:
        widths = [int(TOTAL * n / sum(need)) for n in need]
    tblPr = t._tbl.tblPr
    for old in tblPr.findall(qn("w:tblCellMar")):
        tblPr.remove(old)
    mar = OxmlElement("w:tblCellMar")
    for side in ("left", "right"):
        el = OxmlElement("w:" + side); el.set(qn("w:w"), "60"); el.set(qn("w:type"), "dxa"); mar.append(el)
    look = tblPr.find(qn("w:tblLook"))
    (look.addprevious(mar) if look is not None else tblPr.append(mar))
    lay = tblPr.find(qn("w:tblLayout"))
    if lay is not None:
        lay.set(qn("w:type"), "fixed")
    for g, w in zip(t._tbl.tblGrid.findall(qn("w:gridCol")), widths):
        g.set(qn("w:w"), str(w))
    for row in t.rows:
        for i, c in enumerate(row.cells[:ncol]):
            tcPr = c._tc.get_or_add_tcPr()
            tcW = tcPr.find(qn("w:tcW"))
            if tcW is None:
                tcW = OxmlElement("w:tcW"); tcPr.insert(0, tcW)
            tcW.set(qn("w:w"), str(widths[i])); tcW.set(qn("w:type"), "dxa")
            for shd in tcPr.findall(qn("w:shd")):
                tcPr.remove(shd); n_shade += 1
            for p in c.paragraphs:
                p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
                p.paragraph_format.space_after = Pt(2); p.paragraph_format.space_before = Pt(2)
                p.alignment = 0
                for r in p.runs:
                    r.font.size = Pt(size)
doc.save(dst)

text = "\n".join(p.text for p in doc.paragraphs)
left = sorted(set(re.findall(r"\[[^\]\d][^\]]{8,}\]", text)))
words = lambda s: len(s.split())  # noqa: E731
m = re.search(r"Background:(.*?)Keywords:", text, re.S)
print("removed header lines:", removed)
print("paragraphs", len(doc.paragraphs), "tables", len(doc.tables), "shaded cells cleared", n_shade)
print("abstract words:", words(m.group(1)) + 1 if m else "abstract not found")
print("working text still in square brackets:", left if left else "none")
for lab in ("Working title", "Alternative"):
    if lab in text:
        print("still present on the title page:", lab)
