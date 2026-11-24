"""Enriched per-paper reading digests (v2) for manual extraction.
Usage: python build_digests_v2.py <project_root>
"""
import re, sys, io, csv, zipfile, subprocess, xml.etree.ElementTree as ET
from pathlib import Path
import pandas as pd
R = Path(sys.argv[1]) / '03_extraction'; OUT = R / 'digests2'; OUT.mkdir(exist_ok=True)
x = pd.read_excel(R / 'ai_prefill_PRAD.xlsx', sheet_name='逐篇', dtype=str)
todo = x[x['状态'].isin({'未找到系数，需人工', '系数可能在补充材料'})]
G = re.compile(r'^[A-Z][A-Z0-9]*(?:[-.][A-Z0-9]+)*$')
def T(e): return re.sub(r'\s+', ' ', ''.join(e.itertext())).strip() if e is not None else ''
def rows_md(rows, maxr=30, maxc=8):
    rows = [[str(c).replace('|', '/').strip()[:22] for c in r[:maxc]] for r in rows if any(str(c).strip() for c in r)]
    return '\n'.join('    ' + ' | '.join(r) for r in rows[:maxr]) + (f'\n    … ({len(rows)} rows total)' if len(rows) > maxr else '')
def gene_rows(rows): return sum(1 for r in rows if any(G.match(str(c).strip()) for c in r[:3]))
SENT = re.compile(r'(\b\d{1,2}|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve|thirteen|fourteen|fifteen)[- ](gene|lncRNA|mRNA|IRGP|gene pair)s?\b|signature (?:consist|compos|includ|contain)|following (?:genes|formula)|genes? (?:were|was) (?:finally )?(?:identified|selected|retained|screened)|risk ?score\s*[=:(]|coefficient|LASSO|multivariate Cox|β', re.I)
for _, r in todo.iterrows():
    p = r['PMCID']; root = ET.fromstring((R / 'fulltext_xml' / f'{p}.xml').read_text(encoding='utf-8'))
    body = T(root.find('.//body'))
    out = [f"# {r['PMID']} {p} | {r['标题']} | {r['期刊']} {r['年份']} | {r['状态']}"]
    sents = re.split(r'(?<=[.;])\s+(?=[A-Z(])', body)
    hits = [s for s in sents if SENT.search(s)]
    out.append('## 关键句'); out += ['- ' + s[:500] for s in hits[:18]]
    out.append('## 表格')
    for tw in root.iter('table-wrap'):
        tab = tw.find('.//table'); cap = T(tw.find('caption'))
        rows = [[T(c) for c in tr if c.tag in ('td', 'th')] for tr in tab.iter('tr')] if tab is not None else []
        if gene_rows(rows) >= 2 or re.search(r'coef|lasso|signature|risk|multivariate', cap, re.I):
            out.append(f"**{T(tw.find('label'))}** {cap[:150]}\n" + rows_md(rows, 40))
    out.append('## 图注')
    for fg in root.iter('fig'):
        g = fg.find('.//graphic'); href = g.get('{http://www.w3.org/1999/xlink}href') if g is not None else ''
        out.append(f"- {T(fg.find('label'))} [{href}] {T(fg.find('caption'))[:260]}")
    out.append('## 补充材料')
    d = R / 'fulltext' / p
    for f in sorted(d.iterdir()) if d.exists() else []:
        s = f.suffix.lower()
        try:
            if s in ('.xlsx', '.xls', '.xlsm'):
                for n, df in pd.read_excel(f, sheet_name=None, header=None, dtype=str).items():
                    rows = df.fillna('').values.tolist()
                    if gene_rows(rows) >= 2: out.append(f'**{f.name}/{n}** ({len(rows)} rows)\n' + rows_md(rows))
                    else: out.append(f'- {f.name}/{n}: {len(rows)} rows, no gene-like column; header: ' + ' | '.join(str(c)[:15] for c in (rows[0] if rows else []))[:150])
            elif s in ('.csv', '.tsv'):
                rows = list(csv.reader(io.StringIO(f.read_text(errors='replace')), delimiter='\t' if s == '.tsv' else ','))
                out.append(f'**{f.name}** ({len(rows)} rows)\n' + rows_md(rows) if gene_rows(rows) >= 2 else f'- {f.name}: {len(rows)} rows, no gene-like column')
            elif s == '.docx':
                ns = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
                X = ET.fromstring(zipfile.ZipFile(f).read('word/document.xml'))
                txt = ' '.join(''.join(t.text or '' for t in pp.iter(ns + 't')) for pp in X.iter(ns + 'p'))
                tb = list(X.iter(ns + 'tbl'))
                out.append(f'- {f.name}: {len(tb)} tables; text: {txt[:300]}')
                for k, t in enumerate(tb, 1):
                    rows = [[''.join(n.text or '' for n in tc.iter(ns + 't')) for tc in tr.iter(ns + 'tc')] for tr in t.iter(ns + 'tr')]
                    if gene_rows(rows) >= 2: out.append(f'**{f.name}/tbl{k}**\n' + rows_md(rows))
            elif s == '.pdf':
                t = subprocess.run(['pdftotext', '-layout', str(f), '-'], capture_output=True, text=True).stdout
                ls = [l for l in t.splitlines() if re.search(r'coef|β|beta|lasso', l, re.I) or re.match(r'\s*[A-Z][A-Z0-9-]{1,12}\s+[-−]?\d+\.\d+', l)]
                out.append(f'- {f.name} (pdf, {len(t)} chars){": coef-like lines:" if ls else ""}\n' + '\n'.join('    ' + l.strip()[:120] for l in ls[:30]))
            elif s == '.zip':
                out.append(f'- {f.name} zip: ' + ', '.join(zipfile.ZipFile(f).namelist()[:15]))
            elif s in ('.doc',):
                out.append(f'- {f.name}: legacy .doc')
            else:
                out.append(f'- {f.name}')
        except Exception as e:
            out.append(f'- {f.name}: unreadable ({type(e).__name__})')
    (OUT / f"{r['PMID']}.md").write_text('\n'.join(out), encoding='utf-8')
print(len(list(OUT.glob('*.md'))))
