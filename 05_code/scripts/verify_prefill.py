"""Check every pre-extracted gene and coefficient appears literally in the source full text; write results back."""
import re, xml.etree.ElementTree as ET
from pathlib import Path
import pandas as pd
ROOT = Path(r'D:\article\UroSig-Reeval'); XD = ROOT / '03_extraction' / 'fulltext_xml'; F = ROOT / '03_extraction' / 'ai_prefill_PRAD.xlsx'
sheets = pd.read_excel(F, sheet_name=None, dtype=str)
x = sheets['逐篇']
norm = lambda s: s.translate(str.maketrans({'\u2212': '-', '\u2013': '-', '\u00a0': ' ', '\u2009': ' ', '\u202f': ' '}))
res = []
for i, r in x.iterrows():
    if r['状态'] != '已预提取，待人工核对':
        res.append(''); continue
    root = ET.fromstring((XD / f"{r['PMCID']}.xml").read_text(encoding='utf-8'))
    txt = norm(''.join(root.find('.//body').itertext()))
    compact = re.sub(r'\s', '', txt)
    genes, coefs = r['基因'].split(';'), r['系数'].split(';')
    miss_g = [g for g in genes if g not in txt]
    def coef_found(c):
        v = float(c)
        cands = {c, c.lstrip('-'), f'{abs(v):g}', f'{abs(v)}'}
        for dp in range(1, 8):
            cands.add(f'{abs(v):.{dp}f}')
        return any(k in compact for k in cands if k)
    miss_c = [c for c in coefs if not coef_found(c)]
    res.append('通过' if not miss_g and not miss_c else f"未通过：基因缺 {miss_g[:5]}；系数缺 {miss_c[:5]}")
x['逐字核对'] = res
cols = list(x.columns); cols.insert(cols.index('状态') + 1, cols.pop(cols.index('逐字核对'))); x = x[cols]
sheets['逐篇'] = x
with pd.ExcelWriter(F, engine='openpyxl') as w:
    for k, v in sheets.items():
        v.to_excel(w, sheet_name=k, index=False)
y = x[x['状态'] == '已预提取，待人工核对']
print('checked', len(y), '| pass', int((y['逐字核对'] == '通过').sum()))
for _, r in y[y['逐字核对'] != '通过'].iterrows():
    print(' ', r['PMID'], r['逐字核对'][:160].encode('ascii', 'backslashreplace').decode())
