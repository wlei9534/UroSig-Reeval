import re, xml.etree.ElementTree as ET, pandas as pd
from pathlib import Path
XD = Path(r'D:\article\UroSig-Reeval\03_extraction\fulltext_xml')
x = pd.read_excel(r'D:\article\UroSig-Reeval\03_extraction\ai_prefill_PRAD.xlsx', sheet_name='逐篇', dtype=str)
m = x[x['状态'] == '未找到系数，需人工']
n = 0; mism = 0
for _, r in m.iterrows():
    root = ET.fromstring((XD / (r['PMCID'] + '.xml')).read_text(encoding='utf-8'))
    body = ''.join(root.find('.//body').itertext())
    t = root.find('.//article-title'); t = ''.join(t.itertext()) if t is not None else ''
    pm = [e.text for e in root.iter('article-id') if e.get('pub-id-type') == 'pmid']
    if pm and pm[0] != r['PMID']: mism += 1
    if not re.search(r'signature|risk model|risk score|prognostic model|classifier', body, re.I):
        n += 1
        if n <= 6: print(r['PMCID'], '| pmid in xml', pm[:1], 'csv', r['PMID'], '| body', len(body), '| xml:', t[:70], '| csv:', r['标题'][:70])
print('no-signature-words', n, '| PMID mismatches among missing', mism)
