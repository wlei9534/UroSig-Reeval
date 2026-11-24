import re, sys, collections
import xml.etree.ElementTree as ET
from pathlib import Path
import pandas as pd
ROOT = Path(r'D:\article\UroSig-Reeval'); XD = ROOT / '03_extraction' / 'fulltext_xml'
x = pd.read_excel(ROOT / '03_extraction' / 'ai_prefill_PRAD.xlsx', sheet_name='逐篇', dtype=str)
miss = x[x['状态'] == '未找到系数，需人工']
c = collections.Counter(); ex = collections.defaultdict(list)
for _, r in miss.iterrows():
    j = (XD / f"{r['PMCID']}.xml").read_text(encoding='utf-8'); root = ET.fromstring(j)
    body = ''.join(root.find('.//body').itertext())
    tabs = [' | '.join(''.join(th.itertext()) for th in tw.iter('th'))[:150] for tw in root.iter('table-wrap')]
    flags = {
      'coef_word_in_body': bool(re.search(r'coefficient', body, re.I)),
      'formula_sum_generic': bool(re.search(r'(Σ|∑|sum of)', body)),
      'score_eq': bool(re.search(r'(risk ?score|riskscore)\s*[=:]', body, re.I)),
      'supp_mention_near_coef': bool(re.search(r'coefficient[^.]{0,150}(supplement|additional file|table s\d)', body, re.I) or re.search(r'(supplement|additional file|table s\d)[^.]{0,150}coefficient', body, re.I)),
      'table_with_coef_header': any(re.search(r'coef|β|beta', t, re.I) for t in tabs),
      'any_table': len(tabs) > 0,
      'not_a_signature_paper': not re.search(r'signature|risk model|risk score|prognostic model|classifier', body, re.I),
      'gene_pair': bool(re.search(r'gene pair', body, re.I)),
      'ml_model': bool(re.search(r'random survival forest|xgboost|deep learning|neural network|support vector|gbm|coxboost', body, re.I)),
    }
    for k, v in flags.items():
        if v: c[k] += 1; ex[k].append(r['PMCID'])
    if flags['score_eq']:
        m = re.search(r'(risk ?score|riskscore)\s*[=:][^\n]{0,220}', body, re.I); ex['score_eq_text'].append(m.group(0)[:220] if m else '')
    if flags['table_with_coef_header']:
        ex['coef_tables'].append([t for t in tabs if re.search(r'coef|β|beta', t, re.I)][0])
print('n missing', len(miss)); [print(f'  {k}: {v}') for k, v in c.most_common()]
print('--- score= snippets'); [print('  ', s) for s in ex['score_eq_text'][:12]]
print('--- coef table headers'); [print('  ', s) for s in ex['coef_tables'][:12]]
