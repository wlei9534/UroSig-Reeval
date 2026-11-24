"""Re-derive each record's own DOI / PMCID from PubMed esummary JSON (fixes reference-list contamination).
esummary JSON is downloaded beforehand with curl into 01_search/esummary_raw/."""
import json
from pathlib import Path
import pandas as pd
ROOT = Path(r'D:\article\UroSig-Reeval')
d = pd.read_csv(ROOT / '01_search' / 'pubmed_dedup_2026-10-26.csv', dtype=str)
got = {}
for f in (ROOT / '01_search' / 'esummary_raw').glob('es_*.json'):
    j = json.loads(f.read_text(encoding='utf-8'))
    for u in j['result']['uids']:
        aid = {x['idtype']: x['value'] for x in j['result'][u].get('articleids', [])}
        pmc = aid.get('pmc', '') or aid.get('pmcid', '').replace('pmc-id: ', '').split(';')[0]
        got[u] = {'doi_ok': aid.get('doi', ''), 'pmcid_ok': pmc}
fix = pd.DataFrame.from_dict(got, orient='index'); fix.index.name = 'pmid'
d = d.merge(fix, left_on='pmid', right_index=True, how='left')
chg_doi = int((d['doi'].fillna('') != d['doi_ok'].fillna('')).sum()); chg_pmc = int((d['pmcid'].fillna('') != d['pmcid_ok'].fillna('')).sum())
d['doi'] = d['doi_ok']; d['pmcid'] = d['pmcid_ok'].replace('', pd.NA)
d = d.drop(columns=['doi_ok', 'pmcid_ok'])
d.to_csv(ROOT / '01_search' / 'pubmed_dedup_2026-10-26_idfixed.csv', index=False)
d[['pmid', 'doi', 'pmcid']].to_csv(ROOT / '01_search' / 'article_ids_esummary_2026-10-26.csv', index=False)
print(f'records {len(d)} | esummary hits {len(got)} | DOI changed {chg_doi} | PMCID changed {chg_pmc} | PMCID present {int(d.pmcid.notna().sum())} | DOI present {int((d.doi.fillna("") != "").sum())}')
