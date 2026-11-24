import gzip, sys; sys.path.insert(0, r'D:\article\UroSig-Reeval\05_code')
import pandas as pd
from urosig.geo import read_series_matrix
R = r'D:\article\UroSig-Reeval\04_cohorts\raw\GSE54460'
_, ph = read_series_matrix(R + r'\GSE54460_series_matrix.txt.gz')
with gzip.open(R + r'\GSE54460_FPKM-genes-TopHat2-106samples-12-4-13.txt.gz', 'rt') as f:
    hdr = [next(f).rstrip('\n').split('\t') for _ in range(14)]
meta = pd.DataFrame({r[0].strip(): r[2:] for r in hdr[1:14] if r[0].strip()}, index=hdr[0][2:])
ext = meta['External ID'].str.strip(); key = ph['external id'].str.strip()
un = [e for e in ext if not ((key == e) | key.str.endswith(e)).any()]
print('unmatched FPKM ext ids:', un)
print('GEO ext ids not hit:', [k for k in key if not any(k == e or k.endswith(e) for e in ext)])
# fingerprint by column order
m = meta.reset_index(drop=True); p = ph.reset_index()
f = lambda s: pd.to_numeric(s.astype(str).str.replace(',', ''), errors='coerce')
cmp = pd.DataFrame({'fu_fpkm': f(m['Months total F/U']), 'fu_geo': f(p['months total f/u']), 'bcr_fpkm': m['BCR'], 'bcr_geo': p['bcr'],
                    'psa_fpkm': f(m['prepsa']), 'psa_geo': f(p['prepsa']), 'reads_fpkm': f(m['no. of uniquely mapped reads']), 'reads_geo': f(p['no of uniquely mapped reads'])})
ok_reads = (cmp.reads_fpkm == cmp.reads_geo)
print('column-order agreement: reads', ok_reads.sum(), '/', len(cmp), '| bcr', (cmp.bcr_fpkm.str.strip() == cmp.bcr_geo.str.strip()).sum(),
      '| f/u within 0.5', ((cmp.fu_fpkm - cmp.fu_geo).abs() < 0.5).sum())
# exact match on uniquely mapped reads (unique fingerprint)
rg = f(p.set_index('index')['no of uniquely mapped reads'])
print('reads unique in GEO:', rg.is_unique, '| FPKM cols matched by reads:', f(m['no. of uniquely mapped reads']).isin(rg).sum())
print(cmp[~ok_reads].head(12).to_string())
