import gzip, sys; sys.path.insert(0, r'D:\article\UroSig-Reeval\05_code')
import pandas as pd
from urosig.geo import read_series_matrix
R = r'D:\article\UroSig-Reeval\04_cohorts\raw\GSE54460'
_, ph = read_series_matrix(R + r'\GSE54460_series_matrix.txt.gz')
with gzip.open(R + r'\GSE54460_FPKM-genes-TopHat2-106samples-12-4-13.txt.gz', 'rt') as f:
    hdr = [next(f).rstrip('\n').split('\t') for _ in range(14)]
meta = pd.DataFrame({r[0].strip(): r[2:] for r in hdr[1:14] if r[0].strip()}, index=hdr[0][2:])
num = lambda s: pd.to_numeric(s.astype(str).str.replace('"', '').str.replace(',', '').str.strip(), errors='coerce')
key = ph['external id'].str.strip()
for col, e in meta['External ID'].str.strip().items():
    hit = key.index[(key == e) | key.str.endswith(e)]
    if len(hit) > 1:
        print(col, e, 'fpkm reads', num(pd.Series([meta.loc[col, 'no. of uniquely mapped reads']]))[0], 'fu', meta.loc[col,'Months total F/U'],
              '| geo:', [(h, ph.loc[h,'title'], ph.loc[h,'no of uniquely mapped reads'], ph.loc[h,'months total f/u']) for h in hit])
