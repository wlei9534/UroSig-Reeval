import gzip, sys, itertools, re
sys.path.insert(0, r'D:\article\UroSig-Reeval\05_code')
import numpy as np, pandas as pd
from urosig.geo import read_series_matrix
R = r'D:\article\UroSig-Reeval\04_cohorts\raw'
with gzip.open(R + r'\GSE54460\GSE54460_FPKM-genes-TopHat2-106samples-12-4-13.txt.gz', 'rt') as f:
    head = [next(f) for _ in range(4)]
for l in head: print('FPKM|', l[:300].rstrip())
_, ph = read_series_matrix(R + r'\GSE54460\GSE54460_series_matrix.txt.gz')
print('54460 pheno cols', list(ph.columns)[:20]); print(ph[['title']].head(3).to_string())
for g in ['GSE70769', 'GSE70768', 'GSE116918']:
    e, p = read_series_matrix(fr'{R}\{g}\{g}_series_matrix.txt.gz')
    v = e.values.ravel(); v = v[~np.isnan(v)]
    print(g, e.shape, 'ids', list(e.index[:3]), 'q', np.round(np.percentile(v, [0, 1, 50, 99, 100]), 2), 'nan', int(np.isnan(e.values).sum()))
for name in ['GPL10558.annot.gz']:
    with gzip.open(R + r'\platforms\\' + name, 'rt', errors='replace') as f:
        for l in f:
            if l.startswith('ID\t'): print(name, l[:400]); print(next(f)[:300]); break
with gzip.open(R + r'\platforms\GPL25318_family.soft.gz', 'rt', errors='replace') as f:
    for l in f:
        if l.startswith('!platform_table_begin'):
            print('GPL25318 hdr', next(f)[:500]); print(next(f)[:400]); print(next(f)[:400]); break
