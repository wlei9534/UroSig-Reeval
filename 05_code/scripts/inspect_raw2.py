import gzip, sys
sys.path.insert(0, r'D:\article\UroSig-Reeval\05_code')
import pandas as pd
from urosig.geo import read_series_matrix
R = r'D:\article\UroSig-Reeval\04_cohorts\raw'
with gzip.open(R + r'\GSE54460\GSE54460_FPKM-genes-TopHat2-106samples-12-4-13.txt.gz', 'rt') as f:
    rows = [next(f).split('\t') for _ in range(40)]
for i, r in enumerate(rows): print(i, '|', r[0][:40], '|', r[1][:25], '|', r[2][:25], '|', len(r))
_, ph = read_series_matrix(R + r'\GSE54460\GSE54460_series_matrix.txt.gz')
print(ph[['title', 'cgsr id', 'external id', 'bcr', 'months to bcr', 'months total f/u']].head(8).to_string())
print('dup external ids:', ph['external id'].duplicated().sum(), 'dup cgsr', ph['cgsr id'].duplicated().sum())
e, _ = read_series_matrix(R + r'\GSE116918\GSE116918_series_matrix.txt.gz')
ids = set()
with gzip.open(R + r'\platforms\GPL25318_family.soft.gz', 'rt', errors='replace') as f:
    on = False
    for l in f:
        if l.startswith('!platform_table_begin'): on = True; next(f); continue
        if l.startswith('!platform_table_end'): break
        if on: ids.add(l.split('\t', 1)[0])
print('GPL25318 table ids', len(ids), 'matrix ids', len(e.index), 'overlap', len(ids & set(e.index)), list(ids)[:3])
