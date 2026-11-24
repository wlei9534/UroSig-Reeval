import sys; sys.path.insert(0, r'D:\article\UroSig-Reeval\05_code')
from urosig.geo import read_series_matrix
R = r'D:\article\UroSig-Reeval\04_cohorts\raw'
for g, cols in {'GSE54460': ['tissue','bcr','pstage','gleason score','tgleason','sms','prepsa','collection site'],
                'GSE70769': ['tumour gleason','biochemical relapse (bcr)','pathology stage','positive surgical margins (psm)','extra capsular extension (ece)','psa at diag','clinical stage'],
                'GSE70768': ['sample type','tumour gleason','biochemical relapse (bcr)','pathology stage','positive surgical margins (psm)','psa at diag'],
                'GSE116918': ['t-stage','gleason grade','psa (ng/ml)','tissue','bcr event (1=yes, 0=no)']}.items():
    _, p = read_series_matrix(fr'{R}\{g}\{g}_series_matrix.txt.gz')
    if g == 'GSE70768': p = p[p['sample type'] == 'Tumour']
    print('=====', g, len(p))
    for c in cols:
        vc = p[c].value_counts(dropna=False)
        print(f'  {c} [{len(vc)}]:', dict(vc.head(14)) if len(vc) <= 30 or c in ('pstage','pathology stage') else f'numeric-ish e.g. {list(vc.index[:5])}')
