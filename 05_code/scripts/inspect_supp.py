import json, gzip, sys; sys.path.insert(0, r'D:\article\UroSig-Reeval\05_code')
import pandas as pd
from urosig.geo import read_series_matrix
R = r'D:\article\UroSig-Reeval\04_cohorts\raw'
P = pd.DataFrame(json.load(open(R + r'\MSKCC\cbioportal_patient_clinical.json'))).pivot(index='patientId', columns='clinicalAttributeId', values='value')
S = pd.DataFrame(json.load(open(R + r'\MSKCC\cbioportal_sample_clinical.json'))).pivot(index='sampleId', columns='clinicalAttributeId', values='value')
print('cBio patients', len(P), 'samples', len(S), '| sample types', S['SAMPLE_TYPE'].value_counts().to_dict() if 'SAMPLE_TYPE' in S else '')
print('DFS_STATUS', P['DFS_STATUS'].value_counts().to_dict(), '| DFS_MONTHS non-null', int(P['DFS_MONTHS'].notna().sum()))
print('ids', list(S.index[:3]), list(P.index[:3]))
e, ph = read_series_matrix(R + r'\MSKCC\GSE21032-GPL10264_series_matrix.txt.gz')
print('GEO gene-level', e.shape, list(ph['title'][:5]))
print('pheno cols', list(ph.columns))
for c in ph.columns:
    if c not in ('title',): print('  ', c, dict(ph[c].value_counts().head(6)))
c = pd.read_csv(R + r'\TCGA-PRAD\TCGA-PRAD.clinical.tsv.gz', sep='\t', low_memory=False)
print('TCGA clin', c.shape)
print([x for x in c.columns if any(k in x.lower() for k in ['biochem', 'recur', 'psa', 'gleason', 'margin', 'pathologic_t', 'follow', 'days_to', 'new_tumor', 'primary_pattern', 'secondary_pattern'])])
sv = pd.read_csv(R + r'\TCGA-PRAD\TCGA-PRAD.survival.tsv.gz', sep='\t'); print('TCGA surv', sv.shape, list(sv.columns))
