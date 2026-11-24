import json, sys, subprocess; sys.path.insert(0, r'D:\article\UroSig-Reeval\05_code')
import pandas as pd
from urosig.geo import read_series_matrix
R = r'D:\article\UroSig-Reeval\04_cohorts\raw'
P = pd.DataFrame(json.load(open(R + r'\MSKCC\cbioportal_patient_clinical.json'))).pivot(index='patientId', columns='clinicalAttributeId', values='value')
S = pd.DataFrame(json.load(open(R + r'\MSKCC\cbioportal_sample_clinical.json')))
s2p = S.drop_duplicates('sampleId').set_index('sampleId')['patientId']
_, ph = read_series_matrix(R + r'\MSKCC\GSE21032-GPL10264_series_matrix.txt.gz')
prim = ph[ph['tumor type'] == 'Primary tumor']
pid = prim['sample id'].map(lambda x: s2p.get(x, x))
d = P.reindex(pid.values)
ok = d['DFS_MONTHS'].notna() & d['DFS_STATUS'].notna()
print('MSKCC primaries with expr', len(prim), '| linked to cBio patient', int(pid.isin(P.index).sum()), '| with DFS', int(ok.sum()), '| events', int(d.loc[ok, 'DFS_STATUS'].str.startswith('1').sum()))
for st in ['prad_tcga_pan_can_atlas_2018', 'prad_tcga']:
    out = subprocess.run(['curl.exe', '-s', f'https://www.cbioportal.org/api/studies/{st}/clinical-attributes?projection=SUMMARY'], capture_output=True, text=True).stdout
    try:
        a = json.loads(out); print(st, [x['clinicalAttributeId'] for x in a if any(k in x['clinicalAttributeId'] for k in ['DFS', 'PFS', 'BCR', 'BIOCHEM', 'PSA', 'GLEASON', 'MARGIN', 'PATH_T', 'DSS'])])
    except Exception as e: print(st, 'ERR', out[:200])
