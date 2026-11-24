"""QC of processed cohorts: label sanity (known risk factors vs BCR), clinical completeness, gene coverage."""
import json, sys
from pathlib import Path
ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else r'D:\article\UroSig-Reeval'); sys.path.insert(0, str(ROOT / '05_code'))
import numpy as np, pandas as pd
from urosig.metrics import harrell_c, hr_per_sd
P = ROOT / '04_cohorts' / 'processed'; OUT = ROOT / '06_results' / 'qc'; OUT.mkdir(parents=True, exist_ok=True)
PROLIF = ['MKI67', 'TOP2A', 'CCNB1', 'BIRC5', 'CDK1', 'BUB1', 'CENPF', 'KIF20A', 'RRM2', 'TPX2', 'ASPM', 'PTTG1']  # generic cell-cycle genes, QC only
rows, genesets = [], {}
for c in ['GSE54460', 'GSE70769', 'GSE70768', 'MSKCC_GSE21032', 'GSE116918', 'TCGA-PRAD']:
    e = pd.read_csv(P / c / 'expr_log2.csv.gz', index_col=0); k = pd.read_csv(P / c / 'clinical.csv', index_col=0)
    genesets[c] = set(e.index)
    d = k.dropna(subset=['time_bcr', 'bcr'])
    z = e.loc[[g for g in PROLIF if g in e.index], d.index]; z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1), axis=0)
    prolif = z.mean(axis=0)
    r = {'cohort': c, 'n': len(d), 'events': int(d.bcr.sum()), 'prolif_genes': z.shape[0]}
    for name, x in [('prolif', prolif), ('gleason_sum', d['gleason_sum']), ('log_psa', d['log_psa'])]:
        ok = x.notna()
        if ok.sum() > 20 and d.loc[ok, 'bcr'].sum() >= 5:
            h = hr_per_sd(d.loc[ok, 'time_bcr'].values, d.loc[ok, 'bcr'].values, x[ok].values)
            r[f'{name}_C'] = round(harrell_c(d.loc[ok, 'time_bcr'], d.loc[ok, 'bcr'], x[ok]), 3)
            r[f'{name}_HR'] = f"{h['hr']:.2f} ({h['lo']:.2f}-{h['hi']:.2f})"
    for col in ['psa', 'gleason_sum', 'isup', 'pt_cat', 'margin']:
        r[f'miss_{col}'] = int(d[col].isna().sum()) if col in d else 'n/a'
    rows.append(r)
q = pd.DataFrame(rows); q.to_csv(OUT / 'qc_label_sanity.csv', index=False)
print(q.to_string(index=False))
common = set.intersection(*genesets.values())
print('genes common to all cohorts:', len(common), '| per cohort:', {c: len(s) for c, s in genesets.items()})
