"""Recover the 17-lncRNA multivariate Cox coefficients of PMID 35414328 from the authors' own supplementary data.

Usage: python refit_35414328.py <supp_dir>   (supp_dir = 03_extraction/fulltext_manual/supp_35414328)

The paper gives no coefficients. Supplementary file 5 (lassoSigExp.xlsx) holds the training-set expression of the
17 signature lncRNAs (Ensembl IDs) with Time/Status (319 patients, 58 events); Supplementary file 4 gives the
authors' univariate Cox coefficient for each candidate lncRNA by name.
Step 1: ENSG -> paper name, by refitting univariate Cox per ENSG and matching to file 4's coefficients.
Step 2: refit the multivariate Cox on all 17 (the model the paper describes).
Check : the median of the refitted risk score must equal the threshold printed in the paper's Fig 3c (-2.5775).
"""
import sys
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

warnings.filterwarnings("ignore")
D = sys.argv[1].rstrip("/") + "/"
d = pd.read_excel(D + "Supplementaryfile5.lassoSigExp.xlsx", header=0)
u = pd.read_excel(D + "Supplementaryfile4.ARg_lncRNAsunvarCox.results.xlsx", header=0)

name = {}
for e in d.columns[2:]:
    c = CoxPHFitter().fit(d[["Time", "Status", e]], "Time", "Status").params_[e]
    j = (u["coef"] - c).abs().idxmin()
    name[e] = u.loc[j, "symbol"]
    print(f"{e}\t{name[e]}\tunivariate refit {c:.6f} vs paper {u.loc[j, 'coef']:.6f}")

m = CoxPHFitter().fit(d, "Time", "Status")
lp = d[m.params_.index].values @ m.params_.values
print("n", len(d), "events", int(d["Status"].sum()), "| median risk score", round(float(np.median(lp)), 4),
      "(paper Fig 3c threshold: -2.5775)")
print("genes:", ";".join(name[e] for e in m.params_.index))
print("coefs:", ";".join(f"{v:.6f}" for v in m.params_.values))
print("ensg :", ";".join(m.params_.index))
