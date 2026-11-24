"""Original modelling endpoint of each locked signature, for sensitivity analysis S6 (no outcome data read).

Usage: python classify_endpoint.py <project_root>
Source: the 'endpoint' field the blind reviewers recorded from each paper (03_extraction/review/out_*.jsonl).
Rule (fixed before any evaluation was run): the endpoint term that appears FIRST in the recorded text is taken as
the endpoint the model was built on.
  BCR/RFS class : BCR, BCRFS, BRFS/bRFS, BFS, BFFS, biochemical, PSA-free, RFS, DFS, recurrence/relapse, 复发
  other         : OS, overall, PFS, PFI, 进展, MFS, metastasis, survival without qualifier, 未明确
32264916_OS is OS by construction.
Writes 06_results/eligibility/original_endpoint_class.csv.
"""
import glob
import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(sys.argv[1])
X = ROOT / "03_extraction"
BCR = r"BCR|BRFS|bRFS|\bBFS\b|BFFS|biochemical|PSA-free|\bRFS\b|\bDFS\b|recurrence|relapse|复发"
OTHER = r"\bOS\b|overall|\bPFS\b|\bPFI\b|进展|\bMFS\b|metasta|未明确"

rev = {}
for f in sorted(glob.glob(str(X / "review" / "out_*.jsonl"))):
    for line in open(f, encoding="utf-8"):
        r = json.loads(line)
        if "pmid" in r and "endpoint" in r:
            rev.setdefault(r["pmid"], r)

sig = pd.read_csv(X / "locked" / "signatures_locked_2026-11-01.csv", dtype=str).fillna("")
sig = sig[sig["分析角色"] != "excluded"]
rows = []
for _, s in sig.iterrows():
    txt = "OS (32264916_OS model)" if s["sig_id"].endswith("_OS") else rev.get(s["PMID"], {}).get("endpoint", "")
    mb, mo = re.search(BCR, txt, re.I), re.search(OTHER, txt, re.I)
    if mb and (not mo or mb.start() < mo.start()):
        cls = "BCR/RFS"
    elif mo:
        cls = "other"
    else:
        cls = "unclear"
    rows.append({"sig_id": s["sig_id"], "PMID": s["PMID"], "分析角色": s["分析角色"], "原文终点（复核员记录）": txt,
                 "终点类别": cls})
out = pd.DataFrame(rows)
out.to_csv(ROOT / "06_results" / "eligibility" / "original_endpoint_class.csv", index=False, encoding="utf-8-sig")
print(out.groupby(["分析角色", "终点类别"]).size().to_string())
