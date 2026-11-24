"""Input/output locations per extraction lock. UROSIG_LOCK=v1 reproduces the 2026-11-01 results; default is v2.

v1: 292 keyword-prioritised papers, locked 2026-11-01 (results in 06_results/evaluation).
v2: v1 + 300 papers from the dual title/abstract screening, locked 2026-11-14 (results in 06_results/evaluation_v2).
"""
import os
from pathlib import Path

LOCK = os.environ.get("UROSIG_LOCK", "v2")
_V = {
    "v1": dict(sig="03_extraction/locked/signatures_locked_2026-11-01.csv",
               terms="03_extraction/locked/signature_terms_locked_2026-11-01.csv",
               dirs="06_results/eligibility/equal_weight_directions.csv",
               leak="06_results/eligibility/leakage_matrix.csv",
               cov="06_results/coverage/gene_coverage_by_signature.csv",
               missing="06_results/coverage/missing_terms.csv",
               elig="06_results/eligibility/eligibility_signature_by_cohort.csv",
               elig_summary="06_results/eligibility/eligibility_summary.csv",
               endpoint="06_results/eligibility/original_endpoint_class.csv",
               im=["03_extraction/review/out_IM*.jsonl"],
               ev="06_results/evaluation"),
    "v2": dict(sig="03_extraction/locked/signatures_locked_v2_2026-11-14.csv",
               terms="03_extraction/locked/signature_terms_locked_v2_2026-11-14.csv",
               dirs="06_results/eligibility/equal_weight_directions_v2.csv",
               leak="06_results/eligibility/leakage_matrix_v2.csv",
               cov="06_results/coverage/gene_coverage_by_signature_v2.csv",
               missing="06_results/coverage/missing_terms_v2.csv",
               elig="06_results/eligibility/eligibility_signature_by_cohort_v2.csv",
               elig_summary="06_results/eligibility/eligibility_summary_v2.csv",
               endpoint="06_results/eligibility/original_endpoint_class_v2.csv",
               im=["03_extraction/review/out_IM*.jsonl", "03_extraction/v2/im/out_IM*.jsonl"],
               ev="06_results/evaluation_v2"),
}[LOCK]


def p(root, key):
    return Path(root) / _V[key]


def im_globs(root):
    return [str(Path(root) / g) for g in _V["im"]]
