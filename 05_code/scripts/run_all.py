"""Run the primary analysis and sensitivity analyses S1-S10 of SAP v2 on this machine, in parallel shards.

Usage: python run_all.py <project_root> [--jobs 14] [--only primary,S3,...]
Each evaluation run is split into --jobs shards (seeds are per signature x cohort, so results do not depend on
the split), merged, then pooled. Progress: 06_results/evaluation/run_all.log; per-shard output: <tag>/logs/.
"""
import argparse
import os
import subprocess
import sys
import time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from urosig.paths import p as P, im_globs, _V as PATHS  # noqa: E402


ap = argparse.ArgumentParser()
ap.add_argument("root")
ap.add_argument("--jobs", type=int, default=14)
ap.add_argument("--only", default="")
ap.add_argument("--pool-only", action="store_true", help="re-run only the pooling steps on existing runs")
ap.add_argument("--nboot", default="1000")
ap.add_argument("--nrand", default="1000")
ap.add_argument("--force", action="store_true", help="re-run evaluations even if merged results or complete shards exist")
a = ap.parse_args()
ROOT = Path(a.root).resolve()
S = ROOT / "05_code" / "scripts"
EV = P(ROOT, "ev")
EV.mkdir(parents=True, exist_ok=True)
LOG = EV / "run_all.log"
PY = sys.executable
MAIN4 = "GSE54460,GSE70769,GSE70768,MSKCC_GSE21032"
ENV = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONIOENCODING="utf-8")


def log(msg):
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    print(line, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def n_expected(cohorts, roles):
    """Pairs the primary-type runs must produce: eligible (cov 0.8) signature x cohort in the active lock."""
    import pandas as pd
    e = pd.read_csv(P(ROOT, "elig"), dtype={"sig_id": str})
    return int((e["可评估_0.8"] & e["cohort"].isin(cohorts.split(",")) & e["分析角色"].isin(roles.split(","))).sum())


def evaluate(tag, extra, expect=0):
    if a.pool_only:
        return
    d = EV / tag
    (d / "logs").mkdir(parents=True, exist_ok=True)
    if (d / "pairs.csv").exists() and not a.force:
        log(f"[{tag}] already merged, skipped (use --force to redo)")
        return
    t0 = time.time()
    if all((d / f"pairs_shard{i:02d}of{a.jobs:02d}.csv").exists() and (d / f"run_info_shard{i:02d}of{a.jobs:02d}.json").exists()
           for i in range(a.jobs)) and not a.force:
        log(f"[{tag}] all {a.jobs} shards already present, merging only")
        procs = None
    else:
        procs = []
    for i in range(a.jobs if procs is not None else 0):
        out = open(d / "logs" / f"shard{i:02d}.log", "w", encoding="utf-8")
        cmd = [PY, str(S / "run_evaluation.py"), str(ROOT), "--tag", tag, "--shard", f"{i}/{a.jobs}"] + extra
        procs.append((subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT, env=ENV), out))
    if procs:
        log(f"[{tag}] started {a.jobs} shards: {' '.join(extra)}")
    bad = []
    for i, (p, out) in enumerate(procs or []):
        p.wait()
        out.close()
        if p.returncode:
            bad.append(i)
    if bad:
        log(f"[{tag}] FAILED shards {bad}; see {d / 'logs'}")
        raise SystemExit(1)
    r = subprocess.run([PY, str(S / "merge_shards.py"), str(ROOT), "--tag", tag] + (["--expect", str(expect)] if expect else []),
                       capture_output=True, text=True, env=ENV)
    log(f"[{tag}] done in {time.time() - t0:.0f}s; {r.stdout.strip()} {r.stderr.strip()[-300:]}")
    if r.returncode:
        raise SystemExit(1)


def pool(tag, src, *filters):
    r = subprocess.run([PY, str(S / "pool_results.py"), str(ROOT), "--tag", tag, "--src", src] + list(filters),
                       capture_output=True, text=True, env=ENV)
    log(f"[pool {tag}] rc={r.returncode} {r.stdout.strip()[-1500:]} {r.stderr.strip()[-500:]}")


def exploratory():
    r = subprocess.run([PY, str(S / "exploratory.py"), str(ROOT), "--jobs", str(a.jobs), "--nboot", a.nboot],
                       capture_output=True, text=True, env=ENV)
    log(f"[E1E2] rc={r.returncode} {r.stdout.strip()[-1500:]} {r.stderr.strip()[-800:]}")


BASE = ["--cov", "0.8", "--nboot", a.nboot, "--nrand", a.nrand, "--seed", "20261101"]
PLAN = [
    ("primary", lambda: (evaluate("run_primary", BASE + ["--mode", "z", "--cohorts", MAIN4, "--roles", "main"], n_expected(MAIN4, "main")),
                         pool("primary", "run_primary"))),
    ("S1", lambda: pool("S1_cov100", "run_primary", "--cov", "1.0")),
    ("S2", lambda: (evaluate("run_raw", BASE + ["--mode", "raw", "--cohorts", MAIN4, "--roles", "main"], n_expected(MAIN4, "main")),
                    pool("S2a_raw_log2", "run_raw"),
                    evaluate("run_tcga", BASE + ["--mode", "tcga", "--cohorts", MAIN4, "--roles", "main"]),
                    pool("S2b_tcga_sd", "run_tcga"))),
    ("S3", lambda: (evaluate("run_116918", BASE + ["--mode", "z", "--cohorts", "GSE116918", "--roles", "main"]),
                    pool("S3_plus_GSE116918", "run_primary,run_116918"))),
    ("S4", lambda: pool("S4_without_MSKCC", "run_primary", "--cohorts", "GSE54460,GSE70769,GSE70768")),
    ("S5", lambda: pool("S5_GSE54460_only", "run_primary", "--cohorts", "GSE54460")),
    ("S6", lambda: pool("S6_endpoint_BCR_RFS", "run_primary", "--endpoint", "BCR/RFS")),
    ("S7", lambda: (evaluate("run_sens_roles", BASE + ["--mode", "z", "--cohorts", MAIN4, "--roles", "sens_figure,sens_hr,sens_lasso"]),
                    pool("S7_plus_low_precision", "run_primary,run_sens_roles"))),
    ("S8", lambda: (evaluate("run_equal", BASE + ["--mode", "z", "--cohorts", MAIN4, "--roles", "equal_weight,equal_signed"]),
                    pool("S8_equal_weight", "run_equal"))),
    ("S9", lambda: evaluate("run_full_formula", BASE + ["--mode", "z", "--cohorts", MAIN4, "--roles", "main",
                                                        "--full-formula", "--only", "34002865,39832058"])),
    ("S10", lambda: (evaluate("run_drop_repeats", BASE + ["--mode", "z", "--cohorts", "GSE54460", "--roles", "main", "--drop-repeats"]),
                     pool("S10_GSE54460_no_repeats", "run_primary,run_drop_repeats"))),
    ("E3", lambda: (evaluate("run_met_116918", BASE + ["--mode", "z", "--cohorts", "GSE116918", "--roles", "main", "--endpoint", "met"]),
                    pool("E3_met_GSE116918", "run_met_116918"))),
    ("E12", lambda: None if a.pool_only else exploratory()),
]
only = set(a.only.split(",")) if a.only else None
log(f"run_all start: jobs={a.jobs} python={sys.version.split()[0]} only={a.only or 'all'}")
for name, fn in PLAN:
    if (only and name not in only) or (not only and name.startswith("E")):
        continue
    fn()
log("run_all finished")
