"""Exploratory analyses E1 and E2 (pre-specified in decision_log 2026-11-12; not for conclusions).

Usage: python exploratory.py <project_root> [--jobs 14] [--nboot 1000] [--seed 20261101]
E1  paired bootstrap of (signature C - grade C) on the same resampled patients, per signature x cohort of the primary
    analysis; pooled per signature by random effects (REML + Hartung-Knapp) on the difference scale.
E2  proliferation benchmark meta-PCNA (Venet et al. 2011, PMID 22028643): the 1% of genes most positively correlated
    (Pearson) with PCNA in TCGA-PRAD (reference cohort, outcome not used), PCNA itself excluded, candidates restricted
    to genes measured in all four primary cohorts. Score = mean within-cohort z of those genes. Reports its C per cohort,
    paired-bootstrap (signature C - meta-PCNA C), and Spearman correlation of each signature score with meta-PCNA.
Outputs: 06_results/evaluation/exploratory/{metapcna_genes.csv, metapcna_cohorts.csv, pairs_E1E2.csv, signatures_E1E2.csv, summary_E1E2.json}
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from urosig.paths import p as P, im_globs, _V as PATHS  # noqa: E402

MAIN4 = ["GSE54460", "GSE70769", "GSE70768", "MSKCC_GSE21032"]


def read_expr(root, c):
    with gzip.open(Path(root) / "04_cohorts/processed" / c / "expr_log2.csv.gz", "rt", encoding="utf-8") as f:
        return pd.read_csv(f, index_col=0)


def load(root, c):
    expr = read_expr(root, c)
    clin = pd.read_csv(Path(root) / "04_cohorts/processed" / c / "clinical.csv", index_col=0)
    clin = clin[clin["time_bcr"].notna() & clin["bcr"].notna()]
    clin = clin.loc[clin.index.intersection(expr.columns)]
    return expr[clin.index], clin


def zmean(expr, genes):
    x = expr.loc[[g for g in genes if g in expr.index]].astype(float)
    z = x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1, ddof=1).replace(0, np.nan), axis=0).fillna(0)
    return z.mean(axis=0)


def harrell(t, e, s):
    from lifelines.utils import concordance_index
    return float(concordance_index(t, -np.asarray(s, float), e))


def work(job):
    root, c, sids, nboot, seed = job
    sys.path.insert(0, str(Path(root) / "05_code"))
    from urosig.sigscore import score_signature
    from scipy.stats import spearmanr
    terms = pd.read_csv(P(root, "terms"), dtype=str).fillna("")
    dirs_df = pd.read_csv(P(root, "dirs"), dtype={"sig_id": str})
    mp_genes = pd.read_csv(P(root, "ev") / "exploratory/metapcna_genes.csv")["gene"].tolist()
    expr, clin = load(root, c)
    mp = zmean(expr, mp_genes).reindex(clin.index).to_numpy(float)
    T = clin["time_bcr"].to_numpy(float); E = clin["bcr"].to_numpy(int)
    if clin["isup"].notna().mean() >= 0.9:
        G = clin["isup"].astype(float).to_numpy()
    else:
        g = clin["gleason_sum"].astype(float)
        G = pd.Series(np.select([g <= 6, g == 7], [1.0, 2.0], 3.0), index=clin.index).where(g.notna()).to_numpy(float)
    out = []
    for sid in sids:
        tt = terms[terms["sig_id"] == sid]
        d = dict(zip(dirs_df.loc[dirs_df["sig_id"] == sid, "基因（原文写法）"], dirs_df.loc[dirs_df["sig_id"] == sid, "方向"]))
        sc = score_signature(sid, c, tt, expr, "z", d, None)
        s = sc.values.reindex(clin.index).to_numpy(float)
        ok = ~np.isnan(s) & ~np.isnan(G)
        t, e, x, gr, m = T[ok], E[ok], s[ok], G[ok], mp[ok]
        rng = np.random.default_rng(zlib.crc32(f"{seed}|E1E2|{sid}|{c}".encode()))
        c_s, c_g, c_m = harrell(t, e, x), harrell(t, e, gr), harrell(t, e, m)
        dg, dm = [], []
        n = len(t)
        for _ in range(nboot):
            i = rng.integers(0, n, n)
            if e[i].sum() < 2:
                continue
            cs = harrell(t[i], e[i], x[i])
            dg.append(cs - harrell(t[i], e[i], gr[i]))
            dm.append(cs - harrell(t[i], e[i], m[i]))
        dg, dm = np.array(dg), np.array(dm)
        rho = spearmanr(s[~np.isnan(s)], mp[~np.isnan(s)]).statistic
        out.append({"sig_id": sid, "cohort": c, "n": n, "events": int(e.sum()), "c_sig": c_s, "c_grade": c_g, "c_metapcna": c_m,
                    "d_grade": c_s - c_g, "d_grade_se": dg.std(ddof=1), "d_grade_lo": np.percentile(dg, 2.5), "d_grade_hi": np.percentile(dg, 97.5),
                    "d_mp": c_s - c_m, "d_mp_se": dm.std(ddof=1), "d_mp_lo": np.percentile(dm, 2.5), "d_mp_hi": np.percentile(dm, 97.5),
                    "rho_sig_metapcna": float(rho), "nboot_used": len(dg)})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--jobs", type=int, default=14)
    ap.add_argument("--nboot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261101)
    a = ap.parse_args()
    root = Path(a.root).resolve()
    sys.path.insert(0, str(root / "05_code"))
    from urosig.meta import random_effects
    OUT = P(root, "ev") / "exploratory"
    OUT.mkdir(parents=True, exist_ok=True)

    # meta-PCNA genes from TCGA-PRAD (no outcome used)
    tc = read_expr(root, "TCGA-PRAD").astype(float)
    common = set(tc.index)
    for c in MAIN4:
        common &= set(read_expr(root, c).index)
    tc = tc.loc[sorted(common)]
    tc = tc[tc.std(axis=1) > 0]
    r = tc.T.corrwith(tc.loc["PCNA"]).drop("PCNA").sort_values(ascending=False)
    k = int(round(0.01 * len(r)))
    pd.DataFrame({"gene": r.index[:k], "pearson_r_with_PCNA_TCGA": r.values[:k]}).to_csv(OUT / "metapcna_genes.csv", index=False)
    print(f"meta-PCNA: {k} genes (top 1% of {len(r)} candidate genes; r from {r.iloc[k-1]:.3f} to {r.iloc[0]:.3f})", flush=True)

    pairs = pd.read_csv(P(root, "ev") / "run_primary/pairs.csv", dtype={"sig_id": str})
    pairs = pairs[pairs["status"] == "ok"]
    jobs = []
    for c in MAIN4:
        sids = pairs.loc[pairs["cohort"] == c, "sig_id"].tolist()
        for i in range(0, len(sids), 8):
            jobs.append((str(root), c, sids[i:i + 8], a.nboot, a.seed))
    rows = []
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for k_, res in enumerate(ex.map(work, jobs)):
            rows += res
            print(f"{k_ + 1}/{len(jobs)} chunks", flush=True)
    p = pd.DataFrame(rows).sort_values(["cohort", "sig_id"])
    p.to_csv(OUT / "pairs_E1E2.csv", index=False, encoding="utf-8-sig")

    # meta-PCNA C per cohort (same patients as the grade benchmark; identical for all pairs of a cohort when no NaN)
    coh = p.groupby("cohort").agg(n=("n", "max"), events=("events", "max"), c_metapcna=("c_metapcna", "median"),
                                  c_grade=("c_grade", "median"), c_sig_median=("c_sig", "median")).reset_index()
    coh.to_csv(OUT / "metapcna_cohorts.csv", index=False, encoding="utf-8-sig")

    srows = []
    for sid, g in p.groupby("sig_id"):
        r_ = {"sig_id": sid, "k": len(g)}
        for key in ("d_grade", "d_mp"):
            y, se = g[key].to_numpy(float), g[key + "_se"].to_numpy(float)
            if len(g) > 1:
                re = random_effects(y, se ** 2)
                est, lo, hi = re["est"], re["lo"], re["hi"]
            else:
                est, lo, hi = y[0], g[key + "_lo"].iloc[0], g[key + "_hi"].iloc[0]
            r_.update({key: est, key + "_lo": lo, key + "_hi": hi})
        r_["median_rho_metapcna"] = float(g["rho_sig_metapcna"].median())
        srows.append(r_)
    s = pd.DataFrame(srows)
    s.to_csv(OUT / "signatures_E1E2.csv", index=False, encoding="utf-8-sig")
    from scipy.stats import spearmanr
    rho_c = spearmanr(p["rho_sig_metapcna"].abs(), p["c_sig"]).statistic
    summ = {"n_signatures": len(s), "n_pairs": len(p), "metapcna_n_genes": k,
            "E1_grade": {"median_pooled_dC": float(s["d_grade"].median()), "n_ci_above_0": int((s["d_grade_lo"] > 0).sum()),
                         "n_ci_below_0": int((s["d_grade_hi"] < 0).sum()), "n_ci_includes_0": int(((s["d_grade_lo"] <= 0) & (s["d_grade_hi"] >= 0)).sum())},
            "E2_metapcna": {"cohort_C": dict(zip(coh["cohort"], coh["c_metapcna"].round(4))),
                            "median_pooled_dC_sig_minus_metapcna": float(s["d_mp"].median()),
                            "n_ci_above_0": int((s["d_mp_lo"] > 0).sum()), "n_ci_below_0": int((s["d_mp_hi"] < 0).sum()),
                            "n_point_above_0": int((s["d_mp"] > 0).sum()),
                            "median_abs_rho_sig_metapcna": float(p["rho_sig_metapcna"].abs().median()),
                            "spearman_absrho_vs_C_pairs": float(rho_c)}}
    (OUT / "summary_E1E2.json").write_text(json.dumps(summ, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps(summ, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
