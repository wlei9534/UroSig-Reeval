"""Build harmonised prostate-cancer cohorts: genes x samples (log2) + clinical table.

Usage: python build_prad_cohorts.py <project_root>
Outputs: 04_cohorts/processed/<cohort>/{expr_log2.csv.gz, clinical.csv, qc.json} and
         04_cohorts/processed/cohort_summary.csv
"""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else r"D:\article\UroSig-Reeval")
sys.path.insert(0, str(ROOT / "05_code"))
from urosig.clinical import bcr_time, isup_group, parse_gleason, parse_pt, to_float  # noqa: E402
from urosig.geo import collapse_to_genes, read_series_matrix  # noqa: E402
from urosig.genes import build_symbol_map, harmonise_index  # noqa: E402
from urosig.preprocess import fpkm_to_log2tpm, quantile_normalise  # noqa: E402

RAW, OUT = ROOT / "04_cohorts" / "raw", ROOT / "04_cohorts" / "processed"
HGNC = pd.read_csv(RAW / "reference" / "hgnc_complete_set_2026-10-26.txt", sep="\t", dtype=str,
                   usecols=["symbol", "status", "prev_symbol", "alias_symbol", "refseq_accession", "ensembl_gene_id"])
HGNC = HGNC[HGNC["status"] == "Approved"]
SYMMAP = build_symbol_map(HGNC)


def gleason_cols(raw: pd.Series) -> pd.DataFrame:
    pats = raw.map(parse_gleason)
    p = pats.map(lambda x: x[0] if x else np.nan)
    s = pats.map(lambda x: x[1] if x else np.nan)
    isup = [isup_group(int(a), int(b)) if pd.notna(a) else np.nan for a, b in zip(p, s)]
    return pd.DataFrame({"gleason_primary": p, "gleason_secondary": s, "gleason_sum": p + s, "isup": isup}, index=raw.index)


def finish(name, expr, clin, notes, setting):
    clin = clin.copy()
    clin["cohort"], clin["setting"] = name, setting
    clin["log_psa"] = np.log(clin["psa"].where(clin["psa"] > 0))
    clin["isup_cat"] = clin["isup"].map(lambda v: np.nan if pd.isna(v) else ("4-5" if v >= 4 else str(int(v))))
    expr = expr[clin.index]
    d = OUT / name
    d.mkdir(parents=True, exist_ok=True)
    expr.round(4).to_csv(d / "expr_log2.csv.gz", compression="gzip")
    clin.to_csv(d / "clinical.csv")
    usable = clin.dropna(subset=["time_bcr", "bcr"])
    v = expr.to_numpy().ravel()
    v = v[~np.isnan(v)]
    qc = {"cohort": name, "setting": setting, "samples": int(expr.shape[1]), "genes": int(expr.shape[0]),
          "bcr_usable": int(len(usable)), "bcr_events": int(usable["bcr"].sum()),
          "median_followup_censored": float(usable.loc[usable.bcr == 0, "time_bcr"].median()),
          "expr_quantiles_0_50_100": [round(float(x), 2) for x in np.percentile(v, [0, 50, 100])],
          "missing_clinical": {c: int(clin[c].isna().sum()) for c in ["psa", "isup", "gleason_sum", "pt_cat", "margin"] if c in clin},
          "notes": notes}
    (d / "qc.json").write_text(json.dumps(qc, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps({k: qc[k] for k in ["cohort", "samples", "genes", "bcr_usable", "bcr_events"]}, ensure_ascii=False))
    return qc


def illumina_map():
    with gzip.open(RAW / "platforms" / "GPL10558.annot.gz", "rt", errors="replace") as f:
        lines = f.readlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("ID\t"))
    end = next((i for i, l in enumerate(lines) if l.startswith("!platform_table_end")), len(lines))
    from io import StringIO
    a = pd.read_csv(StringIO("".join(lines[start:end])), sep="\t", dtype=str, usecols=["ID", "Gene symbol"])
    a = a[a["Gene symbol"].notna() & ~a["Gene symbol"].str.contains("///", na=False)]
    return a.set_index("ID")["Gene symbol"]


def almac_map():
    rows, on = [], False
    with gzip.open(RAW / "platforms" / "GPL25318_family.soft.gz", "rt", errors="replace") as f:
        for l in f:
            if l.startswith("!platform_table_begin"):
                hdr = next(f).rstrip("\n").split("\t")
                iid, isym = hdr.index("ID"), hdr.index("Gene Symbol")
                on = True
                continue
            if l.startswith("!platform_table_end"):
                break
            if on:
                c = l.rstrip("\n").split("\t")
                rows.append((c[iid], c[isym] if len(c) > isym else ""))
    a = pd.DataFrame(rows, columns=["ID", "sym"]).set_index("ID")["sym"]
    return a[(a != "---") & (a != "") & ~a.str.contains("///")]


def build_illumina(gse, tumour_only=False):
    expr, ph = read_series_matrix(RAW / gse / f"{gse}_series_matrix.txt.gz")
    notes = ["series-matrix values (submitter-normalised, log2); re-quantile-normalised within cohort"]
    if tumour_only:
        ph = ph[ph["sample type"] == "Tumour"]
        notes.append("benign samples removed")
    expr = quantile_normalise(expr[ph.index])
    g = collapse_to_genes(expr, illumina_map())
    g, rep = harmonise_index(g, SYMMAP)
    notes.append(f"probe->gene GPL10558 annot; HGNC harmonisation {rep}")
    T, E = bcr_time(ph["biochemical relapse (bcr)"], ph["time to bcr (months)"], ph["total follow up (months)"], ("Y",), ("N",))
    clin = pd.DataFrame({"patient_id": ph.index, "time_bcr": T, "bcr": E, "psa": to_float(ph["psa at diag"]),
                         "pt_cat": ph["pathology stage"].map(parse_pt),
                         "margin": ph["positive surgical margins (psm)"].str.upper().map({"Y": 1.0, "N": 0.0})}, index=ph.index)
    clin = clin.join(gleason_cols(ph["tumour gleason"]))
    return finish(gse, g, clin, notes, "radical prostatectomy")


def build_54460():
    _, ph = read_series_matrix(RAW / "GSE54460" / "GSE54460_series_matrix.txt.gz")
    with gzip.open(RAW / "GSE54460" / "GSE54460_FPKM-genes-TopHat2-106samples-12-4-13.txt.gz", "rt") as f:
        lines = f.readlines()
    hdr = [l.rstrip("\n").split("\t") for l in lines[:16]]
    meta = pd.DataFrame({r[0].strip(): r[2:] for r in hdr[1:14] if r[0].strip()}, index=hdr[0][2:])
    body = [l.rstrip("\n").split("\t") for l in lines[16:]]
    fp = pd.DataFrame([b[2:] for b in body], index=[b[1] for b in body], columns=hdr[0][2:]).apply(pd.to_numeric, errors="coerce")
    # link FPKM columns to GSM via External ID; GEO 'external id' may carry a site prefix
    ext = meta["External ID"].str.strip()
    key = ph["external id"].str.strip()
    num = lambda s: pd.to_numeric(s.astype(str).str.replace('"', "").str.replace(",", "").str.strip(), errors="coerce")
    reads_f, reads_g = num(meta["no. of uniquely mapped reads"]), num(ph["no of uniquely mapped reads"])
    gsm_of = {}
    for col, e in ext.items():
        hit = key.index[(key == e) | key.str.endswith(e)]
        if len(hit) > 1:  # repeat libraries of one patient: disambiguate by read count, then follow-up
            by_reads = [h for h in hit if reads_g[h] == reads_f[col]]
            hit = by_reads if len(by_reads) == 1 else [
                h for h in (by_reads or hit) if abs(num(pd.Series([ph.loc[h, "months total f/u"]]))[0] - num(pd.Series([meta.loc[col, "Months total F/U"]]))[0]) < 1]
        if len(hit) == 1:
            gsm_of[col] = hit[0]
    if len(gsm_of) != len(ext) or len(set(gsm_of.values())) != len(ext):
        raise RuntimeError(f"GSE54460 sample mapping incomplete: {len(gsm_of)}/{len(ext)}")
    # concordance check: follow-up and BCR carried in both files must agree
    m2 = meta.rename(index=gsm_of)
    fu_ok = ((num(m2["Months total F/U"]) - num(ph.loc[m2.index, "months total f/u"])).abs() < 1).mean()
    bcr_ok = (m2["BCR"].str.strip() == ph.loc[m2.index, "bcr"].str.strip()).mean()
    print(f"GSE54460 mapping: {len(gsm_of)} columns; follow-up agreement {fu_ok:.0%}, BCR agreement {bcr_ok:.0%}")
    if bcr_ok < 1 or fu_ok < 0.95:
        raise RuntimeError("GSE54460 FPKM/GEO clinical fields disagree after mapping")
    fp = fp.rename(columns=gsm_of)
    tpm = fpkm_to_log2tpm(fp.fillna(0))
    tpm = tpm.assign(_m=tpm.mean(axis=1)).sort_values("_m", ascending=False)
    tpm = tpm[~tpm.index.duplicated()].drop(columns="_m")
    g, rep = harmonise_index(tpm, SYMMAP)
    # one library per patient: keep the original library (GEO title without "-2");
    # the repeat entries carry different follow-up in GEO, so these patients are flagged for sensitivity analysis
    dup_pat = set(ph.loc[ph["external id"].duplicated(keep=False), "external id"])
    ph = ph[~ph["title"].str.endswith("-2")]
    if ph["external id"].duplicated().any():
        raise RuntimeError("GSE54460: duplicate patients remain after dedup")
    T, E = bcr_time(ph["bcr"], ph["months to bcr"], ph["months total f/u"], ("1",), ("0",))
    clin = pd.DataFrame({"patient_id": ph["external id"], "time_bcr": T, "bcr": E, "psa": to_float(ph["prepsa"]),
                         "pt_cat": ph["pstage"].map(parse_pt),
                         "margin": ph["sms"].str.lower().map({"positive": 1.0, "negative": 0.0}),
                         "site": ph["collection site"],
                         "flag_inconsistent_repeat": ph["external id"].isin(dup_pat).astype(int)}, index=ph.index)
    clin = clin.join(gleason_cols(ph["gleason score"]))
    notes = ["FPKM (TopHat2) -> TPM -> log2(TPM+1)", f"HGNC harmonisation {rep}",
             f"6 patients with 2 libraries: kept original library (title without -2); repeat entries have different follow-up in GEO -> flag_inconsistent_repeat=1 for sensitivity analysis; {len(ph)} patients",
             "time = months to BCR for events, months total f/u for censored",
             "development cohort of the 24-gene panel in PMID 24713434: exclude that panel here (leakage)"]
    return finish("GSE54460", g, clin, notes, "radical prostatectomy (FFPE, 3 sites)")


def build_116918():
    expr, ph = read_series_matrix(RAW / "GSE116918" / "GSE116918_series_matrix.txt.gz")
    expr = quantile_normalise(expr)
    g = collapse_to_genes(expr, almac_map())
    g, rep = harmonise_index(g, SYMMAP)
    T, E = bcr_time(ph["bcr event (1=yes, 0=no)"], ph["follow-up time (bcr, months)"], ph["follow-up time (bcr, months)"], ("1",), ("0",))
    Tm, Em = bcr_time(ph["met event (1=yes, 0=no)"], ph["follow-up time (met, months)"], ph["follow-up time (met, months)"], ("1",), ("0",))
    gs = to_float(ph["gleason grade"])
    clin = pd.DataFrame({"patient_id": ph.index, "time_bcr": T, "bcr": E, "time_met": Tm, "met": Em,
                         "psa": to_float(ph["psa (ng/ml)"]), "ct_stage": ph["t-stage"].where(ph["t-stage"] != "Unknown"),
                         "gleason_primary": np.nan, "gleason_secondary": np.nan, "gleason_sum": gs,
                         "isup": gs.map(lambda v: np.nan if pd.isna(v) else (1 if v <= 6 else (np.nan if v == 7 else (4 if v == 8 else 5)))),
                         "pt_cat": np.nan, "margin": np.nan}, index=ph.index)
    notes = ["Almac PCa DSA, series-matrix RMA values; re-quantile-normalised", f"probe->gene GPL25318 soft table; HGNC {rep}",
             "Gleason sum only: ISUP 2 vs 3 not separable (isup NaN for sum 7); use gleason_sum in clinical model",
             "radiotherapy biopsies: clinical T stage, no pT/margin"]
    return finish("GSE116918", g, clin, notes, "radical radiotherapy (diagnostic biopsy)")


def _cbio(path):
    d = pd.DataFrame(json.load(open(path, encoding="utf-8")))
    key = "sampleId" if "sampleId" in d and d["sampleId"].notna().any() and "SAMPLE" in path.name.upper() else "patientId"
    return d.pivot_table(index=key, columns="clinicalAttributeId", values="value", aggfunc="first"), d


def build_mskcc():
    expr, ph = read_series_matrix(RAW / "MSKCC" / "GSE21032-GPL10264_series_matrix.txt.gz")
    ph = ph[ph["tumor type"] == "Primary tumor"]
    rows, on = [], False
    with gzip.open(RAW / "platforms" / "GPL10264_family.soft.gz", "rt", errors="replace") as f:
        for l in f:
            if l.startswith("!platform_table_begin"):
                next(f); on = True; continue
            if l.startswith("!platform_table_end"):
                break
            if on:
                c = l.rstrip("\n").split("\t")
                rows.append((c[0], c[1].split(".")[0] if len(c) > 1 else ""))
    acc = pd.DataFrame(rows, columns=["ID", "acc"]).set_index("ID")["acc"]
    ref2sym = {}
    for sym, refs in zip(HGNC["symbol"], HGNC["refseq_accession"].fillna("")):
        for r in filter(None, refs.split("|")):
            ref2sym.setdefault(r.split(".")[0], sym)
    probe2gene = acc.map(ref2sym).dropna()
    expr = expr[ph.index]
    if np.nanmax(expr.values) > 100:  # GEO values on linear scale
        expr = np.log2(expr.clip(lower=1))
    expr = quantile_normalise(expr)
    g = collapse_to_genes(expr, probe2gene)
    g, rep = harmonise_index(g, SYMMAP)
    P, _ = _cbio(RAW / "MSKCC" / "cbioportal_patient_clinical.json")
    Sraw = pd.DataFrame(json.load(open(RAW / "MSKCC" / "cbioportal_sample_clinical.json", encoding="utf-8")))
    S = Sraw.pivot_table(index="sampleId", columns="clinicalAttributeId", values="value", aggfunc="first")
    s2p = Sraw.drop_duplicates("sampleId").set_index("sampleId")["patientId"]
    sid = ph["sample id"]
    pid = sid.map(lambda x: s2p.get(x, x))
    pc = P.reindex(pid.values).set_index(ph.index)
    sc = S.reindex(sid.values).set_index(ph.index)
    ev = pc["DFS_STATUS"].str.startswith("1").map({True: "1", False: "0"})
    T, E = bcr_time(ev, pc["DFS_MONTHS"], pc["DFS_MONTHS"], ("1",), ("0",))
    pats = (sc["GLEASON_SCORE_1"].astype(str) + "+" + sc["GLEASON_SCORE_2"].astype(str)) if {"GLEASON_SCORE_1", "GLEASON_SCORE_2"} <= set(sc.columns) else ph["biopsy_gleason_grade"]
    clin = pd.DataFrame({"patient_id": pid, "time_bcr": T, "bcr": E, "psa": np.nan,
                         "pt_cat": ph["pathological_stage"].map(parse_pt), "margin": np.nan}, index=ph.index)
    clin = clin.join(gleason_cols(pats))
    notes = [f"GEO GPL10264 gene-level exon array (series matrix, linear scale -> log2), re-quantile-normalised; probe->RefSeq->HGNC: {len(probe2gene)}/{len(acc)} probes mapped; {rep}",
             "outcome = cBioPortal prad_mskcc DFS_MONTHS/DFS_STATUS; equivalence with BCR definition of the original study 未核验",
             "Gleason patterns from cBioPortal sample attributes; PSA and margin not available -> clinical model without PSA"]
    return finish("MSKCC_GSE21032", g, clin, notes, "radical prostatectomy")


def build_tcga():
    x = pd.read_csv(RAW / "TCGA-PRAD" / "TCGA-PRAD.star_tpm.tsv.gz", sep="\t", index_col=0)
    x.index = x.index.str.split(".").str[0]
    ens2sym = HGNC.dropna(subset=["ensembl_gene_id"]).set_index("ensembl_gene_id")["symbol"]
    ens2sym = ens2sym[~ens2sym.index.duplicated(keep=False)]  # drop ambiguous Ensembl IDs
    x = x[~x.index.duplicated()]  # PAR_Y copies share the stripped Ensembl ID
    prim = [c for c in x.columns if c[13:15] == "01"]
    x = x[prim]
    x = x.loc[:, ~x.columns.str[:12].duplicated()]
    x.columns = x.columns.str[:12]
    g = x[x.index.isin(ens2sym.index)].copy()
    g.index = ens2sym.reindex(g.index).values
    g = g.assign(_m=g.mean(axis=1)).sort_values("_m", ascending=False)
    g = g[~g.index.duplicated()].drop(columns="_m")
    P, _ = _cbio(RAW / "TCGA-PRAD" / "cbioportal_prad_tcga_patient_clinical.json")
    gdc = pd.read_csv(RAW / "TCGA-PRAD" / "TCGA-PRAD.clinical.tsv.gz", sep="\t", low_memory=False)
    gdc["pid"] = gdc.iloc[:, 0].astype(str).str[:12]
    fu = pd.to_numeric(gdc.groupby("pid")["days_to_last_follow_up.diagnoses"].max(), errors="coerce")
    P = P.reindex(g.columns)
    ev = P["BIOCHEMICAL_RECURRENCE_INDICATOR"].map({"YES": "1", "NO": "0"})
    T, E = bcr_time(ev, to_float(P["DAYS_TO_BIOCHEMICAL_RECURRENCE_FIRST"]) / 30.44, fu.reindex(g.columns) / 30.44, ("1",), ("0",))
    pats = P["GLEASON_PATTERN_PRIMARY"].astype(str) + "+" + P["GLEASON_PATTERN_SECONDARY"].astype(str)
    clin = pd.DataFrame({"patient_id": g.columns, "time_bcr": T.values, "bcr": E.values, "psa": to_float(P["PSA_MOST_RECENT_RESULTS"]).values,
                         "pt_cat": P["PATH_T_STAGE"].map(parse_pt).values,
                         "margin": P["RESIDUAL_TUMOR"].map({"R0": 0.0, "R1": 1.0, "R2": 1.0}).values}, index=g.columns)
    clin = clin.join(gleason_cols(pats.set_axis(g.columns)))
    notes = ["GDC hub STAR TPM (as distributed, log2 scale), primary tumours (-01), one aliquot per patient; Ensembl->HGNC",
             "BCR: cBioPortal prad_tcga indicator + days to first BCR; censored time = GDC days_to_last_follow_up",
             "PSA field = most recent PSA (not preoperative) -> do not use as baseline covariate without check",
             "IN-SAMPLE REFERENCE ONLY: training cohort of most published signatures; never enters the main analysis"]
    return finish("TCGA-PRAD", g, clin, notes, "radical prostatectomy (in-sample reference)")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    only = sys.argv[2].split(",") if len(sys.argv) > 2 else None
    builders = {"GSE54460": build_54460, "GSE70769": lambda: build_illumina("GSE70769"),
                "GSE70768": lambda: build_illumina("GSE70768", tumour_only=True), "GSE116918": build_116918,
                "MSKCC_GSE21032": build_mskcc, "TCGA-PRAD": build_tcga}
    for k, f in builders.items():
        if only is None or k in only:
            f()
    qcs = [json.loads((OUT / k / "qc.json").read_text(encoding="utf-8")) for k in builders if (OUT / k / "qc.json").exists()]
    pd.DataFrame([{k: v for k, v in q.items() if k not in ("notes", "missing_clinical")} for q in qcs]).to_csv(OUT / "cohort_summary.csv", index=False)
