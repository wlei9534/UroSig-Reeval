"""Readers for GEO series matrix files (expression + sample annotations)."""
from __future__ import annotations

import csv
import gzip
import io
from pathlib import Path

import numpy as np
import pandas as pd


def _open(path):
    path = Path(path)
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else open(path, encoding="utf-8", errors="replace")


def _cells(line):
    return next(csv.reader([line.rstrip("\n")], delimiter="\t"))


def read_series_matrix(path):
    """Return (expr: probes x samples float DataFrame, pheno: samples x fields DataFrame).

    Characteristics 'key: value' are split into columns named by key; title and
    source fields are kept as 'title' and 'source'.
    """
    meta, table, in_table = {}, [], False
    accessions = None
    char_rows = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                table.append(line)
                continue
            if not line.startswith("!Sample_"):
                continue
            c = _cells(line)
            key, vals = c[0][len("!Sample_"):], c[1:]
            if key == "geo_accession":
                accessions = vals
            elif key.startswith("characteristics"):
                char_rows.append(vals)
            elif key in ("title", "source_name_ch1"):
                meta["title" if key == "title" else "source"] = vals
    pheno = pd.DataFrame(index=accessions)
    for k, v in meta.items():
        pheno[k] = v
    for row in char_rows:
        for acc, cell in zip(accessions, row):
            if ":" in cell:
                k, v = cell.split(":", 1)
                pheno.loc[acc, k.strip()] = v.strip()
    if table:
        expr = pd.read_csv(io.StringIO("".join(table)), sep="\t", index_col=0, na_values=["null", "NA", ""])
        expr.index = expr.index.astype(str)
        expr = expr.apply(pd.to_numeric, errors="coerce").astype(float)
    else:
        expr = pd.DataFrame(columns=accessions, dtype=float)
    return expr, pheno


def collapse_to_genes(expr: pd.DataFrame, probe2gene: pd.Series) -> pd.DataFrame:
    """Map probes to HGNC symbols; keep the probe with the highest mean expression per gene."""
    g = probe2gene.reindex(expr.index)
    keep = g.notna() & (g.astype(str).str.len() > 0)
    e = expr[keep].copy()
    e["_gene"] = g[keep].values
    e["_mean"] = e.drop(columns="_gene").mean(axis=1)
    e = e.sort_values("_mean", ascending=False).drop_duplicates("_gene")
    return e.set_index("_gene").drop(columns="_mean")
