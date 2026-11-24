"""Fetch open-access full text (JATS XML) and supplementary files from Europe PMC for the v2 full-text set.
Usage: python -I fetch_epmc_v2.py <ids.csv with PMID,PMCID> <out_dir>"""
import csv, sys, time, urllib.request, zipfile, io, json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
ids = list(csv.DictReader(open(sys.argv[1], encoding="utf-8")))
OUT = Path(sys.argv[2])
B = "https://www.ebi.ac.uk/europepmc/webservices/rest/"
def get(u, t=240):
    for k in range(3):
        try:
            return urllib.request.urlopen(u, timeout=t).read()
        except Exception as e:
            err = e; time.sleep(5)
    raise err
def one(r):
    pmc = r["PMCID"]; res = {"PMID": r["PMID"], "PMCID": pmc}
    x = OUT / "fulltext_xml" / f"{pmc}.xml"
    try:
        if not x.exists():
            b = get(B + f"{pmc}/fullTextXML"); x.write_bytes(b)
        res["xml"] = x.stat().st_size
    except Exception as e:
        res["xml"] = f"ERR {e}"
    d = OUT / "fulltext" / pmc
    try:
        if not d.exists():
            b = get(B + f"{pmc}/supplementaryFiles", 600)
            d.mkdir(parents=True)
            with zipfile.ZipFile(io.BytesIO(b)) as z:
                for n in z.namelist():
                    p = (d / Path(n).name)
                    if n.endswith("/"):
                        continue
                    p.write_bytes(z.read(n))
        res["supp_files"] = len(list(d.iterdir()))
    except Exception as e:
        res["supp_files"] = f"ERR {e}"
    return res
with ThreadPoolExecutor(4) as ex:
    out = list(ex.map(one, ids))
json.dump(out, open(OUT / "fetch_log.json", "w"), indent=1)
print(sum(isinstance(r["xml"], int) for r in out), "xml ok of", len(out))
