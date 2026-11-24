"""Fetch PMC open-access article files from the PMC OA S3 bucket (pmc-oa-opendata) for the v2 full-text set.
Usage: python -I fetch_pmc_s3_v2.py <ids.csv with PMID,PMCID> <out_dir>
For each PMCID the highest available version's JSON manifest is read; all listed files are downloaded to
<out_dir>/fulltext/<PMCID>/ and the JATS XML is also copied to <out_dir>/fulltext_xml/<PMCID>.xml."""
import csv, json, sys, time, urllib.request, urllib.error, shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
B = "https://pmc-oa-opendata.s3.amazonaws.com/"
ids = list(csv.DictReader(open(sys.argv[1], encoding="utf-8")))
OUT = Path(sys.argv[2])
def get(u, t=120):
    for k in range(3):
        try:
            return urllib.request.urlopen(u, timeout=t).read()
        except urllib.error.HTTPError as e:
            if e.code in (403, 404):
                return None
            err = e
        except Exception as e:
            err = e
        time.sleep(3)
    raise err
def one(r):
    pmc = r["PMCID"]; res = {"PMID": r["PMID"], "PMCID": pmc}
    man, ver = None, None
    for v in range(1, 6):
        j = get(f"{B}{pmc}.{v}/{pmc}.{v}.json")
        if j is None:
            break
        man, ver = json.loads(j), v
    if man is None:
        res["status"] = "not in OA bucket"; return res
    d = OUT / "fulltext" / pmc; d.mkdir(parents=True, exist_ok=True)
    keys = []
    for k, v in man.items():
        for u in (v if isinstance(v, list) else [v]):
            if isinstance(u, str) and u.startswith("s3://pmc-oa-opendata/"):
                keys.append(u.split("s3://pmc-oa-opendata/")[1].split("?")[0])
    n = 0
    for k in sorted(set(keys)):
        p = d / Path(k).name
        if not p.exists():
            b = get(B + k, 300)
            if b is None:
                continue
            p.write_bytes(b)
        n += 1
    xmls = [p for p in d.glob("*.xml")]
    if xmls:
        shutil.copy(xmls[0], OUT / "fulltext_xml" / f"{pmc}.xml")
    res.update(status="ok", version=ver, files=n, xml=bool(xmls), pdf=bool(list(d.glob(f"{pmc}*.pdf"))),
               license=man.get("license_code", ""))
    return res
with ThreadPoolExecutor(8) as ex:
    out = list(ex.map(one, ids))
json.dump(out, open(OUT / "fetch_s3_log.json", "w"), indent=1)
from collections import Counter
print(Counter(r["status"] for r in out), sum(r.get("xml", False) for r in out), "with xml")
