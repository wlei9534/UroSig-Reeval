"""Merge the shard outputs of run_evaluation.py into pairs.csv / random_c.npz / run_info.json.

Usage: python merge_shards.py <project_root> --tag primary [--expect 384]
"""
import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from urosig.paths import p as P, im_globs, _V as PATHS  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("root")
ap.add_argument("--tag", required=True)
ap.add_argument("--expect", type=int, default=0)
a = ap.parse_args()
D = P(a.root, "ev") / a.tag
pf = sorted(glob.glob(str(D / "pairs_shard*.csv")))
infos = [json.loads(Path(f).read_text(encoding="utf-8")) for f in sorted(glob.glob(str(D / "run_info_shard*.json")))]
n = {int(i["args"]["shard"].split("/")[1]) for i in infos}
assert len(n) == 1 and len(pf) == n.pop() == len(infos), f"incomplete shards in {D}: {len(pf)} files"
def _read(f):
    try:
        return pd.read_csv(f, dtype={"sig_id": str})
    except pd.errors.EmptyDataError:  # shard with no pairs
        return pd.DataFrame()


p = pd.concat([_read(f) for f in pf], ignore_index=True)
assert not p.duplicated(["sig_id", "cohort"]).any()
if a.expect:
    assert len(p) == a.expect, f"{len(p)} pairs, expected {a.expect}"
p = p.sort_values(["cohort", "sig_id"]).reset_index(drop=True)
p.to_csv(D / "pairs.csv", index=False, encoding="utf-8-sig")
store = {}
for f in sorted(glob.glob(str(D / "random_c_shard*.npz"))):
    z = np.load(f)
    for k in z.files:
        assert k not in store
        store[k] = z[k]
np.savez_compressed(D / "random_c.npz", **store)
info = dict(infos[0])
info["args"] = {k: v for k, v in info["args"].items() if k != "shard"}
info.update(n_pairs=len(p), n_shards=len(pf), seconds_max_shard=max(i["seconds"] for i in infos),
            seconds_sum=sum(i["seconds"] for i in infos))
(D / "run_info.json").write_text(json.dumps(info, ensure_ascii=False, indent=1), encoding="utf-8")
print(f"{a.tag}: {len(p)} pairs ({(p['status'] == 'ok').sum()} ok), {len(store)} random sets, {len(pf)} shards")
