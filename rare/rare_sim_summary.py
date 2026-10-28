import os, sys, glob, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "dev"))
from methods_v2 import bh  # noqa: E402
out = {}
for od in sorted(glob.glob(os.path.join(HERE, "sim_S*"))):
    tr = np.load(os.path.join(od, "truth.npz")); alt = tr["altB"]; cls = tr["cls"]
    res = {m: dict(fdp=[], pw=[], ss=[]) for m in ("LCP", "PB")}
    for fn in sorted(glob.glob(os.path.join(od, "rep*.npz"))):
        z = np.load(fn)
        for m in res:
            p = z[m].astype(float); r = bh(p); n = r.sum()
            res[m]["fdp"].append((r & ~alt).sum() / max(n, 1)); res[m]["pw"].append((r & alt).sum() / alt.sum()); res[m]["ss"].append((p[cls == "SS"] < 0.05).mean())
    R = len(res["LCP"]["fdp"])
    out[os.path.basename(od)] = {m: {k: [float(np.mean(v)), float(np.std(v, ddof=1) / np.sqrt(R))] for k, v in d.items()} | {"R": R} for m, d in res.items()}
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, "rare_sim_summary.json"), "w"), indent=1)
