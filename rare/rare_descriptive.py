"""Descriptive (not in RARE_PREREG.md decision rule): replication of S1/S2 call sets in Lee 2026 (independent donors, full depth) and
within-donor correlation in the tested cells themselves."""
import os, sys, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "real"))
from methods_v2 import bh  # noqa: E402
z = np.load(os.path.join(HERE, "..", "real", "runs", "lee_full.npz")); G = 417; a, b = np.triu_indices(G, 1)
lee = pd.DataFrame(dict(g1=a, g2=b, lee_LCP=bh(z["LCP"].astype(float)), lee_T=z["LCP_stat"].astype(float))).set_index(["g1", "g2"])
out = {}
for key in ("S1", "S2"):
    P = pd.read_csv(os.path.join(HERE, f"rare_pairs_{key}.csv.gz")).set_index(["g1", "g2"]).join(lee)
    sets = {"PB-only": P[P.PB & ~P.LCP], "Both": P[P.PB & P.LCP], "LCP-only": P[~P.PB & P.LCP], "neither": P[~P.PB & ~P.LCP]}
    r = {}
    for k, X in sets.items():
        sg = np.sign(X.r_PB) if k != "LCP-only" else np.sign(X.T_LCP)
        r[k] = dict(n=int(len(X)), repl_lee=float((X.lee_LCP & (np.sign(X.lee_T) == sg)).mean()) if len(X) else None,
                    med_absrw_own=float(np.nanmedian(np.abs(X.r_w))) if len(X) else None, med_leakshare=float(np.nanmedian(X.share)) if len(X) else None)
    out[key] = r; print(key, json.dumps(r))
json.dump(out, open(os.path.join(HERE, "rare_descriptive.json"), "w"), indent=1)
