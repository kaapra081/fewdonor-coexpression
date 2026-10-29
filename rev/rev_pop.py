"""R1 (and R4 with --real): python rev_pop.py --sub P2 --k0 0 --k1 100"""
import os, sys, json, argparse, time, numpy as np, scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "dev"))
from revm import stats_from_cells, methods  # noqa: E402
from methods_v2 import bh  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--sub", required=True); ap.add_argument("--k0", type=int, default=0); ap.add_argument("--k1", type=int, default=100)
ap.add_argument("--real", action="store_true"); a = ap.parse_args(); si = ["P1", "P2", "P3"].index(a.sub); P = os.path.join(HERE, "..", "pop")
X = sp.load_npz(os.path.join(P, f"pop_{a.sub}_X.npz")).tocsr(); M = np.load(os.path.join(P, f"pop_{a.sub}_meta.npz"), allow_pickle=True)
lib = M["lib"]; don = M["donor"].astype(str); sym = M["sym"].astype(str); iu = np.triu_indices(X.shape[1], 1)
def run(lab, seed):
    ud = np.unique(lab); rows = [np.where(lab == u)[0] for u in ud]
    f, h = stats_from_cells(((X[r].toarray().astype(np.float64), lib[r].astype(np.float64)) for r in rows), seed)
    return methods(f, h, iu)
if a.real:
    o = run(don, 44_000_000)
    Y = {"UTY", "USP9Y", "ZFY", "DDX3Y", "KDM5D", "NLGN4Y", "TTTY14"}; xi = np.where(sym == "XIST")[0]
    res = []
    for k in range(len(iu[0])):
        g1, g2 = sym[iu[0][k]], sym[iu[1][k]]
        if ("XIST" in (g1, g2)) and (g1 in Y or g2 in Y):
            res.append(dict(pair=f"{g1}-{g2}", p_XDS=float(o["XDS"][k]), sign_XDS=int(o["sign_XDS"][k]), p_PB=float(o["PB"][k])))
        elif g1 in Y and g2 in Y:
            res.append(dict(pair=f"{g1}-{g2}", p_XDS=float(o["XDS"][k]), sign_XDS=int(o["sign_XDS"][k]), p_PB=float(o["PB"][k])))
    json.dump(dict(sub=a.sub, calls={m: int(bh(o[m]).sum()) for m in ("PB", "XDS", "XDSM", "KX")}, sexpairs=res), open(os.path.join(HERE, "out", f"real_{a.sub}.json"), "w"), indent=1)
    print(a.sub, "real", {m: int(bh(o[m]).sum()) for m in ("PB", "XDS", "XDSM", "KX")}, len(res), flush=True); sys.exit()
for k in range(a.k0, a.k1):
    fn = os.path.join(HERE, "out", f"R1_{a.sub}_k{k}.json")
    if os.path.exists(fn): continue
    t0 = time.time(); seed = 48_000_000 + 10 * k + si; lab = np.random.default_rng(seed).permutation(don); o = run(lab, seed)
    r = {m: dict(V=int(bh(o[m]).sum()), minp=float(o[m].min()), n3=int((o[m] < 1e-3).sum()), n4=int((o[m] < 1e-4).sum()), n5=int((o[m] < 1e-5).sum()),
                 t05=float((o[m] < 0.05).mean())) for m in ("PB", "XDS", "XDSM", "KX")}
    r["P"] = int(len(iu[0])); json.dump(r, open(fn, "w"))
    print(a.sub, k, {m: r[m]["V"] for m in ("PB", "XDS", "XDSM", "KX")}, round(time.time() - t0), "s", flush=True)
