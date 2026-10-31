"""REV3_PREREG A3 (small-D run-level, clip [0.1,10], seeds 60M) and B3 (clip [0.5,2], P3, seeds 62M). python rev3_A3B3.py --mode A3 --sub P1 --k0 0 --k1 10"""
import os, sys, json, argparse, time, numpy as np, scipy.sparse as sp, torch
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("dev", "dstd", "pop", "rev"): sys.path.insert(0, os.path.join(HERE, "..", p))
from methods_v2 import bh, paule_mandel  # noqa: E402
from dstd import _ms, donor_level_cells, methods_d  # noqa: E402
from pop_core import run_methods, subset  # noqa: E402
from revm import stats_from_cells, methods as rmethods  # noqa: E402
torch.set_num_threads(1)
ap = argparse.ArgumentParser(); ap.add_argument("--mode", required=True); ap.add_argument("--sub", required=True); ap.add_argument("--k0", type=int, default=0); ap.add_argument("--k1", type=int, default=100); a = ap.parse_args()
si = ["P1", "P2", "P3"].index(a.sub); P = os.path.join(HERE, "..", "pop"); OD = os.path.join(HERE, "out"); os.makedirs(OD, exist_ok=True)
lo, hi = (0.1, 10.0) if a.mode == "A3" else (0.5, 2.0); base = 60_000_000 if a.mode == "A3" else 62_000_000
DS = (25, 50, 100) if a.mode == "A3" else (200, 400, None)
X = sp.load_npz(os.path.join(P, f"pop_{a.sub}_X.npz")).tocsr(); M = np.load(os.path.join(P, f"pop_{a.sub}_meta.npz"), allow_pickle=True)
lib = M["lib"]; don = M["donor"].astype(str); iu = np.triu_indices(X.shape[1], 1)
ud = np.unique(don); rows = [np.where(don == u)[0] for u in ud]
st = [_ms(X[r].toarray().astype(np.float64), lib[r].astype(np.float64)) for r in rows]
m = np.array([s_[0] for s_ in st]); ps = np.array([s_[1] for s_ in st]); tau2, mu, _ = paule_mandel(m, ps)
th = mu[None, :] + tau2[None, :] / (tau2[None, :] + ps) * (m - mu[None, :]); f = np.clip(mu[None, :] / np.maximum(th, 1e-30), lo, hi)
idx = np.empty(X.shape[0], dtype=np.int64)
for i, r in enumerate(rows): idx[r] = i
Xs = sp.csr_matrix(X.multiply(f[idx]))
for k in range(a.k0, a.k1):
    fn = os.path.join(OD, f"{a.mode}_{a.sub}_k{k}.json")
    if os.path.exists(fn): continue
    t0 = time.time(); seed = base + 10 * k + si; rng = np.random.default_rng(seed); lab = rng.permutation(don)
    pu = np.unique(lab); prow = [np.where(lab == u)[0] for u in pu]; cells = [(Xs[r].toarray(), lib[r].astype(np.float64)) for r in prow]
    d = donor_level_cells(cells, seed + 7); d["donors"] = pu; fs, hs = stats_from_cells(cells, seed); Dall = len(pu)
    order = rng.permutation(Dall); res = {}
    for D in DS:
        D = D or Dall; sel = np.sort(order[:D]); sd = subset(d, sel)
        o1 = run_methods(sd, iu, seed=seed + D); o2 = methods_d(sd, iu, seed=seed + D + 3)
        fsub = {kk: v[sel] for kk, v in fs.items()}; hsub = [{kk: v[sel] for kk, v in h.items()} for h in hs]
        o3 = rmethods(fsub, hsub, iu)
        pv = dict(PB=o1["PB"], LCP=o1["LCP"], LCPd=o2["LCPd"], XDSperm=o2["LCPXd"], XDS=o3["XDS"], XDSM=o3["XDSM"], KX=o3["KX"])
        res["all" if D == Dall else str(D)] = {mm: dict(V=int(bh(np.asarray(v, float)).sum()), t05=float((np.asarray(v) < 0.05).mean()), minp=float(np.min(v))) for mm, v in pv.items()}
    json.dump(res, open(fn, "w")); print(a.mode, a.sub, k, {D: {mm: r[mm]["V"] for mm in r} for D, r in res.items()}, round(time.time() - t0), "s", flush=True)
