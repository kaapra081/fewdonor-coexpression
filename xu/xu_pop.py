"""E1: python xu_pop.py --sub P2 --k 3"""
import os, sys, argparse, time, numpy as np, scipy.sparse as sp, torch
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("dstd", "pop"): sys.path.insert(0, os.path.join(HERE, "..", p))
sys.path.insert(0, HERE)
from pop_core import run_methods, subset, bh  # noqa: E402
from dstd import donor_level_cells, methods_d  # noqa: E402
from xu import xu  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--sub", required=True); ap.add_argument("--k", type=int, required=True); ap.add_argument("--threads", type=int, default=1)
a = ap.parse_args(); torch.set_num_threads(a.threads); si = ["P1", "P2", "P3"].index(a.sub)
P = os.path.join(HERE, "..", "pop")
X = sp.load_npz(os.path.join(P, f"pop_{a.sub}_X.npz")).tocsr(); M = np.load(os.path.join(P, f"pop_{a.sub}_meta.npz"), allow_pickle=True)
lib = M["lib"]; don = M["donor"].astype(str); G = X.shape[1]; iu = np.triu_indices(G, 1)
od = os.path.join(HERE, "out"); os.makedirs(od, exist_ok=True); t0 = time.time()
seed = 41_000_000 + 10 * a.k + si; rng = np.random.default_rng(seed); lab = rng.permutation(don)
ud = np.unique(lab); rows_of = [np.where(lab == u)[0] for u in ud]
cell = lambda i: (X[rows_of[i]].toarray().astype(np.float64), lib[rows_of[i]].astype(np.float64))
d = donor_level_cells((cell(i) for i in range(len(ud))), seed + 7); d["donors"] = ud; Dall = len(ud)
order = rng.permutation(Dall)
for D in (Dall, 400, 200, 100):
    if D > Dall: continue
    idx = np.sort(order[:D]); sd = subset(d, idx)
    o = {}
    if D >= 200:
        r1 = run_methods(sd, iu, seed=seed + D); r2 = methods_d(sd, iu, seed=seed + D + 3)
        o.update(PB=r1["PB"], LCP=r1["LCP"], XDS=r2["LCPXd"], r_w=r1["r_w"])
    else:
        import pop_core; from methods_v2 import pearson_fisher
        o["PB"] = pearson_fisher(np.log1p(1e6 * sd["m"]), iu)["p"]; o["r_w"] = run_methods.__globals__["np"].zeros(0)
    o["XU"] = xu((cell(i) for i in idx), sd["m"], sd["psi"], sd["kappa"], sd["S"], iu)["p"]
    if D < 200:   # r_w for top-1% from the full-D computation of the same pseudo-donor set is not available; recompute cheaply
        r1 = run_methods.__globals__  # placeholder to keep namespace
        Cw = sd["C"][:, iu[0], iu[1]].astype(np.float64).mean(0)
        sig2 = np.maximum((sd["psi"] - sd["m"] / sd["S"][:, None]) / sd["kappa"][:, None], 0).mean(0)
        with np.errstate(invalid="ignore", divide="ignore"): o["r_w"] = Cw / np.sqrt(sig2[iu[0]] * sig2[iu[1]])
    tag = "all" if D == Dall else D
    np.savez_compressed(os.path.join(od, f"XU_{a.sub}_k{a.k}_D{tag}.npz"), **{m: np.asarray(v, np.float32) for m, v in o.items()})
    print(a.sub, a.k, tag, {m: int(bh(o[m].astype(float)).sum()) for m in ("PB", "LCP", "XDS", "XU") if m in o}, round(time.time() - t0), "s", flush=True)
