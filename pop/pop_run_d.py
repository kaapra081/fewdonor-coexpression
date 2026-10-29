"""POP_PREREG Part A (PB, LCP, LCPX, LCPXS) + DSTD_PREREG E3 (LCPd, LCPXd, PBlin) on identical pseudo-donors.
python pop_run_d.py --sub P2 --k 0 --Ds 1164,800,400,200   (D values processed in the order given; 'all' = every pseudo-donor)"""
import os, sys, argparse, time, numpy as np, scipy.sparse as sp, torch
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "dstd"))
from pop_core import run_methods, subset, bh  # noqa: E402
from dstd import donor_level_cells, methods_d  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--sub", required=True); ap.add_argument("--k", type=int, default=0); ap.add_argument("--Ds", required=True)
ap.add_argument("--threads", type=int, default=1); ap.add_argument("--halves", action="store_true"); a = ap.parse_args(); torch.set_num_threads(a.threads)
si = ["P1", "P2", "P3"].index(a.sub)
X = sp.load_npz(os.path.join(HERE, f"pop_{a.sub}_X.npz")).tocsr(); M = np.load(os.path.join(HERE, f"pop_{a.sub}_meta.npz"), allow_pickle=True)
lib = M["lib"]; don = M["donor"].astype(str); G = X.shape[1]; iu = np.triu_indices(G, 1)
od = os.path.join(HERE, "out"); os.makedirs(od, exist_ok=True); t0 = time.time()
seed = 41_000_000 + 10 * a.k + si; rng = np.random.default_rng(seed); lab = rng.permutation(don)
ud = np.unique(lab); cells = ((X[np.where(lab == u)[0]].toarray().astype(np.float64), lib[np.where(lab == u)[0]].astype(np.float64)) for u in ud)
d = donor_level_cells(cells, seed + 7); d["donors"] = ud; Dall = len(ud)
order = rng.permutation(Dall)
def save(name, o):
    np.savez_compressed(os.path.join(od, name), **{k: np.asarray(v, np.float32) for k, v in o.items()})
    print(name, {m: int(bh(o[m].astype(float)).sum()) for m in ("PB", "LCP", "LCPX", "LCPXS", "LCPd", "LCPXd", "PBlin") if m in o}, round(time.time() - t0), "s", flush=True)
for Ds in a.Ds.split(","):
    D = Dall if Ds == "all" else int(Ds)
    if D > Dall: continue
    sd = subset(d, np.sort(order[:D])); o = run_methods(sd, iu, seed=seed + D); o.update(methods_d(sd, iu, seed=seed + D + 3))
    save(f"AD_{a.sub}_k{a.k}_D{'all' if D == Dall else D}.npz", o)
if a.halves:
    for j, rows in enumerate((np.sort(order[: Dall // 2]), np.sort(order[Dall // 2:]))):
        sd = subset(d, rows); o = run_methods(sd, iu, seed=seed + 99 + j); o.update(methods_d(sd, iu, seed=seed + 199 + j)); save(f"ADhalf_{a.sub}_k{a.k}_h{j}.npz", o)
