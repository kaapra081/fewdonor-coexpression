"""POP_PREREG Part B + DSTD methods on real donors. python pop_real.py --sub P2"""
import os, sys, argparse, time, numpy as np, scipy.sparse as sp, torch
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "dstd"))
from pop_core import run_methods, subset, bh  # noqa: E402
from dstd import donor_level_cells, methods_d  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--sub", required=True); ap.add_argument("--threads", type=int, default=2); a = ap.parse_args(); torch.set_num_threads(a.threads)
si = ["P1", "P2", "P3"].index(a.sub)
X = sp.load_npz(os.path.join(HERE, f"pop_{a.sub}_X.npz")).tocsr(); M = np.load(os.path.join(HERE, f"pop_{a.sub}_meta.npz"), allow_pickle=True)
lib = M["lib"]; don = M["donor"].astype(str); G = X.shape[1]; iu = np.triu_indices(G, 1); od = os.path.join(HERE, "out"); t0 = time.time()
seed = 42_000_000 + si; ud = np.unique(don)
cells = ((X[np.where(don == u)[0]].toarray().astype(np.float64), lib[np.where(don == u)[0]].astype(np.float64)) for u in ud)
d = donor_level_cells(cells, seed + 7); d["donors"] = ud; Dall = len(ud)
def save(name, o):
    np.savez_compressed(os.path.join(od, name), **{k: np.asarray(v, np.float32) for k, v in o.items()})
    print(name, {m: int(bh(o[m].astype(float)).sum()) for m in ("PB", "LCP", "LCPX", "LCPXS", "LCPd", "LCPXd", "PBlin")}, round(time.time() - t0), "s", flush=True)
o = run_methods(d, iu, seed=seed); o.update(methods_d(d, iu, seed=seed + 3)); save(f"B_{a.sub}_full.npz", o)
order = np.random.default_rng(42_000_000).permutation(Dall)
for j, rows in enumerate((np.sort(order[: Dall // 2]), np.sort(order[Dall // 2:]))):
    sd = subset(d, rows); o = run_methods(sd, iu, seed=seed + 1 + j); o.update(methods_d(sd, iu, seed=seed + 11 + j)); save(f"B_{a.sub}_h{j}.npz", o)
