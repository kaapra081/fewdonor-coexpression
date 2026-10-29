"""python pop_run.py --sub P2 --part A --k 0   |   --part B"""
import os, sys, argparse, time, json, numpy as np, scipy.sparse as sp, torch
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from pop_core import donor_level, subset, run_methods, bh  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--sub", required=True); ap.add_argument("--part", required=True); ap.add_argument("--k", type=int, default=0)
ap.add_argument("--threads", type=int, default=2); a = ap.parse_args(); torch.set_num_threads(a.threads)
si = ["P1", "P2", "P3"].index(a.sub)
X = sp.load_npz(os.path.join(HERE, f"pop_{a.sub}_X.npz")).tocsr(); M = np.load(os.path.join(HERE, f"pop_{a.sub}_meta.npz"), allow_pickle=True)
lib = M["lib"]; don = M["donor"].astype(str); G = X.shape[1]; iu = np.triu_indices(G, 1)
od = os.path.join(HERE, "out"); os.makedirs(od, exist_ok=True)
def save(name, o):
    np.savez_compressed(os.path.join(od, name), **{k: np.asarray(v, np.float32) for k, v in o.items()})
    print(name, {m: int(bh(o[m].astype(float)).sum()) for m in ("PB", "LCP", "LCPX", "LCPXS")}, round(time.time() - t0), "s", flush=True)
t0 = time.time()
if a.part == "A":
    seed = 41_000_000 + 10 * a.k + si; rng = np.random.default_rng(seed)
    lab = rng.permutation(don)                                  # pseudo-donors: n_d kept, nuclei random
    d = donor_level(X, lib, lab, split_seed=seed + 7); Dall = len(d["donors"])
    order = rng.permutation(Dall)                                # nested subsets
    for D in (25, 50, 100, 200, 400, 800, Dall):
        if D > Dall: continue
        o = run_methods(subset(d, np.sort(order[:D])), iu, seed=seed + D); save(f"A_{a.sub}_k{a.k}_D{D}.npz", o)
    h = np.sort(order[: Dall // 2]), np.sort(order[Dall // 2:])
    for j, rows in enumerate(h):
        o = run_methods(subset(d, rows), iu, seed=seed + 99 + j); save(f"Ahalf_{a.sub}_k{a.k}_h{j}.npz", o)
else:
    seed = 42_000_000 + si; d = donor_level(X, lib, don, split_seed=seed + 7); Dall = len(d["donors"])
    o = run_methods(d, iu, seed=seed); save(f"B_{a.sub}_full.npz", o)
    order = np.random.default_rng(42_000_000).permutation(Dall)
    for j, rows in enumerate((np.sort(order[: Dall // 2]), np.sort(order[Dall // 2:]))):
        o = run_methods(subset(d, rows), iu, seed=seed + 1 + j); save(f"B_{a.sub}_h{j}.npz", o)
