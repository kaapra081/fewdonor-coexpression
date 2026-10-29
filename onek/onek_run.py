"""ONEK_PREREG: python onek_run.py --typ B1 --k 0   (k = -1: real donors)"""
import os, sys, argparse, time, numpy as np, scipy.sparse as sp, torch
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("dstd", "pop"): sys.path.insert(0, os.path.join(HERE, "..", p))
from pop_core import run_methods, subset, bh  # noqa: E402
from dstd import donor_level_cells, methods_d  # noqa: E402
sys.path.insert(0, os.path.join(HERE, "..", "rev")); from revm import stats_from_cells, methods as rmethods  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--typ", required=True); ap.add_argument("--k", type=int, required=True); a = ap.parse_args(); torch.set_num_threads(1)
ti = ["B1", "B2", "B3"].index(a.typ)
X = sp.load_npz(os.path.join(HERE, f"onek_{a.typ}_X.npz")).tocsr(); M = np.load(os.path.join(HERE, f"onek_{a.typ}_meta.npz"), allow_pickle=True)
lib = M["lib"]; don = M["donor"].astype(str); iu = np.triu_indices(X.shape[1], 1); t0 = time.time()
seed = 45_000_000 + 10 * a.k + ti if a.k >= 0 else 47_000_000 + ti
rng = np.random.default_rng(seed); lab = rng.permutation(don) if a.k >= 0 else don
ud = np.unique(lab); rows = [np.where(lab == u)[0] for u in ud]
d = donor_level_cells(((X[r].toarray().astype(np.float64), lib[r].astype(np.float64)) for r in rows), seed + 7); d["donors"] = ud; Dall = len(ud)
order = rng.permutation(Dall)
for D in ((Dall, 200) if a.k >= 0 else (Dall,)):
    sd = subset(d, np.sort(order[:D])); o = run_methods(sd, iu, seed=seed + D); o.update(methods_d(sd, iu, seed=seed + D + 3))
    sel = np.sort(order[:D]); fs, hs = stats_from_cells(((X[rows[i]].toarray().astype(np.float64), lib[rows[i]].astype(np.float64)) for i in sel), seed + D + 11)
    o["XDSM"] = rmethods(fs, hs, iu)["XDSM"]
    keep = {m: o[m] for m in ("PB", "LCP", "LCPXd", "XDSM", "r_w", "sign_PB")}
    tag = ("all" if D == Dall else str(D)) if a.k >= 0 else "real"
    np.savez_compressed(os.path.join(HERE, "out", f"ON_{a.typ}_k{a.k}_D{tag}.npz"), **{k: np.asarray(v, np.float32) for k, v in keep.items()})
    print(a.typ, a.k, tag, {m: int(bh(o[m].astype(float)).sum()) for m in ("PB", "LCP", "LCPXd", "XDSM")}, {m: round(float((o[m] < 0.05).mean()), 3) for m in ("PB", "LCP", "LCPXd", "XDSM")}, round(time.time() - t0), "s", flush=True)
