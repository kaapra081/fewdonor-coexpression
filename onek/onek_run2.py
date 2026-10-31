"""REV2_PREREG E: python onek_run2.py --typ B1 --k0 3 --k1 33 (full D; PB, PBlin, LCP, XDS-M)"""
import os, sys, json, argparse, time, numpy as np, scipy.sparse as sp, torch
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("dstd", "pop", "rev", "dev"): sys.path.insert(0, os.path.join(HERE, "..", p))
from pop_core import run_methods, bh  # noqa: E402
from dstd import donor_level_cells  # noqa: E402
from revm import stats_from_cells, methods as rmethods  # noqa: E402
from methods_v2 import pearson_fisher  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--typ", required=True); ap.add_argument("--k0", type=int, default=3); ap.add_argument("--k1", type=int, default=33); a = ap.parse_args(); torch.set_num_threads(1)
ti = ["B1", "B2", "B3"].index(a.typ)
X = sp.load_npz(os.path.join(HERE, f"onek_{a.typ}_X.npz")).tocsr(); M = np.load(os.path.join(HERE, f"onek_{a.typ}_meta.npz"), allow_pickle=True)
lib = M["lib"]; don = M["donor"].astype(str); iu = np.triu_indices(X.shape[1], 1)
for k in range(a.k0, a.k1):
    fn = os.path.join(HERE, "out", f"E_{a.typ}_k{k}.json")
    if os.path.exists(fn): continue
    t0 = time.time(); seed = 45_000_000 + 10 * k + ti; lab = np.random.default_rng(seed).permutation(don)
    ud = np.unique(lab); cells = [(X[np.where(lab == u)[0]].toarray().astype(np.float64), lib[np.where(lab == u)[0]].astype(np.float64)) for u in ud]
    d = donor_level_cells(cells, seed + 7); d["donors"] = ud; o = run_methods(d, iu, seed=seed)
    f_, hs = stats_from_cells(cells, seed + 11); r3 = rmethods(f_, hs, iu)
    pv = dict(PB=o["PB"], PBlin=pearson_fisher(d["m"], iu)["p"], LCP=o["LCP"], XDSM=r3["XDSM"])
    res = {m: dict(V=int(bh(np.asarray(v, float)).sum()), t05=float((np.asarray(v) < 0.05).mean())) for m, v in pv.items()}
    json.dump(res, open(fn, "w")); print(a.typ, k, {m: res[m]["V"] for m in res}, round(time.time() - t0), "s", flush=True)
