"""REV2_PREREG F: A5 stress (donor-specific overdispersion shared by genes). python rev2_F.py --design L200"""
import os, sys, inspect, argparse, json, numpy as np, torch
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("dev", "dstd", "pop", "rev"): sys.path.insert(0, os.path.join(HERE, "..", p))
import sim_v2 as S  # noqa: E402
from run_dev import DESIGNS  # noqa: E402
from dstd import donor_level_cells, methods_d  # noqa: E402
from pop_core import run_methods  # noqa: E402
from revm import stats_from_cells, methods as rmethods  # noqa: E402
torch.set_num_threads(1)
src = inspect.getsource(S.simulate)
old = "xi = rng.gamma(phi, 1.0 / phi, size=(n, G))"
assert old in src
src = src.replace(old, 'odf = np.exp(st.get("od_sd", 0.0) * rng.normal()); xi = rng.gamma(phi * odf, 1.0 / (phi * odf), size=(n, G))')
exec(src, S.__dict__)
ap = argparse.ArgumentParser(); ap.add_argument("--design", required=True); a = ap.parse_args()
k = ["L200", "L500"].index(a.design); cfg = dict(DESIGNS["D40"]); cfg.update(D=int(a.design[1:]), nd_median=20, nd_min=4)
pan = S.make_panel(seed=4112026); st = S.setting(**cfg); st["od_sd"] = 0.7; tr = S.truth(pan, st); iu = tr["iu"]
od = os.path.join(HERE, "out", "F_" + a.design); os.makedirs(od, exist_ok=True)
np.savez_compressed(os.path.join(od, "truth.npz"), cls=tr["cls"].astype(str), altB=tr["altB"])
for rep in range(50):
    fn = os.path.join(od, f"rep{rep:04d}.npz")
    if os.path.exists(fn): continue
    seed = 9_900_000 + 100_003 * k + rep; sim = S.simulate(pan, st, seed); rows = [r for r in sim["donor_idx"] if len(r) >= 4]
    cells = [(sim["Y"][r].astype(np.float64), sim["s_hat"][r].astype(np.float64)) for r in rows]
    d = donor_level_cells(cells, seed + 7); r1 = run_methods(d, iu, seed=seed); r2 = methods_d(d, iu, seed=seed + 3)
    f_, hs = stats_from_cells(cells, seed); r3 = rmethods(f_, hs, iu)
    np.savez_compressed(fn, PB=r1["PB"].astype(np.float32), LCP=np.asarray(r1["LCP"], np.float32), LCPd=np.asarray(r2["LCPd"], np.float32),
                        XDS=r3["XDS"].astype(np.float32), XDSM=r3["XDSM"].astype(np.float32), KX=r3["KX"].astype(np.float32))
print(a.design, "done", flush=True)
