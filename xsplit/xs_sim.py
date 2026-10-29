"""python xs_sim.py --design F_D16 --rep0 0 --rep1 100  (same data seeds as fix/fix_sim.py)"""
import os, sys, json, argparse, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "dev")); sys.path.insert(0, HERE)
import sim_v2 as S, methods_v2 as M  # noqa: E402
from run_dev import DESIGNS  # noqa: E402
from lcpx import lcpx  # noqa: E402
import torch; torch.set_num_threads(1)
ORDER = ["F_D16", "F_D24", "F_D40", "F_S1", "F_S2", "F_D16_n100"]
meta = json.load(open(os.path.join(HERE, "..", "rare", "rare_prep_meta.json")))
ap = argparse.ArgumentParser(); ap.add_argument("--design", required=True); ap.add_argument("--rep0", type=int, default=0); ap.add_argument("--rep1", type=int, default=100)
ap.add_argument("--B", type=int, default=100_000)
a = ap.parse_args(); k = ORDER.index(a.design)
if a.design in ("F_S1", "F_S2"):
    nl = np.array(meta[a.design[2:]]["n_list"]); S.make_nd = lambda st, rng: nl.copy()
    cfg = {kk: vv for kk, vv in DESIGNS["D40"].items() if kk not in ("nd_median", "nd_min")}; cfg.update(D=len(nl), nd_median=40, nd_min=4)
elif a.design == "F_D16_n100":
    cfg = dict(DESIGNS["D16"]); cfg.update(nd_median=100)
else:
    cfg = dict(DESIGNS[a.design[2:]])
pan = S.make_panel(seed=4112026); st = S.setting(**cfg); tr = S.truth(pan, st); iu = tr["iu"]
od = os.path.join(HERE, "runs", a.design); os.makedirs(od, exist_ok=True)
ref = np.load(os.path.join(HERE, "..", "fix", "runs", a.design, "truth.npz"))
assert (ref["altB"] == tr["altB"]).all()
for rep in range(a.rep0, a.rep1):
    fn = os.path.join(od, f"rep{rep:04d}.npz")
    if os.path.exists(fn): continue
    t0 = time.time(); seed = 9_000_000 + 100_003 * k + rep
    sim = S.simulate(pan, st, seed)
    rows = [r for r in sim["donor_idx"] if len(r) >= 4]
    cells = [(sim["Y"][r].astype(np.float64), sim["s_hat"][r].astype(np.float64)) for r in rows]
    o = lcpx(cells, iu, B=a.B, seed=seed, paule_mandel=M.paule_mandel)
    np.savez_compressed(fn, LCPX=o["p_LCPX"].astype(np.float32), LCPXS=o["p_LCPXS"].astype(np.float32))
    print(a.design, rep, round(time.time() - t0, 1), flush=True)
