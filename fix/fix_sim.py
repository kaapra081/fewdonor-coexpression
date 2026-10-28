"""FIX_PREREG.md Item 2 evaluation. python fix_sim.py --design F_D16 --rep0 0 --rep1 100"""
import os, sys, json, argparse, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "dev")); sys.path.insert(0, HERE)
import sim_v2 as S, methods_v2 as M  # noqa: E402
from run_dev import DESIGNS  # noqa: E402
from lcp2 import lcp_both  # noqa: E402
import torch; torch.set_num_threads(1)
ORDER = ["F_D16", "F_D24", "F_D40", "F_S1", "F_S2", "F_D16_n100"]
meta = json.load(open(os.path.join(HERE, "..", "rare", "rare_prep_meta.json")))
ap = argparse.ArgumentParser(); ap.add_argument("--design", required=True); ap.add_argument("--rep0", type=int, default=0); ap.add_argument("--rep1", type=int, default=100)
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
np.savez_compressed(os.path.join(od, "truth.npz"), cls=tr["cls"].astype(str), altB=tr["altB"])
for rep in range(a.rep0, a.rep1):
    fn = os.path.join(od, f"rep{rep:04d}.npz")
    if os.path.exists(fn): continue
    t0 = time.time(); seed = 9_000_000 + 100_003 * k + rep
    sim = S.simulate(pan, st, seed); d = M.donor_stats(sim["Y"], sim["s_hat"], sim["donor_idx"])
    rows = [r for r in sim["donor_idx"] if len(r) >= 4]
    cells = [(sim["Y"][r].astype(np.float64), sim["s_hat"][r].astype(np.float64)) for r in rows]
    o = lcp_both(d["m"], d["psi"], d["C"].astype(np.float64), d["kappa"], cells, iu, B=100_000, seed=seed + 5, paule_mandel=M.paule_mandel)
    pb = M.pb_logcpm(d, iu)["p"]
    np.savez_compressed(fn, LCP=o["p_LCP"].astype(np.float32), LCP2=o["p_LCP2"].astype(np.float32), PB=np.asarray(pb, np.float32))
    print(a.design, rep, round(time.time() - t0, 1), flush=True)
