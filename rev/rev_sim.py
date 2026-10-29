"""R2 simulations: python rev_sim.py --design L500  (data seeds identical to XU_PREREG / COUPLE_PREREG / FIX)"""
import os, sys, json, argparse, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
for p in ("dev",): sys.path.insert(0, os.path.join(HERE, "..", p))
import sim_v2 as S  # noqa: E402
from run_dev import DESIGNS  # noqa: E402
from revm import stats_from_cells, methods  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--design", required=True); a = ap.parse_args()
LO = ["L200", "L500", "L1000"]; SO = ["LS200", "LS500", "LS1000"]; FO = ["F_D16", "F_D24", "F_D40", "F_S1", "F_S2", "F_D16_n100"]
pan = S.make_panel(seed=4112026); d = a.design; R = 50
if d in LO: k = LO.index(d); cfg = dict(DESIGNS["D40"]); cfg.update(D=int(d[1:]), nd_median=20, nd_min=4); base = 9_500_000
elif d in SO: k = SO.index(d); cfg = dict(DESIGNS["D40"]); cfg.update(D=int(d[2:]), nd_median=20, nd_min=4, kappa=0.3, se_mult=2.0); base = 9_700_000
else:
    k = FO.index(d); base = 9_000_000; R = 100; meta = json.load(open(os.path.join(HERE, "..", "rare", "rare_prep_meta.json")))
    if d in ("F_S1", "F_S2"):
        nl = np.array(meta[d[2:]]["n_list"]); S.make_nd = lambda st, rng: nl.copy()
        cfg = {kk: vv for kk, vv in DESIGNS["D40"].items() if kk not in ("nd_median", "nd_min")}; cfg.update(D=len(nl), nd_median=40, nd_min=4)
    elif d == "F_D16_n100": cfg = dict(DESIGNS["D16"]); cfg.update(nd_median=100)
    else: cfg = dict(DESIGNS[d[2:]])
st = S.setting(**cfg); tr = S.truth(pan, st); iu = tr["iu"]; od = os.path.join(HERE, "runs", d); os.makedirs(od, exist_ok=True)
np.savez_compressed(os.path.join(od, "truth.npz"), cls=tr["cls"].astype(str), altB=tr["altB"])
for rep in range(R):
    fn = os.path.join(od, f"rep{rep:04d}.npz")
    if os.path.exists(fn): continue
    seed = base + 100_003 * k + rep; sim = S.simulate(pan, st, seed); rows = [r for r in sim["donor_idx"] if len(r) >= 4]
    f, h = stats_from_cells([(sim["Y"][r].astype(np.float64), sim["s_hat"][r].astype(np.float64)) for r in rows], seed)
    o = methods(f, h, iu); np.savez_compressed(fn, **{m: o[m].astype(np.float32) for m in ("XDSM", "KX")})
print(d, "done", flush=True)
