"""E3 (and E4 with --fix): python xu_lsim.py --design L500"""
import os, sys, json, argparse, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("dev", "dstd", "pop"): sys.path.insert(0, os.path.join(HERE, "..", p))
sys.path.insert(0, HERE)
import sim_v2 as S  # noqa: E402
from run_dev import DESIGNS  # noqa: E402
from dstd import donor_level_cells, methods_d  # noqa: E402
from pop_core import run_methods  # noqa: E402
from xu import xu  # noqa: E402
import torch; torch.set_num_threads(1)
ap = argparse.ArgumentParser(); ap.add_argument("--design", required=True); ap.add_argument("--rep0", type=int, default=0); ap.add_argument("--rep1", type=int, default=50)
a = ap.parse_args()
LORDER = ["L200", "L500", "L1000"]; SORDER = ["LS200", "LS500", "LS1000"]; FORDER = ["F_D16", "F_D24", "F_D40", "F_S1", "F_S2", "F_D16_n100"]
pan = S.make_panel(seed=4112026); od = os.path.join(HERE, "runs", a.design); os.makedirs(od, exist_ok=True)
if a.design in LORDER:
    k = LORDER.index(a.design); cfg = dict(DESIGNS["D40"]); cfg.update(D=int(a.design[1:]), nd_median=20, nd_min=4); base = 9_500_000
elif a.design in SORDER:
    k = SORDER.index(a.design); cfg = dict(DESIGNS["D40"]); cfg.update(D=int(a.design[2:]), nd_median=20, nd_min=4, kappa=0.3, se_mult=2.0); base = 9_700_000
else:
    k = FORDER.index(a.design); base = 9_000_000
    meta = json.load(open(os.path.join(HERE, "..", "rare", "rare_prep_meta.json")))
    if a.design in ("F_S1", "F_S2"):
        nl = np.array(meta[a.design[2:]]["n_list"]); S.make_nd = lambda st, rng: nl.copy()
        cfg = {kk: vv for kk, vv in DESIGNS["D40"].items() if kk not in ("nd_median", "nd_min")}; cfg.update(D=len(nl), nd_median=40, nd_min=4)
    elif a.design == "F_D16_n100":
        cfg = dict(DESIGNS["D16"]); cfg.update(nd_median=100)
    else:
        cfg = dict(DESIGNS[a.design[2:]])
st = S.setting(**cfg); tr = S.truth(pan, st); iu = tr["iu"]
np.savez_compressed(os.path.join(od, "truth.npz"), cls=tr["cls"].astype(str), altB=tr["altB"])
for rep in range(a.rep0, a.rep1):
    fn = os.path.join(od, f"rep{rep:04d}.npz")
    if os.path.exists(fn): continue
    t0 = time.time(); seed = base + 100_003 * k + rep
    sim = S.simulate(pan, st, seed); rows = [r for r in sim["donor_idx"] if len(r) >= 4]
    cells = [(sim["Y"][r].astype(np.float64), sim["s_hat"][r].astype(np.float64)) for r in rows]
    o = {}
    x = xu(cells, *[None] * 0, **{}) if False else None
    d = donor_level_cells(cells, seed + 7, need_C=a.design in LORDER + SORDER)
    o["XU"] = xu(cells, d["m"], d["psi"], d["kappa"], d["S"], iu)["p"]
    if a.design in LORDER + SORDER:
        r1 = run_methods(d, iu, seed=seed); r2 = methods_d(d, iu, seed=seed + 3)
        o.update(PB=r1["PB"], LCP=r1["LCP"], XDS=r2["LCPXd"], LCPd=r2["LCPd"])
    np.savez_compressed(fn, **{m: np.asarray(v, np.float32) for m, v in o.items()})
    print(a.design, rep, round(time.time() - t0, 1), flush=True)
