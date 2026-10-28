"""Calibration gate (RARE_PREREG.md): LCP and PB on simulated data with S1's donor number and exact per-donor cell numbers.
python rare_sim.py --set S1 --rep0 0 --rep1 100"""
import os, sys, json, argparse
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "dev"))
import sim_v2 as S, methods_v2 as M  # noqa: E402
from run_dev import DESIGNS  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--set", default="S1"); ap.add_argument("--rep0", type=int, default=0); ap.add_argument("--rep1", type=int, default=100)
a = ap.parse_args()
nl = np.array(json.load(open(os.path.join(HERE, "rare_prep_meta.json")))[a.set]["n_list"])
S.make_nd = lambda st, rng: nl.copy()
cfg = {k: v for k, v in DESIGNS["D40"].items() if k not in ("nd_median", "nd_min")}; cfg["D"] = len(nl)
pan = S.make_panel(); st = S.setting(**{**cfg, "nd_median": 40, "nd_min": 4}); tr = S.truth(pan, st); iu = tr["iu"]
od = os.path.join(HERE, f"sim_{a.set}"); os.makedirs(od, exist_ok=True)
np.savez_compressed(os.path.join(od, "truth.npz"), cls=tr["cls"].astype(str), altB=tr["altB"])
for rep in range(a.rep0, a.rep1):
    fn = os.path.join(od, f"rep{rep:04d}.npz")
    if os.path.exists(fn): continue
    seed = 7_000_000 + rep
    sim = S.simulate(pan, st, seed); d = M.donor_stats(sim["Y"], sim["s_hat"], sim["donor_idx"])
    p_lcp = M.lcsp_test(d, iu, seed=seed, studentize=False)["p"]; p_pb = M.pb_logcpm(d, iu)["p"]
    np.savez_compressed(fn, LCP=np.asarray(p_lcp, np.float32), PB=np.asarray(p_pb, np.float32)); print(rep, flush=True)
