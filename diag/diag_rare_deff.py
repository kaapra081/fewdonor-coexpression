"""Median D_eff/D and Proposition-1 efficiency for the rare-supertype calibration designs (S1, S2), same computation as diag_sim.py."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "dev"))
import sim_v2 as S  # noqa: E402
from run_dev import DESIGNS  # noqa: E402
meta = json.load(open(os.path.join(HERE, "..", "rare", "rare_prep_meta.json"))); out = {}
for key in ("S1", "S2"):
    nl = np.array(meta[key]["n_list"]); S.make_nd = lambda st, rng, nl=nl: nl.copy()
    cfg = {k: v for k, v in DESIGNS["D40"].items() if k not in ("nd_median", "nd_min")}; cfg["D"] = len(nl)
    pan = S.make_panel(); st = S.setting(**{**cfg, "nd_median": 40, "nd_min": 4}); tr = S.truth(pan, st); iu = tr["iu"]; alt = tr["altB"]; sB = np.diag(tr["SB"])
    are, deff = [], []
    for rep in range(10):
        sim = S.simulate(pan, st, 7_000_000 + rep); rho = sB[None, :] / (sB[None, :] + sim["psi_true"])
        pr = rho[:, iu[0][alt]] * rho[:, iu[1][alt]]
        are.append(np.median(np.sqrt(pr).mean(0) ** 2 / pr.mean(0))); deff.append(np.median(pr.mean(0)))
    out[key] = dict(ARE_median=float(np.mean(are)), Deff_over_D_median=float(np.mean(deff))); print(key, out[key])
json.dump(out, open(os.path.join(HERE, "diag_rare_deff.json"), "w"), indent=1)
