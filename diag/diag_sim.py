"""Post hoc diagnostics on simulated designs.
(1) Proposition 1: LCP's asymptotic efficiency relative to the Theorem-1 information bound, ARE_jk = [mean_d sqrt(rho_dj rho_dk)]^2 / mean_d rho_dj rho_dk,
    and D_eff/D = mean_d rho_dj rho_dk, from the true sampling variances (median over between-donor pairs; 10 replicates per design).
(2) Depth-state coupling statistic: pooled within-donor Pearson correlation between log library size and x_ij = Y_ij / s_i, per gene."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "dev"))
import sim_v2 as S  # noqa: E402
from run_dev import DESIGNS, DESIGN_INDEX  # noqa: E402
from run_heldout import HO, HO_INDEX  # noqa: E402


def coupling(Y, s, rows):
    ls = np.log(s); x = Y / s[:, None]
    lc = np.concatenate([ls[r] - ls[r].mean() for r in rows]); xc = np.vstack([x[r] - x[r].mean(0) for r in rows])
    return (lc @ xc) / np.sqrt((lc ** 2).sum() * (xc ** 2).sum(0))


out = {}
pan = S.make_panel(); pan_ho = None
for name, cfg, base, k, panel_seed in [(d, DESIGNS[d], 2_000_000, DESIGN_INDEX[d], None) for d in ("D16", "D24", "D40", "D16_STRESS_depth")] + \
                                   [(d, HO[d], 5_000_000, HO_INDEX[d], 4112026) for d in ("HO_seaad24", "HO_seaad24_cap40", "HO_lee40", "HO_D32")]:
    p = pan if panel_seed is None else (pan_ho if pan_ho is not None else S.make_panel(seed=panel_seed))
    if panel_seed is not None: pan_ho = p
    st = S.setting(**cfg); tr = S.truth(p, st); iu = tr["iu"]; alt = tr["altB"]; sB = np.diag(tr["SB"])
    are, deff, cpl_W, cpl_other = [], [], [], []
    for rep in range(10):
        sim = S.simulate(p, st, base + (100_019 if panel_seed else 100_003) * k + rep)
        rho = sB[None, :] / (sB[None, :] + sim["psi_true"])
        a, b = iu[0][alt], iu[1][alt]
        pr = rho[:, a] * rho[:, b]
        are.append(np.median(np.sqrt(pr).mean(0) ** 2 / pr.mean(0))); deff.append(np.median(pr.mean(0)))
        if rep < 3:
            c = coupling(sim["Y"].astype(float), sim["s_hat"], sim["donor_idx"])
            lab = p["labels"]; cpl_W.append(np.median(c[lab == "W"])); cpl_other.append(np.median(np.abs(c[lab != "W"])))
    out[name] = dict(ARE_median=float(np.mean(are)), Deff_over_D_median=float(np.mean(deff)),
                     coupling_W_median=float(np.mean(cpl_W)), coupling_absother_median=float(np.mean(cpl_other)))
    print(name, {k: round(v, 3) for k, v in out[name].items()}, flush=True)
json.dump(out, open(os.path.join(HERE, "diag_sim.json"), "w"), indent=1)
