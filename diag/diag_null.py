"""Null references for diag_real.py: (a) PC1 share of standardized donor means after independently permuting each gene across donors
(100 permutations); (b) within-donor SD of log library size in simulated null (D16) and stress (D16_STRESS_depth) designs."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "dev"))
from methods_v2 import paule_mandel  # noqa: E402
import sim_v2 as S  # noqa: E402
from run_dev import DESIGNS, DESIGN_INDEX  # noqa: E402
out = {}; rng = np.random.default_rng(1)
for coh in ("seaad", "lee"):
    z = np.load(os.path.join(os.environ.get("FEWDONOR_STATS", os.path.join(HERE, "..", "real")), f"real_stats_{coh}.npz"), allow_pickle=True)
    tau2, mu, _ = paule_mandel(z["m"], z["psi"]); R = (z["m"] - mu) / np.sqrt(tau2[None, :] + z["psi"])
    sh = []
    for _ in range(100):
        P = np.column_stack([rng.permutation(R[:, g]) for g in range(R.shape[1])]); P = P - P.mean(0)
        s = np.linalg.svd(P, compute_uv=False); sh.append(s[0] ** 2 / (s ** 2).sum())
    out[coh] = dict(pc1_null_mean=float(np.mean(sh)), pc1_null_max=float(np.max(sh)))
pan = S.make_panel()
for d in ("D16", "D16_STRESS_depth"):
    sim = S.simulate(pan, S.setting(**DESIGNS[d]), 2_000_000 + 100_003 * DESIGN_INDEX[d])
    ls = np.log(sim["s_hat"]); out[d] = dict(sd_loglib_median=float(np.median([ls[r].std(ddof=1) for r in sim["donor_idx"]])))
print(json.dumps(out)); json.dump(out, open(os.path.join(HERE, "diag_null.json"), "w"), indent=1)
