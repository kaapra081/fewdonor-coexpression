"""Minimal example on simulated data: 40 donors, ~40 cells each."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..")); sys.path.insert(0, os.path.join(HERE, "..", "dev"))
import sim_v2 as S
from lcp import donor_statistics, lcp_test, bh

panel = S.make_panel(); st = S.setting(D=40); tr = S.truth(panel, st)
sim = S.simulate(panel, st, seed=1)
donor = np.concatenate([[d] * len(r) for d, r in enumerate(sim["donor_idx"])])
stats_ = donor_statistics(sim["Y"], sim["s_hat"], donor)        # Y: cells x genes counts; s: total UMI per cell
res = lcp_test(stats_, B=20000, seed=1)
calls = bh(res["p"])
print(f"{calls.sum()} pairs called at FDR 0.05; true between-donor pairs: {tr['altB'].sum()}; "
      f"false calls: {(calls & ~tr['altB']).sum()}")
