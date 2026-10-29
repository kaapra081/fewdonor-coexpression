"""Post hoc (disclosed): XDS on real PsychAD donors with 3 split seeds; Jaccard of BH call sets; sign of PB false calls in POP null."""
import os, sys, json, glob, numpy as np, scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("dstd", "pop", "dev"): sys.path.insert(0, os.path.join(HERE, "..", p))
from dstd import donor_level_cells, methods_d  # noqa: E402
from pop_core import bh  # noqa: E402
out = {}
for sub in ("P1", "P2", "P3"):
    P = os.path.join(HERE, "..", "pop"); X = sp.load_npz(os.path.join(P, f"pop_{sub}_X.npz")).tocsr()
    M = np.load(os.path.join(P, f"pop_{sub}_meta.npz"), allow_pickle=True); lib = M["lib"]; don = M["donor"].astype(str)
    iu = np.triu_indices(X.shape[1], 1); ud = np.unique(don); rows = [np.where(don == u)[0] for u in ud]
    calls = []
    for s in (1, 2, 3):
        d = donor_level_cells(((X[r].toarray().astype(np.float64), lib[r].astype(np.float64)) for r in rows), 44_000_000 + s, need_C=False)
        d["C"] = np.zeros((len(ud), 1, 1), np.float32)
        o = methods_d.__wrapped__(d, iu) if hasattr(methods_d, "__wrapped__") else None
        from dstd import design_R
        from pop_core import normal_perm_p
        RA, _ = design_R(d["mA"], d["psiA"], d["kappaA"], d["SA"]); RB, _ = design_R(d["mB"], d["psiB"], d["kappaB"], d["SB"])
        T = (RA[:, iu[0]] * RB[:, iu[1]]).sum(0); p = normal_perm_p(RA, RB, iu[0], iu[1], T); calls.append(bh(p))
    J = [float((calls[i] & calls[j]).sum() / max((calls[i] | calls[j]).sum(), 1)) for i, j in ((0, 1), (0, 2), (1, 2))]
    out[sub] = dict(n_calls=[int(c.sum()) for c in calls], jaccard=J)
    # sign of PB false calls in the POP null (full D, K = 3)
    pos = []
    for f in sorted(glob.glob(os.path.join(HERE, "..", "pop", "out", f"AD_{sub}_k*_Dall.npz"))):
        z = np.load(f); c = bh(z["PB"].astype(float)); pos.append(float((z["sign_PB"][c] > 0).mean()))
    out[sub]["PB_false_frac_positive"] = pos
    print(sub, out[sub], flush=True)
json.dump(out, open(os.path.join(HERE, "split_results.json"), "w"), indent=1)
