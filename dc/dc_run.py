"""python dc_run.py --sub P1 --k 10   |  --k -1 for real donors"""
import os, sys, argparse, time, numpy as np, scipy.sparse as sp, anndata as ad
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("dstd", "pop", "dev"): sys.path.insert(0, os.path.join(HERE, "..", p))
from dstd import donor_level_cells, design_R  # noqa: E402
from pop_core import std_R, bh  # noqa: E402
from methods_v2 import paule_mandel  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--sub", required=True); ap.add_argument("--k", type=int, required=True); a = ap.parse_args()
si = ["P1", "P2", "P3"].index(a.sub); P = os.path.join(HERE, "..", "pop")
X = sp.load_npz(os.path.join(P, f"pop_{a.sub}_X.npz")).tocsr(); M = np.load(os.path.join(P, f"pop_{a.sub}_meta.npz"), allow_pickle=True)
lib = M["lib"]; don = M["donor"].astype(str); G = X.shape[1]; iu = np.triu_indices(G, 1); A_, B_ = iu
obs = ad.read_h5ad("/Users/kaavinprasanna/Documents/recomb27/data/pfc_myeloid_1469.h5ad", backed="r").obs
adlab = obs.drop_duplicates("donor_id").set_index("donor_id")["AD"].astype(str).to_dict()
keep = np.array([adlab.get(u, "nan") in ("0.0", "1.0") for u in don]); X = X[np.where(keep)[0]]; lib = lib[keep]; don = don[keep]
t0 = time.time()
if a.k >= 0:
    seed = 41_000_000 + 10 * a.k + si; lab = np.random.default_rng(seed).permutation(don)   # slot u keeps n_u; label of u
else:
    seed = 43_000_000 + si; lab = don
ud = np.unique(lab); grp = np.array([adlab[u] == "1.0" for u in ud])
cells = ((X[np.where(lab == u)[0]].toarray().astype(np.float64), lib[np.where(lab == u)[0]].astype(np.float64)) for u in ud)
d = donor_level_cells(cells, seed + 7)
def welch(Q):
    q1, q0 = Q[grp], Q[~grp]
    z = (q1.mean(0) - q0.mean(0)) / np.sqrt(q1.var(0, ddof=1) / len(q1) + q0.var(0, ddof=1) / len(q0)); return 2 * norm.sf(np.abs(z))
out = {}
L = np.log1p(1e6 * d["m"])
def pear(Mx):
    Mc = Mx - Mx.mean(0); cov = Mc.T @ Mc / (len(Mx) - 1); sd = np.sqrt(np.diag(cov)); return np.clip((cov / np.outer(sd, sd))[iu], -0.999999, 0.999999)
r1, r0 = pear(L[grp]), pear(L[~grp]); n1, n0 = grp.sum(), (~grp).sum()
out["PB"] = 2 * norm.sf(np.abs((np.arctanh(r1) - np.arctanh(r0)) / np.sqrt(1 / (n1 - 3) + 1 / (n0 - 3))))
R = std_R(d["m"], d["psi"]); tau2, mu, _ = paule_mandel(d["m"], d["psi"]); v = tau2[None, :] + d["psi"]
Z = np.empty((len(ud), len(A_)), np.float32)
for c0 in range(0, len(A_), 20000):
    aa, bb = A_[c0:c0 + 20000], B_[c0:c0 + 20000]
    Z[:, c0:c0 + 20000] = R[:, aa] * R[:, bb] - d["kappa"][:, None] * d["C"][:, aa, bb] / np.sqrt(v[:, aa] * v[:, bb])
out["LCP"] = welch(Z.astype(np.float64)); del Z
RA, _ = design_R(d["mA"], d["psiA"], d["kappaA"], d["SA"]); RB, _ = design_R(d["mB"], d["psiB"], d["kappaB"], d["SB"])
out["XDS"] = welch(RA[:, A_] * RB[:, B_])
Cw = d["C"][:, A_, B_].astype(np.float64).mean(0); sig2 = np.maximum((d["psi"] - d["m"] / d["S"][:, None]) / d["kappa"][:, None], 0).mean(0)
with np.errstate(invalid="ignore", divide="ignore"): out["r_w"] = Cw / np.sqrt(sig2[A_] * sig2[B_])
out["n1"] = np.array([n1]); out["n0"] = np.array([n0])
tag = f"k{a.k}" if a.k >= 0 else "real"
np.savez_compressed(os.path.join(HERE, "out", f"DC_{a.sub}_{tag}.npz"), **{k: np.asarray(v_, np.float32) for k, v_ in out.items()})
print(a.sub, tag, "AD", n1, "ctl", n0, {m: int(bh(out[m].astype(float)).sum()) for m in ("PB", "LCP", "XDS")},
      {m: round(float((out[m] < 0.05).mean()), 3) for m in ("PB", "LCP", "XDS")}, round(time.time() - t0), "s", flush=True)
