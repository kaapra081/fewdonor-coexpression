"""REV4 Amendment 1, R6: tables for the Mantis donor and gene maps (PsychAD P1). python rev4_R6.py --sub P1 -> out/R6_P1.npz"""
import os, sys, argparse, numpy as np, scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("dev", "dstd", "pop", "rev"): sys.path.insert(0, os.path.join(HERE, "..", p))
from methods_v2 import bh, paule_mandel, pearson_fisher  # noqa: E402
from dstd import _ms, design_R  # noqa: E402
from revm import stats_from_cells, methods as rmethods  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--sub", default="P1"); a = ap.parse_args(); si = ["P1", "P2", "P3"].index(a.sub)
OD = os.environ.get("REV4_OUT", os.path.join(HERE, "out")); os.makedirs(OD, exist_ok=True)
P = os.path.join(HERE, "..", "pop"); X = sp.load_npz(os.path.join(P, f"pop_{a.sub}_X.npz")).tocsr(); M = np.load(os.path.join(P, f"pop_{a.sub}_meta.npz"), allow_pickle=True)
lib = M["lib"].astype(float); don = M["donor"].astype(str); sym = M["sym"].astype(str); G = X.shape[1]; iu = np.triu_indices(G, 1); A_, B_ = iu
# donor-centring exactly as REV2 A / R3
ud = np.unique(don); rows = [np.where(don == u)[0] for u in ud]
st = [_ms(X[r].toarray().astype(np.float64), lib[r]) for r in rows]
m = np.array([s_[0] for s_ in st]); ps = np.array([s_[1] for s_ in st]); tau2, mu, _ = paule_mandel(m, ps)
th = mu[None, :] + tau2[None, :] / (tau2[None, :] + ps) * (m - mu[None, :]); f = np.clip(mu[None, :] / np.maximum(th, 1e-30), 0.1, 10.0)
idx = np.empty(X.shape[0], dtype=np.int64)
for i, r in enumerate(rows): idx[r] = i
Xs = sp.csr_matrix(X.multiply(f[idx]))


def profiles(Xm, labs):
    pu = np.unique(labs); cl = [(Xm[np.where(labs == u)[0]].toarray().astype(np.float64), lib[np.where(labs == u)[0]]) for u in pu]
    fu, hs = stats_from_cells(cl, 44_000_001); L = np.log1p(1e6 * fu["m"]); R, _ = design_R(fu["m"], fu["psi"], fu["kappa"], fu["S"])
    nd = np.array([len(s) for _, s in cl], float); dp = np.array([s.mean() for _, s in cl]); return pu, cl, fu, hs, L, R, nd, dp


out = {}
# R6a-i null pseudo-donors (R3 seed, k = 0)
lab = np.random.default_rng(61_000_000 + si).permutation(don)
pu0, cl0, fu0, hs0, L0, R0, nd0, dp0 = profiles(Xs, lab)
pb0 = pearson_fisher(L0, iu); c0 = bh(pb0["p"])
out.update(null_donor=pu0, null_L=L0, null_R=R0, null_nd=nd0, null_depth=dp0, null_rPB=pb0["est"], null_cPB=c0)
# R6a-ii real donors
pu1, cl1, fu1, hs1, L1, R1, nd1, dp1 = profiles(X, don)
pb1 = pearson_fisher(L1, iu); o1 = rmethods(fu1, hs1, iu); cP = bh(pb1["p"]); cX = bh(o1["XDSM"])
out.update(real_donor=pu1, real_L=L1, real_R=R1, real_nd=nd1, real_depth=dp1, real_rPB=pb1["est"], real_pPB=pb1["p"], real_cPB=cP, real_pXDSM=o1["XDSM"], real_sXDSM=o1["sign_XDSM"], real_cXDSM=cX)
# XDS-M signed statistic magnitude for pair-level comparison: z from two-sided p and sign
from scipy.stats import norm  # noqa: E402
out["real_zXDSM"] = o1["sign_XDSM"] * norm.isf(np.clip(o1["XDSM"], 1e-300, 1) / 2)
# within-donor cell-level correlation r_w (real donors with >= 10 cells; proportions y / lib, as R4)
S = np.zeros(len(A_)); N = np.zeros(len(A_))
for Y, s in cl1:
    if len(s) < 10: continue
    Z = Y / s[:, None]; Zc = Z - Z.mean(0); sd = np.sqrt((Zc ** 2).sum(0)); ok = sd > 0
    r = ((Zc.T @ Zc) / np.maximum(np.outer(sd, sd), 1e-300))[iu]; v = ok[A_] & ok[B_]; S[v] += r[v]; N[v] += 1
out["r_w"] = np.where(N > 0, S / np.maximum(N, 1), np.nan); out["r_w_n"] = N
# gene-level descriptors
out["sym"] = sym; out["mean_umi"] = np.asarray(X.mean(0)).ravel(); out["detect"] = np.asarray((X > 0).mean(0)).ravel()
np.savez_compressed(os.path.join(OD, f"R6_{a.sub}.npz"), **out)
print("R6", a.sub, "null PB calls", int(c0.sum()), "real PB", int(cP.sum()), "real XDSM", int(cX.sum()), "rw donors", int(sum(len(s) >= 10 for _, s in cl1)), flush=True)
