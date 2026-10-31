"""REV3_PREREG B2: covariate-adjusted (PBres) and precision-weighted (PBw) pseudobulk on the donor-centred null. python rev3_B2.py --sub P1"""
import os, sys, json, argparse, numpy as np, scipy.sparse as sp
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("dev", "dstd", "pop"): sys.path.insert(0, os.path.join(HERE, "..", p))
from methods_v2 import bh, paule_mandel, pearson_fisher  # noqa: E402
from dstd import _ms  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--sub", required=True); ap.add_argument("--k0", type=int, default=0); ap.add_argument("--k1", type=int, default=20); a = ap.parse_args()
si = ["P1", "P2", "P3"].index(a.sub); P = os.path.join(HERE, "..", "pop"); OD = os.path.join(HERE, "out"); os.makedirs(OD, exist_ok=True)
X = sp.load_npz(os.path.join(P, f"pop_{a.sub}_X.npz")).tocsr(); M = np.load(os.path.join(P, f"pop_{a.sub}_meta.npz"), allow_pickle=True)
lib = M["lib"].astype(float); don = M["donor"].astype(str); G = X.shape[1]; iu = np.triu_indices(G, 1)
ud = np.unique(don); rows = [np.where(don == u)[0] for u in ud]
st = [_ms(X[r].toarray().astype(np.float64), lib[r]) for r in rows]
m = np.array([s_[0] for s_ in st]); ps = np.array([s_[1] for s_ in st]); tau2, mu, _ = paule_mandel(m, ps)
th = mu[None, :] + tau2[None, :] / (tau2[None, :] + ps) * (m - mu[None, :]); f = np.clip(mu[None, :] / np.maximum(th, 1e-30), 0.1, 10.0)
idx = np.empty(X.shape[0], dtype=np.int64)
for i, r in enumerate(rows): idx[r] = i
Xs = sp.csr_matrix(X.multiply(f[idx]))
def fz(r, neff): return 2 * norm.sf(np.abs(np.arctanh(np.clip(r, -0.999999, 0.999999)) * np.sqrt(neff)))
for k in range(a.k0, a.k1):
    fn = os.path.join(OD, f"B2_{a.sub}_k{k}.json")
    if os.path.exists(fn): continue
    seed = 61_000_000 + 10 * k + si; lab = np.random.default_rng(seed).permutation(don); pu = np.unique(lab)
    mm = []; nd = []; dp = []
    for u in pu:
        r = np.where(lab == u)[0]; Y = np.asarray(Xs[r].sum(0)).ravel(); S = lib[r].sum(); mm.append(Y / S); nd.append(len(r)); dp.append(lib[r].mean())
    L = np.log1p(1e6 * np.array(mm)); nd = np.array(nd, float); dp = np.array(dp); D = len(pu)
    Z = np.column_stack([np.ones(D), np.log(nd), np.log(dp)]); Rres = L - Z @ np.linalg.lstsq(Z, L, rcond=None)[0]
    Rc = Rres - Rres.mean(0); C = Rc.T @ Rc; sd = np.sqrt(np.diag(C)); rres = (C / np.outer(sd, sd))[iu]
    w = nd / nd.sum(); Lw = L - (w[:, None] * L).sum(0); Cw = (Lw * w[:, None]).T @ Lw; sw = np.sqrt(np.diag(Cw)); rw = (Cw / np.outer(sw, sw))[iu]
    neff = nd.sum() ** 2 / (nd ** 2).sum()
    pv = dict(PB=pearson_fisher(L, iu)["p"], PBres=fz(rres, D - 5), PBw=fz(rw, neff - 3))
    res = {m_: dict(V=int(bh(v).sum()), t05=float((v < 0.05).mean())) for m_, v in pv.items()}; res["D"] = D; res["neff"] = float(neff)
    json.dump(res, open(fn, "w")); print("B2", a.sub, k, {m_: res[m_]["V"] for m_ in pv}, flush=True)
