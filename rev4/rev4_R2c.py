"""REV4 Amendment 2 (post hoc) R2c: bank-centred XDS-M and pseudobulk on the R2 reassignments. python rev4_R2c.py --sub P2 --k0 0 --k1 40"""
import os, sys, json, argparse, time, numpy as np, scipy.sparse as sp
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("dev", "dstd", "pop", "rev"): sys.path.insert(0, os.path.join(HERE, "..", p))
from methods_v2 import bh, paule_mandel, pearson_fisher  # noqa: E402
from dstd import _ms, design_R  # noqa: E402
from revm import stats_from_cells  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--sub", required=True); ap.add_argument("--k0", type=int, default=0); ap.add_argument("--k1", type=int, default=40); a = ap.parse_args()
si = ["P1", "P2", "P3"].index(a.sub); OD = os.path.join(HERE, "out"); P = os.path.join(HERE, "..", "pop")
X = sp.load_npz(os.path.join(P, f"pop_{a.sub}_X.npz")).tocsr(); M = np.load(os.path.join(P, f"pop_{a.sub}_meta.npz"), allow_pickle=True)
lib = M["lib"].astype(float); don = M["donor"].astype(str); G = X.shape[1]; iu = np.triu_indices(G, 1); A_, B_ = iu
src = json.load(open(os.path.join(HERE, "donor_source.json"))); bank = np.array([src[d_] for d_ in don])
ud = np.unique(don); rows = [np.where(don == u)[0] for u in ud]   # donor-centring identical to rev4.py centre()
st = [_ms(X[r].toarray().astype(np.float64), lib[r]) for r in rows]
m = np.array([s_[0] for s_ in st]); ps = np.array([s_[1] for s_ in st]); tau2, mu, _ = paule_mandel(m, ps)
th = mu[None, :] + tau2[None, :] / (tau2[None, :] + ps) * (m - mu[None, :]); f = np.clip(mu[None, :] / np.maximum(th, 1e-30), 0.1, 10.0)
idx = np.empty(X.shape[0], dtype=np.int64)
for i, r in enumerate(rows): idx[r] = i
Xs = sp.csr_matrix(X.multiply(f[idx]))


def tt(q):
    t = q.mean(0) / (q.std(0, ddof=1) / np.sqrt(q.shape[0])); return 2 * norm.sf(np.abs(t))


def centre_by(Z, g):
    Z = Z.copy()
    for b in np.unique(g): Z[g == b] -= Z[g == b].mean(0)
    return Z


for k in range(a.k0, a.k1):
    fn = os.path.join(OD, f"R2c_{a.sub}_k{k}.json")
    if os.path.exists(fn): continue
    t0 = time.time(); seed = 64_000_000 + 10 * k + si; rng = np.random.default_rng(seed); lab = don.copy()
    for bk in ("H", "M", "R"):   # identical draw order to rev4.py R2
        ix = np.where(bank == bk)[0]; lab[ix] = rng.permutation(don[ix])
    pu = np.unique(lab); gb = np.array([src[u] for u in pu])
    cells = [(Xs[np.where(lab == u)[0]].toarray(), lib[np.where(lab == u)[0]]) for u in pu]
    fs, hs = stats_from_cells(cells, seed); D = len(pu)
    L = np.log1p(1e6 * fs["m"]); Lb = centre_by(L, gb)
    pv = {"PB": pearson_fisher(L, iu)["p"]}
    Cb = Lb.T @ Lb; sd = np.sqrt(np.diag(Cb)); rb = np.clip((Cb / np.outer(sd, sd))[iu], -0.999999, 0.999999)
    pv["PBb"] = 2 * norm.sf(np.abs(np.arctanh(rb) * np.sqrt(D - 3 - 2)))
    Rs = [(design_R(h["mA"], h["psiA"], h["kappaA"], h["SA"])[0], design_R(h["mB"], h["psiB"], h["kappaB"], h["SB"])[0]) for h in hs]
    pM = np.empty(len(A_)); pMb = np.empty(len(A_))
    for c0 in range(0, len(A_), 10000):
        aa, bb = A_[c0:c0 + 10000], B_[c0:c0 + 10000]
        q = sum(RA[:, aa] * RB[:, bb] for RA, RB in Rs) / len(Rs); pM[c0:c0 + 10000] = tt(q)
        qb = sum(centre_by(RA[:, aa], gb) * centre_by(RB[:, bb], gb) for RA, RB in Rs) / len(Rs); pMb[c0:c0 + 10000] = tt(qb)
    pv["XDSM"] = pM; pv["XDSMb"] = pMb
    res = {mm: dict(V=int(bh(v).sum()), t05=float((v < 0.05).mean())) for mm, v in pv.items()}
    sdL = L.std(0); res["bank_shift_median"] = {b: float(np.median(np.abs(L[gb == b].mean(0) - L.mean(0)) / np.maximum(sdL, 1e-12))) for b in np.unique(gb)}
    res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w"))
    print("R2c", a.sub, k, {mm: res[mm]["V"] for mm in pv}, res["bank_shift_median"], res["sec"], "s", flush=True)
