"""REV2_PREREG D1/D3: python rev2_D.py --sub P1 --k 10 (donor-centred null) | --k -1 (real) | --k -2 (real, split halves)"""
import os, sys, json, argparse, time, numpy as np, scipy.sparse as sp
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("dev", "dstd", "rev"): sys.path.insert(0, os.path.join(HERE, "..", p))
from methods_v2 import bh, paule_mandel  # noqa: E402
from dstd import _ms, design_R  # noqa: E402
from revm import stats_from_cells  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--sub", required=True); ap.add_argument("--k", type=int, required=True); a = ap.parse_args()
si = ["P1", "P2", "P3"].index(a.sub); P = os.path.join(HERE, "..", "pop")
X = sp.load_npz(os.path.join(P, f"pop_{a.sub}_X.npz")).tocsr(); M = np.load(os.path.join(P, f"pop_{a.sub}_meta.npz"), allow_pickle=True)
lib = M["lib"]; don = M["donor"].astype(str); iu = np.triu_indices(X.shape[1], 1); A_, B_ = iu
_aj = os.path.join(HERE, "adlab.json")   # cached donor -> AD label (same values as the h5ad obs column 'AD')
if os.path.exists(_aj): adlab = json.load(open(_aj))
else:
    import anndata as ad
    obs = ad.read_h5ad("/Users/kaavinprasanna/Documents/recomb27/data/pfc_myeloid_1469.h5ad", backed="r").obs
    adlab = obs.drop_duplicates("donor_id").set_index("donor_id")["AD"].astype(str).to_dict()
keep = np.array([adlab.get(u, "nan") in ("0.0", "1.0") for u in don]); X = X[np.where(keep)[0]]; lib = lib[keep]; don = don[keep]
if a.k >= 0:   # donor-centred cells, then pseudo-donors that keep each real donor's cell number and AD label
    ud0 = np.unique(don); rows0 = [np.where(don == u)[0] for u in ud0]
    st = [_ms(X[r].toarray().astype(np.float64), lib[r].astype(np.float64)) for r in rows0]
    m = np.array([s_[0] for s_ in st]); ps = np.array([s_[1] for s_ in st]); tau2, mu, _ = paule_mandel(m, ps)
    th = mu[None, :] + tau2[None, :] / (tau2[None, :] + ps) * (m - mu[None, :]); f = np.clip(mu[None, :] / np.maximum(th, 1e-30), 0.1, 10.0)
    idx = np.empty(X.shape[0], dtype=np.int64)
    for i, r in enumerate(rows0): idx[r] = i
    X = sp.csr_matrix(X.multiply(f[idx])); seed = 41_000_000 + 10 * a.k + si; lab = np.random.default_rng(seed).permutation(don)
else:
    seed = 43_000_000 + si; lab = don


def pear(Mx):
    Mc = Mx - Mx.mean(0); cov = Mc.T @ Mc / (len(Mx) - 1); sd = np.sqrt(np.diag(cov))
    return np.clip((cov / np.outer(sd, sd))[iu], -0.999999, 0.999999)


def welch(q1, q0):
    z = (q1.mean(0) - q0.mean(0)) / np.sqrt(q1.var(0, ddof=1) / len(q1) + q0.var(0, ddof=1) / len(q0))
    return 2 * norm.sf(np.abs(z)), np.sign(z)


def run(Xm, libm, lab_, seed_):
    ud = np.unique(lab_); grp = np.array([adlab[u] == "1.0" for u in ud])
    cells = [(Xm[np.where(lab_ == u)[0]].toarray().astype(np.float64), libm[np.where(lab_ == u)[0]].astype(np.float64)) for u in ud]
    f_, hs = stats_from_cells(cells, seed_)
    L = np.log1p(1e6 * f_["m"]); n1, n0 = grp.sum(), (~grp).sum(); out = {}
    zpb = (np.arctanh(pear(L[grp])) - np.arctanh(pear(L[~grp]))) / np.sqrt(1 / (n1 - 3) + 1 / (n0 - 3))
    out["PB"], out["sPB"] = 2 * norm.sf(np.abs(zpb)), np.sign(zpb)

    def qg(h, g):
        RA, _ = design_R(h["mA"][g], h["psiA"][g], h["kappaA"][g], h["SA"][g]); RB, _ = design_R(h["mB"][g], h["psiB"][g], h["kappaB"][g], h["SB"][g])
        return RA[:, A_] * RB[:, B_]
    out["XDSg"], out["sXDSg"] = welch(qg(hs[0], grp), qg(hs[0], ~grp))
    q1 = sum(qg(h, grp) for h in hs) / len(hs); q0 = sum(qg(h, ~grp) for h in hs) / len(hs)
    out["XDSMg"], out["sXDSMg"] = welch(q1, q0)
    # REV2 Amendment 1 (code audit, before any D computation): unstandardized group-centred cross-split product,
    # whose mean is the between-donor covariance whatever the donor's cell number
    def qu(h, g):
        mA = h["mA"][g]; mB = h["mB"][g]; return (mA - mA.mean(0))[:, A_] * (mB - mB.mean(0))[:, B_]
    u1 = sum(qu(h, grp) for h in hs) / len(hs); u0 = sum(qu(h, ~grp) for h in hs) / len(hs)
    out["XDSMgu"], out["sXDSMgu"] = welch(u1, u0)
    return out, int(n1), int(n0)


t0 = time.time()
if a.k == -2:
    ud = np.unique(don); order = np.random.default_rng(42_000_000).permutation(len(ud)); hv = []
    for j, rows in enumerate((np.sort(order[: len(ud) // 2]), np.sort(order[len(ud) // 2:]))):
        sel = set(ud[rows]); mask = np.array([d_ in sel for d_ in don])
        o, _, _ = run(X[np.where(mask)[0]], lib[mask], don[mask], seed + 1 + j); hv.append(o)
    res = {}
    for m_ in ("PB", "XDSg", "XDSMg", "XDSMgu"):
        c0 = bh(hv[0][m_]); rep = c0 & (hv[1][m_] < 0.05) & (hv[0]["s" + m_] == hv[1]["s" + m_])
        res[m_] = dict(calls_h0=int(c0.sum()), replicated=int(rep.sum()))
    json.dump(res, open(os.path.join(HERE, "out", f"D3_{a.sub}.json"), "w")); print(a.sub, "halves", res, flush=True); sys.exit()
o, n1, n0 = run(X, lib, lab, seed)
tag = f"k{a.k}" if a.k >= 0 else "real"
res = {m_: dict(V=int(bh(o[m_]).sum()), t05=float((o[m_] < 0.05).mean())) for m_ in ("PB", "XDSg", "XDSMg", "XDSMgu")}; res["n1"] = n1; res["n0"] = n0
json.dump(res, open(os.path.join(HERE, "out", f"D1_{a.sub}_{tag}.json"), "w")); print(a.sub, tag, res, round(time.time() - t0), "s", flush=True)
