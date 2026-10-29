"""R6: python rev_centred.py --sub P1 --k 0"""
import os, sys, json, argparse, time, numpy as np, scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
for p in ("dev", "dstd"): sys.path.insert(0, os.path.join(HERE, "..", p))
from revm import stats_from_cells, methods  # noqa: E402
from methods_v2 import bh, paule_mandel  # noqa: E402
from dstd import _ms  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--sub", required=True); ap.add_argument("--k", type=int, required=True); a = ap.parse_args()
si = ["P1", "P2", "P3"].index(a.sub); P = os.path.join(HERE, "..", "pop")
X = sp.load_npz(os.path.join(P, f"pop_{a.sub}_X.npz")).tocsr(); M = np.load(os.path.join(P, f"pop_{a.sub}_meta.npz"), allow_pickle=True)
lib = M["lib"]; don = M["donor"].astype(str); iu = np.triu_indices(X.shape[1], 1); t0 = time.time()
ud = np.unique(don); rows = [np.where(don == u)[0] for u in ud]
st = [_ms(X[r].toarray().astype(np.float64), lib[r].astype(np.float64)) for r in rows]
m = np.array([s_[0] for s_ in st]); ps = np.array([s_[1] for s_ in st]); tau2, mu, _ = paule_mandel(m, ps)
th = mu[None, :] + tau2[None, :] / (tau2[None, :] + ps) * (m - mu[None, :]); f = np.clip(mu[None, :] / np.maximum(th, 1e-30), 0.1, 10.0)
Xc = X.astype(np.float64).tolil(copy=True) if False else None
dense_scale = np.empty(X.shape[0], dtype=np.int64); 
for i, r in enumerate(rows): dense_scale[r] = i
Xs = sp.csr_matrix(X.multiply(f[dense_scale]))   # row-wise gene-specific scaling by the donor factor
seed = 49_000_000 + 10 * a.k + si; lab = np.random.default_rng(seed).permutation(don)
pu = np.unique(lab); prow = [np.where(lab == u)[0] for u in pu]
fs, hs = stats_from_cells(((Xs[r].toarray(), lib[r].astype(np.float64)) for r in prow), seed)
o = methods(fs, hs, iu)
# within-pseudo-donor correlation (for top-1%) from the full-donor covariance, computed on the fly
G = X.shape[1]; Cw = np.zeros((G, G)); sig = np.zeros(G); n_ = 0
for r in prow:
    Y = Xs[r].toarray(); s = lib[r].astype(np.float64); x = Y / s[:, None]; xc = x - x.mean(0); C = xc.T @ xc / (len(r) - 1)
    Cw += C; sig += np.maximum(np.diag(C) - (x / s[:, None]).mean(0), 0); n_ += 1
Cw /= n_; sig /= n_; a_, b_ = iu
with np.errstate(invalid="ignore", divide="ignore"): rw = Cw[a_, b_] / np.sqrt(sig[a_] * sig[b_])
top = np.abs(rw) >= np.nanquantile(np.abs(rw), 0.99)
res = {mm: dict(V=int(bh(o[mm]).sum()), t05=float((o[mm] < 0.05).mean()), top1=float((o[mm][top] < 0.05).mean())) for mm in ("PB", "XDS", "XDSM", "KX")}
res["median_abs_rw"] = float(np.nanmedian(np.abs(rw))); json.dump(res, open(os.path.join(HERE, "out", f"R6_{a.sub}_k{a.k}.json"), "w"))
print(a.sub, a.k, res, round(time.time() - t0), "s", flush=True)
