"""REV3_PREREG A2: tail calibration of PB, XDS and XDS-M (S = 1, 5, 10) by gene-expression stratum. python rev3_A2.py --ds B1 --null 0 --k0 0 --k1 20"""
import os, sys, json, argparse, time, numpy as np, scipy.sparse as sp, torch
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("dev", "dstd", "pop", "rev"): sys.path.insert(0, os.path.join(HERE, "..", p))
from methods_v2 import bh, paule_mandel  # noqa: E402
from dstd import _ms, design_R  # noqa: E402
from revm import stats_from_cells, methods as rmethods  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--ds", required=True); ap.add_argument("--null", type=int, required=True)
ap.add_argument("--k0", type=int, default=0); ap.add_argument("--k1", type=int, default=20); ap.add_argument("--dev", default="cuda:0"); a = ap.parse_args()
dev = torch.device(a.dev if torch.cuda.is_available() else "cpu"); torch.set_num_threads(1)
dsi = {"B1": 1, "B2": 2, "P3": 3, "B1x": 1}[a.ds]; BASE = 71_000_000 if a.ds.endswith("x") else 70_000_000; OD = os.path.join(HERE, "out"); os.makedirs(OD, exist_ok=True)
X = sp.load_npz(os.path.join(HERE, "data", f"a2_{a.ds}_X.npz")).tocsr(); M = np.load(os.path.join(HERE, "data", f"a2_{a.ds}_meta.npz"), allow_pickle=True)
lib = M["lib"].astype(float); don = M["donor"].astype(str); umi = M["umi"]; G = X.shape[1]
edges = [(1, np.inf), (0.3, 1), (0.1, 0.3), (0.03, 0.1), (-1, 0.03)]
gs = np.full(G, 4)
for i, (lo, hi) in enumerate(edges): gs[(umi > lo) & (umi <= hi)] = i
iu = np.triu_indices(G, 1); A_, B_ = iu; ps_ = np.maximum(gs[A_], gs[B_])
THR = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
if a.null == 1:   # donor-centred (REV_PREREG R6): remove each real donor's shrunken effect before pooling
    ud0 = np.unique(don); rows0 = [np.where(don == u)[0] for u in ud0]
    st = [_ms(X[r].toarray().astype(np.float64), lib[r]) for r in rows0]
    m0 = np.array([s_[0] for s_ in st]); p0 = np.array([s_[1] for s_ in st]); tau2, mu, _ = paule_mandel(m0, p0)
    th = mu[None, :] + tau2[None, :] / (tau2[None, :] + p0) * (m0 - mu[None, :]); f = np.clip(mu[None, :] / np.maximum(th, 1e-30), 0.1, 10.0)
    idx = np.empty(X.shape[0], dtype=np.int64)
    for i, r in enumerate(rows0): idx[r] = i
    X = sp.csr_matrix(X.multiply(f[idx])); del st, m0, p0, th, f
T = lambda x: torch.as_tensor(np.ascontiguousarray(x), dtype=torch.float32, device=dev)
def ut(Mx): return Mx[iu[0], iu[1]] if isinstance(Mx, np.ndarray) else Mx.cpu().numpy().astype(np.float64)[iu]
def gpu_stats(L, Rs):
    """L: log1p(1e6 m) [D x G]; Rs: list of (RA, RB) per split. Returns dict of upper-triangle statistics (float64)."""
    D = L.shape[0]; out = {}
    Lt = T(L); Lc = Lt - Lt.mean(0); C = Lc.T @ Lc / (D - 1); sd = torch.sqrt(torch.diag(C)); r = ut(C / torch.outer(sd, sd))
    out["PB"] = 2 * norm.sf(np.abs(np.arctanh(np.clip(r, -0.999999, 0.999999)) * np.sqrt(D - 3))); out["sPB"] = np.sign(r)
    RA, RB = T(Rs[0][0]), T(Rs[0][1]); Tm = ut(RA.T @ RB)
    mu = ut(torch.outer(RA.sum(0), RB.sum(0)) / D); var = ut(torch.outer(((RA - RA.mean(0)) ** 2).sum(0), ((RB - RB.mean(0)) ** 2).sum(0)) / (D - 1))
    sdv = np.sqrt(np.maximum(var, 1e-300)); out["XDS"] = norm.sf((np.abs(Tm) - mu) / sdv) + norm.cdf((-np.abs(Tm) - mu) / sdv); out["sXDS"] = np.sign(Tm)
    Sx = len(Rs); RAs = [T(x[0]) for x in Rs]; RBs = [T(x[1]) for x in Rs]
    for S in (1, 5, 10):
        if S > Sx: continue
        Q1 = sum(RAs[s].T @ RBs[s] for s in range(S)) / (S * D)
        Q2 = torch.zeros((G, G), dtype=torch.float32, device=dev)
        for s in range(S):
            for s2 in range(s, S):
                Q2 += (1.0 if s == s2 else 2.0) * ((RAs[s] * RAs[s2]).T @ (RBs[s] * RBs[s2]))
        Q2 /= (S * S * D); q1 = ut(Q1); q2 = ut(Q2); var = np.maximum((q2 - q1 ** 2) * D / (D - 1), 1e-300)
        t = q1 / np.sqrt(var / D); out[f"XDSM{S}"] = 2 * norm.sf(np.abs(t)); out[f"sXDSM{S}"] = np.sign(t)
    return out
def summarize(o):
    res = {}
    for m in ("PB", "XDS", "XDSM1", "XDSM5", "XDSM10"):
        p = o[m]; sg = o["s" + m]; rj = bh(p); r_ = dict(V=int(rj.sum()), strata={})
        for si in range(5):
            w = ps_ == si; pw = p[w]; r_["strata"][str(si)] = dict(n=int(w.sum()), V=int(bh(pw).sum()) if w.any() else 0, Vall=int(rj[w].sum()),
                exc=[int((pw <= t_).sum()) for t_ in THR], pos=[int(((pw <= t_) & (sg[w] > 0)).sum()) for t_ in THR], t05=float((pw < 0.05).mean()) if w.any() else None)
        res[m] = r_
    return res
for k in range(a.k0, a.k1):
    fn = os.path.join(OD, f"A2_{a.ds}_n{a.null}_k{k}.json")
    if os.path.exists(fn): continue
    t0 = time.time(); seed = BASE + 1_000 * a.null + 10 * k + dsi; rng = np.random.default_rng(seed); lab = rng.permutation(don)
    pu = np.unique(lab); prow = [np.where(lab == u)[0] for u in pu]; cells = [(X[r].toarray().astype(np.float64), lib[r]) for r in prow]
    full, halves = stats_from_cells(cells, seed, S=10); del cells; Dall = len(pu); order = rng.permutation(Dall); res = {"D_all": Dall, "check": {}}
    for D in (200, Dall):
        sel = np.sort(order[:D]); fsub = {kk: v[sel] for kk, v in full.items()}; hsub = [{kk: v[sel] for kk, v in h.items()} for h in halves]
        Rs = [(design_R(h["mA"], h["psiA"], h["kappaA"], h["SA"])[0], design_R(h["mB"], h["psiB"], h["kappaB"], h["SB"])[0]) for h in hsub]
        o = gpu_stats(np.log1p(1e6 * fsub["m"]), Rs)
        # float64 CPU check on 2,000 random pairs (frozen validity rule)
        cr = np.random.default_rng(seed + D); pick = np.sort(cr.choice(len(A_), 2000, replace=False)); iu2 = (A_[pick], B_[pick]); md = 0.0
        for S, key in ((1, "XDSM1"), (5, "XDSM5"), (10, "XDSM10")):
            oc = rmethods(fsub, hsub[:S], iu2)
            for mm, ref in ((key, oc["XDSM"]),) + (((("XDS", oc["XDS"]), ("PB", oc["PB"])) if S == 5 else ())):
                g_ = o[mm][pick]; ok = (ref > 1e-8) & (g_ > 1e-8); md = max(md, float(np.max(np.abs(np.log10(g_[ok]) - np.log10(ref[ok])))) if ok.any() else 0.0)
        res["check"][str(D)] = md; res["all" if D == Dall else str(D)] = summarize(o)
    res["valid"] = all(v <= 0.05 for v in res["check"].values()); res["device"] = str(dev)
    json.dump(res, open(fn, "w")); print(a.ds, a.null, k, "D", Dall, "check", res["check"], "V", {m: res["all"][m]["V"] for m in res["all"]}, round(time.time() - t0), "s", flush=True)
