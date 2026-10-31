"""REV2_PREREG B and C: python rev2_C.py --sub P1"""
import os, sys, json, argparse, numpy as np, scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("dev", "dstd", "rev"): sys.path.insert(0, os.path.join(HERE, "..", p))
from methods_v2 import bh  # noqa: E402
from dstd import design_R  # noqa: E402
from revm import stats_from_cells, methods as rmethods  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--sub", required=True); a = ap.parse_args(); si = ["P1", "P2", "P3"].index(a.sub)
P = os.path.join(HERE, "..", "pop"); X = sp.load_npz(os.path.join(P, f"pop_{a.sub}_X.npz")).tocsr(); M = np.load(os.path.join(P, f"pop_{a.sub}_meta.npz"), allow_pickle=True)
lib = M["lib"]; don = M["donor"].astype(str); sym = M["sym"].astype(str); iu = np.triu_indices(X.shape[1], 1); out = {}
def cells_for(lab):
    ud = np.unique(lab); return [(X[np.where(lab == u)[0]].toarray().astype(np.float64), lib[np.where(lab == u)[0]].astype(np.float64)) for u in ud]
# B: variance heterogeneity in the pooled null (POP seed k = 0)
seed = 41_000_000 + si; lab = np.random.default_rng(seed).permutation(don); cl = cells_for(lab)
f, h = stats_from_cells(cl, seed, S=1); RA, _ = design_R(h[0]["mA"], h[0]["psiA"], h[0]["kappaA"], h[0]["SA"]); n = np.array([len(s) for _, s in cl])
q = np.quantile(n, np.linspace(0, 1, 11)); dec = np.clip(np.searchsorted(q, n, side="right") - 1, 0, 9)
out["A5"] = dict(decile_median_n=[float(np.median(n[dec == i])) for i in range(10)], mean_R2=[float((RA[dec == i] ** 2).mean()) for i in range(10)])
# C: real donors
cl = cells_for(don); calls, pv = [], None
for s_ in (1, 2, 3):
    f, h = stats_from_cells(cl, 44_000_000 + s_); o = rmethods(f, h, iu); calls.append(bh(o["XDSM"]))
    if s_ == 1: pv = o
J = [float((calls[i] & calls[j]).sum() / max((calls[i] | calls[j]).sum(), 1)) for i, j in ((0, 1), (0, 2), (1, 2))]
out["seeds"] = dict(n_calls=[int(c.sum()) for c in calls], jaccard=J, frac_pairs=float(calls[0].mean()))
Y = {"UTY", "USP9Y", "ZFY", "DDX3Y", "KDM5D", "NLGN4Y", "TTTY14"}; sx = []
from methods_v2 import pearson_fisher
pbr = pearson_fisher(np.log1p(1e6 * f["m"]), iu)
for kk in range(len(iu[0])):
    g1, g2 = sym[iu[0][kk]], sym[iu[1][kk]]
    if ("XIST" in (g1, g2) and (g1 in Y or g2 in Y)) or (g1 in Y and g2 in Y):
        sx.append(dict(pair=f"{g1}-{g2}", xist="XIST" in (g1, g2), p_XDSM=float(pv["XDSM"][kk]), s_XDSM=int(pv["sign_XDSM"][kk]), p_PB=float(pbr["p"][kk]), r_PB=float(pbr["est"][kk])))
out["sexpairs"] = sx
ud = np.unique(don); order = np.random.default_rng(42_000_000).permutation(len(ud)); hv = []
for rows in (np.sort(order[: len(ud) // 2]), np.sort(order[len(ud) // 2:])):
    sel = set(ud[rows]); lab_h = np.array([d_ if d_ in sel else "" for d_ in don]); keep = lab_h != ""
    clh = [(X[np.where(don == u)[0]].toarray().astype(np.float64), lib[np.where(don == u)[0]].astype(np.float64)) for u in ud[rows]]
    f, h = stats_from_cells(clh, 42_000_000 + si); hv.append(rmethods(f, h, iu))
c0 = bh(hv[0]["XDSM"]); rep = c0 & (hv[1]["XDSM"] < 0.05) & (hv[1]["sign_XDSM"] == hv[0]["sign_XDSM"])
out["halves"] = dict(calls_h0=int(c0.sum()), replicated=int(rep.sum()))
json.dump(out, open(os.path.join(HERE, "out", f"C_{a.sub}.json"), "w"), indent=1); print(a.sub, {k: v for k, v in out.items() if k != "sexpairs"}, flush=True)
