"""Plan 6e (post hoc, rev5/PLAN6E.md). python rev5_e.py --mode C|D|A|R --sub P1 [--k0 --k1] -> $REV5_OUT/e_<mode>_<sub>[_k<k>].json"""
import os, sys, json, argparse, time, importlib.util, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(); ap.add_argument("--mode", required=True); ap.add_argument("--sub", default="P1")
ap.add_argument("--k0", type=int, default=0); ap.add_argument("--k1", type=int, default=20); a = ap.parse_args()
sys.argv = [sys.argv[0], "--mode", "none", "--sub", a.sub]
spec = importlib.util.spec_from_file_location("r5", os.path.join(HERE, "rev5.py")); r5 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r5)
from methods_v2 import bh  # noqa: E402
from revm import stats_from_cells  # noqa: E402
from dstd import design_R  # noqa: E402
X, lib, don, sym, G, iu, A_, B_, si, elig = r5.X, r5.lib, r5.don, r5.sym, r5.G, r5.iu, r5.A_, r5.B_, r5.si, r5.elig; OD = r5.OD


def hc_splits(hs, nd, W):
    nA = np.floor(nd / 2); nB = nd - nA; out = []
    for h in hs:
        pr = []
        for L_, nh in (("A", nA), ("B", nB)):
            m = h["m" + L_]; _, v = design_R(m, h["psi" + L_], h["kappa" + L_], h["S" + L_])
            Wh = np.column_stack([np.ones(len(nd)), np.log(h["kappa" + L_] * h["S" + L_]), np.log(nh)]); pr.append(r5.resid(r5.resid(m, W) / np.sqrt(v), Wh))
        out.append(tuple(pr))
    return out


def run(cl, hs, W, sel=None):
    nd = np.array([len(s) for _, s in cl], float)
    if sel is not None: nd = nd[sel]; hs = [{k: v[sel] for k, v in h.items()} for h in hs]; W = W[sel]
    return dict(XDSM=r5.xdsm_W(r5.Rsplits(hs), iu)[0], XDSMhc=r5.xdsm_W(hc_splits(hs, nd, W), iu)[0])


def summ(pv, extra=None):
    o = {}
    for m_, v in pv.items():
        v = np.asarray(v, float); c = bh(v); o[m_] = dict(calls=int(c.sum()), t05=float((v < 0.05).mean()))
        if extra is not None: o[m_].update({nm: int((c & msk).sum()) for nm, msk in extra.items()}); o[m_].update({"t05_" + nm: float((v[msk] < 0.05).mean()) for nm, msk in extra.items()})
    return o


if a.mode in ("C", "D", "A"):
    Xs = r5.centre(X, lib, don)
    if a.mode == "C": g200 = np.random.default_rng(76_000_000 + si).choice(elig, 200, replace=False); pm = np.zeros(G, bool); pm[g200] = True; planted = pm[A_] & pm[B_]
    if a.mode == "D": dbin = np.digitize(np.log(lib), np.quantile(np.log(lib), np.linspace(0, 1, 11)[1:-1]))
    for k in range(a.k0, a.k1):
        fn = os.path.join(OD, f"e_{a.mode}_{a.sub}_k{k}.json")
        if os.path.exists(fn): continue
        t0 = time.time(); seed = {"C": 75_000_000, "D": 73_000_000, "A": 61_000_000}[a.mode] + 10 * k + si; rng = np.random.default_rng(seed)
        if a.mode == "D":
            lab = don.copy()
            for b in range(10):
                ix = np.where(dbin == b)[0]; lab[ix] = rng.permutation(don[ix])
        else:
            lab = rng.permutation(don)
        pu, cl = r5.cells_of(Xs, lib, lab); W = r5.covW(pu)
        if a.mode == "C":
            case = W[:, 5] == 1; cl = [((Y * np.where(pm, 1.2, 1.0)[None, :]) if case[i] else Y, s) for i, (Y, s) in enumerate(cl)]
        fs, hs = stats_from_cells(cl, seed)
        res = {"all": summ(run(cl, hs, W), dict(planted=planted) if a.mode == "C" else None)}
        if a.mode == "A":
            order = np.random.default_rng(seed + 1).permutation(len(pu)); res["200"] = summ(run(cl, hs, W, np.sort(order[:200])))
        res["D"] = int(len(pu)); res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w")); print("e", a.mode, a.sub, k, {d: {m: v["calls"] for m, v in r.items()} for d, r in res.items() if isinstance(r, dict)}, res["sec"], "s", flush=True)

elif a.mode == "R":
    t0 = time.time(); pu, cl = r5.cells_of(X, lib, don); W = r5.covW(pu); fs, hs = stats_from_cells(cl, 44_000_001); pv = run(cl, hs, W)
    cX = bh(pv["XDSM"]); cH = bh(pv["XDSMhc"]); g = list(sym); YG = r5.YG
    pairs = {"XIST_Y": np.array([i for i in range(len(A_)) if (g[A_[i]] == "XIST" and g[B_[i]] in YG) or (g[B_[i]] == "XIST" and g[A_[i]] in YG)], int),
             "Y_Y": np.array([i for i in range(len(A_)) if g[A_[i]] in YG and g[B_[i]] in YG], int)}
    res = dict(D=int(len(pu)), XDSM=int(cX.sum()), XDSMhc=int(cH.sum()), keep=float((cX & cH).sum() / max(cX.sum(), 1)), sex={nm: dict(called=int(cH[ix].sum()), n=int(len(ix)), median_p=float(np.median(pv["XDSMhc"][ix]))) for nm, ix in pairs.items()})
    res["sec"] = round(time.time() - t0); json.dump(res, open(os.path.join(OD, f"e_R_{a.sub}.json"), "w")); print("e R", a.sub, res, flush=True)
