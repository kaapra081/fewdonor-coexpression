"""Plan 6b (post hoc, rev5/PLAN6B.md). python rev5_post.py --mode C|D|L|R --sub P1 [--k0 0 --k1 20]  -> $REV5_OUT/post_<mode>_<sub>[_k<k>].json"""
import os, sys, json, argparse, time, importlib.util, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(); ap.add_argument("--mode", required=True); ap.add_argument("--sub", default="P1")
ap.add_argument("--k0", type=int, default=0); ap.add_argument("--k1", type=int, default=20); a = ap.parse_args()
sys.argv = [sys.argv[0], "--mode", "none", "--sub", a.sub]          # rev5.py loads data and helpers, runs nothing for mode none
spec = importlib.util.spec_from_file_location("r5", os.path.join(HERE, "rev5.py")); r5 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r5)
from methods_v2 import bh  # noqa: E402
from dstd import _ms, design_R  # noqa: E402
from revm import stats_from_cells  # noqa: E402
X, lib, don, sym, G, iu, A_, B_, si, pool, elig = r5.X, r5.lib, r5.don, r5.sym, r5.G, r5.iu, r5.A_, r5.B_, r5.si, r5.pool, r5.elig
OD = r5.OD


def ca_splits(hs, W):
    """Adjust first (OLS of each half's ratio means on W), standardize second (design-only v*)."""
    out = []
    for h in hs:
        pr = []
        for L_ in "AB":
            m = h["m" + L_]; _, v = design_R(m, h["psi" + L_], h["kappa" + L_], h["S" + L_]); pr.append(r5.resid(m, W) / np.sqrt(v))
        out.append(tuple(pr))
    return out


def rw_pooled():
    num = np.zeros((G, G)); vx = np.zeros(G)
    for u in np.unique(don):
        r = np.where(don == u)[0]
        if len(r) < 2: continue
        x = X[r].toarray() / lib[r][:, None]; xc = x - x.mean(0); num += xc.T @ xc; vx += (xc ** 2).sum(0)
    return (num / np.sqrt(np.maximum(np.outer(vx, vx), 1e-300)))[iu]


def source_split_hs(cl_src, rng, S=5):
    """cl_src: list of (Y, s, src). Returns kept index and S half-stat dicts with all cells of a source donor in one half."""
    keep, parts = [], []
    for i, (Y, s, src) in enumerate(cl_src):
        us = np.unique(src); n = len(s); sp_ = []
        if len(us) < 2: continue
        for _ in range(S):
            for _t in range(20):
                A = []; cum = 0
                for x in rng.permutation(us):
                    if cum >= n / 2: break
                    A.append(x); cum += int((src == x).sum())
                inA = np.isin(src, A)
                if inA.sum() >= 2 and (~inA).sum() >= 2: break
            else:
                sp_ = None; break
            sp_.append(inA)
        if sp_ is None: continue
        keep.append(i); parts.append(sp_)
    hs = []
    for s_ in range(S):
        h = {k: [] for k in ("mA", "psiA", "kappaA", "SA", "mB", "psiB", "kappaB", "SB")}
        for i, sp_ in zip(keep, parts):
            Y, s, _ = cl_src[i]; inA = sp_[s_]
            for L_, msk in (("A", inA), ("B", ~inA)):
                m2, ps2, k2, S2, _, _ = _ms(Y[msk], s[msk]); h["m" + L_].append(m2); h["psi" + L_].append(ps2); h["kappa" + L_].append(k2); h["S" + L_].append(S2)
        hs.append({k: np.array(v) for k, v in h.items()})
    return np.array(keep), hs


def rep(pv, extra=None):
    o = {}
    for m_, v in pv.items():
        v = np.asarray(v, float); c = bh(v); o[m_] = dict(calls=int(c.sum()), t05=float((v < 0.05).mean()))
        if extra is not None: o[m_].update({nm: int((c & msk).sum()) for nm, msk in extra.items()})
        if extra is not None: o[m_].update({"t05_" + nm: float((v[msk] < 0.05).mean()) for nm, msk in extra.items()})
    return o


if a.mode in ("C", "D", "L"):
    Xs = r5.centre(X, lib, don, exact=(a.mode == "L"))
    if a.mode == "C":
        g200 = np.random.default_rng(76_000_000 + si).choice(elig, 200, replace=False); pm = np.zeros(G, bool); pm[g200] = True; planted = pm[A_] & pm[B_]
    if a.mode == "D":
        dbin = np.digitize(np.log(lib), np.quantile(np.log(lib), np.linspace(0, 1, 11)[1:-1]))
    if a.mode == "L": rw = rw_pooled()
    for k in range(a.k0, a.k1):
        fn = os.path.join(OD, f"post_{a.mode}_{a.sub}_k{k}.json")
        if os.path.exists(fn): continue
        t0 = time.time(); seed = {"D": 73_000_000, "L": 74_000_000, "C": 75_000_000}[a.mode] + 10 * k + si; rng = np.random.default_rng(seed); res = {}
        if a.mode == "D":
            lab = don.copy()
            for b in range(10):
                ix = np.where(dbin == b)[0]; lab[ix] = rng.permutation(don[ix])
        elif a.mode == "L":
            lab = don.copy()
            for p_ in np.unique(pool):
                ix = np.where(pool == p_)[0]; lab[ix] = rng.permutation(don[ix])
        else:
            lab = rng.permutation(don)
        pu, cl = r5.cells_of(Xs, lib, lab)
        if a.mode == "C":
            W = r5.covW(pu); case = W[:, 5] == 1; cl = [((Y * np.where(pm, 1.2, 1.0)[None, :]) if case[i] else Y, s) for i, (Y, s) in enumerate(cl)]
        fs, hs = stats_from_cells(cl, seed); Rs = r5.Rsplits(hs); pX, mqX = r5.xdsm_W(Rs, iu); pv = dict(XDSM=pX)
        if a.mode == "C":
            pv["XDSMc"] = r5.xdsm_W(Rs, iu, W)[0]; pv["XDSMca_ad"] = r5.xdsm_W(ca_splits(hs, W[:, [0, 5, 6]]), iu)[0]; pv["XDSMca"] = r5.xdsm_W(ca_splits(hs, W), iu)[0]
            res = rep(pv, dict(planted=planted, other=~planted))
        elif a.mode == "D":
            nd = np.array([len(s) for _, s in cl], float); Wd = np.column_stack([np.ones(len(cl)), np.log(r5.rhat(cl)), np.log(nd)])
            pv["XDSMd"] = r5.xdsm_W(Rs, iu, Wd)[0]; pv["XDSMda"] = r5.xdsm_W(ca_splits(hs, Wd), iu)[0]; res = rep(pv)
        else:
            cX = bh(pX); sg = np.sign(mqX) == -np.sign(rw); res["XDSM_sign"] = dict(calls=int(cX.sum()), opposite_to_rw=int((cX & sg).sum()), frac_opposite=float((cX & sg).sum() / max(cX.sum(), 1)), null_frac_opposite_all_pairs=float(sg.mean()))
            cl_src = [(Y, s, don[np.where(lab == u)[0]]) for (Y, s), u in zip(cl, pu)]
            kp, hsS = source_split_hs(cl_src, np.random.default_rng(seed + 11)); pv["XDSMsrc"] = r5.xdsm_W(r5.Rsplits(hsS), iu)[0]; res.update(rep(pv)); res["D_src"] = int(len(kp))
        res["D"] = int(len(pu)); res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w"))
        print("post", a.mode, a.sub, k, {m_: res[m_]["calls"] for m_ in pv}, res["sec"], "s", flush=True)

elif a.mode == "R":
    t0 = time.time(); pu, cl = r5.cells_of(X, lib, don); fs, hs = stats_from_cells(cl, 44_000_001); Rs = r5.Rsplits(hs); W = r5.covW(pu); res = dict(D=int(len(pu)))
    pX, _ = r5.xdsm_W(Rs, iu); cX = bh(pX); pA = r5.xdsm_W(ca_splits(hs, W), iu)[0]; cA = bh(pA); pS = r5.xdsm_W(ca_splits(hs, W[:, :2]), iu)[0]
    res["XDSMca"] = dict(calls=int(cA.sum()), keep_of_XDSM=float((cA & cX).sum() / max(cX.sum(), 1)), XDSM=int(cX.sum()))
    g = list(sym); YG = r5.YG; yi = [g.index(y) for y in YG if y in g]; xi = g.index("XIST")
    pairs = {"XIST_Y": np.array([i for i in range(len(A_)) if (g[A_[i]] == "XIST" and g[B_[i]] in YG) or (g[B_[i]] == "XIST" and g[A_[i]] in YG)], int),
             "Y_Y": np.array([i for i in range(len(A_)) if g[A_[i]] in YG and g[B_[i]] in YG], int)}
    def sexrep(p):
        c = bh(p); o = {nm: dict(called=int(c[ix].sum()), n=int(len(ix)), median_p=float(np.median(p[ix]))) for nm, ix in pairs.items()}; o["calls_total"] = int(c.sum()); return o
    res["sex_XDSMca_all"] = sexrep(pA); res["sex_XDSMca_sexonly"] = sexrep(pS)
    # H3: bug-fixed R-4(c)
    male = W[:, 1] == 1; Ysum = np.array([Y[:, yi].sum() for Y, _ in cl]); Xx = np.array([Y[:, xi].sum() for Y, _ in cl]); mism = (np.log((Ysum + 1) / (Xx + 1)) > 0) != male; keep = ~mism
    cl3 = []; nrem = 0
    for i in np.where(keep)[0]:
        Y, s = cl[i]; ys = Y[:, yi].sum(1); xs = Y[:, xi]; bad = ((xs >= 2) & (ys == 0)) if male[i] else ((ys >= 2) & (xs == 0))
        if (~bad).sum() >= 4: cl3.append((Y[~bad], s[~bad])); nrem += int(bad.sum())
        else: cl3.append((Y, s))
    _, hs3 = stats_from_cells(cl3, 44_000_001); Rs3 = r5.Rsplits(hs3); W2 = W[keep]
    res["clean_fixed"] = dict(nuclei_removed=nrem, XDSMc=sexrep(r5.xdsm_W(Rs3, iu, W2)[0]), sexonly=sexrep(r5.xdsm_W(Rs3, iu, W2[:, :2])[0]), XDSMca=sexrep(r5.xdsm_W(ca_splits(hs3, W2), iu)[0]), XDSMca_sexonly=sexrep(r5.xdsm_W(ca_splits(hs3, W2[:, :2]), iu)[0]))
    res["sec"] = round(time.time() - t0); json.dump(res, open(os.path.join(OD, f"post_R_{a.sub}.json"), "w"), indent=1); print("post R", a.sub, json.dumps(res)[:1500], flush=True)
