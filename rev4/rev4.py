"""REV4_PREREG R1-R4. python rev4.py --mode R1 --sub P1 --k0 0 --k1 20
R1: small-D donor-centred null, B = 2e6 permutations   R2: within-brain-bank (Source) donor-centred null, full D
R3: PC-adjusted pseudobulk, donor-centred null (B2 seeds)   R3o: OneK1K pooled null (--sub B1|B2|B3)   R4: real donors (--k0 ignored)"""
import os, sys, json, argparse, time, numpy as np, scipy.sparse as sp
from scipy.stats import norm, spearmanr
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("dev", "dstd", "pop", "rev"): sys.path.insert(0, os.path.join(HERE, "..", p))
from methods_v2 import bh, paule_mandel, pearson_fisher  # noqa: E402
from dstd import _ms, donor_level_cells, methods_d  # noqa: E402
from revm import stats_from_cells, methods as rmethods  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--mode", required=True); ap.add_argument("--sub", required=True)
ap.add_argument("--k0", type=int, default=0); ap.add_argument("--k1", type=int, default=20); a = ap.parse_args()
OD = os.environ.get("REV4_OUT", os.path.join(HERE, "out")); os.makedirs(OD, exist_ok=True)


def load_pop(sub):
    P = os.path.join(HERE, "..", "pop"); X = sp.load_npz(os.path.join(P, f"pop_{sub}_X.npz")).tocsr()
    M = np.load(os.path.join(P, f"pop_{sub}_meta.npz"), allow_pickle=True); return X, M


def centre(X, lib, don, lo=0.1, hi=10.0):   # identical to REV2 A / REV3 A3 / B2 donor-centring
    ud = np.unique(don); rows = [np.where(don == u)[0] for u in ud]
    st = [_ms(X[r].toarray().astype(np.float64), lib[r].astype(np.float64)) for r in rows]
    m = np.array([s_[0] for s_ in st]); ps = np.array([s_[1] for s_ in st]); tau2, mu, _ = paule_mandel(m, ps)
    th = mu[None, :] + tau2[None, :] / (tau2[None, :] + ps) * (m - mu[None, :]); f = np.clip(mu[None, :] / np.maximum(th, 1e-30), lo, hi)
    idx = np.empty(X.shape[0], dtype=np.int64)
    for i, r in enumerate(rows): idx[r] = i
    return sp.csr_matrix(X.multiply(f[idx]))


def summ(pv):
    out = {}
    for mm, v in pv.items():
        v = np.asarray(v, float); out[mm] = dict(V=int(bh(v).sum()), t05=float((v < 0.05).mean()), minp=float(np.min(v)), nfloor=int((v <= 1e-5).sum()))
    return out


def fz(r, neff): return 2 * norm.sf(np.abs(np.arctanh(np.clip(r, -0.999999, 0.999999)) * np.sqrt(max(neff, 1.0))))


def corr_iu(Mx, iu):
    Mc = Mx - Mx.mean(0); C = Mc.T @ Mc; sd = np.sqrt(np.maximum(np.diag(C), 1e-300)); return (C / np.outer(sd, sd))[iu]


def pc_adj(L, k, iu):
    Lc = L - L.mean(0); U, s, Vt = np.linalg.svd(Lc, full_matrices=False); R = Lc - (U[:, :k] * s[:k]) @ Vt[:k]
    return fz(corr_iu(R, iu), L.shape[0] - 3 - k)


def pb_family(mm, nd, iu, with_res=False, dp=None):
    L = np.log1p(1e6 * mm); D = len(L); pv = dict(PB=pearson_fisher(L, iu)["p"])
    if with_res:
        Z = np.column_stack([np.ones(D), np.log(nd), np.log(dp)]); Rres = L - Z @ np.linalg.lstsq(Z, L, rcond=None)[0]; pv["PBres"] = fz(corr_iu(Rres, iu), D - 5)
    w = nd / nd.sum(); Lw = L - (w[:, None] * L).sum(0); Cw = (Lw * w[:, None]).T @ Lw; sw = np.sqrt(np.diag(Cw))
    neff = nd.sum() ** 2 / (nd ** 2).sum(); pv["PBw"] = fz((Cw / np.outer(sw, sw))[iu], neff - 3)
    pv["PC5"] = pc_adj(L, 5, iu); pv["PC10"] = pc_adj(L, 10, iu); return pv, float(neff)


if a.mode in ("R1", "R2"):
    import torch; torch.set_num_threads(1)
    from pop_core import run_methods, subset  # noqa: E402
    si = ["P1", "P2", "P3"].index(a.sub); X, M = load_pop(a.sub); lib = M["lib"]; don = M["donor"].astype(str); iu = np.triu_indices(X.shape[1], 1)
    Xs = centre(X, lib, don)
    if a.mode == "R1": base, DS, B = 63_000_000, (25, 50, 100), int(os.environ.get("REV4_B", 2_000_000))   # REV4_B only for smoke tests
    else:
        base, DS, B = 64_000_000, (None,), 100_000
        src = json.load(open(os.path.join(HERE, "donor_source.json"))); bank = np.array([src[d_] for d_ in don])
    for k in range(a.k0, a.k1):
        fn = os.path.join(OD, f"{a.mode}_{a.sub}_k{k}.json")
        if os.path.exists(fn): continue
        t0 = time.time(); seed = base + 10 * k + si; rng = np.random.default_rng(seed)
        if a.mode == "R1": lab = rng.permutation(don)
        else:
            lab = don.copy()
            for bk in ("H", "M", "R"):
                ix = np.where(bank == bk)[0]; lab[ix] = rng.permutation(don[ix])
        pu = np.unique(lab); prow = [np.where(lab == u)[0] for u in pu]; cells = [(Xs[r].toarray(), lib[r].astype(np.float64)) for r in prow]
        d = donor_level_cells(cells, seed + 7); d["donors"] = pu; fs, hs = stats_from_cells(cells, seed); Dall = len(pu)
        order = rng.permutation(Dall); res = {}
        for D in DS:
            D = D or Dall; sel = np.sort(order[:D]); sd = subset(d, sel)
            o1 = run_methods(sd, iu, B=B, seed=seed + D, skip_x=True); o2 = methods_d(sd, iu, B=B, seed=seed + D + 3)
            fsub = {kk: v[sel] for kk, v in fs.items()}; hsub = [{kk: v[sel] for kk, v in h.items()} for h in hs]
            o3 = rmethods(fsub, hsub, iu)
            pv = dict(PB=o1["PB"], LCP=o1["LCP"], LCPd=o2["LCPd"], XDSperm=o2["LCPXd"], XDS=o3["XDS"], XDSM=o3["XDSM"], KX=o3["KX"])
            if a.mode == "R2": pv.pop("XDSperm")   # D > 100: LCPXd uses the normal reference, identical in kind to XDS
            res["all" if D == Dall else str(D)] = summ(pv)
        res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w"))
        print(a.mode, a.sub, k, {D: {mm: r[mm]["V"] for mm in r} for D, r in res.items() if D != "sec"}, res["sec"], "s", flush=True)

elif a.mode == "R3":
    si = ["P1", "P2", "P3"].index(a.sub); X, M = load_pop(a.sub); lib = M["lib"].astype(float); don = M["donor"].astype(str); iu = np.triu_indices(X.shape[1], 1)
    Xs = centre(X, lib, don)
    for k in range(a.k0, a.k1):
        fn = os.path.join(OD, f"R3_{a.sub}_k{k}.json")
        if os.path.exists(fn): continue
        seed = 61_000_000 + 10 * k + si; lab = np.random.default_rng(seed).permutation(don); pu = np.unique(lab)
        mm, nd, dp = [], [], []
        for u in pu:
            r = np.where(lab == u)[0]; mm.append(np.asarray(Xs[r].sum(0)).ravel() / lib[r].sum()); nd.append(len(r)); dp.append(lib[r].mean())
        pv, neff = pb_family(np.array(mm), np.array(nd, float), iu, with_res=True, dp=np.array(dp))
        res = summ(pv); res["D"] = len(pu); res["neff"] = neff; json.dump(res, open(fn, "w")); print("R3", a.sub, k, {m_: res[m_]["V"] for m_ in pv}, flush=True)

elif a.mode == "R3o":
    ti = ["B1", "B2", "B3"].index(a.sub); OK = os.path.join(HERE, "..", "onek")
    X = sp.load_npz(os.path.join(OK, f"onek_{a.sub}_X.npz")).tocsr(); M = np.load(os.path.join(OK, f"onek_{a.sub}_meta.npz"), allow_pickle=True)
    lib = M["lib"].astype(float); don = M["donor"].astype(str); iu = np.triu_indices(X.shape[1], 1)
    for k in range(a.k0, a.k1):
        fn = os.path.join(OD, f"R3o_{a.sub}_k{k}.json")
        if os.path.exists(fn): continue
        seed = 65_000_000 + 10 * k + ti; lab = np.random.default_rng(seed).permutation(don); pu = np.unique(lab)
        mm, nd = [], []
        for u in pu:
            r = np.where(lab == u)[0]; mm.append(np.asarray(X[r].sum(0)).ravel() / lib[r].sum()); nd.append(len(r))
        pv, neff = pb_family(np.array(mm), np.array(nd, float), iu)
        res = summ(pv); res["D"] = len(pu); res["neff"] = neff; json.dump(res, open(fn, "w")); print("R3o", a.sub, k, {m_: res[m_]["V"] for m_ in pv}, flush=True)

elif a.mode == "R4":
    X, M = load_pop(a.sub); lib = M["lib"].astype(float); don = M["donor"].astype(str); sym = M["sym"].astype(str); iu = np.triu_indices(X.shape[1], 1); A_, B_ = iu
    ud = np.unique(don); rows = [np.where(don == u)[0] for u in ud]; cl = [(X[r].toarray().astype(np.float64), lib[r]) for r in rows]
    f, h = stats_from_cells(cl, 44_000_001); o = rmethods(f, h, iu); D = len(ud)
    pb = pearson_fisher(np.log1p(1e6 * f["m"]), iu); cP = bh(pb["p"]); cX = bh(o["XDSM"])
    C = np.zeros((X.shape[1], X.shape[1]))
    for hh in h:
        mA = hh["mA"] - hh["mA"].mean(0); mB = hh["mB"] - hh["mB"].mean(0); C += mA.T @ mB / (D - 1)
    C = (C + C.T) / (2 * len(h)); dg = np.diag(C)
    with np.errstate(invalid="ignore", divide="ignore"):
        rho = np.where((dg[A_] > 0) & (dg[B_] > 0), C[iu] / np.sqrt(np.abs(dg[A_] * dg[B_])), np.nan)
    S = np.zeros(len(A_)); N = np.zeros(len(A_))
    for Y, s in cl:
        if len(s) < 10: continue
        Z = Y / s[:, None]; Zc = Z - Z.mean(0); sd = np.sqrt((Zc ** 2).sum(0)); ok = sd > 0
        R = (Zc.T @ Zc) / np.maximum(np.outer(sd, sd), 1e-300); r = R[iu]; v = ok[A_] & ok[B_]; S[v] += r[v]; N[v] += 1
    rw = np.where(N > 0, S / np.maximum(N, 1), np.nan); dlt = pb["est"] - rho
    def med(x, msk):
        x = x[msk]; x = x[np.isfinite(x)]; return float(np.median(x)) if len(x) else None
    cls = dict(PBonly=cP & ~cX, XDSMonly=cX & ~cP, shared=cP & cX, neither=~cP & ~cX)
    res = dict(D=int(D), n_pairs=int(len(A_)), calls=dict(PB=int(cP.sum()), XDSM=int(cX.sum())), donors_rw=int(sum(len(s) >= 10 for _, s in cl)))
    res["classes"] = {c: dict(n=int(m_.sum()), med_rw=med(rw, m_), med_rPB=med(pb["est"], m_), med_rho=med(rho, m_), med_rPB_minus_rho=med(dlt, m_),
                              sign_agree=float(np.nanmean((np.sign(pb["est"]) == np.sign(rho))[m_ & np.isfinite(rho)])) if (m_ & np.isfinite(rho)).any() else None)
                      for c, m_ in cls.items()}
    fin = np.isfinite(dlt) & np.isfinite(rw); res["spearman_delta_rw"] = float(spearmanr(dlt[fin], rw[fin]).correlation); res["n_finite"] = int(fin.sum())
    res["spearman_rPB_rho"] = float(spearmanr(pb["est"][np.isfinite(rho)], rho[np.isfinite(rho)]).correlation)
    def hubs(p, top=5000, nh=20):
        o_ = np.argsort(p)[:top]; deg = np.bincount(np.concatenate([A_[o_], B_[o_]]), minlength=X.shape[1]); return set(np.argsort(-deg, kind="stable")[:nh])
    hP, hX = hubs(pb["p"]), hubs(o["XDSM"]); res["hub_overlap_top20"] = len(hP & hX); res["hubs_PB"] = sorted(sym[list(hP)].tolist()); res["hubs_XDSM"] = sorted(sym[list(hX)].tolist())
    top = 5000; oP = set(np.argsort(pb["p"])[:top]); oX = set(np.argsort(o["XDSM"])[:top]); res["top5000_pair_overlap"] = len(oP & oX)
    json.dump(res, open(os.path.join(OD, f"R4_{a.sub}.json"), "w"), indent=1); print("R4", a.sub, json.dumps({k: v for k, v in res.items() if not k.startswith("hubs")}), flush=True)
