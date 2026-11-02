"""Plan 6 (rev5/PLAN6.md). python rev5.py --mode S|D|L|C|R --sub P1 [--size full|200] [--k0 0 --k1 20]
Outputs rev5/out/<mode>_<sub>[_<size>]_k<k>.json (R: R_<sub>.json and R_<sub>.npz). SMOKE=1 offsets seeds by +999 and writes to $REV5_OUT."""
import os, sys, json, argparse, time, numpy as np, scipy.sparse as sp
from scipy.stats import norm, spearmanr
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("dev", "dstd", "pop", "rev"): sys.path.insert(0, os.path.join(HERE, "..", p))
from methods_v2 import bh, paule_mandel, pearson_fisher  # noqa: E402
from dstd import _ms, design_R  # noqa: E402
from revm import stats_from_cells, methods as rmethods  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--mode", required=True); ap.add_argument("--sub", default="P1"); ap.add_argument("--size", default="full")
ap.add_argument("--k0", type=int, default=0); ap.add_argument("--k1", type=int, default=20); a = ap.parse_args()
SMOKE = 999 if os.environ.get("SMOKE") else 0
OD = os.environ.get("REV5_OUT", os.path.join(HERE, "out")); os.makedirs(OD, exist_ok=True)
YG = ("UTY", "USP9Y", "ZFY", "DDX3Y", "KDM5D", "NLGN4Y", "TTTY14"); SEX = set(YG) | {"XIST"}


def load(sub):
    if sub.startswith("P"):
        P = os.path.join(HERE, "..", "pop"); X = sp.load_npz(os.path.join(P, f"pop_{sub}_X.npz")).tocsr(); M = np.load(os.path.join(P, f"pop_{sub}_meta.npz"), allow_pickle=True)
    else:
        P = os.path.join(HERE, "..", "onek"); X = sp.load_npz(os.path.join(P, f"onek_{sub}_X.npz")).tocsr(); M = np.load(os.path.join(P, f"onek_{sub}_meta.npz"), allow_pickle=True)
    return X, M["lib"].astype(float), M["donor"].astype(str), (M["sym"].astype(str) if "sym" in M else np.array([f"g{i}" for i in range(X.shape[1])]))


def centre(X, lib, don, lo=0.1, hi=10.0, exact=False):   # rev4_C.centre; exact=True: no shrinkage (Plan 6 L)
    ud = np.unique(don); rows = [np.where(don == u)[0] for u in ud]
    st = [_ms(X[r].toarray().astype(np.float64), lib[r]) for r in rows]
    m = np.array([s_[0] for s_ in st]); ps = np.array([s_[1] for s_ in st]); tau2, mu, _ = paule_mandel(m, ps)
    th = m if exact else mu[None, :] + tau2[None, :] / (tau2[None, :] + ps) * (m - mu[None, :])
    f = np.clip(mu[None, :] / np.maximum(th, 1e-30), lo, hi)
    idx = np.empty(X.shape[0], dtype=np.int64)
    for i, r in enumerate(rows): idx[r] = i
    return sp.csr_matrix(X.multiply(f[idx]))


def cells_of(Xm, lib, lab, units=None):
    pu = np.unique(lab) if units is None else units
    return pu, [(Xm[np.where(lab == u)[0]].toarray().astype(np.float64), lib[np.where(lab == u)[0]]) for u in pu]


def covW(pu):
    cv = json.load(open(os.path.join(HERE, "..", "rev4", "donor_covariates.json"))); rows = []
    for u in pu:
        c = cv[u]; rows.append([1.0, c["male"], c["age"], c["bank"] == "H", c["bank"] == "R", c["AD"] == "case", c["AD"] == "missing", c["eth"] == "AA", c["eth"] == "HL", c["eth"] == "OT"])
    return np.array(rows, float)


def resid(Z, W): return Z - W @ np.linalg.lstsq(W, Z, rcond=None)[0]


def tt(q): return 2 * norm.sf(np.abs(q.mean(0) / (q.std(0, ddof=1) / np.sqrt(q.shape[0]))))


def Rsplits(hs): return [(design_R(h["mA"], h["psiA"], h["kappaA"], h["SA"])[0], design_R(h["mB"], h["psiB"], h["kappaB"], h["SB"])[0]) for h in hs]


def xdsm_W(Rs, iu, W=None, chunk=10000):
    """XDS-M p-values, optionally with every half's R residualized on donor-level W (XDS-Mc / XDS-Md / XDS-Mad)."""
    A_, B_ = iu; Rw = Rs if W is None else [(resid(RA, W), resid(RB, W)) for RA, RB in Rs]; p = np.empty(len(A_)); mq = np.empty(len(A_))
    for c0 in range(0, len(A_), chunk):
        aa, bb = A_[c0:c0 + chunk], B_[c0:c0 + chunk]; q = sum(RA[:, aa] * RB[:, bb] for RA, RB in Rw) / len(Rw)
        p[c0:c0 + chunk] = tt(q); mq[c0:c0 + chunk] = q.mean(0)
    return p, mq


def xds_r(Rs, a, b, chunk=10000):
    dg = np.mean([(RA * RB).mean(0) for RA, RB in Rs], axis=0); r = np.empty(len(a))
    for c0 in range(0, len(a), chunk):
        aa, bb = a[c0:c0 + chunk], b[c0:c0 + chunk]
        num = np.mean([0.5 * ((RA[:, aa] * RB[:, bb]).mean(0) + (RA[:, bb] * RB[:, aa]).mean(0)) for RA, RB in Rs], axis=0)
        with np.errstate(invalid="ignore", divide="ignore"): rr = num / np.sqrt(dg[aa] * dg[bb])
        rr[(dg[aa] <= 0) | (dg[bb] <= 0)] = np.nan; r[c0:c0 + chunk] = rr
    return r


def fz(r, neff): return 2 * norm.sf(np.abs(np.arctanh(np.clip(r, -0.999999, 0.999999)) * np.sqrt(max(neff, 1.0))))


def corr_iu(Mx, iu):
    Mc = Mx - Mx.mean(0); C = Mc.T @ Mc; sd = np.sqrt(np.maximum(np.diag(C), 1e-300)); return (C / np.outer(sd, sd))[iu]


def pbw(m, nd, iu):
    L = np.log1p(1e6 * m); w = nd / nd.sum(); Lw = L - (w[:, None] * L).sum(0); Cw = (Lw * w[:, None]).T @ Lw; sw = np.sqrt(np.maximum(np.diag(Cw), 1e-300))
    return fz((Cw / np.outer(sw, sw))[iu], nd.sum() ** 2 / (nd ** 2).sum() - 3)


def pb_adj(L, iu, W):
    Rr = resid(L, W); Rr = Rr - Rr.mean(0); C = Rr.T @ Rr; sd = np.sqrt(np.maximum(np.diag(C), 1e-300)); r = np.clip((C / np.outer(sd, sd))[iu], -0.999999, 0.999999)
    return 2 * norm.sf(np.abs(np.arctanh(r) * np.sqrt(L.shape[0] - 3 - (W.shape[1] - 1))))


def rhat(cl): return np.array([(s ** 2).sum() / s.sum() for _, s in cl])


def base(cl, seed, iu):
    fs, hs = stats_from_cells(cl, seed); out = rmethods(fs, hs, iu); nd = np.array([len(s) for _, s in cl], float)
    pv = dict(PB=out["PB"], PBw=pbw(fs["m"], nd, iu), KX=out["KX"], XDSM=out["XDSM"]); return fs, hs, pv, nd


def summ(pv, null=None):
    o = {}
    for m_, v in pv.items():
        v = np.asarray(v, float); c = bh(v); nm = np.ones(len(v), bool) if null is None else null
        o[m_] = dict(V=int((c & nm).sum()), calls=int(c.sum()), t05=float((v[nm] < 0.05).mean()))
    return o


def lib_split_cells(Xm, lib, lab, pool, pu):
    """Per pseudo-donor: cells split by library (A = smaller pool ID). Returns kept units and (YA, sA, YB, sB)."""
    keep, halves = [], []
    for u in pu:
        r = np.where(lab == u)[0]; pl = np.unique(pool[r])
        if len(pl) != 2: continue
        rA, rB = r[pool[r] == pl[0]], r[pool[r] == pl[1]]
        if len(rA) < 2 or len(rB) < 2: continue
        keep.append(u); halves.append((Xm[rA].toarray().astype(np.float64), lib[rA], Xm[rB].toarray().astype(np.float64), lib[rB]))
    return np.array(keep), halves


def lib_R(halves):
    h = {k: [] for k in ("mA", "psiA", "kappaA", "SA", "mB", "psiB", "kappaB", "SB")}
    for YA, sA, YB, sB in halves:
        for lab_, Y, s in (("A", YA, sA), ("B", YB, sB)):
            m2, ps2, k2, S2, _, _ = _ms(Y, s); h["m" + lab_].append(m2); h["psi" + lab_].append(ps2); h["kappa" + lab_].append(k2); h["S" + lab_].append(S2)
    return Rsplits([{k: np.array(v) for k, v in h.items()}])


si = {"P1": 0, "P2": 1, "P3": 2, "B1": 0, "B2": 1, "B3": 2}[a.sub]
X, lib, don, sym = load(a.sub); G = X.shape[1]; iu = np.triu_indices(G, 1); A_, B_ = iu
mu_cell = np.asarray(X.sum(0)).ravel() / X.shape[0]; elig = np.where((mu_cell >= 0.1) & ~np.isin(sym, list(SEX)))[0]
pool = np.load(os.path.join(HERE, f"pool_{a.sub}.npz"), allow_pickle=True)["pool"].astype(str) if a.sub.startswith("P") else None


def fn_of(k): return os.path.join(OD, f"{a.mode}_{a.sub}{'_' + a.size if a.mode == 'S' else ''}_k{k}.json")


if a.mode == "S":
    pu0, cl0 = cells_of(X, lib, don); st0 = [_ms(Y, s) for Y, s in cl0]; m0 = np.array([t[0] for t in st0]); p0 = np.array([t[1] for t in st0])
    tau2, mu0, _ = paule_mandel(m0, p0); s_ = float(np.clip(np.median(np.sqrt(tau2[elig]) / np.maximum(mu0[elig], 1e-30)), 0.1, 0.5)); del cl0, st0
    RHO = (0.1, 0.2, 0.3, 0.4); mods = np.random.default_rng(72_000_000 + si).choice(elig, 40, replace=False).reshape(4, 10)
    gm = -np.ones(G, int)
    for i, mm in enumerate(mods): gm[mm] = i
    true = (gm[A_] >= 0) & (gm[A_] == gm[B_]); null = ~true; tmod = np.where(true, gm[A_], -1)
    rlin = np.array([(np.exp(s_ ** 2 * r) - 1) / (np.exp(s_ ** 2) - 1) for r in RHO]); Xs = centre(X, lib, don)
    for k in range(a.k0, a.k1):
        fn = fn_of(k)
        if os.path.exists(fn): continue
        t0 = time.time(); seed = (70_000_000 if a.size == "full" else 71_000_000) + 10 * k + si + SMOKE; rng = np.random.default_rng(seed)
        lab = rng.permutation(don); pu = np.unique(lab)
        if a.size != "full": pu = pu[np.sort(rng.choice(len(pu), int(a.size), replace=False))]
        Dk = len(pu); z = rng.standard_normal((Dk, 4)); e = rng.standard_normal((Dk, G)); U = s_ * e
        for i in range(4): U[:, mods[i]] = s_ * (np.sqrt(RHO[i]) * z[:, [i]] + np.sqrt(1 - RHO[i]) * e[:, mods[i]])
        F = np.exp(U - s_ ** 2 / 2); _, cl = cells_of(Xs, lib, lab, units=pu); cl = [(Y * F[i][None, :], s) for i, (Y, s) in enumerate(cl)]
        fs, hs, pv, nd = base(cl, seed, iu); res = dict(s=s_, rho=list(RHO), rho_lin=rlin.tolist(), D=Dk, modules=[[str(sym[g]) for g in mm] for mm in mods])
        for m_, v in pv.items():
            c = bh(np.asarray(v, float)); V = int((c & null).sum()); R_ = int(c.sum())
            res[m_] = dict(calls=R_, V=V, FDP=V / max(R_, 1), power=float(c[true].mean()), power_mod=[float(c[tmod == i].mean()) for i in range(4)], t05=float((np.asarray(v)[null] < 0.05).mean()))
        Rs = Rsplits(hs); nidx = np.random.default_rng(seed + 5).choice(np.where(null)[0], 2000, replace=False); tidx = np.where(true)[0]
        L = np.log1p(1e6 * fs["m"]); Lc = L - L.mean(0); sd = np.sqrt((Lc ** 2).sum(0))
        def rpb(ix): return (Lc[:, A_[ix]] * Lc[:, B_[ix]]).sum(0) / (sd[A_[ix]] * sd[B_[ix]])
        xt, xn, pt, pn = xds_r(Rs, A_[tidx], B_[tidx]), xds_r(Rs, A_[nidx], B_[nidx]), rpb(tidx), rpb(nidx)
        res["est"] = dict(XDSr=dict(bias_mod=[float(np.nanmean(xt[tmod[tidx] == i]) - rlin[i]) for i in range(4)], rmse_mod=[float(np.sqrt(np.nanmean((xt[tmod[tidx] == i] - rlin[i]) ** 2))) for i in range(4)],
                                     null_mean=float(np.nanmean(xn)), null_absmean=float(np.nanmean(np.abs(xn))), undefined=int(np.isnan(xt).sum() + np.isnan(xn).sum())),
                          rPB=dict(bias_mod=[float(np.mean(pt[tmod[tidx] == i]) - rlin[i]) for i in range(4)], rmse_mod=[float(np.sqrt(np.mean((pt[tmod[tidx] == i] - rlin[i]) ** 2))) for i in range(4)],
                                   null_mean=float(np.mean(pn)), null_absmean=float(np.mean(np.abs(pn)))))
        res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w")); print("S", a.sub, a.size, k, {m_: (res[m_]["calls"], round(res[m_]["FDP"], 3), round(res[m_]["power"], 2)) for m_ in pv}, res["sec"], "s", flush=True)

elif a.mode in ("D", "L", "C"):
    Xs = centre(X, lib, don, exact=(a.mode == "L"))
    if a.mode == "D":
        edges = np.quantile(np.log(lib), np.linspace(0, 1, 11)[1:-1]); dbin = np.digitize(np.log(lib), edges)
        sd_real = float(np.std(np.log(rhat(cells_of(X, lib, don)[1]))))
        lab0 = np.random.default_rng(61_000_000 + si).permutation(don); sd_plain = float(np.std(np.log(rhat(cells_of(X, lib, lab0)[1]))))
    if a.mode == "C":
        g200 = np.random.default_rng(76_000_000 + si).choice(elig, 200, replace=False); pm = np.zeros(G, bool); pm[g200] = True; planted = pm[A_] & pm[B_]
    for k in range(a.k0, a.k1):
        fn = fn_of(k)
        if os.path.exists(fn): continue
        t0 = time.time(); seed = {"D": 73_000_000, "L": 74_000_000, "C": 75_000_000}[a.mode] + 10 * k + si + SMOKE; rng = np.random.default_rng(seed); res = {}
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
        pu, cl = cells_of(Xs, lib, lab)
        if a.mode == "C":
            W = covW(pu); case = W[:, 5] == 1; cl = [((Y * np.where(pm, 1.2, 1.0)[None, :]) if case[i] else Y, s) for i, (Y, s) in enumerate(cl)]
        fs, hs, pv, nd = base(cl, seed, iu); Rs = Rsplits(hs)
        if a.mode == "D":
            rh = rhat(cl); Wd = np.column_stack([np.ones(len(cl)), np.log(rh), np.log(nd)]); pv["XDSMd"] = xdsm_W(Rs, iu, Wd)[0]
            res["sd_logrhat"] = dict(this=float(np.std(np.log(rh))), real=sd_real, plain_null=sd_plain)
        if a.mode == "L":
            ku, halves = lib_split_cells(Xs, lib, lab, pool, pu); pv["XDSMlib"] = xdsm_W(lib_R(halves), iu)[0]; res["D_lib"] = int(len(ku))
        if a.mode == "C":
            pv["PBc"] = pb_adj(np.log1p(1e6 * fs["m"]), iu, W); pv["XDSMad"] = xdsm_W(Rs, iu, W[:, [0, 5, 6]])[0]; pv["XDSMc"] = xdsm_W(Rs, iu, W)[0]
            for m_, v in pv.items():
                v = np.asarray(v, float); c = bh(v); res[m_] = dict(calls=int(c.sum()), calls_planted=int((c & planted).sum()), calls_other=int((c & ~planted).sum()), t05_planted=float((v[planted] < 0.05).mean()), t05_other=float((v[~planted] < 0.05).mean()))
            res["n_case"] = int(case.sum())
        else:
            res.update(summ(pv))
        res["D"] = int(len(pu)); res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w"))
        print(a.mode, a.sub, k, {m_: res[m_]["calls"] for m_ in pv}, res["sec"], "s", flush=True)

elif a.mode == "R":
    t0 = time.time(); pu, cl = cells_of(X, lib, don); fs, hs, pv, nd = base(cl, 44_000_001, iu); Rs = Rsplits(hs); res = dict(D=int(len(pu)))
    cX = bh(pv["XDSM"]); cP = bh(pv["PB"]); res["calls"] = dict(PB=int(cP.sum()), XDSM=int(cX.sum()), KX=int(bh(pv["KX"]).sum()), PBw=int(bh(pv["PBw"]).sum()))
    if a.sub.startswith("P"):
        # R-2 effect sizes
        xr = xds_r(Rs, A_, B_); L = np.log1p(1e6 * fs["m"]); rP = corr_iu(L, iu)
        num = np.zeros((G, G)); vx = np.zeros(G)
        for Y, s in cl:
            if len(s) < 2: continue
            x = Y / s[:, None]; xc = x - x.mean(0); num += xc.T @ xc; vx += (xc ** 2).sum(0)
        rw = (num / np.sqrt(np.outer(vx, vx)))[iu]; dlt = rP - xr; ok = np.isfinite(dlt) & np.isfinite(rw)
        q = np.quantile(rw[ok], [0.1, 0.9]); res["R2"] = dict(spearman_delta_rw=float(spearmanr(dlt[ok], rw[ok])[0]), spearman_rPB_XDSr=float(spearmanr(rP[ok], xr[ok])[0]),
                                                              median_delta_top_rw=float(np.median(dlt[ok][rw[ok] >= q[1]])), median_delta_bottom_rw=float(np.median(dlt[ok][rw[ok] <= q[0]])),
                                                              median_rPB=float(np.median(rP)), median_XDSr=float(np.nanmedian(xr)), undefined=int((~np.isfinite(xr)).sum()))
        # R-3 XDS-Md and XDS-Mlib
        rh = rhat(cl); Wd = np.column_stack([np.ones(len(cl)), np.log(rh), np.log(nd)]); cD = bh(xdsm_W(Rs, iu, Wd)[0])
        ku, halves = lib_split_cells(X, lib, don, pool, pu); pL, mqL = xdsm_W(lib_R(halves), iu); cL = bh(pL)
        kidx = np.isin(pu, ku); RsK = [(RA[kidx], RB[kidx]) for RA, RB in Rs[:1]]; c1 = bh(xdsm_W(RsK, iu)[0])
        res["R3"] = dict(XDSMd=int(cD.sum()), XDSMd_keep=float((cD & cX).sum() / max(cX.sum(), 1)), XDSMlib=int(cL.sum()), D_lib=int(len(ku)), XDSt1_libdonors=int(c1.sum()),
                         XDSMlib_keep_of_XDSM=float((cL & cX).sum() / max(cX.sum(), 1)), XDSMlib_keep_of_t1=float((cL & c1).sum() / max(c1.sum(), 1)))
        # R-4 sex diagnostics
        g = list(sym); yi = [g.index(y) for y in YG if y in g]; xi = g.index("XIST")
        pairs = {"XIST_Y": np.array([i for i in range(len(A_)) if (g[A_[i]] == "XIST" and g[B_[i]] in YG) or (g[B_[i]] == "XIST" and g[A_[i]] in YG)], int),
                 "Y_Y": np.array([i for i in range(len(A_)) if g[A_[i]] in YG and g[B_[i]] in YG], int)}
        def sexrep(p, mq=None):
            c = bh(p); o = {nm: dict(called=int(c[ix].sum()), n=int(len(ix)), median_p=float(np.median(p[ix]))) for nm, ix in pairs.items()}
            if mq is not None: o["XIST_Y"]["mean_q_sign"] = float(np.sign(np.mean(mq[pairs["XIST_Y"]])))
            o["calls_total"] = int(c.sum()); return o
        W = covW(pu); male = W[:, 1] == 1; Ysum = np.array([Y[:, yi].sum() for Y, _ in cl]); Xx = np.array([Y[:, xi].sum() for Y, _ in cl])
        emale = np.log((Ysum + 1) / (Xx + 1)) > 0; mism = emale != male; R4 = dict(mismatch=int(mism.sum()), mismatch_male_meta=int((mism & male).sum()), mismatch_female_meta=int((mism & ~male).sum()))
        pc, mqc = xdsm_W(Rs, iu, W); R4["all_XDSMc"] = sexrep(pc, mqc); R4["all_sexonly"] = sexrep(xdsm_W(Rs, iu, W[:, :2])[0])
        keep = ~mism; pu2 = pu[keep]; cl2 = [cl[i] for i in np.where(keep)[0]]; _, hs2 = stats_from_cells(cl2, 44_000_001); Rs2 = Rsplits(hs2); W2 = W[keep]
        p2, mq2 = xdsm_W(Rs2, iu, W2); R4["nomismatch_XDSMc"] = sexrep(p2, mq2); R4["nomismatch_sexonly"] = sexrep(xdsm_W(Rs2, iu, W2[:, :2])[0])
        cl3 = []; nrem = 0
        for i in np.where(keep)[0]:
            Y, s = cl[i]; ys = Y[:, yi].sum(1); xs = Y[:, xi]; bad = ((xs >= 2) & (ys == 0)) if male[i] else ((ys >= 2) & (xs == 0)); nrem += int(bad.sum())
            cl3.append((Y[~bad], s[~bad]) if (~bad).sum() >= 2 else (Y, s))
        _, hs3 = stats_from_cells(cl3, 44_000_001); Rs3 = Rsplits(hs3); p3, mq3 = xdsm_W(Rs3, iu, W2); R4["clean_XDSMc"] = sexrep(p3, mq3); R4["clean_sexonly"] = sexrep(xdsm_W(Rs3, iu, W2[:, :2])[0]); R4["nuclei_removed"] = nrem
        for nm, msk in (("male", male & keep), ("female", ~male & keep)):
            ii = np.where(msk)[0]; _, hsS = stats_from_cells([cl[i] for i in ii], 44_000_001); pS, mqS = xdsm_W(Rsplits(hsS), iu); R4["within_" + nm] = sexrep(pS, mqS); R4["within_" + nm]["D"] = int(len(ii))
        R4["lib"] = sexrep(pL, mqL)
        Rf, _ = design_R(fs["m"], fs["psi"], fs["kappa"], fs["S"]); ysc = Rf[:, yi].mean(1); mm_ = male & keep
        R4["LOY_spearman_age_Yscore_males"] = float(spearmanr(W[mm_, 2], ysc[mm_])[0]); R4["LOY_n"] = int(mm_.sum())
        res["R4"] = R4
        np.savez_compressed(os.path.join(OD, f"R_{a.sub}.npz"), xr=xr.astype(np.float32), rP=rP.astype(np.float32), rw=rw.astype(np.float32), pX=pv["XDSM"], pP=pv["PB"])
    res["sec"] = round(time.time() - t0); json.dump(res, open(os.path.join(OD, f"R_{a.sub}.json"), "w"), indent=1); print("R", a.sub, json.dumps(res)[:2000], flush=True)
