"""PLAN7 analyses (lupus/PLAN7.md, sha256 in PLAN7.sha256).
python lupus_run.py --mode A|A2|B|D|R --sub LM|FC [--k0 0 --k1 20]
A: donor-centred null (PB, PBw, KX, XDS-M); A2: same with 100 pseudo-donors; B: SLE-label null (PB-DC, XDS-Mgu);
D: library-preserving null (PB, PBw, KX, XDS-M, XDS-Mlib; PB-DC, XDS-Mgu as D3); R: real donors (C-1, C-2, D-R).
Outputs $LUPUS_OUT/<mode>_<sub>_k<k>.json (R: R_<sub>.json/.npz). SMOKE=1 offsets seeds by +999. Resumable (skips existing outputs)."""
import os, sys, json, argparse, time, numpy as np, scipy.sparse as sp
from scipy.stats import norm, spearmanr
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.environ.get("LUPUS_REPO", os.path.join(HERE, ".."))
for p in ("dev", "dstd", "pop", "rev", "xsplit"): sys.path.insert(0, os.path.join(REPO, p))
from methods_v2 import bh, paule_mandel, pearson_fisher  # noqa: E402
from dstd import _ms, design_R  # noqa: E402
from revm import stats_from_cells  # noqa: E402
from revm import methods as rmethods  # noqa: E402
DATA = os.environ.get("LUPUS_DATA", HERE); OD = os.environ.get("LUPUS_OUT", os.path.join(HERE, "out")); SMOKE = 999 if os.environ.get("SMOKE") else 0
ISG = ["ISG15", "IFI6", "IFI44L", "IFI44", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "IFIT1", "IFIT2", "IFIT3", "IFIT5", "IFI27", "RSAD2", "IFITM1", "IFITM3", "STAT1", "IRF7", "XAF1", "LY6E", "USP18",
       "EPSTI1", "SIGLEC1", "HERC5", "CMPK2", "SAMD9L", "PLSCR1", "IFI35", "BST2", "TRIM22"]
SLE = "systemic lupus erythematosus"


def load(key):
    X = sp.load_npz(os.path.join(DATA, f"lupus_{key}_X.npz")).tocsr().astype(np.float32); M = np.load(os.path.join(DATA, f"lupus_{key}_meta.npz"), allow_pickle=True)
    return X, M["lib"].astype(float), M["donor"].astype(str), M["library"].astype(str), M["sym"].astype(str), M["disease"].astype(str)


def centre(X, lib, don, lo=0.1, hi=10.0, exact=False):   # rev4_C.centre / rev5.centre (exact=True: no shrinkage)
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


def tt(q): return 2 * norm.sf(np.abs(q.mean(0) / (q.std(0, ddof=1) / np.sqrt(q.shape[0]))))


def Rsplits(hs): return [(design_R(h["mA"], h["psiA"], h["kappaA"], h["SA"])[0], design_R(h["mB"], h["psiB"], h["kappaB"], h["SB"])[0]) for h in hs]


def xdsm_W(Rs, iu, chunk=10000):
    A_, B_ = iu; p = np.empty(len(A_)); mq = np.empty(len(A_))
    for c0 in range(0, len(A_), chunk):
        aa, bb = A_[c0:c0 + chunk], B_[c0:c0 + chunk]; q = sum(RA[:, aa] * RB[:, bb] for RA, RB in Rs) / len(Rs)
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


def pbw(m, nd, iu):   # rev4.pb_family PBw
    L = np.log1p(1e6 * m); w = nd / nd.sum(); Lw = L - (w[:, None] * L).sum(0); Cw = (Lw * w[:, None]).T @ Lw; sw = np.sqrt(np.maximum(np.diag(Cw), 1e-300))
    return fz((Cw / np.outer(sw, sw))[iu], nd.sum() ** 2 / (nd ** 2).sum() - 3)


def base(cl, seed, iu):
    fs, hs = stats_from_cells(cl, seed); out = rmethods(fs, hs, iu); nd = np.array([len(s) for _, s in cl], float)
    return fs, hs, dict(PB=out["PB"], PBw=pbw(fs["m"], nd, iu), KX=out["KX"], XDSM=out["XDSM"]), nd


def summ(pv):
    o = {}
    for m_, v in pv.items():
        v = np.asarray(v, float); c = bh(v); o[m_] = dict(calls=int(c.sum()), t05=float((v < 0.05).mean()))
    return o


def pear(Mx, iu):   # rev2_D.pear
    Mc = Mx - Mx.mean(0); cov = Mc.T @ Mc / (len(Mx) - 1); sd = np.sqrt(np.diag(cov)); return np.clip((cov / np.outer(sd, sd))[iu], -0.999999, 0.999999)


def welch(q1, q0):
    z = (q1.mean(0) - q0.mean(0)) / np.sqrt(q1.var(0, ddof=1) / len(q1) + q0.var(0, ddof=1) / len(q0)); return 2 * norm.sf(np.abs(z)), np.sign(z)


def dc_from_stats(fs, hs, grp, iu):
    """rev2_D PB-DC and REV2 Amendment 1 XDS-Mgu; group 1 = grp True (SLE). sign +1: correlation larger in group 1."""
    A_, B_ = iu; L = np.log1p(1e6 * fs["m"]); n1, n0 = int(grp.sum()), int((~grp).sum())
    zpb = (np.arctanh(pear(L[grp], iu)) - np.arctanh(pear(L[~grp], iu))) / np.sqrt(1 / (n1 - 3) + 1 / (n0 - 3))
    out = dict(PBDC=2 * norm.sf(np.abs(zpb)), sPBDC=np.sign(zpb))
    def qu(h, g):
        mA = h["mA"][g]; mB = h["mB"][g]; return (mA - mA.mean(0))[:, A_] * (mB - mB.mean(0))[:, B_]
    u1 = sum(qu(h, grp) for h in hs) / len(hs); u0 = sum(qu(h, ~grp) for h in hs) / len(hs)
    out["XDSMgu"], out["sXDSMgu"] = welch(u1, u0); out["n1"], out["n0"] = n1, n0
    return out


def dc_summ(o):
    return {m_: dict(calls=int(bh(o[m_]).sum()), t05=float((o[m_] < 0.05).mean())) for m_ in ("PBDC", "XDSMgu")}


def lib_halves(cl, lcodes, seed, S=5, tries=100):
    """XDS-Mlib: per split, libraries of a pseudo-donor go to half A or B by fair coins (rng seeded [seed, 17, split, donor index]);
    accepted if both halves hold >= 2 cells (<= `tries` draws); a donor without an accepted draw in some split is dropped from the run."""
    D = len(cl); parts = [[None] * S for _ in range(D)]; ok = np.ones(D, bool)
    for i in range(D):
        lc = lcodes[i]; ul = np.unique(lc)
        if len(ul) < 2: ok[i] = False; continue
        pos = np.searchsorted(ul, lc)
        for s_ in range(S):
            rg = np.random.default_rng([int(seed), 17, s_, i]); found = None
            for _ in range(tries):
                inA = (rg.random(len(ul)) < 0.5)[pos]
                if inA.sum() >= 2 and (~inA).sum() >= 2: found = inA; break
            if found is None: ok[i] = False; break
            parts[i][s_] = found
    kept = np.where(ok)[0]; hs = []
    for s_ in range(S):
        h = {k: [] for k in ("mA", "psiA", "kappaA", "SA", "mB", "psiB", "kappaB", "SB")}
        for i in kept:
            Y, s = cl[i]; inA = parts[i][s_]
            for lab_, msk in (("A", inA), ("B", ~inA)):
                m2, ps2, k2, S2, _, _ = _ms(Y[msk], s[msk]); h["m" + lab_].append(m2); h["psi" + lab_].append(ps2); h["kappa" + lab_].append(k2); h["S" + lab_].append(S2)
        hs.append({k: np.array(v) for k, v in h.items()})
    return kept, hs


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--mode", required=True); ap.add_argument("--sub", required=True); ap.add_argument("--k0", type=int, default=0); ap.add_argument("--k1", type=int, default=20); a = ap.parse_args()
    os.makedirs(OD, exist_ok=True); si = {"LM": 0, "FC": 1}[a.sub]
    X, lib, don, library, sym, disease = load(a.sub); G = X.shape[1]; iu = np.triu_indices(G, 1); A_, B_ = iu
    sle = {d: (t == SLE) for d, t in zip(don, disease)}; isg = np.isin(sym, ISG)
    def fn_of(k): return os.path.join(OD, f"{a.mode}_{a.sub}_k{k}.json")
    if a.mode in ("A", "A2", "B", "D"):
        Xs = centre(X, lib, don, exact=(a.mode == "D")); lcode = np.unique(library, return_inverse=True)[1]
        for k in range(a.k0, a.k1):
            fn = fn_of(k)
            if os.path.exists(fn): continue
            t0 = time.time(); res = {}
            seed = {"A": 80_000_000, "B": 81_000_000, "D": 82_000_000, "A2": 83_000_000}[a.mode] + 10 * k + si + SMOKE; rng = np.random.default_rng(seed)
            if a.mode == "D":
                lab = don.copy()
                for p_ in np.unique(library):
                    ix = np.where(library == p_)[0]; lab[ix] = rng.permutation(don[ix])
            else:
                lab = rng.permutation(don)
            pu = np.unique(lab); units = pu
            if a.mode == "A2": units = pu[np.sort(rng.choice(len(pu), 100, replace=False))]
            rows = [np.where(lab == u)[0] for u in units]; cl = [(Xs[r].toarray().astype(np.float64), lib[r]) for r in rows]; nd = np.array([len(r) for r in rows], float)
            res["D"] = int(len(units)); res["nd_median"] = float(np.median(nd)); res["nd_min"] = int(nd.min())
            if a.mode == "B":
                fs, hs = stats_from_cells(cl, seed); grp = np.array([sle[u] for u in units]); res.update(dc_summ(dc_from_stats(fs, hs, grp, iu))); res["n1"], res["n0"] = int(grp.sum()), int((~grp).sum())
            else:
                fs, hs, pv, nd = base(cl, seed, iu)
                if a.mode == "D":
                    lcl = [lcode[r] for r in rows]; kept, hl = lib_halves(cl, lcl, seed); pv["XDSMlib"] = xdsm_W(Rsplits(hl), iu)[0]; res["D_lib"] = int(len(kept))
                res.update(summ(pv))
                if a.mode == "D":   # D3 (descriptive): SLE-label DC under the library-preserving null
                    grp = np.array([sle[u] for u in units]); res["DC"] = dc_summ(dc_from_stats(fs, hs, grp, iu)); res["DC"]["n1"], res["DC"]["n0"] = int(grp.sum()), int((~grp).sum())
            res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w")); print(a.mode, a.sub, k, {m_: (v["calls"] if isinstance(v, dict) and "calls" in v else v) for m_, v in res.items() if m_ not in ("DC",)}, flush=True)
    elif a.mode == "R":
        fn = os.path.join(OD, f"R_{a.sub}.json")
        if os.path.exists(fn): sys.exit()
        t0 = time.time(); pu, cl = cells_of(X, lib, don); fs, hs, pv, nd = base(cl, 44_000_001, iu); Rs = Rsplits(hs); res = dict(D=int(len(pu)), G=int(G), nd_median=float(np.median(nd)))
        cP, cW, cK, cX = (bh(pv[m_]) for m_ in ("PB", "PBw", "KX", "XDSM")); res["C1"] = dict(calls=dict(PB=int(cP.sum()), PBw=int(cW.sum()), KX=int(cK.sum()), XDSM=int(cX.sum())), n_pairs=int(len(A_)),
              frac_XDSM_in_PB=float((cX & cP).sum() / max(cX.sum(), 1)), frac_PB_in_XDSM=float((cX & cP).sum() / max(cP.sum(), 1)), PB_only=int((cP & ~cX).sum()), XDSM_only=int((cX & ~cP).sum()), both=int((cP & cX).sum()))
        xr = xds_r(Rs, A_, B_); L = np.log1p(1e6 * fs["m"]); rP = corr_iu(L, iu); num = np.zeros((G, G)); vx = np.zeros(G)
        for Y, s in cl:
            if len(s) < 2: continue
            x = Y / s[:, None]; xc = x - x.mean(0); num += xc.T @ xc; vx += (xc ** 2).sum(0)
        rw = (num / np.sqrt(np.outer(vx, vx)))[iu]; dlt = rP - xr; ok = np.isfinite(dlt) & np.isfinite(rw); q = np.quantile(rw[ok], [0.1, 0.9])
        res["C1"]["R2"] = dict(spearman_delta_rw=float(spearmanr(dlt[ok], rw[ok])[0]), spearman_rPB_XDSr=float(spearmanr(rP[ok], xr[ok])[0]), median_delta_top_rw=float(np.median(dlt[ok][rw[ok] >= q[1]])),
                               median_delta_bottom_rw=float(np.median(dlt[ok][rw[ok] <= q[0]])), median_rPB=float(np.median(rP)), median_XDSr=float(np.nanmedian(xr)), undefined=int((~np.isfinite(xr)).sum()))
        grp = np.array([sle[u] for u in pu]); o = dc_from_stats(fs, hs, grp, iu); cD = {m_: bh(o[m_]) for m_ in ("PBDC", "XDSMgu")}; ii = isg[A_] & isg[B_]
        res["C2"] = dict(n1=int(grp.sum()), n0=int((~grp).sum()), calls={m_: int(c.sum()) for m_, c in cD.items()}, overlap=int((cD["PBDC"] & cD["XDSMgu"]).sum()),
                         pos_sign={m_: int((c & (o["s" + m_] > 0)).sum()) for m_, c in cD.items()}, n_isg_pairs=int(ii.sum()), isg_called={m_: int((c & ii).sum()) for m_, c in cD.items()},
                         isg_called_pos={m_: int((c & ii & (o["s" + m_] > 0)).sum()) for m_, c in cD.items()}, median_nd_sle=float(np.median(nd[grp])), median_nd_control=float(np.median(nd[~grp])),
                         t05={m_: float((o[m_] < 0.05).mean()) for m_ in ("PBDC", "XDSMgu")})
        res["C1"]["isg_pairs_XDSM"] = dict(n=int(ii.sum()), called=int((cX & ii).sum()), pos=int((cX & ii & (np.nan_to_num(xr) > 0)).sum()), PB_called=int((cP & ii).sum()))
        lcode = np.unique(library, return_inverse=True)[1]; lcl = [lcode[np.where(don == u)[0]] for u in pu]; kept, hl = lib_halves(cl, lcl, 44_000_001); pL, mqL = xdsm_W(Rsplits(hl), iu); cL = bh(pL)
        kidx = np.zeros(len(pu), bool); kidx[kept] = True; RsK = [(RA[kidx], RB[kidx]) for RA, RB in Rs[:1]]; c1 = bh(xdsm_W(RsK, iu)[0])
        res["DR"] = dict(XDSMlib=int(cL.sum()), D_lib=int(len(kept)), XDSt1_libdonors=int(c1.sum()), XDSMlib_keep_of_XDSM=float((cL & cX).sum() / max(cX.sum(), 1)), XDSMlib_keep_of_t1=float((cL & c1).sum() / max(c1.sum(), 1)),
                         t1_keep_of_XDSM=float((c1 & cX).sum() / max(cX.sum(), 1)), XDSMlib_t05=float((pL < 0.05).mean()))
        np.savez_compressed(os.path.join(OD, f"R_{a.sub}.npz"), xr=xr.astype(np.float32), rP=rP.astype(np.float32), rw=rw.astype(np.float32), pX=pv["XDSM"], pP=pv["PB"], pPBDC=o["PBDC"], pXDSMgu=o["XDSMgu"], sPBDC=np.nan_to_num(o["sPBDC"]).astype(np.int8), sXDSMgu=np.nan_to_num(o["sXDSMgu"]).astype(np.int8))
        res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w"), indent=1); print("R", a.sub, json.dumps(res)[:1800], flush=True)
