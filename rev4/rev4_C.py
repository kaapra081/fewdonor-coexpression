"""REV4 Amendment 3: python rev4_C.py --mode C1|C1real|C2|C3|C4|C5|C6 --sub P1 [--k0 --k1]   (outputs rev4/out/C*_*.json)"""
import os, sys, json, argparse, time, numpy as np, scipy.sparse as sp
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("dev", "dstd", "pop", "rev"): sys.path.insert(0, os.path.join(HERE, "..", p))
from methods_v2 import bh, paule_mandel, pearson_fisher  # noqa: E402
from dstd import _ms, design_R, donor_level_cells, methods_d  # noqa: E402
from revm import stats_from_cells, methods as rmethods  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("--mode", required=True); ap.add_argument("--sub", default="P1")
ap.add_argument("--k0", type=int, default=0); ap.add_argument("--k1", type=int, default=20); a = ap.parse_args()
OD = os.path.join(HERE, "out"); os.makedirs(OD, exist_ok=True)


def load(sub):
    if sub.startswith("P"):
        P = os.path.join(HERE, "..", "pop"); X = sp.load_npz(os.path.join(P, f"pop_{sub}_X.npz")).tocsr(); M = np.load(os.path.join(P, f"pop_{sub}_meta.npz"), allow_pickle=True)
    else:
        P = os.path.join(HERE, "..", "onek"); X = sp.load_npz(os.path.join(P, f"onek_{sub}_X.npz")).tocsr(); M = np.load(os.path.join(P, f"onek_{sub}_meta.npz"), allow_pickle=True)
    return X, M["lib"].astype(float), M["donor"].astype(str), (M["sym"].astype(str) if "sym" in M else None)


def centre(X, lib, don, lo=0.1, hi=10.0):
    ud = np.unique(don); rows = [np.where(don == u)[0] for u in ud]
    st = [_ms(X[r].toarray().astype(np.float64), lib[r]) for r in rows]
    m = np.array([s_[0] for s_ in st]); ps = np.array([s_[1] for s_ in st]); tau2, mu, _ = paule_mandel(m, ps)
    th = mu[None, :] + tau2[None, :] / (tau2[None, :] + ps) * (m - mu[None, :]); f = np.clip(mu[None, :] / np.maximum(th, 1e-30), lo, hi)
    idx = np.empty(X.shape[0], dtype=np.int64)
    for i, r in enumerate(rows): idx[r] = i
    return sp.csr_matrix(X.multiply(f[idx]))


def cells_of(Xm, lib, lab):
    pu = np.unique(lab); return pu, [(Xm[np.where(lab == u)[0]].toarray().astype(np.float64), lib[np.where(lab == u)[0]]) for u in pu]


def covW(pu):
    cv = json.load(open(os.path.join(HERE, "donor_covariates.json"))); rows = []
    for u in pu:
        c = cv[u]; rows.append([1.0, c["male"], c["age"], c["bank"] == "H", c["bank"] == "R", c["AD"] == "case", c["AD"] == "missing", c["eth"] == "AA", c["eth"] == "HL", c["eth"] == "OT"])
    return np.array(rows, float)


def resid(Z, W):
    return Z - W @ np.linalg.lstsq(W, Z, rcond=None)[0]


def tt(q):
    return 2 * norm.sf(np.abs(q.mean(0) / (q.std(0, ddof=1) / np.sqrt(q.shape[0]))))


def xdsm_variants(hs, iu, W=None, groups=None, chunk=10000):
    A_, B_ = iu; Rs = [(design_R(h["mA"], h["psiA"], h["kappaA"], h["SA"])[0], design_R(h["mB"], h["psiB"], h["kappaB"], h["SB"])[0]) for h in hs]
    out = {"XDSM": np.empty(len(A_))}
    if W is not None: out["XDSMc"] = np.empty(len(A_)); RsW = [(resid(RA, W), resid(RB, W)) for RA, RB in Rs]
    if groups is not None:
        out["XDSMb"] = np.empty(len(A_)); Gm = np.array([[g == b for b in np.unique(groups)] for g in groups], float); RsG = [(resid(RA, Gm), resid(RB, Gm)) for RA, RB in Rs]
    for c0 in range(0, len(A_), chunk):
        aa, bb = A_[c0:c0 + chunk], B_[c0:c0 + chunk]; sl = slice(c0, c0 + chunk)
        out["XDSM"][sl] = tt(sum(RA[:, aa] * RB[:, bb] for RA, RB in Rs) / len(Rs))
        if W is not None: out["XDSMc"][sl] = tt(sum(RA[:, aa] * RB[:, bb] for RA, RB in RsW) / len(RsW))
        if groups is not None: out["XDSMb"][sl] = tt(sum(RA[:, aa] * RB[:, bb] for RA, RB in RsG) / len(RsG))
    return out


def pb_adj(L, iu, W):
    Rr = resid(L, W); Rr = Rr - Rr.mean(0); C = Rr.T @ Rr; sd = np.sqrt(np.maximum(np.diag(C), 1e-300)); r = np.clip((C / np.outer(sd, sd))[iu], -0.999999, 0.999999)
    return 2 * norm.sf(np.abs(np.arctanh(r) * np.sqrt(L.shape[0] - 3 - (W.shape[1] - 1))))


def summ(pv): return {m: dict(V=int(bh(np.asarray(v, float)).sum()), t05=float((np.asarray(v) < 0.05).mean())) for m, v in pv.items()}


def top_sets(score, A_, B_, G, top=5000, nh=20):   # score: larger = stronger (|r| for PB, -log p for XDS-M)
    o = np.argsort(-score, kind="stable")[:top]; deg = np.bincount(np.concatenate([A_[o], B_[o]]), minlength=G); return set(o.tolist()), set(np.argsort(-deg, kind="stable")[:nh].tolist())


si = {"P1": 0, "P2": 1, "P3": 2, "B1": 0, "B2": 1, "B3": 2}[a.sub]
X, lib, don, sym = load(a.sub); G = X.shape[1]; iu = np.triu_indices(G, 1); A_, B_ = iu

if a.mode == "C1":
    Xs = centre(X, lib, don)
    for k in range(a.k0, a.k1):
        fn = os.path.join(OD, f"C1_{a.sub}_k{k}.json")
        if os.path.exists(fn): continue
        t0 = time.time(); seed = 66_000_000 + 10 * k + si; lab = np.random.default_rng(seed).permutation(don)
        pu, cl = cells_of(Xs, lib, lab); W = covW(pu); fs, hs = stats_from_cells(cl, seed); L = np.log1p(1e6 * fs["m"])
        pv = xdsm_variants(hs, iu, W=W); pv["PB"] = pearson_fisher(L, iu)["p"]; pv["PBc"] = pb_adj(L, iu, W)
        res = summ(pv); res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w")); print("C1", a.sub, k, {m: res[m]["V"] for m in pv}, res["sec"], "s", flush=True)

elif a.mode == "C1real":
    pu, cl = cells_of(X, lib, don); W = covW(pu); fs, hs = stats_from_cells(cl, 44_000_001)
    pv = xdsm_variants(hs, iu, W=W); cX = bh(pv["XDSM"]); cC = bh(pv["XDSMc"])
    res = dict(D=len(pu), XDSM=int(cX.sum()), XDSMc=int(cC.sum()), both=int((cX & cC).sum()), only_XDSM=int((cX & ~cC).sum()), only_XDSMc=int((cC & ~cX).sum()), p_covariates=int(W.shape[1]))
    json.dump(res, open(os.path.join(OD, f"C1real_{a.sub}.json"), "w")); print("C1real", a.sub, res, flush=True)

elif a.mode == "C2":
    src = json.load(open(os.path.join(HERE, "donor_source.json"))); bank = np.array([src[d_] for d_ in don]); Xs = centre(X, lib, don)
    mu = np.asarray(Xs.sum(0)).ravel() / lib.sum(); gfac = np.ones((3, G)); banks = ("H", "M", "R"); bstat = {}
    for i, b in enumerate(banks):
        ix = np.where(bank == b)[0]; mub = np.asarray(Xs[ix].sum(0)).ravel() / lib[ix].sum(); gfac[i] = np.clip(mu / np.maximum(mub, 1e-30), 0.1, 10.0)
        bstat[b] = dict(cells=int(len(ix)), mean_log_depth=float(np.log(lib[ix]).mean()), Es2_over_Es=float((lib[ix] ** 2).mean() / lib[ix].mean()))
    Xb = sp.csr_matrix(Xs.multiply(gfac[np.array([banks.index(b) for b in bank])]))
    for k in range(a.k0, a.k1):
        fn = os.path.join(OD, f"C2_{a.sub}_k{k}.json")
        if os.path.exists(fn): continue
        t0 = time.time(); seed = 64_000_000 + 10 * k + si; rng = np.random.default_rng(seed); lab = don.copy()
        for bk in ("H", "M", "R"):   # identical draw order to rev4.py R2
            ix = np.where(bank == bk)[0]; lab[ix] = rng.permutation(don[ix])
        pu, cl = cells_of(Xb, lib, lab); gb = np.array([src[u] for u in pu]); fs, hs = stats_from_cells(cl, seed)
        pv = xdsm_variants(hs, iu, groups=gb); pv["PB"] = pearson_fisher(np.log1p(1e6 * fs["m"]), iu)["p"]
        res = summ(pv); res["bank_stats"] = bstat; res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w")); print("C2", a.sub, k, {m: res[m]["V"] for m in pv}, res["sec"], "s", flush=True)

elif a.mode == "C3":
    pu, cl = cells_of(X, lib, don); res = {}
    f1, h1 = stats_from_cells(cl, 44_000_001); f2, h2 = stats_from_cells(cl, 44_000_002)
    nl = lambda p: -np.log(np.maximum(p, 1e-320))
    p1 = rmethods(f1, h1, iu)["XDSM"]; p2 = rmethods(f2, h2, iu)["XDSM"]; rPB = np.abs(pearson_fisher(np.log1p(1e6 * f1["m"]), iu)["est"])
    res["ties"] = dict(XDSM_p_zero=int((p1 <= 1e-320).sum()))
    T1, H1 = top_sets(nl(p1), A_, B_, G); T2, H2 = top_sets(nl(p2), A_, B_, G); TP, HP = top_sets(rPB, A_, B_, G)
    res["XDSM_seed_vs_seed"] = dict(top5000=len(T1 & T2), hubs20=len(H1 & H2)); res["PB_vs_XDSM"] = dict(top5000=len(TP & T1), hubs20=len(HP & H1))
    order = np.random.default_rng(42_000_000).permutation(len(pu)); halves = (np.sort(order[: len(pu) // 2]), np.sort(order[len(pu) // 2:])); TX, HX, TQ, HQ = [], [], [], []
    for j, rows in enumerate(halves):
        clh = [cl[i] for i in rows]; fh, hh = stats_from_cells(clh, 42_000_000 + si)
        t_, h_ = top_sets(nl(rmethods(fh, hh, iu)["XDSM"]), A_, B_, G); TX.append(t_); HX.append(h_)
        t_, h_ = top_sets(np.abs(pearson_fisher(np.log1p(1e6 * fh["m"]), iu)["est"]), A_, B_, G); TQ.append(t_); HQ.append(h_)
    res["XDSM_half_vs_half"] = dict(top5000=len(TX[0] & TX[1]), hubs20=len(HX[0] & HX[1])); res["PB_half_vs_half"] = dict(top5000=len(TQ[0] & TQ[1]), hubs20=len(HQ[0] & HQ[1]))
    if a.sub == "P1" and os.path.exists(os.path.join(OD, "R6_P1.npz")):
        Z = np.load(os.path.join(OD, "R6_P1.npz"), allow_pickle=True); c0 = Z["null_cPB"].astype(bool); cP = Z["real_cPB"].astype(bool); cX = Z["real_cXDSM"].astype(bool)
        pbo = cP & ~cX; t = [[int((pbo & c0).sum()), int((pbo & ~c0).sum())], [int((~pbo & c0).sum()), int((~pbo & ~c0).sum())]]
        res["PBonly_null_OR"] = dict(table=t, OR=float(t[0][0] * t[1][1] / max(t[0][1] * t[1][0], 1)))
    json.dump(res, open(os.path.join(OD, f"C3_{a.sub}.json"), "w"), indent=1); print("C3", a.sub, json.dumps(res), flush=True)

elif a.mode == "C4":
    from pop_core import subset  # noqa: E402
    Xs = centre(X, lib, don)
    for k in range(a.k0, a.k1):
        fn = os.path.join(OD, f"C4_{a.sub}_k{k}.json")
        if os.path.exists(fn): continue
        t0 = time.time(); seed = 50_000_000 + 10 * k + si; rng = np.random.default_rng(seed); lab = rng.permutation(don)   # as REV2 A
        pu = np.unique(lab); prow = [np.where(lab == u)[0] for u in pu]; cells = [(Xs[r].toarray(), lib[r].astype(np.float64)) for r in prow]
        d = donor_level_cells(cells, seed + 7); d["donors"] = pu; order = rng.permutation(len(pu)); D = 200; sd = subset(d, np.sort(order[:D]))
        oN = methods_d(sd, iu, seed=seed + D + 3); oP = methods_d(sd, iu, B=2_000_000, seed=seed + D + 3, perm_max_D=200)
        pv = dict(XDS_normal=oN["LCPXd"], LCPd_normal=oN["LCPd"], XDS_perm=oP["LCPXd"], LCPd_perm=oP["LCPd"])
        res = summ(pv); res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w")); print("C4", a.sub, k, {m: res[m]["V"] for m in pv}, res["sec"], "s", flush=True)

elif a.mode == "C5":
    ud = np.unique(don); num = np.zeros(G); den_x = np.zeros(G); den_l = 0.0; nc = 0
    for u in ud:
        r = np.where(don == u)[0]
        if len(r) < 2: continue
        ls = np.log(lib[r]); lc = ls - ls.mean(); x = X[r].toarray() / lib[r][:, None]; xc = x - x.mean(0)
        num += lc @ xc; den_x += (xc ** 2).sum(0); den_l += (lc ** 2).sum(); nc += len(r)
    c = num / np.sqrt(den_l * np.maximum(den_x, 1e-300)); ok = den_x > 0
    res = dict(coupling_median=float(np.median(c[ok])), frac_abs_gt_0045=float(np.mean(np.abs(c[ok]) > 0.045)), n_genes=int(ok.sum()), n_cells=int(nc), n_donors=int(len(ud)))
    json.dump(res, open(os.path.join(OD, f"C5_{a.sub}.json"), "w")); print("C5", a.sub, res, flush=True)

elif a.mode == "C6":
    Xs = centre(X, lib, don); lab = np.random.default_rng(61_000_000 + si).permutation(don); pu = np.unique(lab)   # first R3 reassignment
    mm = np.array([np.asarray(Xs[np.where(lab == u)[0]].sum(0)).ravel() / lib[np.where(lab == u)[0]].sum() for u in pu])
    c0 = bh(pearson_fisher(np.log1p(1e6 * mm), iu)["p"]); deg = np.bincount(np.concatenate([A_[c0], B_[c0]]), minlength=G)
    res = dict(calls=int(c0.sum()), degree={s_: int(v) for s_, v in zip(sym, deg)})
    json.dump(res, open(os.path.join(OD, f"C6_{a.sub}.json"), "w")); print("C6", a.sub, res["calls"], flush=True)
