"""REV4 Amendment 5 (post hoc). python rev4_A5.py --mode null20|top|sexctl --sub P1|P2|P3
null20: packed BH calls of null pseudobulk for the 20 donor-centred R3 reassignments -> out/A4/A5null20_<s>.npz
top:    indices of the 5,000 strongest pairs (null PB |r| k=0, real PB |r|, real XDS-M p) from out/R6_<s>.npz -> out/A4/A5top_<s>.npz
sexctl: XIST-Y / Y-Y pairs under XDS-M, XDS-Mc (all C1 covariates) and sex-only adjustment on real donors -> out/A5sex_<s>.json"""
import os, sys, json, argparse, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.argv_backup = list(sys.argv)
ap = argparse.ArgumentParser(); ap.add_argument("--mode", required=True); ap.add_argument("--sub", default="P1"); a = ap.parse_args()
sys.argv = [sys.argv[0], "--mode", "none", "--sub", a.sub]   # rev4_C parses its own args at import
import importlib.util  # noqa: E402
spec = importlib.util.spec_from_file_location("rc", os.path.join(HERE, "rev4_C.py")); rc = importlib.util.module_from_spec(spec); spec.loader.exec_module(rc)
OD = os.path.join(HERE, "out"); os.makedirs(os.path.join(OD, "A4"), exist_ok=True); si = rc.si; X, lib, don, sym = rc.X, rc.lib, rc.don, rc.sym; G = rc.G; iu = rc.iu; A_, B_ = iu

if a.mode == "null20":
    Xs = rc.centre(X, lib, don); calls = []
    for k in range(20):
        lab = np.random.default_rng(61_000_000 + 10 * k + si).permutation(don); pu = np.unique(lab)
        mm = np.array([np.asarray(Xs[np.where(lab == u)[0]].sum(0)).ravel() / lib[np.where(lab == u)[0]].sum() for u in pu])
        c0 = rc.bh(rc.pearson_fisher(np.log1p(1e6 * mm), iu)["p"]); calls.append(np.packbits(c0.astype(bool))); print("null20", a.sub, k, int(c0.sum()), flush=True)
    np.savez_compressed(os.path.join(OD, "A4", f"A5null20_{a.sub}.npz"), calls=np.array(calls), n=len(A_), sym=sym)
elif a.mode == "top":
    Z = np.load(os.path.join(OD, f"R6_{a.sub}.npz"), allow_pickle=True)
    o = lambda s: np.argsort(-s, kind="stable")[:5000]
    px = Z["real_pXDSM"]; zx = np.abs(Z["real_zXDSM"]); ordx = np.lexsort((-zx, px))[:5000]   # smallest p, ties by |z|
    np.savez_compressed(os.path.join(OD, "A4", f"A5top_{a.sub}.npz"), null_PB=o(np.abs(Z["null_rPB"])), real_PB=o(np.abs(Z["real_rPB"])), real_XDSM=ordx, sym=Z["sym"])
    print("top", a.sub, "XDS-M p at rank 5000:", float(px[ordx[-1]]), "ties at p<=1e-300:", int((px <= 1e-300).sum()), flush=True)
elif a.mode == "sexctl":
    Y = {"UTY", "USP9Y", "ZFY", "DDX3Y", "KDM5D", "NLGN4Y", "TTTY14"}; g = list(sym)
    pu, cl = rc.cells_of(X, lib, don); W = rc.covW(pu); fs, hs = rc.stats_from_cells(cl, 44_000_001)
    pv = rc.xdsm_variants(hs, iu, W=W); pvs = rc.xdsm_variants(hs, iu, W=W[:, :2])   # intercept + male
    cX, cC, cS = rc.bh(pv["XDSM"]), rc.bh(pv["XDSMc"]), rc.bh(pvs["XDSMc"])
    res = {}
    for nm, test in (("XIST_Y", lambda x, y: ("XIST" in (x, y)) and (x in Y or y in Y)), ("Y_Y", lambda x, y: x in Y and y in Y)):
        idx = np.array([i for i in range(len(A_)) if test(g[A_[i]], g[B_[i]])], int)
        res[nm] = dict(n=int(len(idx)), called_XDSM=int(cX[idx].sum()), called_XDSMc=int(cC[idx].sum()), called_sexonly=int(cS[idx].sum()),
                       median_p_XDSM=float(np.median(pv["XDSM"][idx])) if len(idx) else None, median_p_XDSMc=float(np.median(pv["XDSMc"][idx])) if len(idx) else None,
                       median_p_sexonly=float(np.median(pvs["XDSMc"][idx])) if len(idx) else None)
    res["calls_total"] = dict(XDSM=int(cX.sum()), XDSMc=int(cC.sum()), sexonly=int(cS.sum()))
    json.dump(res, open(os.path.join(OD, f"A5sex_{a.sub}.json"), "w"), indent=1); print("sexctl", a.sub, json.dumps(res), flush=True)
