"""PLAN7 Amendment 1 (post hoc; lupus/AMENDMENT1.md). python lupus_amend1.py --mode E1|E2|ER --sub LM|FC [--k0 0 --k1 20]
Batch-composition covariates for XDS-M and XDS-Mlib under the library-preserving null (E1), the donor-centred null (E2) and the real donors (ER)."""
import os, sys, json, argparse, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import lupus_run as L  # noqa: E402  (functions only; no side effects)
from methods_v2 import bh  # noqa: E402


def resid(Z, W): return Z - W @ np.linalg.lstsq(W, Z, rcond=None)[0]


def covariates(coh_codes, lib_codes, n_coh, n_lib):
    """rows: per pseudo-donor arrays of cohort / library codes of its cells. Returns W_coh (1 + n_coh - 1 columns) and W_pc (1 + 10)."""
    D = len(coh_codes); Fc = np.zeros((D, n_coh)); Fl = np.zeros((D, n_lib))
    for i, (cc, ll) in enumerate(zip(coh_codes, lib_codes)):
        Fc[i] = np.bincount(cc, minlength=n_coh) / len(cc); Fl[i] = np.bincount(ll, minlength=n_lib) / len(ll)
    Wc = np.column_stack([np.ones(D), Fc[:, 1:]]); Fl = Fl - Fl.mean(0); U, s, Vt = np.linalg.svd(Fl, full_matrices=False); Wp = np.column_stack([np.ones(D), U[:, :10]])
    return Wc, Wp


def p_of(Rs, iu, W=None, keep=None):
    if keep is not None: Rs = [(RA[keep], RB[keep]) for RA, RB in Rs]
    if W is not None: Rs = [(resid(RA, W), resid(RB, W)) for RA, RB in Rs]
    return L.xdsm_W(Rs, iu)[0]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--mode", required=True); ap.add_argument("--sub", required=True); ap.add_argument("--k0", type=int, default=0); ap.add_argument("--k1", type=int, default=20); a = ap.parse_args()
    OD = L.OD; os.makedirs(OD, exist_ok=True); si = {"LM": 0, "FC": 1}[a.sub]; SM = L.SMOKE
    X, lib, don, library, sym, disease = L.load(a.sub); G = X.shape[1]; iu = np.triu_indices(G, 1)
    M = np.load(os.path.join(L.DATA, f"lupus_{a.sub}_meta.npz"), allow_pickle=True); cohort = M["cohort"].astype(str)
    lcode = np.unique(library, return_inverse=True)[1]; ccode = np.unique(cohort, return_inverse=True)[1]; n_lib, n_coh = int(lcode.max() + 1), int(ccode.max() + 1)
    Xs = L.centre(X, lib, don, exact=(a.mode == "E1")) if a.mode != "ER" else None

    def variants(cl, rows, hs, seed, with_lib):
        Rs = L.Rsplits(hs); Wc, Wp = covariates([ccode[r] for r in rows], [lcode[r] for r in rows], n_coh, n_lib); pv = dict(XDSM=p_of(Rs, iu), XDSMcoh=p_of(Rs, iu, Wc), XDSMpc=p_of(Rs, iu, Wp)); extra = {}
        if with_lib:
            kept, hl = L.lib_halves(cl, [lcode[r] for r in rows], seed); Rl = L.Rsplits(hl); pv["XDSMlib"] = p_of(Rl, iu); pv["XDSMlibcoh"] = p_of(Rl, iu, Wc[kept]); pv["XDSMlibpc"] = p_of(Rl, iu, Wp[kept]); extra["D_lib"] = int(len(kept))
        return pv, extra

    if a.mode in ("E1", "E2"):
        for k in range(a.k0, a.k1):
            fn = os.path.join(OD, f"{a.mode}_{a.sub}_k{k}.json")
            if os.path.exists(fn): continue
            t0 = time.time(); seed = (84_000_000 if a.mode == "E1" else 80_000_000) + 10 * k + si + SM; rng = np.random.default_rng(seed)
            if a.mode == "E1":
                lab = don.copy()
                for p_ in np.unique(library):
                    ix = np.where(library == p_)[0]; lab[ix] = rng.permutation(don[ix])
            else:
                lab = rng.permutation(don)
            pu = np.unique(lab); rows = [np.where(lab == u)[0] for u in pu]; cl = [(Xs[r].toarray().astype(np.float64), lib[r]) for r in rows]
            fs, hs, pv0, nd = L.base(cl, seed, iu); pv, extra = variants(cl, rows, hs, seed, a.mode == "E1")
            if a.mode == "E1": pv["PB"] = pv0["PB"]
            res = L.summ(pv); res.update(extra); res["D"] = int(len(pu)); res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w"))
            print(a.mode, a.sub, k, {m: res[m]["calls"] for m in pv}, flush=True)
    else:
        fn = os.path.join(OD, f"ER_{a.sub}.json")
        if os.path.exists(fn): sys.exit()
        t0 = time.time(); pu = np.unique(don); rows = [np.where(don == u)[0] for u in pu]; cl = [(X[r].toarray().astype(np.float64), lib[r]) for r in rows]
        fs, hs, pv0, nd = L.base(cl, 44_000_001, iu); pv, extra = variants(cl, rows, hs, 44_000_001, True); c = {m: bh(v) for m, v in pv.items()}; c["PB"] = bh(pv0["PB"])
        res = dict(D=int(len(pu)), calls={m: int(v.sum()) for m, v in c.items()}, keep_of_XDSM={m: float((c[m] & c["XDSM"]).sum() / max(c["XDSM"].sum(), 1)) for m in c}, t05={m: float((v < 0.05).mean()) for m, v in pv.items()}, **extra)
        res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w"), indent=1); print("ER", a.sub, json.dumps(res), flush=True)
