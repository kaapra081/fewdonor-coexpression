"""PLAN7 Amendment 2 (post hoc; lupus/AMENDMENT2.md). python lupus_amend2.py --mode F1|F2|F2B|FR --sub LM|FC [--k0 0 --k1 20]
Pool-composition covariates (24 pools of lanes with identical donor sets; lupus_pools.json) for XDS-M, XDS-Mlib, pseudobulk and XDS-Mgu."""
import os, sys, json, argparse, time, numpy as np
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import lupus_run as L  # noqa: E402
from methods_v2 import bh  # noqa: E402
from lupus_amend1 import resid  # noqa: E402


def pb_adj(Lg, iu, W):   # rev4_C.pb_adj
    Rr = resid(Lg, W); Rr = Rr - Rr.mean(0); C = Rr.T @ Rr; sd = np.sqrt(np.maximum(np.diag(C), 1e-300)); r = np.clip((C / np.outer(sd, sd))[iu], -0.999999, 0.999999)
    return 2 * norm.sf(np.abs(np.arctanh(r) * np.sqrt(Lg.shape[0] - 3 - (W.shape[1] - 1))))


def dc_pool(hs, grp, iu, W):
    A_, B_ = iu
    def qu(h, g):
        mA = resid(h["mA"], W)[g]; mB = resid(h["mB"], W)[g]; return (mA - mA.mean(0))[:, A_] * (mB - mB.mean(0))[:, B_]
    u1 = sum(qu(h, grp) for h in hs) / len(hs); u0 = sum(qu(h, ~grp) for h in hs) / len(hs); return L.welch(u1, u0)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--mode", required=True); ap.add_argument("--sub", required=True); ap.add_argument("--k0", type=int, default=0); ap.add_argument("--k1", type=int, default=20); a = ap.parse_args()
    OD = L.OD; os.makedirs(OD, exist_ok=True); si = {"LM": 0, "FC": 1}[a.sub]; SM = L.SMOKE
    X, lib, don, library, sym, disease = L.load(a.sub); G = X.shape[1]; iu = np.triu_indices(G, 1); A_, B_ = iu
    sle = {d: (t == SLE) for d, t in zip(don, disease)} if False else {d: (t == L.SLE) for d, t in zip(don, disease)}; isg = np.isin(sym, L.ISG); ii = isg[A_] & isg[B_]
    pm = json.load(open(os.path.join(L.DATA, "lupus_pools.json")))["lib2pool"]; pools = sorted(set(pm.values())); pidx = {p: i for i, p in enumerate(pools)}
    lcode = np.unique(library, return_inverse=True)[1]; ulib = np.unique(library); pcode = np.array([pidx[pm[l]] for l in library]); NP = len(pools)
    Xs = L.centre(X, lib, don, exact=(a.mode == "F1")) if a.mode != "FR" else None

    def Wpool(rows):
        F = np.zeros((len(rows), NP))
        for i, r in enumerate(rows): F[i] = np.bincount(pcode[r], minlength=NP) / len(r)
        return np.column_stack([np.ones(len(rows)), F[:, 1:]])

    def run_set(cl, rows, fs, hs, pv0, seed, with_lib, with_dc, grp):
        Rs = L.Rsplits(hs); W = Wpool(rows); pv = dict(XDSM=pv0["XDSM"], XDSMpool=L.xdsm_W([(resid(RA, W), resid(RB, W)) for RA, RB in Rs], iu)[0]); ex = {}
        pv["PB"] = pv0["PB"]; pv["PBpool"] = pb_adj(np.log1p(1e6 * fs["m"]), iu, W)
        if with_lib:
            kept, hl = L.lib_halves(cl, [lcode[r] for r in rows], seed); Rl = L.Rsplits(hl); pv["XDSMlib"] = L.xdsm_W(Rl, iu)[0]
            pv["XDSMlibpool"] = L.xdsm_W([(resid(RA, W[kept]), resid(RB, W[kept])) for RA, RB in Rl], iu)[0]; ex["D_lib"] = int(len(kept))
        dc = None
        if with_dc:
            o = L.dc_from_stats(fs, hs, grp, iu); pg, sg = dc_pool(hs, grp, iu, W); dc = dict(PBDC=o["PBDC"], XDSMgu=o["XDSMgu"], XDSMgupool=pg, sPBDC=o["sPBDC"], sXDSMgu=o["sXDSMgu"], sXDSMgupool=sg)
        return pv, dc, ex

    if a.mode in ("F1", "F2", "F2B"):
        for k in range(a.k0, a.k1):
            fn = os.path.join(OD, f"{a.mode}_{a.sub}_k{k}.json")
            if os.path.exists(fn): continue
            t0 = time.time(); seed = {"F1": 85_000_000, "F2": 80_000_000, "F2B": 81_000_000}[a.mode] + 10 * k + si + SM; rng = np.random.default_rng(seed)
            if a.mode == "F1":
                lab = don.copy()
                for p_ in np.unique(library):
                    ix = np.where(library == p_)[0]; lab[ix] = rng.permutation(don[ix])
            else:
                lab = rng.permutation(don)
            pu = np.unique(lab); rows = [np.where(lab == u)[0] for u in pu]; cl = [(Xs[r].toarray().astype(np.float64), lib[r]) for r in rows]; grp = np.array([sle[u] for u in pu])
            fs, hs, pv0, nd = L.base(cl, seed, iu); pv, dc, ex = run_set(cl, rows, fs, hs, pv0, seed, a.mode == "F1", a.mode in ("F1", "F2B"), grp)
            if a.mode == "F2B": pv = {}
            if a.mode == "F2": pv = {m: pv[m] for m in ("XDSM", "XDSMpool", "PBpool")}
            res = L.summ(pv) if pv else {}
            if dc is not None: res["DC"] = {m: dict(calls=int(bh(dc[m]).sum()), t05=float((dc[m] < 0.05).mean())) for m in ("PBDC", "XDSMgu", "XDSMgupool")}
            if a.mode == "F2B": res["DC"].pop("PBDC")
            res.update(ex); res["D"] = int(len(pu)); res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w"))
            print(a.mode, a.sub, k, {m: v["calls"] for m, v in res.items() if isinstance(v, dict) and "calls" in v}, {m: v["calls"] for m, v in res.get("DC", {}).items()}, flush=True)
    else:
        fn = os.path.join(OD, f"FR_{a.sub}.json")
        if os.path.exists(fn): sys.exit()
        t0 = time.time(); pu = np.unique(don); rows = [np.where(don == u)[0] for u in pu]; cl = [(X[r].toarray().astype(np.float64), lib[r]) for r in rows]; grp = np.array([sle[u] for u in pu])
        fs, hs, pv0, nd = L.base(cl, 44_000_001, iu); pv, dc, ex = run_set(cl, rows, fs, hs, pv0, 44_000_001, True, True, grp)
        c = {m: bh(v) for m, v in pv.items()}; cd = {m: bh(dc[m]) for m in ("PBDC", "XDSMgu", "XDSMgupool")}
        res = dict(D=int(len(pu)), calls={m: int(v.sum()) for m, v in c.items()}, keep_of_XDSM={m: float((c[m] & c["XDSM"]).sum() / max(c["XDSM"].sum(), 1)) for m in c}, t05={m: float((v < 0.05).mean()) for m, v in pv.items()}, **ex)
        res["DC"] = dict(calls={m: int(v.sum()) for m, v in cd.items()}, pos={m: int((cd[m] & (dc["s" + m] > 0)).sum()) for m in cd}, isg_called={m: int((cd[m] & ii).sum()) for m in cd}, isg_called_pos={m: int((cd[m] & ii & (dc["s" + m] > 0)).sum()) for m in cd},
                         overlap_XDSMgu_pool=int((cd["XDSMgu"] & cd["XDSMgupool"]).sum()), keep_of_XDSMgu=float((cd["XDSMgu"] & cd["XDSMgupool"]).sum() / max(cd["XDSMgu"].sum(), 1)), t05={m: float((dc[m] < 0.05).mean()) for m in cd}, n_isg_pairs=int(ii.sum()))
        res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w"), indent=1); print("FR", a.sub, json.dumps(res), flush=True)
