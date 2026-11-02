"""Plan 6d (post hoc, rev5/PLAN6D.md). python rev5_d.py --mode D|A|O|S|R --sub P1|B1 [--k0 --k1] -> $REV5_OUT/d_<mode>_<sub>[_k<k>].json"""
import os, sys, json, argparse, time, importlib.util, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(); ap.add_argument("--mode", required=True); ap.add_argument("--sub", default="P1")
ap.add_argument("--k0", type=int, default=0); ap.add_argument("--k1", type=int, default=20); a = ap.parse_args()
sys.argv = [sys.argv[0], "--mode", "none", "--sub", a.sub]
spec = importlib.util.spec_from_file_location("r5", os.path.join(HERE, "rev5.py")); r5 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r5)
from methods_v2 import bh  # noqa: E402
from revm import stats_from_cells  # noqa: E402
X, lib, don, sym, G, iu, A_, B_, si, elig = r5.X, r5.lib, r5.don, r5.sym, r5.G, r5.iu, r5.A_, r5.B_, r5.si, r5.elig; OD = r5.OD


def variants(cl, hs, sel=None):
    nd = np.array([len(s) for _, s in cl], float)
    if sel is not None: nd = nd[sel]; hs = [{k: v[sel] for k, v in h.items()} for h in hs]
    Rs = r5.Rsplits(hs); nA = np.floor(nd / 2); nB = nd - nA; Rh = []
    for (RA, RB), h in zip(Rs, hs):
        WA = np.column_stack([np.ones(len(nd)), np.log(h["kappaA"] * h["SA"]), np.log(nA)]); WB = np.column_stack([np.ones(len(nd)), np.log(h["kappaB"] * h["SB"]), np.log(nB)])
        Rh.append((r5.resid(RA, WA), r5.resid(RB, WB)))
    return dict(XDSM=r5.xdsm_W(Rs, iu)[0], XDSMh=r5.xdsm_W(Rh, iu)[0])


def summ(pv, null=None, true=None):
    o = {}
    for m_, v in pv.items():
        c = bh(np.asarray(v, float)); nm = np.ones(len(v), bool) if null is None else null
        o[m_] = dict(calls=int(c.sum()), V=int((c & nm).sum()), t05=float((np.asarray(v)[nm] < 0.05).mean()))
        if true is not None: o[m_].update(FDP=float((c & nm).sum() / max(c.sum(), 1)), power=float(c[true].mean()))
    return o


if a.mode in ("A", "O", "D"):
    Xs = X if a.mode == "O" else r5.centre(X, lib, don)
    if a.mode == "D": dbin = np.digitize(np.log(lib), np.quantile(np.log(lib), np.linspace(0, 1, 11)[1:-1]))
    for k in range(a.k0, a.k1):
        fn = os.path.join(OD, f"d_{a.mode}_{a.sub}_k{k}.json")
        if os.path.exists(fn): continue
        t0 = time.time(); seed = {"A": 61_000_000, "O": 65_000_000, "D": 73_000_000}[a.mode] + 10 * k + si; rng = np.random.default_rng(seed)
        if a.mode == "D":
            lab = don.copy()
            for b in range(10):
                ix = np.where(dbin == b)[0]; lab[ix] = rng.permutation(don[ix])
        else:
            lab = rng.permutation(don)
        pu, cl = r5.cells_of(Xs, lib, lab); fs, hs = stats_from_cells(cl, seed); res = {"all": summ(variants(cl, hs))}
        if a.mode == "A":
            order = np.random.default_rng(seed + 1).permutation(len(pu))
            for D in (100, 200): res[str(D)] = summ(variants(cl, hs, np.sort(order[:D])))
        res["D"] = int(len(pu)); res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w")); print("d", a.mode, a.sub, k, {d: {m: v["calls"] for m, v in r.items()} for d, r in res.items() if isinstance(r, dict)}, res["sec"], "s", flush=True)

elif a.mode == "S":
    pu0, cl0 = r5.cells_of(X, lib, don); st0 = [r5._ms(Y, s) for Y, s in cl0]; m0 = np.array([t[0] for t in st0]); p0 = np.array([t[1] for t in st0])
    tau2, mu0, _ = r5.paule_mandel(m0, p0); s_ = float(np.clip(np.median(np.sqrt(tau2[elig]) / np.maximum(mu0[elig], 1e-30)), 0.1, 0.5)); del cl0, st0
    RHO = (0.1, 0.2, 0.3, 0.4); mods = np.random.default_rng(72_000_000 + si).choice(elig, 40, replace=False).reshape(4, 10); gm = -np.ones(G, int)
    for i, mm in enumerate(mods): gm[mm] = i
    true = (gm[A_] >= 0) & (gm[A_] == gm[B_]); null = ~true; Xs = r5.centre(X, lib, don)
    for k in range(a.k0, a.k1):
        fn = os.path.join(OD, f"d_S_{a.sub}_k{k}.json")
        if os.path.exists(fn): continue
        t0 = time.time(); seed = 70_000_000 + 10 * k + si; rng = np.random.default_rng(seed); lab = rng.permutation(don); pu = np.unique(lab)
        Dk = len(pu); z = rng.standard_normal((Dk, 4)); e = rng.standard_normal((Dk, G)); U = s_ * e
        for i in range(4): U[:, mods[i]] = s_ * (np.sqrt(RHO[i]) * z[:, [i]] + np.sqrt(1 - RHO[i]) * e[:, mods[i]])
        F = np.exp(U - s_ ** 2 / 2); _, cl = r5.cells_of(Xs, lib, lab, units=pu); cl = [(Y * F[i][None, :], s) for i, (Y, s) in enumerate(cl)]
        fs, hs = stats_from_cells(cl, seed); res = summ(variants(cl, hs), null, true); res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w"))
        print("d S", a.sub, k, {m: (v["calls"], round(v["FDP"], 3), round(v["power"], 2)) for m, v in res.items() if isinstance(v, dict)}, res["sec"], "s", flush=True)

elif a.mode == "R":
    t0 = time.time(); pu, cl = r5.cells_of(X, lib, don); fs, hs = stats_from_cells(cl, 44_000_001); pv = variants(cl, hs)
    cX = bh(pv["XDSM"]); cH = bh(pv["XDSMh"]); res = dict(D=int(len(pu)), XDSM=int(cX.sum()), XDSMh=int(cH.sum()), keep=float((cX & cH).sum() / max(cX.sum(), 1)), new=int((cH & ~cX).sum()))
    res["sec"] = round(time.time() - t0); json.dump(res, open(os.path.join(OD, f"d_R_{a.sub}.json"), "w")); print("d R", a.sub, res, flush=True)
