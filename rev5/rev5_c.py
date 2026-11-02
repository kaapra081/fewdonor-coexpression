"""Plan 6c (post hoc, rev5/PLAN6C.md). python rev5_c.py --mode A|O|S --sub P1|B1 [--k0 --k1] -> $REV5_OUT/c_<mode>_<sub>_k<k>.json"""
import os, sys, json, argparse, time, importlib.util, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(); ap.add_argument("--mode", required=True); ap.add_argument("--sub", default="P1")
ap.add_argument("--k0", type=int, default=0); ap.add_argument("--k1", type=int, default=20); a = ap.parse_args()
sys.argv = [sys.argv[0], "--mode", "none", "--sub", a.sub]
spec = importlib.util.spec_from_file_location("r5", os.path.join(HERE, "rev5.py")); r5 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r5)
spec2 = importlib.util.spec_from_file_location("r5p", os.path.join(HERE, "rev5_post.py"))
from methods_v2 import bh  # noqa: E402
from revm import stats_from_cells  # noqa: E402
from dstd import design_R  # noqa: E402
X, lib, don, sym, G, iu, A_, B_, si, elig = r5.X, r5.lib, r5.don, r5.sym, r5.G, r5.iu, r5.A_, r5.B_, r5.si, r5.elig; OD = r5.OD


def ca_splits(hs, W):
    out = []
    for h in hs:
        pr = []
        for L_ in "AB":
            m = h["m" + L_]; _, v = design_R(m, h["psi" + L_], h["kappa" + L_], h["S" + L_]); pr.append(r5.resid(m, W) / np.sqrt(v))
        out.append(tuple(pr))
    return out


def variants(cl, hs, sel=None):
    if sel is not None: cl = [cl[i] for i in sel]; hs = [{k: v[sel] for k, v in h.items()} for h in hs]
    nd = np.array([len(s) for _, s in cl], float); Wd = np.column_stack([np.ones(len(cl)), np.log(r5.rhat(cl)), np.log(nd)]); Rs = r5.Rsplits(hs)
    return dict(XDSM=r5.xdsm_W(Rs, iu)[0], XDSMd=r5.xdsm_W(Rs, iu, Wd)[0], XDSMda=r5.xdsm_W(ca_splits(hs, Wd), iu)[0])


def summ(pv, null=None, true=None):
    o = {}
    for m_, v in pv.items():
        c = bh(np.asarray(v, float)); nm = np.ones(len(v), bool) if null is None else null
        o[m_] = dict(calls=int(c.sum()), V=int((c & nm).sum()), t05=float((np.asarray(v)[nm] < 0.05).mean()))
        if true is not None: o[m_].update(FDP=float((c & nm).sum() / max(c.sum(), 1)), power=float(c[true].mean()))
    return o


if a.mode in ("A", "O"):
    Xs = r5.centre(X, lib, don) if a.mode == "A" else X
    for k in range(a.k0, a.k1):
        fn = os.path.join(OD, f"c_{a.mode}_{a.sub}_k{k}.json")
        if os.path.exists(fn): continue
        t0 = time.time(); seed = (61_000_000 if a.mode == "A" else 65_000_000) + 10 * k + si; lab = np.random.default_rng(seed).permutation(don)
        pu, cl = r5.cells_of(Xs, lib, lab); fs, hs = stats_from_cells(cl, seed); res = {"all": summ(variants(cl, hs))}
        if a.mode == "A":
            order = np.random.default_rng(seed + 1).permutation(len(pu))
            for D in (100, 200): res[str(D)] = summ(variants(cl, hs, np.sort(order[:D])))
        res["D"] = int(len(pu)); res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w")); print("c", a.mode, a.sub, k, {d: {m: v["calls"] for m, v in r.items()} for d, r in res.items() if isinstance(r, dict)}, res["sec"], "s", flush=True)

elif a.mode == "S":
    pu0, cl0 = r5.cells_of(X, lib, don); st0 = [r5._ms(Y, s) for Y, s in cl0]; m0 = np.array([t[0] for t in st0]); p0 = np.array([t[1] for t in st0])
    tau2, mu0, _ = r5.paule_mandel(m0, p0); s_ = float(np.clip(np.median(np.sqrt(tau2[elig]) / np.maximum(mu0[elig], 1e-30)), 0.1, 0.5)); del cl0, st0
    RHO = (0.1, 0.2, 0.3, 0.4); mods = np.random.default_rng(72_000_000 + si).choice(elig, 40, replace=False).reshape(4, 10); gm = -np.ones(G, int)
    for i, mm in enumerate(mods): gm[mm] = i
    true = (gm[A_] >= 0) & (gm[A_] == gm[B_]); null = ~true; Xs = r5.centre(X, lib, don)
    for k in range(a.k0, a.k1):
        fn = os.path.join(OD, f"c_S_{a.sub}_k{k}.json")
        if os.path.exists(fn): continue
        t0 = time.time(); seed = 70_000_000 + 10 * k + si; rng = np.random.default_rng(seed); lab = rng.permutation(don); pu = np.unique(lab)
        Dk = len(pu); z = rng.standard_normal((Dk, 4)); e = rng.standard_normal((Dk, G)); U = s_ * e
        for i in range(4): U[:, mods[i]] = s_ * (np.sqrt(RHO[i]) * z[:, [i]] + np.sqrt(1 - RHO[i]) * e[:, mods[i]])
        F = np.exp(U - s_ ** 2 / 2); _, cl = r5.cells_of(Xs, lib, lab, units=pu); cl = [(Y * F[i][None, :], s) for i, (Y, s) in enumerate(cl)]
        fs, hs = stats_from_cells(cl, seed); res = summ(variants(cl, hs), null, true); res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w"))
        print("c S", a.sub, k, {m: (v["calls"], round(v["FDP"], 3), round(v["power"], 2)) for m, v in res.items() if isinstance(v, dict)}, res["sec"], "s", flush=True)
