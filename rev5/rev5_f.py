"""Plan 8 (post hoc, rev5/PLAN8.md, sha256 in PLAN8.sha256).
PsychAD:  python rev5_f.py --mode B|S --sub P1|P2|P3 [--k0 --k1]   B: brain-bank-stratified null (Plan 5 R2 seeds), S: planted signal (Plan 6d S seeds)
Lupus:    python rev5_f.py --mode LA|LD|LR --sub LM|FC [--k0 --k1] LA: donor-centred null, LD: library-preserving null (Plan 7 seeds), LR: real donors
Outputs $REV5_OUT/f_<mode>_<sub>[_k<k>].json. XDS-Mh is the Plan 6d definition (rev5_d.variants, copied verbatim below)."""
import os, sys, json, argparse, time, types, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(); ap.add_argument("--mode", required=True); ap.add_argument("--sub", required=True)
ap.add_argument("--k0", type=int, default=0); ap.add_argument("--k1", type=int, default=20); a = ap.parse_args()
OD = os.environ.get("REV5_OUT", os.path.join(HERE, "out")); os.makedirs(OD, exist_ok=True)


def resid(Z, W): return Z - W @ np.linalg.lstsq(W, Z, rcond=None)[0]


def half_W(nd, h, extra=None):
    """Plan 6d half-design columns [1, log(kappa^h S^h), log n^h] (+ optional donor-level columns, XDS-Mhb)."""
    nA = np.floor(nd / 2); nB = nd - nA
    WA = np.column_stack([np.ones(len(nd)), np.log(h["kappaA"] * h["SA"]), np.log(nA)]); WB = np.column_stack([np.ones(len(nd)), np.log(h["kappaB"] * h["SB"]), np.log(nB)])
    if extra is not None: WA = np.column_stack([WA, extra]); WB = np.column_stack([WB, extra])
    return WA, WB


def bh(p, q=0.05):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p); thr = q * np.arange(1, m + 1) / m; ok = np.nonzero(p[o] <= thr)[0]
    c = np.zeros(m, bool)
    if len(ok): c[o[:ok.max() + 1]] = True
    return c


def matched(p, true, levels=(0.05, 0.10)):
    """Oracle-matched power: longest prefix of the p-value ordering with false-discovery proportion <= L."""
    o = np.argsort(np.asarray(p, float), kind="stable"); t = true[o]; tp = np.cumsum(t); n = np.arange(1, len(t) + 1); fdp = (n - tp) / n; out = {}
    for L in levels:
        ok = np.nonzero(fdp <= L)[0]; out[f"power_at_fdp{L:.2f}"] = float(tp[ok.max()] / t.sum()) if len(ok) else 0.0
    return out


def summ(pv, null=None, true=None):
    o = {}
    for m_, v in pv.items():
        v = np.asarray(v, float); c = bh(v); nm = np.ones(len(v), bool) if null is None else null
        o[m_] = dict(calls=int(c.sum()), V=int((c & nm).sum()), t05=float((v[nm] < 0.05).mean()))
        if true is not None: o[m_].update(FDP=float((c & nm).sum() / max(c.sum(), 1)), power=float(c[true].mean()), **matched(v, true))
    return o


if a.mode in ("B", "S"):
    # rev5.py executed as a module with its per-cell pool line disabled (pool_P*.npz is not used by modes B or S)
    src = open(os.path.join(HERE, "rev5.py")).read(); line = [l for l in src.split("\n") if l.startswith("pool = np.load(")]; assert len(line) == 1
    src = src.replace(line[0], "pool = None"); sys.argv = [os.path.join(HERE, "rev5.py"), "--mode", "none", "--sub", a.sub]
    r5 = types.ModuleType("r5"); r5.__file__ = os.path.join(HERE, "rev5.py"); exec(compile(src, r5.__file__, "exec"), r5.__dict__)
    from revm import stats_from_cells  # noqa: E402  (path set by rev5.py)
    bh = r5.bh
    X, lib, don, G, iu, A_, B_, si, elig = r5.X, r5.lib, r5.don, r5.G, r5.iu, r5.A_, r5.B_, r5.si, r5.elig

    def variants(cl, hs, extra=None):   # rev5_d.variants (XDSM, XDSMh) + XDSMhb when `extra` (bank indicators) is given
        nd = np.array([len(s) for _, s in cl], float); Rs = r5.Rsplits(hs); Rh, Rb = [], []
        for (RA, RB), h in zip(Rs, hs):
            WA, WB = half_W(nd, h); Rh.append((r5.resid(RA, WA), r5.resid(RB, WB)))
            if extra is not None: WA, WB = half_W(nd, h, extra); Rb.append((r5.resid(RA, WA), r5.resid(RB, WB)))
        out = dict(XDSM=r5.xdsm_W(Rs, iu)[0], XDSMh=r5.xdsm_W(Rh, iu)[0])
        if extra is not None: out["XDSMhb"] = r5.xdsm_W(Rb, iu)[0]
        return out

    if a.mode == "B":   # Plan 5 R2 reassignments: donor-centred cells permuted among donors of the same PsychAD Source bank
        srcb = json.load(open(os.path.join(HERE, "..", "rev4", "donor_source.json"))); bank = np.array([srcb[d_] for d_ in don]); Xs = r5.centre(X, lib, don)
        for k in range(a.k0, a.k1):
            fn = os.path.join(OD, f"f_B_{a.sub}_k{k}.json")
            if os.path.exists(fn): continue
            t0 = time.time(); seed = 64_000_000 + 10 * k + si; rng = np.random.default_rng(seed); lab = don.copy()
            for bk in ("H", "M", "R"):
                ix = np.where(bank == bk)[0]; lab[ix] = rng.permutation(don[ix])
            pu = np.unique(lab); prow = [np.where(lab == u)[0] for u in pu]; cl = [(Xs[r].toarray(), lib[r].astype(np.float64)) for r in prow]
            fs, hs = stats_from_cells(cl, seed); bk_u = np.array([srcb[u] for u in pu]); extra = np.column_stack([bk_u == "H", bk_u == "R"]).astype(float)
            res = dict(all=summ(variants(cl, hs, extra)), D=int(len(pu)), sec=round(time.time() - t0)); json.dump(res, open(fn, "w"))
            print("f B", a.sub, k, {m: v["calls"] for m, v in res["all"].items()}, res["sec"], "s", flush=True)

    else:   # Plan 6d S: planted modules on donor-centred cells (code of rev5_d mode S), plus PB and KX from revm.methods
        pu0, cl0 = r5.cells_of(X, lib, don); st0 = [r5._ms(Y, s) for Y, s in cl0]; m0 = np.array([t[0] for t in st0]); p0 = np.array([t[1] for t in st0])
        tau2, mu0, _ = r5.paule_mandel(m0, p0); s_ = float(np.clip(np.median(np.sqrt(tau2[elig]) / np.maximum(mu0[elig], 1e-30)), 0.1, 0.5)); del cl0, st0
        RHO = (0.1, 0.2, 0.3, 0.4); mods = np.random.default_rng(72_000_000 + si).choice(elig, 40, replace=False).reshape(4, 10); gm = -np.ones(G, int)
        for i, mm in enumerate(mods): gm[mm] = i
        true = (gm[A_] >= 0) & (gm[A_] == gm[B_]); null = ~true; Xs = r5.centre(X, lib, don)
        for k in range(a.k0, a.k1):
            fn = os.path.join(OD, f"f_S_{a.sub}_k{k}.json")
            if os.path.exists(fn): continue
            t0 = time.time(); seed = 70_000_000 + 10 * k + si; rng = np.random.default_rng(seed); lab = rng.permutation(don); pu = np.unique(lab)
            Dk = len(pu); z = rng.standard_normal((Dk, 4)); e = rng.standard_normal((Dk, G)); U = s_ * e
            for i in range(4): U[:, mods[i]] = s_ * (np.sqrt(RHO[i]) * z[:, [i]] + np.sqrt(1 - RHO[i]) * e[:, mods[i]])
            F = np.exp(U - s_ ** 2 / 2); _, cl = r5.cells_of(Xs, lib, lab, units=pu); cl = [(Y * F[i][None, :], s) for i, (Y, s) in enumerate(cl)]
            fs, hs = stats_from_cells(cl, seed); pv = variants(cl, hs); o = r5.rmethods(fs, hs, iu); pv["PB"] = o["PB"]; pv["KX"] = o["KX"]
            res = summ(pv, null, true); res["n_true"] = int(true.sum()); res["sec"] = round(time.time() - t0); json.dump(res, open(fn, "w"))
            print("f S", a.sub, k, {m: (v["calls"], round(v["FDP"], 3), round(v["power"], 2), round(v["power_at_fdp0.05"], 2)) for m, v in res.items() if isinstance(v, dict)}, res["sec"], "s", flush=True)

else:   # lupus (Plan 7 data and seeds); lupus_run.py imported as a module (its main block is guarded)
    import importlib.util
    spec = importlib.util.spec_from_file_location("lr", os.path.join(HERE, "..", "lupus", "lupus_run.py")); lr = importlib.util.module_from_spec(spec); spec.loader.exec_module(lr)
    bh = lr.bh
    si = {"LM": 0, "FC": 1}[a.sub]; X, lib, don, library, sym, disease = lr.load(a.sub); G = X.shape[1]; iu = np.triu_indices(G, 1); A_, B_ = iu

    def lvariants(cl, hs):
        nd = np.array([len(s) for _, s in cl], float); Rs = lr.Rsplits(hs); Rh = []
        for (RA, RB), h in zip(Rs, hs):
            WA, WB = half_W(nd, h); Rh.append((resid(RA, WA), resid(RB, WB)))
        pX, mX = lr.xdsm_W(Rs, iu); pH, mH = lr.xdsm_W(Rh, iu)
        return dict(XDSM=pX, XDSMh=pH), dict(XDSM=mX, XDSMh=mH)

    if a.mode in ("LA", "LD"):
        Xs = lr.centre(X, lib, don, exact=(a.mode == "LD"))
        for k in range(a.k0, a.k1):
            fn = os.path.join(OD, f"f_{a.mode}_{a.sub}_k{k}.json")
            if os.path.exists(fn): continue
            t0 = time.time(); seed = {"LA": 80_000_000, "LD": 82_000_000}[a.mode] + 10 * k + si; rng = np.random.default_rng(seed)
            if a.mode == "LD":
                lab = don.copy()
                for p_ in np.unique(library):
                    ix = np.where(library == p_)[0]; lab[ix] = rng.permutation(don[ix])
            else:
                lab = rng.permutation(don)
            pu = np.unique(lab); rows = [np.where(lab == u)[0] for u in pu]; cl = [(Xs[r].toarray().astype(np.float64), lib[r]) for r in rows]
            fs, hs = lr.stats_from_cells(cl, seed); pv, _ = lvariants(cl, hs)
            res = dict(all=summ(pv), D=int(len(pu)), sec=round(time.time() - t0)); json.dump(res, open(fn, "w"))
            print("f", a.mode, a.sub, k, {m: v["calls"] for m, v in res["all"].items()}, res["sec"], "s", flush=True)
    elif a.mode == "LR":
        fn = os.path.join(OD, f"f_LR_{a.sub}.json")
        if not os.path.exists(fn):
            t0 = time.time(); pu, cl = lr.cells_of(X, lib, don); fs, hs = lr.stats_from_cells(cl, 44_000_001); pv, mq = lvariants(cl, hs)
            cX = bh(pv["XDSM"]); cH = bh(pv["XDSMh"]); isg = np.isin(sym, lr.ISG); ii = isg[A_] & isg[B_]
            res = dict(D=int(len(pu)), XDSM=int(cX.sum()), XDSMh=int(cH.sum()), keep=float((cX & cH).sum() / max(cX.sum(), 1)), new=int((cH & ~cX).sum()),
                       isg_pairs=int(ii.sum()), isg_XDSM=int((cX & ii).sum()), isg_XDSMh=int((cH & ii).sum()), isg_XDSMh_pos=int((cH & ii & (mq["XDSMh"] > 0)).sum()), sec=round(time.time() - t0))
            json.dump(res, open(fn, "w")); print("f LR", a.sub, res, flush=True)
