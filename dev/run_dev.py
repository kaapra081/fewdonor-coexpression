"""Run development replicates for one design. Usage:
    python run_dev.py --design D16 --rep0 0 --rep1 200 --out runs_dev [--seedbase 2000000]
Saves one compressed npz per replicate with p-values (float32) for every method and interval bounds.
Designs and seeds are fixed in DEV_PREREG.md (hashed before the formal run)."""
import os, sys, time, argparse, json
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1"); os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sim_v2 as S, methods_v2 as M  # noqa: E402

LN40 = dict(nd_median=40, nd_min=10)
DESIGNS = {
    "D8":  dict(D=8, **LN40),
    "D12": dict(D=12, **LN40),
    "D16": dict(D=16, **LN40),
    "D24": dict(D=24, **LN40),
    "D40": dict(D=40, **LN40),
    "D16_strong": dict(D=16, g=0.45, **LN40),              # secondary: stronger between programs
    "D16_STRESS_depth": dict(D=16, kappa=0.3, se_mult=2.0, **LN40),   # STRESS: depth-state coupling
}
DESIGN_INDEX = {k: i for i, k in enumerate(DESIGNS)}


def one_rep(pan, st, tr, seed):
    iu = tr["iu"]
    t0 = time.time()
    sim = S.simulate(pan, st, seed)
    d = M.donor_stats(sim["Y"], sim["s_hat"], sim["donor_idx"])
    out = {}
    o = M.lcsp_test(d, iu, seed=seed, studentize=False); out["LCP"] = o["p"]; out["LCP_est"] = o["est"]
    t_lcp = time.time() - t0
    out["LCP_stud"] = M.lcsp_test(d, iu, seed=seed + 1, studentize=True)["p"]
    ci = M.eq_interval(d, iu); out["EQ_est"] = ci["est"]; out["EQ_lo"] = ci["lo"]; out["EQ_hi"] = ci["hi"]
    cm = M.cmb_mom(d, iu); out["CMB_mom"] = cm["p"]
    K = d["m"].shape[0]; tq = __import__("scipy").stats.t.ppf(0.975, K - 1)
    out["CMB_lo"] = cm["est"] - tq * cm["se"]; out["CMB_hi"] = cm["est"] + tq * cm["se"]
    out["CMB_naive"] = M.pearson_fisher(d["m"], iu)["p"]
    out["PB_logcpm"] = M.pb_logcpm(d, iu)["p"]
    out["E5"] = M.e5_logstd(sim["Y"], sim["s_hat"], sim["donor_idx"], iu)["p"]
    out["ORACLE"] = M.oracle_perm(sim, iu, seed + 2, studentize=False)["p"]
    out["BOUND"] = M.bound_test(sim, tr, iu, seed + 3, studentize=False)["p"]
    out = {k: np.asarray(v, np.float32) for k, v in out.items()}
    out["runtime_LCP_s"] = np.float32(t_lcp); out["runtime_total_s"] = np.float32(time.time() - t0)
    out["nd"] = sim["nd"].astype(np.int32)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True); ap.add_argument("--rep0", type=int, default=0)
    ap.add_argument("--rep1", type=int, default=200); ap.add_argument("--out", default="runs_dev")
    ap.add_argument("--seedbase", type=int, default=2_000_000)
    a = ap.parse_args()
    pan = S.make_panel(); st = S.setting(**DESIGNS[a.design]); tr = S.truth(pan, st)
    od = os.path.join(a.out, a.design); os.makedirs(od, exist_ok=True)
    json.dump(dict(design=a.design, setting={k: (float(v) if isinstance(v, (int, float)) else v) for k, v in st.items()},
                   seedbase=a.seedbase), open(os.path.join(od, "config.json"), "w"), indent=1)
    np.savez_compressed(os.path.join(od, "truth.npz"), sb=tr["sb"], rB=tr["rB"], cls=tr["cls"].astype(str),
                        altB=tr["altB"], sBdiag=np.diag(tr["SB"]))
    for rep in range(a.rep0, a.rep1):
        fn = os.path.join(od, f"rep{rep:04d}.npz")
        if os.path.exists(fn):
            continue
        seed = a.seedbase + 100_003 * DESIGN_INDEX[a.design] + rep
        res = one_rep(pan, st, tr, seed)
        np.savez_compressed(fn + ".tmp.npz", **res); os.replace(fn + ".tmp.npz", fn)
        print(a.design, rep, round(float(res["runtime_total_s"]), 1), flush=True)


if __name__ == "__main__":
    main()
