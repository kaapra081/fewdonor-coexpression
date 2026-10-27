"""Frozen held-out evaluation (HELDOUT_SPEC.md). New panel (seed 4112026), new designs, new seeds.
python run_heldout.py --design HO_seaad24 --rep0 0 --rep1 200 --out runs_heldout"""
import os, sys, time, argparse, json
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1"); os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sim_v2 as S, methods_v2 as M  # noqa: E402
from run_dev import one_rep  # noqa: E402

PANEL_SEED = 4112026
HO = {
    "HO_seaad24":       dict(D=24, pool="seaad_v3", pool_n="real"),
    "HO_seaad24_cap40": dict(D=24, pool="seaad_v3", nd_median=40, nd_min=10),
    "HO_lee40":         dict(D=40, pool="lee80x1600", nd_median=60, nd_logsd=0.8, nd_min=10),
    "HO_D32":           dict(D=32, nd_median=40, nd_min=10),
    "HO_weak40":        dict(D=40, g=0.20, g_Q=0.20, nd_median=40, nd_min=10),
    "HO_D16_n100":      dict(D=16, nd_median=100, nd_min=10),
}
HO_INDEX = {k: i for i, k in enumerate(HO)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True); ap.add_argument("--rep0", type=int, default=0)
    ap.add_argument("--rep1", type=int, default=200); ap.add_argument("--out", default="runs_heldout")
    ap.add_argument("--seedbase", type=int, default=5_000_000)
    a = ap.parse_args()
    pan = S.make_panel(seed=PANEL_SEED); st = S.setting(**HO[a.design]); tr = S.truth(pan, st); iu = tr["iu"]
    od = os.path.join(a.out, a.design); os.makedirs(od, exist_ok=True)
    json.dump(dict(design=a.design, setting={k: (float(v) if isinstance(v, (int, float)) else v) for k, v in st.items()},
                   panel_seed=PANEL_SEED, seedbase=a.seedbase), open(os.path.join(od, "config.json"), "w"), indent=1)
    np.savez_compressed(os.path.join(od, "truth.npz"), sb=tr["sb"], rB=tr["rB"], cls=tr["cls"].astype(str),
                        altB=tr["altB"], sBdiag=np.diag(tr["SB"]))
    for rep in range(a.rep0, a.rep1):
        fn = os.path.join(od, f"rep{rep:04d}.npz")
        if os.path.exists(fn):
            continue
        seed = a.seedbase + 100_019 * HO_INDEX[a.design] + rep
        res = one_rep(pan, st, tr, seed)
        sim = S.simulate(pan, st, seed)              # identical replicate (same seed) for the intervals
        d = M.donor_stats(sim["Y"], sim["s_hat"], sim["donor_idx"])
        L = d["kappa"][:, None] * d["C"][:, iu[0], iu[1]]
        bt = M.boot_t_interval(d["m"][:, iu[0]], d["m"][:, iu[1]], L, seed=seed + 11)
        mt = sim["mu_true"]
        obt = M.boot_t_interval(mt[:, iu[0]], mt[:, iu[1]], np.zeros_like(L), seed=seed + 12)
        res.update(BT_lo=bt["lo"].astype(np.float32), BT_hi=bt["hi"].astype(np.float32),
                   OBT_lo=obt["lo"].astype(np.float32), OBT_hi=obt["hi"].astype(np.float32))
        np.savez_compressed(fn + ".tmp.npz", **res); os.replace(fn + ".tmp.npz", fn)
        print(a.design, rep, round(float(res["runtime_total_s"]), 1), flush=True)


if __name__ == "__main__":
    main()
