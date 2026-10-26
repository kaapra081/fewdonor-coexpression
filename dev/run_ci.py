"""Change log 1 (DEV_PREREG.md): amended interval runs. Re-simulates exactly the same replicates as run_dev.py
(same seeds) and computes the donor bootstrap-t interval (BOOT_T) for the leak-corrected equal-weight estimator,
plus the same interval on the TRUE donor means (ORACLE_BOOT_T, leak term zero).
python run_ci.py --design D16 --rep0 0 --rep1 200 --out runs_ci"""
import os, sys, time, argparse
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sim_v2 as S, methods_v2 as M  # noqa: E402
from run_dev import DESIGNS, DESIGN_INDEX  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True); ap.add_argument("--rep0", type=int, default=0)
    ap.add_argument("--rep1", type=int, default=200); ap.add_argument("--out", default="runs_ci")
    ap.add_argument("--seedbase", type=int, default=2_000_000)
    a = ap.parse_args()
    pan = S.make_panel(); st = S.setting(**DESIGNS[a.design]); tr = S.truth(pan, st); iu = tr["iu"]
    od = os.path.join(a.out, a.design); os.makedirs(od, exist_ok=True)
    np.savez_compressed(os.path.join(od, "truth.npz"), sb=tr["sb"], cls=tr["cls"].astype(str), altB=tr["altB"])
    for rep in range(a.rep0, a.rep1):
        fn = os.path.join(od, f"rep{rep:04d}.npz")
        if os.path.exists(fn):
            continue
        seed = a.seedbase + 100_003 * DESIGN_INDEX[a.design] + rep
        t0 = time.time()
        sim = S.simulate(pan, st, seed)
        d = M.donor_stats(sim["Y"], sim["s_hat"], sim["donor_idx"])
        L = d["kappa"][:, None] * d["C"][:, iu[0], iu[1]]
        o = M.boot_t_interval(d["m"][:, iu[0]], d["m"][:, iu[1]], L, seed=seed + 11)
        mt = sim["mu_true"]
        orc = M.boot_t_interval(mt[:, iu[0]], mt[:, iu[1]], np.zeros_like(L), seed=seed + 12)
        np.savez_compressed(fn + ".tmp.npz", BT_est=o["est"].astype(np.float32), BT_lo=o["lo"].astype(np.float32),
                            BT_hi=o["hi"].astype(np.float32), OBT_lo=orc["lo"].astype(np.float32),
                            OBT_hi=orc["hi"].astype(np.float32), runtime_s=np.float32(time.time() - t0))
        os.replace(fn + ".tmp.npz", fn)
        print(a.design, rep, round(time.time() - t0, 1), flush=True)


if __name__ == "__main__":
    main()
