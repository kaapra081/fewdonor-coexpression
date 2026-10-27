"""Run LCP and comparators on real per-donor statistics.
python real_run.py --stats real_stats_seaad.npz --out runs/seaad_full.npz [--subset-D 24 --subset-seed 3]"""
import os, sys, argparse, time
os.environ.setdefault("OMP_NUM_THREADS", "2"); os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "dev"))
import methods_v2 as M  # noqa: E402


def load(path, rows=None):
    z = np.load(path, allow_pickle=True)
    d = {k: z[k] for k in z.files}
    if rows is not None:
        for k in ("m", "psi", "C", "kappa", "S", "n", "donors", "sex"):
            d[k] = d[k][rows]
    d["C"] = d["C"].astype(np.float64)
    return d


def run(d, seed):
    G = d["m"].shape[1]; iu = np.triu_indices(G, 1); a, b = iu
    t0 = time.time()
    o = M.lcsp_test(d, iu, seed=seed, studentize=False, batch=40)
    out = dict(LCP=o["p"], LCP_stat=o["stat"], LCP_est=o["est"])
    cm = M.cmb_mom(d, iu); out.update(CMB_mom=cm["p"], BCES_est=cm["est"], BCES_se=cm["se"])
    pb = M.pb_logcpm(d, iu); out.update(PB_logcpm=pb["p"], PB_r=pb["est"])
    out["CMB_naive"] = M.pearson_fisher(d["m"], iu)["p"]
    D = d["m"].shape[0]
    mc = d["m"] - d["m"].mean(0)
    Sjk = (mc[:, a] * mc[:, b]).sum(0) / (D - 1)
    leak = (d["kappa"][:, None] * d["C"][:, a, b]).mean(0)
    out.update(S_pb=Sjk, leak=leak)
    # within-donor (count scale) equal-weight covariance and correlation (Poisson-corrected variances)
    Cw = d["C"][:, a, b].mean(0)
    var_lat = np.maximum(np.einsum("dgg->dg", d["C"]) - 0, 0).mean(0)   # sample variance of x (includes Poisson)
    out["W_cov"] = Cw
    out["runtime_s"] = np.float32(time.time() - t0)
    return {k: np.asarray(v, np.float32) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--subset-D", type=int, default=0); ap.add_argument("--subset-seed", type=int, default=0)
    a = ap.parse_args()
    z = np.load(a.stats, allow_pickle=True); Dfull = len(z["donors"])
    rows = None
    if a.subset_D:
        rows = np.sort(np.random.default_rng(a.subset_seed).choice(Dfull, size=a.subset_D, replace=False))
    d = load(a.stats, rows)
    res = run(d, seed=97 + a.subset_seed)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    np.savez_compressed(a.out, **res, rows=(rows if rows is not None else np.arange(Dfull)), gene_sym=d["gene_sym"],
                        gene_ids=d["gene_ids"])
    print(a.out, "runtime", float(res["runtime_s"]), "LCP BH", int(M.bh(res["LCP"]).sum()), "PB BH", int(M.bh(res["PB_logcpm"]).sum()), flush=True)


if __name__ == "__main__":
    main()
