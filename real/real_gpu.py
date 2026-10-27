"""Real-data analyses on GPU (Kaggle T4). LCP p-values with B donor permutations for ALL pairs (no refinement step),
plus CPU comparators. Writes compact summaries only.
python real_gpu.py --stats-dir . --out real_out --device cuda:0 --part full|cap|sub24|sub40
"""
import os, sys, json, argparse, time
import numpy as np
import torch
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "dev"))
from methods_v2 import paule_mandel, bh, cmb_mom, pb_logcpm, pearson_fisher  # noqa: E402


def load(path, rows=None):
    z = np.load(path, allow_pickle=True)
    d = {k: z[k] for k in z.files}
    if rows is not None:
        for k in ("m", "psi", "C", "kappa", "S", "n", "donors", "sex"):
            d[k] = d[k][rows]
    d["C"] = d["C"].astype(np.float64)
    return d


def lcp_gpu(d, B, seed, device, batch=400):
    m, psi, C, kap = d["m"], d["psi"], d["C"], d["kappa"]
    D, G = m.shape
    tau2, mu, _ = paule_mandel(m, psi)
    v = tau2[None, :] + psi
    R = (m - mu[None, :]) / np.sqrt(v)
    iu = np.triu_indices(G, 1); a, b = iu
    Lt = kap[:, None] * C[:, a, b] / np.sqrt(v[:, a] * v[:, b])
    obs = ((R[:, a] * R[:, b]).sum(0) - Lt.sum(0)) * np.sqrt(D / (D + (Lt ** 2).sum(0)))
    Rt = torch.tensor(R, dtype=torch.float32, device=device)
    ai = torch.tensor(a, device=device); bi = torch.tensor(b, device=device)
    aobs = torch.tensor(np.abs(obs) - 1e-6, dtype=torch.float32, device=device)
    cnt = torch.zeros(len(a), dtype=torch.int64, device=device)
    g = torch.Generator(device=device); g.manual_seed(int(seed))
    done = 0
    while done < B:
        bb = min(batch, B - done)
        perms = torch.argsort(torch.rand(bb, D, generator=g, device=device), dim=1)
        Rp = Rt[perms]                                   # bb x D x G
        T = torch.bmm(Rp.transpose(1, 2), Rt.expand(bb, D, G))   # bb x G x G
        cnt += (T[:, ai, bi].abs() >= aobs).sum(0)
        done += bb
    p = (cnt.cpu().numpy() + 1) / (B + 1)
    w = 1.0 / np.sqrt(v[:, a] * v[:, b])
    est = ((R[:, a] * R[:, b]) - Lt).sum(0) / w.sum(0)
    return dict(p=p, T=obs, est=est, tau2=tau2, iu=iu)


def comparators(d, iu):
    a, b = iu; D = d["m"].shape[0]
    cm = cmb_mom(d, iu)
    pb = pb_logcpm(d, iu)
    nv = pearson_fisher(d["m"], iu)
    mc = d["m"] - d["m"].mean(0)
    S_pb = (mc[:, a] * mc[:, b]).sum(0) / (D - 1)
    leak = (d["kappa"][:, None] * d["C"][:, a, b]).mean(0)
    # within-donor latent correlation (equal donor weights): Poisson-corrected variances from psi
    sig2 = np.maximum((d["psi"] - d["m"] / d["S"][:, None]) / d["kappa"][:, None], 0).mean(0)
    Cw = d["C"][:, a, b].mean(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        r_w = Cw / np.sqrt(sig2[a] * sig2[b])
    return dict(p_MC=cm["p"], est_MC=cm["est"], se_MC=cm["se"], p_PB=pb["p"], r_PB=pb["est"], p_naive=nv["p"],
                S_pb=S_pb, leak=leak, r_w=r_w)


def rho_b(est, tau2, iu):
    a, b = iu
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where((tau2[a] > 0) & (tau2[b] > 0), est / np.sqrt(tau2[a] * tau2[b]), np.nan)


def run_one(path, rows, B, seed, device):
    d = load(path, rows)
    t0 = time.time()
    o = lcp_gpu(d, B, seed, device)
    cmp_ = comparators(d, o["iu"])
    res = dict(p_LCP=o["p"], T_LCP=o["T"], est_LCP=o["est"], **cmp_)
    est_bces = cmp_["est_MC"]
    res["rhoB"] = rho_b(est_bces, o["tau2"], o["iu"])
    res["runtime_s"] = time.time() - t0
    return res, d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats-dir", default="."); ap.add_argument("--out", default="real_out")
    ap.add_argument("--device", default="cuda:0"); ap.add_argument("--part", required=True)
    ap.add_argument("--nsub", type=int, default=20); ap.add_argument("--cohort", default="seaad,lee")
    ap.add_argument("--B-override", type=int, default=0)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    jobs = []
    for coh in a.cohort.split(","):
        f = os.path.join(a.stats_dir, f"real_stats_{coh}.npz")
        Dfull = len(np.load(f, allow_pickle=True)["donors"])
        if a.part == "full":
            jobs.append((f"{coh}_full", f, None, 400_000, 11))
        elif a.part == "cap":
            for r in range(3):
                jobs.append((f"{coh}_cap40_r{r}", os.path.join(a.stats_dir, f"real_stats_{coh}_cap40_r{r}.npz"), None, 100_000, 21 + r))
        elif a.part in ("sub24", "sub40"):
            Dsub = int(a.part[3:])
            for s in range(a.nsub):
                rows = np.sort(np.random.default_rng(1000 * Dsub + s).choice(Dfull, size=Dsub, replace=False))
                jobs.append((f"{coh}_{a.part}_s{s:02d}", f, rows, 100_000, 31 + s))
    for name, f, rows, B, seed in jobs:
        if a.B_override:
            B = a.B_override
        fn = os.path.join(a.out, name + ".npz")
        if os.path.exists(fn):
            continue
        res, d = run_one(f, rows, B, seed, a.device)
        np.savez_compressed(fn, **{k: np.asarray(v, np.float32) for k, v in res.items()},
                            rows=(rows if rows is not None else np.arange(len(d["donors"]))), B=B)
        print(name, "B", B, "runtime", round(res["runtime_s"], 1), "LCP", int(bh(res["p_LCP"]).sum()),
              "PB", int(bh(res["p_PB"]).sum()), "MC", int(bh(res["p_MC"]).sum()), flush=True)


if __name__ == "__main__":
    main()
