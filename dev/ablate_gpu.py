"""POST HOC ablations (not pre-registered; requested by internal review). Same development replicates (same seeds) as
run_dev.py. Every variant uses B = 100,000 GPU permutations for all pairs, no refinement, so the variants are
directly comparable with each other (including 'LCP_gpu', the full method computed the same way).
Variants:
  LCP_gpu      full LCP statistic
  LCP_noleak   LCP without the leak term
  LCP_truev    LCP standardized with the TRUE total variances (Sigma_B,jj + psi_dj) and true means; estimated leak
  LCP_unstd    leak-corrected raw covariance (no standardization), donor permutation
  PB_perm      pseudobulk log-CPM donor means, standardized by sample SD, donor permutation (no leak correction)
python ablate_gpu.py --design D16 --rep0 0 --rep1 200 --out runs_ablate --device cuda:0"""
import os, sys, argparse, time
import numpy as np
try:
    import torch
except ImportError:  # CPU fallback (numpy) when torch is unavailable
    torch = None
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sim_v2 as S, methods_v2 as M  # noqa: E402
from run_dev import DESIGNS, DESIGN_INDEX  # noqa: E402


def perm_p_np(R, obs, iu, B, seed, batch=500):
    D, G = R.shape; a, b = iu
    ao = np.abs(obs) - 1e-6 * np.maximum(1.0, np.abs(obs))
    cnt = np.zeros(len(a), dtype=np.int64); rng = np.random.default_rng(int(seed)); done = 0
    R = R.astype(np.float64)
    while done < B:
        bb = min(batch, B - done)
        perms = np.argsort(rng.random((bb, D)), axis=1)
        T = np.matmul(R[perms].transpose(0, 2, 1), R)
        cnt += (np.abs(T[:, a, b]) >= ao).sum(0); done += bb
    return (cnt + 1) / (B + 1)


def perm_p(R, obs, iu, B, seed, device, batch=500):
    if device == "cpu" or torch is None:
        return perm_p_np(R, obs, iu, B, seed, batch)
    D, G = R.shape
    Rt = torch.tensor(R, dtype=torch.float32, device=device)
    ai = torch.tensor(iu[0], device=device); bi = torch.tensor(iu[1], device=device)
    ao = torch.tensor(np.abs(obs) - 1e-6 * np.maximum(1.0, np.abs(obs)), dtype=torch.float32, device=device)
    cnt = torch.zeros(len(iu[0]), dtype=torch.int64, device=device)
    g = torch.Generator(device=device); g.manual_seed(int(seed) % (2 ** 31))
    done = 0
    while done < B:
        bb = min(batch, B - done)
        perms = torch.argsort(torch.rand(bb, D, generator=g, device=device), dim=1)
        T = torch.bmm(Rt[perms].transpose(1, 2), Rt.expand(bb, D, G))
        cnt += (T[:, ai, bi].abs() >= ao).sum(0); done += bb
    return (cnt.cpu().numpy() + 1) / (B + 1)


def variants(d, sim, tr, iu):
    m, psi, C, kap = d["m"], d["psi"], d["C"], d["kappa"]; D = m.shape[0]; a, b = iu
    tau2, mu, _ = M.paule_mandel(m, psi); v = tau2[None, :] + psi
    out = {}
    def lcp_stat(R, vv, leak=True):
        Lt = kap[:, None] * C[:, a, b] / np.sqrt(vv[:, a] * vv[:, b]) if leak else np.zeros((D, len(a)))
        return ((R[:, a] * R[:, b]).sum(0) - Lt.sum(0)) * np.sqrt(D / (D + (Lt ** 2).sum(0)))
    R = (m - mu) / np.sqrt(v)
    out["LCP_gpu"] = (R, lcp_stat(R, v))
    out["LCP_noleak"] = (R, lcp_stat(R, v, leak=False))
    vt = np.diag(tr["SB"])[None, :] + sim["psi_true"]; mut = np.asarray(S.make_panel()["mu"])
    Rt = (m - mut) / np.sqrt(vt)
    out["LCP_truev"] = (Rt, lcp_stat(Rt, vt))
    # unweighted donor deviations; each gene rescaled by one constant (its SD across donors), which leaves the
    # permutation test unchanged but keeps the statistic away from float32 underflow
    sd = m.std(0, ddof=1); Rc = (m - m.mean(0)) / sd
    Lr = kap[:, None] * C[:, a, b] / (sd[a] * sd[b])
    out["LCP_unstd"] = (Rc, (Rc[:, a] * Rc[:, b]).sum(0) - Lr.sum(0) * (D - 1) / D)
    L = np.log1p(1e6 * m); Rl = (L - L.mean(0)) / L.std(0, ddof=1)
    out["PB_perm"] = (Rl, (Rl[:, a] * Rl[:, b]).sum(0))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True); ap.add_argument("--rep0", type=int, default=0)
    ap.add_argument("--rep1", type=int, default=200); ap.add_argument("--out", default="runs_ablate")
    ap.add_argument("--seedbase", type=int, default=2_000_000); ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--B", type=int, default=100_000)
    a = ap.parse_args()
    pan = S.make_panel(); st = S.setting(**DESIGNS[a.design]); tr = S.truth(pan, st); iu = tr["iu"]
    od = os.path.join(a.out, a.design); os.makedirs(od, exist_ok=True)
    np.savez_compressed(os.path.join(od, "truth.npz"), sb=tr["sb"], cls=tr["cls"].astype(str), altB=tr["altB"])
    for rep in range(a.rep0, a.rep1):
        fn = os.path.join(od, f"rep{rep:04d}.npz")
        if os.path.exists(fn):
            continue
        t0 = time.time()
        seed = a.seedbase + 100_003 * DESIGN_INDEX[a.design] + rep
        sim = S.simulate(pan, st, seed); d = M.donor_stats(sim["Y"], sim["s_hat"], sim["donor_idx"])
        res = {k: perm_p(R, obs, iu, a.B, seed + 5, a.device).astype(np.float32) for k, (R, obs) in variants(d, sim, tr, iu).items()}
        np.savez_compressed(fn, **res)
        print(a.design, rep, round(time.time() - t0, 1), flush=True)


if __name__ == "__main__":
    main()
