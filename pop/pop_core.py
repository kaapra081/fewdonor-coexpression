"""Shared computations for POP_PREREG.md (PB, LCP, LCPX, LCPX-S at population scale)."""
import os, sys, numpy as np, scipy.sparse as sp, torch
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "dev")); sys.path.insert(0, os.path.join(HERE, "..", "xsplit"))
from methods_v2 import paule_mandel, pearson_fisher, bh  # noqa: E402,F401
from lcpx import half_stats  # noqa: E402


def donor_level(X, lib, donor, split_seed):
    """Per donor: full-data m, psi, kappa, S and within-donor covariance C (G x G, float32); and half-split m, psi."""
    ud = np.unique(donor); rng = np.random.default_rng(split_seed)
    out = {k: [] for k in ("m", "psi", "kappa", "S", "C", "mA", "psiA", "mB", "psiB", "n")}
    idx = {u: np.where(donor == u)[0] for u in ud}
    for u in ud:
        r = idx[u]; Y = X[r].toarray().astype(np.float64); s = lib[r].astype(np.float64); n = len(r)
        S = s.sum(); x = Y / s[:, None]; m = Y.sum(0) / S; kap = (s ** 2).sum() / S ** 2; xc = x - x.mean(0)
        sig2 = (xc ** 2).sum(0) / (n - 1) - (x / s[:, None]).mean(0)
        out["m"].append(m); out["psi"].append(np.maximum(kap * np.maximum(sig2, 0) + m / S, 1e-30)); out["kappa"].append(kap); out["S"].append(S)
        out["C"].append((xc.T @ xc / (n - 1)).astype(np.float32)); out["n"].append(n)
        p = rng.permutation(n); A, B = p[: n // 2], p[n // 2:]
        for h, rr in (("A", A), ("B", B)):
            mm, ps = half_stats(Y[rr], s[rr]); out["m" + h].append(mm); out["psi" + h].append(ps)
    d = {k: np.array(v) for k, v in out.items()}; d["donors"] = ud
    return d


def subset(d, rows):
    return {k: (v[rows] if isinstance(v, np.ndarray) and v.shape[:1] == (len(d["donors"]),) else v) for k, v in d.items()}


def std_R(m, psi):
    tau2, mu, _ = paule_mandel(m, psi); return (m - mu[None, :]) / np.sqrt(tau2[None, :] + psi)


def perm_counts(L, Rr, a, b, obs, stud_obs=None, B=100_000, seed=0, batch=200):
    """Two-sided donor-permutation p-values for sum_d L[pi(d),a] Rr[d,b] (and its studentized version)."""
    batch = int(os.environ.get("PERM_BATCH", batch))   # REV4: larger GPU batches; default unchanged
    dv = torch.device(os.environ.get("PERM_DEVICE", "cpu") if torch.cuda.is_available() or os.environ.get("PERM_DEVICE", "cpu") == "cpu" else "cpu")   # optional GPU; permutations still drawn on the CPU generator
    D = L.shape[0]; Lt = torch.tensor(L, dtype=torch.float32, device=dv); Rt = torch.tensor(Rr, dtype=torch.float32, device=dv)
    ai = torch.tensor(a, device=dv); bi = torch.tensor(b, device=dv); o = torch.tensor(np.abs(obs) - 1e-5, dtype=torch.float32, device=dv)
    c = torch.zeros(len(a), dtype=torch.int64, device=dv); cs = torch.zeros(len(a), dtype=torch.int64, device=dv)
    if stud_obs is not None: os_ = torch.tensor(np.abs(stud_obs) - 1e-5, dtype=torch.float32, device=dv); L2 = Lt ** 2; R2t = (Rt ** 2).T.unsqueeze(0)
    RtT = Rt.T.unsqueeze(0); g = torch.Generator(); g.manual_seed(int(seed) % (2 ** 31)); done = 0
    while done < B:
        bb = min(batch, B - done); perms = torch.argsort(torch.rand(bb, D, generator=g), dim=1).to(dv)
        S1 = torch.matmul(RtT, Lt[perms])                      # [bb, G, G]: entry (k, j) = sum_d Rr[d,k] L[pi(d), j]
        S1ab = S1[:, bi, ai]; c += (S1ab.abs() >= o).sum(0)
        if stud_obs is not None:
            S2 = torch.matmul(R2t, L2[perms])[:, bi, ai]
            var = (S2 - S1ab ** 2 / D) / (D - 1); t = (S1ab / D) / torch.sqrt(torch.clamp(var, min=1e-30)) * np.sqrt(D)
            cs += (t.abs() >= os_).sum(0)
        done += bb
    p = (c.cpu().numpy() + 1) / (B + 1)
    return (p, (cs.cpu().numpy() + 1) / (B + 1)) if stud_obs is not None else p


def normal_perm_p(L, Rr, a, b, obs):
    """Normal approximation with the exact permutation mean and variance of sum_d L[pi(d),a] Rr[d,b] (Hoeffding)."""
    D = L.shape[0]; La, Rb = L[:, a], Rr[:, b]
    mu = La.sum(0) * Rb.sum(0) / D
    var = ((La - La.mean(0)) ** 2).sum(0) * ((Rb - Rb.mean(0)) ** 2).sum(0) / (D - 1)
    sd = np.sqrt(np.maximum(var, 1e-300))
    return norm.sf((np.abs(obs) - mu) / sd) + norm.cdf((-np.abs(obs) - mu) / sd)


def run_methods(d, iu, B=100_000, seed=0, chunk=20000, skip_x=False):
    a, b = iu; D = len(d["donors"]); out = {}
    out["PB"] = pearson_fisher(np.log1p(1e6 * d["m"]), iu)["p"]
    R = std_R(d["m"], d["psi"]); tau2, mu, _ = paule_mandel(d["m"], d["psi"]); v = tau2[None, :] + d["psi"]
    obs = np.empty(len(a)); sig2 = np.maximum((d["psi"] - d["m"] / d["S"][:, None]) / d["kappa"][:, None], 0).mean(0)
    r_w = np.empty(len(a))
    for c0 in range(0, len(a), chunk):
        aa, bb = a[c0:c0 + chunk], b[c0:c0 + chunk]
        Cab = d["C"][:, aa, bb].astype(np.float64)
        Lt = d["kappa"][:, None] * Cab / np.sqrt(v[:, aa] * v[:, bb])
        num = (R[:, aa] * R[:, bb]).sum(0) - Lt.sum(0)
        obs[c0:c0 + chunk] = num * np.sqrt(D / (D + (Lt ** 2).sum(0)))
        with np.errstate(invalid="ignore", divide="ignore"):
            r_w[c0:c0 + chunk] = Cab.mean(0) / np.sqrt(sig2[aa] * sig2[bb])
    RA = std_R(d["mA"], d["psiA"]); RB = std_R(d["mB"], d["psiB"]); q = RA[:, a] * RB[:, b]
    T = q.sum(0); t = q.mean(0) / np.maximum(q.std(0, ddof=1), 1e-300) * np.sqrt(D)
    if D <= 100:
        out["LCP"] = perm_counts(R, R, a, b, obs, B=B, seed=seed)
        if skip_x: out["LCPX"] = out["LCPXS"] = None   # REV4 R1: LCPX not reported; skipped to save permutations
        else: out["LCPX"], out["LCPXS"] = perm_counts(RA, RB, a, b, T, stud_obs=t, B=B, seed=seed + 1)
    else:
        out["LCP"] = normal_perm_p(R, R, a, b, obs)
        out["LCPX"] = normal_perm_p(RA, RB, a, b, T)
        out["LCPXS"] = 2 * norm.sf(np.abs(t))
    out["sign_PB"] = np.sign(pearson_fisher(np.log1p(1e6 * d["m"]), iu)["est"]); out["sign_LCP"] = np.sign(obs); out["sign_LCPX"] = np.sign(T)
    out["r_w"] = r_w
    return out
