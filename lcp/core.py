"""Core implementation (numpy). Notation follows the paper.

Input: a cells x genes UMI count matrix Y (dense or scipy.sparse), per-cell library sizes s (total UMI over all
genes), and a donor label per cell. Output: per-pair LCP p-values for H0 'donor effects of genes j and k are
independent' (implies Sigma_B,jk = 0), the unbiased leak-corrected (BCES-type) estimate of Sigma_B,jk and the leak.

Validity: calibrated in simulation from ~24 donors; NOT calibrated in the BH-relevant tail at <= 16 donors
(see paper, Section 6). Intervals are not provided because covariance intervals undercover at these donor numbers.
"""
import numpy as np
import scipy.sparse as sp
from scipy import stats


def donor_statistics(Y, s, donors, min_cells=4):
    """Per-donor count-scale statistics.
    Returns dict with m (D x G ratio means), psi (sampling variances), C (D x G x G within-donor covariance of
    x = Y/s), kappa (D,), S (D,), n (D,), donors (D,)."""
    donors = np.asarray(donors).astype(str)
    s = np.asarray(s, float)
    out = dict(m=[], psi=[], C=[], kappa=[], S=[], n=[], donors=[])
    for d in np.unique(donors):
        r = np.where(donors == d)[0]
        if len(r) < min_cells:
            continue
        Yd = Y[r].toarray() if sp.issparse(Y) else np.asarray(Y[r])
        Yd = Yd.astype(float); sd = s[r]; S = sd.sum(); n = len(r)
        x = Yd / sd[:, None]; m = Yd.sum(0) / S; kap = (sd ** 2).sum() / S ** 2
        xc = x - x.mean(0)
        sig2 = (xc ** 2).sum(0) / (n - 1) - (x / sd[:, None]).mean(0)       # Poisson-corrected latent variance
        out["m"].append(m); out["psi"].append(np.maximum(kap * np.maximum(sig2, 0) + m / S, 1e-30))
        out["C"].append(xc.T @ xc / (n - 1)); out["kappa"].append(kap); out["S"].append(S); out["n"].append(n)
        out["donors"].append(d)
    return {k: np.array(v) for k, v in out.items()}


def paule_mandel(m, psi, iters=60):
    """Per-gene Paule-Mandel between-donor variance tau^2 and mean, given sampling variances psi (D x G)."""
    D = m.shape[0]
    lo = np.zeros(m.shape[1]); hi = np.maximum(m.var(0, ddof=1) * 10, 1e-30)

    def Q(t):
        w = 1.0 / (t[None, :] + psi); mu = (w * m).sum(0) / w.sum(0)
        return (w * (m - mu) ** 2).sum(0)
    q0 = Q(lo)
    for _ in range(iters):
        mid = 0.5 * (lo + hi); big = Q(mid) > D - 1
        lo = np.where(big, mid, lo); hi = np.where(big, hi, mid)
    tau2 = np.where(q0 <= D - 1, 0.0, 0.5 * (lo + hi))
    w = 1.0 / (tau2[None, :] + psi)
    return tau2, (w * m).sum(0) / w.sum(0)


def bh(p, q=0.05):
    p = np.asarray(p, float); ok = np.isfinite(p); rej = np.zeros(len(p), bool)
    if ok.sum() == 0:
        return rej
    ps = np.sort(p[ok]); k = np.where(ps <= q * np.arange(1, len(ps) + 1) / len(ps))[0]
    if len(k):
        rej[ok] = p[ok] <= ps[k.max()]
    return rej


def bces_estimate(stats_, pairs=None):
    """Unbiased leak-corrected estimate S_B,jk - mean_d kappa_d C_d,jk (equal donor weights) and the leak."""
    m, C, kap = stats_["m"], stats_["C"], stats_["kappa"]; D, G = m.shape
    a, b = pairs if pairs is not None else np.triu_indices(G, 1)
    mc = m - m.mean(0)
    S = (mc[:, a] * mc[:, b]).sum(0) / (D - 1)
    leak = (kap[:, None] * C[:, a, b]).mean(0)
    return S - leak, leak


def lcp_test(stats_, pairs=None, B=20000, seed=0, refine_below=5e-4, B_refine=400_000, batch=200):
    """LCP p-values for all requested pairs (default: all j < k).
    Returns dict(p, T, est_bces, leak, pairs)."""
    m, psi, C, kap = stats_["m"], stats_["psi"], stats_["C"], stats_["kappa"]
    D, G = m.shape
    a, b = pairs if pairs is not None else np.triu_indices(G, 1)
    tau2, mu = paule_mandel(m, psi)
    v = tau2[None, :] + psi
    R = (m - mu[None, :]) / np.sqrt(v)
    Lt = kap[:, None] * C[:, a, b] / np.sqrt(v[:, a] * v[:, b])
    T = ((R[:, a] * R[:, b]).sum(0) - Lt.sum(0)) * np.sqrt(D / (D + (Lt ** 2).sum(0)))
    aT = np.abs(T) - 1e-12
    rng = np.random.default_rng(seed); cnt = np.zeros(len(a))
    for s0 in range(0, B, batch):
        bb = min(batch, B - s0)
        perms = np.argsort(rng.random((bb, D)), axis=1)
        Tp = np.matmul(R[perms].transpose(0, 2, 1), R)[:, a, b]
        cnt += (np.abs(Tp) >= aT).sum(0)
    p = (cnt + 1) / (B + 1)
    idx = np.where(p <= refine_below)[0]
    if len(idx) and B_refine:
        PERM = np.argsort(rng.random((B_refine, D)), axis=1).astype(np.int32)
        X = R[:, a[idx]]; Yv = R[:, b[idx]]; c2 = np.zeros(len(idx))
        ch = max(1, int(2e7 // max(1, D * len(idx))))
        for s0 in range(0, B_refine, ch):
            Tp = np.einsum("bdp,dp->bp", X[PERM[s0:s0 + ch]], Yv)
            c2 += (np.abs(Tp) >= aT[idx][None, :]).sum(0)
        p[idx] = (c2 + 1) / (B_refine + 1)
    est, leak = bces_estimate(stats_, (a, b))
    return dict(p=p, T=T, est_bces=est, leak=leak, pairs=(a, b))
