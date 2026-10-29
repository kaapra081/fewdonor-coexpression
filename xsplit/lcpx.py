"""LCPX and LCPX-S (XS_PREREG.md): cross-split donor-permutation tests with no leak estimate."""
import numpy as np
import torch


def half_stats(Y, s):
    S = s.sum(); x = Y / s[:, None]; m = Y.sum(0) / S; kap = (s ** 2).sum() / S ** 2; n = len(s)
    xc = x - x.mean(0); sig2 = (xc ** 2).sum(0) / (n - 1) - (x / s[:, None]).mean(0)
    return m, np.maximum(kap * np.maximum(sig2, 0) + m / S, 1e-30)


def split_R(cells, seed, paule_mandel):
    rng = np.random.default_rng(seed); out = {"A": ([], []), "B": ([], [])}
    for Y, s in cells:
        n = len(s); idx = rng.permutation(n); A, B = idx[: n // 2], idx[n // 2:]
        for h, r in (("A", A), ("B", B)):
            m, p = half_stats(Y[r], s[r]); out[h][0].append(m); out[h][1].append(p)
    R = {}
    for h in ("A", "B"):
        m = np.array(out[h][0]); psi = np.array(out[h][1]); tau2, mu, _ = paule_mandel(m, psi)
        R[h] = (m - mu[None, :]) / np.sqrt(tau2[None, :] + psi)
    return R["A"], R["B"]


def lcpx(cells, iu, B=100_000, seed=0, paule_mandel=None, batch=200, device="cpu"):
    RA, RB = split_R(cells, seed + 7, paule_mandel); D, G = RA.shape; a, b = iu
    q = RA[:, a] * RB[:, b]
    T = q.sum(0); sd = q.std(0, ddof=1); t = q.mean(0) / np.maximum(sd, 1e-300) * np.sqrt(D)
    A_ = torch.tensor(RA, dtype=torch.float32); B_ = torch.tensor(RB, dtype=torch.float32)
    A2 = A_ ** 2; B2 = B_ ** 2
    ai = torch.tensor(a); bi = torch.tensor(b)
    oT = torch.tensor(np.abs(T) - 1e-5, dtype=torch.float32); ot = torch.tensor(np.abs(t) - 1e-5, dtype=torch.float32)
    cT = torch.zeros(len(a), dtype=torch.int64); ct = torch.zeros(len(a), dtype=torch.int64)
    g = torch.Generator(); g.manual_seed(int(seed) % (2 ** 31)); done = 0
    At = A_.T.unsqueeze(0); A2t = A2.T.unsqueeze(0)
    while done < B:
        bb = min(batch, B - done)
        perms = torch.argsort(torch.rand(bb, D, generator=g), dim=1)
        Bp = B_[perms]; S1 = torch.matmul(At, Bp)[:, ai, bi]                    # sum_d RA_dj RB_pi(d)k
        S2 = torch.matmul(A2t, B2[perms])[:, ai, bi]                              # sum_d (RA_dj RB_pi(d)k)^2
        var = (S2 - S1 ** 2 / D) / (D - 1)
        tp = (S1 / D) / torch.sqrt(torch.clamp(var, min=1e-30)) * np.sqrt(D)
        cT += (S1.abs() >= oT).sum(0); ct += (tp.abs() >= ot).sum(0); done += bb
    return dict(p_LCPX=(cT.numpy() + 1) / (B + 1), p_LCPXS=(ct.numpy() + 1) / (B + 1), T=T, t=t)
