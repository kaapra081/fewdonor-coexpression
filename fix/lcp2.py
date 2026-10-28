"""LCP and LCP2 (FIX_PREREG.md Item 2) with shared donor permutations.
LCP2 adds max(sum_d delta_d, 0) to the scaling denominator, where delta_d = jackknife Var(z_d) - jackknife Var(r_dj r_dk)
over the nuclei of donor d (v, mu-hat held fixed), i.e. Var(leak estimate) - 2 Cov(r r, leak estimate)."""
import numpy as np
import torch


def jack_delta(cells, mu, v, a, b, chunk=20000):
    """cells: list of (Y_d [n x G] float array, s_d [n]); mu [G]; v [D x G]. Returns delta [D x P] (one row per donor)."""
    D = len(cells); P = len(a); out = np.zeros((D, P))
    for d, (Y, s) in enumerate(cells):
        n = len(s); S = s.sum(); Q = (s ** 2).sum(); m = Y.sum(0) / S
        x = Y / s[:, None]; xb = x.mean(0)
        mi = (S * m[None, :] - Y) / (S - s)[:, None]
        dev = (mi - mu[None, :]) / np.sqrt(v[d])[None, :]
        ki = (Q - s ** 2) / (S - s) ** 2
        xbi = (n * xb[None, :] - x) / (n - 1)
        for c0 in range(0, P, chunk):
            aa, bb = a[c0:c0 + chunk], b[c0:c0 + chunk]
            rr = dev[:, aa] * dev[:, bb]
            Sxx = (x[:, aa] * x[:, bb]).sum(0)
            Ci = (Sxx[None, :] - x[:, aa] * x[:, bb] - (n - 1) * xbi[:, aa] * xbi[:, bb]) / (n - 2)
            li = ki[:, None] * Ci / np.sqrt(v[d, aa] * v[d, bb])[None, :]
            z = rr - li
            Vz = (n - 1) / n * ((z - z.mean(0)) ** 2).sum(0)
            Vr = (n - 1) / n * ((rr - rr.mean(0)) ** 2).sum(0)
            out[d, c0:c0 + chunk] = Vz - Vr
    return out


def lcp_both(m, psi, C, kap, cells, iu, B=100_000, seed=0, paule_mandel=None, batch=400, device="cpu"):
    D, G = m.shape; a, b = iu
    tau2, mu, _ = paule_mandel(m, psi); v = tau2[None, :] + psi
    R = (m - mu[None, :]) / np.sqrt(v)
    Lt = kap[:, None] * C[:, a, b] / np.sqrt(v[:, a] * v[:, b])
    num = (R[:, a] * R[:, b]).sum(0) - Lt.sum(0)
    obs1 = num * np.sqrt(D / (D + (Lt ** 2).sum(0)))
    dl = jack_delta(cells, mu, v, a, b).sum(0)
    obs2 = num * np.sqrt(D / (D + (Lt ** 2).sum(0) + np.maximum(dl, 0)))
    Rt = torch.tensor(R, dtype=torch.float32, device=device)
    ai = torch.tensor(a, device=device); bi = torch.tensor(b, device=device)
    o1 = torch.tensor(np.abs(obs1) - 1e-6, dtype=torch.float32); o2 = torch.tensor(np.abs(obs2) - 1e-6, dtype=torch.float32)
    c1 = torch.zeros(len(a), dtype=torch.int64); c2 = torch.zeros(len(a), dtype=torch.int64)
    g = torch.Generator(device=device); g.manual_seed(int(seed) % (2 ** 31)); done = 0
    while done < B:
        bb = min(batch, B - done)
        perms = torch.argsort(torch.rand(bb, D, generator=g, device=device), dim=1)
        T = torch.bmm(Rt[perms].transpose(1, 2), Rt.expand(bb, D, G))[:, ai, bi].abs()
        c1 += (T >= o1).sum(0); c2 += (T >= o2).sum(0); done += bb
    return dict(p_LCP=(c1.numpy() + 1) / (B + 1), p_LCP2=(c2.numpy() + 1) / (B + 1), T1=obs1, T2=obs2, delta_sum=dl)
