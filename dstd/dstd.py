"""DSTD_PREREG.md: LCPd, LCPXd (design-only standardization) and PBlin."""
import os, sys, numpy as np
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "pop")); sys.path.insert(0, os.path.join(HERE, "..", "dev"))
from pop_core import perm_counts, normal_perm_p  # noqa: E402
from methods_v2 import paule_mandel, pearson_fisher  # noqa: E402


def _ms(Y, s):
    S = s.sum(); x = Y / s[:, None]; m = Y.sum(0) / S; kap = (s ** 2).sum() / S ** 2; n = len(s); xc = x - x.mean(0)
    sig2 = (xc ** 2).sum(0) / (n - 1) - (x / s[:, None]).mean(0)
    return m, np.maximum(kap * np.maximum(sig2, 0) + m / S, 1e-30), kap, S, xc, n


def donor_level_cells(cells, split_seed, need_C=True):
    """cells: iterable of (Y_d dense [n x G], s_d [n]). Returns full and half-split donor statistics."""
    rng = np.random.default_rng(split_seed); o = {k: [] for k in ("m", "psi", "kappa", "S", "C", "n")}
    for h in ("A", "B"):
        for k in ("m", "psi", "kappa", "S"): o[k + h] = []
    for Y, s in cells:
        m, ps, kap, S, xc, n = _ms(Y, s)
        o["m"].append(m); o["psi"].append(ps); o["kappa"].append(kap); o["S"].append(S); o["n"].append(n)
        if need_C: o["C"].append((xc.T @ xc / (n - 1)).astype(np.float32))
        p = rng.permutation(n)
        for h, rr in (("A", p[: n // 2]), ("B", p[n // 2:])):
            m2, ps2, k2, S2, _, _ = _ms(Y[rr], s[rr])
            o["m" + h].append(m2); o["psi" + h].append(ps2); o["kappa" + h].append(k2); o["S" + h].append(S2)
    d = {k: np.array(v) for k, v in o.items() if len(v)}; d["donors"] = np.arange(len(d["m"]))
    return d


def design_R(m, psi, kap, S):
    mu = (m * S[:, None]).sum(0) / S.sum()
    sig2 = np.maximum(((psi - m / S[:, None]) / kap[:, None]).mean(0), 0)
    ps = kap[:, None] * sig2[None, :] + mu[None, :] / S[:, None]
    tau2, mu2, _ = paule_mandel(m, ps); v = tau2[None, :] + ps
    return (m - mu2[None, :]) / np.sqrt(v), v


def methods_d(d, iu, B=100_000, seed=0, chunk=20000, perm_max_D=100):
    a, b = iu; D = len(d["m"]); out = {}
    out["PBlin"] = pearson_fisher(d["m"], iu)["p"]
    R, v = design_R(d["m"], d["psi"], d["kappa"], d["S"]); obs = np.empty(len(a))
    for c0 in range(0, len(a), chunk):
        aa, bb = a[c0:c0 + chunk], b[c0:c0 + chunk]
        Lt = d["kappa"][:, None] * d["C"][:, aa, bb].astype(np.float64) / np.sqrt(v[:, aa] * v[:, bb])
        obs[c0:c0 + chunk] = ((R[:, aa] * R[:, bb]).sum(0) - Lt.sum(0)) * np.sqrt(D / (D + (Lt ** 2).sum(0)))
    RA, _ = design_R(d["mA"], d["psiA"], d["kappaA"], d["SA"]); RB, _ = design_R(d["mB"], d["psiB"], d["kappaB"], d["SB"])
    T = (RA[:, a] * RB[:, b]).sum(0)
    if D <= perm_max_D:
        out["LCPd"] = perm_counts(R, R, a, b, obs, B=B, seed=seed); out["LCPXd"] = perm_counts(RA, RB, a, b, T, B=B, seed=seed + 1)
    else:
        out["LCPd"] = normal_perm_p(R, R, a, b, obs); out["LCPXd"] = normal_perm_p(RA, RB, a, b, T)
    out["sign_LCPd"] = np.sign(obs)
    return out
