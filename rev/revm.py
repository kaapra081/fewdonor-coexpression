"""REV_PREREG: XDS (normal-moment), XDS-M (multi-split, studentized) and KX (kappa extrapolation) from cell lists."""
import os, sys, numpy as np
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("dstd", "pop", "dev"): sys.path.insert(0, os.path.join(HERE, "..", p))
from dstd import _ms, design_R  # noqa: E402
from pop_core import normal_perm_p  # noqa: E402
from methods_v2 import pearson_fisher  # noqa: E402


def stats_from_cells(cells, seed, S=5):
    """Full-donor stats and S random half-splits (split 0 uses seed+7 to match XDS in the paper)."""
    full = {k: [] for k in ("m", "psi", "kappa", "S")}; halves = [{k: [] for k in ("mA", "psiA", "kappaA", "SA", "mB", "psiB", "kappaB", "SB")} for _ in range(S)]
    rngs = [np.random.default_rng(seed + 7)] + [np.random.default_rng(seed + 100 + s) for s in range(1, S)]
    for Y, s in cells:
        m, ps, kap, Ssum, _, n = _ms(Y, s)
        full["m"].append(m); full["psi"].append(ps); full["kappa"].append(kap); full["S"].append(Ssum)
        for h, rg in zip(halves, rngs):
            p = rg.permutation(n)
            for lab, rr in (("A", p[: n // 2]), ("B", p[n // 2:])):
                m2, ps2, k2, S2, _, _ = _ms(Y[rr], s[rr])
                h["m" + lab].append(m2); h["psi" + lab].append(ps2); h["kappa" + lab].append(k2); h["S" + lab].append(S2)
    full = {k: np.array(v) for k, v in full.items()}; halves = [{k: np.array(v) for k, v in h.items()} for h in halves]
    return full, halves


def methods(full, halves, iu, chunk=10000):
    a, b = iu; D = full["m"].shape[0]; P = len(a); out = {}
    out["PB"] = pearson_fisher(np.log1p(1e6 * full["m"]), iu)["p"]
    Rs = [(design_R(h["mA"], h["psiA"], h["kappaA"], h["SA"])[0], design_R(h["mB"], h["psiB"], h["kappaB"], h["SB"])[0]) for h in halves]
    RA, RB = Rs[0]; T = (RA[:, a] * RB[:, b]).sum(0); out["XDS"] = normal_perm_p(RA, RB, a, b, T); out["sign_XDS"] = np.sign(T)
    R, v = design_R(full["m"], full["psi"], full["kappa"], full["S"])
    pM = np.empty(P); pK = np.empty(P); sM = np.empty(P)
    for c0 in range(0, P, chunk):
        aa, bb = a[c0:c0 + chunk], b[c0:c0 + chunk]
        q = np.zeros((D, len(aa)))
        for RA_, RB_ in Rs: q += RA_[:, aa] * RB_[:, bb]
        q /= len(Rs); t = q.mean(0) / (q.std(0, ddof=1) / np.sqrt(D)); pM[c0:c0 + chunk] = 2 * norm.sf(np.abs(t)); sM[c0:c0 + chunk] = np.sign(t)
        y = R[:, aa] * R[:, bb]; x = full["kappa"][:, None] / np.sqrt(v[:, aa] * v[:, bb])
        xm = x.mean(0); ym = y.mean(0); sxx = ((x - xm) ** 2).sum(0); beta = ((x - xm) * (y - ym)).sum(0) / sxx; alpha = ym - beta * xm
        e = y - alpha - beta * x
        w = 1.0 / D - xm * (x - xm) / sxx            # intercept = sum_d w_d y_d
        se = np.sqrt((w ** 2 * e ** 2).sum(0)); pK[c0:c0 + chunk] = 2 * norm.sf(np.abs(alpha / se))
    out["XDSM"] = pM; out["KX"] = pK; out["sign_XDSM"] = sM
    return out
