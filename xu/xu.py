"""XU_PREREG.md: exact within-donor U-statistic, design-only weights, studentized across donors."""
import os, sys, numpy as np
from scipy.stats import norm, t as tdist
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "dev"))
from methods_v2 import paule_mandel  # noqa: E402


def design_v(m, psi, kap, S):
    mu = (m * S[:, None]).sum(0) / S.sum()
    sig2 = np.maximum(((psi - m / S[:, None]) / kap[:, None]).mean(0), 0)
    ps = kap[:, None] * sig2[None, :] + mu[None, :] / S[:, None]
    tau2, _, _ = paule_mandel(m, ps)
    return tau2[None, :] + ps


def xu(cells, m, psi, kap, S, iu, small_t=True):
    """cells: list/iterator of (Y_d [n x G], s_d [n]) in the same donor order as m (D x G). Returns p, t, est."""
    D, G = m.shape; a, b = iu
    v = design_v(m, psi, kap, S); isd = 1.0 / np.sqrt(v)
    c = (m.sum(0)[None, :] - m) / (D - 1)                      # leave-one-donor-out equal-weight centre
    s1 = np.zeros(len(a)); s2 = np.zeros(len(a))
    for d, (Y, s) in enumerate(cells):
        Yp = Y.sum(0); den = S[d] ** 2 - (s ** 2).sum()
        U = (np.outer(Yp, Yp) - Y.T @ Y) / den                    # theta_dj theta_dk unbiased (off-diagonal)
        Uab = U[a, b]
        q = isd[d, a] * isd[d, b] * (Uab - c[d, b] * m[d, a] - c[d, a] * m[d, b] + c[d, a] * c[d, b])
        s1 += q; s2 += q * q
    mean = s1 / D; var = (s2 - D * mean ** 2) / (D - 1)
    tt = mean / np.sqrt(np.maximum(var, 1e-300) / D)
    p = 2 * (tdist.sf(np.abs(tt), D - 1) if (small_t and D < 100) else norm.sf(np.abs(tt)))
    return dict(p=p, t=tt, est=mean)
