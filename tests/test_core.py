"""Fast checks: leak identity, Paule-Mandel root, null calibration of LCP on a small simulation."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, ".."))
from lcp import donor_statistics, lcp_test, paule_mandel  # noqa: E402


def _sim(D=30, n=60, G=6, rho_w=0.6, seed=0):
    """Genes 0-1 co-vary within donors only (between-null); genes 2-3 co-vary between donors only."""
    rng = np.random.default_rng(seed); Y = []; s = []; don = []
    for d in range(D):
        ub = rng.normal(); mu = np.full(G, 2e-3) * np.exp(0.3 * np.r_[rng.normal(), rng.normal(), ub, ub, rng.normal(), rng.normal()] - 0.045)
        z = rng.normal(size=n); e = rng.normal(size=(n, G))
        w = np.exp(0.4 * np.c_[z, rho_w * z + np.sqrt(1 - rho_w ** 2) * e[:, 1], e[:, 2:]] - 0.08)
        theta = mu * w; sd = np.exp(rng.normal(8.5, 0.5, n))
        Y.append(rng.poisson(sd[:, None] * theta)); s.append(sd); don += [f"d{d}"] * n
    return np.vstack(Y), np.concatenate(s), np.array(don)


def test_pm_root():
    Y, s, don = _sim()
    st = donor_statistics(Y, s, don); tau2, mu = paule_mandel(st["m"], st["psi"])
    r2 = ((st["m"] - mu) ** 2 / (tau2 + st["psi"])).sum(0)
    D = st["m"].shape[0]
    assert np.all((np.abs(r2 - (D - 1)) < 1e-3 * D) | (tau2 == 0))


def test_leak_positive_for_within_pair():
    Y, s, don = _sim()
    st = donor_statistics(Y, s, don)
    leak = (st["kappa"] * st["C"][:, 0, 1]).mean()
    assert leak > 0


def test_signal_detected_and_shapes():
    Y, s, don = _sim(D=40)
    st = donor_statistics(Y, s, don); o = lcp_test(st, B=2000, B_refine=0)
    a, b = o["pairs"]; k = np.where((a == 2) & (b == 3))[0][0]
    assert o["p"].shape == (15,) and o["p"][k] < 0.01


def test_null_pvalues_roughly_uniform():
    ps = []
    for sd in range(40):
        Y, s, don = _sim(D=30, seed=100 + sd)
        st = donor_statistics(Y, s, don); o = lcp_test(st, B=500, B_refine=0, seed=sd)
        a, b = o["pairs"]; ps.append(o["p"][np.where((a == 0) & (b == 1))[0][0]])   # within-only pair
    assert 0.0 <= np.mean(np.array(ps) < 0.1) <= 0.25
