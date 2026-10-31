"""Smoke tests of the paper's pipeline on tiny synthetic data (all of tests/ runs in well under a minute; no data files needed).

The simulation is a miniature of the paper's real-cell null: 80 donors with 5-15 cells each, 30 genes (435 pairs), no between-donor
covariance at all, but 12 genes share a within-donor cell-state factor. Every call is therefore false. Pseudobulk (PB) is expected to
call dozens of the 66 within-program pairs, while XDS, XDS-M and KX are expected to call none. Thresholds sit far from the values
observed over 10 seeds (PB 66-80 BH calls; XDS-M 0 calls, pooled per-pair type-I 0.074), so they test the pipeline, not the exact numbers.
"""
import os, sys
import numpy as np
import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)         # for `from lcp import bh`; the test file must not rely on another test having set this
for p in ("rev", "dev", "dstd", "pop", "xsplit"): sys.path.insert(0, os.path.join(ROOT, p))
pytest.importorskip("torch")            # rev/revm.py imports pop/pop_core.py, which imports torch
from methods_v2 import bh               # noqa: E402
from revm import stats_from_cells, methods   # noqa: E402

G, N_PROG = 30, 12
IU = np.triu_indices(G, 1)
PROG = (IU[0] < N_PROG) & (IU[1] < N_PROG)


def sim_cells(D=80, nd=(5, 15), seed=0, strength=0.8, between=None):
    """List of (counts [n x G], library size [n]) per donor. `between` = (j, k, sd): genes j, k share a donor-level log-rate factor."""
    rng = np.random.default_rng(seed)
    mu = np.exp(rng.normal(np.log(2e-3), 0.5, G)); load = np.zeros(G); load[:N_PROG] = 1.0; cells = []
    for _ in range(D):
        n = int(rng.integers(nd[0], nd[1] + 1)); s = np.exp(rng.normal(np.log(4000), 0.4, n)); z = rng.normal(size=n)
        w = np.exp(strength * z[:, None] * load[None, :] - 0.5 * strength ** 2 * load[None, :])
        if between is not None:
            j, k, sd = between; u = rng.normal() * sd; w[:, j] *= np.exp(u); w[:, k] *= np.exp(u)
        cells.append((rng.poisson(s[:, None] * mu[None, :] * w).astype(np.float64), s))
    return cells


@pytest.fixture(scope="module")
def null_runs():
    out = []
    for seed in range(4):
        f, h = stats_from_cells(sim_cells(seed=seed), seed=seed); out.append(methods(f, h, IU))
    return out


def test_outputs_are_valid_p_values(null_runs):
    for o in null_runs:
        for m in ("PB", "XDS", "XDSM", "KX"):
            assert o[m].shape == (len(IU[0]),) and np.all(np.isfinite(o[m])) and o[m].min() > 0 and o[m].max() <= 1, m
        assert set(np.unique(o["sign_XDS"])) <= {-1.0, 0.0, 1.0}


def test_pipeline_is_deterministic():
    runs = [methods(*stats_from_cells(sim_cells(seed=7), seed=7), IU) for _ in range(2)]
    for m in ("PB", "XDS", "XDSM", "KX"): assert np.array_equal(runs[0][m], runs[1][m])


def test_pseudobulk_fails_and_cross_split_tests_hold_on_a_complete_null(null_runs):
    for o in null_runs:
        assert bh(o["PB"]).sum() >= 20 and (o["PB"][PROG] < 0.05).mean() > 0.9     # leak: within-donor program read as between-donor
        for m in ("XDS", "XDSM", "KX"): assert bh(o[m]).sum() <= 2, m               # complete null: any call is a false call
    for m in ("XDS", "XDSM", "KX"):                                                  # per-pair type-I roughly nominal
        assert 0.02 < np.mean(np.concatenate([o[m] for o in null_runs]) < 0.05) < 0.11, m


def test_planted_between_donor_pair_is_detected():
    k = int(np.where((IU[0] == 20) & (IU[1] == 21))[0][0])
    f, h = stats_from_cells(sim_cells(D=120, nd=(10, 30), seed=3, between=(20, 21, 0.7)), seed=3); o = methods(f, h, IU)
    assert o["XDS"][k] < 1e-6 and o["XDSM"][k] < 0.05     # XDS-M is a t-test on skewed cross-products; much less powerful than XDS here (about 2e-3)


def _bh_reference(p, q):
    """Step-up rule written independently of the repository code: k = max{i : p_(i) <= q i / m}; reject every p <= p_(k) (ties go together)."""
    p = np.asarray(p, float); fin = np.where(np.isfinite(p))[0]; ps = sorted(p[fin]); m = len(ps); rej = np.zeros(len(p), bool); k = 0
    for i, v in enumerate(ps, 1):
        if v <= q * i / m: k = i
    if k: rej[fin] = p[fin] <= ps[k - 1]
    return rej


def test_bh_known_example_and_agreement_with_reference():
    p = [0.001, 0.008, 0.039, 0.041, 0.042, 0.060, 0.074, 0.205, 0.212, 0.216]   # step-up thresholds 0.005, 0.010, 0.015, ... at q = 0.05
    assert bh(p).tolist() == [True, True] + [False] * 8
    assert bh(p, q=0.2).sum() == 7                                               # thresholds 0.02 i: p_(7) = 0.074 <= 0.14, p_(8) = 0.205 > 0.16
    rng = np.random.default_rng(0)
    from lcp import bh as bh_lcp
    for trial in range(200):
        m = int(rng.integers(1, 60)); q = float(rng.choice([0.01, 0.05, 0.2]))
        pv = rng.uniform(size=m) ** rng.uniform(1, 6); pv[rng.random(m) < 0.1] = np.nan; pv = np.round(pv, 3) if trial % 3 == 0 else pv   # NaN and ties
        assert np.array_equal(bh(pv, q), _bh_reference(pv, q)) and np.array_equal(bh_lcp(pv, q), _bh_reference(pv, q))
    assert not bh([np.nan, np.nan]).any() and bh([]).shape == (0,)
