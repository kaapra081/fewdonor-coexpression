"""recomb_v2 development simulator: few donors / few cells, exact truth on the declared scale.

Extends round3/benchmark/bench_sim.py (log-normal programs) with
  (i) gene-specific (idiosyncratic) donor effects  b_dj = exp(eta_j e_dj - eta_j^2/2), e ~ N(0,1),
      so that every gene has between-donor variance and rho_B < 1 for program pairs
      (real SEA-AD median between-donor CV^2 0.036, Lee 0.085: recomb_v2/real_between_variance_calibration.json);
  (ii) a 'strong within-only' program S (larger loading lam_S) that maximises leakage of within-donor
      covariance into donor means (between-null, within strongly non-null);
  (iii) a between-donor HETEROGENEITY NULL program Q: genes in half Q1 get donor effect (1 + g_Q u_d),
      genes in half Q2 get (1 + g_Q u_d s_d), u ~ U(-sqrt3, sqrt3), s_d = +-1 per donor. Cross pairs Q1-Q2 have
      Sigma_B,jk = g_Q^2 E[u^2 s] = 0 exactly but DEPENDENT donor effects (weak null); within-half pairs have
      Sigma_B,jk = mu_j mu_k g_Q^2 (non-null);
  (iv) synthetic few-donor designs: D donors, n_d ~ round(LogNormal(log nd_median, nd_logsd)) clipped at nd_min.

Latent model (cell i, donor d, gene j)
  theta_dij = mu_j * beta_dj * b_dj * omega_dij * xi_dij
  beta_dj  = prod_h exp(g_jh u_dh - g_jh^2/2)         u ~ N(0,1)   programs B, M, O (shared donor factors)
  omega_dij = prod_f exp(l_djf z_dif - l_djf^2/2)      z ~ N(0,1)   programs W, M, O, S (cell factors)
  H program: omega *= (1 + l_dj z_di), z ~ U(-sqrt3,sqrt3), l_dj = lam_H * (+-1 iid per donor, gene)
  xi ~ Gamma(phi_j, 1/phi_j)
  Y ~ Poisson(s_di theta_dij)
Exact truth (j != k), equal donor weights:
  Sigma_B,jk = mu_j mu_k (exp(sum_h g_jh g_kh) - 1)                   (idiosyncratic b independent across genes)
  Sigma_B,jj = mu_j^2 (exp(sum_h g_jh^2 + eta_j^2) - 1)
  Sigma_W,jk = mu_j mu_k exp(sum_h g_jh g_kh) (E_l exp(sum_f l_jf l_kf) - 1)   (0 for H pairs on average)
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                                            # bundled copy (remote runs)
sys.path.append(os.path.join(HERE, "..", "..", "round3", "diagnostics"))
import diag_sim as ds  # noqa: E402  (only FITTED constants are used)

SQ3 = np.sqrt(3.0)
V2_COUNTS = dict(W=10, S=8, B=10, M=8, O=8, H=10, Q=10, N=16)   # G = 80
LOGN_PROGS = ("W", "S", "M", "O")
BETWEEN_PROGS = ("B", "M", "O")


def make_panel(seed=26102026, counts=None, eta_mu=np.log(0.19), eta_sd=0.35):
    """Gene panel. mu_j, phi_j from the SEA-AD fit (diag_sim.FITTED); eta_j (idiosyncratic donor log-SD)
    ~ LogNormal(log 0.19, 0.35): median between-donor log-SD of SEA-AD microglia genes is 0.19."""
    counts = dict(V2_COUNTS if counts is None else counts)
    rng = np.random.default_rng(seed)
    labels = np.array(sum([[k] * n for k, n in counts.items()], []))
    G = len(labels)
    F = ds.FITTED
    mu = np.clip(np.exp(rng.normal(F["mu_p"], F["sigma_p"], G)), 3e-5, 2e-2)
    phi = np.clip(np.exp(rng.normal(F["mu_phi"], F["sigma_phi"], G)), 0.2, 50.0)
    eta = np.exp(rng.normal(eta_mu, eta_sd, G))
    osign = np.zeros(G)
    oi = np.where(labels == "O")[0]
    osign[oi[: len(oi) // 2]] = 1.0
    osign[oi[len(oi) // 2:]] = -1.0
    qhalf = np.zeros(G)
    qi = np.where(labels == "Q")[0]
    qhalf[qi[: len(qi) // 2]] = 1.0
    qhalf[qi[len(qi) // 2:]] = 2.0
    return dict(labels=labels, mu=mu, phi=phi, eta=eta, osign=osign, qhalf=qhalf, G=G, seed=seed, counts=counts)


def setting(**kw):
    st = dict(D=16, nd_median=40, nd_logsd=ds.FITTED["cpd_logsd"], nd_min=10, nd_fixed=None,
              lam=0.30, lam_S=0.60, g=0.30, lam_H=0.30, g_Q=0.30, eta_mult=1.0, het_h=0.0,
              sa_mult=1.0, se_mult=1.0, kappa=0.0)
    st.update(kw)
    return st


def _gload(panel, g):
    lab, os_ = panel["labels"], panel["osign"]
    Gm = {h: np.zeros(panel["G"]) for h in BETWEEN_PROGS}
    Gm["B"][lab == "B"] = g
    Gm["M"][lab == "M"] = g
    Gm["O"][lab == "O"] = g * os_[lab == "O"]
    return Gm


def _lam(st, f):
    return st["lam_S"] if f == "S" else st["lam"]


def truth(panel, st, n_mc=200000, seed=11):
    G, lab, mu = panel["G"], panel["labels"], panel["mu"]
    eta = panel["eta"] * st["eta_mult"]
    Gm = _gload(panel, st["g"])
    GG = sum(np.outer(Gm[h], Gm[h]) for h in BETWEEN_PROGS)
    rng = np.random.default_rng(seed)
    m = rng.uniform(1 - st["het_h"], 1 + st["het_h"], n_mc) if st["het_h"] > 0 else np.ones(1)
    EW = np.zeros((G, G))
    for f in LOGN_PROGS:
        idx = lab == f
        EW[np.ix_(idx, idx)] = float(np.mean(np.exp((_lam(st, f) * m) ** 2))) - 1.0
    SB = np.outer(mu, mu) * (np.exp(GG) - 1.0)
    qh = panel["qhalf"]; gq = st["g_Q"]
    QQ = gq ** 2 * ((qh[:, None] == qh[None, :]) & (qh[:, None] > 0))          # E[(beta_Q,j - 1)(beta_Q,k - 1)]
    # Q genes carry no B/M/O loading, so E[beta_j beta_k] - 1 = QQ for Q pairs; product with idiosyncratic b
    SB = SB + np.outer(mu, mu) * QQ
    dQ = np.where(qh > 0, gq ** 2, 0.0)
    np.fill_diagonal(SB, mu ** 2 * ((1.0 + dQ) * np.exp(np.diag(GG) + eta ** 2) - 1.0))
    SW = np.outer(mu, mu) * np.exp(GG) * EW
    np.fill_diagonal(SW, np.nan)
    iu = np.triu_indices(G, 1)
    li, lk = lab[iu[0]], lab[iu[1]]
    same = li == lk
    cls = np.full(len(li), "NN", dtype=object)
    for P, name in (("W", "WW"), ("S", "SS"), ("B", "BB"), ("M", "MM"), ("H", "HH")):
        cls[same & (li == P)] = name
    o = same & (li == "O")
    osg = panel["osign"]
    cls[o & (osg[iu[0]] == osg[iu[1]])] = "OO_same"
    cls[o & (osg[iu[0]] != osg[iu[1]])] = "OO_opp"
    q = same & (li == "Q")
    cls[q & (qh[iu[0]] == qh[iu[1]])] = "QQ_same"
    cls[q & (qh[iu[0]] != qh[iu[1]])] = "QQ_cross"
    dB = np.diag(SB)
    rB = SB[iu] / np.sqrt(dB[iu[0]] * dB[iu[1]])
    return dict(SB=SB, SW=SW, iu=iu, cls=cls, sb=SB[iu], rB=rB,
                altB=np.abs(SB[iu]) > 0, altW=np.abs(SW[iu]) > 0)


def make_nd(st, rng):
    if st.get("nd_fixed"):
        return np.full(st["D"], int(st["nd_fixed"]))
    x = rng.lognormal(np.log(st["nd_median"]), st["nd_logsd"], st["D"])
    return np.clip(np.round(x), st["nd_min"], None).astype(int)


def simulate(panel, st, seed):
    rng = np.random.default_rng(seed)
    G, lab, mu, phi = panel["G"], panel["labels"], panel["mu"], panel["phi"]
    eta = panel["eta"] * st["eta_mult"]
    D = st["D"]; nd = make_nd(st, rng)
    a_d = rng.normal(0.0, ds.FITTED["sigma_a"] * st["sa_mult"], D)
    Gm = _gload(panel, st["g"])
    u = {h: rng.normal(0.0, 1.0, D) for h in BETWEEN_PROGS}
    qh = panel["qhalf"]; gq = st["g_Q"]
    assert gq * SQ3 < 1
    uq = rng.uniform(-SQ3, SQ3, D); sq = rng.choice([-1.0, 1.0], size=D)
    mag = {f: (rng.uniform(1 - st["het_h"], 1 + st["het_h"], D) if st["het_h"] > 0 else np.ones(D))
           for f in LOGN_PROGS}
    idx = {f: np.where(lab == f)[0] for f in LOGN_PROGS + ("H",)}
    bg = max(1.0 - mu.sum(), 0.5)
    Ys, S, Sh, rows, MU, PSI, r0 = [], [], [], [], [], [], 0
    for d in range(D):
        n = int(nd[d])
        logb = sum(Gm[h] * u[h][d] - 0.5 * Gm[h] ** 2 for h in BETWEEN_PROGS)
        logb = logb + eta * rng.normal(0.0, 1.0, G) - 0.5 * eta ** 2
        beta = np.exp(logb)
        beta = beta * np.where(qh == 1, 1.0 + gq * uq[d], np.where(qh == 2, 1.0 + gq * uq[d] * sq[d], 1.0))
        logw = np.zeros((n, G)); zW = None
        for f in LOGN_PROGS:
            z = rng.normal(0.0, 1.0, n)
            if f == "W":
                zW = z
            l = _lam(st, f) * mag[f][d]
            logw[:, idx[f]] += (l * z - 0.5 * l * l)[:, None]
        omega = np.exp(logw)
        Eo2 = np.ones(G)
        for f in LOGN_PROGS:
            Eo2[idx[f]] = np.exp((_lam(st, f) * mag[f][d]) ** 2)
        Eo2[idx["H"]] = 1.0 + st["lam_H"] ** 2
        zH = rng.uniform(-SQ3, SQ3, n)
        omega[:, idx["H"]] *= 1.0 + np.outer(zH, st["lam_H"] * rng.choice([-1.0, 1.0], size=len(idx["H"])))
        xi = rng.gamma(phi, 1.0 / phi, size=(n, G))
        theta = mu[None, :] * beta[None, :] * omega * xi
        s = np.exp(ds.FITTED["m_L"] + a_d[d] + rng.normal(0.0, ds.FITTED["sigma_e"] * st["se_mult"], n))
        if st["kappa"] != 0.0:
            s = s * np.exp(st["kappa"] * zW)
        Y = rng.poisson(s[:, None] * theta)
        T = rng.poisson(s * bg)
        Ys.append(Y.astype(np.int32)); S.append(s); Sh.append(Y.sum(1) + T)
        rows.append(np.arange(r0, r0 + n)); r0 += n
        MU.append(mu * beta)
        mud = mu * beta; S_ = s.sum(); kap = (s ** 2).sum() / S_ ** 2
        PSI.append(kap * mud ** 2 * (Eo2 * (1.0 + 1.0 / phi) - 1.0) + mud / S_)   # true Var(m_dj | donor d)
    return dict(Y=np.vstack(Ys), s_true=np.concatenate(S), s_hat=np.concatenate(Sh).astype(np.float64),
                donor_idx=rows, nd=nd, mu_true=np.array(MU), psi_true=np.array(PSI))


def numeric_truth_check(panel, st, D=4000, n=20, seed=5):
    st2 = dict(st); st2.update(D=D, nd_fixed=n)
    sim = simulate(panel, st2, seed)
    return np.cov(sim["mu_true"], rowvar=False), truth(panel, st)["SB"]
