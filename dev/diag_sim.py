"""Round-3 diagnostic simulator with EXACT two-level truth on a declared scale.

Scale and estimands (declared before any method is run)
-------------------------------------------------------
X = latent relative expression theta (expected UMI per unit library size) of a randomly drawn cell.
Donor population = the generative distribution of donors below; donors are weighted EQUALLY.
    Sigma_W = E_D[ Cov(theta | D) ]        (population-average within-donor covariance)
    Sigma_B = Cov_D( E[theta | D] )        (between-donor covariance of donor mean expression)
Correlations rho_W, rho_B are the corresponding normalized matrices. Pooled-cell (marginal)
association is a cell-weighted mixture of both and is NOT a target here.

Latent model (cell i of donor d, gene j)
    theta_dij = mu_j * beta_dj * (1 + w_dij) * xi_dij
    beta_dj   = 1 + sum_h g_jh u_dh,        u_dh  iid U(-sqrt3, sqrt3)  (donor factors, var 1)
    w_dij     = sum_f l_djf z_dif,          z_dif iid U(-sqrt3, sqrt3)  (cell factors,  var 1)
    xi_dij    iid Gamma(phi_j, 1/phi_j)     (mean 1; gene-specific extra-Poisson noise)
Bounded factors keep theta > 0 exactly, so every moment below is exact (no truncation).
Programs (one factor each): W within-only; B between-only; M both (same sign); O both with
OPPOSITE between-donor sign across its two halves; H sign-heterogeneous within-donor null:
l_dj = lam_H * s_dj, s_dj iid +-1 per donor AND gene, so every donor has non-zero within-donor
covariance but the population average is exactly 0; N independent genes.
Loading heterogeneity (het='mag'): l_djf = lam * m_df with m_df ~ U(1-h, 1+h) shared by the
program's genes within donor d (same sign, varying magnitude).

Measurement
    Y_dij ~ Poisson(s_di * theta_dij),  s_di = exp(m_L + a_d + e_di) [* exp(kappa * z_di,W)]
    background T_di ~ Poisson(s_di * (1 - sum_j mu_j)); observed total UMI = sum_j Y_dij + T_di.
kappa != 0 couples depth to the W-program cell state (misspecification F: depth not ignorable).

Closed forms (independence of u, z, xi, loading signs/magnitudes)
    E[theta_j | d] = mu_j beta_dj
    Cov(theta_j, theta_k | d) = mu_j mu_k beta_dj beta_dk sum_f l_djf l_dkf            (j != k)
    Var(theta_j | d) = mu_j^2 beta_dj^2 [ (1 + sum_f l_djf^2)(1 + 1/phi_j) - 1 ]
    Sigma_W,jk = mu_j mu_k (1 + sum_h g_jh g_kh) E[sum_f l_djf l_dkf]
    Sigma_W,jj = mu_j^2 (1 + sum_h g_jh^2) [ (1 + E sum_f l_djf^2)(1 + 1/phi_j) - 1 ]
    Sigma_B,jk = mu_j mu_k sum_h g_jh g_kh                                           (all j, k)
Per-donor latent correlation (does not depend on beta):
    rho_d,jk = sum_f l_djf l_dkf / sqrt(A_dj A_dk),  A_dj = (1 + sum_f l_djf^2)(1 + 1/phi_j) - 1
"""
import numpy as np

SQ3 = np.sqrt(3.0)
FITTED = dict(m_L=8.563972204174574, sigma_a=0.23693207451419476, sigma_e=0.5784360094439615,
              mu_p=-8.367791515929262, sigma_p=0.6502578949974087,
              mu_phi=0.659813174941819, sigma_phi=0.7060692696515399,
              cpd_logsd=0.6623538497231674)   # runs/round2_benchmark/fitted_real_data_params.json (SEA-AD v3)
WITHIN_PROGS = ("W", "M", "O", "H")
BETWEEN_PROGS = ("B", "M", "O")
DEFAULT_COUNTS = dict(W=10, B=10, M=10, O=8, H=10, N=16)


def make_panel(seed=20261025, counts=None, mu_shift=0.0, min_mu=3e-5, gamma_noise=True):
    """Fixed gene panel. mu_j ~ LogNormal(mu_p + mu_shift, sigma_p) clipped to [min_mu, 2e-2];
    phi_j ~ LogNormal(mu_phi, sigma_phi) clipped to [0.2, 50] (inf -> Poisson-only latent)."""
    counts = dict(DEFAULT_COUNTS if counts is None else counts)
    rng = np.random.default_rng(seed)
    labels = np.array(sum([[k] * n for k, n in counts.items()], []))
    G = len(labels)
    mu = np.clip(np.exp(rng.normal(FITTED["mu_p"] + mu_shift, FITTED["sigma_p"], G)), min_mu, 2e-2)
    phi = (np.clip(np.exp(rng.normal(FITTED["mu_phi"], FITTED["sigma_phi"], G)), 0.2, 50.0)
           if gamma_noise else np.full(G, np.inf))
    osign = np.zeros(G)
    oi = np.where(labels == "O")[0]
    osign[oi[: len(oi) // 2]] = 1.0
    osign[oi[len(oi) // 2:]] = -1.0
    return dict(labels=labels, mu=mu, phi=phi, osign=osign, G=G, seed=seed, counts=counts)


def base_setting(**kw):
    st = dict(D=40, nd_mode="fixed", nd=200, nd_median=200, nd_logsd=FITTED["cpd_logsd"], nd_min=10,
              lam=0.3, g=0.3, lam_H=0.3, het="none", het_h=0.0, sa_mult=1.0, se_mult=1.0,
              kappa=0.0, latent_only=False, size_factor="estimated")
    st.update(kw)
    lim = max(st["lam"] * (1 + st["het_h"]), st["lam_H"]) * SQ3
    assert lim < 1 and st["g"] * SQ3 < 1, "loadings would allow theta <= 0"
    return st


def _between_loading(panel, g):
    lab, os_ = panel["labels"], panel["osign"]
    Gm = {h: np.zeros(panel["G"]) for h in BETWEEN_PROGS}
    Gm["B"][lab == "B"] = g
    Gm["M"][lab == "M"] = g
    Gm["O"][lab == "O"] = g * os_[lab == "O"]
    return Gm


def truth(panel, st, n_mc=200000, seed=7):
    """Exact Sigma_W, Sigma_B (G x G), rho_W, rho_B, pair classes, and the population mean of
    per-donor Fisher-z latent correlations (zbar_W; Monte Carlo over the loading distribution only,
    exact for het='none')."""
    G, lab, mu, phi = panel["G"], panel["labels"], panel["mu"], panel["phi"]
    lam, g, lamH = st["lam"], st["g"], st["lam_H"]
    inv_phi = np.where(np.isfinite(phi), 1.0 / phi, 0.0)
    Gm = _between_loading(panel, g)
    GG = sum(np.outer(Gm[h], Gm[h]) for h in BETWEEN_PROGS)          # sum_h g_jh g_kh
    Em2 = 1.0 + st["het_h"] ** 2 / 3.0 if st["het"] == "mag" else 1.0  # E[m^2], m ~ U(1-h,1+h)
    LL = np.zeros((G, G)); L2 = np.zeros(G)
    for f in ("W", "M", "O"):
        idx = lab == f
        LL[np.ix_(idx, idx)] += lam ** 2 * Em2
        L2[idx] += lam ** 2 * Em2
    L2[lab == "H"] += lamH ** 2                                          # E[l^2] for H genes
    SW = np.outer(mu, mu) * (1.0 + GG) * LL
    dW = mu ** 2 * (1.0 + np.diag(GG)) * ((1.0 + L2) * (1.0 + inv_phi) - 1.0)
    np.fill_diagonal(SW, dW)
    SB = np.outer(mu, mu) * GG
    with np.errstate(invalid="ignore", divide="ignore"):
        rW = SW / np.sqrt(np.outer(dW, dW))
        dB = np.diag(SB)
        rB = np.where(np.outer(dB, dB) > 0, SB / np.sqrt(np.outer(dB, dB)), 0.0)
    # population mean of per-donor Fisher-z latent correlation (the CS-CORE-per-donor + meta target)
    rng = np.random.default_rng(seed)
    m = rng.uniform(1 - st["het_h"], 1 + st["het_h"], n_mc) if st["het"] == "mag" else np.ones(1)
    A = lambda l2, ip: (1.0 + l2) * (1.0 + ip) - 1.0
    zbar = np.zeros((G, G))
    for f in ("W", "M", "O"):
        idx = np.where(lab == f)[0]
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                j, k = idx[a], idx[b]
                l2 = (lam * m) ** 2
                rd = l2 / np.sqrt(A(l2, inv_phi[j]) * A(l2, inv_phi[k]))
                zbar[j, k] = zbar[k, j] = np.mean(np.arctanh(rd))
    iu = np.triu_indices(G, 1)
    li, lk = lab[iu[0]], lab[iu[1]]
    same = li == lk
    cls = np.full(len(li), "NN", dtype=object)
    for P, name in (("W", "WW"), ("B", "BB"), ("M", "MM"), ("H", "HH")):
        cls[same & (li == P)] = name
    o_pair = same & (li == "O")
    osg = panel["osign"]
    cls[o_pair & (osg[iu[0]] == osg[iu[1]])] = "OO_same"
    cls[o_pair & (osg[iu[0]] != osg[iu[1]])] = "OO_opp"
    return dict(SW=SW, SB=SB, rW=rW, rB=rB, zbarW=zbar, iu=iu, cls=cls,
                altW=np.abs(SW[iu]) > 0, altB=np.abs(SB[iu]) > 0)


def make_nd(st, rng):
    if st["nd_mode"] == "fixed":
        return np.full(st["D"], int(st["nd"]))
    x = rng.lognormal(np.log(st["nd_median"]), st["nd_logsd"], st["D"])
    return np.clip(np.round(x), st["nd_min"], None).astype(int)


def simulate(panel, st, seed):
    """One replicate. Returns counts Y (N x G int32), s_true, s_hat (observed total UMI),
    donor_idx (list of row arrays), theta (N x G float32) and donor-level draws."""
    rng = np.random.default_rng(seed)
    G, lab, mu, phi = panel["G"], panel["labels"], panel["mu"], panel["phi"]
    D = st["D"]; nd = make_nd(st, rng)
    sa, se = FITTED["sigma_a"] * st["sa_mult"], FITTED["sigma_e"] * st["se_mult"]
    a_d = rng.normal(0.0, sa, D)
    Gm = _between_loading(panel, st["g"])
    u = {h: rng.uniform(-SQ3, SQ3, D) for h in BETWEEN_PROGS}
    mag = {f: (rng.uniform(1 - st["het_h"], 1 + st["het_h"], D) if st["het"] == "mag" else np.ones(D))
           for f in ("W", "M", "O")}
    idx = {f: np.where(lab == f)[0] for f in WITHIN_PROGS}
    bg = max(1.0 - mu.sum(), 0.5)
    Ys, Th, S, Sh, rows, r0 = [], [], [], [], [], 0
    for d in range(D):
        n = int(nd[d])
        beta = 1.0 + sum(Gm[h] * u[h][d] for h in BETWEEN_PROGS)
        w = np.zeros((n, G))
        zW = None
        for f in WITHIN_PROGS:
            z = rng.uniform(-SQ3, SQ3, n)
            if f == "W":
                zW = z
            if f == "H":
                load = st["lam_H"] * rng.choice([-1.0, 1.0], size=len(idx[f]))
            else:
                load = np.full(len(idx[f]), st["lam"] * mag[f][d])
            w[:, idx[f]] += np.outer(z, load)
        xi = np.ones((n, G))
        fin = np.isfinite(phi)
        if fin.any():
            xi[:, fin] = rng.gamma(phi[fin], 1.0 / phi[fin], size=(n, int(fin.sum())))
        theta = mu[None, :] * beta[None, :] * (1.0 + w) * xi
        s = np.exp(FITTED["m_L"] + a_d[d] + rng.normal(0.0, se, n))
        if st["kappa"] != 0.0:
            s = s * np.exp(st["kappa"] * zW)
        Y = rng.poisson(s[:, None] * theta)
        T = rng.poisson(s * bg)
        Ys.append(Y.astype(np.int32)); Th.append(theta.astype(np.float32)); S.append(s)
        Sh.append(Y.sum(1) + T); rows.append(np.arange(r0, r0 + n)); r0 += n
    return dict(Y=np.vstack(Ys), theta=np.vstack(Th), s_true=np.concatenate(S),
                s_hat=np.concatenate(Sh).astype(np.float64), donor_idx=rows, nd=nd, a_d=a_d)


def numeric_truth_check(panel, st, D=4000, n=60, seed=99):
    """Independent numerical check of the closed forms: simulate many donors on the latent scale
    and compare empirical E_D[Cov(theta|D)] and Cov_D(E[theta|D]) with truth()."""
    st2 = dict(st); st2["D"] = D; st2["nd_mode"] = "fixed"; st2["nd"] = n
    sim = simulate(panel, st2, seed)
    th = sim["theta"].astype(np.float64)
    G = panel["G"]
    SW = np.zeros((G, G)); means = []
    for r in sim["donor_idx"]:
        x = th[r]; xc = x - x.mean(0)
        SW += xc.T @ xc / (len(r) - 1)
        means.append(x.mean(0))
    SW /= D
    M = np.array(means)
    SB_emp = np.cov(M, rowvar=False) - SW / n          # finite-cell correction (equal n)
    tr = truth(panel, st)
    return dict(SW_emp=SW, SB_emp=SB_emp, SW=tr["SW"], SB=tr["SB"])
