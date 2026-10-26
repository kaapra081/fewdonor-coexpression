"""Candidate between-donor co-expression tests (recomb_v2 development).

Notation (per donor d, gene j): S_d = sum_i s_i, m_dj = sum_i Y_ij / S_d (ratio mean, unbiased for mu_dj),
kappa_d = sum_i s_i^2 / S_d^2, x_ij = Y_ij / s_i, C_d = sample covariance of x (unbiased for Cov_d(theta)
off the diagonal), sig2_dj = Var_i(x_ij) - mean_i(x_ij/s_i) (Poisson-corrected latent variance),
psi_dj = kappa_d sig2_dj + m_dj/S_d  (sampling variance of m_dj), l_djk = kappa_d C_d,jk (sampling covariance:
the 'leak'). Per gene: Paule-Mandel tau2_j (= Sigma_B,jj), v_dj = tau2_j + psi_dj, RE-weighted mean.
"""
import numpy as np
from scipy import stats


def donor_stats(Y, s, donor_idx, min_cells=4, split_seed=None):
    out = dict(m=[], psi=[], C=[], kappa=[], S=[], n=[])
    if split_seed is not None:
        rng = np.random.default_rng(split_seed)
        out.update(mA=[], mB=[], psiA=[], psiB=[])
    for rows in donor_idx:
        if len(rows) < min_cells:
            continue
        sd = s[rows]; Yd = Y[rows].astype(np.float64)
        st = _one(Yd, sd)
        for k in ("m", "psi", "C", "kappa", "S"):
            out[k].append(st[k])
        out["n"].append(len(rows))
        if split_seed is not None:
            perm = rng.permutation(len(rows)); h = len(rows) // 2
            a, b = perm[:h], perm[h:]
            sa, sb = _one(Yd[a], sd[a], cov=False), _one(Yd[b], sd[b], cov=False)
            out["mA"].append(sa["m"]); out["psiA"].append(sa["psi"])
            out["mB"].append(sb["m"]); out["psiB"].append(sb["psi"])
    return {k: np.array(v) for k, v in out.items()}


def _one(Yd, sd, cov=True):
    S = sd.sum(); n = len(sd)
    x = Yd / sd[:, None]
    m = Yd.sum(0) / S
    kap = (sd ** 2).sum() / S ** 2
    xc = x - x.mean(0)
    var = (xc ** 2).sum(0) / (n - 1)
    sig2 = var - (x / sd[:, None]).mean(0)
    psi = kap * np.maximum(sig2, 0.0) + m / S
    psi = np.maximum(psi, 1e-30)
    o = dict(m=m, psi=psi, kappa=kap, S=S)
    if cov:
        o["C"] = xc.T @ xc / (n - 1)
    return o


def paule_mandel(m, psi, iters=60):
    """Vectorised Paule-Mandel tau^2 per gene (columns) with known sampling variances psi (D x G)."""
    D = m.shape[0]
    lo = np.zeros(m.shape[1]); hi = np.maximum(m.var(0, ddof=1) * 10, 1e-30)
    def Q(t):
        w = 1.0 / (t[None, :] + psi)
        mu = (w * m).sum(0) / w.sum(0)
        return (w * (m - mu) ** 2).sum(0)
    q0 = Q(lo)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        big = Q(mid) > D - 1
        lo = np.where(big, mid, lo); hi = np.where(big, hi, mid)
    tau2 = np.where(q0 <= D - 1, 0.0, 0.5 * (lo + hi))
    w = 1.0 / (tau2[None, :] + psi)
    mu = (w * m).sum(0) / w.sum(0)
    return tau2, mu, w / w.sum(0)


def _pairs(G):
    return np.triu_indices(G, 1)


def bh(p, q=0.05):
    p = np.asarray(p, float); ok = np.isfinite(p)
    rej = np.zeros(len(p), bool)
    if ok.sum() == 0:
        return rej
    ps = np.sort(p[ok]); m = len(ps)
    k = np.where(ps <= q * np.arange(1, m + 1) / m)[0]
    if len(k):
        rej[ok] = p[ok] <= ps[k.max()]
    return rej


def cauchy_comb(p1, p2):
    p1 = np.clip(p1, 1e-300, 1 - 1e-16); p2 = np.clip(p2, 1e-300, 1 - 1e-16)
    t = 0.5 * np.tan((0.5 - p1) * np.pi) + 0.5 * np.tan((0.5 - p2) * np.pi)
    return stats.cauchy.sf(t)


# ---------------------------------------------------------------- (B) conditional leak-offset regression
def cond_test(ds_, iu, floor=True, eb=False, direction="both"):
    m, psi, C, kap = ds_["m"], ds_["psi"], ds_["C"], ds_["kappa"]
    D, G = m.shape
    tau2, mu, _ = paule_mandel(m, psi)
    v = tau2[None, :] + psi                                  # D x G total variance of donor means
    a, b = iu
    L = kap[:, None] * C[:, a, b]                            # D x P leak (sampling covariance)
    res = {}
    for name, (j, k) in (("jk", (a, b)), ("kj", (b, a))):
        x = m[:, j] - mu[None, j]                            # centred regressor gene
        z = x / v[:, j]
        y = m[:, k] - L * z                                  # leak offset removed
        u = np.maximum(v[:, k] - L ** 2 / v[:, j], 1e-3 * v[:, k])
        w = 1.0 / u
        sw = w.sum(0); zb = (w * z).sum(0) / sw; yb = (w * y).sum(0) / sw
        Szz = (w * (z - zb) ** 2).sum(0); Szy = (w * (z - zb) * (y - yb)).sum(0)
        beta = Szy / Szz
        rss = (w * (y - yb - beta * (z - zb)) ** 2).sum(0)
        q = rss / (D - 2)
        if floor:
            q = np.maximum(q, 1.0)
        se = np.sqrt(q / Szz)
        t = beta / se
        p = 2 * stats.t.sf(np.abs(t), df=D - 2)
        res[name] = dict(est=beta, se=se, p=p, t=t)
    p = cauchy_comb(res["jk"]["p"], res["kj"]["p"]) if direction == "both" else res[direction]["p"]
    est = 0.5 * (res["jk"]["est"] + res["kj"]["est"])
    se = 0.5 * (res["jk"]["se"] + res["kj"]["se"])
    return dict(est=est, se=se, p=p, df=np.full(len(a), D - 2.0), p_jk=res["jk"]["p"], p_kj=res["kj"]["p"],
                est_jk=res["jk"]["est"], se_jk=res["jk"]["se"], tau2=tau2)


# ---------------------------------------------------------------- (A) leak-corrected weighted score, exact Gaussian null variance
def lcw_test(ds_, iu, weights="opt"):
    m, psi, C, kap = ds_["m"], ds_["psi"], ds_["C"], ds_["kappa"]
    D, G = m.shape
    tau2, mu, al = paule_mandel(m, psi)                       # al: D x G centring weights (sum 1)
    v = tau2[None, :] + psi
    a, b = iu
    L = kap[:, None] * C[:, a, b]                             # D x P
    x = m[:, a] - mu[None, a]; y = m[:, b] - mu[None, b]
    if weights == "opt":
        w = 1.0 / (v[:, a] * v[:, b])
    else:
        w = np.ones_like(x)
    alj, alk = al[:, a], al[:, b]
    A = (alj * alk).sum(0)
    bias = L * (1 - alj - alk) + (alj * alk * L).sum(0)[None, :]
    num = (w * (x * y - bias)).sum(0)
    den = (w * (1 - alj - alk + A[None, :])).sum(0)
    est = num / den
    # exact Gaussian variance of sum_d w_d x_d y_d under H0 (Isserlis), centring included
    P = len(a)
    I = np.eye(D)
    Pa = I[None] - np.transpose(alj)[:, None, :]              # P x D x D : row d -> e_d - alpha
    Pb = I[None] - np.transpose(alk)[:, None, :]
    Vj = np.transpose(v[:, a]); Vk = np.transpose(v[:, b]); Lt = np.transpose(L)
    Sxx = np.einsum("pde,pe,pfe->pdf", Pa, Vj, Pa)
    Syy = np.einsum("pde,pe,pfe->pdf", Pb, Vk, Pb)
    Sxy = np.einsum("pde,pe,pfe->pdf", Pa, Lt, Pb)
    Wt = np.transpose(w)
    # Var(sum_d w_d x_d y_d) = sum_{d,e} w_d w_e (Sxx_de Syy_de + Sxy_de Sxy_ed)
    var = (np.einsum("pd,pe,pde,pde->p", Wt, Wt, Sxx, Syy)
           + np.einsum("pd,pe,pde,ped->p", Wt, Wt, Sxy, Sxy))
    zst = num / np.sqrt(var)
    p = 2 * stats.norm.sf(np.abs(zst))
    return dict(est=est, se=np.sqrt(var) / np.abs(den), p=p, z=zst)


# ---------------------------------------------------------------- (C) cell-split conditional regression
def split_test(ds_, iu, floor=True):
    mA, mB, pA, pB = ds_["mA"], ds_["mB"], ds_["psiA"], ds_["psiB"]
    D = mA.shape[0]
    tA, muA, _ = paule_mandel(mA, pA); tB, muB, _ = paule_mandel(mB, pB)
    vA = tA[None, :] + pA; vB = tB[None, :] + pB
    a, b = iu
    out = []
    for (X, VX, muX, Yv, VY, j, k) in ((mA, vA, muA, mB, vB, a, b), (mB, vB, muB, mA, vA, a, b)):
        z = (X[:, j] - muX[None, j]) / VX[:, j]
        y = Yv[:, k]; w = 1.0 / VY[:, k]
        sw = w.sum(0); zb = (w * z).sum(0) / sw; yb = (w * y).sum(0) / sw
        Szz = (w * (z - zb) ** 2).sum(0); beta = (w * (z - zb) * (y - yb)).sum(0) / Szz
        q = (w * (y - yb - beta * (z - zb)) ** 2).sum(0) / (D - 2)
        if floor:
            q = np.maximum(q, 1.0)
        se = np.sqrt(q / Szz)
        out.append((beta, 2 * stats.t.sf(np.abs(beta / se), D - 2)))
    return dict(est=0.5 * (out[0][0] + out[1][0]), p=cauchy_comb(out[0][1], out[1][1]), p1=out[0][1])


# ---------------------------------------------------------------- comparators
def cmb_mom(ds_, iu):
    M = ds_["m"]; kap = ds_["kappa"]; C = ds_["C"]; K = M.shape[0]
    a, b = iu
    Mj, Mk = M[:, a], M[:, b]
    KC = kap[:, None] * C[:, a, b]
    SP = (Mj * Mk).sum(0); Sj = Mj.sum(0); Sk = Mk.sum(0); SKC = KC.sum(0)
    est = (SP - Sj * Sk / K) / (K - 1) - SKC / K
    loo = np.empty((K, len(a)))
    for e in range(K):
        sp = SP - Mj[e] * Mk[e]; sj = Sj - Mj[e]; sk = Sk - Mk[e]
        loo[e] = (sp - sj * sk / (K - 1)) / (K - 2) - (SKC - KC[e]) / (K - 1)
    se = np.sqrt((K - 1) / K * ((loo - loo.mean(0)) ** 2).sum(0))
    p = 2 * stats.t.sf(np.abs(est / se), df=K - 1)
    return dict(est=est, se=se, p=p)


def pearson_fisher(M, iu):
    K = M.shape[0]
    Mc = M - M.mean(0); cov = Mc.T @ Mc / (K - 1); sd = np.sqrt(np.diag(cov))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = cov[iu] / (sd[iu[0]] * sd[iu[1]])
    p = 2 * stats.norm.sf(np.abs(np.arctanh(np.clip(r, -0.999999, 0.999999)) * np.sqrt(K - 3)))
    return dict(est=r, p=p)


def pb_logcpm(ds_, iu):
    return pearson_fisher(np.log1p(1e6 * ds_["m"]), iu)


# ================================================================ main candidate: LCSP
# Leak-Corrected, precision-Standardized, (optionally Studentized) donor-Permutation test.
def standardize(ds_):
    m, psi = ds_["m"], ds_["psi"]
    tau2, mu, _ = paule_mandel(m, psi)
    v = tau2[None, :] + psi
    return (m - mu[None, :]) / np.sqrt(v), v, tau2, mu


def _perm_stats(R, R2, a, b, perms):
    """For a batch of permutations (bb x D int), return sums S1 = sum_d R[pi(d),a] R[d,b] and
    S2 = sum_d R[pi(d),a]^2 R[d,b]^2 for all pairs (bb x P)."""
    Rp = R[perms]; R2p = R2[perms]
    T1 = np.matmul(Rp.transpose(0, 2, 1), R)[:, a, b]
    T2 = np.matmul(R2p.transpose(0, 2, 1), R2)[:, a, b]
    return T1, T2


def _tstat(S1, S2, D):
    # one-sample t statistic of the D products: mean / (sd / sqrt(D)), computed from sums
    var = np.maximum((S2 - S1 ** 2 / D) / (D - 1), 1e-300)
    return (S1 / D) / np.sqrt(var / D)


def lcsp_test(ds_, iu, B=20000, seed=0, studentize=True, leak=True, refine_thr=5e-4, B2=400_000,
              batch=1000, R_override=None, L_override=None, ci_level=0.95):
    if R_override is None:
        R, v, tau2, mu = standardize(ds_)
    else:
        R, v = R_override
        tau2 = None
    D, G = R.shape
    a, b = iu
    if L_override is not None:
        Lt = L_override
    elif leak:
        Lt = ds_["kappa"][:, None] * ds_["C"][:, a, b] / np.sqrt(v[:, a] * v[:, b])
    else:
        Lt = np.zeros((D, len(a)))
    q = R[:, a] * R[:, b] - Lt                                 # D x P leak-corrected standardized products
    S1o = q.sum(0); S2o = (q ** 2).sum(0)
    if studentize:
        obs = _tstat(S1o, S2o, D)
    else:
        obs = S1o * np.sqrt(D / (D + (Lt ** 2).sum(0)))
    R2 = R ** 2
    rng = np.random.default_rng(seed)
    cnt = np.zeros(len(a))
    aobs = np.abs(obs) - 1e-12
    for s0 in range(0, B, batch):
        bb = min(batch, B - s0)
        perms = np.argsort(rng.random((bb, D)), axis=1)
        T1, T2 = _perm_stats(R, R2, a, b, perms)
        st = _tstat(T1, T2, D) if studentize else T1
        cnt += (np.abs(st) >= aobs).sum(0)
    p = (cnt + 1) / (B + 1)
    idx = np.where(p <= refine_thr)[0]
    if len(idx):
        PERM = np.argsort(rng.random((B2, D)), axis=1).astype(np.int32)
        for i in idx:
            x = R[:, a[i]]; y = R[:, b[i]]
            T1 = (x[PERM] * y[None, :]).sum(1)
            if studentize:
                T2 = ((x ** 2)[PERM] * (y ** 2)[None, :]).sum(1)
                st = _tstat(T1, T2, D)
            else:
                st = T1
            p[i] = ((np.abs(st) >= aobs[i]).sum() + 1) / (B2 + 1)
    # point estimate and t-interval for Sigma_B,jk: q_d estimates Sigma_B,jk * w_d, w_d = 1/sqrt(v_dj v_dk)
    w = 1.0 / np.sqrt(v[:, a] * v[:, b])
    est = S1o / w.sum(0)
    se_q = np.sqrt(np.maximum((S2o - S1o ** 2 / D) / (D - 1), 0) / D)
    tq = stats.t.ppf(0.5 + ci_level / 2, D - 1)
    lo = (S1o / D - tq * se_q) / w.mean(0); hi = (S1o / D + tq * se_q) / w.mean(0)
    return dict(p=p, stat=obs, est=est, lo=lo, hi=hi)


def bound_test(sim, tr, iu, seed, **kw):
    """Noise-limited per-pair BOUND: true donor means + independent Gaussian noise with the TRUE sampling
    variances (no leak), standardized with the TRUE total variances, same permutation calibration.
    Approximates the best power any per-pair test can attain at this noise level."""
    rng = np.random.default_rng(seed + 77)
    mu_t, psi_t = sim["mu_true"], sim["psi_true"]
    mstar = mu_t + np.sqrt(psi_t) * rng.normal(size=mu_t.shape)
    sB = np.diag(tr["SB"])
    v = sB[None, :] + psi_t
    mbar = mu_t.mean(0) * 0 + (mu_t.mean(0))   # population mean unknown in finite D: use sample-free truth below
    R = (mstar - mstar.mean(0)) / np.sqrt(v)
    D = R.shape[0]
    return lcsp_test(None, iu, seed=seed, R_override=(R, v), L_override=np.zeros((D, len(iu[0]))), **kw)


def oracle_perm(sim, iu, seed, **kw):
    M = sim["mu_true"]
    R = (M - M.mean(0)) / M.std(0, ddof=1)
    D = R.shape[0]
    return lcsp_test(None, iu, seed=seed, R_override=(R, np.ones_like(R)), L_override=np.zeros((D, len(iu[0]))), **kw)


def e5_logstd(Y, s, donor_idx, iu):
    """E5 (round-3 comparator): precision-standardized donor means of log1p(CP10k), Pearson + Fisher z."""
    A = np.log1p(Y.astype(np.float64) / s[:, None] * 1e4)
    n = np.array([len(r) for r in donor_idx], float)
    M = np.vstack([A[r].mean(0) for r in donor_idx]); V = np.vstack([A[r].var(0, ddof=1) for r in donor_idx])
    s2w = ((n[:, None] - 1) * V).sum(0) / (n.sum() - len(n))
    wts = 1.0 / (s2w[None, :] / n[:, None]); mm = (wts * M).sum(0) / wts.sum(0)
    s2b = np.maximum(0.0, ((M - M.mean(0)) ** 2).sum(0) / (len(n) - 1) - (s2w[None, :] / n[:, None]).mean(0))
    U = (M - mm) / np.sqrt(s2b[None, :] + s2w[None, :] / n[:, None])
    return pearson_fisher(U, iu)


def eq_interval(ds_, iu, level=0.95):
    """Interval for Sigma_B,jk: equal-weight leak-corrected (BCES-type) unbiased estimator
    est = S_B,jk - mean_d(kappa_d C_d,jk), with its exact Gaussian finite-sample variance (Isserlis, centring
    included) evaluated at plug-in values, and a t_{D-1} reference."""
    m, psi, C, kap = ds_["m"], ds_["psi"], ds_["C"], ds_["kappa"]; D, G = m.shape
    tau2, mu, _ = paule_mandel(m, psi); v = tau2[None, :] + psi
    a, b = iu; L = kap[:, None] * C[:, a, b]
    x = m[:, a] - m[:, a].mean(0); y = m[:, b] - m[:, b].mean(0)
    est = (x * y).sum(0) / (D - 1) - L.mean(0)
    # Var of sum_d x_d y_d with x = P m_j, y = P m_k, P = I - J/D; Cov blocks diag(v_j), diag(v_k), diag(c)
    c = L + est[None, :]
    P_ = np.eye(D) - 1.0 / D
    vj, vk = v[:, a], v[:, b]
    # (P diag(u) P)_de = u_d delta_de - (u_d + u_e)/D + sum(u)/D^2
    def blk(u):
        su = u.sum(0)
        return (np.eye(D)[None] * np.transpose(u)[:, :, None]
                - (np.transpose(u)[:, :, None] + np.transpose(u)[:, None, :]) / D + (su / D ** 2)[:, None, None])
    Sxx, Syy, Sxy = blk(vj), blk(vk), blk(c)
    var = (Sxx * Syy).sum((1, 2)) + (Sxy * np.transpose(Sxy, (0, 2, 1))).sum((1, 2))
    se = np.sqrt(var) / (D - 1)
    tq = stats.t.ppf(0.5 + level / 2, D - 1)
    return dict(est=est, lo=est - tq * se, hi=est + tq * se, se=se)


# ================================================================ amended interval (change log 1, DEV_PREREG.md)
def _bces_parts(m, LC, iu):
    """Per-donor inputs for the equal-weight leak-corrected (BCES-type) estimator."""
    return m[:, iu[0]], m[:, iu[1]], LC


def _bces_est_se(X, Y, L):
    """X, Y: (..., D, P) donor means; L: (..., D, P) leak terms. Returns unbiased estimate
    S_B,jk - mean_d L_d and a sandwich SE (SD of per-donor centred products / sqrt(D))."""
    D = X.shape[-2]
    xc = X - X.mean(-2, keepdims=True); yc = Y - Y.mean(-2, keepdims=True)
    q = xc * yc * D / (D - 1) - L
    est = q.mean(-2)
    se = q.std(-2, ddof=1) / np.sqrt(D)
    return est, se


def boot_t_interval(X, Y, L, B=999, seed=0, level=0.95, min_distinct=4):
    """Donor bootstrap-t interval for Sigma_B,jk. Resamples donors (with their leak terms) with replacement;
    resamples with fewer than `min_distinct` distinct donors are redrawn."""
    D = X.shape[0]
    est, se = _bces_est_se(X, Y, L)
    rng = np.random.default_rng(seed)
    tb = np.empty((B, X.shape[1]))
    i = 0
    while i < B:
        idx = rng.integers(0, D, D)
        if len(np.unique(idx)) < min_distinct:
            continue
        e_b, s_b = _bces_est_se(X[idx], Y[idx], L[idx])
        tb[i] = (e_b - est) / np.maximum(s_b, 1e-300)
        i += 1
    lo_q, hi_q = np.quantile(tb, [0.5 + level / 2, 0.5 - level / 2], axis=0)
    return dict(est=est, se=se, lo=est - lo_q * se, hi=est - hi_q * se)
