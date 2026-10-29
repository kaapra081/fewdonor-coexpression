# XU: exact within-donor U-statistic with design-only weights, studentized across donors. Frozen plan (29 Oct 2026, 10:05 PDT)

## Motivation and what was seen before this plan
XDS (cross-split, design-only standardization) held every null tested (POP_RESULTS.md) but uses half of each donor's
nuclei per gene and loses power. XU uses all ordered pairs of distinct nuclei instead of one split and replaces the donor
permutation by a studentized test, whose variance is estimated from the donors themselves (so within-donor fourth-moment
terms are included). XU was designed after seeing all XS/DSTD/POP results; it is evaluated only on fresh seeds below.

## Definition (per donor d, pair j<k)
- u_djk = (Y_dj+ Y_dk+ - sum_c Y_cj Y_ck) / (S_d^2 - sum_c s_c^2): exact unbiased estimate of theta_dj theta_dk given the donor
  (cells i.i.d. within donor, E[Y_cj | cell] = s_c lambda_cj, independent Poisson-type measurement across genes).
- m_dj = Y_dj+ / S_d; leave-one-donor-out equal-weight centre c_dj = mean_{d' != d} m_d'j.
- q_djk = w_djk (u_djk - c_dk m_dj - c_dj m_dk + c_dj c_dk), w_djk = 1/sqrt(v*_dj v*_dk), with v* = tau2 + psi* the
  design-only variance of DSTD_PREREG (pooled sigma2_j, mu_j; Paule-Mandel tau2_j).
- t_jk = mean_d q / (sd_d(q)/sqrt(D)); two-sided p from N(0,1). (For D < 100, also reported with t_{D-1}; out of scope for decisions.)
- Null tested: Sigma_B,jk = 0 (weak null), given that donor effects are independent of the design.

## Evaluation (fresh seeds only)
E1. PsychAD pseudo-donor null (POP design): P1, P2 with k = 3, 4, 5; P3 with k = 4, 5, 6; D in {100, 200, 400, all}.
    Same methods as before also re-run on these reassignments (PB, LCP, XDS) for a paired comparison.
E2. SEA-AD real-cell null: new reassignments k = 5..9 (seeds 31,000,000 + k), S1 and S2; XU with t_{D-1}.
E3. Large-cohort simulations (new): designs L200, L500, L1000 = D40 settings with D in {200, 500, 1000}, nd_median 20, nd_min 4,
    held-out panel (seed 4112026), data seeds 9,500,000 + 100,003 k + rep, R = 50 each. PB, LCP, XDS, XU with normal reference
    (exact-moment normal approximation for permutation statistics). FDR (BH 0.05), power, weak-null type-I (QQ_cross), MCSE.
E4. The six FIX designs (existing paired data seeds, R = 100), XU with t_{D-1}: reported, no decision.

## Decisions
- "XU holds the real-cell null at population scale": in E1 at every D >= 200 for every subtype, mean BH calls <= 1,
  all-pair type-I <= 0.055 and top-1% type-I <= 0.05 + 3 SE.
- "XU is the recommended large-cohort test": the above AND E3 FDR <= 0.05 + 2 MCSE in L200, L500, L1000 AND XU power >= 1.5 x XDS
  power in L500 and L1000.
- "Weak-null validity": XU weak-null type-I <= 0.05 + 2 MCSE in L200, L500, L1000.
- If XU fails, XDS remains the recommended large-cohort test and the failure is reported.
