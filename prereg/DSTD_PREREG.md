# Design-only standardization (LCPd, LCPXd): frozen plan (28 Oct 2026, 18:15 PDT)

## What was seen before this plan (exploratory; disclosed)
- Interim (~45/100 replicates) XS_PREREG simulation results: LCPX loses most power and is not calibrated at ~15 cells.
- On ONE PsychAD pseudo-donor reassignment (ADAM_GPNMB, POP Part A k = 0, seed 41,000,002), D = 25: LCPX type-I 0.068 / 0.114
  (all / top-1%). Diagnostic on the same reassignment (all 744 pseudo-donors): precision-standardized means R_dj correlate with log n_d
  (median over genes 0.155; half-split 0.201), raw ratio means do not (-0.000), log-CPM means do (0.102). Replacing the donor's own
  sampling-variance estimate by a design-only one removed the correlation (-0.002) and reduced the mean null pair correlation from
  0.040 to 0.008.
- Mechanism (hypothesis): psi_hat_dj uses donor d's own nuclei, so R_dj = (m_dj - mu_j)/sqrt(tau2_j + psi_hat_dj) is a nonlinear
  function of the donor's data whose mean depends on n_d; all genes share n_d, inducing a spurious between-donor covariance that
  does not shrink with D. log-CPM pseudobulk has the analogous Jensen effect.
That reassignment (PsychAD P3, k = 0) is excluded from every LCPd/LCPXd evaluation below.

## Definitions
- Pooled gene parameters: mu_j = sum_d Y_dj / sum_d S_d; sigma2_j = mean_d max((psi_hat_dj - m_dj/S_d)/kappa_d, 0) clipped at 0.
- Design-only sampling variance: psi*_dj = kappa_d sigma2_j + mu_j / S_d (depends on the donor only through kappa_d, S_d).
- tau2_j by Paule-Mandel with psi*; R*_dj = (m_dj - mu'_j)/sqrt(tau2_j + psi*_dj), mu'_j the PM mean.
- **LCPd** = LCP with R*, v* = tau2 + psi* in place of R, v (leak term kappa_d C_d,jk / sqrt(v*_dj v*_dk); same scaling and permutation).
- **LCPXd** = LCPX with R* computed within each half (kappa, S of the half; sigma2_j pooled over that half's donors).
- **PBlin** (diagnostic comparator): Pearson correlation of raw ratio means m (no log), Fisher z.
- p-values as in the corresponding earlier plans (B = 10^5 permutations; for D >= 200 the exact-moment normal approximation).

## Evaluation (all on data not used to define LCPd)
E1. Simulations: the FIX_PREREG designs and seeds (F_D16, F_D24, F_D40, F_S1, F_S2, F_D16_n100; R = 100), paired with LCP, LCP2,
    LCPX, PB already computed. FDR, power (MCSE), weak-null type-I (QQ_cross).
E2. SEA-AD real-cell null: the FIX_PREREG reassignments (S1, S2; k = 0..4).
E3. PsychAD Part A (POP_PREREG design): P1 and P2 k = 0, 1, 2; P3 k = 1, 2, 3.

## Decisions
**Adopt LCPd in place of LCP** if all hold:
 (a) LCPd FDR <= 0.05 + 2 MCSE in F_D24, F_D40, F_S2, F_D16_n100;
 (b) LCPd holds the SEA-AD real-cell null (FIX_PREREG criterion) in S1 and S2;
 (c) LCPd power >= 0.9 x LCP power in F_D24 and F_D40;
 (d) in E3, for every subtype, at D = all: LCPd mean BH calls <= 1, all-pair type-I <= 0.055 and top-1% type-I <= 0.05 + 3 SE.
**Claim "design-induced bias is a distinct artifact"** if, in E3 at D = all, LCP's all-pair type-I exceeds 0.055 in at least two
subtypes while LCPd's does not.
Same decisions are evaluated for LCPXd (as a secondary method). All outcomes are reported.
