# Cross-split LCP (LCPX): frozen plan (written 28 Oct 2026 17:28 PDT, before any LCPX computation)

## Motivation
LCP and LCP2 estimate and subtract the within-donor leak. Both fail the frozen real-cell null for the top-1% within-correlated
pairs (FIX_RESULTS.md), and LCP is not calibrated at 16 donors x ~40 cells or ~15 cells per donor. LCPX avoids estimating the
leak: within each donor the nuclei are split at random into two disjoint halves A and B; gene j's donor mean comes from A and
gene k's from B. Given the donor, the halves are independent, so E[r^A_dj r^B_dk] contains no within-donor term and, under the
strong null, the per-donor product is a product of independent factors whose variance the donor permutation reproduces.

## Definitions
- Split: for each donor with n_d >= 4, a uniformly random subset A of floor(n_d/2) nuclei; B = the rest. Split seed = data seed + 7.
- Per half h in {A, B}: m^h, psi^h exactly as donor_stats (ratio mean, kappa*sigma^2 + m/S); tau2^h, mu^h by Paule-Mandel;
  R^h_dj = (m^h_dj - mu^h_j) / sqrt(tau2^h_j + psi^h_dj).
- Pair (j<k): q_d = R^A_dj R^B_dk.
- **LCPX** (primary): T = sum_d q_d; two-sided donor-permutation p-value, permuting the rows of R^B (B = 10^5).
- **LCPX-S** (secondary, studentized): t = mean_d q_d / (sd_d(q_d)/sqrt(D)); same permutations applied to t (DiCiccio-Romano studentized permutation).
- No tuning after this point. One split per data set (no split aggregation).

## Evaluation (same data as FIX_PREREG Item 2: identical seeds, so the LCP / LCP2 / PB p-values in fix/runs are paired)
- Simulations: F_D16, F_D24, F_D40, F_S1, F_S2, F_D16_n100; R = 100 each; held-out panel; BH q = 0.05.
  Report FDR, power (with MCSE) and weak-null type-I (classes QQ_cross, QQ_same; alpha = 0.05) for LCPX, LCPX-S, LCP, LCP2, PB.
- Real-cell null: the same 10 pseudo-donor reassignments as FIX_PREREG Item 1 (seeds 31,000,000 + k); same criterion
  ("holds" = all-pair type-I <= 0.055 and top-1% within-correlated type-I <= 0.05 + 3 SE).

## Decisions
- **Adopt LCPX as the recommended test** if all hold:
  (a) LCPX FDR <= 0.05 + 2 MCSE in all six designs;
  (b) LCPX holds the real-cell null in S1 and S2;
  (c) LCPX power >= 0.8 x LCP power in F_D24 and F_D40.
  If (a) and (b) hold but (c) fails, LCPX is reported as the calibrated-everywhere but lower-powered option, and LCP stays primary.
  If (a) or (b) fails, that is reported as a negative result.
- **Claim "LCPX-S tests Sigma_B = 0 (weak null)"** only if its weak-null type-I <= 0.05 + 2 MCSE and its FDR <= 0.05 + 2 MCSE in
  F_D24, F_D40, F_S1 and F_S2. F_D16 results are reported either way.

## Amendment 1 (28 Oct 17:40 PDT; before any full run; only a one-replicate smoke test at B = 2000 had been run)
The weak null is the class QQ_cross only. QQ_same pairs have Sigma_B != 0 in the simulation truth (altB = True) and are
alternatives. Weak-null type-I is therefore computed over QQ_cross. Nothing else changes. Original hash in XS_PREREG_original.sha256.
