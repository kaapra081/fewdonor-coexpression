# Two fixes: frozen plan (written 28 Oct 2026, before any of the computations below)

## Item 1. Real-cell null with exact truth (pseudo-donors)

**Data.**
- S1: SEA-AD Micro-PVM_1.
- S2: Micro-PVM_3-SEAAD.
- Both use the same nuclei, donors (>= 4 nuclei) and 417-gene panel as RARE_PREREG.md.

**Null construction.**
- Pool the nuclei of a dataset.
- Randomly permute the donor labels across nuclei. This keeps every donor's nucleus count n_d.
- The pseudo-donors are then exchangeable draws from one pooled population. So Σ_B = 0 for every pair, while every real within-cell co-expression program is kept.
- K = 5 independent permutations per dataset; seeds 31,000,000 + k.

**Tests.** PB, LCP (B = 10^5, no refinement) and LCP2 (Item 2), each with BH at q = 0.05.

**Reported, per method.**
- Mean number of BH calls over the K permutations (every call is a false positive).
- Fraction of permutations with at least one call.
- Type-I error at alpha = 0.05 over all pairs.
- Type-I error at alpha = 0.05 over the top 1% of pairs by within-pseudo-donor correlation.

**Decision.**
- "The leak produces spurious pseudobulk discoveries in real cells at natural cell numbers" = YES for a dataset if PB's mean number of calls is >= 50 and >= 10 × LCP's mean number of calls (LCP's counted as at least 1).
- LCP (and LCP2) "hold the real-cell null" if:
  - type-I over all pairs is <= 0.05 + 0.005; and
  - type-I over the top-1% within-correlated pairs is <= 0.05 + 3·SE, where SE is binomial over pairs × permutations.

## Item 2. LCP2: adding the variance of the per-donor leak estimate (Theorem 2(b))

**Definitions.**
- Per donor d and pair jk, z_d = r_dj r_dk − ℓ̃_d, as in LCP.
- Compute delete-one-nucleus jackknife variances, holding v_dj, v_dk, μ̂ fixed:
  - V_d^z for z_d;
  - V_d^r for r_dj r_dk.
- δ̂_d = V_d^z − V_d^r. This estimates Var(ℓ̃_d) − 2 Cov(r_dj r_dk, ℓ̃_d).

**Statistic.**
- T2 = [Σ_d z_d] · sqrt( D / (D + Σ_d ℓ̃_d² + max(Σ_d δ̂_d, 0)) ).
- The reference distribution is the same donor-permutation distribution as LCP. B = 10^5, no refinement, in both cases, with shared permutations.

**Evaluation.**
- New seeds: 9,000,000 + 100,003·k + rep.
- The held-out gene panel (seed 4112026) is used; the development panel is not.
- R = 100 per design.
- Designs:
  - F_D16 (16 donors, ~40 cells);
  - F_D24;
  - F_D40;
  - F_S1 (79 donors, S1's exact cell numbers);
  - F_S2 (83 donors, S2's exact cell numbers);
  - F_D16_n100.
- The non-cell-number settings are those of the development designs.

**Adopt LCP2 in place of LCP** if all three hold:
- (a) LCP2's FDR <= 0.05 + 2·MCSE in F_D24, F_D40, F_S2 and F_D16_n100;
- (b) LCP2's FDR <= 0.05 + 2·MCSE in at least one of F_D16 and F_S1, the two designs where LCP failed;
- (c) LCP2's power >= 0.9 × LCP's power in F_D24 and F_D40.

**Otherwise.** LCP stays the method. LCP2 is reported as a tested and rejected variant.

**Both outcomes are reported.**
