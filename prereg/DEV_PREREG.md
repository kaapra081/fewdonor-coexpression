# Development pre-registration and GO/NO-GO gates (recomb_v2)

Written Mon 26 Oct 2026, 11:40 PDT, **before any formal development replicate was run**. Its sha256 is recorded in `DEV_PREREG.sha256` and in project memory. Any later change is logged at the end with a timestamp and reason, before the affected results are inspected.

## 0. What has already been seen (disclosure)

Exploratory runs (seed ranges 7.0e6, 9.0e6–9.3e6; R = 20–40 per setting) were used to develop the method. They showed: normal-theory between-donor tests (including Pearson on true donor means) are anti-conservative in the BH-relevant tails at 16 donors; a leak-corrected, precision-standardized donor-permutation test is not; studentizing the permutation statistic restores weak-null validity but removes most power at 16 donors; Wald-type 95% intervals for Σ_B,jk undercover (0.81–0.92) at 16 donors with log-normal donor effects, *including on the true donor means* (0.84–0.87). Those runs are not used for any gate. The formal runs below use disjoint seeds.

## 1. Frozen code

`recomb_v2/dev/` (sha256): sim_v2.py 3755e662…; methods_v2.py 6e5d8484…; run_dev.py ddef72f1…; diag_sim.py (copy of round-3, constants only) e03468fd…. Full hashes in `DEV_PREREG.sha256`.

## 2. Estimand and hypotheses

Σ_B,jk = Cov_D(E[θ_j | D], E[θ_k | D]) on the latent relative-expression scale, equal donor weights (round3/estimands_and_method.md, Part A). Tested hypothesis of the primary test: **H0: the donor-level effects of genes j and k are independent** (implies Σ_B,jk = 0). The weak null Σ_B,jk = 0 with dependent donor effects is a declared stress null (QQ_cross).

## 3. Methods (frozen)

- **LCP (primary test, proposed).** Per gene: Paule–Mandel τ̂²_j from donor ratio means m_dj with count-derived sampling variances ψ̂_dj = κ_d σ̂²_dj + m_dj/S_d; standardized deviations r_dj = (m_dj − μ̂_j)/√(τ̂²_j + ψ̂_dj). Per pair: T = [Σ_d (r_dj r_dk − ℓ̃_d)] · √(D / (D + Σ_d ℓ̃_d²)) with leak ℓ̃_d = κ_d Ĉ_d,jk / √(v̂_dj v̂_dk). Two-sided p-value from 20,000 random donor permutations of gene j's standardized vector (refined with 400,000 permutations when p ≤ 5·10⁻⁴). BH at q = 0.05 over all pairs.
- **EQ interval (primary interval, proposed).** Equal-weight leak-corrected unbiased estimator Σ̂_B,jk = S_B,jk − mean_d κ_d Ĉ_d,jk, its exact Gaussian finite-sample variance (Isserlis; centring included) at plug-in values, t_{D−1}, 95%.
- **Comparators:** CMB_mom (round-3 moment correction + donor jackknife, t_{D−1}) and its jackknife interval; CMB_naive (count-scale donor-mean Pearson, Fisher z); PB_logcpm (log-of-sums pseudobulk Pearson); E5 (precision-standardized log1p donor means); LCP_stud (studentized variant, reported); **ORACLE** (same permutation test on the true donor means μ_dj, standardized by their sample SD); **BOUND** (same permutation test on true donor means plus independent Gaussian noise with the *true* sampling variances, no leak, standardized with the *true* total variances: the noise-limited per-pair benchmark).

## 4. Designs (R = 200 replicates each; seeds 2,000,000 + 100,003·design_index + rep)

Panel: sim_v2.make_panel(seed 26102026), G = 80 genes, 3,160 pairs; programs W10 (within-only, λ 0.30), S8 (strong within-only, λ 0.60), B10 (between-only, g 0.30), M8 (both), O8 (both, opposite signs), H10 (sign-heterogeneous within-null), Q10 (between heterogeneity: Q1–Q2 pairs are the weak null), N16; log-normal programs; gene-specific donor effects η_j ~ LogNormal(log 0.19, 0.35) (SEA-AD median between-donor log-SD 0.19). Cells per donor n_d = round(LogNormal(log 40, 0.66)), min 10; per-cell depth as in round 3 (SEA-AD fit).

| Design | D | Other | Role |
|---|---|---|---|
| D8, D12, D16, D24 | 8, 12, 16, 24 | — | brief's few-donor ladder |
| D40 | 40 | — | rare-cell design (amended gate 3) |
| D16_strong | 16 | g = 0.45 | secondary effect size |
| D16_STRESS_depth | 16 | depth–state coupling κ = 0.3, within-donor depth SD × 2 | **stress** (not gated) |

Pair classes: between-null = WW, SS, HH, NN (development nulls, gated), QQ_cross (weak null, stress, reported); between-alternative = BB, MM, OO_same, OO_opp, QQ_same.

## 5. Metrics (all with Monte Carlo standard errors across replicates)

- Type-I per null class at α = 0.05 (unadjusted p), MCSE = SD over replicates of the per-replicate class rate / √R.
- FDR = mean over replicates of the BH false discovery proportion (false = any pair with Σ_B,jk = 0, QQ_cross included), MCSE = SD/√R.
- Power = mean over replicates of the BH discovery fraction among Σ_B,jk ≠ 0 pairs, MCSE = SD/√R.
- Coverage of 95% intervals per pair class (pooled over pairs; MCSE from per-replicate class coverage).
- Runtime per replicate; peak memory recorded separately.

## 6. GO/NO-GO gates (Wed 28 Oct 18:00 PDT)

GO only if all four hold.

1. **Novelty.** A concrete new estimator, test or theorem that survives the skeptical-reviewer paragraph of novelty_v2.md (judged in writing in GO_NOGO.md).
2. **Validity.** For LCP in every non-stress design (D8, D12, D16, D24, D40, D16_strong): FDR ≤ 0.05 + 2·MCSE, and type-I ≤ 0.05 + 3·MCSE in every development null class (WW, SS, HH, NN).
3. **Usefulness (amended 26 Oct with the user's approval, before formal runs).**
   - (a) LCP BH power ≥ 0.8 × BOUND power at D16 and at D40. If BOUND power < 0.05 in a design, that design's criterion is instead LCP power ≥ BOUND power − 2·MCSE of the paired difference.
   - (b) LCP BH power ≥ 0.5 × ORACLE power at D40.
   - The brief's original gate (D16: LCP power ≥ 0.5 or ≥ 0.5 × ORACLE) is reported alongside, labelled "original gate 3".
4. **Intervals.** EQ 95% interval coverage ≥ 0.93 in every pair class of every non-stress design.

If GO: an abstract matching what the paper will claim, without held-out numbers. If NO-GO: stated plainly with the numbers.

## 7. Change log

(none)
