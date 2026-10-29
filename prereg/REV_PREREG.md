# Analyses responding to the internal review: frozen plan (29 Oct 2026, before any of the computations below)

## R1. Run-level error in the complete real-cell null (PsychAD)
Under a complete null every BH call is false, so FDR = P(V >= 1). 100 new pseudo-donor reassignments per subtype (P1-P3, full D;
seeds 48,000,000 + 10 k + subtype index, k = 0..99). Methods: PB, XDS (as in the paper), XDS-M and KX (below).
Per run: number of BH calls V, minimum p, and counts of p below 1e-3, 1e-4, 1e-5 (expected 89, 8.9, 0.89 for ~89,000 pairs).
Decision "XDS controls the FDR at population scale" = YES for a subtype if the fraction of runs with V >= 1 is <= 0.05 + 2 SE (binomial).
Same decision for XDS-M and KX.

## R2. Two baselines/variants (definitions frozen here)
- XDS-M: S = 5 independent random splits per donor (split seeds data seed + 100 + s); each split standardized as XDS (design-only,
  pooled within the half); per-donor q_d = mean over splits of R^A_dj R^B_dk; one-sample t test of (q_d) across donors, N(0,1) reference.
- KX (kappa-extrapolation): y_d = R*_dj R*_dk with full-donor design-only standardization; x_d = kappa_d / sqrt(v*_dj v*_dk);
  OLS of y on (1, x); test the intercept with the HC0 sandwich standard error, N(0,1) reference.
Evaluation: simulations L200/L500/L1000, LS200/LS500/LS1000 (same data seeds as XU_PREREG/COUPLE_PREREG, R = 50) and the six small
designs (R = 100): FDR, power (MCSE). SEA-AD real-cell null, new reassignments k = 10..14 (seeds 31,000,000 + k).

## R3. Differential co-expression with group-specific centring
XDS-DCg: as XDS-DC (DC_PREREG) but mu, sigma2 and tau2 are estimated within each group (AD, control), so group mean shifts in
expression do not enter the products. Re-run on the DC_PREREG null reassignments (k = 10, 11, 12) and on real donors.
Decision "XDS-DCg holds" as in DC_PREREG. Real-donor calls reported descriptively, with the overlap with the pooled-centre XDS-DC calls.

## R4. Positive control on real donors (descriptive)
XDS on real PsychAD donors (split seed 44,000,001): p-value and sign for all pairs between XIST and Y-linked genes present in the panel.
Expected: negative between-donor association (sex), detected by XDS.

## R5. Size of pseudobulk false edges (descriptive)
From the POP null at full D (K = 3): |r| of PB false calls (back-transformed from the Fisher-z p-value); median and counts above 0.1, 0.2.

## Amendment R6 (added before R6 was computed; R1/R2 runs were in progress, R3-R5 done)
Donor-centred pseudo-donor null, to bound the pooling inflation (pooling raised the median |within-donor r| by 37-46%).
For each real donor d and gene j: Paule-Mandel mu_j, tau2_j from (m_dj, psi_hat_dj); shrunk donor effect
theta~_dj = mu_j + tau2_j/(tau2_j + psi_hat_dj) (m_dj - mu_j); factor f_dj = mu_j / theta~_dj clipped to [0.1, 10];
every count of donor d's cells is multiplied by f_dj (library sizes unchanged). Then pool and reassign as in POP (K = 3,
seeds 49,000,000 + 10 k + subtype index), full D. Methods PB, XDS, XDS-M, KX. Report BH calls, per-pair type-I and top-1% type-I
(by within-pseudo-donor correlation), and the median |within-pseudo-donor r| relative to real donors. Descriptive; no decision.
