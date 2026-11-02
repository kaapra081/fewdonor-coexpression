# Plan 6 (2 Nov 2026, ~17:00 UTC). Written after internal review round 8 (review_round8.md, items M1, M2, M4, M5, M6, M10) and
# before any computation below. Code: rev5/rev5.py (smoke tests use seeds offset by +999 and write to a scratch directory; their
# outputs are used only to check that the code runs). Everything here is pre-specified relative to its own results; all of it was
# designed after the earlier results were known, and the paper says so.
# Data: PsychAD PFC myeloid subtypes P1-P3 (pop/pop_<s>_X.npz, as in every earlier plan), OneK1K B1-B3 (onek/), per-cell PsychAD
# library (poolID) from the same h5ad (rev5/pool_<s>.npz; rows in the same order as pop_<s>_meta.npz, checked donor by donor).
# Methods: PB (log-CPM Pearson, Fisher z), PBw (cell-number-weighted pseudobulk, as REV4 R3), KX, XDS-M (5 splits), all as in
# revm.py / rev4.py. BH at 5% over all panel pairs unless stated. "Donor-centred" = rev4_C.centre (shrunken donor-mean removal).

S. Planted between-donor covariance on real cells (M1; positive-control null).
 Per subtype: genes with pooled mean >= 0.1 UMI per cell, excluding XIST and the 7 Y genes, are eligible; 40 are drawn
 (seed 72,000,000 + si) into 4 modules of 10 with latent correlations rho = 0.1, 0.2, 0.3, 0.4 (180 true pairs).
 Between-donor log-SD s = median over eligible genes of tau_j / mu_j (Paule-Mandel on the real donors' ratio means), clipped to
 [0.1, 0.5]. Run k (k = 0..19) at size full and at size 200: seed 70,000,000 (full) or 71,000,000 (200) + 10 k + si; donor-centred
 cells are dealt to pseudo-donors with the real cell numbers (permutation of donor labels); for size 200, 200 pseudo-donors are
 drawn at random. Pseudo-donor d draws z_dk ~ N(0,1) (k = 1..4) and e_dj ~ N(0,1) (every gene); u_dj = s (sqrt(rho_k) z_dk +
 sqrt(1-rho_k) e_dj) for module genes, s e_dj otherwise; every cell of d has its gene-j counts multiplied by exp(u_dj - s^2/2).
 Truth: true pairs = within-module pairs; their true between-donor correlation is rho_lin = (exp(s^2 rho)-1)/(exp(s^2)-1).
 Per method and run: number of BH calls, FDP = false calls / max(calls, 1), power (overall and per module), per-pair type-I error
 at 0.05 over the null pairs. FDR = mean FDP over the 20 runs, MCSE = SD/sqrt(20).
 Effect size: XDS-r_jk = qbar_jk / sqrt(qbar_jj qbar_kk), with qbar_jk the mean over donors and splits of
 (R^A_dj R^B_dk + R^A_dk R^B_dj)/2 and qbar_jj the mean of R^A_dj R^B_dj (cross-half reliability); undefined if a qbar_jj <= 0.
 Compared with the pseudobulk correlation r_PB (log-CPM) for the 180 true pairs and 2,000 random null pairs per run:
 mean bias against rho_lin per module, RMSE, and mean |estimate| on null pairs.
 Rules. (S-a) XDS-M controls FDR with signal present if FDR <= 0.05 + 2 MCSE at both sizes in all 3 subtypes.
 (S-b) PB fails with signal present if its FDR at full size is >= 0.20 in all 3 subtypes. PBw and KX are judged by the bound of
 (S-a). (S-c) XDS-r is called approximately unbiased if |mean bias| <= 0.05 for every module at full size in all 3 subtypes.
 Power is descriptive.

D. Depth-matched null (M5; Lemma 1(ii), donor-specific depth distributions).
 Per subtype, k = 0..19, seed 73,000,000 + 10 k + si, full size. Cells are binned into 10 equal-count bins of log library size;
 within each bin the donor labels of the donor-centred cells are permuted. Each pseudo-donor therefore has its real donor's number
 of cells in every depth bin (its depth distribution) but cells from random donors. Methods: PB, KX, XDS-M and XDS-Md, which is
 XDS-M with every half's R residualized on the donor design covariates [1, log rhat_d, log n_d], rhat_d = sum s^2 / sum s over the
 donor's cells (a function of the design only). Reported: BH calls, per-pair type-I, runs with >= 1 call; SD of log rhat_d across
 pseudo-donors here, across real donors, and across the pseudo-donors of the plain donor-centred null (first reassignment).
 Rules. (D-a) XDS-M holds if runs with >= 1 false call <= 2 of 20 and mean per-pair type-I <= 0.055 in all 3 subtypes.
 (D-b) If D-a fails in any subtype, XDS-Md is evaluated by the same rule and is recommended only if it holds in all 3.

L. Library-preserving null (M4; PsychAD poolID; most donors are split across two libraries).
 Per subtype, k = 0..19, seed 74,000,000 + 10 k + si, full size. Cells are centred by exact donor-mean removal (factor
 clip(mu_j / m_dj, 0.1, 10), no shrinkage, so that residual donor biology is not carried into pseudo-donors); within each library
 the donor labels are permuted, so every pseudo-donor has its real donor's libraries and per-library cell numbers, and residual
 library effects are kept. Methods: PB, KX, XDS-M, and XDS-Mlib, which splits each donor's cells by library instead of at random
 (A = the library with the lexicographically smaller ID), uses only donors with >= 2 cells in each of exactly two libraries, and
 applies the XDS-M one-sample t test to the single split. Reported as for D, plus the number of donors used by XDS-Mlib.
 Rules. (L-a) XDS-M holds by the rule of D-a. (L-b) XDS-Mlib holds by the same rule. Descriptive otherwise.

C. Planted-confounder null (M6).
 Per subtype, k = 0..19, seed 75,000,000 + 10 k + si, full size; donor-centred cells dealt to pseudo-donors labelled with real donor
 IDs (as C1), which carry that donor's covariates. 200 eligible genes (as S; seed 76,000,000 + si) have their counts multiplied by
 1.2 in every pseudo-donor whose donor is an AD case. The 19,900 planted pairs are associated between donors marginally and null
 given AD status. Methods: PB, PBc (all C1 covariates), XDS-M, XDS-Mad (adjusted for [1, AD case, AD missing]), XDS-Mc (all C1
 covariates). Reported: calls among planted and among other pairs, per-pair type-I among planted pairs.
 Rules. (C-a) Adjustment works if XDS-Mc and XDS-Mad each have runs with >= 1 false call (any pair) <= 2 of 20 and per-pair type-I
 among planted pairs <= 0.055 in all 3 subtypes. (C-b) The confounder is detectable if XDS-M's mean calls among planted pairs
 exceed 1,000 in all 3 subtypes. Descriptive otherwise.

R. Real donors (M2, M4, M6, M10). Split seed 44,000,001 as C1real; all real donors.
 (R-1) OneK1K B1-B3 (ONEK-4): BH calls of PB and XDS-M. Descriptive.
 (R-2) P1-P3: XDS-r and r_PB for all pairs; pooled within-donor correlation r_w of cell ratios; Spearman correlation across pairs
 of (r_PB - XDS-r) with r_w, and the median of (r_PB - XDS-r) in the top and bottom deciles of r_w. Descriptive; the prediction
 from Proposition 1 is a positive correlation.
 (R-3) P1-P3: XDS-Md (as D) and XDS-Mlib (as L) on real donors: calls and the fraction of XDS-M calls they keep; for XDS-Mlib also
 against the single-split t test on random split 0 (XDS-t1), its equal-power comparator. Descriptive.
 (R-4) Sex control diagnostics, P1-P3, the 7 XIST-Y and 21 Y-Y pairs. (a) Expression sex per donor: male if log((Y+1)/(XIST+1)) > 0,
 with Y the summed UMI of the 7 Y genes and XIST the XIST UMI over the donor's cells; number of donors whose expression sex differs
 from the metadata. (b) XDS-Mc and sex-only XDS-M after excluding those donors. (c) Additionally removing discordant nuclei (in
 metadata-male donors: XIST >= 2 and Y = 0; in metadata-female donors: Y >= 2 and XIST = 0). (d) XDS-M within male donors and within
 female donors separately (BH over all pairs within the stratum). (e) XDS-Mlib for the sex pairs. (f) Among male donors, Spearman
 correlation of age with the donor's mean full-data design-standardized Y-gene score (loss of chromosome Y increases with age).
 For each: number of the 7 and 21 pairs called and their median p; sign of the mean cross-product for XIST-Y. Descriptive.
