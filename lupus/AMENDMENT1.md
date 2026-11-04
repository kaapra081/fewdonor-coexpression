# Plan 7 Amendment 1 (2 Nov 2026, ~18:30 UTC). POST HOC: written after the Plan 7 results of sections A-D were known.
# Plan 7 (sha256 2868352532bac24df94aae87f9583dcb3b21e2d9abbb2339b4d7f842e47954c2) is unchanged. What was seen before this amendment: all
# Plan 7 results. In the library-preserving null D, XDS-M and the pre-specified library-split variant XDS-Mlib both fail in classical
# monocytes (about 48,700 and 34,700 BH calls per run), pseudobulk fails more, and in the few-cell type XDS-Mlib is close to the bound
# (0.1 calls per run, per-pair type-I 0.057). A library-preserving permutation keeps every donor's library composition, so any library or
# processing-cohort effect on gene means is a real between-pseudo-donor covariance, not a sampling artifact. This amendment asks, as a
# diagnostic and not as a replacement for Plan 7, whether donor-level covariates for the batch composition remove it. Nothing in Plan 7
# is re-judged by this amendment. Code: lupus/lupus_amend1.py; functions as in lupus_run.py.

Covariates (per pseudo-donor u, from the cells that carry its label in the run):
 W_coh = [1, fraction of u's cells in Processing_Cohort 2, 3, 4] (cohort 1 is the reference; 4 columns).
 W_pc  = [1, the top 10 principal components of the centred D x 88 matrix of u's library-composition fractions] (11 columns; the PCA is
         recomputed in every run from the pseudo-donors of that run).
Variants (every half's R of XDS-M is residualized on the covariates before the cross-product; p-values from the XDS-M one-sample t test):
 XDS-Mcoh   = XDS-M with random within-donor cell splits, W_coh.
 XDS-Mpc    = XDS-M with random splits, W_pc.
 XDS-Mlibcoh = XDS-Mlib (Plan 7 D procedure: libraries of a donor assigned to halves by coin flips, generator seeded by the run seed of the run it is in) with W_coh.
 XDS-Mlibpc = XDS-Mlib with W_pc.
Runs, per cell type (si = 0 LM, 1 FC):
 E1 (library-preserving null): k = 0..19, seed 84,000,000 + 10 k + si; construction exactly as Plan 7 D (exact centring; donor labels permuted
    within library_uuid). Methods: PB, XDS-M, XDS-Mlib (for reference, same run) and the four variants.
 E2 (validity check): the donor-centred null of Plan 7 A with its own seeds (80,000,000 + 10 k + si); methods: XDS-M, XDS-Mcoh, XDS-Mpc.
 ER (real donors, split seed 44,000,001): calls of XDS-M and the four variants, and the fraction of XDS-M calls that each variant keeps.
Rules (the 'holds' rule is that of Plan 7: runs with >= 1 call <= 2 of 20 and mean per-pair type-I <= 0.055).
 (E-a) A variant restores control in the cohort if it holds in E1 in both cell types. (E-b) A variant is valid only if it holds in E2 in both
 cell types; a variant that fails E2 is reported as invalid whatever it does in E1. (E-c) ER is descriptive.
All results of this amendment are post hoc and are labelled so in the report; they are not evidence for the pre-specified rules of Plan 7.
