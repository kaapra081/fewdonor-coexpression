# Plan 6b (2 Nov 2026, ~17:25 UTC). POST HOC. Written after the Plan 6 results were seen (job 48a2a282, all 66 tasks exit 0):
# D-a failed (XDS-M false calls in 20/20 depth-matched runs in P1 and P2) and D-b held (XDS-Md); L-a and L-b failed (XDS-M and
# XDS-Mlib false calls in 19-20/20 library-preserving runs); C-a failed in P1 and P2 (XDS-Mc / XDS-Mad: 20/20 runs with false calls,
# per-pair type-I among planted pairs 0.26-0.37); S-a, S-b, S-c held. Everything below is post hoc and will be labelled so.
# Code: rev5/rev5_post.py (imports rev5.py unchanged for data loading and helpers).

H1. Adjust first, standardize second (XDS-Mca). Diagnosis of C-a: XDS-Mc residualizes the design-standardized R on covariates,
    but a covariate effect that is multiplicative on expression becomes, after standardization, delta mu_j / sqrt(v*_dj), which
    varies with the donor's design; its residual is shared by genes and gives a positive cross-product. XDS-Mca instead
    residualizes each half's ratio means m^h (OLS across donors, per gene) on W and then divides by sqrt(v*^h_dj) from design_R.
    Its mean is zero under the conditional null whatever the design, because E[m^h_dj | W, design] is linear in W for indicator
    covariates and v* is a function of the design. Evaluated on: (a) the planted-confounder null of Plan 6 C (same seeds; AD-only
    and all-covariate W), rule as C-a; (b) the depth-matched null of Plan 6 D with W = [1, log rhat_d, log n_d], rule as D-a;
    (c) real donors (split seed 44,000,001): total calls and the fraction of XDS-M calls kept, all-covariate W; and the 7 XIST-Y and
    21 Y-Y pairs under all-covariate and sex-only W.
H2. Library-preserving null diagnosis. Under exact donor centring each real donor's cells sum to a fixed total, so two disjoint
    subsets of one donor's cells have means with covariance -Sigma^W/N (sampling without replacement). Within-library permutation
    puts many cells of the same real donor into both halves of a pseudo-donor; this predicts false calls whose sign is opposite to
    the within-donor correlation r_w. Diagnostics on the Plan 6 L null (same seeds): (a) the sign of XDS-M BH calls relative to
    r_w (fraction with sign(t) = -sign(r_w)); (b) source-split XDS-M: each pseudo-donor's cells are split by SOURCE donor (all
    cells of one real donor go to the same half; source donors are dealt at random, greedily balancing cell numbers; 5 splits;
    pseudo-donors with < 2 cells in a half are dropped); library effects remain shared by both halves, the finite-population term
    does not. If source-split XDS-M holds (rule of D-a) while XDS-M fails, the L failure is attributed to the null's construction
    rather than to library effects; otherwise library effects are not excluded.
H3. Bug fix of Plan 6 R-4(c): donors left with < 4 cells after removing discordant nuclei keep all their cells (the pre-specified
    code kept donors with 2-3 cells, whose halves of 1 cell gave undefined variances and NaN p-values for every pair in P1).
