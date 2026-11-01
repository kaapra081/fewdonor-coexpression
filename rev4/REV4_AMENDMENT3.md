# REV4 Amendment 3 (2 Nov 2026, ~00:00 UTC). Written after an internal round-5 review of the rev4 draft and after the R1 (two subtypes), R2, R2c,
# R3-R6 results had been seen; before any computation below.

C1. Covariate-adjusted XDS-M (XDS-Mc). Donor covariates W = intercept, male sex, age in years (development_stage; "90 year-old and over" = 90),
    brain bank (H, R versus M), AD status (case, missing versus control), ethnicity (African American, Hispanic, other or unknown versus European
    American). In every split the design-standardized half-donor values R^A, R^B are residualized on W by least squares before the per-donor products
    are formed and averaged over the 5 splits; one-sample t with a normal reference. PBc: log-CPM residualized on W, Pearson, Fisher z with D - 3 - p.
    (a) Null: donor-centred cells, full random reassignment (as R3), each pseudo-donor carrying the covariates of the real donor whose label it gets;
    20 reassignments per subtype, seeds 66,000,000 + 10 k + subtype index; methods XDS-M, XDS-Mc, PB, PBc. Rule: XDS-Mc holds if at most 2 of 20
    runs make a BH call, in each subtype. (b) Real donors (split seed 44,000,001): BH calls of XDS-M and XDS-Mc and their overlap; descriptive.
C2. Bank-mean-removed bank-stratified null. After donor-centring, every cell of brain bank b is multiplied by mu_j / mu_bj (ratio of the pooled ratio
    mean of gene j over all donor-centred cells to that over bank b's cells, clipped to [0.1, 10]); then the R2 reassignments (seeds 64,000,000 +
    10 k + si, 40 per subtype). Methods PB, XDS-M, XDS-Mb. Interpretation stated now: if XDS-M is within the n = 40 bound (<= 4/40) in every subtype,
    the R2 failures are attributed to residual bank means; otherwise donor-to-donor depth-distribution differences (Lemma 1, case ii) are implicated.
    Also reported: per bank the mean of log depth and E[s^2]/E[s] over cells.
C3. Noise baselines for the real-donor overlaps (R4): for P1-P3, overlap of the 5,000 strongest pairs and of the 20 highest-degree genes (degree
    within those 5,000 pairs) between (i) XDS-M with split seeds 44,000,001 and 44,000,002, (ii) XDS-M on two random donor halves (seed 42,000,000,
    as REV2 C), (iii) pseudobulk on the same halves. For P1, the odds ratio that a real pseudobulk-only call is a false call in the R6 null.
C4. Permutation XDS at 200 donors: P1, REV2 A reassignments k = 0..19 (seeds 50,000,000 + 10 k), nested D = 200 subset exactly as REV2 A;
    XDS and LCPd with 2,000,000 donor permutations, and the normal-moment XDS on the same runs. Rule: permutation XDS holds if <= 2/20 runs make a call.
C5. Depth-state coupling statistic of diag_real.py (pooled within-donor Pearson correlation of log library size with Y/lib) for PsychAD P1-P3 and
    OneK1K B1-B3 (all donors with >= 2 cells); median over genes and fraction with |c| > 0.045. Replaces the "Lee et al. fresh-tissue" value, which
    came from an 80-donor PsychAD microglia subset, not a separate cohort.
C6. Mantis gene map, robustness: (a) eta-squared of null false-edge degree (and of PB-only and XDS-M degree) across the Mantis leaf clusters versus
    10,000 random partitions of the same sizes (gene-label permutation), P1. (b) Per-gene null false-edge degree in the first R3 reassignment of P2 and
    P3 for genes in the P1 map; Spearman correlation with P1 and the rank of the heat-shock leaf (heat-shock genes were in the panel a priori as
    within-cell-program controls). Descriptive.
C7. REV3 A3 coincidence (XDS-M 49/100 runs at 25 donors in every subtype): run-index overlap across subtypes (computed from the REV3 per-run files).
