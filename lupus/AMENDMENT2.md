# Plan 7 Amendment 2 (2 Nov 2026, ~19:20 UTC). POST HOC: written after the Plan 7 results and the results of Amendment 1 were known.
# Plan 7 (sha256 2868352532bac24df94aae87f9583dcb3b21e2d9abbb2339b4d7f842e47954c2) and Amendment 1 (sha256 8cc98ba32e83f77540d9ff49e83b384ea851bf9016b64771d23096924e5610a7)
# are unchanged. What was seen before this amendment: all Plan 7 results; in Amendment 1, composition covariates for the 4 processing cohorts
# and the top-10 PCs of the 88-library composition did not restore control of XDS-M or XDS-Mlib in classical monocytes under the library-preserving
# null (best: XDS-Mlibcoh, 9,258 mean BH calls per run), while XDS-Mlibcoh held in the few-cell type. An inspection of the library-by-donor
# incidence of the dataset (all 1,263,676 cells; counts only) then showed that the 88 libraries consist of 24 groups of libraries with an
# identical donor set (20 groups of 4 libraries, one of 3, two of 2, one of 1): the lanes of multiplexed pools of 7-19 donors. 180 of 261 donors
# lie in one pool, 62 in two, 10 in three, 9 in four; each pool lies in one processing cohort. Pool-level effects on gene means are therefore shared by
# all lanes of a pool and by every half of a lane-split, and the 4 cohort fractions of Amendment 1 are too coarse to absorb them. This amendment
# asks whether donor-level pool composition removes the failure. All of it is post hoc and none of it re-judges Plan 7.
# Code: lupus/lupus_amend2.py; pool map: lupus/lupus_pools.json (library_uuid -> pool, derived as stated above from the h5ad obs).

Covariates: W_pool = [1, fraction of the pseudo-donor's cells in pools 2..24] (24 columns; pool numbering by smallest library_uuid; pool 1 is the reference).
Variants:
 XDS-Mpool     = XDS-M (random within-donor cell splits) with every half's R residualized on W_pool.
 XDS-Mlibpool  = XDS-Mlib (Plan 7 D procedure; donors without an accepted library split dropped) with W_pool (rows of the kept donors).
 PBpool        = pseudobulk log-CPM residualized on W_pool before Pearson, Fisher z with n - 3 - 23 (rev4_C.pb_adj).
 XDS-Mgupool   = XDS-Mgu (Plan 7 B) after residualizing the half means mA, mB (all donors jointly) on W_pool, then group-centring and the Welch test.
Runs, per cell type (si = 0 LM, 1 FC), 20 reassignments each:
 F1 (library-preserving null): seed 85,000,000 + 10 k + si, construction as Plan 7 D (exact centring; donor labels permuted within library_uuid).
    Co-expression methods: PB, PBpool, XDS-M, XDS-Mlib, XDS-Mpool, XDS-Mlibpool. Differential co-expression with the real SLE labels carried by the pseudo-donors:
    PB-DC, XDS-Mgu, XDS-Mgupool.
 F2 (validity, co-expression): donor-centred null of Plan 7 A, seeds 80,000,000 + 10 k + si; methods XDS-M, XDS-Mpool, PBpool.
 F2B (validity, differential co-expression): SLE-label null of Plan 7 B, seeds 81,000,000 + 10 k + si; methods XDS-Mgu, XDS-Mgupool.
 FR (real donors, split seed 44,000,001): BH calls of PB, PBpool, XDS-M, XDS-Mpool, XDS-Mlibpool and the fraction of XDS-M calls each keeps; SLE vs control calls of PB-DC,
    XDS-Mgu, XDS-Mgupool, their overlaps, and for the 496 ISG-ISG pairs the number called with positive sign (SLE larger).
Rules ('holds' as in Plan 7: runs with >= 1 call <= 2 of 20 and mean per-pair type-I <= 0.055).
 (F-a) XDS-Mpool and XDS-Mlibpool each restore control in the cohort if they hold in F1 in both cell types; (F-b) a variant is valid only if it holds in F2
 in both cell types; (F-c) XDS-Mgupool restores control for differential co-expression if it holds in the DC part of F1 in both cell types and is valid only if
 it holds in F2B in both. (F-d) FR is descriptive. A variant that fails F-b is reported as invalid whatever it does in F1.
