# Plan 7 (2 Nov 2026, ~17:40 UTC): second-disease replication in the Perez et al. 2022 lupus PBMC cohort.
# Written before any expression value of this dataset was read and before any analysis below. What was read before this plan:
# Census obs (cell type, donor, disease, sex, ethnicity, raw_sum) and the h5ad obs columns library_uuid, sample_uuid,
# Processing_Cohort, disease_state (HTTP range reads of the obs group only; no expression), used only to count cells, donors and
# libraries so that the rules below can be applied. Like Plan 6 it is designed after the earlier results were known, and the paper says so.
# Code: lupus/lupus_prep.py (Census extraction), lupus/lupus_run.py (analyses); methods are the repo's own (rev/revm.py stats_from_cells and
# methods; dstd/dstd.py design_R; dev/methods_v2.py bh, pearson_fisher, paule_mandel; rev4/rev4.py pb_family -> PBw; rev4/rev4_C.py centre and
# cells_of; rev2/rev2_D.py PB-DC and XDS-Mgu; rev5/rev5.py base, xdsm_W, Rsplits, xds_r, lib_split_cells idea). Smoke tests use seeds +999,
# run one reassignment on the real panel and write to a scratch directory; their outputs are used only to check that the code runs.
# Methods and conventions are unchanged from Plan 6: PB (log1p(1e6 m) Pearson, Fisher z), PBw (cell-number-weighted pseudobulk), KX, XDS-M
# (5 random within-donor half splits, design-standardized, one-sample t across donors), BH at 5% over all panel pairs unless stated.
# "Donor-centred" = rev4_C.centre (shrunken donor-mean removal, factor clip(mu_j / theta_dj, 0.1, 10)); "exact centring" = the same
# with theta_dj = the donor's own ratio mean (no shrinkage), as Plan 6 L.

0. Data.
 Dataset: CELLxGENE collection 436154da-bcf1-4130-9c8b-120ff9a888f2 (Perez et al. 2022 Science, doi 10.1126/science.abf1970), dataset
 218acb0f-9f2f-4f76-b90b-15a4b7c7f629 ("multiplexed scRNA-seq of 1.2 million PBMCs from adult lupus samples"), dataset version
 4118e166-34f5-4c1f-9eed-c64b90a3dace; CELLxGENE Census release 2025-11-08 (the release used for OneK1K). 1,263,676 cells, all 10x 3' v2,
 all primary; 261 donors = 162 SLE ('systemic lupus erythematosus') + 99 'normal' (donor-level, constant within donor); 244 female, 17 male.
 Row order of the h5ad obs equals the Census soma_joinid order (checked: donor, cell type, disease, sex and observation_joinid agree on
 all 1,263,676 rows). Batch variables from the h5ad obs: library_uuid (88 libraries; 7-19 donors each; every donor is in 4-16 libraries),
 sample_uuid (274 samples; 11 donors have 2-3 samples), Processing_Cohort (4 cohorts; each library lies in exactly one; cohort 1 has no SLE donor).
 Library size lib = Census obs raw_sum (sum of raw UMI over all measured genes of the cell).
 Cell types. LM = Census cell_type 'classical monocyte' (CD14+ monocytes; the type with a strong interferon program).
 FC (few cells per donor) is chosen by the rule: among Census cell types with at least 200 donors having >= 4 cells, the type with the smallest
 median cells per donor (over donors with >= 4 cells). By the obs counts this is 'plasmacytoid dendritic cell' (235 donors, 5,181 cells, median 15
 per donor; the next smallest is 'lymphocyte', 258 donors, median 21). LM: 261 donors, 307,429 cells, median 1,005.
 Donor inclusion. A donor is used for a cell type if it has >= 4 cells of that type, counting all its samples. One unit per donor_id: the cells of
 all samples of a donor are pooled (11 donors have 2-3 samples; visit and sample are not modelled). SLE label = Census disease.
 Gene panel (the PsychAD rule of pop/pop_prep.py and real/real_prep.py): cells of the included donors; cps_g = mean over these cells of
 count_g / lib * 1e4 over all Census genes; top 400 genes by cps; then every forced gene absent from the top 400 with cps > 0 is added; panel
 sorted by Census gene index (so no data-driven ordering). Built separately for LM and FC. Adaptations that must change, and nothing else:
 (i) cps is computed exactly over all included cells by streaming the Census X 'raw' layer in cell blocks (no subsample); (ii) genes are Census
 features, symbol = feature_name, first match for the forced symbols; (iii) the forced list keeps the heat-shock and sex-linked genes of FORCE
 (HSPA1A HSPA1B HSP90AA1 HSP90AB1 HSPH1 HSPD1 HSPE1 DNAJB1 DNAJA1 HSPA8 HSPB1 XIST UTY USP9Y ZFY DDX3Y KDM5D NLGN4Y TTTY14 KDM5C KDM6A) and
 replaces the four brain perivascular-macrophage markers (CD163 MRC1 LYVE1 F13A1), which have no role in blood, by the interferon-stimulated
 gene list ISG = ISG15 IFI6 IFI44L IFI44 MX1 MX2 OAS1 OAS2 OAS3 OASL IFIT1 IFIT2 IFIT3 IFIT5 IFI27 RSAD2 IFITM1 IFITM3 STAT1 IRF7 XAF1 LY6E USP18
 EPSTI1 SIGLEC1 HERC5 CMPK2 SAMD9L PLSCR1 IFI35 BST2 TRIM22 (fixed here, before any expression was read; also used to label ISG-ISG pairs in C).
 Provenance reported: Census release, dataset ids, cell and donor counts per group, panel size, forced genes present, ISG genes in each panel.

A. Donor-centred null (positive: pseudobulk fails with many donors; this is the main replication of the paper's key result).
 Per cell type (si = 0 for LM, 1 for FC), k = 0..19, seed 80,000,000 + 10 k + si, full size. Donor-centred cells; the donor-label vector is
 permuted over all cells, so every pseudo-donor has a real donor's number of cells and Sigma_B = 0. Methods: PB, PBw, KX, XDS-M. Reported per
 method: BH calls (all false), per-pair type-I (fraction of pairs with p < 0.05), runs with >= 1 call, over the 20 runs; D, n_d summary.
 Descriptive size check (A2): the same with 100 pseudo-donors drawn at random (seed 83,000,000 + 10 k + si), PB and XDS-M; ratio of mean PB calls
 at full size to those at 100.
 Rules. A method 'holds' in a cell type if runs with >= 1 false call <= 2 of 20 and mean per-pair type-I <= 0.055; otherwise it 'fails'.
 (A-a) XDS-M holds in the cohort if it holds in both cell types. (A-b) PB fails in a cell type if it does not hold. (A-c) The lupus data
 replicate the key result if XDS-M holds in the cohort and PB fails in at least one cell type; they replicate it fully if PB fails in both.
 (A-d) Stated expectation from the theory (leak ~ 1/n_d): PB's mean per-pair type-I is larger in FC than in LM. If PB holds in LM it is
 reported as such. PBw and KX are judged by the same 'holds' rule; descriptive otherwise.

B. SLE-label null for false differential co-expression (replicates rev2_D D1 with SLE labels).
 Per cell type, k = 0..19, seed 81,000,000 + 10 k + si. Donor-centred cells are dealt to pseudo-donors that keep each real donor's cell number
 and SLE label (donor-label permutation; pseudo-donor u carries real donor u's label), so the true between-group difference in
 between-donor co-expression is zero and the group sizes and group cell numbers are the real ones. Methods: PB-DC (Fisher-z difference of
 log-CPM Pearson correlations, normal reference) and XDS-Mgu (unstandardized group-centred cross-split products averaged over 5 splits, Welch t),
 as in rev2_D. Reported: BH calls per run, per-pair type-I, runs with >= 1 call; for the real donors (C-2) the real calls; the ratio
 rho_PB = (mean over the 20 null runs of PB-DC calls) / (real-donor PB-DC calls) and rho_XDS likewise (undefined if the real count is 0).
 Rules. (B-a) XDS-Mgu holds in a cell type by the A rule (runs >= 1 call <= 2 of 20 and mean per-pair type-I <= 0.055) in both cell types.
 (B-b) PB-DC fails in a cell type if it does not hold. (B-c) The real PB-DC calls are called not interpretable as differential co-expression
 in a cell type if rho_PB >= 0.5.

C. Real donors (all donors of the cell type, split seed 44,000,001).
 (C-1) BH calls of PB, PBw, KX, XDS-M; fraction of XDS-M calls also called by PB and the reverse; the Proposition 1 check of Plan 6 R-2:
 XDS-r and r_PB for all pairs, pooled within-donor correlation r_w of cell ratios, Spearman correlation across pairs of (r_PB - XDS-r) with r_w,
 and the median of (r_PB - XDS-r) in the top and bottom deciles of r_w; the prediction from Proposition 1 is a positive Spearman correlation.
 (C-2) SLE vs control: BH calls of PB-DC and XDS-Mgu, sign (correlation in SLE larger than in control = positive), overlap of the two call
 sets, and, for pairs with both genes in ISG, the number of such pairs, the number called by each method and the fraction of those calls with
 positive sign. Descriptive. The cohort-by-disease table (cohort 1 has no SLE donor) is reported as a caveat on every real-donor DC result.

D. Library-preserving null (Plan 6 L with library_uuid as the pool).
 Per cell type, k = 0..19, seed 82,000,000 + 10 k + si, full size. Exact centring; within each library_uuid the donor labels of the cells are
 permuted, so every pseudo-donor has its real donor's libraries and per-library cell numbers and residual library effects are kept.
 Methods: PB, KX, XDS-M and XDS-Mlib. XDS-Mlib (adapted because every donor lies in 4-16 libraries, not two): for each of 5 splits the libraries of
 the pseudo-donor are assigned to half A or half B by independent fair coins (generator seeded by the run seed, split and donor index); a draw
 is accepted if both halves hold >= 2 cells (up to 100 draws); donors with no accepted draw in some split are dropped from that run (their
 number is reported); R of each half is design_R over the kept donors; the p-value is the XDS-M one-sample t test on the split-averaged
 cross-product. Descriptive extra (D3): PB-DC and XDS-Mgu SLE-label calls under this null (pseudo-donor u carries real donor u's label;
 real donors' libraries and cohorts stay attached to their label, so a cohort-label confounding is kept).
 Rules. (D-a) XDS-M holds in the cohort by the A rule. (D-b) XDS-Mlib holds in the cohort by the same rule. (D-c) If D-a fails, XDS-Mlib is
 recommended only if D-b holds. PB, KX descriptive (holds or fails by the same rule).
 Real donors (D-R): XDS-Mlib calls, number of donors used, fraction of XDS-M calls it keeps, and the single-split t test on a random split
 of the same donors (XDS-t1) as its equal-power comparator, with the fraction of XDS-t1 calls kept.

E. Reporting and limits. All pre-specified results are reported against these rules whatever the outcome. Analyses of Plan 6 that are not run
 here (planted covariance S, depth-matched null D, planted confounder C, sex diagnostics R-4) are not part of this plan. Any further analysis is
 labelled post hoc. Deviations in data scope or number of reassignments are declared with magnitude in LUPUS_RESULTS.md. Compute: Modal CPU
 image im-vNcuP00ryKr0nZlE9sgcfw, results under ./out/, each job resumable by skipping existing outputs.
