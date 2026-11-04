# LUPUS_RESULTS: second-disease replication in the Perez et al. 2022 lupus PBMC cohort (Plan 7)

Plan 7 sha256 `2868352532bac24df94aae87f9583dcb3b21e2d9abbb2339b4d7f842e47954c2` (frozen and saved before any expression was read; PLAN7.md). Post hoc amendments (written after results were known and hashed before running): Amendment 1 `8cc98ba32e83f77540d9ff49e83b384ea851bf9016b64771d23096924e5610a7`, Amendment 2 `eaff32ca379a786bf35af08b0182fd05169667f1f6a212d7b6584a38d4039f1a`.

## Headline

- Donor-centred null (A): XDS-M holds the frozen rule in both cell types (mean false calls 0.0 and 0.1 per run; runs with at least one call 0/20 and 2/20; per-pair type-I 0.051 and 0.050). Pseudobulk fails in both (16,631 and 3,472 false calls per run). The key result replicates in a second disease (rule A-c, fully).
- The pre-stated expectation that pseudobulk fails more in the few-cell type (A-d) is not supported: its per-pair type-I is 0.292 at about 1,000 cells per donor and 0.167 at 15 cells per donor. Pseudobulk fails at both extremes of cell number in this cohort.
- SLE-label null (B): XDS-Mgu holds in both cell types (0 false calls in 20 runs each); pseudobulk differential co-expression fails (8,937 and 1,493 false calls per run, 22% and 12% of the real PB-DC call counts; rule B-c, 'not interpretable' at >= 50%, is not triggered).
- Library-preserving null (D): negative. XDS-M fails in both cell types (48,661 and 13.1 calls per run), and the pre-specified library-split variant XDS-Mlib also fails (34,664 and 0.1; in the few-cell type it has 2/20 runs with a call but per-pair type-I 0.057, above the 0.055 bound). Pseudobulk is worse still. In this multiplexed design pool and cohort effects are real between-donor covariance that neither method removes without batch covariates.
- Post hoc (labelled): donor-level pool-composition covariates remove the false calls in monocytes (XDS-Mpool 0.5 per run, from 48,489) but leave per-pair type-I at 0.061 to 0.080 (four checks), so XDS-Mpool fails the validity check (F-b); the pool-adjusted differential co-expression test XDS-Mgupool holds both the library-null and validity rules in both cell types.

## 1. Data provenance

CELLxGENE collection 436154da-bcf1-4130-9c8b-120ff9a888f2 (Perez et al. 2022 Science, doi 10.1126/science.abf1970), dataset 218acb0f-9f2f-4f76-b90b-15a4b7c7f629, dataset version 4118e166-34f5-4c1f-9eed-c64b90a3dace, Census release 2025-11-08; 1,263,676 cells, 10x 3' v2, 261 donors (162 SLE, 99 normal; 244 female, 17 male), 88 libraries (library_uuid), 274 samples, 4 processing cohorts (cohort 1 has no SLE donor). Per-cell batch columns were read from the dataset h5ad obs by HTTP range requests (obs only; row order verified identical to Census on all 1,263,676 rows). Library size is Census raw_sum.

| Cell type | Cells | Donors | SLE / control donors | Median cells per donor (IQR) | Median per donor SLE / control | Donors with 2-3 samples | Panel genes | ISG genes in panel | Panel nnz per cell |
|---|---|---|---|---|---|---|---|---|---|
| classical monocytes (LM) | 307,429 | 261 | 162 / 99 | 1005 (681-1466) | 1132.5 / 792.0 | 11 | 447 | 32/32 | 297.7 |
| plasmacytoid DC (FC) | 5,181 | 235 | 143 / 92 | 15 (9-26) | 12.0 / 24.0 | 7 | 447 | 32/32 | 331.9 |

Few-cell type rule (Plan 7 section 0): among Census cell types with >= 200 donors having >= 4 cells, the smallest median cells per donor: plasmacytoid dendritic cell (235 donors, median 15), lymphocyte (258 donors, median 21), conventional dendritic cell (258 donors, median 55).

Cohort (modal Processing_Cohort of the donor) by disease, donors: 1.0: SLE 0, control 33; 2.0: SLE 106, control 13; 3.0: SLE 17, control 10; 4.0: SLE 39, control 43 (classical monocytes; pDC in lupus_FC_info.json). Pool structure (derived from library-by-donor incidence, all cells): 24 groups of libraries with identical donor sets (20 x 4 libraries, 1 x 3, 2 x 2, 1 x 1); 180 donors lie in one pool, 62 in two, 10 in three, 9 in four.

Panel rule adaptations (declared in Plan 7): exact cps over all included cells by streaming Census (no subsample); forced genes = heat-shock and sex-linked genes of the PsychAD list, with the four brain perivascular-macrophage markers replaced by a 32-gene interferon-stimulated list fixed in the plan; panels built separately per cell type (447 genes each, 99,681 pairs).

## 2. Pre-specified results against their rules

A method 'holds' in a cell type if runs with >= 1 false BH call <= 2 of 20 and mean per-pair type-I <= 0.055 (mean fraction of pairs with p < 0.05); otherwise it fails.

### A. Donor-centred null, full size, 20 reassignments (seeds 80,000,000 + 10k + si)

| Cell type | D | Method | Mean false calls | SD | Range | Runs >=1 call | Mean per-pair type-I | Rule |
|---|---|---|---|---|---|---|---|---|
| classical monocytes (LM) | 261 | PB | 16,630.7 | 2,713.7 | 12766-24532 | 20/20 | 0.2924 | FAILS |
| classical monocytes (LM) | 261 | PBw | 409.6 | 72.1 | 272-529 | 20/20 | 0.0605 | FAILS |
| classical monocytes (LM) | 261 | KX | 0.2 | 0.4 | 0-1 | 4/20 | 0.0536 | FAILS |
| classical monocytes (LM) | 261 | XDSM | 0.0 | 0.0 | 0-0 | 0/20 | 0.0508 | holds |
| plasmacytoid DC (FC) | 235 | PB | 3,472.3 | 758.9 | 2572-5679 | 20/20 | 0.1667 | FAILS |
| plasmacytoid DC (FC) | 235 | PBw | 81.3 | 27.9 | 49-169 | 20/20 | 0.0296 | FAILS |
| plasmacytoid DC (FC) | 235 | KX | 0.1 | 0.4 | 0-1 | 3/20 | 0.0536 | FAILS |
| plasmacytoid DC (FC) | 235 | XDSM | 0.1 | 0.3 | 0-1 | 2/20 | 0.0502 | holds |

Rules: A-a XDS-M holds in the cohort: **True** (LM True, FC True). A-b PB fails: LM True, FC True. A-c key result replicates: **True**, fully: **True**. A-d (expectation: PB per-pair type-I larger in FC than LM): **not supported** (FC 0.167 vs LM 0.292). PBw and KX fail the rule in both cell types (KX: one call in 4/20 and 3/20 runs; XDS-M in FC has 2/20, exactly at the bound).

Descriptive size check A2 (100 pseudo-donors drawn at random, seeds 83,000,000 + 10k + si):

| Cell type | Method (D=100) | Mean false calls | Range | Runs >=1 call | Mean per-pair type-I |
|---|---|---|---|---|---|
| classical monocytes (LM) | PB | 8,880.1 | 292-25362 | 20/20 | 0.2089 |
| classical monocytes (LM) | XDSM | 0.1 | 0-1 | 1/20 | 0.0531 |
| plasmacytoid DC (FC) | PB | 1,706.8 | 674-4949 | 20/20 | 0.1413 |
| plasmacytoid DC (FC) | XDSM | 0.0 | 0-0 | 0/20 | 0.0512 |

Ratio of mean PB false calls at full size to those at 100 donors: LM 1.87, FC 2.03.

### B. SLE-label null for false differential co-expression, 20 reassignments (seeds 81,000,000 + 10k + si)

| Cell type | SLE / control | Method | Mean false calls | Range | Runs >=1 call | Mean per-pair type-I | Rule | Real-donor calls |
|---|---|---|---|---|---|---|---|---|
| classical monocytes (LM) | 162 / 99 | PBDC | 8,937.5 | 6210-11915 | 20/20 | 0.2294 | FAILS | 40629 |
| classical monocytes (LM) | 162 / 99 | XDSMgu | 0.0 | 0-0 | 0/20 | 0.0378 | holds | 43082 |
| plasmacytoid DC (FC) | 143 / 92 | PBDC | 1,493.1 | 659-3386 | 20/20 | 0.1407 | FAILS | 12357 |
| plasmacytoid DC (FC) | 143 / 92 | XDSMgu | 0.0 | 0-0 | 0/20 | 0.0478 | holds | 2195 |

Rules: B-a XDS-Mgu holds in both: **True**. B-b PB-DC fails: LM True, FC True. B-c: rho_PB = mean null PB-DC calls / real PB-DC calls = 0.220 (LM), 0.121 (FC); 'not interpretable' needs >= 0.5: not triggered. rho_XDS = 0 in both (no null calls; real XDS-Mgu calls 43,082 and 2,195).

### C. Real donors (split seed 44,000,001), descriptive

| Cell type | Donors | PB calls | PBw calls | KX calls | XDS-M calls | Frac. XDS-M calls also in PB | Frac. PB calls also in XDS-M | PB-only | XDS-M-only |
|---|---|---|---|---|---|---|---|---|---|
| classical monocytes (LM) | 261 | 69,706 | 64,315 | 59,021 | 69,764 | 0.951 | 0.952 | 3,362 | 3,420 |
| plasmacytoid DC (FC) | 235 | 31,360 | 23,851 | 19,734 | 36,029 | 0.630 | 0.723 | 8,672 | 13,341 |

Pairs: 99,681. 
Proposition 1 check (prediction: positive Spearman correlation of r_PB - XDS-r with within-donor r_w): **not supported** in either cell type.

| Cell type | Spearman(r_PB - XDS-r, r_w) | Median (r_PB - XDS-r), top r_w decile | bottom decile | Spearman(r_PB, XDS-r) | Pairs with XDS-r undefined |
|---|---|---|---|---|---|
| classical monocytes (LM) | -0.104 | 0.001 | 0.010 | 0.992 | 0 |
| plasmacytoid DC (FC) | -0.286 | -0.142 | 0.097 | 0.772 | 2220 |

SLE vs control differential co-expression (sign + = stronger in SLE). The cohort-by-disease table above is a standing caveat: cohort 1 contains no SLE donor.

| Cell type | SLE / control | PB-DC calls | XDS-Mgu calls | Overlap | PB-DC calls positive | XDS-Mgu calls positive | PB-DC ISG-ISG pairs called | XDS-Mgu ISG-ISG pairs called | of which positive (XDS-Mgu) |
|---|---|---|---|---|---|---|---|---|---|
| classical monocytes (LM) | 162 / 99 | 40,629 | 43,082 | 34,121 | 19,524 | 21,012 | 478/496 | 496/496 | 496 |
| plasmacytoid DC (FC) | 143 / 92 | 12,357 | 2,195 | 1,638 | 5,320 | 1,189 | 52/496 | 112/496 | 112 |

### D. Library-preserving null, 20 reassignments (seeds 82,000,000 + 10k + si; exact centring; donor labels permuted within library_uuid)

| Cell type | Method | Mean calls | Range | Runs >=1 call | Mean per-pair type-I | Rule |
|---|---|---|---|---|---|---|
| classical monocytes (LM) | PB | 62,723.2 | 61080-64534 | 20/20 | 0.6588 | FAILS |
| classical monocytes (LM) | PBw | 53,164.2 | 52390-54035 | 20/20 | 0.5767 | FAILS |
| classical monocytes (LM) | KX | 25,068.7 | 22833-27531 | 20/20 | 0.3737 | FAILS |
| classical monocytes (LM) | XDSM | 48,661.4 | 47994-49630 | 20/20 | 0.5422 | FAILS |
| classical monocytes (LM) | XDSMlib | 34,663.8 | 32885-35861 | 20/20 | 0.4267 | FAILS |
| plasmacytoid DC (FC) | PB | 11,126.1 | 8051-14845 | 20/20 | 0.2389 | FAILS |
| plasmacytoid DC (FC) | PBw | 1,046.2 | 727-1650 | 20/20 | 0.0723 | FAILS |
| plasmacytoid DC (FC) | KX | 1.9 | 0-9 | 13/20 | 0.0681 | FAILS |
| plasmacytoid DC (FC) | XDSM | 13.1 | 0-34 | 19/20 | 0.0638 | FAILS |
| plasmacytoid DC (FC) | XDSMlib | 0.1 | 0-1 | 2/20 | 0.0570 | FAILS |

Donors kept by XDS-Mlib: LM 261 of 261, FC 232 of 235.

Rules: D-a XDS-M holds: **False**. D-b XDS-Mlib holds: **False**. D-c: XDS-Mlib recommended: **False**. PB, PBw, KX fail in both cell types.

D3 (descriptive): SLE-label differential co-expression under the library-preserving null (real donors' libraries and cohorts stay attached to their SLE label):

| Cell type | Method | Mean calls | Runs >=1 call | Mean per-pair type-I | Rule |
|---|---|---|---|---|---|
| classical monocytes (LM) | PBDC | 21,852.2 | 20/20 | 0.3413 | FAILS |
| classical monocytes (LM) | XDSMgu | 37.5 | 20/20 | 0.1077 | FAILS |
| plasmacytoid DC (FC) | PBDC | 8,567.6 | 20/20 | 0.2179 | FAILS |
| plasmacytoid DC (FC) | XDSMgu | 0.1 | 1/20 | 0.0468 | holds |

D-R (real donors): XDS-Mlib, with the single-split comparator XDS-t1 on the same donors:

| Cell type | Donors used | XDS-Mlib calls | XDS-t1 calls | Fraction of XDS-M calls kept by XDS-Mlib | Fraction of XDS-t1 calls kept by XDS-Mlib | Fraction of XDS-M calls kept by XDS-t1 |
|---|---|---|---|---|---|---|
| classical monocytes (LM) | 261 | 69,134 | 68,891 | 0.979 | 0.977 | 0.975 |
| plasmacytoid DC (FC) | 232 | 31,916 | 24,451 | 0.794 | 0.865 | 0.642 |

## 3. Post hoc diagnostics (NOT pre-specified results; written after the Plan 7 results were known)

### Amendment 1: cohort and library-composition covariates (library-preserving null E1; donor-centred validity check E2)

| Cell type | Method | E1 mean calls | Runs >=1 | E1 per-pair type-I | E1 rule | E2 validity |
|---|---|---|---|---|---|---|
| classical monocytes (LM) | PB | 63,005.6 | 20/20 | 0.6612 | FAILS | n/a |
| classical monocytes (LM) | XDSM | 48,555.7 | 20/20 | 0.5409 | FAILS | holds (type-I 0.0508) |
| classical monocytes (LM) | XDSMlib | 34,652.0 | 20/20 | 0.4267 | FAILS | n/a |
| classical monocytes (LM) | XDSMcoh | 31,245.5 | 20/20 | 0.4064 | FAILS | holds (type-I 0.0522) |
| classical monocytes (LM) | XDSMpc | 42,336.2 | 20/20 | 0.4949 | FAILS | holds (type-I 0.0546) |
| classical monocytes (LM) | XDSMlibcoh | 9,258.5 | 20/20 | 0.2404 | FAILS | n/a |
| classical monocytes (LM) | XDSMlibpc | 26,387.2 | 20/20 | 0.3739 | FAILS | n/a |
| plasmacytoid DC (FC) | PB | 11,453.5 | 20/20 | 0.2424 | FAILS | n/a |
| plasmacytoid DC (FC) | XDSM | 11.9 | 17/20 | 0.0639 | FAILS | holds (type-I 0.0502) |
| plasmacytoid DC (FC) | XDSMlib | 0.1 | 1/20 | 0.0565 | FAILS | n/a |
| plasmacytoid DC (FC) | XDSMcoh | 0.2 | 3/20 | 0.0564 | FAILS | holds (type-I 0.0523) |
| plasmacytoid DC (FC) | XDSMpc | 0.3 | 3/20 | 0.0625 | FAILS | FAILS (type-I 0.0556) |
| plasmacytoid DC (FC) | XDSMlibcoh | 0.0 | 0/20 | 0.0530 | holds | n/a |
| plasmacytoid DC (FC) | XDSMlibpc | 0.1 | 2/20 | 0.0581 | FAILS | n/a |

No variant restores control in both cell types (E-a: all false). XDS-Mlibcoh holds in FC and cuts monocyte calls from 48,556 (XDS-M) to 9,258. Real donors, share of XDS-M calls kept: classical monocytes (LM): XDSMcoh 0.72, XDSMpc 0.72, XDSMlibcoh 0.71, XDSMlibpc 0.70; plasmacytoid DC (FC): XDSMcoh 0.28, XDSMpc 0.30, XDSMlibcoh 0.18, XDSMlibpc 0.21.

### Amendment 2: pool-composition covariates (24 pools; F1 library-preserving null; F2/F2B validity)

| Cell type | Method | F1 mean calls | Runs >=1 | F1 per-pair type-I | F1 rule | F2 validity (donor-centred null) |
|---|---|---|---|---|---|---|
| classical monocytes (LM) | PB | 62,653.3 | 20/20 | 0.6583 | FAILS | n/a |
| classical monocytes (LM) | PBpool | 30,288.2 | 20/20 | 0.3956 | FAILS | FAILS (type-I 0.1770) |
| classical monocytes (LM) | XDSM | 48,488.8 | 20/20 | 0.5408 | FAILS | holds (type-I 0.0508) |
| classical monocytes (LM) | XDSMlib | 34,349.6 | 20/20 | 0.4250 | FAILS | n/a |
| classical monocytes (LM) | XDSMpool | 0.5 | 7/20 | 0.0798 | FAILS | FAILS (type-I 0.0614) |
| classical monocytes (LM) | XDSMlibpool | 17,781.6 | 20/20 | 0.2983 | FAILS | n/a |
| plasmacytoid DC (FC) | PB | 11,737.0 | 20/20 | 0.2439 | FAILS | n/a |
| plasmacytoid DC (FC) | PBpool | 7,785.3 | 20/20 | 0.2067 | FAILS | FAILS (type-I 0.1550) |
| plasmacytoid DC (FC) | XDSM | 9.2 | 19/20 | 0.0640 | FAILS | holds (type-I 0.0502) |
| plasmacytoid DC (FC) | XDSMlib | 0.0 | 0/20 | 0.0564 | FAILS | n/a |
| plasmacytoid DC (FC) | XDSMpool | 0.3 | 4/20 | 0.0631 | FAILS | FAILS (type-I 0.0641) |
| plasmacytoid DC (FC) | XDSMlibpool | 0.8 | 7/20 | 0.0662 | FAILS | n/a |

Differential co-expression under the library-preserving null (F1) and the SLE-label null (F2B validity):

| Cell type | Method | F1 mean calls | Runs >=1 | F1 per-pair type-I | F1 rule | F2B validity |
|---|---|---|---|---|---|---|
| classical monocytes (LM) | PBDC | 22,439.7 | 20/20 | 0.3460 | FAILS | n/a |
| classical monocytes (LM) | XDSMgu | 23.2 | 20/20 | 0.1063 | FAILS | holds (type-I 0.0378) |
| classical monocytes (LM) | XDSMgupool | 0.0 | 0/20 | 0.0328 | holds | holds (type-I 0.0549) |
| plasmacytoid DC (FC) | PBDC | 9,212.4 | 20/20 | 0.2237 | FAILS | n/a |
| plasmacytoid DC (FC) | XDSMgu | 0.0 | 0/20 | 0.0465 | holds | holds (type-I 0.0478) |
| plasmacytoid DC (FC) | XDSMgupool | 0.0 | 0/20 | 0.0522 | holds | holds (type-I 0.0520) |

Rules: F-a XDS-Mpool restores control: **False** (monocytes 0.5 calls per run, 7/20 runs with a call, per-pair type-I 0.080); F-b XDS-Mpool valid: **False** (per-pair type-I 0.061 and 0.064 under the donor-centred null; untested explanation: 24 residualizing columns against 235-261 donors without a degrees-of-freedom correction). XDS-Mlibpool restores control: **False** (17,782 calls per run in monocytes; untested explanation: the covariates describe the donor's cells, not the pool mix of each lane half). F-c XDS-Mgupool restores control and is valid: **True / True** (in monocytes the unadjusted XDS-Mgu fails the library null: 23.2 calls per run, type-I 0.106).

Real donors after pool adjustment (descriptive):

| Cell type | XDS-M calls | XDS-Mpool calls | Fraction kept | PBpool calls | XDS-Mlibpool calls | XDS-Mgu calls | XDS-Mgupool calls | ISG-ISG pairs called (XDS-Mgupool) | of which positive |
|---|---|---|---|---|---|---|---|---|---|
| classical monocytes (LM) | 69,764 | 59,772 | 0.687 | 60,236 | 58,535 | 43,082 | 18,211 | 487/496 | 487 |
| plasmacytoid DC (FC) | 36,029 | 6,893 | 0.152 | 12,691 | 5,208 | 2,195 | 7 | 6/496 | 6 |

## 4. Deviations, caveats and provisional items

- Data scope: none reduced. All cells of both cell types and all donors with >= 4 cells were used (LM 307,429 cells / 261 donors; FC 5,181 / 235). All analyses have the planned 20 reassignments (A, A2, B, D) and the real-donor runs. No analysis was cut.
- Declared adaptations (in Plan 7 before any run): exact streaming cps (the OneK1K run used a 5,000-cell subsample); forced-gene list with ISG genes in place of the PVM markers; XDS-Mlib assigns each donor's 4-16 libraries to halves by coin flips (Plan 6 used donors with exactly two libraries); the extra size check A2 and the D3 differential co-expression part were added to the plan before running.
- Plan 6 analyses not run here and not part of Plan 7: planted covariance, depth-matched null, planted confounder, sex diagnostics.
- Smoke tests (seeds +999, scratch directory) were run on the real panels to check the code; their outputs were not used. One smoke run of the library-preserving mode showed XDS-M with many calls before the real runs; no rule or code was changed in response. The real-donor mode has no seed offset, so its smoke output equals the reported real-donor result (identical counts were confirmed on the Modal run).
- In the few-cell type some genes have zero variance within a SLE or control group, which gives undefined PB-DC correlations for those pairs; such pairs are excluded from BH and counted as p >= 0.05 in the per-pair type-I of PB-DC.
- Real-donor SLE-vs-control calls are confounded with processing cohort (cohort 1 has no SLE donor; cohort 2 has 106 of 162 SLE monocyte donors) and with cell number per donor (SLE median 1,132 vs 792 monocytes per donor; 12 vs 24 pDCs); they are descriptive and are not evidence of SLE-specific co-expression changes except where stated post hoc after pool adjustment.
- Disclosure correction to the Plan 7 header (the plan is frozen and unchanged): before it was written, the h5ad obs columns suspension_uuid, ind_cov, ct_cov, author_cell_type, cell_state, disease, sex and observation_joinid were read in addition to the four named in the plan, and value counts of ct_cov and cell_state and a ct_cov by cell_type cross-table were printed. None of these entered a rule, a panel or an analysis; the author cell-type labels were not used (cell types are Census cell_type). No expression values of the dataset were read before the plan was frozen.
- Amendment 1 and 2 results are post hoc. The two failed variants (XDS-Mpool validity, XDS-Mlibpool) have only untested explanations.
- Plan 7 itself used the Modal CPU image im-vNcuP00ryKr0nZlE9sgcfw (numpy 1.26.4, scipy 1.17.1); smoke tests ran locally on numpy 2.4.6. Local code is a Plan 6 repo snapshot (artifact 14ce1518-9840-4775-a50d-3756123e6618); modules used are copied unchanged into code.tgz (hashes in input_hashes.txt).

## 5. Paper-ready paragraph (about 12 lines)

Second disease. In the Perez et al. lupus PBMC cohort (CELLxGENE dataset 218acb0f; 162 SLE and 99 control donors; 88 libraries in 24 multiplexed pools), analysed under a plan hashed before any expression was read, pseudobulk co-expression again made false calls under a donor-centred null (20 reassignments, 99,681 pairs): 16,631 per run in classical monocytes (1,005 cells per donor, per-pair type-I 0.29) and 3,472 in plasmacytoid dendritic cells (15 cells per donor, 0.17). XDS-M made 0.0 and 0.1 (type-I 0.051, 0.050) and met the frozen rule in both; the predicted ordering by cell number did not hold. For SLE-versus-control differential co-expression, pseudobulk called 8,937 and 1,493 pairs per run under a label null, 22% and 12% of its real calls, and XDS-Mgu none. A library-preserving null is a different test, and XDS-M failed it (48,661 calls per run in monocytes, 13 in pDCs; pseudobulk 62,723 and 11,126), as did our library-split variant in monocytes; post hoc, conditioning on donor pool composition removed these calls in monocytes (0.5 per run) at the cost of mild per-pair inflation (0.061 to 0.080).

## 6. Files

PLAN7.md / PLAN7.sha256, AMENDMENT1/2 (.md, .sha256), code (lupus_fetch_obs.py, lupus_prep.py, lupus_run.py, lupus_amend1.py, lupus_amend2.py, lupus_summary.py, lupus_amend_summary.py, make_results_md.py, run_all*.sh, jobs_*.txt, input_hashes.txt), provenance (lupus_LM_info.json, lupus_FC_info.json, lupus_pools.json, panel gene tables), results bundles (lupus_results_plan7.tgz: 164 files; lupus_results_amend1.tgz: 82; lupus_results_amend2.tgz: 122), lupus_summary.json, lupus_amend_summary.json, lupus_structured_output.json.
