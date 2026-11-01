# Fourth-round (external AI review) analyses: frozen plan (1 Nov 2026, ~19:30 UTC, before any computation below)

Context: an external AI program-committee review of the full draft (1 Nov 2026). It pointed out that with B = 10^5 permutations the smallest
permutation p-value is about 1e-5, while BH over ~88,000 pairs needs p <= 5.7e-7 for a first call, so permutation-based tests (XDS for D <= 100,
LCP and LCPd for D <= 100) can make a BH call only if >= 18 pairs tie at the floor. Their small-cohort "no false call" results are therefore
not evidence of calibration. Nothing below has been computed. Run-level rule as before: a method holds in a cell if the fraction of runs with
>= 1 BH call (q = 0.05) <= 0.05 + 2 SE (binomial, p = 0.05): n = 20 -> <= 0.147 (2 of 20); n = 40 -> <= 0.119 (4 of 40).

R1. Permutation resolution at small D. PsychAD P1-P3, donor-centred null exactly as REV3 A3 (clip [0.1, 10]), nested D = 25, 50, 100;
    20 reassignments per subtype (seeds 63,000,000 + 10 k + subtype index). Permutation-based tests (LCP, LCPd, XDS-perm) with
    B = 2,000,000 permutations (smallest p = 5e-7 < q/m); continuous-reference tests (XDS normal-moment, XDS-M, KX) recomputed on the same runs.
    Also recorded per run and method: number of pairs with p <= 1e-5 (the B = 1e5 floor). Decision per method, subtype, D by the run-level rule.
    If compute does not allow 20 runs per subtype by Tue 3 Nov 12:00 UTC, report the runs completed (stated as such).

R2. Batch-stratified null (Lemma 1 failure mode: donors differ in depth distribution). PsychAD P1-P3, cells donor-centred as REV2 A,
    then pseudo-donor labels permuted only among cells of the same brain bank (obs column 'Source': M, H, R), so each pseudo-donor keeps a real
    donor's cell number and inherits its bank's depth and composition distribution; donor-centring removes donor (and bank) mean effects.
    40 reassignments per subtype at full D (seeds 64,000,000 + 10 k + subtype index). Methods PB, LCP, LCPd, XDS, XDS-M, KX (normal references).
    Also reported: median library size and median cells per donor by bank. Decision: XDS-M holds by the run-level rule (n = 40) in each subtype.

R3. More pseudobulk baselines. PsychAD P1-P3, donor-centred null, 20 reassignments at full D with the REV3 B2 seeds (61,000,000 + 10 k + si):
    PC-adjusted pseudobulk (remove the top k = 5 and k = 10 principal components of the centred donor x gene log1p(1e6 m) matrix, Pearson,
    Fisher z with D - 3 - k), alongside PB, PBres and PBw. OneK1K pooled null, 10 reassignments per cell type (seeds 65,000,000 + 10 k + type
    index): PB, PBw, PC-adjusted k = 5, 10. Descriptive (mean BH calls, runs with >= 1 call, per-pair type-I); no decision.

R4. Real-donor consequence. PsychAD P1-P3 real donors: PB and XDS-M (S = 5) BH calls and their overlap; a between-donor correlation estimate
    from unstandardized cross-split moments averaged over the 5 splits, rho_jk = C_jk / sqrt(C_jj C_kk) with C from (m^A - mean)(m^B - mean);
    the mean within-donor cell-level correlation r_w of each pair (normalized expression, donors with >= 10 cells). Reported: number of
    pseudobulk-only, XDS-M-only and shared calls; median r_w and median (r_PB - rho) in each class; Spearman correlation of (r_PB - rho) with r_w;
    overlap of the 20 highest-degree genes under PB and XDS-M at matched numbers of edges (top 5,000 pairs by p-value). Descriptive, post hoc in
    spirit (no decision), split seeds 44,000,001.

R5. Differential co-expression null. REV2 D1 donor-centred null with 17 more reassignments per subtype (k = 13..29; same seed formula),
    giving 20 in total; methods PB, XDS-g, XDS-Mg, XDS-Mgu. Decision: XDS-Mgu holds by the run-level rule (n = 20) in each subtype.
