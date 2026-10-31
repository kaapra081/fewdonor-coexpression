# Third-round review analyses: frozen plan (31 Oct 2026, ~20:15 UTC, before any computation below)

Context: internal review round 3 (review_round3.md). Interim REV2-A results (P3 100 runs, P2 40 runs) were seen before this plan;
nothing below has been computed. All runs use real public data; seeds are new and listed. Run-level rule as in REV2:
a method "holds" in a cell if the fraction of runs with >= 1 BH call (q = 0.05) <= 0.05 + 2 SE (binomial, SE at p = 0.05).

A2. Tail calibration by gene-expression stratum (validated gene range).
   Data: OneK1K (CELLxGENE Census 2025-11-08) transitional stage B cells (B1) and naive B cells (B2), donors with >= 4 cells;
   genes ranked 1-2,000 by mean CP10k on a 5,000-cell random subsample (seed 46,000,000, as ONEK_PREREG; disclosed: ranking subsample).
   If the PsychAD h5ad becomes readable, ADAM_GPNMB (P3) with the same 2,000-gene rule is added as A2-P3 (same seeds offset +3).
   Gene stratum by mean raw UMI per cell over all cells of the type: s1 > 1, s2 0.3-1, s3 0.1-0.3, s4 0.03-0.1, s5 < 0.03.
   Pair stratum = stratum of the lower-expressed gene of the pair.
   Nulls: pooled reassignment (POP construction) and donor-centred reassignment (REV_PREREG R6 construction, clip [0.1, 10]).
   20 reassignments per data set and null (seeds 70,000,000 + 1,000 * null(0 pooled, 1 centred) + 10 k + dataset(1 = B1, 2 = B2, 3 = P3)),
   at D = 200 (random subset, rng of the same seed) and D = all.
   Methods: PB (Fisher-z of log1p(1e6 m) Pearson), XDS (split 0, Hoeffding-moment normal reference), XDS-M with S = 1, 5, 10 splits
   (one-sample t of per-donor mean cross-products, normal reference; splits as revm.stats_from_cells: split 0 seed + 7, split s seed + 100 + s).
   Implementation: exact matrix-product form of the same statistics on GPU (float32), checked in every run against the float64 CPU code
   (revm.methods) on 2,000 random pairs; a run is invalid if max |log10 p_gpu - log10 p_cpu| > 0.05 for p > 1e-8.
   Outputs per run, method, D, stratum: number of pairs; number with p <= 1e-3, 1e-4, 1e-5, 1e-6, 1e-7 (and how many positive);
   BH calls (q = 0.05) over all pairs and within stratum.
   Decision (per data set, null, D, method, stratum): "tail-calibrated" if the pooled-over-runs observed/expected ratio at p <= 1e-4
   and at p <= 1e-5 are both <= 1.5 (strata with expected count < 20 at 1e-5 are "not assessable" at that threshold), AND the
   fraction of runs with >= 1 within-stratum BH call holds by the run-level rule (n = 20: <= 0.147).
   The validated gene range for XDS-M (S = 5) is the set of strata that are tail-calibrated in every data set, null and D.

A3. Run-level error at small D. PsychAD P1-P3, donor-centred null exactly as REV2-A; 100 reassignments per subtype
   (seeds 60,000,000 + 10 k + subtype index); nested random subsets D = 25, 50, 100. Methods PB, LCP, LCPd, XDS (permutation
   p-values as in DSTD for D <= 100: LCPXd), XDS-M, KX. Decision per method, subtype, D by the run-level rule (n = 100: <= 0.0936).
   If the compute budget does not allow 100 runs per subtype by Tue 3 Nov 12:00 UTC, report the runs completed (stated as such).

B2. Stronger pseudobulk baselines. PsychAD P1-P3, donor-centred null (seeds 61,000,000 + 10 k + subtype index), 20 reassignments,
   D = all: (i) PBres: residuals of each gene's log1p(1e6 m_d) on [1, log n_d, log mean depth_d], then Pearson/Fisher-z (df D - 5);
   (ii) PBw: weighted Pearson of log1p(1e6 m_d) with weights n_d (Fisher-z with Kish effective D). Report mean BH calls and
   all-pair type-I. Decision: "covariate adjustment does not fix pseudobulk" if mean calls >= 100 for both variants in all subtypes.

B3. Donor-centring sensitivity. P3, donor-centred null as REV2-A but shrinkage-factor clip [0.5, 2] instead of [0.1, 10];
   40 reassignments (seeds 62,000,000 + 10 k + 2); D = 200, 400, all; methods PB, LCPd, XDS, XDS-M, KX. Descriptive comparison with
   REV2-A P3 (no decision).

Priority on the shared 4-CPU host: REV2-A (P2, P1) > A2 (GPU) > REV2 C, D, E > B2 > A3 > B3 > REV2 D2, F > REV2-A P3 raw-file rerun.
