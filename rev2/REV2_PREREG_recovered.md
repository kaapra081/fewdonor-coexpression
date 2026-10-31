<!-- Recovered from the session transcript because the canonical file (~/Documents/recomb27/recomb_v2/rev2/REV2_PREREG.md, sha256 1a7e4de1...) is temporarily unreadable. Content should be verbatim but whitespace may differ, so this copy does not reproduce the frozen hash. -->
# Second-round review analyses: frozen plan (29 Oct 2026, 14:30 PDT, before any computation below)

A. Donor-centred null, run-level error, by cohort size. PsychAD P1-P3; donor effects removed before pooling exactly as REV_PREREG R6;
   100 reassignments per subtype (seeds 50,000,000 + 10 k + subtype index, k = 0..99); nested random subsets D in {200, 400, all}.
   Methods PB, LCP, LCPd, XDS (single split, normal-moment reference), XDS-M, KX. Per run: BH calls V.
   Decision per method, subtype and D: \"holds\" if the fraction of runs with V >= 1 is <= 0.05 + 2 SE (binomial, n = 100).
   Primary question: does LCPd hold where XDS-M holds? If LCPd holds at every D and subtype, it is reported as a co-recommended test.
B. Variance-heterogeneity diagnostic (A5), descriptive: in the pooled PsychAD null (POP seeds, k = 0), mean over genes of (R^A_dj)^2 for
   XDS's design-only standardization, by decile of n_d.
C. Real donors with XDS-M (P1-P3): BH calls; replication in random donor halves (seed 42,000,000 halves as POP); three split seeds
   (call counts, Jaccard); XIST-Y and Y-Y pairs for PB and XDS-M.
D. Differential co-expression.
   D1. Donor-centred DC null (DC_PREREG design on the donor-centred cells; k = 10, 11, 12): PB-DC, XDS-DCg, XDS-M-DCg
       (XDS-M with group-specific centring: per-donor mean over 5 splits of group-centred cross-products; Welch z).
   D2. True-DC simulation: group 0 = D40 settings with D = 300, nd_median 15; group 1 = D16_strong between-program strength (g = 0.45)
       with D = 300, nd_median 40; held-out panel; seeds 9,800,000 + rep, R = 50. True DC pairs: |SB_1 - SB_0| > 0 by the analytic truth.
       Report FDR and power (BH 0.05) of PB-DC, XDS-DCg and XDS-M-DCg.
   D3. Split-half replication of real-donor XDS-DCg and XDS-M-DCg calls (P1; random donor halves, seed 42,000,000).
E. OneK1K: 30 further pooled reassignments per cell type (k = 3..32), full D; PB, PBlin, LCP, XDS-M; run-level fraction for XDS-M over
   all 33 runs; PBlin versus PB to test the transformation-bias explanation (descriptive).
F. A5 stress simulation: L200 and L500 with a donor-specific overdispersion factor shared by all genes (per-donor gamma shape
   phi * exp(0.7 z_d), z_d ~ N(0,1)); seeds 9,900,000 + 100,003 k + rep, R = 50; PB, LCP, LCPd, XDS, XDS-M, KX FDR and power.
   No decision; reported as the limit of (A5).
G. Oracle (matched-FDP) power from existing simulation p-values: per replicate, the largest top-k set with FDP <= 0.05; descriptive.
