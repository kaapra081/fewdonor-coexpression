# REV3_PREREG Amendment 1 (31 Oct 2026, ~20:40 UTC; before any A2x computation)

Reason: in the frozen A2 data, all 2,000 top-ranked OneK1K transitional B cell genes have mean > 0.1 UMI/cell (strata 265/625/1,110/0/0)
and naive B cells have only 50 genes in 0.03-0.1, so the low-expression strata named in the plan are not assessable. Only the first A2
run per null had been seen (B1, D = all: XDS-M S = 5 made 0 calls in both nulls; PB about 1.0-1.15 million of 2.0 million pairs).

A2x: OneK1K B1, genes ranked 2,001-3,000 by the same CP10k ranking on the same 5,000-cell subsample (seed 46,000,000); pairs among these
1,000 genes; same strata, methods, GPU implementation and float64 check as A2; 10 reassignments per null
(seeds 71,000,000 + 1,000 * null + 10 k + 1); D = 200 and all. Same decision rule as A2 with n = 10 runs (run-level bound 0.05 + 2 SE = 0.188).
