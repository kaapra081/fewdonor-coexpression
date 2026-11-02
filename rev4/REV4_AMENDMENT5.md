# REV4 Amendment 5 (2 Nov 2026, ~02:25 UTC). Written after internal review round 6 (review_round6.md, items N1, N2, N4) and after the
# Amendment 4 results had been seen; before any computation below. All analyses are post hoc and will be labelled so.
# Also recorded: REV4 R1 used B = 2,000,000 donor permutations (rev4.py default; REV4_B unset in the running processes, checked
# 02:27 UTC from /proc/<pid>/environ on the compute host).

E1. Robustness of the function-modularity result (Amendment 4).
 (a) Expression-stratified relabelling. Kept genes are binned into 10 equal-count bins of log mean UMI per cell (per subtype); the
     10,000 relabellings (seed 20,261,105) permute cluster labels only within bins. Statistics as Amendment 4 D3, for every edge set,
     leaf partition and level-4 partition.
 (b) All 20 donor-centred null reassignments per subtype (R3 seeds 61,000,000 + 10 k + si, k = 0..19): null pseudobulk BH calls;
     stratified Q z-score (1,000 stratified relabellings per reassignment, seed 20,261,106 + k); report mean and SD of z over the 20 and
     the number of reassignments with stratified p < 0.01.
 (c) Equal-size edge sets: the 5,000 strongest pairs of (i) null pseudobulk (|r|, reassignment k = 0), (ii) real pseudobulk (|r|),
     (iii) real XDS-M (smallest p, ties broken by |z|); stratified Q z for each.
 (d) Second, non-Mantis annotation: GO biological-process terms of each gene's reviewed human UniProt entry (retrieved now; non-coding
     genes have none). Terms annotated to 2..50 of the kept genes are used. Statistic: fraction of edges whose two genes share at least one
     such term; compared with the stratified relabelling of gene identities (10,000; seed 20,261,107); ratio to the permutation median
     and p, for every edge set of (a) and (c).
 (e) Partition check: fraction of within-leaf pairs (all pairs, not edges) sharing a GO term as in (d), against all pairs, with the
     stratified relabelling p.
 Rule stated now. The claim that pseudobulk's false edges are function-coherent stays in the abstract and contributions only if
 (a) null pseudobulk stratified p < 0.01 on the leaf partition in >= 2 of 3 subtypes AND (b) mean stratified z over the 20 reassignments
 > 2.33 in >= 2 of 3 subtypes. Otherwise it is removed from the abstract and contributions and Section 5.2 reports the result as not
 robust. If (d) also gives p < 0.01 for null pseudobulk edges in >= 2 subtypes, the text may say "by Mantis clusters and by GO";
 otherwise only the Mantis partition is named. The comparison with real edge sets is made only at equal size (c).
E2. Positive control for covariate adjustment (C1): on real donors (split seed 44,000,001), P1-P3, the 7 XIST-Y and 21 Y-Y pairs:
    number called (BH over all pairs) and median p by XDS-M, XDS-Mc (all covariates of C1) and XDS-M adjusted for sex only. Expected:
    called by XDS-M, not by the sex-adjusted versions. Descriptive.
E3. For P2 and P3, the odds ratio that a real pseudobulk-only call is a false call in the first null reassignment (as C3 for P1). Descriptive.
