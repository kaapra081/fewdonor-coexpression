# REV4 Amendment 4 (2 Nov 2026, ~01:20 UTC). Written after C1 (P1), C2 (P1 30 runs, P3), C3, C5 and C6 had been seen, in particular after
# seeing that the heat-shock leaf's top rank in Homeo_CECR2 did NOT replicate in PVM_CD163 / ADAM_GPNMB (ranks 21 and 17 of 49).
# Before any computation below. All analyses below are post hoc and will be labelled so.

Question: are pseudobulk's false edges organized by gene function in the same way as real between-donor co-expression, i.e. would a
functional-coherence check (as commonly used to validate co-expression networks) distinguish them?

D1. Function partition: the Mantis gene map of the 424 Homeo_CECR2 genes (rev4/mantis/gene_map_mantis.csv), built from UniProt/NCBI function
    text only. Primary: the 49 leaf clusters. Secondary: Mantis level 4 (10 clusters); the 19 genes without a level-4 node keep their level-3
    cluster. For PVM_CD163 and ADAM_GPNMB only genes present in the Homeo_CECR2 map are used (others dropped, with their edges).
D2. Edge sets per subtype (P1, P2, P3), from rev4_R6.py --sub <s> (null = first R3 donor-centred reassignment, seed 61,000,000 + si; real
    donors with split seed 44,000,001; BH 5%): (i) null pseudobulk calls (all false), (ii) real pseudobulk calls, (iii) real XDS-M calls,
    (iv) pseudobulk-only = (ii) and not (iii), (v) XDS-M-only = (iii) and not (ii). P2 and P3 are run now with the unchanged script.
D3. Statistics per edge set and partition: number of edges m; fraction within cluster; Newman modularity Q = sum_c [e_c/m - (d_c / 2m)^2]
    (configuration-model expectation, i.e. adjusted for gene degree); gene-label permutation p-values for both (10,000 permutations of the
    cluster labels over genes, cluster sizes preserved, seed 20,261,103; p = (1 + #{perm >= obs}) / 10,001); Q as a z-score against the
    permutation distribution.
D4. Expression check (P1): per-gene degree (null false, pseudobulk-only, XDS-M) residualized by least squares on log mean UMI per cell and
    detection fraction; eta-squared of the residual across the 49 leaves against 10,000 label permutations (seed 20,261,104).
D5. Interpretation stated now. "False edges are function-coherent" if the null pseudobulk edge set has permutation p < 0.01 for Q on the leaf
    partition in at least 2 of 3 subtypes; "not function-coherent" otherwise. Either way we report Q for real XDS-M and pseudobulk edges for
    comparison, and say plainly that a functional-coherence check cannot certify an edge set as between-donor signal if the false edges pass it.
