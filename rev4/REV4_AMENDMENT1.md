# REV4 Amendment 1 (1 Nov 2026, ~19:55 UTC; before any R6 computation): R6, explanatory maps in Mantis

Motivation: show directly what a pseudobulk co-expression network is made of. Descriptive only; no pass/fail decision. Mantis (Kellis Lab
semantic-mapping platform) is used to lay out and explore the maps; every number reported is computed outside Mantis from the tables uploaded.

R6a. Donor maps. (i) PsychAD P1 donor-centred null, R3 reassignment k = 0 (seed 61,000,000): 1,054 pseudo-donors that keep real cell numbers
     and have no between-donor co-expression by construction. (ii) The 1,054 real P1 donors. For each, two profiles over the 424 genes: the centred
     pseudobulk log1p(1e6 m) profile and the design-standardized residual profile R (design_R, as used by XDS/XDS-M). Map layout: precomputed 2-D UMAP
     (umap-learn, n_neighbors 15, min_dist 0.1, random_state 0) of each profile; overlays: cells per donor n_d, mean UMI depth, brain bank (real only),
     profile norm. Reported: Spearman(profile norm, n_d); cross-validated R^2 (5-fold) of log n_d from a 15-nearest-neighbour regression on each profile.
R6b. Pair level. Spearman correlation, over all gene pairs, between the null pseudobulk correlation r_PB^null (R6a-i) and the mean within-donor cell-level
     correlation r_w of real donors (donors with >= 10 cells, as R4); the same for the real-donor pseudobulk correlation and for the XDS-M statistic.
R6c. Gene map. P1 genes. Layout: Mantis semantic embedding of the gene's reviewed human UniProt protein name and function text (no result enters the text).
     Overlays: real-donor PB degree, XDS-M degree, PB-only degree (BH q = 0.05, as R4), null PB false degree (R6a-i), mean UMI per cell, detection fraction,
     mean |r_w|. Reported: Spearman(null false degree, PB-only degree), Spearman(null false degree, mean |r_w|), Spearman(null false degree, mean UMI);
     with Mantis's automatic clusters of the gene map (membership read back from Mantis), Kruskal-Wallis tests of PB-only degree and XDS-M degree across
     clusters, and the clusters with the highest median PB-only degree. Genes without a reviewed UniProt entry are labelled by symbol only and reported.
