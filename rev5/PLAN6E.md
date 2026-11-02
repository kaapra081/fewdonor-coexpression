# Plan 6e (2 Nov 2026, ~18:00 UTC). POST HOC, after Plan 6d (XDS-Mh met the 6d rule and becomes the recommended test).
# The recommended recipe with covariates combines the two post hoc corrections exactly as the paper's Algorithm 1 will state:
# XDS-Mhc: for every split and half h, (1) residualize the half ratio means m^h on the donor covariates W (OLS across donors),
# (2) divide by sqrt(v*^h) (design-only), (3) residualize the result on the half's own design W^h = [1, log rhat^h, log n^h],
# (4) cross-products, average over the 5 splits, one-sample t test. Code: rev5/rev5_e.py.
Evaluations (same seeds): (a) planted-confounder null (Plan 6 C), W = all C1 covariates; (b) depth-matched null (Plan 6 D) with the
same W attached to the pseudo-donors' labels; (c) donor-centred null (Plan 6c a), full size and D = 200, same W; (d) real donors
P1-P3: calls, fraction of XDS-M calls kept, and the 7 XIST-Y and 21 Y-Y pairs.
Rule: XDS-Mhc is stated as validated with covariates if <= 2 of 20 runs with >= 1 false call and mean per-pair type-I <= 0.055 in
(a) (any pair; type-I over planted pairs), (b) and (c) in all 3 subtypes; otherwise the paper says the combination is untested/failed.
