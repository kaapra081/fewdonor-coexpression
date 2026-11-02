# Plan 6c (2 Nov 2026, ~17:45 UTC). POST HOC (written after Plan 6 results; before any Plan 6b or 6c result was seen except the P3
# smoke run of 6b, k = 0). Purpose: decide whether a design-adjusted XDS-M should replace XDS-M as the recommended test, given the
# Plan 6 D-a failure. Variants: XDS-Md (Plan 6 D: R residualized on [1, log rhat_d, log n_d]) and XDS-Mda (Plan 6b H1: half means
# residualized on the same W first, then design-standardized). Code: rev5/rev5_c.py.
(a) Donor-centred null at full size, the 20 REV4-R3 reassignments (seed 61,000,000 + 10 k + si): BH calls, per-pair type-I of
    XDS-M, XDS-Md, XDS-Mda. (b) Subsets of the same reassignments of D = 100 and D = 200 pseudo-donors (drawn with seed + 1).
(c) OneK1K pooled null, B1-B3, the 10 REV4-R3o reassignments (seed 65,000,000 + 10 k + ti), raw cells.
(d) Planted-signal null of Plan 6 S at full size (same seeds): FDR, MCSE and power of XDS-Md and XDS-Mda.
Rule. XDS-Mda is recommended if it has <= 2 of 20 runs (<= 1 of 10 in (c)) with >= 1 false call and mean per-pair type-I <= 0.055 in
Plan 6 D (via 6b H1b), (a), (b) at D = 200, (c), in every subtype/cell type, and FDR <= 0.05 + 2 MCSE in (d). Otherwise XDS-Md if it
meets the same rule; otherwise XDS-M remains recommended with the depth caveat. D = 100 is reported descriptively.
