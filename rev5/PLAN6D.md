# Plan 6d (2 Nov 2026, ~17:50 UTC). POST HOC, written after the Plan 6c results: by the 6c rule neither XDS-Md nor XDS-Mda is
# recommended (XDS-Md: 8/10 OneK1K-B1 pooled-null runs with false calls; XDS-Mda: B1-B3 and D = 200 failures), so XDS-M stays the
# recommended test unless the variant below passes. Diagnosis: XDS-Md and XDS-Mda adjust both halves on rhat_d computed from ALL of
# the donor's cells; with depth-state coupling, the depths of half-A cells carry information on half-A states, so a full-donor
# regressor couples the two halves. Variant XDS-Mh: each half's R^h is residualized on its OWN design, W^h_d = [1, log rhat^h_d,
# log n^h_d], rhat^h_d = sum_{c in h} s_c^2 / sum_{c in h} s_c (= kappa^h_d S^h_d), for every split. It still removes the term
# shared by both halves through F_d (Lemma 1(ii)) but uses no information from the other half. This is the third and last variant
# tried; all three are reported. Code: rev5/rev5_d.py.
Evaluations (same seeds as before): (a) depth-matched null (Plan 6 D); (b) donor-centred null at full size and D = 200, 100
(Plan 6c a, b); (c) OneK1K pooled null B1-B3 (Plan 6c c); (d) planted-signal null at full size (Plan 6 S); (e) real donors P1-P3,
B1-B3: calls and the fraction of XDS-M calls kept.
Rule (as Plan 6c): XDS-Mh is recommended if <= 2 of 20 runs (<= 1 of 10 in (c)) with >= 1 false call and mean per-pair type-I
<= 0.055 in (a), (b) at full size and D = 200, and (c), in every subtype/cell type, and FDR <= 0.05 + 2 MCSE in (d). Otherwise
XDS-M remains the recommended test and donor-to-donor depth heterogeneity is stated as an unresolved limitation.
