# REV2_PREREG Amendment 1 (31 Oct 2026, ~20:35 UTC; written before any REV2 part D computation)

Reason: the internal code audit (synthetic data) showed that XDSg and XDS-Mg compare means of *design-standardized*
cross-split products. When the between-donor covariance of a pair is non-zero, that mean is Sigma^B_jk / sqrt(v_j v_k), and v depends on
cells per donor, so two groups with equal covariance but different cell numbers differ in mean product (synthetic: mean q 0.258 at 15
cells/donor vs 0.402 at 40). XDSg/XDS-Mg therefore remain valid tests of "no between-donor covariance in either group" but are not tests of
equal covariance when groups differ in cell number, which AD cases and controls do.

Change: add one method to D1 (null and real), D2 (true-DC simulation) and D3 (split halves): XDS-Mgu = Welch test of the group difference
in the mean of the UNSTANDARDIZED group-centred cross-split product (m^A_dj - mean_g m^A_j)(m^B_dk - mean_g m^B_k), averaged over the same
S = 5 splits. Its expectation within a group is the between-donor covariance whatever the cell number (disjoint halves remove the leak).
All frozen methods (PB, XDSg, XDS-Mg) are still run and reported. Decision rules for XDS-Mgu are the same as those frozen for XDS-Mg.
