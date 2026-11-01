# REV4 Amendment 2 (1 Nov 2026, ~20:50 UTC). POST HOC, written after seeing partial R2 results (XDS-M 8/40 runs with a false call in PVM_CD163)

Observation that motivates it: the R2 null is not exactly null if brain banks differ in residual mean expression after donor-centring. The centring
factor uses the shrunken donor effect, so for donors with few cells most of a bank-level shift survives; pseudo-donors of the same bank then share a
bank-level mean, which is real between-pseudo-donor covariance (rank <= 2 over three banks) that a valid test should detect. R2 therefore cannot by
itself separate miscalibration from residual bank structure.

R2c (diagnostic, post hoc, no decision rule): same 40 R2 reassignments per subtype (seeds 64,000,000 + 10 k + si). Methods: XDS-M as in R2; XDS-Mb =
XDS-M with the design-standardized half-donor values R^A, R^B centred within brain bank before forming the per-donor products (as XDS-Mg centres within
AD group); PB and PBb (log-CPM centred within bank, Fisher z with D - 3 - 2). Reported: runs with >= 1 BH call and mean calls, per method and subtype;
also per bank the mean of the donor-centred log-CPM shift (median over genes of |bank mean - overall mean| / overall SD of pseudo-donor log-CPM).
Interpretation rule stated in advance of R2c: if XDS-Mb has at most the bound's number of failing runs where XDS-M fails, the R2 failures are attributed
to residual bank structure (and the paper recommends centring within batch); otherwise they are reported as miscalibration under bank-level depth
differences.
