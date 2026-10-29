# Differential co-expression (AD vs control) confounded by cell number: frozen plan (29 Oct 2026, 11:40 PDT)

Seen before this plan: only the per-donor nucleus counts by AD status (Mann-Whitney p 8e-5 to 5e-8 in the three subtypes;
AD-vs-control medians 21 vs 36 Homeo_CECR2, 23 vs 18 PVM_CD163, 17 vs 10 ADAM_GPNMB). No co-expression by AD status computed.

Data: PsychAD subtypes P1-P3 (POP_PREREG), donors with AD label 0 or 1 (CELLxGENE column `AD`) and >= 4 nuclei.

A. Null. K = 3 pseudo-donor reassignments per subtype (seeds 41,000,000 + 10 k + subtype index, k = 10, 11, 12): nuclei pooled and
   reassigned with every real donor's nucleus count kept; each pseudo-donor carries the AD label of the real donor whose count it
   has. Between-donor covariance is zero in both groups, so every differential co-expression (DC) call is false.
B. Real donors: the same tests on real donors.

Tests per pair (BH q = 0.05 over all pairs):
- PB-DC: Fisher z difference of log-CPM pseudobulk Pearson correlations, z = (atanh r1 - atanh r0)/sqrt(1/(n1-3) + 1/(n0-3)).
- LCP-DC: Welch z for the difference in group means of the LCP per-donor terms z_d = R_dj R_dk - Lt_d (paper's LCP standardization).
- XDS-DC: Welch z for the difference in group means of q_d = R*A_dj R*B_dk (XDS halves; mu, v* pooled over both groups).

Decisions:
- "Pseudobulk DC is confounded by cell number" = YES if PB-DC mean false calls >= 100 in at least two subtypes.
- "XDS-DC holds" = YES if in every subtype XDS-DC mean false calls <= 1 and all-pair type-I <= 0.055.
Real-donor results are reported descriptively (number of DC calls per method, overlap).
