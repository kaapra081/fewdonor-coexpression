# Results for XS_PREREG, POP_PREREG and DSTD_PREREG (29 Oct 2026, 00:05 PDT)

Plan hashes: xsplit/XS_PREREG.sha256 (with Amendment 1), pop/POP_PREREG.sha256, dstd/DSTD_PREREG.sha256.
Naming used in the abstract: XDS = LCPXd in code (cross-split, design-only standardization).

## XS_PREREG (LCPX, cross-split with each donor's own variance) — NEGATIVE
- FDR (R = 100): F_D16 0.028, F_D16_n100 0.048, F_D24 0.087, F_D40 0.050, F_S1 0.191, F_S2 0.067 -> (a) fails.
- SEA-AD real-cell null: type-I 0.110 / 0.066 (top 1%: 0.153 / 0.074) -> (b) fails. Power 0.19 vs LCP 0.50 at D40 -> (c) fails.
- Decision: not adopted; weak-null claim for LCPX-S fails.

## DSTD_PREREG
- LCPd (design-only standardization, no split): FDR <= 0.05 + 2 MCSE in all six designs (F_S1 0.049); SEA-AD real-cell null FAILS
  (top-1% type-I 0.197 / 0.141); power 0.88x / 0.81x of LCP at D40 / D24. Not adopted. PsychAD E3 fails (mean false calls 38-1,218).
- LCPXd = XDS: FDR 0.025-0.049 in all six designs; SEA-AD real-cell null HOLDS (S1: 0 calls, type-I 0.050, top-1% 0.043;
  S2: 0 calls, 0.050, 0.050); PsychAD E3 HOLDS in all three subtypes (below); criterion (c) FAILS (power 0.018 vs 0.122 at D24,
  0.156 vs 0.501 at D40, 0.665 vs 0.929 at F_S2). By the frozen rule XDS does not replace LCP; it is the only method that holds
  every null tested.
- Claim "design-induced bias is a distinct artifact" as frozen (LCP vs LCPd): FAILS (LCPd all-pair type-I 0.058/0.059/0.118 > 0.055).
  The analogous LCPX vs LCPXd contrast (same split, standardization only differs): LCPX 15,538-32,801 false calls vs LCPXd 0-1.

## POP_PREREG (PsychAD, CELLxGENE 09642770..., 289,279 nuclei)
Subtypes: P1 Homeo_CECR2 (1,054 donors >= 4 nuclei, median 28), P2 PVM_CD163 (1,164, median 19), P3 ADAM_GPNMB (744, median 12).
Part A pseudo-donor null at full D (K = 3; P3 uses k = 1..3, k = 0 was the exploratory reassignment):
| subtype | PB calls | PB type-I | LCP calls | LCP type-I | XDS calls | XDS type-I | XDS top-1% |
| P1 | 19,158 | 0.335 | 4,986 | 0.201 | 0 | 0.050 | 0.053 |
| P2 | 31,617 | 0.441 | 12,880 | 0.301 | 0.67 | 0.050 | 0.051 |
| P3 | 24,544 | 0.378 | 9,125 | 0.265 | 0 | 0.048 | 0.047 |
At D = 200: PB 7,036 / 10,013 / 10,443; LCP 137 / 145 / 566; XDS 1 / 0 / 0.
Replicated false calls across pseudo-donor halves: PB 3,798-12,827; LCP 502-2,285; XDS 0 in all 9.
- "LCPX holds at population scale": NO. "Replication does not protect pseudobulk": YES (all >= 50).
- "PB false discoveries grow with cohort size" (frozen: calls at full D >= 3x calls at D = 100 and >= 100; k = 0 for P1/P2, k = 1 for P3): YES in all three.
  PB calls at D = 25 / 50 / 100 / full: P1 1,778 / 2,578 / 3,929 / 19,267; P2 1,063 / 3,320 / 6,288 / 31,571; P3 1,409 / 1,803 / 4,031 / 24,962 (ratio full/100: 4.9, 5.0, 6.2).
  LCP calls at D <= 100: 0 except 19 (P2, D = 100); at full D 4,610-13,707. XDS: 0 at every D.
Part B real donors (full D): XDS calls 32,790 / 49,852 / 55,021 (37 / 56 / 63% of pairs); PB-only calls 13,307 / 10,086 / 14,643;
replication in independent donor halves: XDS 79 / 82 / 68%, PB 73 / 80 / 62%.
Caveat: pooling real donors into pseudo-donors turns real between-donor covariance into within-pseudo-donor covariance, so this
null's within-donor structure is somewhat stronger than in real donors.
