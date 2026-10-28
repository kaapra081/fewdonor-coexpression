# Rare-supertype analysis: outcomes under the frozen rule (RARE_PREREG.md, commit 38765b5, sha256 4e389922…)

## Calibration gate (S1)

- Simulation matched to S1's exact cell numbers: D = 79, median n_d = 15, R = 100.
- LCP: FDR 0.075 (MCSE 0.004); power 0.472; type-I for within-only pairs 0.086.
- PB: FDR 0.549; power 0.399; type-I for within-only pairs 0.450.
- **The gate FAILS** (0.075 > 0.05 + 2·MCSE). By the rule, the S1 Q1 and Q2 results are descriptive.

## S1 (Micro-PVM_1, primary): D = 79, nuclei = 1361, median n_d = 15

**Calls.** PB 9970, LCP 8698; PB-only 3926, both 6044, LCP-only 2654.

**Q1 = NO** (cond b False, cond c True).
- Median |r_w| in REF cells: PB-only 0.050, both 0.058, LCP-only 0.045.
- Mann–Whitney p: vs both 1; vs LCP-only 1.25e-07.
- REF replication: PB-only 0.39, both 0.63, LCP-only 0.51.
- Median leak share: PB-only 0.130, both 0.049.

**Q2 = YES.** Top-20 hub Jaccard 0.33; random removal median 0.67, 5th percentile 0.48.
- Hubs dropped: IGSF21, LAPTM5, SLC25A37, SFMBT2, SORL1, FGD4, SPRED1, SNX29, SSH2, PTPN2.
- Hubs added: SLC11A1, DOCK4, ADAM28, FTH1, ATM, CHST11, HIF1A, MAF, EPB41L3, RAB31.

**Descriptive, not in the rule.**
- Replication in Lee 2026: PB-only 0.42, both 0.53, LCP-only 0.42, neither 0.29.
- Median |r_w| in the tested cells themselves (not independent of the test): PB-only 0.106, both 0.097, LCP-only 0.071.

## S2 (Micro-PVM_3, secondary): D = 83, nuclei = 4518, median n_d = 43

**Calls.** PB 23215, LCP 21427; PB-only 4176, both 19039, LCP-only 2388.

**Q1 = NO** (cond b False, cond c True).
- Median |r_w| in REF cells: PB-only 0.047, both 0.058, LCP-only 0.045.
- Mann–Whitney p: vs both 1; vs LCP-only 0.00176.
- REF replication: PB-only 0.45, both 0.81, LCP-only 0.67.
- Median leak share: PB-only 0.042, both 0.025.

**Q2 = NO.** Top-20 hub Jaccard 0.67; random removal median 0.74, 5th percentile 0.60.
- Hubs dropped: MEF2C, SHTN1, CYFIP1, NRIP1.
- Hubs added: DOCK10, APBB1IP, A2M, PTPN2.

**Descriptive, not in the rule.**
- Replication in Lee 2026: PB-only 0.45, both 0.62, LCP-only 0.44, neither 0.30.
- Median |r_w| in the tested cells themselves (not independent of the test): PB-only 0.072, both 0.080, LCP-only 0.063.

## S2 calibration simulation (not in the plan)

- R = 100.
- LCP: FDR 0.051 (MCSE 0.003); type-I for within-only pairs 0.060.
- PB: FDR 0.169; type-I for within-only pairs 0.382.

## Deviations

- The S2 calibration simulation was added after the plan was frozen; it is labelled descriptive.
- No other deviations.
