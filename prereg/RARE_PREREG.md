# Rare-supertype analysis: frozen decision rule (written before any rare-supertype statistic was computed)

Date: 28 Oct 2026. Data: SEA-AD MTG Microglia-PVM, CELLxGENE dataset c76098ba-eed3-45b1-98f2-96fcac55ed18
(version ab53c08f-03a3-4e0e-973a-44beabd78f60), 10x 3' v3 nuclei, ADNC != Reference (84 donors, 34,794 nuclei).
Only per-donor cell counts by supertype have been looked at so far (to choose the datasets below).

## Datasets
- **S1 (primary): Micro-PVM_1.** 1,367 nuclei; 79 donors with >= 4 nuclei. Few cells per donor occur naturally here; nothing is down-sampled.
- **S2 (secondary): Micro-PVM_3-SEAAD.** 4,518 nuclei; 83 donors with >= 4 nuclei; median about 44 per donor.
- **REF: Micro-PVM_2.** 19,130 nuclei from the same 84 donors; median 174 per donor. These are disjoint cells, used only as an independent reference.

Donors with fewer than 4 nuclei of a supertype are excluded from that supertype's analysis.

## Genes
- Use the 417 genes of the main analysis (fixed earlier). Within a dataset, drop any gene with total UMI count < 20, and all pairs that involve it.
- Library size is the UMI total over all 35,483 genes, as in the main analysis.

## Tests
- **LCP:** the frozen method; B = 10^5 donor permutations for all pairs, no refinement; seed 28102026.
- **PB:** pseudobulk log-CPM Pearson, normal theory.
- Both use BH at q = 0.05.

Call sets in S1 and S2:
- PB-only = called by PB and not by LCP;
- Both = called by both;
- LCP-only = called by LCP and not by PB.

## Calibration gate
- Simulate the S1 design: D = number of S1 donors; per-donor cell numbers equal to S1's observed n_d; fitted depths; development gene-panel settings; R = 100; seeds 7,000,000 + rep.
- If LCP's FDR > 0.05 + 2·MCSE, the S1 results are reported as "LCP not calibrated at these cell numbers", and Q1/Q2 are interpreted only descriptively.

## Q1. Are PB-only pairs within-donor programs?
Within-donor reference: the within-donor latent correlation r_w,ref of each pair, computed in REF cells (disjoint from S1/S2, same donors).

Q1 = **YES** if all of the following hold. Otherwise NO (or "not assessable" if (a) fails):
- (a) |PB-only| >= 50.
- (b) Median |r_w,ref| of PB-only pairs exceeds that of Both pairs and that of LCP-only pairs. Each comparison is a one-sided Mann–Whitney test with p < 0.01. If LCP-only has < 50 pairs, the comparison with it is dropped.
- (c) Replication as between-donor signal in REF: the fraction of PB-only pairs that are REF LCP calls of the same sign is at least 0.10 lower than the same fraction for Both pairs.

## Q2. Does removing the PB-only pairs change a downstream result?
- Downstream result: the top-20 hub genes, by degree, of the BH-called pseudobulk network.
- Compare the PB network with the PB network after removing the PB-only pairs. Use the Jaccard index J of the two top-20 hub sets.
- Control: remove the same number of PB edges at random (200 draws) and compute J each time.

Q2 = **YES** if both hold:
- J < 0.5;
- J is below the 5th percentile of the random-removal distribution.

Ties in degree are broken by gene index.

## Reporting
- All outcomes are reported, whichever way they go: Q1, Q2 and the calibration gate, for S1 (primary) and S2 (secondary).
- Nothing is re-run with different choices. Any deviation is listed explicitly.
