# CRITERIA: every frozen decision criterion and its outcome

Built from `criteria_inputs.tgz` (prereg/, rev2/, rev3/, rev4/, results/, paper2/). Every number below is read from a file in that bundle; the file and key are given in the Source column. Where a result file stores a decision flag, the verdict was recomputed from the underlying numbers and asserted equal to the flag. Plans are in the order FIX, XS, POP, DSTD, XU, COUPLE, DC, REV, ONEK, REV2, REV3, REV4.

## Legend

- **PASS / FAIL**: the frozen rule was met / not met. For rules evaluated per subtype, cohort size or method cell, PASS means every cell met it.
- **MIXED**: the rule was evaluated in several cells and some met it, some did not (counts in Evidence). This label is an addition to the PASS / FAIL / NOT EVALUATED / DESCRIPTIVE set, used only where a single PASS or FAIL would hide the split.
- **DESCRIPTIVE**: the plan reports the quantity without a frozen threshold (no criterion). **NOT EVALUATED**: planned, no results in the bundle. **NOT DETERMINABLE FROM BUNDLE**: the outcome cannot be read from the supplied files.
- **(partial)**: based on partial runs (REV4 R1: 13 of 20 runs in P1, 12 in P2, none in P3; REV4 R2c: 36 of 40 runs in P1). For R1, FAIL (partial) cannot change with more runs (the bound is at most 2 runs for n = 12, 13 and 20, and every failing cell already has at least 3 runs with a call); PASS and MIXED (partial) are provisional.
- **Run-level rule**: a method holds in a cell if the fraction of runs with at least one BH call (q = 0.05) is at most 0.05 + 2 SE, SE binomial at p = 0.05 (REV3 and REV4 plans state this; it gives at most 9 of 100, 4 of 40, 2 of 20, 1 of 10, 4 of 33 runs).
- **MCSE**: Monte Carlo standard error as stored in the result files; **top-1%**: type-I over the 1% of pairs most correlated within pseudo-donors; **SE** in top-1% bounds is the stored binomial SE (pairs x reassignments).

## Summary

75 criteria rows. DESCRIPTIVE: 21; FAIL: 19; PASS: 13; MIXED: 12; MIXED (partial): 4; NOT DETERMINABLE FROM BUNDLE: 2; PASS (partial): 2; FAIL (partial): 1; NOT EVALUATED: 1.

| Plan | Rows | PASS | FAIL | MIXED | DESCRIPTIVE | NOT EVALUATED | NOT DETERMINABLE FROM BUNDLE | PASS (partial) | FAIL (partial) | MIXED (partial) |
|---|---|---|---|---|---|---|---|---|---|---|
| FIX | 4 | 1 | 3 |  |  |  |  |  |  |  |
| XS | 2 |  | 2 |  |  |  |  |  |  |  |
| POP | 4 | 2 | 1 |  | 1 |  |  |  |  |  |
| DSTD | 3 |  | 3 |  |  |  |  |  |  |  |
| XU | 4 | 1 | 2 |  | 1 |  |  |  |  |  |
| COUPLE | 1 |  |  |  | 1 |  |  |  |  |  |
| DC | 3 | 2 |  |  | 1 |  |  |  |  |  |
| REV | 9 | 2 | 1 | 1 | 5 |  |  |  |  |  |
| ONEK | 4 | 1 |  | 2 |  |  | 1 |  |  |  |
| REV2 | 13 | 1 | 3 | 1 | 7 |  | 1 |  |  |  |
| REV3 | 15 | 2 | 4 | 7 | 2 |  |  |  |  |  |
| REV4 | 13 | 1 |  | 1 | 3 | 1 |  | 2 | 1 | 4 |

## FIX_PREREG (real-cell pseudo-donor null in SEA-AD; LCP2)

- Files: `prereg/FIX_PREREG.md`, `prereg/FIX_RESULTS.md`
- Hash: prereg/FIX_PREREG.sha256 matches
- Scope in bundle: Item 1: K = 5 reassignments x 2 supertypes; Item 2: R = 100 x 6 designs.

| ID | Criterion | Outcome | Evidence | Source (file: key) |
|---|---|---|---|---|
| FIX-1a | Leak produces spurious pseudobulk discoveries in real cells: PB mean calls ≥50 and ≥10× LCP’s (LCP counted as ≥1). | **PASS** | PB mean false calls 3,713.8 (S1), 2,713.8 (S2); LCP 0.0, 0.0. Decision YES in both. | `results/fix/fix_results.json: null.{S1,S2}.PB.mean_calls, .LCP.mean_calls, .decision_leak_real` |
| FIX-1b | LCP holds the real-cell null: all-pair type-I ≤0.055 and top-1% within-correlated type-I ≤0.05+3 SE. | **FAIL** | All-pair type-I 0.081 / 0.061 (S1/S2; not met); top-1% 0.092 / 0.099 vs bound 0.0599 (SE 0.0033): not met. 0 BH calls in all 5 reassignments. | `results/fix/fix_results.json: null.{S1,S2}.LCP.typeI, .typeI_top1, .typeI_top1_se, .LCP_holds` |
| FIX-1c | LCP2 holds the real-cell null: all-pair type-I ≤0.055 and top-1% within-correlated type-I ≤0.05+3 SE. | **FAIL** | All-pair type-I 0.046 / 0.050 (S1/S2; met); top-1% 0.068 / 0.085 vs bound 0.0599 (SE 0.0033): not met. 0 BH calls in all 5 reassignments. | `results/fix/fix_results.json: null.{S1,S2}.LCP2.typeI, .typeI_top1, .typeI_top1_se, .LCP2_holds` |
| FIX-2 | Adopt LCP2 if (a) FDR ≤0.05+2 MCSE in F_D24, F_D40, F_S2, F_D16_n100; (b) also in F_D16 or F_S1; (c) power ≥0.9× LCP (F_D24, F_D40). | **FAIL** | (a) met: LCP2 FDR 0.033, 0.041, 0.042, 0.039. (b) met (F_D16 0.052, F_S1 0.025). (c) not met: power ratio 0.87 at F_D24, 0.96 at F_D40. Not adopted. | `results/fix/fix_results.json: sim.<design>.{LCP,LCP2}.{FDR,FDR_mcse,power}; LCP2_decision.{a,b,c,adopt}` |

Notes and deviations:

- No deviation from the plan (FIX_RESULTS.md, Deviations). The pooling caveat (null_caveat.json) is a descriptive addition.

## XS_PREREG (LCPX, cross-split LCP) with Amendment 1

- Files: `prereg/XS_PREREG.md`, `prereg/POP_RESULTS.md (XS section)`
- Hash: prereg/XS_PREREG.sha256 matches the amended file; XS_PREREG_original.sha256 hashes the pre-amendment file
- Scope in bundle: R = 100 x 6 designs; SEA-AD null K = 5 x 2.

| ID | Criterion | Outcome | Evidence | Source (file: key) |
|---|---|---|---|---|
| XS-1 | Adopt LCPX if (a) FDR ≤0.05+2 MCSE in all six designs; (b) holds SEA-AD null in S1, S2; (c) power ≥0.8× LCP in F_D24, F_D40. | **FAIL** | (a) not met: LCPX FDR 0.191 in F_S1 (bound 0.078) and 0.067 in F_S2 (0.055). (b) not met: type-I 0.110/0.066, top-1% 0.153/0.074. (c) not met: power ratio 0.16 (F_D24), 0.38 (F_D40). | `results/xsplit/xs_results.json: sim.<design>.{LCPX,LCP}; null.{S1,S2}.LCPX; decision.{a,b,c,adopt}` |
| XS-2 | Claim LCPX-S tests Σ_B=0 only if weak-null type-I and FDR ≤0.05+2 MCSE in F_D24, F_D40, F_S1, F_S2 (weak null = QQ_cross). | **FAIL** | LCPX-S F_S1: FDR 0.106 (bound 0.081), weak-null 0.070 (0.062); F_S2 FDR 0.0557 (bound 0.0556). F_D24, F_D40 within bounds. Claim not made. | `results/xsplit/xs_results.json: sim.<design>.LCPXS.{FDR,FDR_mcse,WN,WN_mcse}; decision.weak_null_claim` |

Notes and deviations:

- Amendment 1 (before any full run; only a B = 2000 smoke test had been run) restricted the weak null to QQ_cross; the table uses the amended definition. F_S2 LCPX-S FDR exceeds its bound by 0.00005 (0.0557 vs 0.0556).

## POP_PREREG (PsychAD population-scale null and real donors)

- Files: `prereg/POP_PREREG.md`, `prereg/POP_RESULTS.md`
- Hash: prereg/POP_PREREG.sha256 matches (path pop/POP_PREREG.md)
- Scope in bundle: Part A: K = 3 per subtype at D = 25, 50, 100 (K = 1 only), 200, 400, 800, all.

| ID | Criterion | Outcome | Evidence | Source (file: key) |
|---|---|---|---|---|
| POP-1 | PB false discoveries grow with cohort size: PB mean calls at D=all ≥3× those at D=100 and ≥100 (per subtype). | **PASS** | YES in P1, P2, P3. Calls at D=all 19,158/31,617/24,544; ratio to D=100 4.9/5.0/6.2 (same reassignment; with K=3 means 4.9/5.0/6.1). D=100 point has K=1. | `results/pop/pop_results.json: A.<P>.all.PB.{calls,mean_calls}, A.<P>.100.PB.calls; prereg/POP_RESULTS.md` |
| POP-2 | LCPX holds the real-cell null at population scale: at every D, all-pair type-I ≤0.055, top-1% ≤0.05+3 SE, mean calls ≤1. | **FAIL** | NO in all three subtypes; criterion fails at every D. LCPX mean false calls at D=all 25,901/32,801/15,538; all-pair type-I 0.416/0.497/0.354; already 0.064/0.076/0.073 at D=25. The stated prediction (LCPX near 0 at every D) was not borne out. | `results/pop/pop_results.json: A.<P>.<D>.LCPX.{mean_calls,typeI,typeI_top1,se_top1,holds}` |
| POP-3 | Replication does not protect pseudobulk: PB replicated false calls ≥50 in the Part A null at D=half (per subtype). | **PASS** | YES in all three. Mean replicated false calls PB 3,876/12,515/7,876 (minimum single run 3,798); LCP 594/2,130/1,212; XDS 0 in all 9. | `results/pop/pop_results.json: rep_null.<P>.{PB,LCP,LCPXd}.{replicated_false,mean}` |
| POP-4 | Part B, real donors at full D: calls, overlaps, replication in donor halves. “Results are reported whatever they are” (no threshold). | **DESCRIPTIVE** | PB calls 53,263/57,280/37,325; XDS 49,852/55,021/32,790; PB-only 13,307/10,086/14,643. Half-to-half replication XDS 79/82/68%, PB 73/80/62%. | `results/pop/pop_results.json: B.<P>.{PB,LCPXd,PB_only,rep_PB,rep_LCPXd}` |

Notes and deviations:

- POP_RESULTS.md lists the real-donor XDS call counts as 32,790 / 49,852 / 55,021, which is not P1/P2/P3 order; the table uses pop_results.json (P1 49,852; P2 55,021; P3 32,790). Deviation (stated in appendix): the D<=100 points were run after the larger ones and with one reassignment per subtype (K = 1, not 3). P3 uses reassignments k = 1..3 (k = 0 was the exploratory reassignment excluded in DSTD_PREREG). POP-1 ratios use the same reassignment (k = 0 for P1, P2; k = 1 for P3), as in POP_RESULTS.md.

## DSTD_PREREG (LCPd, LCPXd = XDS)

- Files: `prereg/DSTD_PREREG.md`, `prereg/POP_RESULTS.md (DSTD section)`
- Hash: prereg/DSTD_PREREG.sha256 matches
- Scope in bundle: E1: R = 100 x 6 designs; E2: K = 5 x 2; E3: K = 3 x 3 (P1, P2 k = 0-2; P3 k = 1-3).

| ID | Criterion | Outcome | Evidence | Source (file: key) |
|---|---|---|---|---|
| DSTD-1 | Adopt LCPd over LCP if (a) FDR bound in four designs; (b) holds SEA-AD null; (c) power ≥0.9× LCP; (d) holds PsychAD null at D=all. | **FAIL** | (a) met (FDR 0.037–0.062). (b) not met: top-1% type-I 0.197/0.141. (c) not met: power ratio 0.81 (D24), 0.88 (D40). (d) not met: PsychAD mean false calls 57.7/38.3/1,218; type-I 0.058/0.059/0.118. | `results/dstd/d_results.json: sim, null, decision_LCPd.{a,b,c}; (d) computed from results/pop/pop_results.json A.<P>.all.LCPd` |
| DSTD-2 | Same adoption rule (a)–(d) for LCPXd (= XDS) as secondary method. | **FAIL** | (a) met (FDR 0.025–0.049, all six designs). (b) met: SEA-AD type-I 0.050/0.050, 0 calls. (c) not met: power ratio 0.14 (D24), 0.31 (D40). (d) met: mean calls 0.0/0.67/0.0, type-I 0.050/0.050/0.048. XDS does not replace LCP by the frozen rule. | `results/dstd/d_results.json: sim.<design>.LCPXd, null.<S>.LCPXd, decision_LCPXd.{a,b,c}; (d) from results/pop/pop_results.json A.<P>.all.LCPXd` |
| DSTD-3 | Claim “design-induced bias is a distinct artifact” if LCP all-pair type-I &gt;0.055 in ≥2 subtypes at D=all while LCPd’s does not. | **FAIL** | LCP type-I 0.201/0.301/0.265 (&gt;0.055 in 3 of 3) but LCPd 0.058/0.059/0.118 also &gt;0.055 in 3 of 3. Claim not supported as frozen. | `results/pop/pop_results.json: A.<P>.all.{LCP,LCPd}.typeI` |

Notes and deviations:

- Written after one exploratory reassignment (PsychAD P3, k = 0), excluded from every LCPd/LCPXd evaluation; interim XS simulation results (~45/100 replicates) were also seen before the plan. Criterion (d) is not stored as a flag in d_results.json; it was recomputed here from pop_results.json (same K = 3 reassignments).

## XU_PREREG (exact within-donor U-statistic, studentized across donors)

- Files: `prereg/XU_PREREG.md`
- Hash: prereg/XU_PREREG.sha256 matches
- Scope in bundle: E1: K = 3 x 3 subtypes, D in {100, 200, 400, all}; E2: K = 5 x 2; E3: R = 50 x 3 designs; E4: R = 100 x 6 designs.

| ID | Criterion | Outcome | Evidence | Source (file: key) |
|---|---|---|---|---|
| XU-1 | XU holds the real-cell null at population scale: at every D≥200, every subtype, mean calls ≤1, all-pair type-I ≤0.055, top-1% ≤0.05+3 SE. | **FAIL** | Fails in all 9 subtype×D cells. XU mean false calls at D=200: 235.7/73.3/391.0; at D=all 19,934/14,258/12,073. All-pair type-I at D=200 0.125/0.106/0.127. | `results/xu/xu_results.json: E1.<P>.<D>.XU.{mean_calls,typeI,typeI_top1,holds}; decision.holds_pop` |
| XU-2 | XU is the recommended large-cohort test: XU-1 AND E3 FDR ≤0.05+2 MCSE in L200, L500, L1000 AND power ≥1.5× XDS in L500, L1000. | **FAIL** | XU-1 not met. E3 FDR within bound (0.042/0.045/0.048). Power ratio XU/XDS 1.01 (L500), 1.00 (L1000), below 1.5 (XDS power 0.99, 1.00). XDS remains the large-cohort test, as the plan prescribes. | `results/xu/xu_results.json: E3.<L>.{XU,XDS}.{FDR,FDR_mcse,power}; decision.recommended` |
| XU-3 | Weak-null validity: XU weak-null type-I ≤0.05+2 MCSE in L200, L500, L1000. | **PASS** | XU weak-null type-I 0.052/0.045/0.029 vs bounds 0.072/0.067/0.060. | `results/xu/xu_results.json: E3.<L>.XU.{WN,WN_mcse}; decision.weak_null` |
| XU-4 | E2 (SEA-AD null, XU with t_{D-1}) and E4 (six small designs): “reported, no decision.” | **DESCRIPTIVE** | E2 XU mean false calls 0.2 (S1), 61.4 (S2); all-pair type-I 0.079/0.109; top-1% 0.102/0.119. E4 results are in the bundle (E4.&lt;design&gt;.XU), not tabulated here. | `results/xu/xu_results.json: E2.{S1,S2}; E4` |

Notes and deviations:

- Designed after all XS/DSTD/POP results were seen; evaluated on fresh seeds only. XU-3 passes but does not make XU recommended.

## COUPLE_PREREG (depth-state coupling stress test)

- Files: `prereg/COUPLE_PREREG.md`
- Hash: prereg/COUPLE_PREREG.sha256 matches
- Scope in bundle: R = 50 x 3 designs (LS200, LS500, LS1000).

| ID | Criterion | Outcome | Evidence | Source (file: key) |
|---|---|---|---|---|
| COUPLE-1 | No adoption decision. Prediction: PB and LCP FDR far above 0.05; XDS FDR inflated less than LCP and XU. | **DESCRIPTIVE** | Prediction consistent with results. FDR at LS200/LS500/LS1000: PB 0.524/0.605/0.718; LCP 0.104/0.143/0.191; XU 0.060/0.104/0.223; XDS 0.059/0.047/0.056 (below LCP in all three, below XU at LS500, LS1000, equal at LS200: 0.0594 vs 0.0596). | `results/xu/couple_results.json: <LS>.<method>.FDR` |

## DC_PREREG (AD differential co-expression)

- Files: `prereg/DC_PREREG.md`
- Hash: prereg/DC_PREREG.sha256 matches
- Scope in bundle: K = 3 null reassignments x 3 subtypes (k = 10, 11, 12).

| ID | Criterion | Outcome | Evidence | Source (file: key) |
|---|---|---|---|---|
| DC-1 | Pseudobulk DC is confounded by cell number: PB-DC mean false calls ≥100 in at least two subtypes. | **PASS** | YES in 3 of 3: mean false calls 5,188/3,134/1,002 (range 718–5,634); per-pair type-I 0.13–0.21. | `results/dc/dc_summary.json: <P>.PB.null_calls, .null_typeI` |
| DC-2 | XDS-DC holds: in every subtype mean false calls ≤1 and all-pair type-I ≤0.055. | **PASS** | YES. XDS-DC false calls 0 in all 9 runs; type-I 0.048–0.052 (max 0.0518). LCP-DC (not decided): 0–6 calls, type-I 0.048–0.056. | `results/dc/dc_summary.json: <P>.XDS.null_calls, .null_typeI; <P>.LCP` |
| DC-3 | Real-donor results reported descriptively (number of DC calls per method, overlap). | **DESCRIPTIVE** | Real-donor DC calls P1/P2/P3: PB 5,806/2,441/2,561; LCP 2,319/485/1; XDS 816/19/0. Overlap between methods not tabulated in the bundle. | `results/dc/dc_summary.json: <P>.{PB,LCP,XDS}.real_calls` |

Notes and deviations:

- The plan states only that real-donor results are reported descriptively; group sizes n1/n0 per subtype are in dc_summary.json.

## REV_PREREG (internal review round 1; R1–R5, amendment R6)

- Files: `prereg/REV_PREREG.md`
- Hash: prereg/REV_PREREG.sha256 matches
- Scope in bundle: R1: 100 runs x 3 subtypes; R2: SEA-AD k = 10..14; R3: K = 3; R6: K = 3 x 3.

| ID | Criterion | Outcome | Evidence | Source (file: key) |
|---|---|---|---|---|
| REV-R1a | XDS controls FDR at population scale (per subtype): fraction of 100 runs with ≥1 BH call ≤0.05+2 SE (0.0936; ≤9 runs). | **FAIL** | Fails in 3 of 3: runs with a call 11, 20, 11 of 100 (P1, P2, P3). Context: PB 100, 100, 100 of 100 (mean calls 18,703/32,128/23,407). | `results/rev/out/R1_<P>_k0..k99.json: XDS.V, PB.V (count of runs with V>=1)` |
| REV-R1b | Same decision for XDS-M (≤9 of 100 runs with a call, per subtype). | **PASS** | Holds in 3 of 3: runs with a call 4, 5, 2 of 100 (P1, P2, P3). | `results/rev/out/R1_<P>_k*.json: XDSM.V` |
| REV-R1c | Same decision for KX (≤9 of 100 runs with a call, per subtype). | **MIXED** | Holds in P1, P2 (4, 8 of 100); fails in P3 (36 of 100; max 25 calls in one run). | `results/rev/out/R1_<P>_k*.json: KX.V` |
| REV-R1d | Per run, counts of p&lt;1e-3, 1e-4, 1e-5 recorded (expected 89, 8.9, 0.89 for ~89,000 pairs); no threshold. | **DESCRIPTIVE** | XDS-M mean counts per run: p&lt;1e-3 85.0–88.8, p&lt;1e-4 8.3–9.1, p&lt;1e-5 0.79–0.89 (range over subtypes; pairs = 89,676). XDS: 89.0–89.5, 8.9–9.8, 0.86–1.31. | `results/rev/out/R1_<P>_k*.json: <method>.{n3,n4,n5}, P` |
| REV-R2 | R2 defines XDS-M (S=5) and KX; evaluated in L/LS and six small designs and the SEA-AD null (k=10..14). No decision rule. | **DESCRIPTIVE** | XDS-M FDR 0.043–0.054 in the six large designs; 0.135–0.265 in F_D16, F_D24, F_D16_n100. KX FDR 0.268–0.897 in F_D16, F_D24, F_D40, F_D16_n100. SEA-AD top-1% type-I (S1/S2): XDS-M 0.070/0.043, KX 0.055/0.061. | `results/rev/rev_sim_results.json: <design>.{XDSM,KX}.FDR; paper2/seanull_summary.json: XDSM, KX` |
| REV-R3 | XDS-DCg holds “as in DC_PREREG”: per subtype mean false calls ≤1 and all-pair type-I ≤0.055. | **PASS** | YES. False calls in the 3 runs per subtype: P1 0, 0, 0; P2 0, 0, 1; P3 0, 0, 0; type-I 0.048–0.052. Real donors: 117/0/0 calls; overlap with pooled-centre XDS-DC calls not determinable from bundle. | `paper2/dcg_summary.json: <P>.null_calls, .null_typeI, .real_calls; results/rev/rev_dc_r5.json: <P>.XDSg_real, .XDSg_in_pooled` |
| REV-R4 | Positive control (descriptive): XDS on real donors, XIST–Y-linked pairs; “Expected: negative between-donor association (sex), detected by XDS.” | **DESCRIPTIVE** | Expectation met: 21 of 21 XIST–Y pairs negative with p ≤ 5e-40 (3 subtypes × 7 pairs). Y–Y pairs (reported in the paper, not in this plan): 63 of 63 positive, max p 3e-08. | `results/rev/out/real_<P>.json: sexpairs[].{pair,p_XDS,sign_XDS}` |
| REV-R5 | Size of pseudobulk false edges (descriptive): absolute r of PB false calls in the POP null at full D; median and counts above 0.1, 0.2. | **DESCRIPTIVE** | Median absolute r 0.10–0.13; share &gt;0.1: 55–84%; share &gt;0.2: 4–14% (P1/P2/P3, n = 57,475/94,852/73,631 false calls). | `results/rev/rev_dc_r5.json: <P>.false_r_median, .false_r_gt01, .false_r_gt02, .false_r_n` |
| REV-R6 | Amendment R6: donor-centred null to bound pooling inflation (K=3, full D); BH calls, type-I, median within-donor r vs real. Descriptive; no decision. | **DESCRIPTIVE** | PB mean calls 8,708–22,872 per subtype; XDS 0–3 per run (type-I 0.049–0.052). Pooling raised the median absolute within-pseudo-donor r by 37–46% over real donors. | `results/rev/out/R6_<P>_k0..k2.json; results/rev/pooling.json: <P>.ratio_median` |

Notes and deviations:

- R1 decision rule: SE is binomial at p = 0.05 (confirmed by the later plans, “SE at p = 0.05”; the rule gives ≤9 of 100). Under the alternative reading (SE at the observed fraction) XDS in P1 and P3 (11 of 100) would be borderline; the paper’s reading is the one used. REV_PREREG states no explicit decision for PB, for R2 or for R4–R6.

## ONEK_PREREG (OneK1K blood) with Amendment 1

- Files: `prereg/ONEK_PREREG.md`
- Hash: prereg/ONEK_PREREG.sha256 matches (hashes the file including the amendment)
- Scope in bundle: K = 3 reassignments x 3 cell types, D = 200 and all.

| ID | Criterion | Outcome | Evidence | Source (file: key) |
|---|---|---|---|---|
| ONEK-1 | PB false discoveries grow with cohort size (per cell type): PB mean calls at D=all ≥2× those at D=200 and ≥100. | **MIXED** | YES in 2 of 3. Ratio all/200: 2.25 (transitional B; 24,641 calls at D=all), 1.68 (naive B; not met), 8.51 (naive CD4 T). Median cells per donor 20/56/244. | `paper2/onek_summary.json: <B>.all.PB.mean_calls, <B>.200.PB.mean_calls; results/onek/onek_<B>_info.json` |
| ONEK-2 | XDS holds (per cell type) at both D: mean calls ≤1, all-pair type-I ≤0.055, top-1% type-I ≤0.05+3 SE. | **PASS** | Holds in 3 of 3. XDS mean calls ≤1.0; all-pair type-I 0.049–0.050; top-1% ≤0.052 (bound 0.063). | `paper2/onek_summary.json: <B>.{all,200}.LCPXd.{mean_calls,typeI,top1,holds}` |
| ONEK-3 | Amendment 1: XDS-M holds (per cell type), same rule at both D. | **MIXED** | Holds in 2 of 3 (naive B, naive CD4 T). Fails in transitional B at D=200: top-1% type-I 0.079 &gt; 0.063 (all-pair 0.053, 0 calls). Run-level frequencies over K=3 (too few for a decision): 0 of 18 XDS-M runs had a call across all cells. | `paper2/onek_summary.json: <B>.{all,200}.XDSM.{typeI,top1,holds,calls}` |
| ONEK-4 | Real donors: XDS and PB calls, PB-only count (descriptive). | **NOT DETERMINABLE FROM BUNDLE** | No OneK1K real-donor results (XDS/PB calls or PB-only counts) are in the bundle. | `(none in bundle)` |

Notes and deviations:

- Top-1% SE for ONEK-2/3 is not stored in the summary; recomputed as binomial over 798 top pairs × K = 3 reassignments (SE 0.0045, bound 0.0634); the recomputed verdicts agree with the stored holds flags in all 12 cells.

## REV2_PREREG (internal review round 2) with Amendment 1

- Files: `rev2/REV2_PREREG.md`, `rev2/REV2_AMENDMENT1.md`
- Hash: no .sha256 file in bundle; sha256 of bundle files: PREREG 1a7e4de14a1f…, AMENDMENT1 4f42dbb837b0…
- Scope in bundle: A: 100 runs x 3 subtypes x D in {200, 400, all}; E: 30 runs x 3 cell types.

| ID | Criterion | Outcome | Evidence | Source (file: key) |
|---|---|---|---|---|
| REV2-A1 | Part A, per method/subtype/D (200, 400, all): “holds” if fraction of 100 runs with V≥1 ≤0.05+2 SE (≤9 runs). XDS-M. | **MIXED** | Holds in 8 of 9 cells; fails P3 at D=400 (10 of 100 runs). Runs with a call across cells: 0–10 of 100. | `rev2/interim/REV2_SUMMARY_latest.json: A.<P>.<D>.XDSM.{any,n}` |
| REV2-A2 | Same rule, KX. | **PASS** | Holds in 9 of 9 cells; runs with a call 0–9 of 100. | `rev2/interim/REV2_SUMMARY_latest.json: A.<P>.<D>.KX.{any,n}` |
| REV2-A3 | Same rule, XDS (single split, normal-moment reference). | **FAIL** | Holds in 0 of 9 cells; runs with a call 11–42 of 100. | `rev2/interim/REV2_SUMMARY_latest.json: A.<P>.<D>.XDS.{any,n}` |
| REV2-A4 | Same rule, PB and LCP. | **FAIL** | PB and LCP hold in 0 of 9 cells (100 of 100 runs with a call in every cell). Mean BH calls PB 2,513–23,779; LCP 114–8,346. | `rev2/interim/REV2_SUMMARY_latest.json: A.<P>.<D>.{PB,LCP}.{any,meanV}` |
| REV2-A5 | Primary question: does LCPd hold where XDS-M holds? If LCPd holds at every D and subtype it is reported as co-recommended. | **FAIL** | LCPd holds in 3 of 9 cells (all in P1: 8, 8, 6 of 100 runs); fails P2 (21, 21, 42) and P3 (17, 13, 20). Of XDS-M’s 8 holding cells LCPd holds in 3. Not co-recommended. | `rev2/interim/REV2_SUMMARY_latest.json: A.<P>.<D>.LCPd.{any,n}` |
| REV2-B | Part B, descriptive: mean over genes of squared design-standardized half-donor values (R^A) by decile of n_d (A5 diagnostic). | **DESCRIPTIVE** | Mean (R^A)^2 by decile 0.89–1.00 across subtypes (P1 0.916–0.951, P2 0.895–0.996, P3 0.896–0.959). | `rev2/interim/REV2_SUMMARY_latest.json: C.<P>.A5.{decile_median_n,mean_R2}` |
| REV2-C | Part C, descriptive: real donors with XDS-M (calls, split-half replication, three split seeds with Jaccard, XIST–Y and Y–Y pairs). | **DESCRIPTIVE** | XDS-M calls 56,075/60,272/37,866 (63/68/43% of pairs); Jaccard between split seeds 0.86–0.95; half-to-half replication 74–85%. Sex pairs: XDS-M XIST–Y 21 of 21, Y–Y 63 of 63 significant with the expected sign. | `rev2/interim/REV2_SUMMARY_latest.json: C.<P>.{seeds,halves,sex_XDSM}` |
| REV2-D1 | Part D1, donor-centred DC null (k=10–12): PB-DC, XDS-DCg, XDS-Mg, plus XDS-Mgu (Amendment 1: “same decision rules as XDS-Mg”). No threshold stated. | **DESCRIPTIVE** | PB-DC false calls 995–5,753 per run (9 runs). XDS-DCg 0 calls in all 9 runs; XDS-Mg 1 call in total (max 1); XDS-Mgu 0 in all 9; type-I 0.048–0.053. Real donors (P1/P2/P3): PB 5,806/2,441/2,561; XDS-DCg 117/0/0; XDS-Mg 465/0/0; XDS-Mgu 351/341/0. | `rev2/interim/REV2_SUMMARY_latest.json: D1.<P>.{k10,k11,k12,real}.<method>.{V,t05}` |
| REV2-D2 | Part D2, true-DC simulation (R=50; 101 true-DC pairs of 3,160): FDR and power of PB-DC, XDS-DCg, XDS-Mg (+XDS-Mgu). Reported; no threshold. | **DESCRIPTIVE** | FDR/power: PB-DC 0.485/0.962; XDS-DCg 0.077/0.094; XDS-Mg 0.071/0.264; XDS-Mgu 0.031/0.280 (n_alt 101, n_pairs 3,160, R = 50). | `rev2/interim/REV2_SUMMARY_latest.json: D2.<method>.{FDR,power}, n_alt, n_pairs` |
| REV2-D3 | Part D3, descriptive: split-half replication of real-donor XDS-DCg, XDS-Mg and XDS-Mgu calls (P1). | **DESCRIPTIVE** | Replicated / calls in half 1 (P1): PB-DC 438/3,663; XDS-DCg 0/0; XDS-Mg 2/4; XDS-Mgu 37/50. | `rev2/interim/REV2_SUMMARY_latest.json: D3.P1.<method>.{calls_h0,replicated}` |
| REV2-E | Part E: 30 further OneK1K reassignments per cell type; run-level fraction for XDS-M over all 33 runs; PBlin vs PB (descriptive). No threshold stated. | **DESCRIPTIVE** | XDS-M runs with a call: 4/33, 1/33, 1/33 (transitional B, naive B, naive CD4 T); all within the Part-A rule for n=33 (≤4 of 33) — the paper’s “within the bound”. PB 30/30 in all; PBlin mean calls 19,700/13,894/16,827 vs PB 22,889/22,948/20,035. | `rev2/interim/REV2_SUMMARY_latest.json: E.<B>.{PB,PBlin,LCP,XDSM}.{any,n,meanV}; paper2/onek_summary.json: <B>.all.XDSM.calls (3 earlier runs)` |
| REV2-F | Part F: A5 overdispersion stress simulation (L200, L500; R=50): FDR and power. “No decision.” | **DESCRIPTIVE** | FDR at L200/L500: XDS-M 0.043/0.047; KX 0.038/0.038; LCPd 0.067/0.071; XDS 0.088/0.084; LCP 0.117/0.220; PB 0.540/0.640. Power XDS-M 0.83/0.99, KX 0.17/0.67. | `rev2/interim/REV2_SUMMARY_latest.json: F.<L>.<method>.{FDR,power}` |
| REV2-G | Part G: oracle (matched-FDP) power from existing simulation p-values, per replicate (descriptive). | **NOT DETERMINABLE FROM BUNDLE** | No matched-FDP / oracle-power result is in the bundle (the “oracle” hits in results/dev and results/heldout refer to the DEV/HELDOUT ORACLE test, not this analysis). | `(none in bundle)` |

Notes and deviations:

- Spearman correlation of decile rank with mean (R^A)^2, computed here (n = 10 deciles, not in the bundle): −0.70/−0.65/−0.22 for P1/P2/P3 — a small drift in magnitude (range 0.895–0.996); the paper’s statement of no trend refers to magnitude.
- The only REV2 summary in the bundle is rev2/interim/REV2_SUMMARY_latest.json (file name says “interim”); it holds complete counts for every part (A: 100 runs in each of 9 cells; E: 30 per cell type; D1 k = 10–12; D2/F R = 50) and reproduces the rev2nums.tex macros, but the bundle does not show whether a later “final” file differs.
- REV2_PREREG states no decision rule for D1, D2, D3, E or F; Amendment 1 says XDS-Mgu follows “the same decision rules as those frozen for XDS-Mg”, but no rule for XDS-Mg is written in REV2_PREREG (REV_PREREG R3 refers to the DC_PREREG rule: mean false calls ≤1, type-I ≤0.055, which the D1 null runs would meet for XDS-DCg, XDS-Mg and XDS-Mgu). Part E: the paper judges XDS-M “within the bound” over 33 runs using the Part-A rule (0.05+2 SE, ≤4 of 33); REV2_PREREG does not state this threshold for E.

## REV3_PREREG (internal review round 3) with Amendment 1

- Files: `rev3/REV3_PREREG.md`, `rev3/REV3_AMENDMENT1.md`
- Hash: no .sha256 file in bundle; sha256 of bundle files: PREREG 314e01284fb9…, AMENDMENT1 14b1eb304ef7…
- Scope in bundle: A2: 20 runs per null (B1, B2, P3); A2x: 10; A3: 100 x 3 subtypes; B2: 20 x 3; B3: 40.

| ID | Criterion | Outcome | Evidence | Source (file: key) |
|---|---|---|---|---|
| REV3-A2v | A2 validity: a run is invalid if the largest GPU–CPU difference in log10 p exceeds 0.05 for p&gt;1e-8 (2,000 random pairs). | **PASS** | 0 invalid runs in all 8 data set×null blocks (140 runs); largest check discrepancy 4.7e-05. | `rev3/results/A2_SUMMARY.json: <block>.{runs,invalid,maxcheck}` |
| REV3-A2a | A2 tail calibration rule (obs/exp ≤1.5 at p≤1e-4 and 1e-5; run-level ≤0.147), PB. | **FAIL** | PB tail-calibrated in 0 of 36 combinations (B1/B2 genes &gt;0.1 UMI/cell 0/24; B2 0.03–0.1 UMI/cell (50 genes; 1e-5 not assessable) 0/4; P3 0/8). | `rev3/results/A2_SUMMARY.json: <block>.<D>.PB.<stratum>.{calibrated,ratio,run_ok}` |
| REV3-A2b | Same rule, XDS (single split 0, Hoeffding-moment normal reference). | **MIXED** | B1/B2 genes &gt;0.1 UMI/cell 16/24; B2 0.03–0.1 UMI/cell (50 genes; 1e-5 not assessable) 3/4; P3 0/8. | `rev3/results/A2_SUMMARY.json: <block>.<D>.XDS.<stratum>.calibrated` |
| REV3-A2c | Same rule, XDS-M with S=5 splits (the validated-range method). | **MIXED** | B1/B2 genes &gt;0.1 UMI/cell 23/24; B2 0.03–0.1 UMI/cell (50 genes; 1e-5 not assessable) 4/4; P3 8/8. One failure: B1 donor-centred null, D=200, &gt;1 UMI/cell (3 of 20 runs with a call; bound 2). | `rev3/results/A2_SUMMARY.json: <block>.<D>.XDSM5.<stratum>.{calibrated,run_ok,anyV,runs,ratio,exp}` |
| REV3-A2d | Same rule, XDS-M with S=1 and S=10 splits. | **MIXED** | S=1: B1/B2 genes &gt;0.1 UMI/cell 23/24; B2 0.03–0.1 UMI/cell (50 genes; 1e-5 not assessable) 3/4; P3 8/8. S=10: B1/B2 genes &gt;0.1 UMI/cell 22/24; B2 0.03–0.1 UMI/cell (50 genes; 1e-5 not assessable) 4/4; P3 8/8. | `rev3/results/A2_SUMMARY.json: <block>.<D>.{XDSM1,XDSM10}.<stratum>.calibrated` |
| REV3-A2r | Validated gene range for XDS-M (S=5): strata tail-calibrated in every data set, null and D. | **DESCRIPTIVE** | Derived here from the flags (UMI/cell): calibrated in every combination at 0.3–1 (12/12) and 0.1–0.3 (12/12); not at &gt;1 (11/12) or 0.03–0.1 (5/8). Strict range: 0.1–1 UMI/cell; the appendix gives counts, not the range. | `rev3/results/A2_SUMMARY.json: <block>.<D>.XDSM5.<stratum>.calibrated (stratum 0..3 = >1, 0.3–1, 0.1–0.3, 0.03–0.1 UMI/cell)` |
| REV3-A2x | Amendment 1: same decision rule on OneK1K B1 genes ranked 2,001–3,000, n=10 runs (bound 0.188, ≤1 of 10 runs). | **MIXED** | Genes 0.1–0.3 UMI/cell (640 genes): XDS-M S=5 4/4, S=1 4/4, S=10 4/4, XDS 0/4, PB 0/4. Genes 0.03–0.1 UMI/cell (360 genes): XDS-M S=5 1/4, S=1 4/4, S=10 3/4, XDS 0/4, PB 0/4. (Combinations = null × D.) | `rev3/results/A2_SUMMARY.json: B1x_n0, B1x_n1 .<D>.<method>.<stratum>.calibrated` |
| REV3-A3a | A3 run-level rule (≤9 of 100 runs with a call; D=25, 50, 100; per subtype), PB. | **FAIL** | Holds in 0 of 9 cells; 100 of 100 runs with a call in every cell. | `rev3/results/A3B3_SUMMARY.json: A3.<P>.<D>.PB.{any,n}` |
| REV3-A3b | Same rule, LCP. | **MIXED** | Holds in 6 of 9 cells (D=25, 50 all subtypes); fails at D=100 in all three. Runs with a call (P1, P2, P3): D=25: 2, 2, 2; D=50: 0, 0, 0; D=100: 38, 10, 46. | `rev3/results/A3B3_SUMMARY.json: A3.<P>.<D>.LCP.{any,n}` |
| REV3-A3c | Same rule, LCPd. | **MIXED** | Holds in 7 of 9 cells; fails at D=25 in P2, P3. Runs with a call: D=25: 2, 24, 12; D=50: 0, 0, 0; D=100: 0, 0, 0. | `rev3/results/A3B3_SUMMARY.json: A3.<P>.<D>.LCPd.{any,n}` |
| REV3-A3d | Same rule, XDS with permutation p-values (LCPXd, as in DSTD; B=1e5). | **PASS** | Holds in 9 of 9 cells; 0 of 900 runs with a call (D=25: 0, 0, 0; D=50: 0, 0, 0; D=100: 0, 0, 0). Caveat (REV4 plan): with B=1e5 the smallest p (~1e-5) cannot reach the BH threshold for ~88,000 pairs. Single-split XDS with normal reference (also in file) holds in 4 of 9. | `rev3/results/A3B3_SUMMARY.json: A3.<P>.<D>.XDSperm.{any,n} (and .XDS)` |
| REV3-A3e | Same rule, XDS-M (S=5). | **MIXED** | Holds in 4 of 9 cells: D=25 none of 3, D=50 one of 3 (P1), D=100 three of 3. Runs with a call: D=25: 49, 49, 49; D=50: 9, 17, 13; D=100: 3, 6, 6. | `rev3/results/A3B3_SUMMARY.json: A3.<P>.<D>.XDSM.{any,n}` |
| REV3-A3f | Same rule, KX. | **FAIL** | Holds in 0 of 9 cells. Runs with a call: D=25: 100, 100, 100; D=50: 88, 78, 71; D=100: 25, 23, 19. | `rev3/results/A3B3_SUMMARY.json: A3.<P>.<D>.KX.{any,n}` |
| REV3-B2 | B2: “covariate adjustment does not fix pseudobulk” if mean calls ≥100 for both PBres and PBw in all subtypes (20 runs, D=all). | **FAIL** | Criterion not met (claim not supported). PBres mean calls 8,576/20,722/6,483 (≥100 in all); PBw 11.5/4,820/11.2 (&lt;100 in P1, P3). Both made ≥1 call in 20 of 20 runs in every subtype; PBw all-pair type-I 0.019/0.143/0.020. | `rev3/results/B2_summary.json: <P>.{PB,PBres,PBw}.{meanV,any,n,t05}` |
| REV3-B3 | B3: donor-centring sensitivity (clip [0.5, 2]; P3; 40 reassignments; D=200, 400, all). Descriptive comparison; no decision. | **DESCRIPTIVE** | Runs with a call (of 40) at D=200/400/all: XDS-M 1/2/1; LCPd 7/36/40; XDS 6/4/8; KX 3/2/6; PB and LCP 40/40/40. | `rev3/results/A3B3_SUMMARY.json: B3.P3.<D>.<method>.{any,n}` |

Notes and deviations:

- The appendix’s “1 of 4 for 0.03–0.1 UMI per cell” (macro AtwoLoDec*) counts only the A2x block (B1, genes 2,001–3,000). The original A2 block B2 also has a 0.03–0.1 stratum (50 genes), where XDS-M (S=5) is flagged calibrated in 4 of 4 combinations, but its expected count at 1e-5 is 19.7 (<20, not assessable at that threshold); pooling both blocks gives 5 of 8. The A2 criterion for the macro “above 0.1 UMI per cell” (23 of 24) covers B1 and B2 strata >0.1 only; A2x 0.1–0.3 UMI/cell (4 of 4) and P3 (8 of 8) are counted separately.
- A3 “XDS” in the frozen plan is the permutation version (LCPXd); the summary file also carries a single-split XDS with normal reference (reported above for context).

## REV4_PREREG (external review) with Amendments 1 and 2

- Files: `rev4/REV4_PREREG.md`, `rev4/REV4_AMENDMENT1.md`, `rev4/REV4_AMENDMENT2.md`
- Hash: no .sha256 file in bundle; sha256 of bundle files: PREREG 4c7d55279263…, AMENDMENT1 92c885d008f9…, AMENDMENT2 07edeeedb9cc…
- Scope in bundle: R1 PARTIAL (13, 12, 0 of 20 runs); R2 40 per subtype (R2c P1: 36); R3 20/10; R5 20; R6 descriptive.

| ID | Criterion | Outcome | Evidence | Source (file: key) |
|---|---|---|---|---|
| REV4-R1a | R1 (all R1 rows): “holds” if ≤2 of 20 runs have a BH call, per subtype and D=25, 50, 100. PB and LCP. | **FAIL (partial)** | PB: holds in 0 of 6 cells (25 of 25, 25 of 25, 25 of 25 runs with a call at D=25, 50, 100). LCP: 0 of 6; D=25: 19/25; D=50: 11/25; D=100: 25/25. Failures cannot reverse with more runs. | `rev4/results/REV4_SUMMARY.json: R1.<P>.<D>.{PB,LCP}.{any,n}` |
| REV4-R1b | R1, LCPd. | **MIXED (partial)** | Holds in 3 of 6 evaluated cells (P1, P2 × D=25, 50, 100); runs with a call (pooled P1+P2, 25 runs): D=25: 20/25; D=50: 5/25; D=100: 3/25. Fails at D=25 (both) and at D=50 in P2 (3 of 12); holds at D=100 (1 of 13, 2 of 12). | `rev4/results/REV4_SUMMARY.json: R1.<P>.<D>.LCPd.{any,n}` |
| REV4-R1c | R1, XDS with permutation p-values (B=2,000,000 in the plan). | **PASS (partial)** | Holds in 6 of 6 evaluated cells (P1, P2 × D=25, 50, 100); runs with a call (pooled P1+P2, 25 runs): D=25: 1/25; D=50: 1/25; D=100: 2/25. Provisional: P3 not run; B actually used is not recorded in the summary. | `rev4/results/REV4_SUMMARY.json: R1.<P>.<D>.XDSperm.{any,n}` |
| REV4-R1d | R1, XDS (single split, normal-moment reference, recomputed on the same runs). | **MIXED (partial)** | Holds in 5 of 6 evaluated cells (P1, P2 × D=25, 50, 100); runs with a call (pooled P1+P2, 25 runs): D=25: 0/25; D=50: 1/25; D=100: 5/25. Fails at D=100 in P2 (4 of 12). | `rev4/results/REV4_SUMMARY.json: R1.<P>.<D>.XDS.{any,n}` |
| REV4-R1e | R1, XDS-M (S=5). | **MIXED (partial)** | Holds in 4 of 6 evaluated cells (P1, P2 × D=25, 50, 100); runs with a call (pooled P1+P2, 25 runs): D=25: 13/25; D=50: 2/25; D=100: 0/25. Fails at D=25 (6 of 13, 7 of 12); holds at D=50, 100. | `rev4/results/REV4_SUMMARY.json: R1.<P>.<D>.XDSM.{any,n}` |
| REV4-R1f | R1, KX. | **MIXED (partial)** | Holds in 2 of 6 evaluated cells (P1, P2 × D=25, 50, 100); runs with a call (pooled P1+P2, 25 runs): D=25: 25/25; D=50: 20/25; D=100: 3/25. Fails at D=25 and D=50; holds at D=100 (2 of 13, 1 of 12). | `rev4/results/REV4_SUMMARY.json: R1.<P>.<D>.KX.{any,n}` |
| REV4-R1g | R1 in P3 (ADAM_GPNMB): 20 runs planned at D=25, 50, 100. | **NOT EVALUATED** | 0 of 20 runs in the bundle (R1.P3.n = 0). P1: 13 of 20 runs, P2: 12 of 20 (partial). Plan: “report the runs completed (stated as such)” if compute does not allow 20. | `rev4/results/REV4_SUMMARY.json: R1.<P>.n` |
| REV4-R2 | R2 decision: XDS-M holds by the run-level rule (n=40: ≤4 of 40 runs with a call) in each subtype (brain-bank-stratified null). | **MIXED** | Holds in 1 of 3 (P3: 2/40); fails P1 15/40 and P2 8/40. Other methods (P1, P2, P3; same rule, not a named decision): PB 40/40, 40/40, 40/40; LCP 40/40, 40/40, 40/40; LCPd 36/40, 29/40, 11/40; XDS 8/40, 12/40, 7/40; KX 23/40, 2/40, 0/40. | `rev4/results/REV4_SUMMARY.json: R2.<P>.<method>.{any,n}` |
| REV4-R2c | Amendment 2 (post hoc; rule fixed before R2c): if XDS-Mb fails ≤ the bound’s runs where XDS-M fails, attribute R2 failures to residual bank structure. | **PASS (partial)** | Rule met. XDS-M fails in P1 (14/36 runs) and P2 (8/40); XDS-Mb 4/36 (P1; exactly at the 4-run bound), 1/40 (P2), 2/40 (P3). Provisional: P1 has 36 of 40 runs. PB and bank-centred PB (PBb): 40/40 runs with a call in P2, P3 and 36/36 in P1. | `rev4/results/REV4_SUMMARY.json: R2c.<P>.{XDSM,XDSMb,PB,PBb}.{any,n}` |
| REV4-R3 | R3: PC-adjusted (k=5, 10) and weighted pseudobulk baselines, PsychAD (20 runs) and OneK1K (10 runs). Descriptive; no decision. | **DESCRIPTIVE** | Every variant made ≥1 call in every run (PsychAD 20/20; OneK1K 10/10). Mean calls PsychAD: PBres 6,483–20,722; PBw 11–4,820; PC5/PC10 4,226–9,222. OneK1K: PB 18,989–23,114; PBw 1,096–2,590; PC 4,084–17,974. | `rev4/results/REV4_SUMMARY.json: R3.<P>.<method>.{any,n,meanV}; R3o.<B>.<method>` |
| REV4-R4 | R4: real-donor consequences on PsychAD: PB vs XDS-M calls and overlap, r_w and (r_PB − rho) by class, hub overlap. Descriptive; no decision. | **DESCRIPTIVE** | Calls PB 53,263/57,280/37,325; XDS-M 56,075/60,272/37,866. PB-only 9,914/7,100/11,357; XDS-M-only 12,726/10,092/11,898; shared 43,349/50,180/25,968. Top-20 hub overlap 11–14; top-5,000-pair overlap 53–70%; Spearman(r_PB−rho, r_w) −0.47 to −0.30. | `rev4/results/REV4_SUMMARY.json: R4.<P>.{calls,classes,hub_overlap_top20,top5000_pair_overlap,spearman_delta_rw}` |
| REV4-R5 | R5 decision: XDS-Mgu holds by the run-level rule (n=20: ≤2 of 20 runs with a call) in each subtype (DC donor-centred null). | **PASS** | Holds in 3 of 3: XDS-Mgu 0/20, 0/20, 0/20 runs with a call. Same rule on others (not named decisions): XDS-Mg 1, 4, 2 of 20 (holds in 2 of 3; fails P2); XDS-g 0, 1, 0; PB 20, 20, 20 (mean calls 4,893/3,598/1,075). | `rev4/results/REV4_SUMMARY.json: R5.<P>.{XDSMgu,XDSMg,XDSg,PB}.{any,n,meanV}` |
| REV4-R6 | Amendment 1 (R6a–c, Mantis maps): donor maps, pair-level and gene-level Spearman statistics. “Descriptive only; no pass/fail decision.” | **DESCRIPTIVE** | P1 profile norm vs n_d (Spearman): null PB −0.98, design-standardized −0.03; real PB −0.81, standardized −0.26. Pairs vs r_w: null PB 0.26, real PB 0.50, real XDS-M 0.51. Genes: null false degree vs PB-only degree 0.10. Mantis-cluster Kruskal–Wallis tests: not determinable from bundle. | `rev4/results/REV4_SUMMARY.json: R6.{null_L,null_R,real_L,real_R,pairs,genes}` |

Notes and deviations:

- R1 is PARTIAL: 13 (P1), 12 (P2) and 0 (P3) of the planned 20 runs per subtype are in the summary (25 of 60). The appendix list item for REV4 does not say so. Verdicts marked partial (table: *; Markdown: (partial)) rest on these runs. Because the bound is at most 2 runs for n = 12, 13 and 20 alike, a cell with ≥3 runs with a call has failed for good (FAIL partial is final for those cells); “holds” cells (PASS and MIXED partial) are provisional.
- Discrepancy for the authors: the appendix says “LCP and LCPd fail at 25–100 donors”. In the R1 data LCPd holds at D=100 in both evaluated subtypes (3/25 runs) and in P1 at D=50; it fails at D=25 (20/25 runs). The statement is accurate for LCP only.
- R1 plan: permutation tests with B = 2,000,000; the summary does not record B, so it cannot be confirmed from the bundle. R2 plan item “median library size and cells per donor by bank” is not in the summary; R2c bank shifts are in R2c.<P>.bank_shift_median. R6c Mantis-cluster results are referenced by a macro (GeneMapAppendix, in appendix.tex) that is not defined in the supplied macro files.
- R2c P1 has 36 of 40 reassignments in the bundle; the XDS-M comparison in R2c P1 (14/36) therefore differs slightly from R2 P1 (15/40).

## Verification performed

- Macro cross-check: 139 values from `nums.tex`, `rev2nums.tex` and `rev4nums.tex` were recomputed from the result files; all 139 agree (listing at the end of this file). Six further recomputed values have no macro of that name; the four that quote paper text (power ratio 0.87, DC type-I 0.048–0.052, OneK1K growth ratios 2.2/1.7/8.5 and PB calls 19,787–24,641) agree, and the two XDS-M tail-count ranges have no counterpart in the paper.
- Stored decision flags were re-derived from the underlying numbers using the frozen rules and agree in every case: `A2_SUMMARY.json` `calibrated` and `run_ok` (220 of 220 entries), `pop_results.json` `holds` (100 of 100 cells), `xu_results.json` E1 `holds` (12 of 12), `xs_results.json` and `d_results.json` null `holds` (6 of 6), `fix_results.json` LCP/LCP2 flags and adoption decision, `onek_summary.json` `holds` (12 of 12 cells), `xs_results.json` and `d_results.json` decision flags.
- Run-level verdicts for REV-R1, REV2-A, REV3-A3, REV4-R1/R2/R2c/R5 were counted from the per-run files or per-cell `any`/`n` counts; the counts reproduce the paper macros (`RoneXDSruns`, `RtwoHoldsXDSM`, `AthreeHoldsXDSMFifty`, `RfBankXDSM`, `RfDCguHolds`, and others).
- Not verified: the LaTeX file was not compiled (no TeX installation in the build sandbox; the tectonic package bundle could not be fetched). A static lint confirmed balanced braces and environments, four cells per row, ASCII-only content, and a whitelist of commands. It needs `\usepackage{booktabs,array}` and uses no paper macros.

## Plan hash check

The paper states each plan was hashed before its computations; a hash shows the file was not edited afterwards, not when it was written. Result of comparing bundled plan files with bundled `.sha256` files:

| Plan | File | sha256 (first 16) | Stored hash file |
|---|---|---|---|
| FIX | `prereg/FIX_PREREG.md` | `3b76d65ccf0a30c6` | match |
| XS | `prereg/XS_PREREG.md` | `6022eabd47aa38f5` | match |
| POP | `prereg/POP_PREREG.md` | `083e485480676b19` | match |
| DSTD | `prereg/DSTD_PREREG.md` | `485fcfbf1f438407` | match |
| XU | `prereg/XU_PREREG.md` | `15d353ed3df6fa31` | match |
| COUPLE | `prereg/COUPLE_PREREG.md` | `7c5aae02a193d631` | match |
| DC | `prereg/DC_PREREG.md` | `4d086ce2f1719972` | match |
| REV | `prereg/REV_PREREG.md` | `b79c1ae0bf2ffc5a` | match |
| ONEK | `prereg/ONEK_PREREG.md` | `cbba3c0d5704e060` | match |
| REV2 | `rev2/REV2_PREREG.md` | `1a7e4de14a1f72fb` | no stored .sha256 in bundle |
| REV2 | `rev2/REV2_AMENDMENT1.md` | `4f42dbb837b01472` | no stored .sha256 in bundle |
| REV3 | `rev3/REV3_PREREG.md` | `314e01284fb9e186` | no stored .sha256 in bundle |
| REV3 | `rev3/REV3_AMENDMENT1.md` | `14b1eb304ef76f33` | no stored .sha256 in bundle |
| REV4 | `rev4/REV4_PREREG.md` | `4c7d552792632789` | no stored .sha256 in bundle |
| REV4 | `rev4/REV4_AMENDMENT1.md` | `92c885d008f9bf66` | no stored .sha256 in bundle |
| REV4 | `rev4/REV4_AMENDMENT2.md` | `07edeeedb9cc7510` | no stored .sha256 in bundle |

## Flags for the authors

1. **REV4 / appendix wording on LCPd.** The appendix list says LCP and LCPd fail at 25–100 donors. In the R1 data LCPd fails at D=25 (20 of 25 runs with a call) but holds at D=100 in both evaluated subtypes (3 of 25 runs) and at D=50 in P1 (REV4-R1b; macros `RfLCPdFifty`, `RfLCPdHundred`). The statement is accurate for LCP only.
2. **REV4 R1 is partial (13, 12 and 0 of 20 runs in P1, P2, P3).** The appendix list item for REV4 does not say so. All “holds” verdicts for R1 (including permutation-based XDS holding) are provisional; P3 has no R1 result (REV4-R1g).
3. **REV3 A2 counts.** The appendix figure “1 of 4 for 0.03–0.1 UMI per cell” (macros `AtwoLoDec*`) counts only the A2x block (OneK1K B1, genes 2,001–3,000). The A2 block B2 has a 50-gene 0.03–0.1 stratum flagged calibrated for XDS-M (S=5) in 4 of 4 combinations, but not assessable at 1e-5 (expected count 19.7); counting both blocks gives 5 of 8 (REV3 note).
4. **REV3 validated gene range.** The frozen definition (strata calibrated in every data set, null and D) gives 0.1–1 UMI per cell for XDS-M (S=5), because the >1 UMI stratum fails once (B1, donor-centred null, D=200: 3 of 20 runs). The paper reports 23 of 24 and 1 of 4, not the range itself (REV3-A2r, derived here).
5. **REV2 thresholds.** REV2_PREREG states no threshold for parts D1, D2, D3, E and F, and no rule for XDS-Mg although Amendment 1 refers to one. The paper's “XDS-M within the bound over 33 runs” uses the Part-A rule (at most 4 of 33); the values (4, 1, 1 of 33) satisfy it.
6. **REV2 summary file is named “interim”** (`rev2/interim/REV2_SUMMARY_latest.json`). It is complete for every part and reproduces all `rev2nums.tex` macros compared here, but the bundle does not show whether a later final file exists.
7. **POP_RESULTS.md ordering.** It lists real-donor XDS calls as 32,790 / 49,852 / 55,021; `pop_results.json` gives P1 49,852, P2 55,021, P3 32,790. The table uses the JSON.
8. **DSTD criterion (d)** (PsychAD null at D=all) is not stored as a flag in `d_results.json`; it was recomputed from `pop_results.json` (same K = 3 reassignments).
9. **XS-2 is a near miss:** LCPX-S FDR in F_S2 is 0.0557 against a bound of 0.0556. The claim fails independently in F_S1 (FDR and weak-null type-I).
10. **Run-level SE reading.** The REV plans define SE at p = 0.05 (REV3 and REV4 state it; REV_PREREG and REV2_PREREG say only “binomial”). With the SE taken at the observed fraction instead, XDS in REV-R1a (11 of 100 in P1 and P3) would be borderline. The paper's reading is used throughout.

## Not determinable from the bundle

- **ONEK-4** real-donor OneK1K calls and PB-only counts: no result file in the bundle.
- **REV2-G** oracle (matched-FDP) power: no result file in the bundle.
- **REV4-R6c** Kruskal–Wallis tests across Mantis clusters and the top clusters by PB-only degree: not in `REV4_SUMMARY.json` (`R6.genes` holds only the Spearman statistics); `appendix.tex` calls a macro `\GeneMapAppendix` that is not defined in the supplied macro files.
- **REV4-R1 permutation count:** the plan specifies B = 2,000,000 for LCP, LCPd and XDS-perm; the summary does not record the B used.
- **REV-R3 overlap** of real-donor XDS-DCg calls with the pooled-centre XDS-DC calls: `rev_dc_r5.json` has `XDSg_in_pooled` without a definition.
- **XU E4** (six small designs) and **DC-3** method overlaps: E4 values are in `xu_results.json` but no criterion applies; overlaps are not tabulated.
- **Stored hashes for REV2, REV3, REV4 plans and amendments:** not in the bundle (the paper's `verify_hashes.py` reads `rev3/*.sha256`, which was not bundled); the sha256 of each bundled file is listed above.

## Plans in the bundle not covered by the table

DEV_PREREG, HELDOUT_SPEC and RARE_PREREG all contain explicit criteria; none is mentioned in the appendix list or in the bundled macro files, so they are **not reported in the paper**. Their outcomes were not re-derived here; the repository memos state:

- **DEV_PREREG** (GO/NO-GO gates 1–4): per `prereg/GO_NOGO.md`, gate 2 (validity) FAIL at D = 8, 12, 16; gate 3 (amended) PASS at D16 with D40 not yet scored; gate 4 (amended) FAIL in 15 of 60 class×design cells; gate 1 (novelty) not passed; recommendation NO-GO (written 26 Oct 2026).
- **HELDOUT_SPEC**: validity claim (D ≥ 24 designs: FDR ≤ 0.05 + 2 MCSE, type-I ≤ 0.05 + 3 MCSE), power claim (≥ 0.8× BOUND), few-donor check and weak-null/interval reporting. Results exist in `results/heldout/` and were not scored here.
- **RARE_PREREG**: calibration gate, Q1 and Q2. Per `prereg/RARE_RESULTS.md`: calibration gate FAILS in S1 (LCP FDR 0.075 > 0.05 + 2 MCSE); Q1 = NO in S1 and S2; Q2 = YES in S1, NO in S2 (S1 results descriptive because the gate failed).

## Cross-check of the paper macros against the bundle

139 macro values (nums.tex, rev2nums.tex, rev4nums.tex) were compared with the value recomputed from the result files for this table: 139 agree after normalising formatting, 0 differ.

<details><summary>All compared macros</summary>

| Macro | In macro file | Recomputed |
|---|---|---|
| `PBnullCalls` | `19{,}158--31{,}617` | `19,158–31,617` |
| `PBrepFalse` | `3{,}798--12{,}827` | `3,798–12,827` |
| `RealPBrep` | `62--80\%` | `62–80%` |
| `RealXDSrep` | `68--82\%` | `68–82%` |
| `XDSnullTI` | `0.048--0.050` | `0.048–0.050` |
| `LCPnullTI` | `0.20--0.30` | `0.20–0.30` |
| `LCPnullCalls` | `4{,}986--12{,}880` | `4,986–12,880` |
| `XUnullCalls` | `12{,}073--19{,}934` | `12,073–19,934` |
| `CouplePB` | `0.524, 0.605, 0.718` | `0.524, 0.605, 0.718` |
| `CoupleLCP` | `0.104, 0.143, 0.191` | `0.104, 0.143, 0.191` |
| `CoupleXDS` | `0.059, 0.047, 0.056` | `0.059, 0.047, 0.056` |
| `CoupleXU` | `0.060, 0.104, 0.223` | `0.060, 0.104, 0.223` |
| `RoneXDSruns` | `11--20` | `11–20` |
| `RoneXDSMruns` | `2--5` | `2–5` |
| `RoneKXruns` | `4--36` | `4–36` |
| `XDSMsmallFDR` | `0.135--0.265` | `0.135–0.265` |
| `KXsmallFDR` | `0.27--0.90` | `0.27–0.90` |
| `SeaXDSMTop` | `0.043--0.070` | `0.043–0.070` |
| `SeaKXTop` | `0.055--0.061` | `0.055–0.061` |
| `XDSMLFDR` | `0.043--0.054` | `0.043–0.054` |
| `DCgTI` | `0.048--0.052` | `0.048–0.052` |
| `DCgReal` | `117, 0, 0` | `117, 0, 0` |
| `PosXYmaxp` | `5\times10^{-40}` | `5e-40` |
| `PosYYmaxp` | `3\times10^{-8}` | `3e-08` |
| `FalseRmed` | `0.10--0.13` | `0.10–0.13` |
| `FalseRgtOne` | `55--84\%` | `55–84%` |
| `FalseRgtTwo` | `4--14\%` | `4–14%` |
| `RsixPBcalls` | `8{,}708--22{,}872` | `8,708–22,872` |
| `PoolRatio` | `37--46\%` | `37–46%` |
| `RtwoHoldsXDSM` | `8 of 9` | `8 of 9` |
| `RtwoHoldsKX` | `9 of 9` | `9 of 9` |
| `RtwoHoldsLCPd` | `3 of 9` | `3 of 9` |
| `RtwoHoldsXDS` | `0 of 9` | `0 of 9` |
| `RtwoAnyXDSM` | `41/900` | `41/900` |
| `RtwoAnyXDS` | `260/900` | `260/900` |
| `RtwoAnyLCPd` | `156/900` | `156/900` |
| `RtwoAnyKX` | `51/900` | `51/900` |
| `RtwoRangeXDS` | `11--42` | `11–42` |
| `RtwoRangeXDSM` | `0--10` | `0–10` |
| `AfiveRtwo` | `0.89--1.00` | `0.89–1.00` |
| `RealXDSMjac` | `0.86--0.95` | `0.86–0.95` |
| `RealXDSMrep` | `74--85\%` | `74–85%` |
| `RealXDSMcalls` | `56,075, 60,272, 37,866` | `56,075, 60,272, 37,866` |
| `DCdcPB` | `995--5,753` | `995–5,753` |
| `DCdcPBTI` | `0.13--0.21` | `0.13–0.21` |
| `DCdcXDSTI` | `0.048--0.053` | `0.048–0.053` |
| `DCrealPB` | `5,806, 2,441, 2,561` | `5,806, 2,441, 2,561` |
| `DCrealM` | `465, 0, 0` | `465, 0, 0` |
| `DCrealU` | `351, 341, 0` | `351, 341, 0` |
| `DCrealG` | `117, 0, 0` | `117, 0, 0` |
| `DCtwoFDRPB` | `0.485` | `0.485` |
| `DCtwoPowPB` | `0.962` | `0.962` |
| `DCtwoFDRM` | `0.071` | `0.071` |
| `DCtwoPowM` | `0.264` | `0.264` |
| `DCtwoFDRU` | `0.031` | `0.031` |
| `DCtwoPowU` | `0.280` | `0.280` |
| `DCtwoFDRG` | `0.077` | `0.077` |
| `DCtwoPowG` | `0.094` | `0.094` |
| `OnekEanyXDSM` | `4/33, 1/33, 1/33` | `4/33, 1/33, 1/33` |
| `OnekEPBlin` | `19,700, 13,894, 16,827` | `19,700, 13,894, 16,827` |
| `OnekEPB` | `22,889, 22,948, 20,035` | `22,889, 22,948, 20,035` |
| `FtwoFDRXDSM` | `0.043` | `0.043` |
| `FfiveFDRXDSM` | `0.047` | `0.047` |
| `FtwoFDRLCPd` | `0.067` | `0.067` |
| `FtwoFDRXDS` | `0.088` | `0.088` |
| `FtwoPowKX` | `0.17` | `0.17` |
| `FfivePowKX` | `0.67` | `0.67` |
| `AtwoHiDecM` | `23 of 24` | `23 of 24` |
| `AtwoHiDecS` | `16 of 24` | `16 of 24` |
| `AtwoHiDecMone` | `23 of 24` | `23 of 24` |
| `AtwoHiDecMten` | `22 of 24` | `22 of 24` |
| `AtwoHiDecPB` | `0 of 24` | `0 of 24` |
| `AtwoLoDecM` | `1 of 4` | `1 of 4` |
| `AtwoLoDecS` | `0 of 4` | `0 of 4` |
| `AtwoLoDecMone` | `4 of 4` | `4 of 4` |
| `AtwoLoDecMten` | `3 of 4` | `3 of 4` |
| `AtwoLoDecPB` | `0 of 4` | `0 of 4` |
| `AtwoPthreeDecM` | `8 of 8` | `8 of 8` |
| `AtwoPthreeDecS` | `0 of 8` | `0 of 8` |
| `AtwoPthreeDecPB` | `0 of 8` | `0 of 8` |
| `AthreeXDSMtwentyfive` | `49\%` | `49%` |
| `AthreeXDSMfifty` | `9--17\%` | `9–17%` |
| `AthreeXDSMhundred` | `3--6\%` | `3–6%` |
| `AthreeLCPhundred` | `10--46\%` | `10–46%` |
| `AthreeLCPdtwentyfive` | `2--24\%` | `2–24%` |
| `AthreeLCPsmall` | `0--2\%` | `0–2%` |
| `BtwoPBw` | `11--4,820` | `11–4,820` |
| `BtwoPBres` | `6,483--20,722` | `6,483–20,722` |
| `BthreeXDSM` | `1/40, 2/40, 1/40` | `1/40, 2/40, 1/40` |
| `BthreeLCPd` | `7/40, 36/40, 40/40` | `7/40, 36/40, 40/40` |
| `RfLCPTwentyFive` | `19/25` | `19/25` |
| `RfLCPdTwentyFive` | `20/25` | `20/25` |
| `RfXDSpermTwentyFive` | `1/25` | `1/25` |
| `RfXDSMTwentyFive` | `13/25` | `13/25` |
| `RfKXTwentyFive` | `25/25` | `25/25` |
| `RfLCPFifty` | `11/25` | `11/25` |
| `RfLCPdFifty` | `5/25` | `5/25` |
| `RfXDSpermFifty` | `1/25` | `1/25` |
| `RfXDSMFifty` | `2/25` | `2/25` |
| `RfKXFifty` | `20/25` | `20/25` |
| `RfLCPHundred` | `25/25` | `25/25` |
| `RfLCPdHundred` | `3/25` | `3/25` |
| `RfXDSpermHundred` | `2/25` | `2/25` |
| `RfXDSMHundred` | `0/25` | `0/25` |
| `RfKXHundred` | `3/25` | `3/25` |
| `RfBankXDSM` | `15/40, 8/40, 2/40` | `15/40, 8/40, 2/40` |
| `RfBankPB` | `40/40, 40/40, 40/40` | `40/40, 40/40, 40/40` |
| `RfBankLCP` | `40/40, 40/40, 40/40` | `40/40, 40/40, 40/40` |
| `RfBankLCPd` | `36/40, 29/40, 11/40` | `36/40, 29/40, 11/40` |
| `RfBankXDS` | `8/40, 12/40, 7/40` | `8/40, 12/40, 7/40` |
| `RfBankKX` | `23/40, 2/40, 0/40` | `23/40, 2/40, 0/40` |
| `RfBankXDSMb` | `4/36, 1/40, 2/40` | `4/36, 1/40, 2/40` |
| `RfBankXDSMc` | `14/36, 8/40, 2/40` | `14/36, 8/40, 2/40` |
| `RfBankXDSMbHolds` | `3` | `3` |
| `RfPCadj` | `4{,}226--9{,}222` | `4,226–9,222` |
| `RfOnekPB` | `18{,}989--23{,}114` | `18,989–23,114` |
| `RfOnekPBw` | `1{,}096--2{,}590` | `1,096–2,590` |
| `RfOnekPC` | `4{,}084--17{,}974` | `4,084–17,974` |
| `RfRealPB` | `53{,}263, 57{,}280, 37{,}325` | `53,263, 57,280, 37,325` |
| `RfRealXDSM` | `56{,}075, 60{,}272, 37{,}866` | `56,075, 60,272, 37,866` |
| `RfRealPBonly` | `9{,}914, 7{,}100, 11{,}357` | `9,914, 7,100, 11,357` |
| `RfRealXonly` | `12{,}726, 10{,}092, 11{,}898` | `12,726, 10,092, 11,898` |
| `RfRealShared` | `43{,}349, 50{,}180, 25{,}968` | `43,349, 50,180, 25,968` |
| `RfHubOverlap` | `11--14` | `11–14` |
| `RfTopOverlap` | `53--70\%` | `53–70%` |
| `RfDCguHolds` | `3` | `3` |
| `RfDCgHolds` | `2` | `2` |
| `RfDCgu` | `0, 0, 0` | `0, 0, 0` |
| `RfDCg` | `1, 4, 2` | `1, 4, 2` |
| `RfDCPB` | `1{,}075--4{,}893` | `1,075–4,893` |
| `RsixNormPB` | `$-$0.98` | `−0.98` |
| `RsixNormR` | `$-$0.03` | `−0.03` |
| `RsixRealNormPB` | `$-$0.81` | `−0.81` |
| `RsixRealNormR` | `$-$0.26` | `−0.26` |
| `RsixPairNull` | `0.26` | `0.26` |
| `RsixPairRealPB` | `0.50` | `0.50` |
| `RsixPairRealX` | `0.51` | `0.51` |
| `RsixPairs` | `89{,}676` | `89,676` |
| `RsixNullCalls` | `12{,}109` | `12,109` |

</details>
