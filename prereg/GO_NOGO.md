# GO / NO-GO for RECOMB 2027 (recomb_v2)

Written Mon 26 Oct 2026, ~19:30 PDT, from the formal development runs scored against DEV_PREREG.md (sha256 010de043… for the original file, a860bfa4… after change 1; pushed to the private GitHub repo as commits d00b9c2 and eee7366 before any formal result was read). The 24-donor, 40-donor, strong-effect and stress designs of the *test* runs are still running locally (the Kaggle copy was lost when its container restarted). The outcome below does not depend on them, because gates 2 and 4 already fail in completed designs.

## Headline

**NO-GO.** The proposed test (LCP) is not valid in the BH-relevant tail at 8–16 donors. The pre-registered FDR criterion fails at D = 8, 12 and 16. The amended interval gate also fails. I recommend not registering a RECOMB 2027 abstract for this method.

## Gate by gate (R = 200 replicates per design; MCSE in brackets)

**Gate 2, validity: FAIL.** The criterion is FDR ≤ 0.05 + 2·MCSE and type-I ≤ 0.05 + 3·MCSE in WW, SS, HH and NN.

| Design | LCP FDR | LCP power | Worst null type-I at 0.05 |
|---|---|---|---|
| D8 | 0.876 (0.015) | 0.005 | SS 0.065 (0.004) |
| D12 | 0.246 (0.028) | 0.003 | SS 0.053 (0.003) |
| D16 | 0.089 (0.017) | 0.015 | WW 0.052 (0.002) |

- Type-I at α = 0.05 is close to nominal from 12 donors up. The extreme tail is not.
- The rate of null p < 10⁻³, relative to nominal, is 2.9–19 at D8, 1.3–2.9 at D12 and 1.1–1.6 at D16.
- Most false discoveries come from independent (NN) pairs, so the excess is not only leakage.
- The likely cause is that donors are not exchangeable in shape. A donor with few cells has heavy-tailed noise in both genes at once, and the identity pairing multiplies that noise together, whereas random pairings dilute it. Standardization equalizes variances but not tails.
- The ORACLE test (true donor means) is valid: FDR 0.045–0.047 at D12–D16. So the failure comes from measurement noise and estimated standardization, not from the permutation idea itself.
- For context, every comparator other than CMB_mom also fails here. At D16: pseudobulk FDR 0.22, E5 0.14, CMB_naive 0.19. CMB_mom is valid but has power 0.000.

**Gate 3, usefulness (amended): PASS at D16; D40 not yet scored.**

- At D16, BOUND power is 0.016, below 0.05, so the fallback criterion applies: LCP power must be at least BOUND power minus 2·MCSE. LCP reaches 0.015 against BOUND 0.016, so it passes.
- Theorem 1 predicts that no per-pair test can beat BOUND. ORACLE power at D16 is 0.21, so the original gate 3 fails, as Theorem 1 says it must.

**Gate 4, intervals (amended): FAIL.** The amended criterion is BOOT_T coverage ≥ ORACLE_BOOT_T coverage − 0.02 in every class.

- It fails in 15 of the 60 class × design cells of the six non-stress designs. Examples:
  - D24 OO_same −0.042;
  - D8 QQ_same −0.036;
  - D24 OO_opp −0.032;
  - D40 OO_same −0.029;
  - D24 SS −0.028;
  - D16_strong BB −0.028;
  - D16 OO_same −0.027.
- The original gate (≥ 0.93) fails nearly everywhere. EQ covers 0.86–0.90 for between-donor pairs at D8–D16, and CMB_jack 0.83–0.87.
- Even the oracle bootstrap-t interval on the true donor means covers only 0.86–0.95 for between-donor pairs (0.89 for independent pairs at D16).

**Gate 1, novelty: does not pass on present evidence.**

- novelty_v2.md shows that the moment correction (BCES 1996), weighted debiasing (dIVW) and permutation testing all exist already.
- The new pieces were: the count-derived per-donor leak (Lemma 1); a test calibrated at few donors; the per-pair information bound with D_eff = Σ_d ρ_dj ρ_dk (Theorem 1). The second has now failed its own validity gate.
- Lemma 1 and Theorem 1 are correct and useful. A skeptical reviewer would still call them an adaptation of attenuation theory and of BCES, not a RECOMB-level methods advance.

## What was learned (all from development data; true as stated)

1. **Leak formula (proved).** Donor-mean (pseudobulk) covariance targets Σ_B + κ̄ Σ_W, with κ̄ ≈ mean 1/n_d. This explains pseudobulk FDR of 0.22 at 16 donors.
2. **Information bound (proved for the Gaussian working model).**
   - Any per-pair between-donor test with D donors has, locally, the power of an oracle with D_eff = Σ_d ρ_dj ρ_dk donors.
   - At about 40 cells per donor, D_eff ≈ 0.3·D.
   - Development data fit this: at D16, ORACLE power is 0.21 while BOUND is 0.016 and LCP 0.015.
   - Few-donor, few-cell designs cannot answer pair-level between-donor questions, whatever the method. Adding donors helps much more than adding cells.
3. **Calibration fails at the tail.** At 8–16 donors, every between-donor test that uses noisy donor means is anti-conservative in the BH-relevant tail, ours included; only CMB_mom is not, and it has no power.
4. **Intervals fail even for the oracle.** Linear-scale between-donor covariance intervals undercover under skewed donor effects even when computed on the true donor means, at 16–84 donors.

These four results, together with the round-3 benchmark, make a coherent negative/limits paper. That paper suits a journal, not the RECOMB main track, which requires a new method that works.

## Options for you

1. **Do not register (recommended).** Consolidate the round-3 benchmark and these limits results into one journal paper, for example in Genome Biology or Bioinformatics.
2. **Register an abstract anyway.** It would have to describe a limits/benchmark paper, which RECOMB's call explicitly deprioritizes ("application of existing methods"). The abstract must closely match the final paper and must not be a placeholder. I don't recommend this.
3. **Continue developing the test for the journal paper.** Candidate repairs are a cell-split leak-free version, or restricting permutations to donors with similar cell counts. These are not for 4 Nov.

Nothing has been registered, submitted or made public.

## Update, Mon 26 Oct 22:45 PDT: all development designs finished, and your decision

**Final development numbers for LCP** (R = 200 per design; MCSE in brackets):

| Design | FDR | Power | BOUND power | ORACLE power | Pseudobulk FDR / power | Moment-corrected power |
|---|---|---|---|---|---|---|
| D24 | 0.046 (0.006) | 0.105 | 0.108 | 0.531 | 0.140 / 0.168 | 0.000 |
| D40 | 0.050 (0.002) | 0.463 | 0.458 | 0.851 | 0.129 / 0.494 | 0.010 |
| D16_strong | 0.065 (0.006) | 0.187 | 0.186 | 0.578 | 0.138 / 0.321 | 0.000 |
| D16 stress | 0.230 (0.026) | 0.012 | 0.008 | 0.224 | 0.256 / 0.034 | 0.000 |

**Gates on the final numbers:**
- **Gate 2 (validity):** passes at D24 and D40. It fails at D8, D12, D16 and D16_strong (0.065 > 0.05 + 2·0.006).
- **Gate 3 (amended):** passes. LCP power is at least 0.98 × BOUND everywhere, and 0.54 × ORACLE at D40.
- **Gate 4:** fails.
- **Overall:** NO-GO under the pre-registered rule, unchanged.

**Your decision (26 Oct, ~21:00 PDT):** proceed to a RECOMB submission anyway ("we need recomb… finish it off"). The paper is scoped to what the evidence supports:
1. the leak lemma;
2. the information bound;
3. a test calibrated from ~24 donors that reaches the bound;
4. explicit negative results at ≤16 donors and for intervals.

The frozen held-out evaluation (HELDOUT_SPEC.md) tests the ≥24-donor validity claim on new panels, seeds and real-donor designs. Any failure there will be reported in the paper.

## Held-out evaluation (frozen spec HELDOUT_SPEC.md, commit c986bab), completed 27 Oct 2026

R = 200 per design, panel seed 4112026, seeds 5,000,000 + 100,019·k + rep. HO_lee40, HO_seaad24, HO_D32 and HO_weak40 were
partly run on Kaggle; that machine disconnected before results were downloaded, so all six designs were (re)run in full on
the local machine with the same frozen code and seeds. Summary: dev/ho_results/summary_all.csv.

| Design | LCP FDR (MCSE) | LCP power / BOUND power | PB FDR |
|---|---|---|---|
| HO_seaad24 | 0.046 (0.002) | 0.554 / 0.553 | 0.082 |
| HO_seaad24_cap40 | 0.041 (0.004) | 0.122 / 0.117 | 0.170 |
| HO_lee40 | 0.050 (0.002) | 0.559 / 0.547 | 0.190 |
| HO_D32 | 0.048 (0.004) | 0.307 / 0.300 | 0.139 |
| HO_weak40 | 0.048 (0.009) | 0.025 / 0.029 | 0.299 |
| HO_D16_n100 (few-donor check, no claim) | 0.040 (0.005) | 0.077 / 0.085 | 0.133 |

- Validity claim (FDR ≤ 0.05 + 2·MCSE; type-I ≤ 0.05 + 3·MCSE in WW, SS, HH, NN): PASSES in all five D ≥ 24 designs.
- Power claim (≥ 0.8 × BOUND): PASSES in all five (ratios 0.89–1.04).
- Weak null (QQ_cross), reported, no claim: LCP type-I 0.062–0.099 (anti-conservative, as in development).
- Pseudobulk has higher power than LCP in HO_seaad24 (0.659 vs 0.554) where donors have hundreds of cells and the leak is small.
