# Response to the critical abstract review (28 Oct 2026)

Each point was checked against the result files at commit 0b728fb, before any change was made.

Verdict labels used below:
- **Correct**: the reviewer's statement matches the files.
- **Partly correct**: the substance is right but some detail is off.

The fixes are in abstract v4 and the revised paper.

## A. Factual errors

| # | Verdict | Evidence | Fix |
|---|---|---|---|
| 1 | Correct | `run_heldout.py`: only `HO_seaad24` uses `pool_n="real"`. `HO_seaad24_cap40` draws cells from LogNormal(median 40); `HO_lee40` from LogNormal(median 60). | The abstract no longer describes the designs' origin. The paper says: one design uses real cell numbers and depths, two use real depths only, three use fitted depths. |
| 2 | Correct | `deff_real.json`: median D_eff/D at 40 cells is 0.527 (SEA-AD) and 0.532 (Lee); noise-to-signal 0.30 / 0.31 on the 417 genes. `real_between_variance_calibration.json`: noise-to-signal 0.98 / 0.69 on the 2,007 genes. | The abstract gives D_eff/D 0.90–0.98 at real cell numbers and 0.53 at 40 cells, both for the 417 genes. Section 4 reports both gene sets explicitly: 16 donors are worth 8.4–8.5 oracle donors for the 417 genes and 4.1–5.6 for the 2,007 genes. The earlier "four to five" was also wrong for Lee (5.6). |
| 3 | Correct | `summary_all.csv`: MC FDR 0.000 with power 0.000 at D = 8–16. Ablation: no-leak LCP FDR 0.045 at D = 16. NMO FDR 0.065 at D = 16. | Claim removed. The paper now states which variants kept FDR near 0.05 at 16 donors and why: MC makes no discoveries; no-leak LCP loses within-only protection (type-I 0.108). |
| 4 | Correct | Theorem 1 is stated for a bivariate Gaussian on donor means with known Ψ_d, at σ = 0 with no leak. | The abstract now says "a test that uses only the two genes' donor means … with known sampling variances and no leak". Same wording in the theorem and intro. The section is renamed "What a per-pair test of donor means can detect". New title: "…a leak-corrected permutation test and the limits of per-pair inference". |

## B. Narrative

| # | Verdict | Evidence | Fix |
|---|---|---|---|
| 5 | Partly correct | No-leak FDR is 0.045 / 0.043 / 0.058 at D = 16 / 24 / 40. **The reviewer missed that 0.058 at D = 40 fails FDR control**: MCSE 0.0023, about 3.4 MCSE above 0.05. Full LCP gives 0.049. The leak term's effect on FDR grows with D, because within-only pairs become detectable as donors are added. | Abstract and paper say: standardization and permutation remove most of pseudobulk's excess; the leak term protects within-only pairs (0.048–0.054 with it vs 0.11–0.21 without) and is needed for FDR control at 40 donors. The paper also reports the D = 16 reversal: dropping the leak term lowers FDR there. |
| 6 | Correct | `bound_test` in `methods_v2.py` is LCP's statistic on true means plus noise, with true variances and no leak. | Renamed NMO (noise-matched oracle) in text, tables and figures, with a note that the plan and code call it BOUND. **New post hoc computation** (`diag/diag_sim.py`): Proposition 1's asymptotic efficiency of LCP relative to the Theorem-1 bound is 0.963–0.997 across our designs. Median D_eff/D in the development designs is 0.62–0.63, against 0.53 in the real data at 40 cells. The 1.04 ratio is now described as LCP vs NMO, not vs the bound. |
| 7 | Correct | `HO_seaad24`: PB FDR 0.082, power 0.659; LCP 0.046, power 0.554. The 0.30 upper end comes from `HO_weak40`, where PB power is 0.045. | Both facts are stated in the held-out paragraph. The abstract no longer quotes the 0.08–0.30 range. |
| 8 | Correct, except one rounding | Dev D24 power is 0.1052 and held-out cap40 is 0.122, so the range is 0.11–0.12. The review wrote 0.10. | Abstract: "Power is low with few cells (0.11–0.12 …)". The paper gives 0.105 and 0.463. |
| 9 | Correct | Weak-null type-I: dev 0.061–0.083 at all D (0.072–0.077 at D = 24–40); held-out 0.062–0.099. | Added to the abstract (0.06–0.10) and to Limitations. |
| 10 | Correct | `real_summarize.py`: the share is computed over PB-called pairs. | The metric is now defined in the text, which also says that for within-only pairs the leak is the whole signal. The abstract says "among pseudobulk calls". |

## C. Theory

| # | Verdict | Evidence | Fix |
|---|---|---|---|
| 11 | Correct, and it goes further | The proof outline left γ_d uncorrected. It also claimed the effect of nuisance estimation is O(D^-1/2), but ψ̂_d and κ_d Ĉ_d are per-donor quantities with O(n_d^-1/2) error. A further term (the variance of the leak estimate and its covariance with r r) was also missing. | Theorem 2(b) is now stated for D → ∞ **and** min n_d → ∞. With n_d fixed, the null variance is 1 + O(mean 1/n_d), which the permutation reference does not reproduce. The proof outline is corrected. The introduction says "as donors and cells per donor grow". The D = 16 ablation reversal is linked to this term. |
| 12 | Correct | The body already credited BCES and Muthén; the abstract did not. | The abstract now says "a count-scale instance of the classical multilevel (BCES-type) measurement-error decomposition". Intro point (i) cites Muthén 1994 and Akritas & Bershady 1996 directly. |

## D. Assumptions on real data (new post hoc analyses, `diag/`)

13. **Depth–state coupling: correct, and the finding is unfavourable.**
    - **Statistic:** the pooled within-donor correlation between log library size and library-normalized expression. It has expectation zero under Lemma 1.
    - **Values:**
      - SEA-AD: median 0.046 (|·| 0.054).
      - Lee: median −0.007 (|·| 0.029).
      - The coupled genes of the stress design have median 0.045.
      - |c| exceeds that value for 58% (SEA-AD) and 33% (Lee) of the 417 genes, and exceeds 0.1 for 12–13%.
    - **Difference from the stress design:** the within-donor spread of log depth (0.57–0.61) matches the null simulation (0.56), not the stress design (1.17).
    - **Reported as:** depth–state coupling of stress-design size is common in both cohorts, so the stress FDR of 0.22 is a warning for the real-data calls. Added to the abstract in one clause.
14. **Low-rank donor factor: correct.**
    - PC1 of the standardized donor means carries 21% (SEA-AD) and 28% (Lee) of the variance, against at most 2.6% after per-gene permutation.
    - It is not explained by mean library size, cell number, sex or SEA-AD mitochondrial fraction (|r| ≤ 0.19).
    - The text now attributes the 40–49% call rate to a likely shared donor factor, biological or technical. "Technical" is added in the abstract. No adjustment was made, because it changes the estimand.

## E. Strategic

15. **Correct, and not fixable in text.** The Discussion now says plainly that we have not yet shown a real-data setting where LCP changes a conclusion. The most useful next analysis is a naturally rare SEA-AD cell type (84 donors, tens of cells each). It needs a new public download. The abstract does not promise it.

## F. Presentation
- The abstract is about 370 words, not 300, because it keeps every material negative result.
- The clause "held-out designs frozen after development" is included.
