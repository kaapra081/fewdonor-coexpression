# Plan 8 (4 Nov 2026, ~18:50 UTC). POST HOC, written after an external review of the submitted draft and BEFORE any of the computations below.

Reviewer points addressed: (1) XDS-Mh, the recommended test, was developed on PsychAD and OneK1K and was never run on a held-out cohort
or on the brain-bank-stratified null where XDS-M failed; (2) power comparisons were reported at unequal achieved FDR. Code: `rev5/rev5_f.py`.
XDS-Mh is exactly the Plan 6d definition (`rev5_d.variants`): in each of the 5 splits, each half's design-standardized values are
residualized across donors on [1, log(kappa^h S^h), log n^h] before the cross-products are formed. Nothing below changes XDS-Mh.

## 8A. Held-out cohort: lupus (Perez et al. 2022; Plan 7 data, unchanged)
XDS-Mh was never computed on the lupus data. Same 20 reassignments (seeds) as Plan 7 for each null; LM (classical monocytes) and FC (pDC).
- A: donor-centred null (Plan 7 mode A seeds). XDS-M is recomputed as a reproduction check against Plan 7.
  **Rule (as Plan 6c/6d):** XDS-Mh holds in a cell type if <= 2 of 20 runs have >= 1 BH call (q = 0.05) AND mean per-pair type-I at
  0.05 is <= 0.055. Reported for LM and FC separately; "held in the held-out cohort" requires both.
- D: library-preserving null (Plan 7 mode D seeds). Same rule, reported; XDS-M failed it in Plan 7 (multiplexed pools), and XDS-Mh does
  not model pools, so the expectation is failure in LM. Whatever the outcome, it is reported.
- R: real donors (split seed 44,000,001). Descriptive: XDS-Mh calls, fraction of XDS-M calls kept, ISG-pair calls.

## 8B. PsychAD brain-bank-stratified null (Plan 5 R2; same 40 reassignments per subtype, seeds 64,000,000 + 10k + subtype index)
- XDS-M recomputed (reproduction check against rev4/out/R2_*.json), XDS-Mh, and XDS-Mhb = XDS-Mh with the pseudo-donor's bank
  indicators (H, R; M reference) added to each half's residualization columns.
  **Rule (as Plan 5 R2):** a test holds in a subtype if <= 4 of 40 runs have >= 1 BH call AND mean per-pair type-I <= 0.055.
  **Interpretation fixed now:** if XDS-Mh holds in all 3 subtypes, the paper states it held the bank null. If XDS-Mh fails anywhere and
  XDS-Mhb holds in all 3, the paper recommends adding known batch indicators (XDS-Mhb) and says XDS-Mh alone does not handle bank
  structure. If both fail, the paper states that neither handles bank structure and keeps bank structure as an open limitation.

## 8C. Power at matched achieved FDR (PsychAD planted signal; Plan 6d S runs: same seeds 70,000,000 + 10k + subtype index, 20 runs/subtype)
- Methods: PB, KX, XDS-M, XDS-Mh. Per run: (i) BH at q = 0.05: FDP and power (reproduction check against rev5/out/d_S_*.json for
  XDS-M, XDS-Mh); (ii) oracle-matched power: pairs ordered by p-value, the longest prefix whose false-discovery proportion is <= L,
  power = true pairs in that prefix / 180, for L = 0.05 and 0.10. Reported as mean and SD over runs per subtype. Descriptive (no rule);
  it answers whether pseudobulk's ranking, rather than its calibration, is the problem.

## Reporting
All outcomes go into the paper (main text if they change a claim, appendix otherwise) and the repository (rev5/out/f_*.json,
rev5/out/plan8_summary.json), labelled post hoc. No other analysis is run under this plan.
