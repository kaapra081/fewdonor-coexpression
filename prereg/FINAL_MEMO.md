# Final decision and status (28 Oct 2026)

## Decision
Submit to RECOMB 2027 as a **limits-and-design-rule paper**. Use abstract v5 and the paper as it now stands.

Reasons:
- Registering is free and does not commit us to anything.
- The paper is rigorous and every claim is now accurate.
- If RECOMB rejects it, the same manuscript fits a statistics or methods journal (Biostatistics, Genome Biology) with little change.

Do not reframe it as a "better network method" paper; the real-data evidence does not support that.

Expected RECOMB outcome: borderline to weak reject. The weak points are novelty (the leak decomposition is classical) and the missing proof that the pairs only pseudobulk calls are within-donor programs. The strong points are rigour, frozen held-out and rare-cell analyses, and honest negative results.

## What the paper claims (every number generated from repository files)
- **The leak.** Pseudobulk covariance = between-donor covariance + κ̄ × within-donor covariance. It is estimable per donor (count-scale BCES/Muthén, credited).
- **The bound.** For a per-pair test of donor means, power is capped by D_eff = Σ ρρ. The median D_eff/D is 0.53 at 40 cells for the 417 analysed genes. LCP's asymptotic efficiency relative to the bound is 0.925–0.997 across designs.
- **Where LCP is calibrated (design map, Table 3):**
  - Every design with ≥ 24 donors and ≈ 40+ cells per donor (FDR 0.041–0.051).
  - At 16 donors only with ~100 cells.
  - Not at 15 cells per donor, even with 79 donors (0.075).
  - Pseudobulk is never calibrated (0.08–0.55).
- **Rare SEA-AD supertype (79 donors, median 15 nuclei; decision rule frozen before the run):**
  - 39% of pseudobulk calls are not LCP calls.
  - Removing them changes 10 of the top-20 network hubs (Q2 = YES).
  - The rule's test that they are within-donor programs fails (Q1 = NO).
  - The calibration gate fails, so these results are descriptive.
- **Negative results, all stated:**
  - ≤ 16 donors with 40 cells.
  - Depth–state coupling (FDR 0.22). Coupling of similar size is present in both real cohorts.
  - The weak null (0.06–0.10).
  - Intervals undercover, even on true donor means.
  - A strong low-rank donor factor is present in the real data.

## Files
- Paper: main.pdf. The main text ends on page 9 of the 10 allowed; the appendix is unlimited.
- Abstract: abstract.txt v5.
- Repository: https://github.com/kaapra081/fewdonor-coexpression, private until you make it public.
  - prereg/ holds every analysis plan, its hash and its outcome (DEV_PREREG, HELDOUT_SPEC, RARE_PREREG, RARE_RESULTS, GO_NOGO, REVIEW_RESPONSE).
  - `python paper/make_results.py && python paper/make_rare.py` regenerates every number from results/.

## What you must do
1. Register on EasyChair by **Wed 28 Oct 23:59 AoE** (= Thu 29 Oct 04:59 PDT), with the title and abstract from abstract.txt v5 and your author names.
2. Add the author names and affiliations to paper/main.tex (line 24).
3. Make the repository public.
4. Upload the PDF by **Wed 4 Nov 23:59 AoE**.
