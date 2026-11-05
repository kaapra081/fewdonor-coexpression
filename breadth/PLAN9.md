# Plan 9 (5 Nov 2026, ~05:15 UTC). POST HOC breadth check, written after an external review asked for more independent cohorts, and BEFORE any data below were read.

Question: do the main null results (pseudobulk fails; XDS-M and XDS-Mh hold) carry over to cohorts never used in this project?
Code: `breadth/breadth_job.py` (selection, extraction, analysis in one Modal job). Nothing in XDS-M or XDS-Mh changes.

## 9.1 Selection (fixed now; metadata only)
CELLxGENE Census 2025-11-08, Homo sapiens, is_primary_data == True. A setting is a (dataset, cell_type) pair. Eligible donors are donors with
>= 4 cells of that cell type in that dataset. A setting is eligible if it has >= 200 eligible donors and the cell type is not 'unknown'.
Excluded: every dataset whose collection is one already used in this project (PsychAD, OneK1K/Yazar 2022, Perez 2022 lupus, SEA-AD), matched
by collection name or dataset title containing 'PsychAD', 'OneK1K', 'Yazar', 'lupus', 'SEA-AD' or 'Seattle Alzheimer' (case-insensitive);
excluded collections are listed in the output. Ranking: eligible settings by number of eligible donors, descending. At most ONE setting per
collection (its top-ranked one) so that settings are independent cohorts. The top 8 collections are analysed (fewer if fewer exist).
Per-donor cap: at most 200 cells per donor (seeded random subset, seed 91,000,000), to bound extraction time; stated as a deviation from
the earlier cohorts, which used all cells.

## 9.2 Extraction
Raw counts (X 'raw'); library size = Census raw_sum. Gene panel: the 400 genes with the highest mean counts per 10^4 over the selected cells
(the PsychAD/lupus rule without forced genes), about 79,800 pairs.

## 9.3 Analysis (per setting)
Donor-centred null exactly as Plan 7 mode A (Paule-Mandel-shrunken donor centring, cells permuted across all donors), 20 reassignments,
seeds 90,000,000 + 10k + setting index. Methods: PB, KX, XDS-M, XDS-Mh (Plan 6d definition). BH at q = 0.05.

## 9.4 Rules (fixed now; all results reported whatever they are)
- A method HOLDS in a setting if <= 2 of 20 runs have >= 1 BH call and mean per-pair type-I at 0.05 is <= 0.055 (Plans 6c, 6d, 7, 8).
- Predictions: pseudobulk fails in every setting; XDS-M and XDS-Mh each hold in every setting. The paper reports, for each method, the
  number of settings held, with per-setting runs-with-a-call, mean calls and type-I, plus cells per donor and number of donors.
- Cutoff: if results are not complete by 09:00 UTC on 5 Nov 2026 they are not added to the submitted paper; they still go into the
  repository.
