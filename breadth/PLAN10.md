# Plan 10 (5 Nov 2026, ~06:50 UTC). POST HOC, larger breadth check, written after Plan 9's results were seen and BEFORE any Plan 10 data were read.

Same code, null and rules as Plan 9 (breadth/PLAN9.md, 9.2-9.4) with these changes:
- Eligible settings: (dataset, cell type) with >= 100 donors having >= 4 cells (Plan 9: 200); cell type not 'unknown'.
- Excluded collections: as Plan 9, plus PsychAD's collection by id (84ce6837-548d-4a1f-919f-0bc0d9a3952f), which Plan 9's keywords missed.
- Up to 3 cell types per collection (the 3 with most eligible donors); ranked by eligible donors; at most 24 settings.
- No per-donor cell cap: all cells of eligible donors are used.
- Seeds 93,000,000 + 10k + setting index.
- Rule per setting as Plan 9: a method holds if <= 2 of 20 runs have >= 1 BH call and mean per-pair type-I <= 0.055.
- Settings run on two machines (Modal, Colab) from opposite ends of the ranking; a setting computed on both must agree exactly
  in call counts (reported). Settings not complete by 09:45 UTC on 5 Nov 2026 are reported as not run (no other selection).
