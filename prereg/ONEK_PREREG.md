# Second tissue: OneK1K peripheral blood. Frozen plan (29 Oct 2026, 11:50 PDT; no OneK1K data or metadata inspected yet)
Data: CELLxGENE dataset 3faad104-2ab8-4434-816d-474d8d2641db ("Single-cell eQTL mapping identifies cell type specific genetic control
of autoimmune disease"), queried through the CELLxGENE Census (latest stable release; version recorded). Raw counts; library size =
total counts of the cell over all genes in the dataset.
Cell types (column cell_type), chosen from per-donor cell counts only: among cell types with >= 500 donors having >= 4 cells, the three
whose median cells per donor (over those donors) is closest to 20, 60 and 200 (one type per target, no type used twice).
Genes: the 400 with the highest mean counts-per-10k in that cell type (no forced controls).
Null: K = 3 donor-label reassignments (seeds 45,000,000 + 10 k + type index), real cell counts kept; nested subsets D = 200 and all.
Methods: PB, LCP, XDS exactly as in the paper (normal approximation with exact permutation moments for D >= 200).
Decisions (per cell type): "pseudobulk false discoveries grow with cohort size" = YES if PB mean calls at D = all >= 2 x at D = 200 and
>= 100; "XDS holds" = YES if at both D: mean calls <= 1, all-pair type-I <= 0.055, top-1% type-I <= 0.05 + 3 SE.
Real donors: XDS and PB calls, PB-only count (descriptive).

## Amendment 1 (29 Oct 13:20 PDT, before any OneK1K test was computed; only the B1/B2 fetches had finished)
(i) Genes are ranked on a 5,000-cell random subsample per cell type (seed 46,000,000) rather than all cells, to limit transfer.
(ii) XDS-M (REV_PREREG R2) is added, and a fourth decision: "XDS-M holds" = YES if at both D: mean calls <= 1, all-pair type-I <= 0.055,
top-1% type-I <= 0.05 + 3 SE. Run-level false-call frequencies over K = 3 runs are reported but are too few for a decision.
