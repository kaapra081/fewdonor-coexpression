# Population-scale analysis: frozen plan (28 Oct 2026, 17:42 PDT; only cell metadata (obs) of the dataset has been inspected)

## Data
PsychAD frozen-PFC myeloid nuclei, "Frozen in silico sorted myeloid single nuclei RNA-seq from human prefrontal cortex",
CELLxGENE dataset 09642770-9a4c-4fe6-9f45-67e1249adc33, version 44e42474-0297-463b-92b0-4b47decf699e
(file sha256 12911b2f61c54c5780a44c70aa8e70f7c57a2be4226d30061b454df6b1c586bc). 289,279 nuclei, 1,469 donors, 10x 3' v3.
Raw counts (raw.X if present, otherwise X if integer-valued). Library size = total raw counts of the nucleus over all genes.

Subtypes (column `subtype`), chosen from metadata to span cells per donor:
- P1 = Homeo_CECR2 (56,911 nuclei; ~39 per donor), P2 = PVM_CD163 (37,365; ~25), P3 = ADAM_GPNMB (16,966; ~12).
- Donors with >= 4 nuclei of the subtype.
Genes per subtype: the 400 genes with the highest mean counts-per-10k in that subtype, plus the 17 SEA-AD control genes when
present and expressed (as RARE_PREREG). Selection uses means only.

## Methods
PB (Pearson correlation of donor log-CPM means, t-test), LCP (as in the paper), LCPX and LCPX-S (XS_PREREG.md). LCP2 omitted (cost).
p-values: B = 10^5 donor permutations when D <= 100. When D >= 200, the permutation p-value is replaced by the two-sided normal
approximation with the exact permutation mean and variance of the statistic (Hoeffding's combinatorial CLT); for LCPX-S, by N(0,1).
BH at q = 0.05 over all pairs.

## Part A. Real-cell null with exact truth, as a function of cohort size
For each subtype: pool its nuclei over the included donors and permute donor labels across nuclei (n_d kept), K = 3
reassignments (seeds 41,000,000 + 10 k + subtype index). In each reassignment draw nested random pseudo-donor subsets of
size D in {25, 50, 100, 200, 400, 800, all}. Every call is false. Report, per method and D: mean BH calls, per-pair type-I at
0.05 (all pairs; top-1% within-correlated pairs).
**Prediction (stated before running):** PB false calls increase with D (fixed leak bias, shrinking standard error); LCPX
false calls stay at or near 0 and its type-I near 0.05 at every D, including top-1% pairs.
**Decision "pseudobulk false discoveries grow with cohort size" = YES** for a subtype if PB mean calls at D = all exceed
those at D = 100 by a factor >= 3 and are >= 100.
**Decision "LCPX holds the real-cell null at population scale" = YES** for a subtype if at every D LCPX all-pair type-I
<= 0.055 and top-1% type-I <= 0.05 + 3 SE, and mean BH calls <= 1.

## Part B. Real donors
Full D per subtype: BH calls of PB, LCP, LCPX; overlaps; for PB-only calls (PB yes, LCPX no) vs shared calls, the median
absolute within-donor correlation.
Replication: split donors at random into two disjoint halves (seed 42,000,000); a call in half 1 replicates if the same pair has
p < 0.05 with the same sign in half 2. Report replication rates of PB-only and of shared half-1 calls, and the same quantity in
the Part A null at matched D (pseudo-donor halves).
**Decision "replication does not protect pseudobulk against the leak" = YES** if, in the Part A null at D = half, PB's
replicated false calls are >= 50 for a subtype.
Results are reported whatever they are.
