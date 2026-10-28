# fewdonor-coexpression

Code for **"Between-donor gene co-expression from few donors: a leak-corrected permutation test and a bound on what any test can detect"** (RECOMB 2027 submission).

In multi-donor single-cell RNA-seq, two genes can co-vary *within* donors (cell-state programs) or *between* donors (donor-level variation). This repository provides:

- **LCP** (`lcp/`), a leak-corrected, precision-standardized donor-permutation test of between-donor co-expression on UMI counts. It removes the finite-cell "leak" of within-donor covariance into donor means (pseudobulk covariance targets Σ_B + κ̄ Σ_W).
- The simulator with exact analytic truth, the pre-registered development and held-out evaluations (`dev/`, `prereg/`), the real-data analyses of SEA-AD and Lee et al. 2026 microglia (`real/`), and scripts that regenerate every figure, table and number in the paper (`paper/`).

**Scope and limits (please read).**
- LCP controls the FDR in simulation from **~24 donors** upward, where its power matches a noise-limited oracle benchmark.
- At **≤16 donors** no method based on observed donor means that we tested, LCP included, is calibrated in the tail used by Benjamini–Hochberg.
- LCP tests *independence* of the two genes' donor effects. If donor effects are dependent but uncorrelated, it is not calibrated.
- Covariance intervals are not provided: they undercover at these donor numbers even when computed from the true donor means.

## Install

```bash
conda env create -f environment.yml && conda activate fewdonor-coexpression
```

numpy, scipy, pandas, matplotlib and anndata are required. torch is optional and is only used for GPU permutations in `real/real_gpu.py`.

## Use

```python
from lcp import donor_statistics, lcp_test, bh
st  = donor_statistics(Y, s, donor)   # Y: cells x genes UMI counts (dense or sparse), s: total UMI per cell, donor: labels
res = lcp_test(st, B=20000)           # all gene pairs j<k
calls = bh(res["p"])                  # BH at q = 0.05
res["est_bces"], res["leak"]          # leak-corrected between-donor covariance estimate, and the leak
```

`python examples/quickstart.py` runs the method on one simulated 40-donor data set.

## Reproduce

| Step | Command |
|---|---|
| Tests | `pytest -q tests` |
| Development simulations (R = 200 per design) | `python dev/run_dev.py --design D16 --out runs_dev` (designs D8, D12, D16, D24, D40, D16_strong, D16_STRESS_depth), then `python dev/run_ci.py ...` for intervals and `python dev/analyze_dev.py --runs runs_dev --out results/dev` |
| Held-out evaluation | `python dev/run_heldout.py --design HO_seaad24 --out runs_heldout` (6 designs; see `prereg/HELDOUT_SPEC.md`) |
| Data | `bash data/download.sh` (checksums verified), then `python data/make_seaad_v3.py` and `python data/p1_prepare_subset.py --h5ad data/lee2026_fresh_myeloid.h5ad --out data/P1_lee2026_microglia_80x1600.h5ad` |
| Real-data statistics | `FEWDONOR_DATA=data/ python real/real_prep.py` |
| Split-half leak experiment | `FEWDONOR_DATA=data/ python real/real_prep_split.py`, then `python real/real_gpu.py --part split` |
| Real-data tests | `python real/real_gpu.py --part full` (also `cap`, `sub24`, `sub40`; `--device cpu` works but is slower), then `python real/real_summarize.py` |
| Post hoc ablations (not pre-registered) | `python dev/ablate_gpu.py --design D16 --out runs_ablate --device cuda:0` (or `--device cpu`; designs D16, D24, D40), then `python dev/analyze_ablate.py --runs runs_ablate --out results/ablate` |
| Figures, tables, numbers | `python paper/make_figs.py && python paper/make_results.py`, then compile `paper/main.tex` |

Seeds are fixed in the scripts. Development seeds are 2,000,000 + 100,003·k + rep; held-out seeds are 5,000,000 + 100,019·k + rep, with gene panel seed 4112026.

## Pre-registration

`prereg/` holds the development gates (`DEV_PREREG.md`), the frozen held-out specification and their SHA-256 hashes, each committed before the corresponding runs. It also holds the decision record (`GO_NOGO.md`). The decision record documents two criteria amended before the affected results were read, and reports that the development gates for validity at ≤16 donors and for intervals were **not** met.

## Data

Both datasets are public CELLxGENE datasets:
- SEA-AD MTG Microglia-PVM, dataset `c76098ba-eed3-45b1-98f2-96fcac55ed18`, version `ab53c08f-03a3-4e0e-973a-44beabd78f60` (Gabitto et al. 2024, doi:10.1038/s41593-024-01774-5).
- Lee et al. 2026 fresh ex vivo myeloid, dataset `440da0a0-f321-4f40-9080-bb04031681e5`, version `21a0c6d5-658c-434e-a479-1e51504aec5c` (doi:10.1038/s41588-026-02716-6).

`dev/empirical_designs.npz` contains only per-cell library sizes of the real donors. It is used to give held-out simulations real cell numbers and depths.

## License

MIT (see `LICENSE`). No third-party method code is redistributed.
