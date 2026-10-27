# Held-out evaluation specification (recomb_v2) — frozen before any held-out replicate is run

Written Mon 26 Oct 2026, ~21:40 PDT. sha256 in HELDOUT_SPEC.sha256 and pushed to the private repo before launch.

## Separation from development
- Development (DEV_PREREG.md): panel seed 26102026; synthetic designs D8–D40; seeds 2,000,000 + 100,003·k + rep.
- Held-out: **new panel** (seed 4112026: new μ_j, φ_j, η_j draws), **new designs** (four of six use real per-cell library sizes of real donors from SEA-AD v3 or Lee 2026, round3/empirical_designs.npz; one uses a weaker effect size), **new seeds** 5,000,000 + 100,019·k + rep. Smoke tests used seeds 7.1e6 and were not inspected beyond runtime.
- The method is the frozen development method LCP (DEV_PREREG.md §3). Two implementation-only changes since the development freeze: (i) sim_v2.simulate gained an optional real-donor library-size pool (no effect on development designs); (ii) the refinement step of lcsp_test was vectorised; on development replicate D16/rep0000 it reproduces the stored LCP and ORACLE p-values exactly (max |Δp| 3·10⁻⁸, float32 storage).

## Designs (R = 200 each)
| ID | D | Cells per donor | Depths | Other |
|---|---|---|---|---|
| HO_seaad24 | 24 random SEA-AD donors | real (18–1,068) | real | — |
| HO_seaad24_cap40 | 24 SEA-AD donors | LogNormal(median 40, 0.66), min 10, cells drawn from the donor's real depths | real | rare-cell analogue |
| HO_lee40 | 40 Lee 2026 donors | LogNormal(median 60, 0.8), min 10 | real (Lee) | — |
| HO_D32 | 32 | LogNormal(40, 0.66) | fitted | — |
| HO_weak40 | 40 | LogNormal(40, 0.66) | fitted | g = g_Q = 0.20 (weak) |
| HO_D16_n100 | 16 | LogNormal(100, 0.66) | fitted | **few-donor check** |

## Pre-specified reading
- **Validity claim (D ≥ 24 designs: all except HO_D16_n100).** LCP FDR ≤ 0.05 + 2·MCSE and type-I at 0.05 ≤ 0.05 + 3·MCSE in WW, SS, HH, NN. A failure in any of these designs is reported as a failure of the paper's validity claim.
- **Power claim.** LCP BH power ≥ 0.8 × BOUND power in each D ≥ 24 design; reported with PB_logcpm, E5, CMB_naive, CMB_mom and ORACLE.
- **Few-donor check (HO_D16_n100).** Reported as is; the paper's claim is that calibration is not guaranteed at ≤ 16 donors.
- **Weak null (QQ_cross)** and **intervals** (BOOT_T vs ORACLE_BOOT_T): reported, no claim.
