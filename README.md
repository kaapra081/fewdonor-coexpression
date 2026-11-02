# fewdonor-coexpression

Code, frozen analysis plans and result files for **"Why pseudobulk gene co-expression fails in large single-cell cohorts, and how to fix it"** (RECOMB 2027 submission; manuscript source in `paper2/`, branch `rev2-wip`).

In multi-donor single-cell RNA-seq, a gene pair can co-vary *within* donors (cell-state programs) or *between* donors (donor-level variation). The paper tests between-donor co-expression with
pseudobulk (PB, Pearson correlation of log-CPM donor means) and with several corrections, and evaluates them on **real-cell nulls**: nuclei are reassigned at random to pseudo-donors that keep the real
per-donor cell numbers, so the true between-donor covariance is exactly zero and every Benjamini-Hochberg (BH) call is a false call. Method names in the paper and in the code:

| Paper | Result key in the code | Where it is computed |
|---|---|---|
| PB, PBlin (Pearson correlation of log-CPM, respectively of raw ratio, donor means; Fisher z) | `PB`, `PBlin` | `pop/pop_core.py` (`run_methods`), `dstd/dstd.py` (`methods_d`); `dev/methods_v2.py` (`pearson_fisher`) |
| LCP (leak-corrected, self-standardized permutation test) | `LCP` | `pop/pop_core.py` (`run_methods`); earlier implementation `dev/methods_v2.py` (`lcsp_test`) and `lcp/core.py` |
| LCPd (leak-corrected, design-only standardization) | `LCPd` | `dstd/dstd.py` (`methods_d`) |
| XDS (cross-split, design-standardized donor-permutation test, normal-moment reference) | `LCPXd` (pooled-null summaries), `XDS` (REV code) | `dstd/dstd.py` (`methods_d`); `rev/revm.py` (`methods`) with `pop/pop_core.py` (`normal_perm_p`) |
| XDS-M (5 random splits, studentized across donors) | `XDSM` | `rev/revm.py` (`methods`) |
| KX (kappa extrapolation) | `KX` | `rev/revm.py` (`methods`) |
| LCPX, XU, LCP2 (reported negative results) | `LCPX`, `XU`, `LCP2` | `xsplit/`, `xu/`, `fix/` |

The repository also contains the earlier paper on the leak-corrected test (`paper/`, `lcp/`, `dev/`, `real/`, `rare/`, `fix/`); parts of it feed `paper2/` (SEA-AD nulls, depth-state coupling diagnostic) and are listed below.

## Install and check (no data needed; about 1 minute)

```bash
conda env create -f environment.yml && conda activate fewdonor-coexpression
pytest -q tests                      # 13 tests, about 20 s: LCP core, XDS/XDS-M/KX/BH on a synthetic null, make_nums.py layouts
python prereg/verify_hashes.py       # checks every frozen-plan SHA-256 against the file it names (see "Frozen plans")
python paper2/make_nums.py           # regenerates nums.tex and tables from the committed result files; prints what it could not regenerate
python paper2/make_nums.py --audit   # lists every input it reads: found / missing, file counts and sizes
python rev2/interim_to_summary.py && python paper2/make_rev2.py rev2/interim/REV2_SUMMARY_interim.json   # run-level table from the committed interim snapshot
```

`torch` is required (CPU is enough for the tests); a GPU only speeds up the permutation runs. The versions the tests were last run with are listed in `environment.yml`.

**What `make_nums.py` does and does not regenerate in a clean clone.** It reads committed JSON in `results/` (or the authors' working layout, `<dir>/` next to the scripts). The per-run arrays `*.npz` are gitignored and are not in the repository; their compact summaries
(`paper2/{onek,dc,dcg,seanull}_summary.json`) stand in for them. Where neither exists the script prints a WARNING that names the missing file, keeps the committed value of the affected macros (listed as STALE in `paper2/nums_provenance.json`) and does not rewrite the
affected table. On this branch that applies to `\DCgNull`, `\DCgTI`, `\DCgReal` and `table_ablate.tex` until `paper2/dcg_summary.json` and `paper2/seanull_summary.json` have been produced on a machine that has the `.npz` files (run `make_nums.py` there and commit those two files).
`python paper2/make_nums.py --strict` exits with status 2 if anything was skipped. `nums_provenance.json` records, for every macro, the block that set it, and lists the hand-typed constants (currently only `\NnucPsy`).

## Data sources

All data are public. None is stored in the repository (`*.h5ad` is gitignored). Dataset and version identifiers are CELLxGENE Discover identifiers; a version is downloaded as `https://datasets.cellxgene.cziscience.com/<version>.h5ad`.

| Data | Identifier | Used for | Fetch / prepare |
|---|---|---|---|
| PsychAD prefrontal cortex microglia and perivascular macrophages (1,469 donors, 10x 3' v3; raw counts in `raw.X`) | dataset `09642770-9a4c-4fe6-9f45-67e1249adc33`, version `44e42474-0297-463b-92b0-4b47decf699e`; local file name `pfc_myeloid_1469.h5ad` (sha256 `12911b2f61c54c5780a44c70aa8e70f7c57a2be4226d30061b454df6b1c586bc` as recorded by the authors; not re-verified for this README) | pooled and donor-centred nulls, real donors, differential co-expression; subtypes P1 Homeo_CECR2, P2 PVM_CD163, P3 ADAM_GPNMB (obs `subtype`; donors with at least 4 nuclei) | download as above; `python pop/pop_prep.py P1` (then P2, P3) extracts the 417-gene panel to `pop/pop_P*_X.npz`, `pop_P*_meta.npz`, `pop_P*_info.json` |
| SEA-AD middle temporal gyrus microglia-PVM | dataset `c76098ba-eed3-45b1-98f2-96fcac55ed18`, version `ab53c08f-03a3-4e0e-973a-44beabd78f60` (`data/checksums.sha256`) | real-cell nulls S1 = Micro-PVM_1 (79 donors), S2 = Micro-PVM_3-SEAAD (83 donors); depth-state coupling | `bash data/download.sh`, `python data/make_seaad_v3.py` (84-donor 10x v3 subset), `python rare/rare_prep.py` |
| OneK1K blood (Yazar et al. 2022) | CELLxGENE dataset `3faad104-2ab8-4434-816d-474d8d2641db`, queried through CELLxGENE Census release `2025-11-08`; cell types transitional stage B cell (B1), naive B cell (B2), naive thymus-derived CD4-positive alpha-beta T cell (B3) | replication of the null in 951 / 972 / 981 donors | `python onek/onek_fetch.py` and `onek/onek_fetch_b3.py` (400 genes); `python rev3/rev3_fetch.py` (2,000 genes, B1 and B2). Needs `cellxgene-census`, `tiledbsoma` and network access to the Census S3 bucket |
| Lee et al. 2026 fresh ex vivo myeloid | dataset `440da0a0-f321-4f40-9080-bb04031681e5`, version `21a0c6d5-658c-434e-a479-1e51504aec5c` (`data/checksums.sha256`) | depth-state coupling diagnostic (`\CouplingFracLee`) and the earlier paper | `bash data/download.sh`, `python data/p1_prepare_subset.py ...` |

Four scripts hard-code the authors' data path instead of reading `FEWDONOR_DATA`: `pop/pop_prep.py` (line 5), `dc/dc_run.py` (13), `rev/dcg_run.py` (13) and `rev2/rev2_D.py` (18); edit that path before running them. Several others default to the same absolute path when `FEWDONOR_DATA` is unset.

## Directory map

Scripts write next to themselves (`<dir>/` or `<dir>/out/`); the committed copies of their summary files are under `results/<dir>/` (for example `pop/pop_results.json` is committed as `results/pop/pop_results.json`). `paper2/make_nums.py` accepts both locations.

| Directory | Contents | Produces (committed under `results/`) | Feeds in `paper2/` |
|---|---|---|---|
| `lcp/`, `dev/` | LCP package; simulator `sim_v2.py`; test statistics `methods_v2.py` (`bh`, `paule_mandel`, `pearson_fisher`, LCP); development and held-out runners | `results/dev`, `results/heldout`, `results/ablate` (earlier paper) | `bh` for all BH calls |
| `pop/` | PsychAD pooled null (pseudo-donors, real cell numbers): `pop_prep.py`, `pop_run_d.py` (PB, LCP, LCPd, XDS, PBlin; D = 25 ... all), `pop_real.py` (real donors), `pop_summary.py` | `results/pop/pop_results.json`, `pop_P*_info.json`, `pop_summary.txt` | Sec. 5.2 (Table 2, Fig. 2a), Sec. 5.5; `\DPone ... \PBrepFalse`, `\RealXDSfrac` |
| `dstd/`, `xsplit/`, `xu/`, `fix/` | LCPd/PBlin; LCPX; XU and coupling stress; LCP2: SEA-AD nulls (`*_null.py`) and simulations (`*_sim.py`, `xu_lsim.py`), summaries (`*_summary.py`) | `results/dstd/d_results.json`, `results/xsplit/xs_results.json`, `results/xu/xu_results.json`, `results/fix/fix_results.json`; `results/xu/couple_results.json` (no committed generator) | Sec. 5.1 (Table 1), Sec. 5.3 (Table 3); `\LPBFDRa ... \SeaPBcalls`, `\Couple*` |
| `rev/` | `revm.py` (XDS, XDS-M, KX from cell lists); `rev_pop.py` (R1: 100 reassignments per subtype; `--real`: positive control), `rev_centred.py` (R6 donor-centred), `rev_sim.py` (R2 simulations), `sea_null.py` (SEA-AD XDS-M/KX), `dcg_run.py` (group-centred differential co-expression) | `results/rev/out/R1_*.json` (300), `R6_*.json` (9), `real_P*.json` (3), `results/rev/rev_sim_results.json`, `pooling.json`, `rev_dc_r5.json` (the last three have no committed generator) | Sec. 5.1-5.5, Table 2 and 3 columns for XDS-M and KX; `\Rone*`, `\Rsix*`, `\Pos*`, `\XDSM*`, `\DCg*` |
| `dc/`, `onek/` | AD vs control differential co-expression null (`dc_run.py`); OneK1K runs (`onek_run.py`, `onek_run2.py`) | `results/dc/dc_summary.json`, `results/onek/onek_B*_info.json` (summaries in `paper2/dc_summary.json`, `paper2/onek_summary.json`) | Sec. 5.4 (`\DCtext`), OneK1K paragraph (`\OnekText`, `\OnekAbs`) |
| `split/`, `diag/`, `real/`, `rare/` | split-seed dependence of XDS; post hoc diagnostics; earlier real-data analyses and rare-supertype statistics | `results/split`, `results/diag`, `results/real`, `results/rare` | Appendix (split seeds), `\Coupling*`, `\PBposFrac` |
| `rev2/` | REV2 jobs (`jobs.txt`, `run_all.sh`): `rev2_A.py` (donor-centred run-level null at D = 200, 400, all), `rev2_C.py`, `rev2_D.py`, `rev2_D2.py`, `rev2_F.py`, OneK1K `onek_run2.py`; `rev2_summary.py` -> `rev2/out/REV2_SUMMARY.json` | `rev2/interim/A_interim_0615UTC.txt` (06:15 UTC snapshot: P2 = 40 runs, P3 = 100, P1 absent) | `paper2/make_rev2.py` -> `table_runlevel.tex`, `rev2nums.tex`, `rev2_macros.json` |
| `rev3/` | REV3 jobs: `rev3_A2.py` (tail calibration by gene-expression stratum), `rev3_A3B3.py` (small-D run level; clipping sensitivity), `rev3_B2.py` (adjusted and weighted pseudobulk), `rev3_fetch.py` | none committed yet (outputs go to `rev3/out/*.json`) | - |
| `paper2/` | manuscript (`main.tex`, `appendix.tex`, `abstract_body.tex`, `refs.bib`), generated `nums.tex`, `table_*.tex`, `*_summary.json`, `nums_provenance.json`; generators `make_nums.py`, `make_rev2.py`, `make_fig.py` | | |
| `rev4/` | REV4 jobs: `rev4.py --mode R1|R2|R3|R3o|R4` (small-D null with 2e6 permutations; brain-bank-stratified null; PC-adjusted and weighted pseudobulk; OneK1K pseudobulk variants; real-donor consequences), `rev4_R2c.py` (post hoc bank-centred diagnostic), `rev4_R6.py`, `rev4_R6sum.py`, `rev4_R6pc1.py`, `rev4_R6hs.py`, `rev4_R6clusters.py` (Mantis map tables and statistics), `rev4_C.py` (Amendment 3: C1 covariates, C2 bank-mean-removed null, C3 overlap baselines, C4 permutation XDS at 200 donors, C5 coupling, C6 map robustness), `rev4_A4.py` (Amendment 4: function modularity), `rev4_A5.py` + `rev4_A5local.py` (Amendment 5: stratified/20-run/equal-size/GO modularity, sex positive control); `mantis/` holds the uploaded map CSVs, the Mantis read-back (`gene_map_mantis.csv`, `gene_clusters.csv`) and `mantis_space_manifest.json` | `rev4/out/*.json`, `rev4/out/R6/`, `rev4/results/REV4_SUMMARY.json` | Sec. 5.1-5.3, Fig. 3, Table 2, App. F |
| `prereg/`, `rev2/`, `rev3/` | frozen plans, their `.sha256` files, decision records; `prereg/verify_hashes.py` | | Appendix E |
| `tests/` | `test_core.py` (LCP), `test_xds_smoke.py` (XDS/XDS-M/KX/BH on a synthetic complete null), `test_make_nums.py` (layouts, fallbacks, round trip) | | |
| `paper/` | the earlier LCP paper and its generators (`make_figs.py`, `make_results.py`) | | |

Not reproducible from this repository alone: `paper2/make_fig.py` needs the `apply_figure_style` and `panel_letter` helpers of the authors' plotting environment and defaults `RV2` to the authors' machine; `paper2/figures/fig_schematic.png` has no generator;
`main.tex` includes `figures/fig_growth.pdf` and `figures/fig_schematic.pdf`, which must be committed next to the PNGs (`.gitignore` now allows `paper2/figures/*.pdf`) before the manuscript compiles from a clean clone.

## Frozen plans

Each plan is a Markdown file with a SHA-256 file next to it (`python prereg/verify_hashes.py` recomputes them; as run on this branch every plan hash matches the file it names except the two superseded versions noted). The hashes are self-generated; they show that a plan was not edited after hashing, not when it was written. The only external date is the git commit time (author-set), and for several plans the plan and its results first appear in the same commit, so git does not show that the plan came first.

| Plan | File and hash (first 12 hex) | First commit (Oct 2026, PDT) | Results first committed | Note |
|---|---|---|---|---|
| DEV (development gates) | `prereg/DEV_PREREG.md`; `DEV_PREREG_change1.sha256` a860bfa48ce7 matches; `DEV_PREREG.sha256` 010de04311d3 hashes the version before change 1 | d00b9c2 (26, 11:37); change 1 eee7366 (26, 13:06) | d258bc7 (26, 22:44) | code hashes in both files partly differ from the current `dev/*.py` (edited after the plan) |
| HELDOUT_SPEC | `prereg/HELDOUT_SPEC.md`, 6c7b54122e2b (match; code hashes match) | c986bab (26, 21:12) | 0b728fb (27, 19:28) | |
| RARE | `prereg/RARE_PREREG.md`, 4e389922d98e | 38765b5 (28, 00:01) | 5afb31b (28, 00:41) | |
| FIX | `prereg/FIX_PREREG.md`, 3b76d65ccf0a | 5f0a800 (28, 10:21) | be75558 (28, 15:49) | |
| XS | `prereg/XS_PREREG.md`, 6022eabd47aa (amended); `XS_PREREG_original.sha256` 247328aedadd is the version before Amendment 1 | 22671cd (29, 00:04) | same commit | LCPX is a negative result |
| DSTD | `prereg/DSTD_PREREG.md`, 485fcfbf1f43 | 22671cd (29, 00:04) | same commit | states the disclosed exploratory diagnostic behind `\Rcorrln`, `\Mcorrln` |
| POP | `prereg/POP_PREREG.md`, 083e48548067 (the `.sha256` names it `pop/POP_PREREG.md`, so `sha256sum -c` from the root fails; use `verify_hashes.py`) | 22671cd (29, 00:04) | same commit | |
| XU, COUPLE, DC | `prereg/{XU,COUPLE,DC}_PREREG.md`, 15d353ed3df6, 7c5aae02a193, 4d086ce2f171 | da40d8a (29, 11:29) | same commit | DC states the nucleus-count medians used in `\DCtext` |
| REV, ONEK | `prereg/REV_PREREG.md`, `prereg/ONEK_PREREG.md`, b79c1ae0bf2f, cbba3c0d5704 | b38f509 (29, 14:20) | same commit | |
| REV2 | `rev2/REV2_PREREG_recovered.md`: a transcript-recovered copy that **does not reproduce the frozen hash**; the canonical file (sha256 beginning `1a7e4de1`, per the recovered copy's header) is not in the repository | recovered copy 84308ee (31, 12:35) | interim snapshot 04b4610 (31, 12:35) | the plan's date line reads 29 Oct 14:30 |
| REV3 | `rev3/REV3_PREREG.md`, 314e01284fb9 | 0e8314a (31, 13:04); scripts 351af77 (31, 13:11) | none yet | written after the REV2 interim results were seen (stated in the plan) |
| REV4 | `rev4/REV4_PREREG.md`, 4c7d552792 (Amendment 1 92c885d008, Amendment 2 07edeeedb9, post hoc; Amendment 3 4363134a; Amendments 4 76f1a8d5 and 5 a501dc2b, post hoc) | 74ba4e5, 9a17b84, c378efc (1 Nov) | `rev4/results/REV4_SUMMARY.json` | written after an external AI review; Amendment 2 written after partial R2 results (stated in it) |

Seeds are fixed in the scripts: development 2,000,000 + 100,003 k + rep, held-out 5,000,000 + 100,019 k + rep, simulation gene panel 4112026; `rev/rev_sim.py` 9,500,000 (L designs), 9,700,000 (LS), 9,000,000 (F); `rev/rev_pop.py` R1 48,000,000 + 10 k + subtype index (real donors 44,000,000);
`rev/rev_centred.py` 49,000,000 + 10 k + i; `dc/dc_run.py` and `rev/dcg_run.py` 41,000,000 + 10 k + i (real 43,000,000 + i); `onek/onek_run.py` 45,000,000 + 10 k + cell-type index (real 47,000,000 + index); `rev/sea_null.py` 31,000,000 + k (+1000 for S2); REV2-A 50,000,000 + 10 k + i; REV3 A3 60,000,000, B2 61,000,000, B3 62,000,000, A2 70,000,000 + ... (see plan).

## Regenerating the results (needs the raw data; not executed for this README)

Only the commands in "Install and check" were run when this README was written. The commands below were assembled from the scripts' argument parsers and the frozen plans (reassignment indices `k` and seeds are those stated in the plans);
they need the data above, hours of CPU or GPU time, and they write into the script directories. Run them from the repository root.

```bash
export FEWDONOR_DATA=data/                                   # where the .h5ad files are
# 1. panels
for s in P1 P2 P3; do python pop/pop_prep.py $s; done        # PsychAD 417-gene panel (edit DATA in the script first)
python onek/onek_fetch.py && python onek/onek_fetch_b3.py    # OneK1K, Census 2025-11-08
python rare/rare_prep.py                                     # SEA-AD S1, S2
# 2. pooled PsychAD null (POP Part A; D in 25, 50, 100, 200, 400, 800, all) and real donors (Part B)   -> Sec. 5.2, 5.5
for s in P1 P2; do for k in 0 1 2; do python pop/pop_run_d.py --sub $s --k $k --halves --Ds all,800,400,200,100,50,25; done; done
for k in 1 2 3; do python pop/pop_run_d.py --sub P3 --k $k --halves --Ds all,800,400,200,100,50,25; done   # P3 k = 0 is the exploratory reassignment excluded by DSTD_PREREG
for s in P1 P2 P3; do python pop/pop_real.py --sub $s; done
python pop/pop_summary.py                                    # -> pop/pop_results.json; commit as results/pop/pop_results.json
# 3. run-level null R1, donor-centred null R6, positive control (REV_PREREG)
for s in P1 P2 P3; do python rev/rev_pop.py --sub $s --k0 0 --k1 100; python rev/rev_pop.py --sub $s --real; for k in 0 1 2; do python rev/rev_centred.py --sub $s --k $k; done; done
# 4. SEA-AD real-cell nulls (S1, S2) and simulations   -> Sec. 5.1, 5.3
for st in S1 S2; do
  for k in 0 1 2 3 4;   do python fix/fix_null.py --set $st --k $k; python dstd/d_null.py --set $st --k $k; python xsplit/xs_null.py --set $st --k $k; done
  for k in 5 6 7 8 9;   do python xu/xu_null.py --set $st --k $k; done
  for k in 10 11 12 13 14; do python rev/sea_null.py --set $st --k $k; done
done
for d in F_D16 F_D24 F_D40 F_S1 F_S2 F_D16_n100; do python fix/fix_sim.py --design $d; python xsplit/xs_sim.py --design $d; python dstd/d_sim.py --design $d; done
for d in L200 L500 L1000 LS200 LS500 LS1000 F_D16 F_D24 F_D40 F_S1 F_S2 F_D16_n100; do python rev/rev_sim.py --design $d; done
for d in L200 L500 L1000 LS200 LS500 LS1000 F_D16 F_D24 F_D40 F_S1 F_S2 F_D16_n100; do python xu/xu_lsim.py --design $d; done
for s in P1 P2; do for k in 3 4 5; do python xu/xu_pop.py --sub $s --k $k; done; done; for k in 4 5 6; do python xu/xu_pop.py --sub P3 --k $k; done
python fix/fix_summary.py; python dstd/d_summary.py; python xsplit/xs_summary.py; python xu/xu_summary.py
python split/split_stab.py                                   # split-seed dependence (appendix)
# 5. differential co-expression (DC_PREREG: k = 10, 11, 12; real donors k = -1) and OneK1K (ONEK_PREREG: K = 3)   -> Sec. 5.4
for s in P1 P2 P3; do for k in 10 11 12 -1; do python dc/dc_run.py --sub $s --k $k; python rev/dcg_run.py --sub $s --k $k; done; done
for t in B1 B2 B3; do for k in 0 1 2 -1; do python onek/onek_run.py --typ $t --k $k; done; done
# 6. REV2 (jobs.txt expects the layout of the upload bundle: copy run_all.sh and jobs.txt to the repository root first) and REV3
cp rev2/run_all.sh rev2/jobs.txt . && NPAR=3 bash run_all.sh && python rev2/rev2_summary.py        # -> rev2/out/REV2_SUMMARY.json
python paper2/make_rev2.py rev2/out/REV2_SUMMARY.json        # pass the file explicitly: the script's default name is lower case
python rev3/rev3_fetch.py                                    # 2,000-gene OneK1K panels for REV3 A2
python rev3/rev3_A2.py --ds B1 --null 0 --k0 0 --k1 20      # also --ds B2 and --null 1; A3/B3: rev3_A3B3.py --mode A3|B3 --sub P1 ...; B2: rev3_B2.py --sub P1
# 6b. REV4 (GPU recommended for R1: PERM_DEVICE=cuda PERM_BATCH=2000)
python rev4/rev4.py --mode R1 --sub P1 --k0 0 --k1 20   # also P2, P3; --mode R2 (k 0-40), R3 (k 0-20), R3o --sub B1|B2|B3 (k 0-10), R4
python rev4/rev4_R2c.py --sub P1                         # post hoc bank-centred diagnostic (Amendment 2)
for s in P1 P2 P3; do python rev4/rev4_C.py --mode C1 --sub $s; done   # Amendment 3 (also C1real, C2, C3, C4, C5, C6)
for s in P1 P2 P3; do python rev4/rev4_R6.py --sub $s; python rev4/rev4_A5.py --mode top --sub $s; python rev4/rev4_A5.py --mode sexctl --sub $s; python rev4/rev4_A5.py --mode null20 --sub $s; done
python rev4/rev4_A4.py; python rev4/rev4_A5local.py      # Amendments 4-5 (need the compact A4_<s>.npz extracted from R6_<s>.npz)
python rev4/rev4_R6.py --sub P1 && python rev4/rev4_R6sum.py && python rev4/rev4_R6pca.py && python rev4/rev4_R6pc1.py && python rev4/rev4_R6hs.py && python rev4/rev4_R6clusters.py
python paper2/make_rev4.py                               # rev4nums.tex, table_pbfix/calib/R1/R2.tex
# 7. diagnostics, paper numbers, tables, figure
python real/real_prep.py && python diag/diag_real.py         # depth-state coupling (results/diag/diag_real.json)
python paper2/make_nums.py                                   # writes dcg_summary.json and seanull_summary.json when the .npz are present: commit them
python paper2/make_fig.py                                    # needs the authors' figure-style helpers (see above)
cd paper2 && tectonic main.tex                               # or pdflatex/bibtex; needs figures/*.pdf
```

## Known gaps

What cannot be regenerated from a clean clone of this repository, and why:

- **Raw per-run arrays** (`*.npz`, gitignored): `rev/out/DCg_*`, `rev/out/seanull_*`, `onek/out/ON_*`, `dc/out/DC_*` and every other `*/out/*.npz`. Compact summaries replace them for `onek`, `dc` and `r1`; `paper2/dcg_summary.json` and `paper2/seanull_summary.json` do not exist yet, so `\DCgNull`, `\DCgTI`, `\DCgReal` and `table_ablate.tex` keep their committed values (flagged by `make_nums.py`).
- **Missing generators:** `results/rev/pooling.json`, `results/rev/rev_dc_r5.json`, `results/rev/rev_sim_results.json` and `results/xu/couple_results.json` are committed, but no script in the repository writes them.
- **REV2 and REV3 outputs** (`rev2/out/*.json`, `onek/out/E_*.json`, `rev3/out/*.json`, `REV2_SUMMARY.json`) are not committed; only the 06:15 UTC interim snapshot is (`rev2/interim/`).
- **Inputs derived from the public data** (`pop/pop_P*_X.npz`, `pop_P*_meta.npz`, `onek_B*_X/meta.npz`, `rare_stats_*.npz`, `real_stats_*.npz`) are not committed; `pop_prep.py`, `onek_fetch*.py`, `rare_prep.py` and `real_prep.py` rebuild them from the public files.
- **Hand-typed constant:** `\NnucPsy` (the number of nuclei in the PsychAD file) is not recorded in any committed file.
- **REV2 plan:** the canonical `REV2_PREREG.md` is not in the repository (see "Frozen plans").
- **Code defects found while writing this README:** `rev2/rev2_C.py` reads `sign_XDSM`, which `rev/revm.py` (`methods`) does not return; `rev2/run_all.sh` assumes it was copied to the repository root; `paper2/make_rev2.py` defaults to `rev2_summary.json` while `rev2_summary.py` writes `REV2_SUMMARY.json`; four scripts hard-code a data path (see "Data sources").

## Earlier analyses (LCP paper, `paper/`)

These produce the SEA-AD and Lee 2026 inputs used above and the earlier manuscript.

| Step | Command |
|---|---|
| Development simulations (R = 200 per design) | `python dev/run_dev.py --design D16 --out runs_dev` (designs D8, D12, D16, D24, D40, D16_strong, D16_STRESS_depth), then `python dev/run_ci.py ...` and `python dev/analyze_dev.py --runs runs_dev --out results/dev` |
| Held-out evaluation | `python dev/run_heldout.py --design HO_seaad24 --out runs_heldout` (6 designs; `prereg/HELDOUT_SPEC.md`) |
| Real-data statistics | `FEWDONOR_DATA=data/ python real/real_prep.py`; split-half: `real/real_prep_split.py`, `python real/real_gpu.py --part split`; tests: `python real/real_gpu.py --part full` (also `cap`, `sub24`, `sub40`), `python real/real_summarize.py` |
| Post hoc ablations | `python dev/ablate_gpu.py --design D16 --out runs_ablate --device cuda:0` (or `cpu`; D16, D24, D40), `python dev/analyze_ablate.py --runs runs_ablate --out results/ablate` |
| Post hoc diagnostics | `python diag/diag_sim.py`, `FEWDONOR_DATA=data/ python diag/diag_real.py`, `python diag/diag_null.py` (outputs in `results/diag/`) |
| Rare-supertype analysis | `python rare/rare_run.py`, `python rare/rare_sim.py --set S1`, `python rare/rare_sim_summary.py`, `python rare/rare_descriptive.py` (`prereg/RARE_RESULTS.md`) |
| Figures and numbers of the earlier paper | `python paper/make_figs.py && python paper/make_results.py`, then compile `paper/main.tex` |

The earlier paper's statements about LCP (for example calibration "from about 24 donors") concern its simulations. For real-cell nulls at hundreds of donors the paper in `paper2/` reports different outcomes for LCP; use `results/pop/pop_results.json` and `results/fix/fix_results.json` as the record. `dev/empirical_designs.npz` holds only per-cell library sizes of the real donors.

## License

MIT (see `LICENSE`). No third-party method code is redistributed.
