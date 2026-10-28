# Outcomes of FIX_PREREG.md (frozen at repo commit 5f0a800, sha256 3b76d65c…)

## Item 1: real-cell pseudo-donor null (K = 5 reassignments per supertype)

| Supertype | Method | Mean BH calls (all false) | Range | Per-pair type-I (α = 0.05) | Type-I, top 1% within-correlated pairs |
|---|---|---|---|---|---|
| S1 | PB | 3713.8 | 2976–4601 | 0.176 | 0.674 |
| S1 | LCP | 0.0 | 0–0 | 0.081 | 0.092 |
| S1 | LCP2 | 0.0 | 0–0 | 0.046 | 0.068 |
| S2 | PB | 2713.8 | 1630–3669 | 0.174 | 0.607 |
| S2 | LCP | 0.0 | 0–0 | 0.061 | 0.099 |
| S2 | LCP2 | 0.0 | 0–0 | 0.050 | 0.085 |

Decisions:
- **S1:** leak produces spurious PB discoveries in real cells = **YES**. LCP holds the null = False; LCP2 holds the null = False.
- **S2:** leak produces spurious PB discoveries in real cells = **YES**. LCP holds the null = False; LCP2 holds the null = False.
- "Holds" fails on the per-pair type-I criteria. No method other than PB made any BH discovery.

## Item 2: LCP2 (R = 100 per design; new seeds; held-out panel)

| Design | LCP FDR (MCSE) | LCP2 FDR (MCSE) | PB FDR | LCP power | LCP2 power |
|---|---|---|---|---|---|
| F_D16 | 0.110 (0.025) | 0.052 (0.018) | 0.216 | 0.014 | 0.009 |
| F_D16_n100 | 0.039 (0.008) | 0.039 (0.009) | 0.122 | 0.087 | 0.085 |
| F_D24 | 0.050 (0.007) | 0.033 (0.006) | 0.174 | 0.122 | 0.106 |
| F_D40 | 0.051 (0.004) | 0.041 (0.003) | 0.133 | 0.501 | 0.482 |
| F_S1 | 0.077 (0.004) | 0.025 (0.003) | 0.594 | 0.399 | 0.298 |
| F_S2 | 0.051 (0.002) | 0.042 (0.002) | 0.179 | 0.929 | 0.926 |

Criteria: (a) True, (b) True, (c) False.

**Adopt LCP2 = False.** Criterion (c) fails: at 24 donors LCP2 has 0.87 of LCP's power, below the required 0.9. LCP2 is reported as a calibrated but lower-powered variant; recommending it for few-cell designs is post hoc.

## Deviations

- None from the plan.
- The pooling caveat is a descriptive addition: within-pseudo-donor correlation is somewhat stronger than within real donors; see null_caveat.json.
