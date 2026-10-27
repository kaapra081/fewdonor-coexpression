# method_v2: between-donor co-expression with few donors

Draft, Mon 26 Oct 2026. Status of each result is marked: **[proved]** (proof below), **[sketch]** (proof outline; conditions stated, not fully verified), **[empirical]** (simulation only), **[prior]** (known; cited).

## 1. Setting and assumptions

Donors d = 1..D; cells i = 1..n_d of one cell type; genes j. Latent relative expression θ_i ∈ R^G_+ (expected UMI per unit library size). Equal donor weights.

- **A1 (donors).** The donor-level cell-state distributions F_d are drawn independently from a donor population P. The design (n_d and the library sizes s_d1..s_dn_d) is ancillary: independent of the F_d and θ (it may differ arbitrarily between donors).
- **A2 (cells).** Given F_d and the design, θ_d1, …, θ_dn_d are iid draws from F_d.
- **A3 (counts).** Y_ij | θ_i, s_i ~ Poisson(s_i θ_ij), independently across genes.

Estimands (round3/estimands_and_method.md, Part A): μ_dj = E_{F_d}[θ_j]; Σ_B,jk = Cov_P(μ_dj, μ_dk); Σ_W,jk = E_P[Cov_{F_d}(θ_j, θ_k)].

Per-donor statistics: S_d = Σ_i s_i; ratio mean m_dj = Σ_i Y_ij / S_d; κ_d = Σ_i s_i² / S_d² (= 1/n_d at equal depth); x_ij = Y_ij / s_i; Ĉ_d = sample covariance of the x_i (divisor n_d − 1); σ̂²_dj = Var̂_i(x_ij) − mean_i(x_ij / s_i); ψ̂_dj = κ_d max(σ̂²_dj, 0) + m_dj / S_d.

## 2. The leak

**Lemma 1 [proved].** Under A1–A3, for j ≠ k,

  E[m_dj | F_d, design] = μ_dj,  Cov(m_dj, m_dk | F_d, design) = κ_d Cov_{F_d}(θ_j, θ_k),

and therefore, over the donor draw, Cov(m_dj, m_dk | design) = Σ_B,jk + κ_d Σ_W,jk. Moreover E[Ĉ_d,jk | F_d, design] = Cov_{F_d}(θ_j, θ_k).

*Proof.* Given F_d and s, E[x_ij | θ_i] = θ_ij, so E[m_dj | F_d, s] = Σ_i s_i μ_dj / S_d = μ_dj. For j ≠ k, Cov(x_ij, x_ik | θ_i, s_i) = 0 (independent Poisson draws), so Cov(x_ij, x_ik | F_d, s) = Cov_{F_d}(θ_j, θ_k), and cells are independent given F_d; hence Cov(m_dj, m_dk | F_d, s) = Σ_i s_i² Cov_{F_d}(θ_j, θ_k)/S_d². The law of total covariance over F_d ~ P gives the second display. Ĉ_d,jk is the usual unbiased sample covariance of n_d iid vectors x_i whose off-diagonal covariance is Cov_{F_d}(θ_j, θ_k). ∎

**Consequence.** Any donor-mean (pseudobulk) covariance targets Σ_B,jk + κ̄ Σ_W,jk, with κ̄ the average of κ_d (≈ mean 1/n_d). Pairs with Σ_B,jk = 0 and Σ_W,jk ≠ 0 (within-only co-expression) are called by pseudobulk correlation at a rate that grows with κ̄ Σ_W,jk / √(Var); this is the "leak" behind round-3's pseudobulk FDR of 0.14–0.21 in few-donor designs. The leak is known in form and estimable per donor: ℓ̂_d,jk = κ_d Ĉ_d,jk.

**Proposition 1 [proved; prior].** For fixed centring weights α, β (Σα_d = Σβ_d = 1) and fixed weights w_d, with r_dj = m_dj − Σ_e α_e m_ej and r_dk = m_dk − Σ_e β_e m_ek,

  Σ̂_B(w) = Σ_d w_d [r_dj r_dk − ℓ̂_d(1 − α_d − β_d) − Σ_e α_e β_e ℓ̂_e] / Σ_d w_d (1 − α_d − β_d + Σ_e α_e β_e)

is unbiased for Σ_B,jk. With w ≡ 1 and α = β = 1/D it reduces to S_B,jk − mean_d κ_d Ĉ_d,jk (round-3 CMB_mom).

*Proof.* Write c_d = Cov(m_dj, m_dk) = Σ_B,jk + λ_d with λ_d = κ_d Σ_W,jk (Lemma 1); donors are independent, so E[r_dj r_dk] = c_d(1 − α_d − β_d) + Σ_e α_e β_e c_e. Substitute c_d, use E[ℓ̂_d] = λ_d, and sum. ∎

*Prior art.* The equal-weight form is the BCES moment correction (Akritas & Bershady 1996) with a count-derived, per-donor off-diagonal error covariance; the weighted debiased form parallels dIVW (Ye et al. 2021). We claim neither; the count-derived per-donor leak (Lemma 1) is the single-cell-specific ingredient.

## 3. What any per-pair test can achieve

Working model **G**: m_d ~ N(μ̄, V_d) independently, V_d = Σ_B + Ψ_d with Ψ_d (sampling covariance) known; for the pair (j,k), write a_d = Σ_B,jj + ψ_dj, b_d = Σ_B,kk + ψ_dk, c_d = σ + ℓ_d, σ = Σ_B,jk.

**Theorem 1 (information bound) [proved for G].** The Fisher information for σ contributed by donor d is

  I_d(σ) = (a_d b_d + c_d²) / (a_d b_d − c_d²)².

At σ = 0 with no leak (ℓ_d = 0) the information about σ is orthogonal to (μ̄, Σ_B,jj, Σ_B,kk), so the efficient information is I = Σ_d 1/(a_d b_d). An oracle that observes μ_d has I* = D/(Σ_B,jj Σ_B,kk). Hence

  I / I* = (1/D) Σ_d ρ_dj ρ_dk,  ρ_dj = Σ_B,jj / (Σ_B,jj + ψ_dj)  (donor-specific reliability),

and for local alternatives σ = h/√D, no asymptotically level-α test based on genes j and k alone has power exceeding that of the oracle test with **D_eff = Σ_d ρ_dj ρ_dk** donors.

*Proof.* For a bivariate normal with covariance V = [[a, c], [c, b]], ∂V/∂c = E = [[0,1],[1,0]] and I_cc = ½ tr(V⁻¹ E V⁻¹ E). With V⁻¹ = Δ⁻¹[[b, −c], [−c, a]], Δ = ab − c², one has V⁻¹E = Δ⁻¹[[−c, b], [a, −c]] and tr((V⁻¹E)²) = 2(ab + c²)/Δ², giving I_d. The cross-information with a is ½ tr(V⁻¹ E V⁻¹ E_aa) = −c b /Δ²... ∝ c, which vanishes at c = 0 (likewise for b); the mean parameters are orthogonal to covariance parameters in the normal family. Information adds over independent donors. The oracle case is ψ ≡ 0. The power statement is the standard local asymptotic normality bound (Le Cam's third lemma; e.g. van der Vaart 1998, ch. 15) applied to the one-dimensional parameter σ with the nuisance parameters orthogonal at the null. ∎

**Corollary 1 (design) [proved for G].** With ψ_dj ≈ c_j / n_d (c_j = within-donor variance plus the Poisson term at unit depth), ρ_dj = 1/(1 + c_j/(Σ_B,jj n_d)). Doubling cells per donor raises D_eff by at most a factor 2 and saturates at D; adding donors raises D_eff linearly. In the 16-donor × ~40-cell development design the median noise-to-signal ratio ψ/Σ_B,jj is 0.73 (SEA-AD at 40 cells: 0.98), so D_eff ≈ 0.3·D ≈ 5: the best per-pair test at 16 donors has roughly the power of the oracle at 5 donors. **[empirical, development seeds]** This is why the brief's original gate 3 cannot be met by any per-pair method; it motivated the amended gate.

With a leak, I_d(0) = (a b + ℓ²)/(a b − ℓ²)² ≥ 1/(ab): a *known* leak does not cost information; an unmodelled leak costs validity (Lemma 1).

## 4. The proposed test: LCP (leak-corrected, precision-standardized donor permutation)

1. Per gene j: Paule–Mandel τ̂²_j, i.e. the root of Σ_d (m_dj − μ̂_j)²/(τ² + ψ̂_dj) = D − 1 with μ̂_j the corresponding inverse-variance mean (τ̂² = 0 if no root). v̂_dj = τ̂²_j + ψ̂_dj; r_dj = (m_dj − μ̂_j)/√v̂_dj.
2. Per pair: ℓ̃_d = κ_d Ĉ_d,jk / √(v̂_dj v̂_dk); T = [Σ_d (r_dj r_dk − ℓ̃_d)] · √(D / (D + Σ_d ℓ̃_d²)).
3. Reference: T_π = Σ_d r_π(d)j r_dk over uniformly random permutations π (20,000; refined to 4·10⁵ when p ≤ 5·10⁻⁴); p = (1 + #{|T_π| ≥ |T|})/(1 + B). BH over pairs.

Design choices, each tied to a failure seen in development:
- *Leak correction inside the statistic* (Lemma 1). Without it, within-only pairs are called (SS type-I 0.10 vs 0.058 with it, exploratory).
- *Standardization by model-based total variance.* It equalizes donor variances so that the permuted pairings are approximately exchangeable; the Paule–Mandel root makes Σ_d r_dj² = D − 1 exactly, so the permutation variance is fixed: Var_π(T_π) = (Σ_d r_dj² − D r̄_j²)(Σ_d r_dk² − D r̄_k²)/(D − 1) ≈ D − 1.
- *Permutation, not normal theory.* At 16 donors, normal-theory references — including Pearson on the true donor means — are anti-conservative in the BH-relevant tail (null p < 10⁻⁴ about 3× too often for the oracle, 13–23× for weighted regression versions), because a sum of 16 products of skewed donor effects is far from normal. **[empirical]**
- *Not studentized.* The studentized permutation statistic (DiCiccio & Romano 2017) protects the weak null but removed nearly all power at 16 donors (oracle BH power 0.003 vs 0.18 unstudentized, exploratory). We therefore test the independence null (below) and report the weak null as a stress case.

**Theorem 2 (validity of LCP) [sketch].** Let H0: the donor-level vectors (μ_dj)_d and (μ_dk)_d are independent (implies Σ_B,jk = 0; Σ_W,jk arbitrary).
(a) *Exact case.* If the leak is absent (ℓ ≡ 0, e.g. j and k measured in disjoint cell halves) and the standardized deviations r_·j are exchangeable across donors (e.g. equal designs and known nuisance parameters), the permutation p-value is exactly valid: P(p ≤ u) ≤ u for all u.
(b) *Asymptotic case.* Under A1–A3, H0, bounded designs (κ_d ≤ c/n_min), uniformly bounded (4+δ) moments, and consistent (μ̂_j, τ̂²_j), as D → ∞: T/√D ⇒ N(0, 1) and the permutation distribution of T_π/√D converges in probability to N(0, 1); hence the LCP test is asymptotically exact.
*Outline.* (a) Under H0 with ℓ = 0 the two standardized vectors are independent and one is exchangeable, so T is distributed as T_π. (b) Under H0, E[r_dj r_dk | design] = ℓ_d/√(v_dj v_dk) + o(1) by Lemma 1, so the summands r_dj r_dk − ℓ̃_d are independent, mean ≈ 0, with Var = 1 + ℓ̃_d² + γ_d where γ_d is a fourth-cumulant term of the within-donor sampling error of order 1/n_d (zero under Gaussian sampling error). Lyapunov's CLT gives T/√D ⇒ N(0, 1) after the √(D/(D + Σℓ̃²)) scaling (exact for γ ≡ 0). Hoeffding's combinatorial CLT applies to Σ_d r_π(d)j r_dk because Σ r_dj² = D − 1 and max_d r_dj²/D → 0 in probability, giving the same limit. Remaining items to verify: the effect of estimating (μ̂, τ̂²) on the variance of T (each enters only through one gene's vector, so it is fixed under the permutation of that gene's vector; its effect on the observed T is O_p(D^{-1/2})), and the γ_d term at n_d ≈ 10–40.
(c) *Weak null.* If Σ_B,jk = 0 but the donor effects are dependent, Var(r_dj r_dk) ≠ E r_dj² E r_dk² and the test is not asymptotically exact (the DiCiccio–Romano phenomenon). Development includes such a stress null (QQ_cross).

**Proposition 2 (efficiency of LCP) [sketch].** Under G with local alternatives, the noncentrality of T/√D is h·D⁻¹ Σ_d √(ρ_dj ρ_dk) / √(Σ_B,jj Σ_B,kk); the asymptotic relative efficiency against the bound of Theorem 1 is

  ARE = [D⁻¹ Σ_d √(ρ_dj ρ_dk)]² / [D⁻¹ Σ_d ρ_dj ρ_dk] ≤ 1,

with equality iff ρ_dj ρ_dk is constant across donors (Cauchy–Schwarz). The loss comes from giving every standardized donor equal weight, the price of exchangeability under permutation; score-optimal weights 1/(a_d b_d) make permuted pairings non-exchangeable.

## 5. Intervals (amended, change 1 of DEV_PREREG.md)

BOOT_T: equal-weight unbiased estimator (Proposition 1), sandwich SE from the per-donor leak-corrected products, donor bootstrap-t. **Negative result [empirical]:** with log-normal donor effects, linear-scale Σ_B,jk intervals undercover at every D we examined, *even on the true donor means* (normal theory 0.85–0.88 at D = 16–84; bootstrap-t 0.87–0.94 at D = 16, 0.91 at D = 40), because the per-donor products are heavy-tailed. The paper therefore claims calibrated tests, and reports intervals with their coverage relative to the oracle; it does not claim calibrated intervals below ~40 donors.

## 6. What is and is not new (honest summary; see novelty_v2.md)

Not new: the moment correction (BCES 1996), weighted debiasing (dIVW), permutation tests of association, studentized permutation (DiCiccio & Romano 2017), attenuation by measurement error (Spearman; Fuller 1987), Fisher information of a covariance. Candidate contribution: (i) the count-derived, per-donor leak and its explicit form Σ_B + κ̄Σ_W for pseudobulk (Lemma 1); (ii) a few-donor test that is calibrated in the BH-relevant tails where all normal-theory versions are not (LCP; Theorem 2); (iii) the per-pair information bound with D_eff = Σ ρ_dj ρ_dk, which quantifies what few-donor designs can and cannot detect and gives a cells-versus-donors design rule (Theorem 1, Corollary 1); (iv) the negative result for intervals. Whether (i)–(iv) clear RECOMB's bar is assessed in GO_NOGO.md.

## References (to be verified in novelty_v2_refs.csv before the paper)

Akritas & Bershady 1996 (ApJ 470:706, doi 10.1086/177901); Ye, Shao & Kang 2021 (dIVW, Ann Stat, doi in novelty_v2_refs.csv); DiCiccio & Romano 2017 (JASA, doi 10.1080/01621459.2016.1202117); Hoeffding 1951 (Ann Math Stat, combinatorial CLT); van der Vaart 1998 (Asymptotic Statistics, CUP); Paule & Mandel 1982; Fuller 1987; Spearman 1904.
