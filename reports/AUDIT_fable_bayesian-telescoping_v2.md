# Bayesian Telescoping V2 — second-vendor audit

**Manuscript:** `inbox/bayesian-telescoping/v2/bayesian_telescoping_v2.tex` ("Telescoping Bayesian Inference with Approximate Forward Operators: Corrected Convergence Theory…", V2, September 2026)
**Driver:** `v2/code/verify_bt_v2.py` — rerun by this auditor (exit 0, `ALL ASSERTIONS PASS`), output matches the shipped `verify_output.txt`.
**Comparison materials:** V1 (`paper.txt`) and the first-vendor V1 audit (`AUDIT_bayesian-telescoping_2026-09-22.md`).
**Auditor:** independent second-vendor adversarial pass (fable), 2026-09-22. Every numeral below was re-derived with the auditor's **own** code (`/tmp/audit_v2/check1–4.py`, independent grids, wider scans, higher precision, integrated — not peak-height — masses), not merely by rerunning the author's driver. Citation metadata checked against arXiv/Project Euclid.

---

## Verdict

**ACCEPT WITH MINOR REVISIONS.** V2 does what a correction paper should do. All three headline corrections check out under independent verification:

1. **MLMC regimes** — thresholds γ = β+2, log case at γ = β+2, case-(iii) exponent −2−(β+2−γ)/k, and the cutoff L ≍ ε^{−1/k} are all **correct**; V1's ε^{−1/(k+1)} cutoff was indeed wrong. I re-derived and numerically confirmed each (measured cost slopes −2.539/−3.251/−2.540 vs predicted −2.500/−3.250/−2.500).
2. **Bimodality** — the linear-obstruction theorem is correct, the non-injectivity mechanism is correct, and Table 1 reproduces exactly in my independent implementation (secondary mode −2.02 with mass ≈ 0.09–0.11 at n=1; −4.71 with mass ≈ 2×10⁻¹³ at n=2; strictly unimodal for n ≥ 3 on a scan of λ ∈ [−40, 15]). V2's "coarse-level transient" **refutes** V1's "persistent at n = 64".
3. **Heavy tails** — the demotion to a conjecture is warranted: my own quadrature (different grid, different n-range, plus a Cauchy prior V2 didn't test) gives slope ≈ −2.93 for Gaussian, t₁₀, t₅, t₃, t₂, and t₁ alike — the full −(k+1) rate, no degradation, contradicting V1's ρ(α) formula.

The residual problems are localized and fixable: the **fold-condition threshold as stated is inconsistent with the paper's own Table 1** (it uses the data-implied λ, under which n = 1 should already be unimodal — it isn't); the shipped **driver's case-(iii) assertion silently tests the superseded ε^{−1/(k+1)} scaling**, not the Theorem 4.1 statement; a **change-log line contradicts the corrected theorem** (/(k+1) vs /k); and a handful of sign/typo/labelling slips listed under "New-error hunt". None of these undermines the mathematics; all should be fixed before acceptance.

---

## Claim-by-claim

| # | Claim (V2) | Verdict | Evidence (auditor's own code unless noted) |
|---|---|---|---|
| 1 | Lemma 2.1: Σ_{n≥N} n^{−α} ≤ N^{−α} + N^{1−α}/(α−1) ≤ C_α N^{1−α}; tail = Θ(N^{1−α}); head Σ_{n<N} n^{−α} = Θ(1) for α>1 | **CORRECT** | mpmath ζ-function check, α ∈ {1.5, 2, 3}, N ∈ {10, 100, 1000}: sum always in [N^{1−α}/(α−1), N^{−α}+N^{1−α}/(α−1)]. Head sums → ζ(α) < ∞. V1's reversed inequality confirmed false at every (α,N) tested (e.g. α=2, N=10: 0.10517 > 0.10000). |
| 2 | Thm 4.1: regimes split at γ = β+2; log² at boundary; case (iii) cost O(ε^{−2−(β+2−γ)/k}) | **CORRECT** | Σ n^{(β−γ)/2} converges iff (β−γ)/2 < −1 ⇔ γ > β+2 (elementary + numeric). Empirical cost slopes match −2−(β+2−γ)/k to <2% (see MLMC section). Boundary case: cost·ε² / (log ε / k)² → ≈ 1.05. |
| 3 | Cutoff remark: level-L bias = O(L^{−k}), so L ≍ ε^{−1/k}; V1's ε^{−1/(k+1)} wrong | **CORRECT** | Re-derived: bias ≤ Σ_{n≥L} Cn^{−(k+1)} = Θ(L^{−k}/k) by Lemma 2.1 with α = k+1. Setting Θ(L^{−k}) = ε forces L = Θ(ε^{−1/k}). V1's L = ε^{−1/(k+1)} gives bias Θ(ε^{k/(k+1)}) ≫ ε — insufficient. Numerically confirmed for k=1,2,3. **But see New-error hunt #2: the shipped driver's case-(iii) test does not test this.** |
| 4 | Thm 5.1: affine Fₙ + Gaussian prior + Gaussian noise ⇒ Gaussian (unimodal) posterior at every level | **CORRECT** (one phrasing quibble) | −Φₙ is a concave quadratic; exp(−Φₙ)dμ₀ is Gaussian; standard conjugacy. Verified in closed form. Quibble: the "in particular no choice of discretisation of a linear forward operator" clause silently assumes the discretisation itself is affine (an adaptive-step or nonlinear scheme applied to a linear ODE need not be); the theorem's hypothesis is fine, the corollary phrasing over-reaches slightly. |
| 5 | §5.2: z* = −1.5960716, R(z*) = 0.2703948; T(z) increasing on z<2; Table 1 | **CORRECT numerically** (criterion misstated — see Bimodality) | My root-find: z* = −1.59607163798, R(z*) = 0.270394765205 (matches to all printed digits). R′ has exactly one real root; R has no real roots (so R^n > 0, injectivity governed by monotonicity alone). Table 1 fully reproduced independently, incl. integrated masses. |
| 6 | §6: no rate degradation for heavy tails; V1 Table 9 unreproducible; ρ(α) unproved | **CORRECT** | Independent quadrature ([−60,60], 480k points, n = 40…160 in steps of 4): slopes −2.928 for Gaussian/t₁₀/t₅/t₃/t₂ **and Cauchy (t₁)** — flat contradiction of ρ(α)(k+1) = −2.667/−2.25/−1.5/0. Prop 6.1 is a valid sufficient condition (see Heavy tails). |
| 7 | Prop 3.3 (credible intervals, CDF version) correct; moment version false | **CORRECT** | Counterexample confirmed: N(0,1) vs Logistic(0, √3/π) have identical mean 0 and variance 1, both log-concave, but q₀.₉ = 1.28155 vs 1.21139 (Δ ≈ 0.070). Take μₙ ≡ logistic, μ = normal: moments converge exactly, quantile error is a constant. |
| 8 | §7 TKF cost Θ(T²), no asymptotic advantage | **CORRECT** | Σ_{t≤10⁴} ⌈√t⌉² = 5.067×10⁷ ≈ T²/2 = 5.0×10⁷ ≫ T^{3/2} = 10⁶. Matches driver's 5.07×10⁷. |
| 9 | Citations: O'Brien/Moores arXiv:2502.11510; Agapiou–Castillo Ann. Statist. 52(4) 1433–1459, 10.1214/24-AOS2397 | **CORRECT** | arXiv:2502.11510 = "Here Be Dragons: Bimodal posteriors arise from numerical integration error in longitudinal models", Tess O'Brien, Matthew T. Moores, et al. (2025). Project Euclid confirms Agapiou & Castillo, "Heavy-tailed Bayesian nonparametric adaptation", Ann. Statist. 52(4):1433–1459, Aug 2024, DOI 10.1214/24-AOS2397. Both fix V1's mis-citations. |
| 10 | Hellinger normalisation (§2) | **CORRECT** | ∫(√f−√g)² convention gives range [0,√2] and d_H = √(2(1−ρ)); the ½∫ convention gives [0,1] and √(1−ρ). Both stated correctly; V1's mixing acknowledged. |
| 11 | Increment identity (Prop 3.1) & Cor 3.2 | **CORRECT** | Elementary expansion; re-derived. Cor 3.2 needs the (implicit) growth bound ‖Fₙ(x)‖ ≤ C(1+‖x‖^p) — used but never displayed in V2; state it. |

---

## MLMC threshold+cutoff

**Thresholds (Thm 4.1 i–iii): verified.** With s = (β−γ)/2, Σ_{n<L} n^s = O(1) iff s < −1 ⇔ γ > β+2; = Θ(log L) at s = −1 ⇔ γ = β+2; = Θ(L^{s+1}) = Θ(L^{(β+2−γ)/2}) for s > −1. This is calculus; the proof in the paper is correct. V1's γ > β+1 threshold contained divergent sums in its "O(1)" case (Σ_{n<10⁶} n^{−3/4} ≈ 123 and growing — reproduced).

**Cutoff: V2 is right, V1 was wrong.** The increment bound is Cn^{−(k+1)}; the level-L bias is the *tail sum*, which is Θ(L^{−k}) (Lemma 2.1, α = k+1 — independently confirmed to within 1.4% of L^{−k}/k for k = 1,2,3). Bias ≤ ε therefore requires L ≍ ε^{−1/k}. V1's L ≍ ε^{−1/(k+1)} yields bias Θ(ε^{k/(k+1)}) ≫ ε and was internally inconsistent with V1's own bias corollary — V2's diagnosis is accurate.

**Case-(iii) exponent: verified empirically with my own code.** Using the *correct* cutoff L = ⌈(1/(kε))^{1/k}⌉ and cost = ε^{−2}(Σ_{n<L} n^s)²:

| (β, γ, k) | measured d log Cost / d log ε | predicted −2−(β+2−γ)/k |
|---|---|---|
| (1, 2, 2) | −2.539 | −2.500 |
| (2, 1.5, 2) | −3.251 | −3.250 |
| (1, 2.5, 1) | −2.540 | −2.500 |

Boundary γ = β+2: cost·ε²/(k^{−1}log ε)² → 1.05–1.10 over ε = 10⁻⁴…10⁻⁸, consistent with O(ε^{−2}(log ε)²). **Theorem 4.1 is correct as stated.**

**However — two related defects (see New-error hunt):** (a) the shipped driver's case-(iii) assertion uses `L = eps**(-1/(k+1))` and `pred = eps**(-2-(beta+2-gamma)/(k+1))` — i.e. it verifies the *superseded V1 scaling pair*, not Theorem 4.1, while its printed label claims the /k exponent; (b) the change log states the corrected exponent as "−2−(β+2−γ)/(k+1)", contradicting the theorem's /k. Neither invalidates the theorem (which I verified independently), but a paper whose reproducibility statement is "asserts every claimed sign and threshold" must actually assert the stated theorem.

---

## Bimodality

**Theorem 5.1 (obstruction): correct.** Affine Fₙ ⇒ −Φₙ concave quadratic ⇒ Gaussian posterior at every level ⇒ log-concave, unimodal. Verified in closed form (posterior N(0.526012, 0.023121) for the driver's test point — matches).

**Stability-function facts: all reproduced to full precision.** R′(z) = 1+z+z²/2+z³/6 has exactly one real root z* = −1.59607163798 (the other two roots are complex −0.702±1.807i); R(z*) = 0.270394765205; R has **no real roots**, so R(λh)^n > 0 always and injectivity is governed solely by R's monotonicity; T′(z) = (1−z/2)^{−2} > 0.

**Table 1: fully reproduced with an independent implementation** (4×10⁶-point scan over λ ∈ [−40, 15], integrated masses by splitting at inter-mode minima rather than peak heights):

| n | RK4 modes (mine) | secondary relative mass (mine) | paper |
|---|---|---|---|
| 1 | −1.0198, −2.0160 | 9.4×10⁻² (integrated); 1.07×10⁻¹ (peak height) | −1.02, −2.02, "mass 0.11" |
| 2 | −1.0007, −4.7139 | 2.5×10⁻¹³ (integrated); 1.9×10⁻¹³ (height) | −1.00, −4.71, "2×10⁻¹³" |
| 3–64 | single mode ≈ −1.000 | — | 1 mode (n≥4 rows) |

Exact and trapezoidal posteriors: unimodal at every n tested (n = 1 trapezoidal mode at −0.928 → −1.000 as n grows). At n = 3 I additionally did a high-precision derivative sign-scan around the secondary-branch preimage λ ≈ −7.70: **no critical point exists** (the prior slope dominates; log-posterior there is ≈ −86 below the main mode). So "observable bimodality is a coarse-level phenomenon" is **correct**, and it **refutes V1's "persistent at n = 64"** — V1's persistence claim came from code that hard-coded a bimodal mixture prior (established in the V1 audit; V2's account of this is accurate).

**But the fold condition as stated is wrong — inconsistent with the paper's own Table 1.** Eq. (5.3)/(nstar) says the map is injective (bimodality impossible) when |λ|T ≤ |z*|, and the abstract says the pathology "disappears once n > |λ|T/z*". Take the paper's own experiment: λ* = −1, T = 1, so |λ|T = 1 < 1.5961 — the criterion **as literally stated predicts unimodality at every n ≥ 1, including n = 1**, where the paper's own Table 1 (and my reproduction) shows a secondary mode of mass ≈ 0.1. The resolution is that the relevant λ in the criterion is not the data-implied value but the largest |λ| carrying appreciable **prior** mass: the fold at level n sits at λ_fold = z*n/T, and the n = 1 fold (−1.596) lies only 1.19 prior-sd from the prior mean, with the secondary preimage at −2.09 (2.2 prior-sd, still supported), while for n = 2 the secondary preimage is at −4.91 (7.8 sd) and for n = 3 at −7.70 (13.4 sd) — hence the 10⁻¹³ and then nothing. Using the prior-edge scale (m₀ − 3s₀ = −2.5), the threshold n < 2.5T/|z*| ≈ 1.57 predicts "bimodal only at n = 1", exactly matching Table 1. The mechanism and numbers are right; **the displayed inequality needs the prior-supported λ-range, not the point value** — as written it contradicts the table one inch below it.

Three smaller slips in the same subsection: (i) "injective only on {λh ≤ z*}" has the inequality **backwards** (R is increasing on z ≥ z*, decreasing on z ≤ z*; injectivity holds on λh ≥ z*); (ii) "the fold sits at λ = ±z*n/T" — there is no fold on the positive side (R is strictly increasing for z > z*); (iii) the abstract writes z* ≈ 1.5961 (positive) while §5.2 defines z* = −1.5960716.

**Trapezoidal "can never fold": over-claimed.** T(z) is negative for z < −2, so for **even** n the n-step map λ ↦ T(λh)^n is *not* injective on ℝ: at n = 2, T = 1 the equation T(λ/2)² = e⁻¹ has a second preimage at λ = −16.332 (my high-precision computation). It generates no posterior mode (relative mass ~10⁻⁶⁸⁹ — I checked; and the second preimage needs |λ| > 2n/T, receding with n), so Table 1 is unaffected, but the correct statement is "T is strictly monotone where it is positive / on the prior-supported region", which is in fact what the driver's own assertion string says. The paper text should match.

---

## Heavy tails

**Prop 6.1 (sufficient condition): valid.** Given |ΔΦₙ| ≤ C(1+‖x‖^{2p})n^{−(k+1)} and the uniform exponential moment condition (2.1), the bound E[|e^{−ΔΦ}−1|] ≤ Cn^{−(k+1)}E[(1+‖x‖^{2p})e^{C(1+‖x‖^{2p})n^{−(k+1)}}] closes, and Zₙ/Z_{n+1} = 1 + O(n^{−(k+1)}) follows the same way. The proof is a sketch that leans on V1's Theorem-6.1 split, but the steps shown are correct and standard.

**Numerics: reproduced independently and extended.** With my own quadrature (grid [−60,60] × 480 001 points, n = 40…160 step 4 — different from the driver's 60…198 step 2) on F_n(θ) = θ(1−n⁻²), y = 1.5, σ = 0.1:

| prior | Gaussian | t₁₀ | t₅ | t₃ | t₂ | t₁ (Cauchy, my addition) |
|---|---|---|---|---|---|---|
| my slope | −2.928 | −2.928 | −2.928 | −2.928 | −2.928 | −2.928 |

Full rate −(k+1) = −3 everywhere, **no degradation even at Cauchy tails**, against V1's predicted −2.667/−2.25/−1.5/0. So: V1's Table 9 (−2.61, −2.19, −1.47) is not reproducible from V1's own model, and ρ(α) = (α−2)/(α−1) is unsupported. **Verdict: V2's §6 conclusion is sound.**

One structural observation V2 gets right and should perhaps make louder: in this model the *likelihood* is Gaussian in θ with coefficient 1−n⁻² bounded away from 0, so the **posterior** has Gaussian tails regardless of the prior — condition (2.1) actually *holds* here, Prop 6.1 applies, and the experiment could never have exhibited a penalty. V2's conjecture paragraph does acknowledge this ("which the model above does not satisfy" for the unbounded-increment requirement), and the diagnosis "the sufficient condition fails ⇒ the rate degrades is the invalid inference in V1" is exactly right.

Two cosmetic discrepancies: manuscript table says measured slope −2.988; the shipped driver prints −2.973 (my run confirms −2.973). Both round to −3 and nothing depends on the third digit, but the paper claims the driver "regenerates every numeral". And the driver's t₂ line prints "V1's rho(alpha)*(k+1)=−3.000" where ρ(2) = 0 should give 0.000 — a Python falsy-zero bug (`-3.0*rho if rho else -3.0`); the manuscript's own table correctly shows 0.000.

---

## Credible intervals

**The V2 statement (Prop 3.3) is correct and the V1 moment-convergence version is genuinely false.**

- *Counterexample to the moment version (confirmed with my own computation):* N(0,1) and Logistic(0, √3/π) both have mean 0, variance 1, and are log-concave, yet their 0.9-quantiles differ by 0.0702 (1.28155 vs 1.21139). Setting μₙ ≡ Logistic for all n and μ = N(0,1), the first two moments converge trivially (they are equal) while the quantile error never decreases. So convergence of the first two moments cannot control quantiles, exactly as V2's remark asserts, and V1's Corollary 7.2 was false as stated.
- *The CDF version:* with (i) a uniform density floor c > 0 near q_α and (ii) |Fₙ(q_α) − α| = O(δₙ), the mean-value argument gives |q_α^{(n)} − q_α| ≤ O(δₙ)/c. Correct. Minor polish: the proof's "once |q_α^{(n)} − q_α| < δ" needs the one-line bootstrap (|Fₙ(q_α)−α| ≤ cδ ⇒ the quantile is within δ) to avoid apparent circularity; and hypothesis (i) is honestly flagged as an assumption on the approximate posteriors, automatic in the log-concave (affine/Gaussian) case. Good.

---

## New-error hunt

Ordered by importance. (None overturns a theorem; #1 and #2 are the ones a referee should insist on.)

1. **Fold-condition threshold contradicts the paper's own Table 1** (§5.2, eq. (nstar), abstract, and the "correct diagnostic" remark). With the paper's own parameters |λ|T = 1 < |z*| = 1.596, the stated criterion predicts unimodality already at n = 1 — but Table 1 (correctly) reports an observable second mode of mass ~0.1 at n = 1. The criterion must be stated with the prior-supported λ-range (e.g. |m₀| + 3s₀ = 2.5 ⇒ bimodal only for n < 1.57, matching the table), not the data-implied point value. As printed, the paper's flagship "practitioner check" fails on the paper's flagship example.
2. **The shipped driver does not verify Theorem 4.1(iii) as stated.** `verify_bt_v2.py` §2 uses `L = eps**(-1/(k+1))` and `pred = eps**(-2-(beta+2-gamma)/(k+1))` — the *superseded* V1 cutoff/exponent pair — under a printed label claiming the /k exponent, and passes only because the two errors cancel (the pair is internally consistent) and the tolerance is a full order of magnitude. The manuscript's claim that the driver "asserts every claimed sign and threshold" is therefore not true for case (iii). (I verified the /k statement independently; it is correct — the defect is in the driver, not the theorem.)
3. **Change-log contradicts the corrected theorem:** it records the corrected exponent as "−2−(β+2−γ)/(k+1)", while Theorem 4.1(iii) and the cutoff remark (correctly) use /k. Stale line from the pre-cutoff-fix draft.
4. **"Injective only on {λh ≤ z*}" is the wrong direction** (should be λh ≥ z*: R decreases left of z*, increases right of it); and **"the fold sits at λ = ±z*n/T"** — there is no positive-side fold. Abstract also gives z* with the wrong sign (≈ 1.5961 vs defined −1.5960716).
5. **Trapezoidal "can never fold" is false as a statement about the n-step map:** for even n, T(λh)^n is non-injective via the negative branch of T (second preimage at λ = −16.332 for n = 2, T = 1, y = e⁻¹; verified at 60-digit precision). Its posterior contribution is ~10⁻⁶⁸⁹ and recedes as |λ| > 2n/T, so all conclusions stand, but the text should say "strictly increasing where positive / injective on the prior-supported region" (as the driver's own assertion string already does).
6. **"Mass" in Table 1 is relative peak height, not posterior mass.** The driver reports exp(Δlog-posterior) at the peak; the table's caption promises "relative posterior mass". Integrated masses are 0.094 (n=1) and 2.5×10⁻¹³ (n=2) vs the printed 0.11 and 1.9×10⁻¹³ — same order, different quantity. Relabel or integrate.
7. **Driver/manuscript numeral mismatch:** heavy-tail slope −2.988 (manuscript) vs −2.973 (shipped driver output and my rerun) — the "regenerates every numeral" claim is not literally met. Also the stale, truncated `verify_output_final.txt` sits next to `verify_output.txt` in `code/`; delete one.
8. **Driver falsy-zero bug** (cosmetic): the t₂ row prints "V1's rho(alpha)*(k+1) = −3.000" instead of 0.000 (`-3.0*rho if rho else -3.0` with rho = 0). The manuscript's table has the correct 0.000.
9. **Thm 5.1 "in particular" clause over-reaches:** "no choice of discretisation of a linear forward operator can induce bimodality" tacitly assumes the discretisation is itself affine; a nonlinear scheme (e.g. adaptive stepping with x-dependent steps) applied to a linear problem is not covered by the hypothesis. One qualifying word ("affine/linear discretisation") fixes it.
10. **Cor 3.2 uses an undisplayed growth assumption** ‖Fₙ(x)‖ ≤ C(1+‖x‖^p) (needed to turn ‖y − Fₙ(x)‖ into (1+‖y‖)(1+‖x‖^p)); V2's §1.3 states only the increment bound. State it.
11. Table 1 text: "lies outside the prior-supported region for n ≥ 4" — in fact the secondary critical point is already gone at n = 3 (no stationary point exists near the second preimage λ ≈ −7.70; verified by derivative sign-scan). Understatement, not error; but say n ≥ 3.

Checked and NOT problems: Lemma 2.1's constant C_α = 1 + 1/(α−1) ✓; the head-sum Θ(1) note ✓; the boundary-case (log ε)² cost ✓; the Hessian formula (5.1) including the middle term ✓ (matches the correction demanded in the V1 audit); the Hellinger dual-convention statement ✓; Θ(T²) TKF ✓; both repaired citations ✓; the increment identity ✓; the claim that R has an interior minimum (not merely a critical point) ✓ (R″(z*) > 0; R > 0 on ℝ).

---

## What is sound

- **Lemma 2.1** (tail bound, both directions, head-sum note) — verified independently at high precision; the V1 inequality-direction error is correctly diagnosed and fixed.
- **Theorem 4.1** (MLMC regimes) — thresholds, boundary log², and case-(iii) exponent all confirmed; the cutoff correction L ≍ ε^{−1/k} is right and V1's ε^{−1/(k+1)} was wrong, exactly as V2 says.
- **Theorem 5.1** (linear obstruction) — correct, and genuinely the sharpest correction in the paper: it pins the pathology on nonlinearity in the parameter.
- **§5.2 numerics** — z*, R(z*), Table 1 (mode locations and masses), trapezoidal/exact unimodality: all reproduced by an independent implementation with a wider scan and integrated masses. The refutation of V1's "persistent at n = 64" stands.
- **§6** — Prop 6.1 is a valid sufficient condition; the no-degradation finding is robust (holds even for a Cauchy prior); the demotion of ρ(α) to a conjecture, with honest discussion of why the V1 experiment could never have tested it, is the right call.
- **Prop 3.3** — correct as restated; the falsity of the moment-convergence version is confirmed by explicit counterexample.
- **§7 TKF** — Θ(T²) accounting correct; advantage claim properly withdrawn.
- **Bibliography** — the two previously broken references now point to the right works (verified against arXiv and Project Euclid).
- The change-log/retraction framing is accurate to what the V1 audit found (invalid Hessian argument, hard-coded Experiment 5, off-by-one regimes, mis-citations), with one stale line (New-error #3).

## Recommended repairs

1. **Fix the fold criterion** (New-error #1): restate eq. (nstar), the abstract clause, and the "correct diagnostic" remark in terms of the prior-supported λ-range (fold at λ = z*n/T is relevant iff |z*|n/T lies within the region of non-negligible prior mass), and reconcile with Table 1's n = 1 row. Also fix the inequality direction "{λh ≤ z*}", the spurious "±", and the abstract's sign of z*.
2. **Fix the driver's case-(iii) assertion** to use L ≍ ε^{−1/k} and exponent −2−(β+2−γ)/k as the theorem states (and tighten the order-of-magnitude tolerance); regenerate `verify_output.txt`; delete the stale `verify_output_final.txt`.
3. **Fix the change-log exponent** to −2−(β+2−γ)/k.
4. **Weaken the trapezoidal claim** to "strictly increasing where positive; injective on the prior-supported region for all n" (New-error #5).
5. **Relabel or recompute "mass"** in Table 1 (integrated mass 0.094 / 2.5×10⁻¹³, or caption "relative peak height"); change "for n ≥ 4" to "for n ≥ 3".
6. **Sync the heavy-tail slope numeral** (−2.988 vs −2.973) with the shipped driver, and fix the t₂ print bug.
7. Add the missing growth assumption ‖Fₙ(x)‖ ≤ C(1+‖x‖^p) to §1.3/Cor 3.2; add the one-line bootstrap in Prop 3.3's proof; qualify Thm 5.1's "in particular" clause with "affine discretisation".

**Bottom line:** the three headline corrections are real and verified; the paper's honest downgrade of its former claims is accurate; remaining defects are statement-level (one self-contradicting threshold, one driver assertion that tests the wrong formula, and typo-grade slips) and are all fixable in a minor revision.
