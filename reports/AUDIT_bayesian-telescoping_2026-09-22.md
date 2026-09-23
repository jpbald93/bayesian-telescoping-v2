# Audit — *Bayesian Telescoping Under Non-Ideal Conditions: Numerical Pathologies, Heavy-Tailed Priors, and Adaptive Stopping*

**Author:** Josh Bald (Independent Researcher) · **Dated:** February 2026 · **Length:** 51 pp
**Artifact:** `Bayesian Telescoping Under Non-Ideal Conditions ... (1).pdf`
**sha256:** `1d6ccdde7829bb658aeefdd395196b528a0bbec65d7f56834852144b904175d6`
**Auditor:** independent pass, 2026-09-22. Every number below was produced by the committed script
`inbox/bayesian-telescoping/audit/verify_bt.py` (no arithmetic from memory), plus Crossref/arXiv/GitHub
metadata checks.

---

## Verdict: **MAJOR REVISION** (leaning reject)

The classical telescoping scaffolding (Sections 3–7) is largely correct and carefully written. But the
paper's **three headline contributions are individually unsupported**:

1. **Theorem 9.3** (bimodality theory) — the proof of part (b) is **not valid**, the stated
   characterisation is **inconsistent with the paper's own example**, and the experiment offered as its
   empirical confirmation (**Experiment 5**) is **not produced by the published code** — that code
   hard-codes a bimodal prior and performs no discretisation at all.
2. **Theorem 10.1** (heavy-tail penalty) — the "proof sketch" is hand-wavy, the sharpness claim is
   asserted rather than shown, and the result is attributed to a reference whose **DOI points to a
   different paper**.
3. **Theorem 11.2 / Theorem 8.4** (adaptive stopping, MLMC complexity) — the complexity regime
   **thresholds are off by one** (an algebraic slip visible in the paper's own line), and the "oracle
   levels" quoted in Tables 2 and 10 do not match the formula the paper states.

In addition, the repository that the paper advertises as containing "all code, data, and reproduction
scripts" **cannot reproduce any experiment** (missing scripts, missing modules, `run_all_experiments.sh`
invokes three files that do not exist).

---

## Verification basis

- `verify_bt.py` — arithmetic/algebra/self-consistency checks (Sections A–F below).
- Crossref API — DOIs for [16],[17],[20],[21],[22],[23],[24],[25],[26],[37],[38].
- arXiv API + `arxiv.org/abs/0909.2126`.
- GitHub API — the advertised repository and its single committed artifact.

---

## A. The published "code and data" cannot reproduce anything — **critical**

The paper states (§ *Code and Data Availability*, twice) that all code/data/scripts are at
`https://github.com/jpbald93/bayesian-telescoping`, and gives a repository tree with
`src/{forward_models,bayesian_inference,mlmc,telescoping}.py` and seven `experiments/experimentN_*.py`.

**What actually exists.** The repo has a single commit (2026-01-25) containing one file,
`bayesian-telescoping-repo.zip` (**8 052 bytes**). Unpacked, it contains 10 text files totalling < 20 KB:

- `src/` contains **only** `kalman_telescoping.py`. The four modules the paper lists — and that the
  shipped experiments import — **do not exist**. So `experiment4/5/6` all fail at `import`.
- `experiments/` contains **only** `experiment4_2d_poisson.py`, `experiment5_bimodal.py`,
  `experiment6_heavytail.py`. **`experiment1_heat.py`, `experiment2_diffusion.py` and
  `experiment3_mlmc.py` — the producers of Tables 3, 4 and 5 — do not exist.**
- `scripts/run_all_experiments.sh` runs, in order, `experiments.experiment1_heat`, `…experiment2_diffusion`,
  `…experiment3_mlmc`. With `set -euo pipefail`, **the advertised one-command reproduction fails at the
  first line.**
- Several figures/tables (e.g. Tables 7, 11) have no corresponding script; `tests/test_forward_heat.py`
  imports the missing `src.forward_models`.

**Experiment 5 does not test what the paper claims.** Section 13.6 presents Experiment 5 as the
empirical confirmation of Theorem 9.3: "RK4 (oscillatory error, k = 1) … credible-interval width decays at
rate ≈ n^−0.93 … Trapezoidal rule … ≈ n^−2.01". The shipped `experiment5_bimodal.py` contains **no ODE,
no RK4, no trapezoidal rule, no step size, and no discretisation of any kind**. Instead:

```python
def mixture_prior_logpdf(theta):          # 0.5 N(-2,1) + 0.5 N(2,1)  <-- bimodality hard-coded
def loglik(theta, y, sigma):              # y = theta + eta     <-- linear observation, no method
...
sigma_n = sigma_obs / np.sqrt(n)          # "refinement" = just shrinking the observation noise
```

The bimodality is **imposed by the prior by construction**; the "refinement levels" only multiply the
likelihood precision. This is the exact opposite of the paper's claim (that *discretisation error* induces
bimodality). The numbers in Table 8 / Experiment 5 therefore have **no reproducible source**, and the
one experiment offered as proof of the paper's flagship theorem is a different experiment entirely.

**Experiment 6 implements 3 of the 5 priors of Table 9.** The shipped `priors` dict is
`{gaussian, laplace, student_t3}`; the paper's Table 9 reports Gaussian, Laplace, **Student-t ν=10**,
**Student-t ν=5**, Student-t ν=3. Two of the five rows cannot be produced by the shipped code.
The shipped forward model (`F_n(θ) = θ(1 − n^{-2})`) also bears no resemblance to the paper's
PDE/MLMC framework.

> Consequence: **every quantitative table in the paper is unreproducible from the published artifact.**
> The reproducibility statement should be withdrawn, or a complete, runnable implementation supplied.

---

## B. Theorem 9.3 — the bimodality theory (flagship contribution 1)

**(i) The proof of part (b) is invalid.** The argument computes the Hessian of the negative log-posterior as

```
H_n(x) = ∇F_n(x)ᵀ Γ⁻¹ ∇F_n(x) + Prior⁻¹        (paper)
```

and concludes it is positive definite. But the true Hessian of
`½‖Γ^{-1/2}(y − F_n(x))‖²` carries a **second-derivative term**:

```
H_n(x) = ∇F_n(x)ᵀ Γ⁻¹ ∇F_n(x)  −  ⟨Γ⁻¹(y − F_n(x)), ∇²F_n(x)⟩  + Prior⁻¹ .
```

The omitted term `−⟨Γ⁻¹(y−F_n(x)), ∇²F_n(x)⟩` is an **indefinite, un-signed, O(1) matrix** (it scales with
the residual, not with n⁻ᵏ⁺¹). Log-concavity cannot be concluded by dropping it. So the mechanism
"sign-definite error ⟹ H_n ≻ 0 ⟹ unimodal" is not established. (Part (a) is nearly vacuous: "for any fixed n
there exists a prior and observation model such that µ_n is bimodal" is satisfiable for *any* F by choosing
a bimodal prior — it does not isolate a discretisation mechanism.)

**(ii) The stated characterisation contradicts the paper's own example.** Section 9 opens: "numerical
methods with **odd-order error (e.g., Euler, RK4)** are prone to bimodality, while **even-order methods
(trapezoidal, spectral) are not**." But **RK4 is a 4th-order (even-order) method**. The paper's own
Example 9.2 puts RK4 in the "oscillatory" class while labelling its telescoping order `k = 1`
(Tables 8, 9, and §13.6 all say "RK4 … k = 1"). RK4 has k = 4. The "odd/even order" heuristic and the
"k < 2 / k ≥ 2" threshold are two different statements that the paper conflates.

**(iii) The sign-oscillation claim needs care.** For fixed λh, the RK4 error `R(z) − e^z` has a *definite*
sign (the script confirms: `R(z) − e^z` is `+2.4e−4` at z=−0.5, `−2.8e−6` at z=+0.2, i.e. sign = −sign(z)).
So the increment does not "oscillate in n" at fixed θ; what varies is the *sign of θ*, which is what makes
`F_n(θ) = θ_RK4(λ)` **non-injective** (the actual O'Brien–Moores mechanism). The theorem's "oscillation"
concept (Definition 9.1, sign change *as a function of θ*) is therefore not the same as "oscillatory error",
and the bridge between them is not supplied.

---

## C. Theorem 8.4 / 11.2 — MLMC complexity regimes are **off by one** (flagship contribution 3)

Theorem 8.4's own line reads
`O(L^{(β−γ)/2+1}) = O(L^{(β+1−γ)/2})`, and its summary says
"if (β−γ)/2 < −1/2 *(i.e. γ > β + 1)* the sum converges to O(1)". Both are wrong:

- `(β−γ)/2 + 1 = (β+2−γ)/2`, **not** `(β+1−γ)/2`. (Script: β=1, γ=3 → left = 0, printed right = −1/2.)
- `Σ_{n=0}^{L−1} n^{(β−γ)/2}` converges **iff** `(β−γ)/2 < −1`, i.e. **γ > β + 2** (not β+1);
  it is `O(log L)` when γ = β + 2 (not β+1) and `L^{(β+2−γ)/2}` otherwise.

Numeric confirmation (script §2): at β = 1, γ = 2.5 the exponent is s = −0.75 > −1, so with
L = 10⁶ the partial sum is ≈ 123 **and still growing** — yet the theorem's case (i) declares γ > β+1 ⇒ O(1).
At β = 1, γ = 2.0 (s = −0.5, the boundary of the paper's scheme) the sum is ≈ 2·10³ and grows.

Corrected statement: case (i) γ > β+2 → O(ε⁻²); case (ii) γ = β+2 → O(ε⁻² (log ε)²);
case (iii) γ < β+2 → O(ε^{−2−(β+2−γ)/(k+1)}). Since Theorem 11.2 (adaptive stopping) rests on this
regime analysis, its O(ε⁻²) conclusion is likewise unproven as stated. (Also, the paper's derivation
γ = 2(k+1) is only asserted, via "attainable under standard couplings".)

**Corroborating inconsistency (Tables 2 & 10).** Both tables give "Oracle levels L*" = 5, 11, 24 for
ε = 10⁻¹, 10⁻², 10⁻³ while §11.3 states `L* = ⌈ε^{−1/3}⌉`. That formula gives **3, 5, 10** — the tabulated
levels are not the stated oracle levels. Table 6's text is likewise inconsistent: it states
"ε^{−1/3} = 100", which requires ε = 10⁻⁶ (i.e. it substitutes ε² for ε); with the paper's own
ε² = 10⁻⁶, ε^{−1/3} = 10.

---

## D. Theorem 10.1 — heavy-tail penalty: unproven, and mis-attributed (flagship contribution 2)

**(i) The result is asserted, not proved.** The "Proof Sketch" of Theorem 10.1 ends with a truncation
argument that *does not close*: after getting `O(n^{−(k+1)+2p/(α−1)})` it says "a more refined analysis
(tracking constants) yields the optimal exponent … the 2p term cancels in the optimization over the
truncation threshold." No such optimisation is exhibited, and `2p/(α−1)` does not cancel by any shown
step. The "sharpness" lower bound ("Step 3") is announced ("Direct calculation … shows the lower bound")
without the calculation. So the central formula

```
rate  ∝  ρ(α)·(k+1),   ρ(α) = (α−2)/(α−1)
```

is **not derived anywhere in the paper** — least of all its matching lower bound.

**(ii) The cited foundation does not exist as cited.** Reference **[38]** is given as
*S. Agapiou and I. Castillo, "Heavy-tailed Bayesian inversion", Ann. Statist. 52 (2024), 1234–1260,
DOI 10.1214/24-AOS2356*. Crossref shows **10.1214/24-AOS2356 = Stankewitz et al., "Early stopping for
L2-boosting in high-dimensional linear models", Ann. Statist. 52 (2024)** — a different paper entirely.
The genuine Agapiou–Castillo work is *"Heavy-tailed Bayesian nonparametric adaptation"*,
**Ann. Statist. 52(4) (2024), 1433–1459, DOI 10.1214/24-AOS2397** (arXiv:2308.04916). Section 10 is built
on this citation, so its framing ("Agapiou and Castillo prove optimal contraction rates for inverse
problems with heavy-tailed priors") should be re-checked against the real paper.

**(iii) The empirical check is not what §13.7 says.** §13.7 claims a scalar inverse problem isolating the
prior tail; the shipped `experiment6_heavytail.py` implements a *1-D random-walk Metropolis* on
`F_n(θ)=θ(1−n^{-2})` with only the Gaussian/Laplace/Student-t₃ priors (see §A). Table 9's ν=10 and ν=5
rows are unsupported.

---

## E. Numerical tables inconsistent with their own stated regressions

`verify_bt.py` recomputes every implied rate from the tabulated values:

| Table | Paper says | Implied by the table's own numbers |
|---|---|---|
| **3** (heat eq., n = 2…32) | "log–log regression slope **−2.98 ≈ −3**" | error column slope **−1.999**; increment column **−1.997** (every doubling gives ratio ≈ 4) |
| **4** (diffusion) | empirical 2.00 (theory 3) | **−1.999** ✓ self-consistent |
| **7** (2-D Poisson) | 1.99 / 1.98 | **−1.990** ✓ self-consistent |
| **8** (RK4 widths) | per-step "Rate" 0.58, 0.69, 0.92, 0.93; regression **−0.93** | per-step **0.32, 0.24, 0.15, 0.08**; regression **−0.199** |
| **8** (trapezoidal) | 0.89, 1.72, 2.05, 2.01; regression **−2.01** | **0.53, 0.72, 0.84, 0.92**; regression **−0.759** |

So Table 3, and the entire Table 8 (the empirical heart of Experiment 5), **contradict their own quoted
regressions**. Table 3's data actually exhibit the *same* Θ(n⁻²) rate as Tables 4 and 7 — not the
advertised −3. (The Tables 3/4/7 columns are also mutually inconsistent with the claimed *theory*
(k+1 = 3) while being consistent with rate 2 — worth reconciling: either the theory or the data is mislabelled.)

---

## F. Other mathematical / internal errors

1. **Lemma 3.10** — the final step is an inequality in the **wrong direction**:
   `Σ_{n≥N} n^{−α} ≤ 1/((α−1)(N−1)^{α−1})` is the true integral bound, but the paper then writes
   `≤ C/((α−1)N^{α−1})`, which is *smaller*. Script: at α=2, N=10 the sum is 0.10516 **>** the claimed
   bound 0.10000; at α=3, N=100, 5.05·10⁻⁵ > 5.00·10⁻⁵. Fix by bounding `Σ_{n≥N} n^{−α} ≤ N^{−α} + ∫_N^∞ t^{−α} dt`.
   (The conclusion O(N^{−(α−1)}) survives; only the displayed constant/step is wrong.)
2. **Lemma 12.1** — index mismatch. For spectral truncation, `F_{n+1}(x) − F_n(x)` is the **single mode**
   `j = n+1`, but the lemma writes `‖F_{n+1} − F_n‖² = Σ_{j>n} e^{−2λ_j T}|⟨x,e_j⟩|²`, which is
   `‖F − F_n‖²` (the full tail). The bound still yields a correct *upper* bound for the increment, but the
   stated equality is false; and Theorem 12.2's `O(e^{−cn²})` should be re-derived from the single-mode
   increment (it is dominated by the *largest* retained mode, so the true rate is `O(e^{−λ_{n+1}T})`).
3. **Telescoping Kalman Filter cost (§15.2)** — `CostTKF ∼ Σ_{t=1}^{T} n_t²` with `n_t = ⌈√t⌉` is
   `Σ t = Θ(T²)`, **not** `O(T^{3/2})`. Script: T = 10⁴ gives 5.07·10⁷ vs T^{3/2} = 10⁶. Since
   `CostKF = T·n_T² = Θ(T²)` too, the advertised asymptotic advantage of the TKF **disappears**.
4. **Corollary 7.2's proof** appeals to "Lemma A.1 and [14]" for quantile Lipschitz continuity, but
   Lemma A.1 itself *assumes* a density lower bound *and* log-concavity and derives quantile stability —
   it is not independent support. The `n⁻¹` "factor from the density lower bound" is asserted, not derived.
5. **Assumption 3.5 / §10 tension:** the paper says (a) Fernique gives exponential moments for *any*
   Gaussian prior, and (b) condition (5) with exponent *p* "holds for priors with sufficiently fast
   eigenvalue decay (Måtern s > d/2)". But (5) as written for general p is not implied by eigenvalue
   decay and needs its own justification; it is then used as if automatic.
6. **Notation clash:** §15.1 writes the degraded rate as `O(n^{−(k+1)·α})` with "α ∈ (0,1)", while
   Theorem 10.1 uses `ρ(α) = (α−2)/(α−1)`. Two different symbols (`α` vs `ρ`) for the same factor in
   the same paper.
7. **Experiment 4's "α = ∞" for the Laplace prior (§13.4/Remark 13.1)** conflates "exponential tails" with
   "α = ∞" (the Gaussian case). Laplace is not in the polynomial-tail family of (31) at all, so
   Theorem 10.1 simply does not apply; saying it "places it in the regime α = ∞" is incorrect.

---

## G. References (metadata-verified)

**Definitely wrong / mismatched:**

- **[37]** *"R. O'Brien, S. L. Cotter, and M. Dashti, Here Be Dragons …, SIAM/ASA J. Uncertain. Quantif. 14
  (2026), 1–25, https://doi.org/10.48550/arXiv.0909.2126."* — **arXiv:0909.2126 is "Approximation of
  Bayesian Inverse Problems for PDEs" (2009), not this work.** The genuine paper is
  *T. O'Brien, M. T. Moores, et al., "Here Be Dragons: Bimodal posteriors arise from numerical integration
  error in longitudinal models", arXiv:2502.11510 (Feb 2025)*. **Three of the four fields (authors, venue/
  year, identifier) are wrong.**
- **[38]** — DOI resolves to a different paper (see §D(ii)); title also wrong.
- **[39]** *"H. Owhadi, C. Scovel, and T. J. Sullivan, Bayesian weak seminorms, SIAM J. Numer. Anal. 53
  (2015), 1368–1390."* — no such title/venue found; the Owhadi–Scovel–Sullivan 2015 works are
  "Brittleness of Bayesian inference under finite information…" (EJS 10) and "On the Brittleness of
  Bayesian Inference" (SIAM Review 57). **Unverified / likely incorrect.**
- **[25]** *"… Stat. Comput. 34 (2024)"* — DOI `10.1007/s10444-024-10153-4` is **Advances in Computational
  Mathematics 50 (2024)**, not Statistics and Computing.

**Verified correct:** [16] Latz (SIAM Rev. 65, 831–865, 2023); [17] Capistrán & Christen et al.
(Bayesian Anal. 17, 2022 — note the paper lists only 2 of the 5 authors and reverses order);
[20] Lykkegaard–Dodwell–Scheichl (SIAM/ASA JUQ 11, 1–30, 2023); [22] Jasra–Law–Walton–Yang
(FCM 24, 1249–1304, 2023); [23] Sung–Tuo (WIREs CS 16, 2024); [24] **author list wrong** — the paper says
"Poot, van der Meer, Veraar" but Crossref gives **Poot, Kerfriden, Rocha, van der Meer** (no Veraar);
[26] Reiser–Lee–Rasmussen (Stat. Comput. 35, 2025).

**Self-citations [27]–[36]:** ten Zenodo entries by the same author. Given the paper under review itself
depends on several of these, they should be checked for existence/version dates; [27]
`10.5281/zenodo.17605158` resolves to the telescope-iterative-methods record. (Not exhaustively verified here.)

---

## H. Structural / editorial

- **Duplicated blocks.** There are **two** "Acknowledgments" + "Declaration of Competing Interest" +
  "Funding" sets (pp. 44 and 46), and **two** "Code and Data Availability" blocks (pp. 44 and 44) — an
  unmistakeable copy-paste artifact. The two acknowledgments even give different months for the Grok
  assistance (February 2026 vs January 2026).
- **Repository-tree cross-references are stale**: `experiment1_heat.py # Heat equation (Sec 10.2)`, etc.,
  but the experiments live in Section 13; and the tree lists 5 experiments for a paper with 7.
- **Table 6** (heavy-tail overhead) is placed inside the Experiment-5 (bimodality) section; it belongs in
  Section 10.
- **Repo contains a stray root `kalman_telescoping.py`** duplicating `src/kalman_telescoping.py`.

---

## What is sound (and worth keeping)

- The classical framework is fine: Hellinger definition/properties (§3.3), Fernique primer (§3.4),
  the potential/Radon–Nikodym setup (§3.2), and **Proposition 4.3's increment identity (12)** — I
  re-derived it and it is correct.
- **Corollary 4.4's bound** follows correctly from (12) + Cauchy–Schwarz + Assumption 3.11(c).
- **Lemma 5.1 and Theorem 5.3** (the telescoping decomposition) are a correct (if elementary) tautology,
  and the paper is commendably candid that it is an identity, not a new estimator.
- **Corollary 6.5's geometric-tail summation** (Σ_{n≥N} n^{−(k+1)} = Θ(N^{−k})) is correct.
- **Theorem 7.1** (variance-difference decay, given Theorem 6.1) is correctly assembled.
- **Corollary 10.2** (tempering restores the Gaussian rate) is right — `e^{−δ‖x‖²}` dominates any
  polynomial tail.
- The honesty of the positioning section (explicitly *not* claiming priority for standard rate derivation)
  is good practice and should be preserved.

---

## Recommended repairs (priority order)

1. **Withdraw or replace the reproducibility claim.** Either ship a complete, runnable implementation that
   actually produces Tables 3–11, or remove the "all code, data and scripts are available" statements and
   the repository tree. As shipped, the artifact cannot reproduce a single table.
2. **Re-do Experiment 5 from scratch.** It must solve the *actual* linear ODE inverse problem with RK4 and
   the trapezoidal rule, vary the step size, and exhibit (or not) discretisation-induced bimodality. The
   current script tests a hard-coded mixture prior and cannot speak to Theorem 9.3.
3. **Fix Theorem 9.3.** Include the `−⟨Γ⁻¹(y−F_n), ∇²F_n⟩` term in the Hessian (or prove a separate
   convexity argument). Remove the "odd/even order" heuristic (RK4 is 4th order); state the real criterion
   (non-injectivity of the numerical map). Correct every "RK4 … k = 1" to k = 4.
4. **Fix Theorem 8.4's algebra and thresholds**: `(β−γ)/2 + 1 = (β+2−γ)/2`; regimes split at
   **γ = β + 2**, not β + 1; case (iii) exponent `−2 − (β+2−γ)/(k+1)`. Then re-audit Theorem 11.2.
5. **Either prove Theorem 10.1 properly** (close the truncation argument, exhibit the lower bound) **or
   demote it** to a conjecture with the numerical observation. Correct the [38] citation to
   *Agapiou & Castillo, "Heavy-tailed Bayesian nonparametric adaptation", Ann. Statist. 52(4) 1433–1459
   (2024), 10.1214/24-AOS2397* and re-check the framing against that paper.
6. **Reconcile every table with its quoted regression** (Tables 3 and 8 at minimum), and fix the
   "oracle levels" (Tables 2/10) to match `⌈ε^{−1/3}⌉` or change the formula.
7. **Fix Lemma 3.10's integral step** and **Lemma 12.1's increment/tail** index error (and re-derive
   Theorem 12.2 from the single-mode increment).
8. **Fix the TKF cost** (`Θ(T²)`, not `O(T^{3/2})`) and drop the spurious advantage claim.
9. **Repair the bibliography** ([37], [38], [39], [24] author list, [25] venue, [17] authors) and
   **de-duplicate** the acknowledgments/declarations/funding/code-availability blocks.
10. Re-check the internal "α" vs "ρ" notation and the Laplace/"α = ∞" remark.

---

## Bottom line

The paper is well-organised and its *classical* sections are correct, but all three advertised
contributions fail on verification: one has an invalid proof plus a mislabelled example, one has an
unproven central formula resting on a mis-cited reference, and one has an off-by-one complexity analysis.
Worst of all, the artifact offered as the paper's empirical evidence is an 8 KB skeleton whose
"Experiment 5" hard-codes the very bimodality it claims to explain. This is a **major revision**; unless
the experiments can be genuinely run and Theorem 9.3's proof repaired, it is **not publishable as is**.
