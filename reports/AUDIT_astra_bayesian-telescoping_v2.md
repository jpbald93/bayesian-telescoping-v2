# Bayesian Telescoping V2 — adversarial audit

## Verdict

**MAJOR REVISION. V2 makes several correct and necessary retractions, but is not mathematically sound as submitted.** Its corrected power-sum threshold, affine-Gaussian posterior theorem, heavy-tail numerical rates, TKF cost, and principal citation corrections are sound. However, the replacement multimodality account contains new false statements, **Table 1 suppresses a real RK4 mode at n=2**, and the retained credible-interval proposition is false. The MLMC exponent is correct conditional on the stipulated cutoff, but that cutoff does not follow from the paper's own bias bound.

I ran the supplied driver unchanged to completion: **exit status 0; `ALL ASSERTIONS PASS`**. This does not validate the manuscript: its mode finder discards maxima more than 14 log-density units below the largest mode and scans only [-8,6]. Independent derivative-root finding and adaptive quadrature expose what those tests miss.

**Audit artifacts, all in `Math Audit/reports/`:**
- `astra_v2_driver_output.txt`: my fresh run of the author's driver.
- `astra_v2_independent_checks.py`: my independent NumPy/SciPy/SymPy/mpmath checks; does not import the driver.
- `astra_v2_independent_output.txt`: successful independent run, including root locations, quadrature, rates, symbolic counterexamples, and cost sums.

The independent posterior calculation uses analytic log-density derivatives with bracketing, followed by adaptive integration, not the author's grid-mode or rectangular-quadrature code. Heavy-tail means use infinite-domain adaptive quadrature after centering and scaling by the likelihood. Symbolic calculations supply an explicit counterexample to the quantile proposition.

## Claim-by-claim

| Item | Finding |
|---|---|
| Lemma 2.1 upper tail inequality | **CORRECT.** V1's reversed comparison is correctly identified. |
| Tail asymptotic | **CORRECT for the tail**, but the displayed concluding sum in V2 is `sum_{n<N}`, which is false; it must be `sum_{n>=N}`. |
| Theorem 4.1 threshold gamma versus beta+2 | **CORRECT** under the displayed continuous-allocation model. |
| Case-(iii) exponent | **CORRECT conditional on L ~ epsilon^(-1/(k+1))**; not the RMSE complexity implied by the preceding O(L^-k) bias bound. |
| Theorem 5.1 affine approximate maps + Gaussian prior/noise | **CORRECT**, with the usual existence/nondegeneracy qualifications. Its extrapolation to arbitrary discretizations of a linear exact map is too broad. |
| RK4 critical point and minimum | **CORRECT:** z* = -1.5960716379833215; R(z*) = 0.2703947652051846. |
| RK4 is bimodal only at n=1 in Table 1 | **FALSE as a count of mathematical local modes.** There are two at n=2, one extremely small. |
| RK4 fold criterion / eventual global injectivity | **FALSE as written.** Missing n, wrong branch inequality, and pointwise convergence mistaken for global injectivity. |
| T'(z)>0 on z<2 | **CORRECT for the one-step rational map.** The conclusions about T(lambda/n)^n and the full Gaussian-prior posterior do not follow. |
| Heavy-tail slopes -2.988 for all five priors | **CORRECT for the tested model and regression window.** Independently reproduced. |
| Proposition 6.1 sufficient exponential-moment condition | **CORRECT given the appropriate polynomial potential-increment bound**, eventually in n; the preceding corollary needs repair. |
| V1 heavy-tail formula unsupported / V1 Table 9 not reproduced | **CORRECT.** Faster observed rates refute a compulsory/sharp penalty, not merely a weaker big-O upper bound. |
| TKF sum Theta(T^2) | **CORRECT.** |
| O'Brien–Moores and Agapiou–Castillo citation corrections | **CORRECT identifiers and principal metadata.** |
| Credible-interval Proposition 3.5 | **FALSE.** Two moments do not control quantiles, even for uniformly positive strictly log-concave densities. |

## The MLMC threshold

Starting from independent level samples, the variance constraint is sum V_n/N_n <= epsilon^2. Minimizing sum N_n C_n with real positive N_n gives

`N_n = epsilon^(-2) sqrt(V_n/C_n) S_L`, where `S_L = sum_{n<L} sqrt(V_n C_n)`.

Consequently the optimal relaxed sampling cost is `epsilon^(-2) S_L^2`. For two-sided power laws, `S_L ~ sum n^s`, with `s=(beta-gamma)/2`. Thus:

- gamma > beta+2: S_L stays bounded;
- gamma = beta+2: S_L = Theta(log L);
- gamma < beta+2: S_L = Theta(L^((beta+2-gamma)/2)).

These are the **right corrections** to V1. Independent partial sums are:

| beta=1; gamma | L=10^3 | L=10^4 | L=10^5 | L=10^6 |
|---|---:|---:|---:|---:|
| 2.5 | 19.04956 | 36.55821 | 67.68980 | 123.04981 |
| 3 | 7.48447 | 9.78751 | 12.09014 | 14.39273 |
| 4 | 2.54911 | 2.59237 | 2.60605 | 2.61038 |

The first grows asymptotically as 4 L^(1/4); the last converges to zeta(3/2). Squaring and substituting the manuscript's **assumed** `L ~ epsilon^(-1/(k+1))` gives its case-(iii) exponent exactly.

**Unrepaired bias mismatch.** Theorem 3.3 bounds consecutive increments by O(n^(-(k+1))) and the remaining bias by O(N^-k), correctly losing one power on summation. To guarantee bias O(epsilon), that result requires `L ~ epsilon^(-1/k)`, not `epsilon^(-1/(k+1))`. With only that bias result, the case-(iii) complexity becomes

`O(epsilon^(-2-(beta+2-gamma)/k))`.

A direct O(L^(-(k+1))) bias assumption or a genuinely faster-bias extrapolated estimator would justify the stated cutoff, but neither is supplied. This issue does **not** invalidate the corrected summability threshold.

**Two other qualifications.**
1. The text starts with V_n=O(n^-gamma), C_n=O(n^beta), then switches to two-sided asymptotic comparability. Upper bounds justify upper cost bounds, not exact regime necessity or equivalent power laws for the actual V_n C_n. State Theta assumptions or phrase everything as bounds.
2. Integer sample allocations generally add `sum_{n<L} C_n` to the relaxed bound. Requiring at least one sample per level costs Theta(L^(beta+1)) when C_n ~ n^beta and beta>-1. For example, with k=1, beta=4, gamma=7 and the manuscript's own cutoff, this floor costs Theta(epsilon^-2.5), despite its nominal O(epsilon^-2) first regime. Either explicitly restrict the theorem to the continuous-allocation cost model or include the floor and conditions under which it is negligible.

The paragraph calling n^-3/4 the “V1 boundary” also mislabels its example: for beta=1 the V1 boundary is gamma=2, giving n^-1/2. Gamma=2.5 lies inside V1's incorrectly claimed bounded-sum regime.

## Bimodality

### What the affine theorem actually proves

For F_n(x)=A_n x+b_n in finite dimensions, completing the square gives

`Sigma_n^(-1)=Sigma_0^(-1)+A_n^T Gamma^(-1) A_n`,

`m_n=Sigma_n[ Sigma_0^(-1)m_0 + A_n^T Gamma^(-1)(y-b_n) ]`.

This is Gaussian and, for positive-definite covariance, strictly log-concave with a unique mode. **The stated affine hypothesis and proof are sound.** A linear ODE in its state is not necessarily a linear map of its unknown parameter; V2 correctly explains that distinction.

But “no choice of discretisation of a linear forward operator” must read **“no affine discretisation, with Gaussian prior and noise.”** An arbitrary consistent approximant to a linear exact map need not be affine. For example, F(x)=x and F_n(x)=x+x^2/n^2 converge pointwise and satisfy polynomial increment bounds. With prior N(0,1), y=0 and sigma=.05, the n=1 posterior log-score is

`-x(800 x^2 +1200 x +401)`.

It has maxima at approximately -0.9974873734 and 0, separated by a minimum at -0.5025126266. This is a counterexample to the overbroad sentence, not to the affine theorem.

### RK4: correct constants, incorrect domain logic

Independent high-precision root solving gives

`z* = -1.59607163798332152311280541440`,

`R(z*) = 0.270394765205184607962459613383`.

Moreover `R''(z)=1+z+z^2/2=((z+1)^2+1)/2>0`. Therefore R' has exactly one real zero, R is positive, and for every finite n

`d/dlambda R(lambda T/n)^n = T R(lambda T/n)^(n-1) R'(lambda T/n)`.

The n-step map decreases to its minimum at `lambda=n z*/T` and then increases. It is injective separately on **both** branches, not “only” on the branch lambda h<=z*. It remains non-injective on the entire real line, and even on the entire negative half-line, **for every finite n**. The fold moves left; it does not cease to exist.

For a restricted parameter interval [-B,0], an appropriate sufficient no-fold condition is

`B T/n <= |z*|`, equivalently `n >= B T/|z*|`.

V2's displayed `|lambda| T <= |z*|` omits n and is not a criterion for injectivity on the prior support. Its own n=1 example has the true lambda=-1 satisfying that inequality, yet is bimodal because the prior also covers the other branch. A full-support Gaussian has no finite B.

Non-injectivity does not automatically make the posterior bimodal. The likelihood has two maxima only if the datum is above the numerical minimum R(z*)^n; multiplication by a prior can remove a secondary maximum. Conversely, monotonicity of a nonlinear map alone does not imply posterior unimodality: with the strictly increasing exact map exp(lambda), prior N(-8,1), y=exp(-1) and sigma=.05, my score solver finds stationary points -7.9480539921, -3.3564273771, -1.1605794313, the first and last being maxima.

### Independent mode calculation: Table 1 is not literally correct

I solved the posterior score without the author's height threshold. For the stated Gaussian prior N(-1,.5^2), sigma=.05 and y=exp(-1):

| n | RK4 local maxima | Secondary peak's log-height relative to principal peak |
|---:|---|---:|
| 1 | -2.016036888839; -1.019817987897 | -2.232828220 |
| 2 | **-4.713878692266**; -1.000740067526 | **-29.316283331** |
| 4 | -1.000037362728 | — |
| 8 | -1.000002102853 | — |
| 16 | -1.000000124743 | — |
| 64 | -1.000000000469 | — |

At n=2 the intervening local minimum is -4.042113236371. The secondary peak is inside the driver's scan interval but fails its explicit `lp[i] > -14` test. Its relative height is about 1.85e-13. Independent adaptive quadrature gives probability left of the intervening minimum of **2.47567e-13**, versus **0.0941099** at n=1. It is negligible in practice, but mathematically a mode.

Thus **“a practically important coarse-grid transient” is defensible for this example**. “Bimodal only at n=1,” “disappears for n>=2,” and the literal mode count are false unless an explicit significant-mode definition replaces the ordinary one. The sampled fine levels support absence of a secondary RK4 posterior maximum there; they do not prove a general eventual-unimodality theorem.

### Trapezoidal: one-step derivative does not settle the n-step map

`T'(z)=1/(1-z/2)^2>0` is correct on each side of the pole. But

`d/dlambda T(lambda Tobs/n)^n = Tobs T(lambda Tobs/n)^(n-1) T'(lambda Tobs/n)`.

For even n, it changes sign at lambda Tobs/n=-2: the map folds on z<2 because T(z)<0 for z<-2. The negative-parameter no-fold interval is z in [-2,0], not all z<2. For example, n=2 and the stated datum give two likelihood roots, approximately -0.9797 and -16.33195; the strong prior can suppress the far likelihood branch as a posterior mode.

There is a separate support issue: the stated Gaussian prior also covers the region beyond the pole lambda=2n/Tobs. The likelihood tends to zero at the pole and the posterior tends to zero at positive infinity, with a positive continuous density between. Hence there is a local maximum to the right of that pole. My search finds the following additional trapezoidal posterior maxima:

| n | Far-right maximum | Relative log-height |
|---:|---:|---:|
| 1 | 10.160494620741 | -939.569 |
| 2 | 19.302009857044 | -1585.671 |
| 4 | 46.030222489931 | -7170.896 |
| 8 | 120.738813512164 | -42668.163 |
| 16 | 342.960242112399 | -313504.733 |
| 64 | 3284.157618351844 | -25892432.828 |

These are utterly negligible but refute “one mathematical mode” with an unrestricted Gaussian prior. Assigning any value to the map at the single pole does not change the posterior measure. If the author intends lambda<0 or another stable parameter domain, specify and use the corresponding **truncated** prior. That removes the positive-pole objection, but not the even-n negative-branch issue or the RK4 n=2 maximum.

## Heavy tails

### The numerical correction is right

For F_n(theta)=a_n theta, a_n=1-n^-2, y=1.5 and sigma=.1, independent infinite-domain adaptive quadrature gives these slopes of log consecutive differences versus log n, using n=60,...,238:

| Prior | Independent slope |
|---|---:|
| Gaussian | -2.988031713 |
| t_10 | -2.988035585 |
| t_5 | -2.988037506 |
| t_3 | -2.988038823 |
| t_2 | -2.988039627 |

All round to **-2.988**, as reported. For the Gaussian prior the exact mean is `a_n y/(a_n^2+sigma^2)`; it agrees with my quadrature to floating-point accuracy.

There is a straightforward explanation beyond the regression. For a near 1, the Gaussian likelihood makes the numerator and denominator defining the posterior mean smooth functions of a, including for Student-t priors. If m'(1) is nonzero, as it is in these examples, `m(a_{n+1})-m(a_n) ~ 2 m'(1)n^-3`. The tail of the prior does not impose the proposed penalty in this identifiable linear Gaussian-likelihood model. The first finite-n correction explains why the fitted slope is slightly above -3.

The conclusion **“V1's experiment does not support its reported compulsory/sharp heavy-tail penalty” is correct**. A faster n^-3 rate still satisfies a weaker upper bound O(n^-r), r<3; therefore these observations do not by themselves disprove such an upper bound or a worst-case rate over a larger model class. They do refute claiming that this model exhibits that degradation.

### Proposition 6.1 is a valid sufficient condition, with explicit qualifications

Assume the actual bound `|Delta Phi_n| <= C(1+||x||^(2p)) n^(-(k+1))` and `sup_{n>=n0} E_mu_n exp(a ||x||^(2p)) < infinity` for a>0. For sufficiently large n, C n^(-(k+1)) <= a/2. Absorb the polynomial prefactor into the remaining half of the exponential moment. Then

`E_mu_n |exp(-Delta Phi_n)-1| = O(n^(-(k+1)))`.

Also `Z_{n+1}/Z_n = E_mu_n exp(-Delta Phi_n) = 1+O(n^(-(k+1)))`; inversion gives the ratio used in V2. Normalizing the likelihood ratio yields the same total-variation/ bounded-test-function rate. This repairs the proof's implicit large-n step. Finite initial levels do not matter for the asymptotic conclusion.

**Corrections needed around that result:**
- “Indeed for any non-Gaussian tail” is false. For p=1, a non-Gaussian density proportional to exp(-x^4) satisfies the condition. More generally many sufficiently light non-Gaussian tails do.
- A heavy-tailed **prior** need not yield a heavy-tailed **posterior**. For this example and n>=2, a_n>=3/4; Gaussian likelihood suppression gives a uniform quadratic exponential moment for all the tested Student-t posteriors. At n=1, F_1=0 and the posterior equals the prior, so the manuscript's strict supremum over n>=1 fails for Student-t priors. An eventual-n condition is the natural formulation.
- The tested observable theta is unbounded, whereas Theorem 3.3 covers bounded Lipschitz observables. The experiment is still valid, but its justification is Gaussian likelihood domination, not automatic application of that theorem.
- Define the tail index consistently. V1 defined a survival-tail index, for which t_nu has alpha=nu. V2 instead writes tails “~||x||^-alpha,” ambiguously suggesting density tails; a t_nu density decays as |x|^(-(nu+1)). Also t_2 is a boundary illustration, outside V1's alpha>2 hypothesis.
- The conjecture's use of Theta for arbitrary observables is inappropriate: a constant observable has zero increments, and symmetry can cause cancellation. Formulate degradation as failure of a specified O bound or as a worst-case statement with nondegeneracy assumptions. Its proposed necessary moment failure remains explicitly speculative, not established.

## New-error hunt

### 1. False quantile proposition — substantive

Proposition 3.5 claims quantile convergence at the rate of the first two moments, assuming unimodality and a uniform local density lower bound. These hypotheses are insufficient, even with strict log-concavity.

Here is an explicit counterexample on [-1,1]:

`p_t(x)=3/5-(3/10)x^2+t P_3(x)`, where `P_3(x)=(5x^3-3x)/2`, and `t_n=.01/n`.

Orthogonality gives, exactly,

`integral p_t=1`, `E_t X=0`, `E_t X^2=7/25`, independently of t.

For |t|<=.01, p_t>=.29 and `p_t''=-.6+15tx<0`. Hence each density is strictly concave, unimodal, and strictly log-concave on its support. Nevertheless

`F_t(0)=1/2+t/8`,

so its median is `q(t)=-(5/24)t+O(t^2)`. With t_n=.01/n, medians differ from the limiting median by Theta(1/n), although the moment errors are identically zero and hence O(n^-2). Numerically n q_n tends to -0.0020833333333.

The required repair is **uniform CDF control** `sup_s |F_n(s)-F(s)|=O(delta_n)` plus a density lower bound near the limiting quantile; total variation of that order suffices. The exponential-moment proof can potentially supply total variation directly. Convergence of two moments cannot replace it. Log-concavity also does not automatically provide a uniform-in-n lower density bound near a fixed quantile.

### 2. Incorrect polynomial power in Corollary 3.2

From the correct increment identity and growth `||F_n(x)|| <= C(1+||x||^p)`, the cross term has degree **2p**, not p. V2 puts the 2p power only on the smaller n^(-2(k+1)) term.

Use F_n(x)=(1-n^-2)x, y=0, Gamma=1, p=1, k=2. Then

`|Delta Phi_n| = .5[(1-(n+1)^-2)^2-(1-n^-2)^2] x^2 ~ 2x^2 n^-3`.

The printed bound would control this by `C(1+|x|)n^-3 + C(1+x^2)n^-6`, impossible uniformly in x,n. Setting x=n^3 makes the ratio to that expression with C=1 diverge like n^3. Replace the first polynomial factor by 1+||x||^(2p), allowing appropriate y-dependent constants. This repaired estimate is exactly what Proposition 6.1's proof uses.

### 3. Hellinger normalization is inconsistent

The definition `d_H^2=(1/2) integral (sqrt(p)-sqrt(q))^2` implies

`d_H=sqrt(1-rho)` and `0<=d_H<=1`,

not `sqrt(2(1-rho))`. The stated upper bound sqrt(2) is a loose bound rather than a false inequality, but it signals the same normalization mix-up. For disjoint unit masses, the defined distance is 1, while the claimed formula gives sqrt(2). Choose one convention consistently.

### 4. Tail lemma prints the wrong sum

The correct sandwich is

`N^(1-alpha)/(alpha-1) <= sum_{n>=N} n^-alpha <= N^-alpha + N^(1-alpha)/(alpha-1)`.

This proves the intended Theta result. But the lemma actually concludes `sum_{n<N} n^-alpha = Theta(N^(1-alpha))`. Its left side tends to zeta(alpha)>0 and is Theta(1), not a vanishing quantity. This is a small edit with a real mathematical consequence in the current statement. At alpha=2,N=10 the actual tail is 0.10516633568168575, between 0.1 and 0.11, independently confirming the V1 correction.

### 5. Hessian is attributed to the wrong density

V2 calls its expression with `+Prior^(-1)` the Hessian of `-log(dmu_n/dmu_0)`. That derivative is the likelihood potential plus a constant, so **no prior precision belongs in it**. The displayed formula, including prior precision, is instead the Hessian of the negative log posterior **Lebesgue density**, for a Gaussian prior in finite dimensions. The warning about the missing residual-times-second-derivative term is correct; fix the measure being differentiated.

### 6. Some assumptions are referred to but not stated

Theorem 3.3 invokes “assumed growth” without specifying it, and convergence to the named exact F needs an explicit assumption F_n -> F. Summable increments alone converge to some limit, not necessarily the named exact model. State k>0, uniform growth, consistency, and the norm used for the observation residual.

For the advertised infinite-dimensional Hilbert observation setting, Gamma^(-1/2) need not be defined on a generic Gaussian observation: Gaussian noise is generally not in its Cameron–Martin space. The literal squared-residual potential therefore needs finite-dimensional observations or a properly formulated infinite-dimensional likelihood. This is a scope/assumption issue, not an objection to the finite-dimensional computations.

### 7. The verification claim exceeds what the driver does

The driver passes, but it does not establish every mathematical statement or even reproduce every displayed number:
- mode counts include an undisclosed relative-height cutoff and a finite search interval;
- `rho if rho else ...` prints **-3.000 for t_2**, whereas the manuscript's displayed formal prediction is **0.000**;
- **no TKF sum is computed in the shipped driver**, despite the sentence saying that it reports that sum;
- it tests the tail upper bound using finite truncated sums and does not assert V1's failed lower comparison as claimed;
- testing whether a sum at one L is below 50 is not a proof of asymptotic convergence;
- no credible-interval proposition, general injectivity criterion, or sufficient-moment proof is checked.

Keep the passing driver, but describe it as selected numerical checks, remove the blanket “all mathematical statements ... verified against the driver,” and test analytic mode derivatives without silently excluding genuine extrema.

## What is sound

1. **The V1 tail-integral direction error is real**, and V2's new upper bound is correct.
2. **The increment identity itself is exact.**
3. **Full bounded-expectation increment rate and variance-difference rate** follow from the corrected polynomial estimate and stated uniform exponential moments. For bounded Lipschitz phi, phi^2 is also bounded Lipschitz, so the variance result follows immediately.
4. **The beta+2 MLMC threshold is the right correction**, not beta+1. The conditional case-(iii) algebra is right.
5. **Affine Gaussian conditioning is Gaussian.** That obstruction is a useful correction when kept within its hypotheses.
6. **The RK4 critical point and minimum are correct**, and the n=1 significant bimodality is independently reproduced. The sizeable spurious mass becomes negligible by n=2 in this example.
7. **All five heavy-tail experimental slopes are correct**, and withdrawing the unsupported universal/sharp penalty is warranted.
8. **TKF cost is Theta(T^2).** Since t <= ceil(sqrt(t))^2 <= t+2sqrt(t)+1, the sum is `T(T+1)/2+O(T^(3/2))`. My exact value at T=10,000 is **50,666,650**, matching the manuscript's rounded 5.07e7. The final-level fixed schedule has the same quadratic order.
9. **Citation corrections checked independently:**
   - [arXiv:2502.11510](https://arxiv.org/abs/2502.11510): *Here Be Dragons: Bimodal posteriors arise from numerical integration error in longitudinal models*, **Tess O'Brien, Matthew T. Moores, David Warton, Daniel Falster**; first submitted 2025. V2's shortened author list and identifier are right.
   - [10.1214/24-AOS2397](https://doi.org/10.1214/24-AOS2397): **Sergios Agapiou and Ismaël Castillo**, *Heavy-tailed Bayesian nonparametric adaptation*, **Annals of Statistics 52(4), 1433–1459 (2024)**. Crossref confirms title/authors/volume/issue; Project Euclid search metadata and Castillo's publication page confirm pages.
   - Crossref confirms **10.1214/24-AOS2356** instead belongs to *Early stopping for L2-boosting in high-dimensional linear models*, by **Bernhard Stankewitz**. V2's “Stankewitz et al.” should be singular; the central DOI mismatch diagnosis is correct.

## Recommended repairs

1. **Replace Proposition 3.5:** use CDF/total-variation convergence plus a local density lower bound; remove the two-moment argument.
2. **Correct and rescope Section 5.2:** state the two RK4 monotonicity branches, retain the factor n, distinguish bounded-support injectivity from pointwise consistency, and do not assert global disappearance of the fold.
3. **Repair Table 1:** either count all mathematical maxima (including RK4 n=2 and the trapezoidal positive-pole branch under the actual Gaussian prior), or explicitly report significant modes with a declared threshold and domain. Distinguish mathematical modes from practically relevant posterior mass.
4. **Restrict the affine obstruction to affine approximate maps**, Gaussian prior and Gaussian noise. Do not infer posterior unimodality from injectivity of an arbitrary nonlinear map.
5. **Separate the weak-bias exponent from the consecutive-increment exponent in MLMC.** Use a symbol q for bias O(L^-q), then report case (iii) with denominator q. Here the preceding theorem supplies q=k. Address integer allocation or clearly label the relaxed model.
6. **Fix Corollary 3.2's polynomial power, the tail index typo, Hellinger normalization, and Hessian reference measure.** State the missing consistency and growth hypotheses.
7. **Keep the independently confirmed heavy-tail experiment**, but clarify prior versus posterior tails, eventual exponential moments, the unbounded observable, and what faster convergence does and does not refute. Remove the assertion that all non-Gaussian tails fail exponential integrability.
8. **Amend the driver/reproducibility claims:** remove the false t_2 printout, compute the TKF sum if claimed, and add counterexample-sensitive tests. A passing driver is evidence only for what it actually tests.

**Bottom line:** V2 is directionally a genuine correction of V1, particularly on MLMC power sums, TKF accounting, and the absent heavy-tail penalty. But its replacement bimodality claims and retained quantile proposition need substantive repair before the paper can fairly describe the corrected picture as sound.
