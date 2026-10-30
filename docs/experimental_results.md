# 5.4 Experimental Results

*(Drop-in replacement/extension for Section 5.4. The original subsection described the experimental design only; the experiment has now been performed. All code, raw outputs, and figures are reproducible from `dp_simulation.py` with a fixed random seed.)*

## 5.4.1 Setup

We simulated the scenario specified in the design: a synthetic population of N = 1,000,000 users, of whom exactly 30% (π = 0.30) use the target feature. For every combination of privacy budget ε ∈ {0.1, 0.5, 1.0, 2.0} and sample size n ∈ {1,000; 10,000; 100,000}, we drew n users without replacement, collected one randomized-response report from each, computed the unbiased estimate π̂, and repeated the procedure 100 times. We report the mean absolute error (MAE), the median and 90th percentile of the absolute error, the empirical standard deviation of π̂, the average width of the 95% confidence interval, and the CI coverage (the fraction of runs whose interval contains π).

The randomized-response mechanism is the one described in Sections 3.7 and 4.1: each user answers truthfully with probability p and otherwise flips a fair coin. Under this mechanism ε = ln((1+p)/(1−p)), equivalently p = (e^ε − 1)/(e^ε + 1), and the estimator is the one given in Section 4.1, π̂ = (ŷ − (1−p)/2)/p, where ŷ is the observed fraction of "yes" reports. (See the erratum note on the ε–p relation printed in Section 4.1.) Confidence intervals use the normal approximation Var(π̂) ≈ ŷ(1−ŷ)/(n·p²).

Estimates were **not clipped** to [0, 1]: clipping would improve the error metrics of the weakest conditions at the cost of biasing the estimator, so we report the raw unbiased estimates and separately record how often they fall outside [0, 1].

## 5.4.2 Results

| ε | p | n | MAE | Median AE | P90 AE | Mean 95% CI width | Coverage | π̂ outside [0,1] |
|---|---|---|---|---|---|---|---|---|
| 0.1 | 0.050 | 1,000 | 0.2218 | 0.2002 | 0.4603 | 1.2399 | 0.96 | 14% |
| 0.1 | 0.050 | 10,000 | 0.0858 | 0.0771 | 0.1665 | 0.3922 | 0.96 | 0% |
| 0.1 | 0.050 | 100,000 | 0.0216 | 0.0181 | 0.0412 | 0.1240 | 0.97 | 0% |
| 0.5 | 0.245 | 1,000 | 0.0441 | 0.0368 | 0.0903 | 0.2517 | 0.97 | 0% |
| 0.5 | 0.245 | 10,000 | 0.0166 | 0.0143 | 0.0319 | 0.0796 | 0.96 | 0% |
| 0.5 | 0.245 | 100,000 | 0.0051 | 0.0042 | 0.0099 | 0.0252 | 0.98 | 0% |
| 1.0 | 0.462 | 1,000 | 0.0252 | 0.0207 | 0.0519 | 0.1316 | 0.95 | 0% |
| 1.0 | 0.462 | 10,000 | 0.0081 | 0.0077 | 0.0163 | 0.0417 | 0.99 | 0% |
| 1.0 | 0.462 | 100,000 | 0.0024 | 0.0020 | 0.0045 | 0.0132 | 0.97 | 0% |
| 2.0 | 0.762 | 1,000 | 0.0151 | 0.0131 | 0.0308 | 0.0775 | 0.96 | 0% |
| 2.0 | 0.762 | 10,000 | 0.0046 | 0.0039 | 0.0095 | 0.0245 | 0.97 | 0% |
| 2.0 | 0.762 | 100,000 | 0.0016 | 0.0014 | 0.0033 | 0.0078 | 0.95 | 0% |

**(1) Error decreases monotonically in ε — confirmed.** For every sample size, the MAE strictly decreases as ε grows (e.g., at n = 10,000: 0.0858 → 0.0166 → 0.0081 → 0.0046). The decrease tracks the 1/ε reference closely in the small-ε regime, where p = (e^ε − 1)/(e^ε + 1) ≈ ε/2, and flattens slightly for ε ≥ 1 because p grows sublinearly in ε — exactly the behavior the randomized-response variance formula predicts (Figure 2).

**(2) Error decreases at the √n rate — confirmed.** Fitting a log-log slope of MAE against n gives −0.51, −0.47, −0.51, and −0.48 for ε = 0.1, 0.5, 1.0, and 2.0 respectively, in close agreement with the theoretical −1/2 (Figure 1). Quadrupling the sample roughly halves the error, as claimed in Section 5.3.

**(3) Large samples make small budgets usable — confirmed with one qualification.** At n = 100,000, even ε = 0.1 attains a typical error of about 2 percentage points (MAE 0.0216, P90 0.0412), inside a ±5-point tolerance, and ε = 0.5 reaches sub-percentage-point accuracy (MAE 0.0051). The qualification is interval width: at ε = 0.1 the average 95% CI half-width is ±6.2 points even at n = 100,000, so analyses that must act on the interval rather than the point estimate still need either a larger budget or more samples. The opposite corner of the grid quantifies the warning in Section 5.3: at ε = 0.1 with n = 1,000 the mechanism is effectively unusable — the average CI is wider than the entire [0, 1] range and 14% of the raw estimates fall outside [0, 1].

**Coverage and theory agreement.** Empirical CI coverage stays at 0.95–0.99 across all twelve conditions, consistent with the nominal 95% level. The empirical standard deviation of π̂ matches the theoretical standard error computed from q̄(1−q̄)/(n·p²), where q̄ = π(1+p)/2 + (1−π)(1−p)/2, to within sampling noise (100 runs). A practical note: the mixture-only variance formula, [q₁(1−q₁)π + q₀(1−q₀)(1−π)]/(n·p²), underestimates the total variance when ε is large (at ε = 2.0 it predicts SE 0.0014 against an observed 0.0019–0.0020), because it omits the variance of the sample's true composition; the q̄(1−q̄) form includes it and should be preferred.

**Comparison with the illustrative table in Section 5.2 (n = 10,000).** The measured errors are consistent with the table's qualitative bands: ε = 0.5 gives typical errors of 1.4–3.2 points (median–P90), matching "moderate (±1–3%)"; ε = 1.0 gives 0.8–1.6 points, matching "small (about ±1%)" as a statement about typical error (the 95% CI half-width is ±2.1 points); ε = 0.1 gives an MAE of 8.6 points, matching "large (several % or more)." The table's values can now be cited as measured rather than illustrative.

## 5.4.3 Supplementary Experiment: Laplace Mechanism for a Bounded Mean

To check the standard-error claim of Section 4.2, we ran the same grid for the Laplace mechanism: per-user values on [a, b] = [0, 1] drawn from a fixed synthetic population with true mean 0.300, each user reporting their value plus Laplace noise of scale (b−a)/ε, and the server averaging the noisy reports.

The measured standard deviation of the estimate matches √(σ²_pop/n + 2(b−a)²/(ε²n)), which for small ε is dominated by the noise term (b−a)·√2/(ε√n). For example, at ε = 0.1 and n = 1,000 the empirical SD is 0.424 against 0.447 predicted, while the expression in Section 4.2, SE ≈ (b−a)/(ε√n), predicts 0.316. The approximation in Section 4.2 is therefore correct in form (∝ 1/(ε√n)) but understates the constant by the factor √2, the standard deviation of a unit-scale Laplace variable; we recommend restating it as SE ≈ (b−a)·√2/(ε·√n).

## 5.4.4 Reproducibility

All results were produced by a single script (`dp_simulation.py`, NumPy/pandas/Matplotlib, fixed seed 20260604): 100 repetitions per condition, sampling without replacement from the fixed population, and no post-hoc selection of runs. Condition-level aggregates are in `results.csv` and `laplace_results.csv`; Figures 1 and 2 are `fig1_mae_vs_n.png` and `fig2_mae_vs_eps.png`.

---

## Proposed revision to Section 7, Limitation 3

> **Before:** This study focuses on design and qualitative/quantitative analysis. A large-scale implementation under real traffic and a load evaluation remain future work.
>
> **After:** The proposed trade-off analysis has been validated by Monte Carlo simulation (Section 5.4): the predicted 1/ε and 1/√n error behavior, the confidence-interval coverage, and the illustrative accuracy bands of Section 5.2 were all confirmed on a synthetic population of one million users. A large-scale implementation under real traffic — including client-side deployment, the budget-manager service, and load evaluation — remains future work, as does validation on real behavioral data, whose correlation structure a synthetic Bernoulli population does not capture.

---

## Erratum for Section 4.1 (draft)

> **Erratum.** Section 4.1 states the relation ε = ln(p/(1−p)); this relation belongs to the direct-perturbation variant of randomized response (truthful with probability p, inverted otherwise), not to the mechanism described in this paper. For the forced-response mechanism used here (truthful with probability p, otherwise a fair coin flip) — which is the mechanism assumed by the estimator π̂ = (ŷ − (1−p)/2)/p and by the scenario in Section 3.7 — the correct relation is ε = ln((1+p)/(1−p)), equivalently p = (e^ε − 1)/(e^ε + 1). All analyses and experiments in this paper use the forced-response mechanism with this corrected relation.
