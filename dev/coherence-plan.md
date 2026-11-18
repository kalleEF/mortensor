# Coherent forecasting: ideas from bamp in mortensor

Plan of 2026-11-18, after the review of the bamp development branch
"coherent_forecasting" (commit 4363a78 of 2026-07-13) and its design note
`docs/coherent-forecasting.md`. (A first reading used commit c4c3678 of
2026-06-01 of the same lineage; the findings were checked again against
4363a78.) The review for its author lists the findings on the branch code
itself. This plan takes (A) the branch's models into mortensor and (B) the
ideas of the design note into the coherent methods mortensor already has.

## Rules

- Each method is a published method, or its help page and the
  methods-sources article label it "package-specific (reference unknown)".
- New methods are opt-in; defaults and existing results stay as they are.
- Tests check new methods against the published method or a reference
  implementation, and against exact identities where they hold.
- An item enters only when it adds something mortensor lacks.

## Decisions of the owner (2026-11-18)

1. The bamp branch is not released and its release date is unknown, so its
   models are rewritten in mortensor (no wrapper). Each rewrite carries a
   note: when a bamp release contains the model, mortensor calls bamp and
   deprecates its own sampler.
2. The multi-cause model is rewritten as well, labelled package-specific
   (no published source for the combination was found).
3. The coherent model allows three or more strata, labelled
   package-specific beyond the two populations of Cairns et al. (2011).
4. Correlated deviations in `li_lee()`: the audit decides (B3). Result
   (`dev/audit-li-lee-innovations.md`): no change; joint innovations leave
   the scores unchanged.

## Work packages

| Package | Items | Content |
|---|---|---|
| C1 | B1, B2, A4, B4 | `mt_backtest(multivariate = TRUE)`: energy and variogram scores of all strata jointly; `mt_tune()` and `tuned()` with gap and joint targets; FAQ on rates and probabilities for microsimulation; table of the weights of the common trend of each coherent model |
| C2 | B3 | Audit of the correlation of the stratum deviations of `li_lee()` on `ons_data` and an HMD extract; option for joint innovations if the audit shows a gain |
| C3 | B5 | Sex by cause forecasts from existing tools in the coherence vignette |
| C4 | A1 | `apc_bayes(coherent = TRUE)`: native rewrite of `bamp_coherent()` with the fixes of the review |
| C5 | A2 | Native rewrite of `bamp_multicause()` |

### A1. Coherent joint APC model (C4)

`apc_bayes(coherent = TRUE, common_by = , reference = , cohort_deviation =
)`: shared period and cohort random walks, stratum age effects and levels,
stationary AR(1) deviations of the strata from the shared period effect
with a sampled coefficient, for S >= 2 strata per group.

Published anchor: Cairns et al. (2011), two populations, Bayesian
single-stage fit, common period trend and mean-reverting AR(1) spread with
level zero; Section 4.5 (equal populations: the spread sums to zero, the
default) and Section 4.4 (dominant population: `reference`). Differences,
labelled package-specific: binomial-logit likelihood with bamp's
random-walk smoothing priors (the paper: Poisson, random walk with drift,
AR(2) cohort processes); one shared cohort effect (the paper's alternative
hypothesis, Section 4.3) with optional AR(1) cohort deviations; independent
innovations of trend and deviations; priors; S >= 3.

Fixes of the review in the rewrite: the same constraint for every S (the
branch: a time constraint for S = 2, the sum over strata for S >= 3); rho
sampled for every S (the branch: fixed for S >= 3); the period slope
pinned only when all effects are RW2; exposures not rounded; an exact
Metropolis-Hastings sampler (the branch approximates the Polya-Gamma draws
by normal draws). Forecast innovations use the scaled precision of the
prior, as the branch does since commit 08c7a18.

Tests: constraints hold in every draw; the AR(1) and random-walk forecast
variances match their formulas; the sampler recovers simulated parameters;
with fixed hyperparameters, its posterior means agree with importance
sampling from the Gaussian approximation. `dev/check-coherent-apc.R`
compares the posterior with the branch implementation.

### A2. Multi-cause model (C5)

Native rewrite of `bamp_multicause()`: the deaths of each cell split by a
stick-breaking sequence of binomial-logit APC models of the cause shares
(Linderman, Johnson and Adams 2015 for the representation), period (and
cohort) effects of the shares coupled across causes by a Wishart prior,
all-cause rates from an APC model of the total, cause rates as shares
times the total. Package-specific (reference unknown).

### Not taken over

- `reconcile_apc()`: reconciles means and quantiles separately; mortensor
  reconciles each draw (`reconciled()`, `mt_reconcile()`).
- `bamp_strata()`: shares of events without exposures; mortensor covers
  totals and shares on the rate scale (`coda_causes()`, `relative_model()`,
  `mt_reconcile(method = "proportional")`).
- Forecast options of the branch (damped trends, variance damping, trend
  windows, changepoints, exogenous age drifts): unpublished extensions of
  the forecast, outside this plan.

## Sources

- Cairns AJG, Blake D, Dowd K, Coughlan GD, Khalaf-Allah M (2011). Bayesian
  stochastic mortality modelling for two populations. ASTIN Bulletin
  41(1): 29-59. doi:10.2143/AST.41.1.2084385
- Gneiting T, Raftery AE (2007). Strictly proper scoring rules, prediction,
  and estimation. JASA 102(477): 359-378. doi:10.1198/016214506000001437
- Linderman SW, Johnson MJ, Adams RP (2015). Dependent multinomial models
  made easy: stick-breaking with the Polya-gamma augmentation.
  arXiv:1506.05843
- Li N, Lee R (2005). Coherent mortality forecasts for a group of
  populations. Demography 42(3): 575-594. doi:10.1353/dem.2005.0021
- Scheuerer M, Hamill TM (2015). Variogram-based proper scoring rules for
  probabilistic forecasts of multivariate quantities. Monthly Weather
  Review 143(4): 1321-1334. doi:10.1175/MWR-D-14-00269.1
