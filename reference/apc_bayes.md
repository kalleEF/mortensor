# Bayesian age-period-cohort model (bamp)

Wrapper around the bamp package (Schmid and Held 2007): binomial
likelihood for deaths given exposures, logit link, random walk priors of
first or second order for age, period and cohort effects, optional
overdispersion, MCMC with four chains. Forecasts extrapolate period and
cohort effects as random walks for every posterior draw
([`bamp::predict_apc()`](https://volkerschmid.github.io/bamp/reference/predict_apc.html)),
so the draws include parameter uncertainty.

## Usage

``` r
apc_bayes(
  age = "rw2",
  period = "rw1",
  cohort = "rw1",
  overdispersion = FALSE,
  mcmc = mcmc_control(),
  resolution = c("observed", "grid"),
  shape_years = 5,
  max_refits = 1,
  coherent = FALSE,
  common_by = NULL,
  reference = NULL,
  cohort_deviation = FALSE,
  seed = NULL
)
```

## Arguments

- age, period, cohort:

  Prior for each effect: `"rw1"` or `"rw2"`. The defaults use
  first-order random walks for period and cohort: with second-order
  random walks the extrapolated trends can make forecast intervals span
  orders of magnitude. First-order random walks carry the last period
  effect forward without a trend, so forecasts are biased when mortality
  keeps declining. The vignette "Bayesian age-period-cohort models with
  bamp" compares both.

- overdispersion:

  Include cell-level overdispersion.

- mcmc:

  MCMC settings from [`mcmc_control()`](mcmc_control.md).

- resolution:

  `"observed"` or `"grid"`.

- shape_years:

  Number of recent years used for the shape factors.

- max_refits:

  Refits with longer chains if not converged.

- coherent:

  Fit the strata of each group jointly with shared period and cohort
  effects and mean-reverting stratum deviations (see the section
  "Coherent model"); needs at least two strata per group.

- common_by:

  Stratum variable(s) defining groups with their own common trend
  (`coherent = TRUE`), for example `"sex"` for the deprivation groups of
  each sex. `NULL`: one group of all strata.

- reference:

  With `coherent = TRUE`: `NULL` (the deviations sum to zero over the
  strata of a group) or a named list of values of stratum variables that
  selects one reference stratum in every group, for example
  `list(quintile = "England")`, whose deviation is zero.

- cohort_deviation:

  With `coherent = TRUE`: add AR(1) deviations of the strata along the
  cohorts.

- seed:

  Seed of the native sampler (`coherent = TRUE`); `NULL` draws one from
  R's random number generator, so that
  [`set.seed()`](https://rdrr.io/r/base/Random.html) makes fits
  reproducible.

## Value

A model specification.

## Details

The wrapper makes bamp behave like every other model:

- **Resolution.** With `resolution = "observed"` (default) the model is
  fitted where the deaths were observed: 5-year age groups with annual
  periods use `periods_per_agegroup = 5` (Knorr-Held and Rainer 2001).
  Exposures are aggregated to the death cells. Age groups must have
  equal widths; an open last group is treated as one more group of that
  width.

- **Grid.** Draws for age groups are mapped to the modelling grid with
  age-specific shape factors (single-year rate relative to its group
  rate, exposure-weighted to 1 within each group) estimated from the
  preprocessed training data of the last `shape_years` years. Use
  [`step_ungroup()`](mt_prep.md) in `prep` when deaths are grouped.

- **Convergence.**
  [`bamp::checkConvergence()`](https://volkerschmid.github.io/bamp/reference/checkConvergence.html)
  is called explicitly; a stratum that has not converged is refitted
  with twice as many iterations up to `max_refits` times, then reported
  as failed in the diagnostics.

- **Integers.** bamp needs whole numbers; deaths and exposures are
  rounded and the largest relative rounding change is recorded.

Heterogeneity variants (`"rw1+het"`, `"rw2+het"`) are not offered
because
[`bamp::predict_apc()`](https://volkerschmid.github.io/bamp/reference/predict_apc.html)
does not forecast them in bamp 3.0.0.

## Coherent model (`coherent = TRUE`)

Fits the strata of each group of `common_by` (all strata without it) in
one model, so that their forecasts do not diverge:

`logit p[s, i, j] = mu + a[s] + theta[s, i] + phi[j] + psi[k] + d[s, j]`

for stratum `s`, age group `i`, year `j` and cohort `k`, with binomial
deaths given the exposures. The period effect `phi` and the cohort
effect `psi` are shared and carry the common trend; their random-walk
priors (orders `period`, `cohort`) are scaled to unit generalised
variance (Sørbye and Rue 2014) and forecast as random walks. Each
stratum has its own age effect `theta[s, ]` (random walk of order `age`,
one precision for all strata) and level `a[s]`. The deviations `d[s, ]`
of the strata from the shared period effect are stationary AR(1) series
with a common coefficient `rho` (uniform prior on 0 to 1, sampled) and
innovation precision; they sum to zero over the strata in every year,
or, with `reference`, the deviation of the reference stratum is zero.
Their forecasts revert towards zero, so the gaps between strata stay
bounded. `cohort_deviation = TRUE` adds deviations `c[s, k]` of the same
form along the cohorts. The precisions have Gamma priors with shape 1
and rates 0.5 (age) and 5e-4 (period, cohort, deviations). For the
deviations the branch uses the rate 0.05, which dominates typical data:
the posterior of their precision and the width of the forecast gaps are
then set by the prior.

**Sources.** The model is the coherent model of the development branch
"coherent_forecasting" of bamp (`bamp_coherent()`, by Chris Kypridemos),
rewritten in this package because that branch is not released; when a
bamp release contains it, `apc_bayes(coherent = TRUE)` will call bamp
and the native sampler will be deprecated. Its published anchor is the
two-population Bayesian APC model of Cairns et al. (2011): a common
period trend and a mean-reverting AR(1) spread with level zero,
population-specific age effects, a single-stage Bayesian fit, and two
parameterisations, the spread between equal populations (their Section
4.5, the default here) and the spread from a dominant population
(Section 4.4, `reference`). It differs from that paper, and these parts
are package-specific (reference unknown): the binomial-logit likelihood
and the random-walk smoothing priors of bamp (Schmid and Held 2007) in
place of the Poisson likelihood with a random walk with drift and AR(2)
cohort processes; one cohort effect shared by all strata (the
alternative cohort hypothesis of their Section 4.3), with optional AR(1)
cohort deviations in place of an AR(2) cohort spread; independent
innovations of the common trend and the deviations (the paper:
correlated); the priors; and more than two strata.

**Sampler.** Four chains of a Metropolis-Hastings sampler: all effects
are proposed jointly from a Gaussian approximation of their conditional
posterior (one Newton step from the current state, with sparse Cholesky
factors; Rue and Held 2005, Section 4.4), the precisions are drawn from
their Gamma conditionals, `rho` by a random-walk Metropolis step on the
logit scale whose step size adapts during the burn-in only. Convergence
is checked with the Gelman-Rubin statistic of the hyperparameters, the
last period effect and the level of every stratum in the last year
(`rhat` below 1.1); a group that has not converged is refitted with
twice as many iterations up to `max_refits` times. Overdispersion is not
available in the coherent model.

Draws are posterior draws of the rates. Predictive draws
(`mt_forecast(predictive = TRUE)`) add Poisson variation of the deaths
given the rates of each draw and the exposures, as for the other models
of expected rates; no other observation noise is added. bamp's binomial
likelihood and its overdispersion term (`overdispersion = TRUE`) are not
carried into the predictive draws: the overdispersion only enters the
posterior of the rates. If intervals are too narrow on real data, that
is a property of the model (see [`recalibrated()`](recalibrated.md),
[`ensemble()`](ensemble.md), or overdispersed models such as
`stanmomo(family = "nb")`).

## References

Knorr-Held L, Rainer E (2001). Prognosis of lung cancer mortality in
West Germany: a case study in Bayesian prediction. *Biostatistics* 2(1):
109-129.

Schmid VJ, Held L (2007). Bayesian age-period-cohort modeling and
prediction: BAMP. *Journal of Statistical Software* 21(8): 1-15.
[doi:10.18637/jss.v021.i08](https://doi.org/10.18637/jss.v021.i08)

Cairns AJG, Blake D, Dowd K, Coughlan GD, Khalaf-Allah M (2011).
Bayesian stochastic mortality modelling for two populations. *ASTIN
Bulletin* 41(1): 29-59.
[doi:10.2143/AST.41.1.2084385](https://doi.org/10.2143/AST.41.1.2084385)

Gelman A, Rubin DB (1992). Inference from iterative simulation using
multiple sequences. *Statistical Science* 7(4): 457-472.
[doi:10.1214/ss/1177011136](https://doi.org/10.1214/ss/1177011136)

Rue H, Held L (2005). *Gaussian Markov Random Fields: Theory and
Applications*. Chapman and Hall/CRC, Boca Raton.

Sørbye SH, Rue H (2014). Scaling intrinsic Gaussian Markov random field
priors in spatial modelling. *Spatial Statistics* 8: 39-51.
[doi:10.1016/j.spasta.2013.06.004](https://doi.org/10.1016/j.spasta.2013.06.004)

## Examples

``` r
# \donttest{
if (requireNamespace("bamp", quietly = TRUE)) {
  x <- mt_example(years = 2005:2019)
  fit <- mt_fit(apc_bayes(mcmc = mcmc_control("test")), x, prep = mt_prep(step_ungroup()))
  fc <- mt_forecast(fit, h = 5, n_draws = 200, seed = 1)
}
# }
# coherent model of the two sexes (native sampler, no bamp needed)
# \donttest{
x <- mt_example(years = 2005:2019)
fit <- mt_fit(apc_bayes(coherent = TRUE, mcmc = mcmc_control("test"), seed = 1), x,
              prep = mt_prep(step_ungroup()))
# }
```
