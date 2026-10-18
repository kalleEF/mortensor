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
  max_refits = 1
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

## Examples

``` r
# \donttest{
if (requireNamespace("bamp", quietly = TRUE)) {
  x <- mt_example(years = 2005:2019)
  fit <- mt_fit(apc_bayes(mcmc = mcmc_control("test")), x, prep = mt_prep(step_ungroup()))
  fc <- mt_forecast(fit, h = 5, n_draws = 200, seed = 1)
}
# }
```
