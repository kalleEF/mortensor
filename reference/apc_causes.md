# Multi-cause age-period-cohort model with coupled cause trends

Forecasts the deaths of several causes together with all causes, so that
the cause-specific rates add up to the all-cause rate in every draw and
the trends of the causes move together as estimated. Package-specific
(reference unknown): no published source of this combination is known to
us.

## Usage

``` r
apc_causes(
  var = "cause",
  total_level = "all",
  residual = "auto",
  age = "rw2",
  period = "rw1",
  cohort = "rw1",
  mcmc = mcmc_control(),
  resolution = c("observed", "grid"),
  shape_years = 5,
  max_refits = 1,
  seed = NULL
)
```

## Arguments

- var:

  Stratum variable with the causes.

- total_level:

  Level of `var` for all causes.

- residual:

  Add a residual cause, all deaths not covered by the listed causes:
  `TRUE`, `FALSE`, or `"auto"` (added when the listed causes leave more
  than 0.1% of all-cause deaths unexplained in the median cell). The
  residual cause is modelled but not returned.

- age, period, cohort:

  Prior for each effect: `"rw1"` or `"rw2"`. The defaults use
  first-order random walks for period and cohort: with second-order
  random walks the extrapolated trends can make forecast intervals span
  orders of magnitude. First-order random walks carry the last period
  effect forward without a trend, so forecasts are biased when mortality
  keeps declining. The vignette "Bayesian age-period-cohort models with
  bamp" compares both.

- mcmc:

  MCMC settings from [`mcmc_control()`](mcmc_control.md).

- resolution:

  `"observed"` or `"grid"`.

- shape_years:

  Number of recent years used for the shape factors.

- max_refits:

  Refits with longer chains if not converged.

- seed:

  Seed of the native sampler (`coherent = TRUE`); `NULL` draws one from
  R's random number generator, so that
  [`set.seed()`](https://rdrr.io/r/base/Random.html) makes fits
  reproducible.

## Value

A model specification.

## Details

**Model.** In each group of strata that differ only in `var` (for
example the causes of one sex), the all-cause deaths follow a
binomial-logit age-period-cohort model with random-walk priors (as
[`apc_bayes()`](apc_bayes.md)). The cause deaths split the all-cause
deaths of each cell by stick-breaking: a sequence of binomial-logit
models of the share of each cause among the deaths not yet assigned,
`logit pi[c, i, j] = mu[c] + theta[c, i] + phi[c, j] + psi[c, k]`, the
causes ordered by their number of deaths and the last cause taking the
rest (Linderman, Johnson and Adams 2015 for the stick-breaking
representation of a multinomial). The period effects of the shares have
a joint random-walk prior whose precision across causes is a matrix with
a Wishart prior, so that period trends of different causes can be
correlated, also negatively (one cause replacing another); the cohort
effects likewise. Cause rates are shares times the all-cause rate: they
add up, and hazards built from them are additive. Within an age group,
every cause takes the age pattern of the all-cause rate on the modelling
grid.

**Sources.** The model is the multi-cause model of the development
branch "coherent_forecasting" of bamp (`bamp_multicause()`, by Chris
Kypridemos), rewritten here because the branch is not released; when a
bamp release contains it, `apc_causes()` will call bamp and the native
sampler will be deprecated. Changes: the exact sampler of
`apc_bayes(coherent = TRUE)`; one zero-slope constraint on the period
effect of each share, only when all effects are second-order random
walks (the branch also pins the cohort slope, which removes an
identified trend); the all-cause model fitted with the same sampler, its
period and cohort precisions with Gamma priors (Wishart priors of
dimension 1); the latent-factor coupling of the branch is not taken
over. Priors: Gamma(1, 0.5) for the precisions of the age effects,
Wishart with `C + 2` degrees of freedom and scale `(1e-4 I)^-1` for the
period and cohort precisions of `C` fields, as in the branch.

**Sampler and convergence** as [`apc_bayes()`](apc_bayes.md) with
`coherent = TRUE`: four chains, the Gelman-Rubin statistic of the
precisions of the age effects, the diagonals of the Wishart precisions,
the last period effect and the level of every field in the last year.

## References

Linderman SW, Johnson MJ, Adams RP (2015). Dependent multinomial models
made easy: stick-breaking with the Pólya-gamma augmentation.
arXiv:1506.05843.

Schmid VJ, Held L (2007). Bayesian age-period-cohort modeling and
prediction: BAMP. *Journal of Statistical Software* 21(8): 1-15.
[doi:10.18637/jss.v021.i08](https://doi.org/10.18637/jss.v021.i08)

## See also

[`coda_causes()`](coda_causes.md) for a published model of causes of
death, and [`reconciled()`](reconciled.md) with a hierarchy of causes.

## Examples

``` r
# \donttest{
x <- mt_example_causes(years = 2005:2019)
fit <- mt_fit(apc_causes(mcmc = mcmc_control("test"), seed = 1), x)
fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1)
# }
```
