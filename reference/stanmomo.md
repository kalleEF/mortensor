# Bayesian stochastic mortality models in Stan (StanMoMo)

Wrapper around StanMoMo (Barigou, Goffard, Loisel and Salhi 2023):
Bayesian versions of the Lee-Carter (`"lc"`), Renshaw-Haberman (`"rh"`),
age-period-cohort (`"apc"`), Cairns-Blake-Dowd (`"cbd"`) and M6 (`"m6"`)
models, fitted with Stan to the deaths and exposures of each stratum
with a Poisson or negative binomial likelihood (log link). Period
indexes follow random walks with drift (CBD and M6: correlated), cohort
effects AR(2) processes.

## Usage

``` r
stanmomo(
  model = c("lc", "rh", "apc", "cbd", "m6"),
  family = c("poisson", "nb"),
  h_max = 50,
  chains = 4,
  iter = 2000,
  cores = 1,
  seed = NULL
)
```

## Arguments

- model:

  `"lc"`, `"rh"`, `"apc"`, `"cbd"` or `"m6"`.

- family:

  `"poisson"` or `"nb"` (negative binomial).

- h_max:

  Longest forecast horizon the fit supports (at least 2).

- chains, iter:

  Number of Markov chains and iterations per chain (half of them
  warm-up).

- cores:

  Cores for the chains.

- seed:

  Seed of the sampler; `NULL` draws one from R's random number
  generator, so that [`set.seed()`](https://rdrr.io/r/base/Random.html)
  makes fits reproducible.

## Value

A model specification.

## Details

StanMoMo forecasts inside the fit, so the model is fitted with forecasts
up to `h_max` years. Draws are posterior draws of the expected rates,
rebuilt from the posterior parameters and forecast indexes; they carry
the uncertainty of the parameters and of the future indexes. Predictive
draws add Poisson and extra-Poisson variation of the observed deaths as
for other models (for the negative binomial family this takes the place
of its overdispersion). Fitted rates are posterior means.

The fit warns when the chains have not converged (largest R-hat of the
model parameters above 1.05) or had divergent transitions; then use more
`iter` or `chains`. Cohort models need single years of age. StanMoMo
needs whole numbers: deaths and exposures are rounded (a warning reports
rounding of exposures by more than 1%). Fits take from seconds to
minutes per stratum; `cores` runs the chains in parallel.

## References

Barigou K, Goffard P-O, Loisel S, Salhi Y (2023). Bayesian model
averaging for mortality forecasting using leave-future-out validation.
*International Journal of Forecasting* 39(2): 674-690.
[doi:10.1016/j.ijforecast.2022.01.011](https://doi.org/10.1016/j.ijforecast.2022.01.011)

## See also

[`stmomo()`](stmomo.md) for the same models fitted by maximum
likelihood.

## Examples

``` r
# \donttest{
if (requireNamespace("StanMoMo", quietly = TRUE)) {
  x <- mt_example(deaths = "single", years = 2000:2019)
  fit <- mt_fit(stanmomo("lc", chains = 1, iter = 400, h_max = 10), x)
  fc <- mt_forecast(fit, h = 10, n_draws = 100, seed = 1)
}
#> Warning: StanMoMo Lee-Carter: warnings in 2 stratum/strata, for example female: StanMoMo chains have not converged (largest R-hat 1.07); use more `iter` or `chains`. (see the diagnostics of the fit).
# }
```
