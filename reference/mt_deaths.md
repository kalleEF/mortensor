# Expected deaths from draws and future exposures

Multiplies rate draws by exposures and optionally adds Poisson noise, so
that the result describes observable death counts.

## Usage

``` r
mt_deaths(draws, exposure = NULL, noise = TRUE, seed = NULL)
```

## Arguments

- draws:

  An [mt_draws](mt_draws.md) object.

- exposure:

  Future exposures: an [mt_data](mt_data.md) object with exposures (for
  example a population projection converted with
  [`mt_as_exposure()`](mt_as_exposure.md)), or an array
  `[age, period, stratum]` or matrix `[age, stratum]` (used for every
  period) with dimnames, as [`mt_array()`](mt_array.md) returns. Ages,
  years and strata are matched by their labels. Default: the exposures
  of the last training year.

- noise:

  Add Poisson noise (draws of expected rates only; ignored for
  predictive draws, which contain it).

- seed:

  Random seed.

## Value

An array `[age, period, stratum, draw]` of deaths.

## Details

Predictive draws (`mt_forecast(predictive = TRUE)`) already contain the
variation of the observed deaths (Poisson given the rates, or negative
binomial for `stanmomo(family = "nb")`): for them no noise is added,
whatever `noise`, and the result is their deaths. Give the exposures
that were used for the forecast (its `exposure`); other exposures
rescale the deaths and the result is then not a count.

## Examples

``` r
x <- mt_example(deaths = "single")
fc <- mt_forecast(mt_fit(rwd(), x), h = 2, n_draws = 20, seed = 1)
dim(mt_deaths(fc))
#> [1] 71  2  2 20
```
