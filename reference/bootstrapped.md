# Bootstrap uncertainty of preprocessing and parameters

Wraps a model so that its forecast draws include the uncertainty of the
estimated parameters and, optionally, of the preprocessing (smoothing,
ungrouping). The model is refitted on `n_boot` resampled datasets and
the forecast draws are pooled across refits:

## Usage

``` r
bootstrapped(spec, n_boot = 50, type = c("grid", "data"), seed = 1)
```

## Arguments

- spec:

  A model specification.

- n_boot:

  Number of bootstrap refits.

- type:

  `"grid"` or `"data"`.

- seed:

  Random seed for the resampling.

## Value

A model specification.

## Details

- `type = "grid"`: deaths on the grid (after preprocessing) are drawn
  from Poisson distributions with the fitted rates of the model (or the
  preprocessed rates when the model has no fitted rates); the model is
  refitted without preprocessing. This is the semi-parametric bootstrap
  of Brouhns et al. (2005), generalised to any model.

- `type = "data"`: the observed deaths, at their original resolution,
  are drawn from Poisson distributions with the observed counts as
  means, and preprocessing and model are refitted. This propagates the
  uncertainty of ungrouping and smoothing, for any smoother.

## References

Brouhns N, Denuit M, Van Keilegom I (2005). Bootstrapping the Poisson
log-bilinear model for mortality forecasting. *Scandinavian Actuarial
Journal* 2005(3): 212-224.

## Examples

``` r
x <- mt_example(deaths = "single", years = 2000:2019)
fit <- mt_fit(bootstrapped(rwd(), n_boot = 5), x)
fc <- mt_forecast(fit, h = 3, n_draws = 50, seed = 1)
```
