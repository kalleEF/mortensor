# Fit a model to mortality data

Applies the preprocessing pipeline, takes deaths and exposures on the
grid and fits the model per stratum (or jointly). A failure in one
stratum is recorded in the diagnostics and does not stop the other
strata.

## Usage

``` r
mt_fit(spec, data, prep = NULL, cores = 1, ...)
```

## Arguments

- spec:

  A model specification (`mt_spec`).

- data:

  An [mt_data](mt_data.md) object.

- prep:

  An [mt_prep](mt_prep.md) pipeline or `NULL`.

- cores:

  Number of cores for per-stratum fits (forked processes; 1 on Windows).
  Each stratum is fitted with its own random number stream, started from
  the current generator, so fits that use random numbers (such as
  `lee_carter(bootstrap = TRUE)`) are reproduced by
  [`set.seed()`](https://rdrr.io/r/base/Random.html) with any number of
  cores.

- ...:

  Passed to methods.

## Value

An object of class `mt_fit`.

## Details

Smooth inside the fit (`prep = mt_prep(step_smooth(...))`), not before
it: data returned by [`mt_smooth()`](mt_smooth.md) carry no raw deaths,
so residuals and observation noise would refer to the smoothed rates and
the intervals would be too narrow. Such data give a warning.

## Examples

``` r
x <- mt_example(deaths = "single")
fit <- mt_fit(lee_carter(), x)
fit
#> <mt_fit> Lee-Carter; 2 strata; ages 30 to 100+; periods 1990-2019
```
