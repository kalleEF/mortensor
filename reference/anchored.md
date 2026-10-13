# Control the jump-off of any model

Forecasts of most models start from their fitted rates in the last
observed year. Where the fit misses the level of the last years (for
example after a recent change of trend), the forecast starts with a step
away from the observations; where the observations are noisy, starting
from them carries the noise into every forecast year. `anchored()` wraps
any model specification and moves the start of its forecasts towards the
observed level, per age and stratum:

## Usage

``` r
anchored(
  spec,
  to = c("observed", "smoothed"),
  years = 3,
  decay = 10,
  shrink = TRUE
)
```

## Arguments

- spec:

  A model specification.

- to:

  `"observed"` or `"smoothed"`.

- years:

  Number of last training years whose deaths define the offset.

- decay:

  Years over which the offset fades out (`Inf`: constant).

- shrink:

  Shrink offsets of cells with few deaths.

## Value

A model specification. The fitted object has the offsets in
`$anchor$offset` (`[age, stratum]`, log scale).

## Details

1.  the offset is the log ratio of observed to fitted deaths, summed
    over the last `years` training years (`to = "observed"`: the rates
    before smoothing; `to = "smoothed"`: the rates after preprocessing,
    which differs when the pipeline smooths);

2.  offsets of cells with few deaths are shrunk towards the
    deaths-weighted mean offset of their stratum (empirical Bayes, with
    the Poisson variance of each offset; `shrink = TRUE`);

3.  the forecast draws are multiplied by `exp(offset * w(h))`, where the
    weight `w(h)` falls linearly from 1 at the first forecast year to 0
    after `decay` years (`decay = Inf` keeps the offset, `decay = 0`
    leaves the model's own jump-off).

`years = 1, decay = Inf, shrink = FALSE` is the observed jump-off of
Lee-Carter models (`lee_carter(jump_off = "observed")`). Models that
have no fitted rates (random walks, time series per age) start from the
observed rates already and are left unchanged; in an
[`ensemble()`](ensemble.md) each draw is anchored with the offset of the
component it comes from. For [`reconciled()`](reconciled.md) models the
bottom strata are anchored and the aggregates recomputed, so the draws
stay coherent. [`mt_jump_off()`](mt_jump_off.md) reports the step
between the last observed or fitted year and the first forecast year.

In backtests on England & Wales men (origins 1996-2001, h = 10),
anchoring Lee-Carter to the last 3 years (`decay = 10`) lowered the CRPS
of the first two forecast years from 0.040 to 0.033 and over all
horizons from 0.0687 to 0.0648; the observed jump-off of the last year
alone gave 0.0583. On the ONS deprivation data (h = 5), where the last
observed rates of the quintiles are noisier, anchoring Li-Lee to 3 years
improved the CRPS slightly (0.0658 against 0.0666), while the observed
jump-off of one year made it worse (0.0708).

## Examples

``` r
x <- mt_example(deaths = "single", years = 1995:2019)
fit <- mt_fit(anchored(lee_carter(), years = 3, decay = 10), x)
summary(as.vector(fit$anchor$offset))
#>       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#> -0.0002439 -0.0002439  0.0003177  0.0003177  0.0008793  0.0008793 
fc <- mt_forecast(fit, h = 10, n_draws = 100, seed = 1)
mt_jump_off(fc)
#>       sex  label observed   fitted forecast step_observed step_fitted
#>    <char> <char>    <num>    <num>    <num>         <num>       <num>
#> 1: female female 55.95597 55.95828 56.05399    0.09802129  0.09571211
#> 2:   male   male 51.16230 51.17137 51.26540    0.10310001  0.09403114
```
