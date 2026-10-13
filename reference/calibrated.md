# Calibrate forecast distributions on recent years

Wraps any model specification so that the spread (and optionally the
trend) of its forecast draws is recalibrated on the last `years`
observed years of the training data. Each time the model is fitted (also
at every origin of an outer [`mt_backtest()`](mt_backtest.md)):

## Usage

``` r
calibrated(
  spec,
  years = 10,
  n_origins = 5,
  scale = c("horizon", "constant"),
  bias = FALSE,
  target = c("crps", "interval"),
  level = 0.95,
  by = NULL,
  n_draws = 100,
  seed = 1
)
```

## Arguments

- spec:

  A model specification.

- years:

  Length of the calibration window (the last observed years).

- n_origins:

  Number of origins in the window (evenly spaced from `years` to 2 years
  before the end of the training data).

- scale:

  `"horizon"` (spread factor \\s_1 h^{\gamma}\\) or `"constant"`.

- bias:

  Estimate a correction of the drift (log rate per year).

- target:

  `"crps"` or `"interval"` (interval score at `level`).

- level:

  Level of the central interval for `target = "interval"`.

- by:

  Stratum variables with separate calibration parameters.

- n_draws:

  Draws per inner forecast.

- seed:

  Random seed of the inner forecasts and of the simulated observation
  noise.

## Value

A model specification. The fitted object has a `calibration` element: a
`data.table` with one row per group (`scale_h1`, `slope`, `drift`, the
number of scored cells, and CRPS and coverage of the inner forecasts
before and after calibration).

## Details

1.  the base model is refitted at `n_origins` origins inside the
    calibration window (training data up to the origin) and forecasts
    the remaining years of the window;

2.  the draws are transformed per cell on the log scale, \$\$\log m^\* =
    \tilde m + b\\h + s(h)\\(\log m - \tilde m),\$\$ where \\\tilde m\\
    is the median of the draws in the cell, \\h\\ the horizon, \\s(h) =
    s_1 h^{\gamma}\\ a spread factor (`scale = "horizon"`; \\\gamma =
    0\\ for `scale = "constant"`), constant beyond the longest horizon
    of the window, and \\b\\ a correction of the drift per year
    (`bias = TRUE`);

3.  \\s_1\\, \\\gamma\\ and \\b\\ minimise the mean CRPS (or the
    interval score, `target`) of the log death rates in the window, with
    Poisson (and the model's extra-Poisson) variation of the observed
    deaths, as in [`mt_backtest()`](mt_backtest.md). Draws of an
    ensemble get the observation noise of the component they come from.

The final model is fitted on all training data and its forecasts are
transformed with the estimated parameters. This is minimum-CRPS
estimation of a location-scale recalibration (Gneiting et al. 2005)
applied to simulated paths: ranks within each cell are kept, so paths
stay paths, and derived quantities such as life expectancy can be
computed from the calibrated draws.

The wrapper works with any specification, including
[`ensemble()`](ensemble.md), [`bootstrapped()`](bootstrapped.md) and
joint models. For [`reconciled()`](reconciled.md) models the bottom
strata are calibrated and the aggregates recomputed from them, so the
draws stay coherent. With `by`, separate parameters are estimated for
groups of strata (for example `by = "sex"`); strata outside these groups
(aggregates) use the pooled parameters. Groups with different factors
scale the common variation of joint models (such as the common trend of
[`li_lee()`](li_lee.md)) differently, so it no longer cancels in
contrasts between the groups, and intervals of gaps between them widen;
use pooled parameters when contrasts between the groups are the target.

Calibration needs `years` more training years than the base model. A
spread factor above 1 means that the base model was too narrow in the
window; it widens future intervals in proportion. The parameters are
estimated from few overlapping forecasts (with the defaults, 5 origins
and 30 forecast years, horizons 1 to 10), so they carry their own
uncertainty, and they cannot anticipate changes after the window.
Compare calibrated and uncalibrated models in an outer backtest.

In backtests on England & Wales men (HMD via StMoMo, origins 1996-2001,
h = 10), calibration raised the 95% coverage of
[`lee_carter()`](lee_carter.md) from 0.75 to 0.88 and lowered its CRPS
by 4% (0.0687 to 0.0657); `bias = TRUE` gave the same. On the ONS
deprivation data (6-year window), the coverage of Lee-Carter and Li-Lee
rose with `target = "interval"` (0.83 to 0.90 and 0.84 to 0.89) at
unchanged CRPS, and `bias = TRUE` made the CRPS worse (0.069 to 0.082).
An equal-weight [`ensemble()`](ensemble.md) of Lee-Carter and
[`rwd()`](rwd.md) was already calibrated; calibration then narrowed it
(95% coverage 0.93 to 0.85 on ONS, 0.95 to 0.90 on England & Wales) at
the same CRPS (see the vignette "Forecasting mortality by deprivation
quintile").

`target = "crps"` (default) balances calibration and sharpness;
`target = "interval"` minimises the interval score of the central
`level` interval (Gneiting and Raftery 2007), whose optimum is the
correct coverage of that interval.

## References

Gneiting T, Raftery AE, Westveld AH, Goldman T (2005). Calibrated
probabilistic forecasting using ensemble model output statistics and
minimum CRPS estimation. *Monthly Weather Review* 133(5): 1098-1118.
[doi:10.1175/MWR2904.1](https://doi.org/10.1175/MWR2904.1)

Gneiting T, Raftery AE (2007). Strictly proper scoring rules,
prediction, and estimation. *Journal of the American Statistical
Association* 102(477): 359-378.
[doi:10.1198/016214506000001437](https://doi.org/10.1198/016214506000001437)

## Examples

``` r
x <- mt_example(deaths = "single", years = 1990:2019)
fit <- mt_fit(calibrated(lee_carter(), years = 8, n_origins = 3, n_draws = 50), x)
fit$calibration
#>     group scale_h1      slope drift cells crps_before crps_after cov80_before
#>    <char>    <num>      <num> <num> <int>       <num>      <num>        <num>
#> 1:    all 1.244915 -0.9998819     0  2414  0.03398005 0.03389239    0.7692626
#>    cov80_after cov95_before cov95_after
#>          <num>        <num>       <num>
#> 1:   0.7228666    0.9092792   0.8719967
fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1)
```
