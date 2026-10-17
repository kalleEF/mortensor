# Backtesting and scores

How good is a forecast? The future is unknown, but the past can be used
as a test: fit a model to the years up to some year in the past (the
*origin*), forecast the following years, and compare the forecast with
what was observed. Repeating this from several origins gives a
*rolling-origin backtest* (Tashman 2000). This vignette explains the
design and the scores, and compares three models on one data set.

## Data and models

Deaths and populations in England by sex and deprivation quintile,
2001-2019 (ONS). Three models, chosen as examples:

- a **random walk with drift** per age ([`rwd()`](../reference/rwd.md)):
  each log rate continues its average change of the past. It is a
  standard benchmark in forecasting: a model that does not beat it adds
  little;
- **Lee-Carter** ([`lee_carter()`](../reference/lee_carter.md)), one
  model per stratum;
- **Li-Lee** (`li_lee(common_by = "sex")`), quintiles of each sex pooled
  around a common trend.

``` r

library(mortensor)
library(data.table)
x <- mt_data(deaths = ons_data, population = ons_data, by = c("sex", "imd_decile"), age = "age_group",
             value = c(deaths = "all_causes", population = "population"))
x <- mt_as_exposure(x, type = "midyear")
x <- mt_group_levels(x, "imd_decile", function(d) paste0("Q", ceiling(as.numeric(d) / 2)), name = "quintile")
models <- list("random walk" = rwd(), "Lee-Carter" = lee_carter(), "Li-Lee" = li_lee(common_by = "sex"))
```

## The design

Seven origins, 2008 to 2014, each forecasting five years:

``` r

mt_plot_backtest_design(x, origins = 2008:2014, h = 5)
```

![](backtesting_files/figure-html/unnamed-chunk-2-1.png)

Every model, and every preprocessing step such as smoothing, is fitted
anew at each origin from the training years only, so no information from
the test years leaks into the forecast. Each test year is forecast
several times from different origins and horizons, which averages over
the particular events of a single period.

``` r

bt <- mt_compare(models, x, origins = 2008:2014, h = 5, n_draws = 200, quantities = "life_expectancy")
```

The forecasts are *predictive*: they describe the observed death rates
of the test years, including the Poisson noise of the counts, because
the observed rates are what they are compared with.

## Scores

### Interval coverage

A 95% interval should contain the observed value in about 95% of the
cases. The *empirical coverage* is the share of observed values inside
the interval: below the nominal level, the intervals are too narrow;
above it, too wide.

``` r

mt_plot_coverage(bt)
```

![](backtesting_files/figure-html/unnamed-chunk-4-1.png)

Coverage alone can be gamed: a very wide interval always covers. It is
read together with a score that also rewards narrow intervals.

### CRPS

The continuous ranked probability score (CRPS) compares the whole
forecast distribution with the observed value: it is the integral of the
squared difference between the forecast’s cumulative distribution and a
step at the observed value (Gneiting and Raftery 2007). The figure
shades the difference between the two curves; squared and summed over
the axis, it gives the CRPS:

![](backtesting_files/figure-html/unnamed-chunk-5-1.png)

A forecast that puts its probability close to the observed value has a
small area. The CRPS punishes both a distribution that is off target and
one that is wider than necessary. For a forecast without uncertainty it
equals the absolute error, so it is in the unit of the forecast
quantity: for log death rates, a CRPS of 0.05 means errors of roughly 5%
of the rate; for life expectancy it is in years. Lower is better.

The *weighted interval score* (WIS, Bracher et al. 2021) approximates
the CRPS from a few central intervals; it is reported as well and ranks
models in the same way.

### Point errors

`mae` and `rmse` are the mean absolute and the root mean squared error
of the forecast median: they ignore the intervals. `mase` divides the
absolute error by the mean absolute year-to-year change of the training
data (Hyndman and Koehler 2006): below 1, the forecast errs less than
the typical change of a single year.

### PIT

The probability integral transform (PIT) is the share of forecast draws
below the observed value. For calibrated forecasts it is spread evenly
between 0 and 1, so its histogram is flat. A U shape means that the
observed values often fall in the tails, outside intervals that are too
narrow; a slope means that the forecasts are systematically too high or
too low (Gneiting, Balabdaoui and Raftery 2007).

``` r

mt_plot_pit(bt)
```

![](backtesting_files/figure-html/unnamed-chunk-6-1.png)

All three histograms rise to the right: the observed rates were more
often above the forecast median than below it. The forecasts were too
optimistic, because the decline of mortality in England slowed after
2011. Lee-Carter and Li-Lee also pile up at both ends (a U shape): their
intervals are too narrow.

## Comparing the models

``` r

mt_score_table(bt, metrics = c("crps", "rmse", "cov80", "cov95"))
#>          model   crps  rmse cov80 cov95 crps_le rmse_le cov80_le cov95_le
#>         <fctr>  <num> <num> <num> <num>   <num>   <num>    <num>    <num>
#> 1: random walk 0.0850 0.180 0.803 0.915   0.293   0.503    0.603    0.820
#> 2:  Lee-Carter 0.0678 0.140 0.654 0.829   0.292   0.507    0.626    0.823
#> 3:      Li-Lee 0.0666 0.132 0.673 0.847   0.336   0.601    0.614    0.811
```

The columns without suffix score the log death rates of all cells; the
columns with `_le` score life expectancy at birth (in years), computed
from the same forecasts. Rates and life expectancy can rank models
differently, because life expectancy depends on how the errors of the
ages add up: for the rates Li-Lee is the most accurate (CRPS 0.067) and
the random walk the least (0.085), for life expectancy the random walk
and Lee-Carter are equal (0.29 years) and Li-Lee is the least accurate
(0.34 years). Choose the target that the analysis is about.

``` r

mt_plot_scores(bt)
```

![](backtesting_files/figure-html/unnamed-chunk-8-1.png)

All three models get worse with the horizon, as they should. The random
walk has wide intervals: its 95% coverage is close to nominal (0.92),
and its 50% and 80% intervals cover about as often as they should. The
95% intervals of Lee-Carter and Li-Lee cover 0.83 and 0.85 of the
observed rates, fewer at longer horizons: they describe the uncertainty
of a fixed trend, while the trend of English mortality changed after
2011. “Tuning and recalibration” shows how to widen intervals to the
errors of recent forecasts, and “Choosing and combining models” how
mixtures of models do.

## How far to trust a backtest

- Seven origins of one country are few, and the test years overlap; a
  difference of a few percent in the CRPS can be noise.
- A backtest scores the past. If the future differs from the test period
  (a pandemic, a policy change), so will the accuracy.
- Compare models on the same origins, horizons and cells, as
  [`mt_compare()`](../reference/mt_compare.md) does. Suppressed cells
  that were estimated are not scored.

## References

Bracher J, Ray EL, Gneiting T, Reich NG (2021). Evaluating epidemic
forecasts in an interval format. *PLoS Computational Biology* 17(2):
e1008618.

Gneiting T, Balabdaoui F, Raftery AE (2007). Probabilistic forecasts,
calibration and sharpness. *Journal of the Royal Statistical Society
Series B* 69(2): 243-268.

Gneiting T, Raftery AE (2007). Strictly proper scoring rules,
prediction, and estimation. *Journal of the American Statistical
Association* 102(477): 359-378.

Hyndman RJ, Koehler AB (2006). Another look at measures of forecast
accuracy. *International Journal of Forecasting* 22(4): 679-688.

Tashman LJ (2000). Out-of-sample tests of forecasting accuracy: an
analysis and review. *International Journal of Forecasting* 16(4):
437-450.
