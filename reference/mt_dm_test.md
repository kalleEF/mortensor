# Test differences in forecast accuracy between models (Diebold-Mariano)

Tests whether models of a backtest differ in their expected score with
the Diebold-Mariano test (Diebold and Mariano 1995). By default the
long-run variance is the equal-weighted cosine (EWC) estimate with
fixed-smoothing critical values from a t distribution, as Lazarus,
Lewis, Stock and Watson (2018) recommend for tests of a single
restriction; the Diebold-Mariano statistic is such a test, of the mean
loss difference. The fixed-smoothing test of Coroneo and Iacone (2020)
for the Diebold-Mariano test, with a Daniell-kernel periodogram
estimate, and the version of
[`forecast::dm.test()`](https://pkg.robjhyndman.com/forecast/reference/dm.test.html),
with the small-sample correction of Harvey, Leybourne and Newbold
(1997), are options. It complements [`mt_score_table()`](mt_compare.md),
which shows the mean scores: the test says whether a difference is
larger than the variation between forecast origins.

## Usage

``` r
mt_dm_test(
  backtest,
  ref = NULL,
  metrics = "crps",
  targets = NULL,
  by = NULL,
  alternative = c("two.sided", "less", "greater"),
  varestimator = c("ewc", "daniell", "acf", "bartlett")
)
```

## Arguments

- backtest:

  An [`mt_backtest()`](mt_backtest.md) or
  [`mt_compare()`](mt_compare.md) result with at least two models.
  Models are identified by `model`, and by `"model (prep)"` when a model
  was backtested with several preprocessing pipelines.

- ref:

  Reference model (default: the first model). Every other model is
  compared with it.

- metrics:

  Scores to compare: `"crps"`, `"wis"`, `"ae"` (absolute error of the
  median) and `"se"` (squared error of the median). CRPS and WIS compare
  whole predictive distributions; `"ae"` and `"se"` the medians only.

- targets:

  Targets to test (default: all in the backtest; see `quantities` and
  `contrast` in [`mt_backtest()`](mt_backtest.md)).

- by:

  Optional stratum columns (for example `"sex"`): the losses per origin
  and the tests are then computed separately for each of their values.

- alternative:

  `"two.sided"`, `"less"` (the model has a lower expected score than
  `ref`, that is, it is more accurate) or `"greater"`.

- varestimator:

  `"ewc"` (equal-weighted cosine estimate with fixed-smoothing t
  critical values, default), `"daniell"` (Daniell-kernel periodogram
  estimate with fixed-smoothing t critical values, Coroneo and Iacone
  2020), `"acf"` (sum of the autocovariances, as
  [`forecast::dm.test()`](https://pkg.robjhyndman.com/forecast/reference/dm.test.html))
  or `"bartlett"` (Bartlett weights, as
  [`forecast::dm.test()`](https://pkg.robjhyndman.com/forecast/reference/dm.test.html)).
  See Variance estimators.

## Value

A `data.table` with one row per model (other than `ref`), value of `by`,
target, metric and horizon: `model`, `ref`, the `by` columns, `target`,
`metric`, `horizon` (a factor: the horizons, then `"all"`), `n_origins`
(number of losses tested), `lag` (largest autocovariance lag used; `NA`
for `"ewc"` and `"daniell"`), `df` (degrees of freedom of the t
distribution), `score` and `score_ref` (mean loss over the origins),
`mean_diff` (`score - score_ref`; negative when the model is more
accurate), `statistic` and `p_value`.

## Details

The test uses one loss per forecast origin, the time index of a
rolling-origin backtest. Scores of the cells of one origin and horizon
(ages, strata) concern the same test year and are strongly dependent, so
they are averaged first:

- per horizon `h`: the mean score over the cells of each origin at that
  horizon gives one loss per origin and model. The differences between
  the losses of a model and of `ref` form a series over the origins. The
  errors of forecasts `h` years ahead from consecutive origins build up
  over overlapping years, so the differences are autocorrelated.

- overall (`horizon = "all"`): the mean of these per-horizon losses over
  all horizons gives one loss per origin. Only origins scored at every
  horizon enter.

Only cells scored for both models enter (the same origin, horizon,
stratum, age group, year and target). Scores are not negative, so the
loss is the score itself.

## Variance estimators

With `d_t` the loss differences of the `T` origins (`t = 1, ..., T`, in
the order of the origins) and `dbar` their mean:

- `"ewc"` (default):
  `Lambda_j = sqrt(2 / T) * sum_t (d_t - dbar) * cos(pi * j * (t - 1/2) / T)`
  for `j = 1, ..., nu`, the long-run variance
  `Omega = mean(Lambda_j^2)`, the statistic
  `sqrt(T) * dbar / sqrt(Omega)` and p-values from the t distribution
  with `nu` degrees of freedom, `nu = 0.4 * T^(2/3)` cosine terms as
  Lazarus et al. (2018) recommend for the EWC test, rounded down and at
  least 1 (1 for up to 11 origins, 2 for 12 to 20, 3 for 21 to 31).
  There is no further small-sample correction. The estimate does not use
  the horizon, so there is a test for every horizon with at least 2
  origins; `lag` is `NA`.

- `"daniell"`: the weighted periodogram estimate with the Daniell kernel
  of Coroneo and Iacone (2020, Eq. 8), `Omega = (1 / m) * sum_j I_j` for
  the first `m` Fourier frequencies `lambda_j = 2 * pi * j / T`, with
  `I_j = |sum_t d_t exp(-i lambda_j t)|^2 / T`, and p-values from the t
  distribution with `2 m` degrees of freedom (their Eq. 9).
  `m = floor(T^(1/3))`, at least 1, as in their simulations and
  application (Section 5: the best combination of size and power in a
  preliminary study): 1 for up to 7 origins, 2 for 8 to 26. `lag` is
  `NA`. Coroneo and Iacone (2020) also propose the Bartlett-kernel
  estimate with bandwidth `floor(T^(1/2))` and fixed-b critical values
  given as a polynomial for two quantiles, which give no p-values; it is
  not implemented. They do not use the equal-weighted cosine estimate.

- `"acf"`: the test of `forecast::dm.test(e1, e2, h = h, power = 1)`
  with `e1` and `e2` the losses per origin: the variance is the sum of
  the autocovariances up to lag `h - 1`, the statistic has the
  correction of Harvey et al. (1997), and p-values come from the t
  distribution with `T - 1` degrees of freedom. With origins spaced `s`
  years apart (for example `origins = seq(1981, 2001, by = 10)`), test
  years overlap only up to lag `ceiling(h / s) - 1`, which is used
  instead; with unequal spacing the smallest gap counts. When `T` is not
  larger than the number of lags needed (`T <= h` for consecutive
  origins), or when the variance estimate is not positive, there is no
  test (`statistic` and `p_value` are `NA`, with a warning of class
  `mt_warning_dm`);
  [`forecast::dm.test()`](https://pkg.robjhyndman.com/forecast/reference/dm.test.html)
  falls back to lag 0 in the second case, which rejects a true null
  hypothesis about half the time.

- `"bartlett"`: as `"acf"` with Bartlett weights `1 - k / h` for the
  autocovariances at lags `k = 1, ..., h - 1`, as in
  [`forecast::dm.test()`](https://pkg.robjhyndman.com/forecast/reference/dm.test.html).

## Size and power with few origins

Backtests have few origins, and with consecutive origins the test years
of horizons above 1 overlap. In a simulation under the null hypothesis
(20,000 series of loss differences that are sums of `h` consecutive
independent normal errors, a moving average of order `h - 1` as for
forecasts `h` years ahead from consecutive origins), the `"acf"` test at
the nominal 5% level rejected in 10% (`h = 2`), 16% (`h = 3`) and 19%
(`h = 5`) of the samples with 7 origins, and in 16% to 17% with 15 or 16
origins at `h = 5`; `"bartlett"` in 10% to 15% with 7 origins. Their
p-values are too small with few overlapping origins. The fallback to lag
0 of
[`forecast::dm.test()`](https://pkg.robjhyndman.com/forecast/reference/dm.test.html),
needed in 4% of the samples with 7 origins at `h = 3` and in 22% at
`h = 5`, rejected in 39% to 55% of them, so that the test with the
fallback rejected in 27% at `h = 5` with 7 origins. The EWC test
rejected in 4.9% to 7.2% of the samples with 7 origins and in 4.8% to
7.0% with 15 or 16 origins (`h` from 1 to 5). Its size is nearly
correct, but with few origins it has few degrees of freedom and low
power: with 7 origins (`nu = 1`) it finds only large and consistent
differences. The `"daniell"` test of Coroneo and Iacone (2020) has more
degrees of freedom (2 with 7 origins, 4 with 8 to 26) and rejected in
4.9%, 6.0%, 7.6%, 10.5% and 15.3% of the samples with 7 origins at `h`
from 1 to 5, and in 4.8% to 13.8% with 15 or 16 origins: its size is
close to 5% for `h` up to 2, where test years overlap little, and too
large beyond. Coroneo and Iacone (2020) studied one-step forecasts with
40 or 120 origins. More origins, or origins spaced `h` years apart,
which avoid overlapping test years, give more informative tests. The
p-values are not adjusted for the number of tests (models, targets,
horizons); see
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

## References

Coroneo L, Iacone F (2020). Comparing predictive accuracy in small
samples using fixed-smoothing asymptotics. *Journal of Applied
Econometrics* 35(4): 391-409.
[doi:10.1002/jae.2756](https://doi.org/10.1002/jae.2756)

Diebold FX, Mariano RS (1995). Comparing predictive accuracy. *Journal
of Business & Economic Statistics* 13(3): 253-263.

Harvey D, Leybourne S, Newbold P (1997). Testing the equality of
prediction mean squared errors. *International Journal of Forecasting*
13(2): 281-291.

Lazarus E, Lewis DJ, Stock JH, Watson MW (2018). HAR inference:
recommendations for practice. *Journal of Business & Economic
Statistics* 36(4): 541-559.

## See also

[`mt_compare()`](mt_compare.md), [`mt_score_table()`](mt_compare.md),
[`mt_score()`](mt_score.md).

## Examples

``` r
x <- mt_example(deaths = "single", years = 1990:2012)
bt <- mt_compare(list(rw = rwd(), lc = lee_carter()), x, origins = 2001:2007, h = 3, n_draws = 30)
mt_dm_test(bt)
#>     model    ref target metric horizon n_origins   lag    df      score
#>    <char> <char> <char> <char>  <fctr>     <int> <int> <num>      <num>
#> 1:     lc     rw   rate   crps       1         7    NA     1 0.03593905
#> 2:     lc     rw   rate   crps       2         7    NA     1 0.03653935
#> 3:     lc     rw   rate   crps       3         7    NA     1 0.03743453
#> 4:     lc     rw   rate   crps     all         7    NA     1 0.03663764
#>     score_ref   mean_diff statistic    p_value
#>         <num>       <num>     <num>      <num>
#> 1: 0.04750121 -0.01156216 -4.218645 0.14817141
#> 2: 0.04969118 -0.01315182 -5.091841 0.12345620
#> 3: 0.05588227 -0.01844774 -6.541768 0.09656861
#> 4: 0.05102489 -0.01438724 -5.300051 0.11872013
# the version of forecast::dm.test()
mt_dm_test(bt, varestimator = "acf")
#>     model    ref target metric horizon n_origins   lag    df      score
#>    <char> <char> <char> <char>  <fctr>     <int> <int> <num>      <num>
#> 1:     lc     rw   rate   crps       1         7     0     6 0.03593905
#> 2:     lc     rw   rate   crps       2         7     1     6 0.03653935
#> 3:     lc     rw   rate   crps       3         7     2     6 0.03743453
#> 4:     lc     rw   rate   crps     all         7     2     6 0.03663764
#>     score_ref   mean_diff statistic      p_value
#>         <num>       <num>     <num>        <num>
#> 1: 0.04750121 -0.01156216 -8.165424 0.0001815895
#> 2: 0.04969118 -0.01315182 -5.469997 0.0015579200
#> 3: 0.05588227 -0.01844774 -7.156147 0.0003756169
#> 4: 0.05102489 -0.01438724 -5.214971 0.0019858300
```
