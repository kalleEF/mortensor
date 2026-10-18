# Test differences in forecast accuracy between models (Diebold-Mariano)

Tests whether models of a backtest differ in their expected score with
the Diebold-Mariano test (Diebold and Mariano 1995) and the small-sample
correction of Harvey, Leybourne and Newbold (1997). It complements
[`mt_score_table()`](mt_compare.md), which shows the mean scores: the
test says whether a difference is larger than the variation between
forecast origins.

## Usage

``` r
mt_dm_test(
  backtest,
  ref = NULL,
  metrics = "crps",
  targets = NULL,
  by = NULL,
  alternative = c("two.sided", "less", "greater"),
  varestimator = c("acf", "bartlett")
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

  `"acf"` (sum of the autocovariances, default) or `"bartlett"`
  (Bartlett weights), as in
  [`forecast::dm.test()`](https://pkg.robjhyndman.com/forecast/reference/dm.test.html).

## Value

A `data.table` with one row per model (other than `ref`), value of `by`,
target, metric and horizon: `model`, `ref`, the `by` columns, `target`,
`metric`, `horizon` (a factor: the horizons, then `"all"`), `n_origins`
(number of losses tested), `lag` (largest autocovariance lag used),
`score` and `score_ref` (mean loss over the origins), `mean_diff`
(`score - score_ref`; negative when the model is more accurate),
`statistic` (Diebold-Mariano statistic with the correction of Harvey et
al.) and `p_value`.

## Details

The test uses one loss per forecast origin, the time index of a
rolling-origin backtest. Scores of the cells of one origin and horizon
(ages, strata) concern the same test year and are strongly dependent, so
they are averaged first:

- per horizon `h`: the mean score over the cells of each origin at that
  horizon gives one loss per origin and model. The differences between
  the losses of a model and of `ref` form a series over the origins,
  tested with the autocovariances up to lag `h - 1` (the errors of
  forecasts `h` years ahead from consecutive origins build up over
  overlapping years).

- overall (`horizon = "all"`): the mean of these per-horizon losses over
  all horizons gives one loss per origin, tested with the largest
  horizon `H` (lags up to `H - 1`). Only origins scored at every horizon
  enter.

Only cells scored for both models enter (the same origin, horizon,
stratum, age group, year and target). The statistic, its variance (sum
of autocovariances, or Bartlett weights with
`varestimator = "bartlett"`), the correction factor and the t
distribution with `T - 1` degrees of freedom (`T` origins) are those of
`forecast::dm.test(e1, e2, h = h, power = 1)` with `e1` and `e2` the
losses per origin (scores are not negative, so the loss is the score
itself). With origins spaced `s` years apart (for example
`origins = seq(1981, 2001, by = 10)`), test years overlap only up to lag
`ceiling(h / s) - 1`, which is used instead; with unequal spacing the
smallest gap counts.

With few origins the test has little power: a backtest with 6 origins
finds only large and consistent differences. When the number of origins
`T` is not larger than the number of lags needed (`T <= h` for
consecutive origins), there is no test (`statistic` and `p_value` are
`NA`, with a warning of class `mt_warning_dm`): lag 0 would ignore the
correlation of overlapping test years and give too small p-values. For
long horizons, use more origins or origins spaced `h` years apart. When
the variance estimate with lags up to `h - 1` is not positive, the test
uses lag 0, as
[`forecast::dm.test()`](https://pkg.robjhyndman.com/forecast/reference/dm.test.html)
does, with the same warning. The p-values are not adjusted for the
number of tests (models, targets, horizons); see
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

## References

Diebold FX, Mariano RS (1995). Comparing predictive accuracy. *Journal
of Business & Economic Statistics* 13(3): 253-263.

Harvey D, Leybourne S, Newbold P (1997). Testing the equality of
prediction mean squared errors. *International Journal of Forecasting*
13(2): 281-291.

## See also

[`mt_compare()`](mt_compare.md), [`mt_score_table()`](mt_compare.md),
[`mt_score()`](mt_score.md).

## Examples

``` r
x <- mt_example(deaths = "single", years = 1990:2012)
bt <- mt_compare(list(rw = rwd(), lc = lee_carter()), x, origins = 2001:2007, h = 3, n_draws = 30)
mt_dm_test(bt)
#>     model    ref target metric horizon n_origins   lag      score  score_ref
#>    <char> <char> <char> <char>  <fctr>     <int> <int>      <num>      <num>
#> 1:     lc     rw   rate   crps       1         7     0 0.03573127 0.04750121
#> 2:     lc     rw   rate   crps       2         7     1 0.03613443 0.04969118
#> 3:     lc     rw   rate   crps       3         7     2 0.03729030 0.05588227
#> 4:     lc     rw   rate   crps     all         7     2 0.03638533 0.05102489
#>      mean_diff statistic      p_value
#>          <num>     <num>        <num>
#> 1: -0.01176994 -8.088596 0.0001913927
#> 2: -0.01355675 -6.981889 0.0004294763
#> 3: -0.01859197 -7.064353 0.0004029532
#> 4: -0.01463955 -5.434188 0.0016111375
```
