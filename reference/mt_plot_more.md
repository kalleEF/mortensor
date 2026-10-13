# Further plots: Tucker factors, cause shares and backtest scores

- `mt_plot_tucker()`: the factor matrices of a [`tucker3()`](tucker3.md)
  fit (age, period and stratum modes), one panel per mode and one line
  per factor; signs are arbitrary.

- `mt_plot_shares()`: shares of each cause in the deaths over the
  selected ages, observed and forecast (median and 95% interval), from a
  forecast with a cause variable among the strata (for example from
  [`cause_shares()`](cause_shares.md)).

- `mt_plot_scores()`: mean score of each model by horizon or age from a
  backtest, or one bar per model (`by = "model"`); `target` selects
  rates or a derived quantity scored by the backtest.

## Usage

``` r
mt_plot_tucker(fit)

mt_plot_shares(
  draws,
  var = NULL,
  total_level = NULL,
  ages = NULL,
  groups = NULL,
  causes = NULL
)

mt_plot_scores(
  backtest,
  metric = "crps",
  by = c("horizon", "age", "model"),
  target = c("rate", "life_expectancy", "asr")
)
```

## Arguments

- fit:

  A fitted [`tucker3()`](tucker3.md) model.

- draws:

  An [mt_draws](mt_draws.md) object whose strata include `var`.

- var:

  Stratum variable with the causes (default: the variable with shared
  exposures in the hierarchy of a reconciled forecast).

- total_level:

  Level of `var` for all causes (left out of the shares; default: its
  total in the hierarchy).

- ages:

  Age labels summed over (default all).

- groups:

  Panels to show: combinations of the other stratum variables joined by
  "." (for example `"female.Q1"`; default all).

- causes:

  Levels of `var` to show (default all; shares are always computed among
  all levels).

- backtest:

  An [`mt_backtest()`](mt_backtest.md) result.

- metric:

  Score column, for example `"crps"`, `"wis"` or `"ae"`.

- by:

  `"horizon"`, `"age"` or `"model"`.

- target:

  `"rate"`, `"life_expectancy"` or `"asr"` (see
  [`mt_backtest()`](mt_backtest.md)).

## Value

A ggplot object.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  x <- mt_example(deaths = "single", years = 2000:2019)
  mt_plot_tucker(mt_fit(tucker3(), x))
  bt <- rbind(mt_backtest(rwd(), x, origins = c(2012, 2014), h = 4, n_draws = 50),
              mt_backtest(lee_carter(), x, origins = c(2012, 2014), h = 4, n_draws = 50))
  mt_plot_scores(bt)
}
```
