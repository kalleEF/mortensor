# Summarise backtest scores

Summarise backtest scores

## Usage

``` r
mt_score(
  backtest,
  by = "model",
  target = c("rate", "life_expectancy", "asr", "life_expectancy_gap", "asr_gap")
)
```

## Arguments

- backtest:

  An [`mt_backtest()`](mt_backtest.md) result (or several combined with
  [`rbind()`](https://rdrr.io/r/base/cbind.html)).

- by:

  Grouping columns, for example `c("model", "horizon")`.

- target:

  `"rate"` (log death rates per cell, default), `"life_expectancy"` or
  `"asr"` (when the backtest scored these quantities, see `quantities`
  in [`mt_backtest()`](mt_backtest.md)), or `"life_expectancy_gap"` or
  `"asr_gap"` (differences between strata, see `contrast` in
  [`mt_backtest()`](mt_backtest.md)).

## Value

A `data.table` with the mean CRPS, WIS and absolute error, the root mean
squared error of the median (`rmse`), MASE (mean of absolute error over
the scale; rates only), empirical coverage of the 50%, 80% and 95%
intervals, mean PIT, the number of scored cells and the `unit` of CRPS,
WIS and absolute error (log rates, or years for life expectancy).
Backtests of the same model with different preprocessing (column `prep`)
are scored separately.

## Examples

``` r
x <- mt_example(deaths = "single", years = 1990:2012)
bt <- mt_backtest(rwd(), x, origins = 2007, h = 5, n_draws = 100)
mt_score(bt, by = c("model", "horizon"))
#>                     model horizon       crps        wis        mae       rmse
#>                    <char>   <num>      <num>      <num>      <num>      <num>
#> 1: Random walk with drift       1 0.04454219 0.03992641 0.06272241 0.09469268
#> 2: Random walk with drift       2 0.04654157 0.04148880 0.06166271 0.09425197
#> 3: Random walk with drift       3 0.05565911 0.04972059 0.07710375 0.11415199
#> 4: Random walk with drift       4 0.05812890 0.05190181 0.07358098 0.11643492
#> 5: Random walk with drift       5 0.06429823 0.05740455 0.08308110 0.12352441
#>        mase     cov50     cov80     cov95  pit_mean     n     unit
#>       <num>     <num>     <num>     <num>     <num> <int>   <char>
#> 1: 1.048261 0.4859155 0.7253521 0.8661972 0.5020423   142 log rate
#> 2: 1.157779 0.6126761 0.8591549 0.9718310 0.5157746   142 log rate
#> 3: 1.386932 0.6478873 0.8943662 0.9647887 0.4927465   142 log rate
#> 4: 1.287454 0.7183099 0.9577465 1.0000000 0.4763380   142 log rate
#> 5: 1.459025 0.7605634 0.9507042 1.0000000 0.4895070   142 log rate
```
