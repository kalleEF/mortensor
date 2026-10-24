# Summarise backtest scores

Summarise backtest scores

## Usage

``` r
mt_score(
  backtest,
  by = "model",
  target = c("rate", "life_expectancy", "asr", "life_expectancy_gap", "asr_gap"),
  measures = NULL
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

- measures:

  Optional further summaries (none by default):

  - `"sape"`: the mean standardised absolute predictive error of Raftery
    et al. (2013), `sqrt(pi / 2) |y - median| / sd` with the standard
    deviation `sd` of the predictive draws. It is scaled so that its
    expected value is 1 when the predictive distribution is correct and
    normal; a mean above 1 indicates intervals that are too narrow,
    below 1 too wide. Raftery et al. define it for life expectancy, so
    it is available for `target = "life_expectancy"` only. Rows whose
    draws have no spread (`sd = 0`) are left out;

  - `"cov90"`: the empirical coverage of the central 90% intervals.
    Raftery et al. (2013) report the coverage of the 80%, 90% and 95%
    intervals with their mean half-widths; with
    `by = c("model", "horizon")` they are given by horizon;

  - `"halfwidth"`: the mean half-widths `hw80`, `hw90`, `hw95` of these
    intervals (Raftery et al. 2013, Table 1);

  - `"revision"`: the mean absolute revision, the absolute change of the
    forecast median of a cell (stratum, age group, year) between
    consecutive origins of the backtest, for every cell forecast from
    both. It is a summary of this package (reference unknown).
    Stoeldraijer et al. (2013) name robustness as a purpose of a
    forecast: small changes in the observations when the data are
    updated should result in only modest changes in the forecasts. A
    revision belongs to the later origin and its horizon; the first
    origin has none. Origins further apart than one year (for example
    `origins = c(2005, 2007)`) give revisions over several years of
    data. The published measure of robustness is `"robustness"`;

  - `"robustness"`: the robustness of Stoeldraijer et al. (2018), "the
    standard deviation (SD) of the increase/decrease" of the forecast of
    one target year over successive forecasts, each with one more year
    of data; a lower standard deviation means a more robust forecast.
    They forecast life expectancy at age 65 in 2020 from the fitting
    periods 1960-2005, ..., 1960-2014. Here, for each cell (stratum, age
    group, target year and target), it is the sample standard deviation
    ([`stats::sd()`](https://rdrr.io/r/stats/sd.html)) of the revisions
    of the forecast median between consecutive origins (as for
    `"revision"`), for cells with at least two revisions (`NA`
    otherwise); the column `robustness` is the mean of these standard
    deviations over the cells of a group, each cell counted once. With
    one stratum, `target = "life_expectancy"` and one target year
    (`by = c("model", "period")`) this is the published measure; the
    mean over several cells is package-specific (reference unknown). The
    revisions of a cell span several origins and horizons, so `by`
    cannot contain `"horizon"` or `"origin"`: group by the target year
    (`"period"`). The target years of a backtest are observed years, so
    the origins and `h` must reach them from every origin: for example
    `origins = 2005:2014` and `h = 10` forecast 2015 from all ten
    origins, one year apart like the fitting periods of the paper. The
    life expectancy of a backtest is that at the first age of the grid,
    so at age 65 for data from age 65.

## Value

A `data.table` with the mean CRPS, WIS and absolute error, the root mean
squared error of the median (`rmse`), MASE (mean of absolute error over
the scale, the mean absolute year-to-year change of the training data,
as the mean absolute scaled error of Hyndman and Koehler 2006; rates
only), empirical coverage of the 50%, 80% and 95% intervals, mean PIT,
the number of scored cells (`n`), the number of them whose observed
deaths were floored at half a death (`n_floored`, see
[`mt_backtest()`](mt_backtest.md); `NA` for targets other than `"rate"`
and for backtest tables without the column `floored`) and the `unit` of
CRPS, WIS and absolute error (log rates, or years for life expectancy).
The scores of floored cells can dominate the means of sparse data; the
derived quantities and coarser age groups avoid the floor. Backtests of
the same model with different preprocessing (column `prep`) are scored
separately. With `measures`, the columns `sape`, `cov90`, `hw80`,
`hw90`, `hw95`, `revision` and `robustness` follow.

## References

Hyndman RJ, Koehler AB (2006). Another look at measures of forecast
accuracy. *International Journal of Forecasting* 22(4): 679-688.

Raftery AE, Chunn JL, Gerland P, Sevcikova H (2013). Bayesian
probabilistic projections of life expectancy for all countries.
*Demography* 50(3): 777-801.
[doi:10.1007/s13524-012-0193-x](https://doi.org/10.1007/s13524-012-0193-x)

Stoeldraijer L, van Duin C, van Wissen L, Janssen F (2013). Impact of
different mortality forecasting methods and explicit assumptions on
projected future life expectancy: the case of the Netherlands.
*Demographic Research* 29(13): 323-354.
[doi:10.4054/DemRes.2013.29.13](https://doi.org/10.4054/DemRes.2013.29.13)

Stoeldraijer L, van Duin C, van Wissen L, Janssen F (2018). Comparing
strategies for matching mortality forecasts to the most recently
observed data: exploring the trade-off between accuracy and robustness.
*Genus* 74: 16.
[doi:10.1186/s41118-018-0041-y](https://doi.org/10.1186/s41118-018-0041-y)

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
#>        mase     cov50     cov80     cov95  pit_mean     n n_floored     unit
#>       <num>     <num>     <num>     <num>     <num> <int>     <int>   <char>
#> 1: 1.048261 0.4859155 0.7253521 0.8661972 0.5020423   142         0 log rate
#> 2: 1.157779 0.6126761 0.8591549 0.9718310 0.5157746   142         0 log rate
#> 3: 1.386932 0.6478873 0.8943662 0.9647887 0.4927465   142         0 log rate
#> 4: 1.287454 0.7183099 0.9577465 1.0000000 0.4763380   142         0 log rate
#> 5: 1.459025 0.7605634 0.9507042 1.0000000 0.4895070   142         0 log rate
# life expectancy: SAPE and 90% coverage by horizon (Raftery et al. 2013),
# revisions between consecutive origins
bl <- mt_backtest(rwd(), x, origins = 2004:2007, h = 5, n_draws = 100,
                  quantities = "life_expectancy")
mt_score(bl, by = c("model", "horizon"), target = "life_expectancy",
         measures = c("sape", "cov90", "revision"))
#>                     model horizon       crps        wis        mae       rmse
#>                    <char>   <num>      <num>      <num>      <num>      <num>
#> 1: Random walk with drift       1 0.05241432 0.04643019 0.06818529 0.09138844
#> 2: Random walk with drift       2 0.05633612 0.05021002 0.07538965 0.08831333
#> 3: Random walk with drift       3 0.06926703 0.06167666 0.09661399 0.10912573
#> 4: Random walk with drift       4 0.08127709 0.07284445 0.09343235 0.11813568
#> 5: Random walk with drift       5 0.09710010 0.08682756 0.12023575 0.16650060
#>     mase cov50 cov80 cov95 pit_mean     n n_floored   unit      sape cov90
#>    <num> <num> <num> <num>    <num> <int>     <int> <char>     <num> <num>
#> 1:    NA 0.625 0.750     1  0.52750     8        NA  years 0.7986443     1
#> 2:    NA 0.625 1.000     1  0.55875     8        NA  years 0.6208891     1
#> 3:    NA 0.750 1.000     1  0.50625     8        NA  years 0.6485447     1
#> 4:    NA 0.750 1.000     1  0.53250     8        NA  years 0.5303456     1
#> 5:    NA 0.750 0.875     1  0.55875     8        NA  years 0.5496038     1
#>      revision
#>         <num>
#> 1: 0.09423038
#> 2: 0.10963108
#> 3: 0.12096127
#> 4: 0.13076723
#> 5:         NA
# robustness (Stoeldraijer et al. 2018) of each target year: 2008 and
# 2009 are forecast from all four origins
mt_score(bl, by = c("model", "period"), target = "life_expectancy",
         measures = "robustness")
#>                     model period       crps        wis        mae       rmse
#>                    <char>  <num>      <num>      <num>      <num>      <num>
#> 1: Random walk with drift   2005 0.03567047 0.03165998 0.05139350 0.06238901
#> 2: Random walk with drift   2006 0.07705623 0.06836422 0.12534838 0.12792663
#> 3: Random walk with drift   2007 0.05311679 0.04720499 0.06603429 0.07483395
#> 4: Random walk with drift   2008 0.05608048 0.05022010 0.05675972 0.07472110
#> 5: Random walk with drift   2009 0.06191777 0.05566069 0.05994773 0.07285361
#> 6: Random walk with drift   2010 0.09545905 0.08509163 0.13337241 0.16187205
#> 7: Random walk with drift   2011 0.08459577 0.07581569 0.10584468 0.13570146
#> 8: Random walk with drift   2012 0.14888367 0.13152272 0.23659860 0.25792984
#>     mase     cov50 cov80 cov95 pit_mean     n n_floored   unit robustness
#>    <num>     <num> <num> <num>    <num> <int>     <int> <char>      <num>
#> 1:    NA 1.0000000   1.0     1  0.45500     2        NA  years         NA
#> 2:    NA 0.2500000   0.5     1  0.85250     4        NA  years         NA
#> 3:    NA 0.6666667   1.0     1  0.49000     6        NA  years  0.1367226
#> 4:    NA 0.8750000   1.0     1  0.49875     8        NA  years  0.1383298
#> 5:    NA 1.0000000   1.0     1  0.49500     8        NA  years  0.1546668
#> 6:    NA 0.3333333   1.0     1  0.48000     6        NA  years  0.2498016
#> 7:    NA 0.7500000   1.0     1  0.54500     4        NA  years         NA
#> 8:    NA 0.5000000   0.5     1  0.60000     2        NA  years         NA
```
