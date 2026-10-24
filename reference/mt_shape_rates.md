# Age shapes within groups for disaggregating forecasts

Rates on a fine age grid that give
[`mt_disaggregate()`](mt_disaggregate.md) the age pattern within each
age group of a forecast. The method (`"last"`) takes the rates of the
last `years` years (deaths and exposures summed), the same for every
forecast year: a standard schedule for proportional allocation within
groups.

## Usage

``` r
mt_shape_rates(data, periods, method = "last", years = 5)
```

## Arguments

- data:

  An [mt_data](mt_data.md) object at the fine age resolution (for
  example ungrouped training data).

- periods:

  Forecast years.

- method:

  `"last"` (the rates of the last `years` years).

- years:

  Number of last observed years used.

## Value

An array `[age, period, stratum]` of rates.

## Examples

``` r
x <- mt_example(years = 2000:2019, deaths = "single")
sh <- mt_shape_rates(x, periods = 2020:2022)
xg <- mt_regroup(x, age = c(seq(30, 85, 5), 90))
fc <- mt_forecast(mt_fit(lee_carter(), xg), h = 3, n_draws = 20, seed = 1)
fs <- mt_disaggregate(fc, mt_grid(x)$age, exposure = mt_array(x, "exposure")[, 20, ],
                      method = "shape", shape_rates = sh)
```
