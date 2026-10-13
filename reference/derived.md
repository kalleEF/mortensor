# Derived quantities from forecast draws or observed data

All quantities are computed path by path, so their intervals reflect the
joint uncertainty across ages and years. For observed data (an
[mt_data](mt_data.md) object) the paths are Poisson draws of the
observed death counts (see [`mt_observed()`](mt_observed.md)), so the
intervals are those of the observed rates.

## Usage

``` r
prob_dying(
  draws,
  from,
  to,
  probs = c(0.025, 0.1, 0.5, 0.9, 0.975),
  paths = FALSE
)

asr(
  draws,
  standard = "esp2013",
  probs = c(0.025, 0.1, 0.5, 0.9, 0.975),
  paths = FALSE,
  per = 1e+05
)

yll(
  draws,
  exposure = NULL,
  reference_ex = NULL,
  probs = c(0.025, 0.1, 0.5, 0.9, 0.975)
)
```

## Arguments

- draws:

  An [mt_draws](mt_draws.md) object; `prob_dying()` and `asr()` also
  take an [mt_data](mt_data.md) object for observed rates (with Poisson
  intervals, see [`mt_observed()`](mt_observed.md)).

- from, to:

  Ages (starts of grid intervals; `to` may be the end of the last closed
  interval).

- probs:

  Quantile levels.

- paths:

  Return the values of every path, an array `[period, stratum, draw]`,
  instead of the summary (for example to compute ratios between strata
  path by path).

- standard:

  The name of a standard population (see
  [`standard_population()`](standard_population.md)) or a data frame
  with `age_lo`, `age_hi`, `weight`.

- per:

  Rates per `per` persons (age-standardised rates).

- exposure:

  Future exposures for `yll()` (see [`mt_deaths()`](mt_deaths.md)).

- reference_ex:

  Reference remaining life expectancy per grid age, or `NULL` for the
  forecast's own.

## Value

A `data.table` with stratum columns, `period`, `mean` and quantiles, or
with `paths = TRUE` an array.

## Details

- `prob_dying()`: probability of dying between ages `from` and `to` (for
  example the SDG indicator 3.4.1 uses 30 to 70), `1 - l(to) / l(from)`
  from the period life table;

- `asr()`: age-standardised death rate per `per` persons (default
  100,000) with weights of a standard population, restricted to the ages
  of the grid and renormalised (columns `ages`, `standard` and `per`
  state this; an ASR for ages 50+ is per 100,000 of the standard
  population aged 50+);

- `yll()`: years of life lost, deaths times a reference remaining life
  expectancy at each age (by default the forecast's own period life
  expectancy at the interval midpoint, approximated by the value at the
  start of the interval).

## Examples

``` r
x <- mt_example(deaths = "single")
fc <- mt_forecast(mt_fit(lee_carter(), x), h = 5, n_draws = 100, seed = 1)
prob_dying(fc, from = 30, to = 70)
#>     stratum    sex  label period       mean     q0.025       q0.1       q0.5
#>       <int> <char> <char>  <num>      <num>      <num>      <num>      <num>
#>  1:       1 female female   2020 0.10272149 0.10029047 0.10123483 0.10277695
#>  2:       1 female female   2021 0.10122205 0.09770717 0.09886535 0.10123556
#>  3:       1 female female   2022 0.09998290 0.09583266 0.09721790 0.10003892
#>  4:       1 female female   2023 0.09874500 0.09373185 0.09532758 0.09863632
#>  5:       1 female female   2024 0.09748585 0.09240964 0.09379196 0.09761577
#>  6:       2   male   male   2020 0.16304307 0.15931262 0.16084192 0.16311840
#>  7:       2   male   male   2021 0.16106775 0.15598329 0.15810202 0.16111791
#>  8:       2   male   male   2022 0.15872591 0.15324347 0.15464245 0.15861433
#>  9:       2   male   male   2023 0.15665954 0.14982229 0.15207823 0.15679195
#> 10:       2   male   male   2024 0.15473156 0.14624807 0.14856198 0.15533294
#>          q0.9    q0.975  from    to
#>         <num>     <num> <num> <num>
#>  1: 0.1040814 0.1052221    30    70
#>  2: 0.1033245 0.1045844    30    70
#>  3: 0.1031898 0.1041283    30    70
#>  4: 0.1022911 0.1034675    30    70
#>  5: 0.1009225 0.1032426    30    70
#>  6: 0.1651576 0.1659756    30    70
#>  7: 0.1639472 0.1660846    30    70
#>  8: 0.1626083 0.1654740    30    70
#>  9: 0.1611459 0.1633972    30    70
#> 10: 0.1593443 0.1627950    30    70
asr(fc)
#>     stratum    sex  label period      mean    q0.025      q0.1      q0.5
#>       <int> <char> <char>  <num>     <num>     <num>     <num>     <num>
#>  1:       1 female female   2020  917.8167  900.9276  907.4972  918.2079
#>  2:       1 female female   2021  907.3934  882.9148  890.9984  907.5024
#>  3:       1 female female   2022  898.7630  869.8038  879.4961  899.1765
#>  4:       1 female female   2023  890.1245  855.0679  866.2651  889.4009
#>  5:       1 female female   2024  881.3295  845.7694  855.4901  882.2764
#>  6:       2   male   male   2020 1532.8800 1504.1236 1515.9179 1533.4636
#>  7:       2   male   male   2021 1517.6510 1478.4268 1494.7835 1518.0457
#>  8:       2   male   male   2022 1499.5829 1457.2561 1468.0689 1498.7367
#>  9:       2   male   male   2023 1483.6304 1430.7873 1448.2452 1484.6712
#> 10:       2   male   male   2024 1468.7305 1403.0915 1421.0269 1473.4035
#>          q0.9    q0.975    ages standard   per
#>         <num>     <num>  <char>   <char> <num>
#>  1:  927.2516  935.1484 30-100+  esp2013 1e+05
#>  2:  922.0057  930.7354 30-100+  esp2013 1e+05
#>  3:  921.0720  927.5761 30-100+  esp2013 1e+05
#>  4:  914.8358  922.9969 30-100+  esp2013 1e+05
#>  5:  905.3257  921.4380 30-100+  esp2013 1e+05
#>  6: 1549.1706 1555.4688 30-100+  esp2013 1e+05
#>  7: 1539.8488 1556.3079 30-100+  esp2013 1e+05
#>  8: 1529.5330 1551.6067 30-100+  esp2013 1e+05
#>  9: 1518.2615 1535.6114 30-100+  esp2013 1e+05
#> 10: 1504.3681 1530.9721 30-100+  esp2013 1e+05
```
