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
  per = 1e+05,
  exposure = NULL,
  age_groups = c("standard", "grid")
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

  Exposures (see [`mt_deaths()`](mt_deaths.md)): for `yll()` the future
  exposures that give the deaths; for `asr()` the exposures that weight
  the rates within each standard age group (default: those of each
  observed year for observed data, of the last training year for
  forecasts).

- age_groups:

  `"standard"` (default: direct standardisation on the age groups of the
  standard population) or `"grid"` (the standard weights shared among
  the grid cells, specific to this package). See Age-standardised rates.

- reference_ex:

  Reference remaining life expectancy at the start of each grid
  interval, or `NULL` for the forecast's own.

## Value

A `data.table` with stratum columns, `period`, `mean` and quantiles, or
with `paths = TRUE` an array.

## Details

- `prob_dying()`: probability of dying between ages `from` and `to` (for
  example the SDG indicator 3.4.1 uses 30 to 70), `1 - l(to) / l(from)`
  from the period life table of [`life_table()`](life_table.md)
  (Preston, Heuveline and Guillot 2001);

- `asr()`: age-standardised death rate per `per` persons (default
  100,000), by direct standardisation with the weights of a standard
  population (see Age-standardised rates);

- `yll()`: years of life lost, the deaths of each age interval times the
  remaining life expectancy at the average age at death in the interval
  (see Years of life lost).

## Age-standardised rates

By default (`age_groups = "standard"`) the rates are first combined to
the age groups of the standard population, as for the European Standard
Population (Eurostat 2013): the rate of a standard age group is the sum
of the deaths divided by the sum of the exposures of the grid cells in
it, that is the mean of their rates weighted by their exposures. The ASR
is the mean of these rates weighted by the standard population. The
exposures are those of each observed year for observed data and, for
forecasts, `exposure` (default: the exposures of the last training year,
as in [`mt_deaths()`](mt_deaths.md)). Rates of a group without exposures
(all zero, or draws that hold none) are averaged with equal weights.

Grid cells wider than the standard groups (for example 0-4 for the
groups 0 and 1-4 of the ESP 2013, or an open interval 90+ for the groups
90-94 and 95+) take the sum of the weights of the standard groups they
cover: rates and weights are combined over the smallest age groups that
are unions both of grid cells and of standard groups. Only the ages of
the grid enter, and the weights are renormalised to them (columns
`ages`, `standard` and `per` state this; an ASR for ages 50+ is per
100,000 of the standard population aged 50+); a standard group that the
grid covers only in part enters with the covered share of its weight.

`age_groups = "grid"` instead shares the weight of each standard group
evenly among the grid cells within it (the weight of the open standard
group equally among the grid cells it covers) and weights the rates of
the grid cells with these shares. This is specific to this package
(reference unknown). Where the grid is finer than the standard groups,
every cell of a group gets the same weight, although the older cells of
a group have higher death rates and, at older ages, fewer persons, so
the ASR is higher than with the standard's age groups: by 1.8% to 3.6%
for French rates on single years of age to 100+ (HMD, 1950-2006, ESP
2013), and by more when the grid extends to higher ages, over which the
weight of the open standard group 95+ is then spread evenly (up to 11%
with single years to 110+). On a grid of the standard's age groups (or
coarser) both give the same value.

## Years of life lost

`yll()` multiplies the deaths of each age interval by the remaining life
expectancy at the average age at death in the interval, `x + a_x`,
interpolated linearly between the values at the start and at the end of
the interval, `e_x + a_x / n_x * (e_(x+n) - e_x)` (as for e-dagger in
[`lifespan_variation()`](lifespan_variation.md)); the last interval uses
its `e_x`. With the life expectancy of the forecast's own period life
table (default) these are the period expected years of life lost of
Murray (1994). With `reference_ex`, the remaining life expectancy of a
standard life table at the start of each grid interval, interpolated in
the same way with the `a_x` of the forecast's life table, they are the
standard expected years of life lost of the same paper. Neither
discounting nor age weighting is applied. The deaths are rates times
`exposure` (default: the exposures of the last training year), without
Poisson noise.

## References

Eurostat (2013). Revision of the European Standard Population: report of
Eurostat's task force. Methodologies and Working papers. Publications
Office of the European Union, Luxembourg.
[doi:10.2785/11470](https://doi.org/10.2785/11470)

Murray CJL (1994). Quantifying the burden of disease: the technical
basis for disability-adjusted life years. *Bulletin of the World Health
Organization* 72(3): 429-445.

Preston SH, Heuveline P, Guillot M (2001). *Demography: Measuring and
Modeling Population Processes*. Blackwell.

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
#>  1:       1 female female   2020  905.4925  888.7426  895.2578  905.8803
#>  2:       1 female female   2021  895.1552  870.8804  878.8962  895.2629
#>  3:       1 female female   2022  886.5967  857.8801  867.4904  887.0061
#>  4:       1 female female   2023  878.0306  843.2699  854.3715  877.3120
#>  5:       1 female female   2024  869.3095  834.0513  843.6885  870.2473
#>  6:       2   male   male   2020 1493.9369 1465.5919 1477.2167 1494.5118
#>  7:       2   male   male   2021 1478.9260 1440.2668 1456.3864 1479.3140
#>  8:       2   male   male   2022 1461.1183 1419.4051 1430.0597 1460.2826
#>  9:       2   male   male   2023 1445.3970 1393.3261 1410.5264 1446.4205
#> 10:       2   male   male   2024 1430.7146 1366.0425 1383.7105 1435.3166
#>          q0.9    q0.975    ages standard   per
#>         <num>     <num>  <char>   <char> <num>
#>  1:  914.8499  922.6825 30-100+  esp2013 1e+05
#>  2:  909.6469  918.3054 30-100+  esp2013 1e+05
#>  3:  908.7208  915.1718 30-100+  esp2013 1e+05
#>  4:  902.5358  910.6300 30-100+  esp2013 1e+05
#>  5:  893.1043  909.0838 30-100+  esp2013 1e+05
#>  6: 1509.9958 1516.2048 30-100+  esp2013 1e+05
#>  7: 1500.8061 1517.0321 30-100+  esp2013 1e+05
#>  8: 1490.6371 1512.3974 30-100+  esp2013 1e+05
#>  9: 1479.5267 1496.6290 30-100+  esp2013 1e+05
#> 10: 1465.8329 1492.0557 30-100+  esp2013 1e+05
yll(fc)
#>     stratum    sex  label period     mean   q0.025     q0.1     q0.5     q0.9
#>       <int> <char> <char>  <num>    <num>    <num>    <num>    <num>    <num>
#>  1:       1 female female   2020 966862.8 954425.2 959277.4 967160.4 973790.6
#>  2:       1 female female   2021 959177.0 941057.0 947068.4 959281.2 969947.5
#>  3:       1 female female   2022 952787.7 931264.8 938508.7 953130.2 969262.7
#>  4:       1 female female   2023 946366.4 920196.6 928613.1 945881.9 964682.1
#>  5:       1 female female   2024 939816.6 913177.2 920514.6 940581.4 957675.0
#>  6:       2   male   male   2020 925711.2 912673.8 918036.0 925984.3 933071.8
#>  7:       2   male   male   2021 918801.6 900941.0 908418.1 919001.6 928868.6
#>  8:       2   male   male   2022 910569.9 891219.7 896191.1 910220.8 924206.4
#>  9:       2   male   male   2023 903279.8 878994.6 887066.8 903799.0 919099.6
#> 10:       2   male   male   2024 896435.5 866117.9 874467.2 898638.9 912785.9
#>       q0.975
#>        <num>
#>  1: 979560.9
#>  2: 976338.2
#>  3: 974027.5
#>  4: 970673.7
#>  5: 969531.2
#>  6: 935906.5
#>  7: 936283.8
#>  8: 934167.0
#>  9: 926954.2
#> 10: 924857.6
```
