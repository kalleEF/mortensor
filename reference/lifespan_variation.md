# Lifespan variation from rates or forecast draws

Indices of the variation of the age at death in the period life table,
conditional on survival to age `at`: how unequal lifespans are, as
opposed to how long they are on average (life expectancy). For forecast
draws they are computed path by path, like
[`life_expectancy()`](life_expectancy.md).

## Usage

``` r
lifespan_variation(
  x,
  measure = "edagger",
  age = NULL,
  at = NULL,
  distribution = c("age_at_death", "remaining_life"),
  probs = c(0.025, 0.1, 0.5, 0.9, 0.975),
  paths = FALSE,
  ...
)
```

## Arguments

- x:

  A matrix or array of rates with age in the first dimension, an
  [mt_draws](mt_draws.md) object, or an [mt_data](mt_data.md) object
  (observed rates with Poisson intervals, see
  [`mt_observed()`](mt_observed.md)).

- measure:

  One or more of `"edagger"`, `"sd"`, `"var"`, `"cv"`, `"entropy"`,
  `"gini"`, `"iqr"` (one only for arrays and `paths = TRUE`).

- age:

  An [mt_axis](mt_axis.md) (taken from `x` for draws).

- at:

  Age on which the measures are conditioned (start of an interval;
  default the first age).

- distribution:

  Denominator of the relative measures `"cv"` and `"gini"`: the mean age
  at death (`"age_at_death"`) or the mean remaining life
  (`"remaining_life"`) at age `at`.

- probs:

  Quantiles reported for draws.

- paths:

  For draws, return the values of every path, an array
  `[period, stratum, draw]`.

- ...:

  Passed to the life-table calculation (`ax`, `open`, `k`), see
  [`life_table()`](life_table.md).

## Value

For arrays, an array without the age dimension. For draws and data, a
`data.table` with the stratum columns, `period`, `mean`, quantiles,
`measure` and `at`, computed path by path.

## Details

The ages at death of the life table are taken as `x + a_x` for the
deaths `d_x` of each interval (the discrete formulas of LifeIneq). With
`l` the survivors and `e` the life expectancy at age `at`,
`y_x = x + a_x` and `p_x = d_x / l(at)` for the intervals from `at`:

- `"edagger"`: life disparity, the average remaining life expectancy
  lost at death (Vaupel 1986), `sum(d_x * ebar_x) / l(at)`, where
  `ebar_x` interpolates life expectancy linearly to the average age at
  death in the interval, `e_x + a_x / n_x * (e_(x+n) - e_x)` (the open
  interval uses its `e_x`);

- `"var"`, `"sd"`: variance and standard deviation of the age at death,
  `sum(p_x * (y_x - at - e)^2)`;

- `"cv"`: coefficient of variation, the standard deviation divided by
  the mean age at death `at + e` (or by `e` with
  `distribution = "remaining_life"`);

- `"entropy"`: the life-table entropy of Keyfitz (1977) and Leser
  (1955), `H = e-dagger / e`, the elasticity of life expectancy to a
  proportional change of the rates at all ages (Goldman and Lord 1986);

- `"gini"`: Gini coefficient of the ages at death, the mean absolute
  difference between two ages at death divided by twice the mean age at
  death `at + e` (Hanada 1983; Shkolnikov, Andreev and Begun 2003), or
  by twice the mean remaining life `e` with
  `distribution = "remaining_life"`;

- `"iqr"`: interquartile range, the distance between the ages by which a
  quarter and three quarters of the survivors at age `at` have died. The
  survival curve is interpolated linearly within closed intervals
  (deaths spread evenly, as with the default `a_x = n / 2`) and
  exponentially in the open interval (with the mean `a_x` of that
  interval).

The absolute measures (`edagger`, `sd`, `iqr`) are in years, `var` in
years squared; `cv`, `entropy` and `gini` have no unit.

The discrete formulas treat everyone dying in an interval as dying at
age `x + a_x`, so `var`, `sd`, `cv` and `gini` leave out the spread of
the ages at death within each interval (about `n^2 / 12` years squared
for deaths spread evenly over `n` years). With single years of age to
110+, the standard deviation of French life tables (HMD) is lower by
less than 0.01 years than with ages at death spread evenly within the
intervals, with 5-year groups to 100+ by up to 0.14 years.

The open age interval matters more for every measure. The life table
assumes a constant death rate in it (with the default closure), while
the death rate keeps rising with age, so the ages at death in the open
interval are spread too widely. Compared with single years of age to
100+, an open interval from 85 gave for French life tables of 1950-2006
an e-dagger and an entropy up to 22% higher (2.1 years for women in
2006), an interquartile range and a Gini coefficient up to 13% lower and
a standard deviation up to 4% lower. With an open interval from 90 the
differences were at most 8% (e-dagger, entropy), 5% (interquartile
range), 4% (Gini coefficient) and 1% (standard deviation); from 95 at
most 2%. They grow as more deaths occur at the highest ages. Use single
years of age up to 100 or more where possible, and compare measures only
between tables with the same open age.

The life table needs an open last age interval (otherwise the deaths do
not add up to the survivors at age `at`). The open interval uses the
life-table closure (`open`, see [`life_table()`](life_table.md)); its
average time lived `a_x` is its life expectancy.

## References

Goldman N, Lord G (1986). A new look at entropy and the life table.
*Demography* 23(2): 275-282.

Hanada K (1983). A formula of Gini's concentration ratio and its
application to life tables. *Journal of the Japan Statistical Society,
Japanese Issue* 13(2): 95-98.

Keyfitz N (1977). Mortality comparisons; the male-female ratio. In
*Applied Mathematical Demography*, 54-76. John Wiley and Sons.

Leser CEV (1955). Variations in mortality and life expectation.
*Population Studies* 9(1): 67-71.

Shkolnikov VM, Andreev EM, Begun AZ (2003). Gini coefficient as a life
table function: computation from discrete data, decomposition of
differences and empirical examples. *Demographic Research* 8: 305-358.

van Raalte A, Riffe T, Dudel C (2024). LifeIneq: Calculate indices of
lifespan variation. R package version 00.04.09.
<https://github.com/alysonvanraalte/LifeIneq>

Vaupel JW (1986). How change in age-specific mortality affects life
expectancy. *Population Studies* 40(1): 147-157.

## Examples

``` r
m <- 5e-5 * exp(0.095 * (0:100))
lifespan_variation(m, "edagger", single_ages(0, 100))
#> [1] 10.4904
x <- mt_example(deaths = "single", years = 2015:2019)
lifespan_variation(x, measure = c("edagger", "gini"), at = 65)
#>     stratum    sex  label period       mean     q0.025       q0.1       q0.5
#>       <int> <char> <char>  <num>      <num>      <num>      <num>      <num>
#>  1:       1 female female   2015 8.32653252 8.28650268 8.30044142 8.32662042
#>  2:       1 female female   2016 8.35511514 8.31683126 8.32937049 8.35521479
#>  3:       1 female female   2017 8.35671269 8.31818244 8.33166253 8.35581992
#>  4:       1 female female   2018 8.32388828 8.28524442 8.29841089 8.32333160
#>  5:       1 female female   2019 8.29707412 8.26148309 8.27296038 8.29672086
#>  6:       2   male   male   2015 7.67576125 7.63989496 7.65244014 7.67576397
#>  7:       2   male   male   2016 7.68757690 7.65231498 7.66476279 7.68698523
#>  8:       2   male   male   2017 7.68981676 7.65538380 7.66816994 7.68969267
#>  9:       2   male   male   2018 7.70039788 7.66787043 7.67889076 7.70072088
#> 10:       2   male   male   2019 7.71825775 7.68950978 7.70009451 7.71772973
#> 11:       1 female female   2015 0.06255155 0.06226856 0.06235711 0.06255368
#> 12:       1 female female   2016 0.06260181 0.06231575 0.06243157 0.06260402
#> 13:       1 female female   2017 0.06260496 0.06232645 0.06241969 0.06260100
#> 14:       1 female female   2018 0.06242314 0.06213256 0.06223741 0.06242601
#> 15:       1 female female   2019 0.06221368 0.06194941 0.06203594 0.06221215
#> 16:       2   male   male   2015 0.05960837 0.05933062 0.05944483 0.05961660
#> 17:       2   male   male   2016 0.05965643 0.05939437 0.05949076 0.05965620
#> 18:       2   male   male   2017 0.05967871 0.05944727 0.05952649 0.05968562
#> 19:       2   male   male   2018 0.05973061 0.05949512 0.05957462 0.05973087
#> 20:       2   male   male   2019 0.05989824 0.05967985 0.05976202 0.05989641
#>     stratum    sex  label period       mean     q0.025       q0.1       q0.5
#>       <int> <char> <char>  <num>      <num>      <num>      <num>      <num>
#>           q0.9     q0.975 measure    at
#>          <num>      <num>  <char> <num>
#>  1: 8.35219154 8.36255724 edagger    65
#>  2: 8.37936211 8.39387005 edagger    65
#>  3: 8.38298598 8.39759399 edagger    65
#>  4: 8.34713620 8.36463024 edagger    65
#>  5: 8.32123878 8.33481790 edagger    65
#>  6: 7.69796135 7.70951549 edagger    65
#>  7: 7.71050649 7.71970902 edagger    65
#>  8: 7.70979783 7.72077743 edagger    65
#>  9: 7.72066148 7.73448400 edagger    65
#> 10: 7.73865513 7.75209848 edagger    65
#> 11: 0.06274850 0.06282756    gini    65
#> 12: 0.06278075 0.06288396    gini    65
#> 13: 0.06279280 0.06288415    gini    65
#> 14: 0.06261877 0.06270857    gini    65
#> 15: 0.06239060 0.06248510    gini    65
#> 16: 0.05977399 0.05984672    gini    65
#> 17: 0.05981362 0.05992077    gini    65
#> 18: 0.05982297 0.05990479    gini    65
#> 19: 0.05987881 0.05997550    gini    65
#> 20: 0.06003873 0.06012775    gini    65
#>           q0.9     q0.975 measure    at
#>          <num>      <num>  <char> <num>
fc <- mt_forecast(mt_fit(lee_carter(), x), h = 5, n_draws = 100, seed = 1)
lifespan_variation(fc, measure = "sd", at = 30)
#>     stratum    sex  label period     mean   q0.025     q0.1     q0.5     q0.9
#>       <int> <char> <char>  <num>    <num>    <num>    <num>    <num>    <num>
#>  1:       1 female female   2020 12.88953 12.86809 12.87580 12.88836 12.90281
#>  2:       1 female female   2021 12.87254 12.84479 12.85170 12.87400 12.89205
#>  3:       1 female female   2022 12.85895 12.82497 12.83542 12.85701 12.88240
#>  4:       1 female female   2023 12.84689 12.81108 12.81878 12.84405 12.87891
#>  5:       1 female female   2024 12.83593 12.80117 12.81222 12.83368 12.86471
#>  6:       2   male   male   2020 12.67651 12.67264 12.67470 12.67650 12.67857
#>  7:       2   male   male   2021 12.66743 12.66247 12.66411 12.66727 12.67026
#>  8:       2   male   male   2022 12.65867 12.65206 12.65519 12.65847 12.66264
#>  9:       2   male   male   2023 12.65073 12.64388 12.64643 12.65076 12.65515
#> 10:       2   male   male   2024 12.64346 12.63536 12.63832 12.64351 12.64795
#>       q0.975 measure    at
#>        <num>  <char> <num>
#>  1: 12.91317      sd    30
#>  2: 12.90314      sd    30
#>  3: 12.90010      sd    30
#>  4: 12.88969      sd    30
#>  5: 12.88436      sd    30
#>  6: 12.67976      sd    30
#>  7: 12.67294      sd    30
#>  8: 12.66529      sd    30
#>  9: 12.65812      sd    30
#> 10: 12.65156      sd    30
```
