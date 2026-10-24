# Pessimism of forecasts by cause of death

Compares a forecast of all-cause mortality with the sum of separate
forecasts of the causes of death, as Wilmoth (1995) does. Wilmoth calls
the first the "total" forecast and the second the "all-cause" forecast;
here they are `"total"` and `"causes"`. The sum of cause forecasts tends
to be higher, because it is dominated over time by the causes that fall
slowest or rise. Every measure is computed path by path:

## Usage

``` r
mt_pessimism(
  x,
  causes = NULL,
  var = "cause",
  total_level = "all",
  by = c("total", "age", "cause"),
  at = NULL,
  periods = NULL,
  h = 10,
  to = NULL,
  pairing = c("independent", "draw"),
  n_draws = 1,
  probs = c(0.025, 0.5, 0.975),
  paths = FALSE,
  seed = 1,
  ...
)
```

## Arguments

- x:

  An [mt_draws](mt_draws.md) object (forecast) or an
  [mt_data](mt_data.md) object (observed deaths by cause, for the
  approximation from two periods).

- causes:

  Optional [mt_draws](mt_draws.md) object with the forecasts of the
  causes (see Forecasts).

- var:

  Stratum variable with the causes.

- total_level:

  Level of `var` for all causes.

- by:

  `"total"` (life expectancy), `"age"` (rates and m-ratios) or `"cause"`
  (contributions of the causes to the m-ratios).

- at:

  Age at which life expectancy is computed (start of an age interval;
  default the first age).

- periods:

  For data: a list of two vectors of observed years, the first and the
  second period (default: the two halves of the observed years).

- h, to:

  For data: the number of forecast years after the last observed year,
  or the last forecast year.

- pairing:

  `"independent"` or `"draw"`: how the draws of `x` and `causes` are
  paired.

- n_draws:

  For data: 1 (the observed deaths) or the number of Poisson draws.

- probs:

  Quantile levels.

- paths:

  Return the value of every path (columns `draw` and `value`) instead of
  the mean and quantiles.

- seed:

  Seed of the shuffling of `causes` and of Poisson draws.

- ...:

  Passed to the life table (`ax`, `open`, `k`, `omega`, see
  [`life_table()`](life_table.md)) for life expectancy.

## Value

A `data.table` with the stratum variables other than `var`, `period`,
for `by = "age"` and `"cause"` also `age_lo` and `age`, for `"cause"`
the column `var` with the cause, then `measure`, for `"total"` the age
`at`, and the `mean` and quantiles over paths (`q0.025`, ...). With
`paths = TRUE`, `draw` and `value` replace the mean and quantiles.

## Details

- `by = "total"`: life expectancy at `at` of the total forecast
  (`measure = "total"`) and of the sum of the causes (`"causes"`), and
  their `"difference"`, total minus causes, the difference
  `e0(Tot) - e0(All)` that Wilmoth decomposes (Section 3, equation 15
  and Table 5). It is positive where the forecast by cause is more
  pessimistic.

- `by = "age"`: the rates of the total forecast (`"total"`) and of the
  sum of the causes (`"causes"`), and their ratio, the m-ratio of
  Wilmoth's equation (6), `M(x, t) = sum_i m_i(x, t) / m(x, t)`
  (`"m_ratio"`), his measure of the pessimism of forecasts by cause. It
  is above 1 where the forecast by cause is more pessimistic.

- `by = "cause"`: the contribution of each cause to the m-ratio minus 1
  (equation 13), `m_i(x, t) / m(x, t) - p_i(x)`, where `p_i(x)` is the
  proportion of the deaths at age `x` that are due to cause `i` in the
  last observed year (equation 5). The contributions of an age add up to
  its m-ratio minus 1 (Table 3).

## Forecasts

For [mt_draws](mt_draws.md), the forecast of all causes comes from `x`
and the forecasts of the causes from `causes`. Without `causes`, `x`
holds both: the strata of `var` at `total_level` and those of the
causes, for example the forecasts of [`mt_fit()`](mt_fit.md) with a
separate model for each stratum. With `causes`, for example the draws of
a model of causes such as [`coda_causes()`](coda_causes.md), `x` holds
the total forecast: its strata at `total_level`, or all its strata if
`var` is not a stratum variable of `x`. The causes are the levels of
`var` other than `total_level` and the groups of a hierarchy stored in
the object (see [`mt_aggregate_strata()`](mt_aggregate_strata.md)), and
the strata of the total and of the causes are matched on the other
stratum variables (for example sex and deprivation quintile). The causes
should cover all deaths: a warning reports observed causes that do not
add up to the observed total in the last training year.

Wilmoth derives equations (6) and (13) for proportional rate of change
models, which project the observed rates of the last year `t0` with
constant annual factors. Equation (6) is defined for any forecast; for
other models, and for models whose forecasts start from fitted rates
(his footnote 2), the m-ratio is still the ratio of the forecasts and
the contributions still add up to the m-ratio minus 1, but the m-ratio
of the first forecast years also contains the jump-off of the models
(see [`mt_jump_off()`](mt_jump_off.md)).

Draws of two objects come from different models and are independent: by
default (`pairing = "independent"`) the draws of `causes` are shuffled
before they are paired with those of `x`, as in
[`mt_decompose()`](mt_decompose.md). Use `pairing = "draw"` when
`causes` is derived from the draws of `x`.

## Approximation from two periods

For an [mt_data](mt_data.md) object, `mt_pessimism()` computes Wilmoth's
approximation from the distributions of deaths by cause in two periods
(Section 3, equations 12 and 14). The rates of each cause and period are
its deaths divided by the exposures, both summed over the years of the
period (Appendix C), and the total rate is the sum of the causes (a
warning reports causes that do not add up to the level `total_level` of
the data, where it exists). Both are projected from the midpoint `t2` of
the second period at the geometric mean rate of change between the
midpoints `t1` and `t2`, the "first-to-last" method of his equations (8)
to (11) applied to periods:
`m(x, t) = m2(x) (m2(x) / m1(x))^((t - t2) / (t2 - t1))`. The m-ratio of
these forecasts is equation (12),
`sum_i p_i (p_i / q_i)^((t - t2) / (t2 - t1))`, with `q_i` and `p_i` the
proportions of deaths by cause in the first and the second period, and
the contributions are equation (14). Wilmoth shows (Appendix B) that
this m-ratio exceeds 1 at every horizon unless the two distributions are
equal. The difference of life expectancies is the quantity that his
Table 5 decomposes by age and cause.

The default `periods` are the first and the second half of the observed
years, the "first-to-second" method that Wilmoth uses for his
decompositions; with an odd number of years the middle year is left out,
so that the halves have the same length (package-specific, reference
unknown). His "first-to-last" method uses the first and the last three
years, for example `periods = list(1951:1953, 1988:1990)`. The forecast
years are those after the last observed year (`h`, `to`).

With `n_draws = 1` (default) the observed deaths are used. With more
draws the deaths are Poisson draws, as in
[`mt_observed()`](mt_observed.md), and the intervals show the sampling
variation of the approximation (package-specific, reference unknown). A
cause without deaths at an age in both periods contributes nothing there
(the dashes of Wilmoth's Table 3). A cause without deaths in the first
period and with deaths in the second would grow without bound: its
deaths in the first period are set to 0.5, as [`rwd()`](rwd.md) does
with zero counts, with a warning (package-specific, reference unknown).

## References

Wilmoth JR (1995). Are mortality projections always more pessimistic
when disaggregated by cause of death? *Mathematical Population Studies*
5(4): 293-319.
[doi:10.1080/08898489509525409](https://doi.org/10.1080/08898489509525409)

## See also

[`mt_coherence()`](mt_coherence.md) for the ratio of summed strata to
aggregates in any hierarchy, [`mt_decompose()`](mt_decompose.md) for
decompositions of life expectancy by age and cause,
[`reconciled()`](reconciled.md) and [`coda_causes()`](coda_causes.md)
for forecasts by cause that add up.

## Examples

``` r
x <- mt_example_causes(years = 2000:2019)
# separate Lee-Carter models for all causes and for each cause
fc <- mt_forecast(mt_fit(lee_carter(), x), h = 10, n_draws = 50, seed = 1)
mt_pessimism(fc, at = 65)[period == 2029]
#>       sex period    measure    at       mean     q0.025       q0.5     q0.975
#>    <char>  <num>     <char> <num>      <num>      <num>      <num>      <num>
#> 1: female   2029      total    65 23.7035780 23.3216321 23.7116535 24.0669705
#> 2: female   2029     causes    65 23.4944717 23.2815900 23.5195917 23.6932909
#> 3: female   2029 difference    65  0.2091063 -0.2813627  0.2067860  0.6321682
#> 4:   male   2029      total    65 19.5776676 19.1667844 19.6262465 19.9755737
#> 5:   male   2029     causes    65 19.4764975 19.3181567 19.4929899 19.6335253
#> 6:   male   2029 difference    65  0.1011701 -0.3036909  0.1087295  0.6226814
m <- mt_pessimism(fc, by = "age")
m[measure == "m_ratio" & period == 2029 & age %in% c("50", "70", "90")]
#>       sex period age_lo    age measure     mean    q0.025     q0.5   q0.975
#>    <char>  <num>  <num> <char>  <char>    <num>     <num>    <num>    <num>
#> 1: female   2029     50     50 m_ratio 1.009184 0.9269285 1.005622 1.102265
#> 2: female   2029     70     70 m_ratio 1.019328 0.9395110 1.016333 1.115691
#> 3: female   2029     90     90 m_ratio 1.028305 0.9614267 1.028123 1.085640
#> 4:   male   2029     50     50 m_ratio 1.026957 0.9366452 1.035116 1.099733
#> 5:   male   2029     70     70 m_ratio 1.017059 0.9580756 1.016272 1.096000
#> 6:   male   2029     90     90 m_ratio 1.010961 0.9581379 1.011347 1.074437
# Wilmoth's approximation from the halves 2000-2009 and 2010-2019
a <- mt_pessimism(x, by = "cause", to = 2039)
a[period == 2039 & age == "70"]
#>       sex period age_lo    age  cause      measure        mean      q0.025
#>    <char>  <num>  <num> <char> <char>       <char>       <num>       <num>
#> 1: female   2039     70     70    cvd contribution -0.10070629 -0.10070629
#> 2: female   2039     70     70 cancer contribution  0.12682346  0.12682346
#> 3: female   2039     70     70  other contribution  0.04443719  0.04443719
#> 4:   male   2039     70     70    cvd contribution -0.10560148 -0.10560148
#> 5:   male   2039     70     70 cancer contribution  0.11115490  0.11115490
#> 6:   male   2039     70     70  other contribution  0.06836902  0.06836902
#>           q0.5      q0.975
#>          <num>       <num>
#> 1: -0.10070629 -0.10070629
#> 2:  0.12682346  0.12682346
#> 3:  0.04443719  0.04443719
#> 4: -0.10560148 -0.10560148
#> 5:  0.11115490  0.11115490
#> 6:  0.06836902  0.06836902
```
