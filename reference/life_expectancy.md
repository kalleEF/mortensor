# Life expectancy from rates or forecast draws

Life expectancy of the period life table of
[`life_table()`](life_table.md) (Preston, Heuveline and Guillot 2001),
for draws path by path.

## Usage

``` r
life_expectancy(
  x,
  age = NULL,
  at = NULL,
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

- age:

  An [mt_axis](mt_axis.md) (taken from `x` for draws).

- at:

  Age at which life expectancy is computed (start of an interval;
  default the first age).

- probs:

  Quantiles reported for draws.

- paths:

  For draws, return the values of every path, an array
  `[period, stratum, draw]`.

- ...:

  Passed to the life-table calculation (`ax`, `open`, `k`, `omega`). A
  rate of 0 in the open interval is replaced as in
  [`life_table()`](life_table.md).

## Value

For arrays, an array without the age dimension. For draws, a
`data.table` with the mean and quantiles per stratum and period,
computed path by path.

## References

Preston SH, Heuveline P, Guillot M (2001). *Demography: Measuring and
Modeling Population Processes*. Blackwell.

## Examples

``` r
x <- mt_example(deaths = "single")
m <- mt_array(x, "deaths") / mt_array(x, "exposure")
life_expectancy(m, mt_grid(x)$age, at = 65)[1:3, ]
#>       stratum
#> period   female     male
#>   1990 21.17552 16.96333
#>   1991 21.28877 16.89039
#>   1992 21.36637 16.94218
```
