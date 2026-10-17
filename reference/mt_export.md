# Export forecasts and result tables

Writes a tidy CSV of forecast draws (their summary, `what = "summary"`,
or all draws in long format, `what = "draws"`) or of any result table,
such as those of [`life_expectancy()`](life_expectancy.md),
[`asr()`](derived.md), [`mt_contrast()`](mt_contrast.md),
[`mt_inequality()`](mt_inequality.md) or [`mt_score()`](mt_score.md).
The file name of a forecast records the model, the grid and the forecast
origin unless `file` is given.

## Usage

``` r
mt_export(
  x,
  dir,
  what = c("summary", "draws"),
  file = NULL,
  probs = c(0.025, 0.1, 0.5, 0.9, 0.975),
  digits = NULL
)
```

## Arguments

- x:

  An [mt_draws](mt_draws.md) object, an
  [`mt_decompose()`](mt_decompose.md) result (its summary is written) or
  a data frame.

- dir:

  Output directory.

- what:

  For draws: `"summary"` or `"draws"`.

- file:

  Optional file name (required for tables other than draws unless the
  default `"mortensor_table.csv"` is wanted).

- probs:

  Quantiles for the summary of draws.

- digits:

  Optional number of significant digits of numeric columns.

## Value

The path of the written file, invisibly.

## Examples

``` r
x <- mt_example(deaths = "single", years = 2005:2019)
fc <- mt_forecast(mt_fit(rwd(), x), h = 2, n_draws = 10, seed = 1)
f <- mt_export(fc, tempdir())
head(read.csv(f), 3)
#>   stratum    sex  label age_lo age period         mean       q0.025
#> 1       1 female female     30  30   2020 0.0004868688 0.0003333668
#> 2       1 female female     31  31   2020 0.0006259915 0.0004828314
#> 3       1 female female     32  32   2020 0.0006665331 0.0005311451
#>           q0.1         q0.5         q0.9       q0.975
#> 1 0.0003674486 0.0004761180 0.0006350217 0.0006703991
#> 2 0.0004908948 0.0005503622 0.0008643498 0.0008654412
#> 3 0.0005358900 0.0006447373 0.0008039168 0.0008249697
f2 <- mt_export(life_expectancy(fc, at = 65), tempdir(), file = "e65.csv", digits = 4)
```
