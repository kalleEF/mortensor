# Disaggregate forecast draws to a finer age grid

Maps draws of age groups to finer ages (for example single years), path
by path, when a model was forecast at the resolution of the observed
deaths.

## Usage

``` r
mt_disaggregate(
  draws,
  age,
  exposure,
  method = "shape",
  shape_rates = NULL,
  omega = 110
)
```

## Arguments

- draws:

  An [mt_draws](mt_draws.md) object on age groups.

- age:

  The finer [mt_axis](mt_axis.md); its cells must nest in the groups.

- exposure:

  Exposures on the finer ages, `[age, stratum]` or
  `[age, period, stratum]`, used as weights.

- method:

  `"shape"` (proportional allocation by a standard schedule).

- shape_rates:

  Rates on the finer ages, `[age, stratum]` or `[age, period, stratum]`
  for the forecast periods (see
  [`mt_shape_rates()`](mt_shape_rates.md)).

- omega:

  Maximum age.

## Value

An `mt_draws` object on `age`.

## Details

The method (`"shape"`) is proportional allocation by a standard
schedule. The rate of each finer age is the group rate of the draw times
a fixed age-specific factor taken from rates at the finer resolution
(`shape_rates`), the factors being scaled to an exposure-weighted mean
of 1 within each group, so that the exposure-weighted mean of the finer
rates equals the group rate of every draw.
[`mt_shape_rates()`](mt_shape_rates.md) estimates the shape from the
last observed years of data at the finer resolution; with grouped
training data, estimate those first with [`mt_ungroup()`](mt_ungroup.md)
(for example `mt_shape_rates(mt_ungroup(x), periods)`). A group whose
shape rates are all 0 gets a flat shape. Missing group rates (for
example of a stratum whose fit failed) give missing rates at the finer
ages of that group.

## Examples

``` r
x <- mt_example(years = 2000:2019, deaths = "single")
xg <- mt_regroup(x, age = c(seq(30, 85, 5), 90))
fc <- mt_forecast(mt_fit(lee_carter(), xg), h = 3, n_draws = 20, seed = 1)
e <- mt_array(x, "exposure")[, 20, ]
sh <- mt_shape_rates(x, periods = 2020:2022)
fs <- mt_disaggregate(fc, mt_grid(x)$age, exposure = e, shape_rates = sh)
dim(fs$rates)
#> [1] 71  3  2 20
```
