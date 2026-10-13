# Disaggregate forecast draws to a finer age grid

Maps draws of age groups to finer ages (for example single years), path
by path, when a model was forecast at the resolution of the observed
deaths. Two approaches:

## Usage

``` r
mt_disaggregate(
  draws,
  age,
  exposure,
  method = c("hyman_kannisto", "natural_gompertz", "shape"),
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

  `"shape"`, `"hyman_kannisto"` or `"natural_gompertz"`.

- shape_rates:

  For `"shape"`: rates on the finer ages, `[age, stratum]` or
  `[age, period, stratum]` for the forecast periods (see
  [`mt_shape_rates()`](mt_shape_rates.md)).

- omega:

  Maximum age.

## Value

An `mt_draws` object on `age`.

## Details

- `"shape"`: age-specific factors relative to the group rate, estimated
  from fine-resolution rates (`shape_rates`), exposure-weighted to 1
  within each group. [`mt_shape_rates()`](mt_shape_rates.md) estimates
  them from training data at the finer resolution, constant (mean of the
  last years) or changing over the forecast years (an age-period GAM
  extrapolated in time);

- `"hyman_kannisto"` or `"natural_gompertz"`: a spline through the group
  rates of each draw (monotone Hyman spline with a Kannisto tail, or a
  natural spline with a log-linear tail), corrected so that the
  exposure-weighted mean of the fine rates equals each group rate.

## Examples

``` r
x <- mt_example(years = 2000:2019, deaths = "single")
xg <- mt_regroup(x, age = c(seq(30, 85, 5), 90))
fc <- mt_forecast(mt_fit(lee_carter(), xg), h = 3, n_draws = 20, seed = 1)
e <- mt_array(x, "exposure")[, 20, ]
fs <- mt_disaggregate(fc, mt_grid(x)$age, exposure = e)
dim(fs$rates)
#> [1] 71  3  2 20
```
