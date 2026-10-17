# Estimate grid-resolution data from grouped observations

Fills the modelling grid (single years of age by default) from deaths
and exposures observed in coarser age groups, using the penalized
composite link model (PCLM; Eilers 2007, Rizzi et al. 2015).

## Usage

``` r
mt_ungroup(
  data,
  age = NULL,
  dims = 1,
  method = c("pclm", "monotone", "hyman_kannisto", "natural_gompertz", "gam"),
  exposure_method = c("pclm", "monotone"),
  lambda = NULL,
  criterion = c("bic", "aic"),
  kr = 2,
  deg = 3,
  pord = 2,
  lambda_range = c(0.01, 1e+08),
  preserve_totals = TRUE,
  tail = c("kannisto", "gompertz"),
  tail_width = 2
)
```

## Arguments

- data:

  An [mt_data](mt_data.md) object.

- age:

  Target age axis ([mt_axis](mt_axis.md)). By default the grid when it
  was set (with [`mt_set_grid()`](mt_set_grid.md) or
  [`step_grid()`](mt_prep.md)) and is finer than the data, otherwise
  single years of age from the lowest age to the start of the open
  interval.

- dims:

  1 or 2.

- method:

  Method for deaths: `"pclm"` (default); `"monotone"` (cumulative deaths
  interpolated by a monotone Hyman spline, counts only, 1D, ignores
  exposures); or rate-based spline methods, which need exposures on the
  grid: `"hyman_kannisto"` (monotone spline of log group rates with a
  Kannisto tail for the open interval), `"natural_gompertz"` (natural
  spline with a log-linear tail) and `"gam"` (2D tensor-product GAM of
  log group rates over age and year with a Kannisto tail, package mgcv).
  All three preserve the exposure-weighted group rates exactly. Median
  absolute log errors against single-year rates (England & Wales men
  1980-2011, 5-year groups to 90-100; and the simulated
  [`mt_example()`](mt_example.md) with deaths to 90+): PCLM 0.015/0.024
  and 0.012/0.006 (ages 30-89/90+), `"hyman_kannisto"` 0.015/0.029 and
  0.013/0.028, `"natural_gompertz"` 0.015/0.027 and 0.013/0.030, `"gam"`
  0.015/0.031 and 0.012/0.031. PCLM is the most accurate in the open
  interval; the spline methods are 5 to 10 times faster.

- exposure_method:

  `"pclm"` or `"monotone"`.

- lambda:

  Fixed smoothing parameter(s); `NULL` selects them by `criterion`.

- criterion:

  `"bic"` or `"aic"`.

- kr:

  Knot spacing in years (age and period).

- deg, pord:

  Degree of the B-splines and order of the difference penalty.

- lambda_range:

  Search range for the smoothing parameters.

- preserve_totals:

  Rescale estimates to observed group totals.

- tail, tail_width:

  For `method = "gam"`: the parametric curve that takes over from the
  GAM above the start of the last age group (`"kannisto"`, logistic, or
  `"gompertz"`, log-linear), fitted to the GAM rates from age 70, and
  the width (years) of the logistic blend between the two.

## Value

An `mt_data` object with deaths (and exposures) on the grid. Diagnostics
per fit are attached as attribute `"ungroup"`: smoothing parameters,
effective dimension, the largest relative difference between fitted and
observed group totals before rescaling, and a flag when a smoothing
parameter is at the boundary of the search range.

## Details

- **Exposures** observed coarser than the grid are ungrouped first, as
  counts (`exposure_method = "pclm"`) or by a monotone spline of the
  cumulative counts (`"monotone"`, Hyman filter), which preserves group
  totals exactly.

- **Deaths** are ungrouped as **rates** when exposures are available on
  the grid: the exposures enter as offset, so cohort-size waves stay in
  the exposures and the estimated rate schedule is smooth. Without
  exposures, deaths are ungrouped as counts per year of age (1D only; 2D
  smoothing of death counts cannot represent cohort waves and is
  refused).

- The open age interval is expanded into single years up to `omega`; for
  rates the exposure of the open interval is spread with
  stationary-population weights implied by the fitted rates.

- `dims = 1` fits each stratum and period separately; `dims = 2` fits a
  smooth age-period surface per stratum (two smoothing parameters).

With `preserve_totals = TRUE`, the estimated deaths within each observed
group are rescaled to the observed total, so aggregation reproduces the
data exactly. Cells already observed at grid resolution are kept.

In backtests use [`step_ungroup()`](mt_prep.md) so that ungrouping, in
particular in 2D, sees only the training window.

## References

Eilers PHC (2007). Ill-posed problems with counts, the composite link
model and penalized likelihood. *Statistical Modelling* 7(3): 239-254.

Rizzi S, Gampe J, Eilers PHC (2015). Efficient estimation of smooth
distributions from coarsely grouped data. *American Journal of
Epidemiology* 182(2): 138-147.
[doi:10.1093/aje/kwv020](https://doi.org/10.1093/aje/kwv020)

## Examples

``` r
x <- mt_example(years = 2010:2014)
u <- mt_ungroup(x)
u
#> <mt_data> 2 strata by sex, omega = 110
#>   deaths     single years 30-99, 100+; periods 2010-2014 (width 1); 710 cells
#>   exposure   single years 30-99, 100+; periods 2010-2014 (width 1); 710 cells
#>   grid       <mt_axis: age, 71 intervals> 30 31 32 33 ... 98 99 100+
#>              <mt_axis: period, 5 intervals> 2010 2011 2012 2013 2014
head(attr(u, "ungroup"))
#>    measure stratum period       method lambda_age lambda_period       ed
#>     <char>   <int> <char>       <char>      <num>         <num>    <num>
#> 1:  deaths       1   2010 pclm_rate_1d  19515.355            NA 6.448289
#> 2:  deaths       1   2011 pclm_rate_1d  27584.245            NA 6.086151
#> 3:  deaths       1   2012 pclm_rate_1d   7094.886            NA 7.545750
#> 4:  deaths       1   2013 pclm_rate_1d   2340.519            NA 8.833281
#> 5:  deaths       1   2014 pclm_rate_1d   1880.462            NA 9.096312
#> 6:  deaths       2   2010 pclm_rate_1d  18039.474            NA 6.841390
#>    max_rel_diff at_bound
#>           <num>   <lgcl>
#> 1:   0.05824939    FALSE
#> 2:   0.06301434    FALSE
#> 3:   0.04597185    FALSE
#> 4:   0.06615491    FALSE
#> 5:   0.02673557    FALSE
#> 6:   0.04517797    FALSE
```
