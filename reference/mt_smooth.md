# Smooth death rates on the grid

Smooths death rates over age, over time or over both, with the exposures
as weights or offset. Deaths are replaced by exposure times the smoothed
rate; exposures are unchanged. Cells with missing deaths do not enter
the fit and stay missing (only the 2D GAM predicts them). Deaths and
exposures must be available on the grid (use
[`mt_ungroup()`](mt_ungroup.md) first for grouped data). In a
forecasting pipeline, use it as a preprocessing step
(`mt_prep(step_smooth(...))`): models are then fitted to the smoothed
rates, while the observation noise of their predictive draws, residuals
and backtest scales come from the rates before smoothing ("smooth input,
raw residuals").

## Usage

``` r
mt_smooth(
  data,
  dims = 1,
  lambda = NULL,
  criterion = c("bic", "aic"),
  kr = 2,
  deg = 3,
  pord = 2,
  lambda_range = c(0.01, 1e+08),
  method = c("pspline", "whittaker", "gam"),
  ages = NULL,
  along = c("age", "period"),
  k = NULL
)
```

## Arguments

- data:

  An [mt_data](mt_data.md) object.

- dims:

  1 (over one direction, see `along`) or 2 (age-period surface;
  `"pspline"` or `"gam"`).

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

- method:

  `"pspline"` (Poisson P-splines, default), `"whittaker"` or `"gam"`
  (package mgcv); see Details. With `dims = 1` and `along = "period"`
  the default is `"whittaker"`.

- ages:

  Range of ages that is smoothed, `c(lower, upper)`; other ages keep
  their observed deaths and do not enter the fit. By default all ages
  except the first year of life when the grid starts with age 0 as its
  own interval (infant mortality does not follow the age pattern of the
  following ages).

- along:

  For `dims = 1`: `"age"` or `"period"` (`"whittaker"` only).

- k:

  Basis dimensions of the tensor-product smooth for `dims = 2`,
  `method = "gam"`: ages and years (default: up to 10 each).

## Value

An `mt_data` object with smoothed deaths; diagnostics as attribute
`"smooth"`.

## Details

Methods and directions:

- `dims = 1, along = "age"` (default): each year separately over age,
  with Poisson P-splines (`"pspline"`), Whittaker smoothing of log rates
  with the deaths as weights (`"whittaker"`), or a Gaussian additive
  model of log rates with the deaths as weights (`"gam"`);

- `dims = 1, along = "period"`: each age separately over time, Whittaker
  smoothing of log rates with the deaths as weights and one smoothing
  parameter per stratum (generalised cross-validation pooled over ages).
  The weights are scaled within each age, so the parameter gives the
  same smoothness at every age; ages with many deaths are not left
  unsmoothed;

- `dims = 2`: the age-period surface, with a Poisson P-spline surface
  (`"pspline"`; Currie, Durban and Eilers 2004) or a quasi-Poisson
  additive model with a tensor-product P-spline `te(age, year)` and the
  exposures as offset (`"gam"`; Wood 2017), smoothing parameters by REML
  and overdispersion estimated, so that year-to-year variation beyond
  Poisson noise is not taken for signal.

Smoothing over time removes year-to-year fluctuations from the fitted
trends and from the starting level of forecasts; the last years are
estimated from one side only. It also removes the year-to-year variation
from which models such as Lee-Carter estimate the uncertainty of their
period index, so their intervals become too narrow unless the spread is
re-estimated against the raw observations
([`recalibrated()`](recalibrated.md)). Smoothing over age alone changes
little for large populations and helps for sparse data such as
cause-specific rates of small regions. Whether smoothing improves a
forecast is a question for a backtest; the vignette "Smoothing death
rates" shows one. Use [`mt_smooth_check()`](mt_smooth_check.md) to see
how much was removed and [`mt_plot_smooth()`](mt_smooth_check.md) to
look at it. The 2D P-spline surface is slow for many ages and years
(minutes per stratum); the 2D GAM takes about a second.

## References

Currie ID, Durban M, Eilers PHC (2004). Smoothing and forecasting
mortality rates. *Statistical Modelling* 4(4): 279-298.

Wood SN (2017). *Generalized Additive Models: An Introduction with R*
(2nd ed.). Chapman and Hall/CRC.

## Examples

``` r
x <- mt_example(years = 2001:2003, deaths = "single")
s <- mt_smooth(x)
head(attr(s, "smooth"))
#>    measure stratum period     method lambda_age lambda_period       ed
#>     <char>   <int> <char>     <char>      <num>         <num>    <num>
#> 1:  deaths       1   2001 pspline_1d   84044.06            NA 5.767171
#> 2:  deaths       1   2002 pspline_1d   55318.96            NA 6.280081
#> 3:  deaths       1   2003 pspline_1d   27572.52            NA 7.243212
#> 4:  deaths       2   2001 pspline_1d  101398.86            NA 5.533108
#> 5:  deaths       2   2002 pspline_1d   23289.26            NA 7.546782
#> 6:  deaths       2   2003 pspline_1d   92069.24            NA 5.638856
#>    max_rel_diff at_bound
#>           <num>   <lgcl>
#> 1:    0.1988679    FALSE
#> 2:    0.2322712    FALSE
#> 3:    0.2057177    FALSE
#> 4:    0.2271092    FALSE
#> 5:    0.2455273    FALSE
#> 6:    0.2102315    FALSE
xt <- mt_example(deaths = "single", years = 2000:2015)
st <- mt_smooth(xt, along = "period", method = "whittaker")
mt_smooth_check(xt, st)
#>    stratum    sex  label    mad_log dispersion  trend_ratio correlation
#>      <int> <char> <char>      <num>      <num>        <num>       <num>
#> 1:       1 female female 0.02113118  0.8846054 0.0013826242   0.9994987
#> 2:       2   male   male 0.02072276  0.8750107 0.0009253008   0.9996417
```
