# Smooth death rates on the grid

Smooths death rates over age, over time or over both, with the exposures
as weights or offset. Deaths are replaced by exposure times the smoothed
rate; exposures are unchanged. Cells with missing deaths do not enter
the fit and stay missing (only the 2D GAM predicts them). Deaths and
exposures must be available on the grid (use
[`mt_ungroup()`](mt_ungroup.md) first for grouped data). In a
forecasting pipeline, use it as a preprocessing step
(`mt_prep(step_smooth(...))`): models are then fitted to the smoothed
rates, while their residuals, observed rates, backtest scales and
jump-off offsets come from the deaths before smoothing ("smooth input,
raw residuals"). Predictive draws add Poisson deaths given the model's
rates, as without smoothing.

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

  Fixed smoothing parameter(s) of the P-splines or the Whittaker
  smoother; `NULL` selects them.

- criterion:

  `"bic"` or `"aic"`.

- kr:

  Knot spacing in years (age and period).

- deg, pord:

  Degree of the B-splines and order of the difference penalty.

- lambda_range:

  Search range of the smoothing parameters of the P-splines and the
  Whittaker smoother (in both directions).

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
`"smooth"`, one row per fit, with the smoothing parameter, the effective
degrees of freedom (`edf`, or `ed` for P-splines) and, for the Whittaker
smoother, `at_bound` when the selected parameter is at an end of
`lambda_range` before the fit reached its limit there (interpolation at
the lower end, a polynomial of degree `pord - 1` at the upper end); such
fits are also reported with a warning.

## Details

Methods and directions:

- `dims = 1, along = "age"` (default): each year separately over age,
  with Poisson P-splines and the exposures as offset (`"pspline"`; the
  one-dimensional case of Currie, Durban and Eilers 2004), the weighted
  Whittaker smoother of log rates with the deaths as weights
  (`"whittaker"`; Whittaker 1923; Eilers 2003), or a quasi-Poisson
  additive model `s(age)` with the exposures as offset (`"gam"`; Wood
  2017), smoothing parameter by REML and overdispersion estimated;

- `dims = 1, along = "period"`: each age separately over time, the
  weighted Whittaker smoother of log rates with the deaths as weights
  (the approximate inverse variances of the log rates);

- `dims = 2`: the age-period surface, with a Poisson P-spline surface
  (`"pspline"`; Currie, Durban and Eilers 2004) or a quasi-Poisson
  additive model with a tensor-product P-spline `te(age, year)` and the
  exposures as offset (`"gam"`; Wood 2017), smoothing parameters by REML
  and overdispersion estimated, so that year-to-year variation beyond
  Poisson noise is not taken for signal.

With `"pspline"`, an open age interval of the grid is expanded into
single years up to `omega` and its exposure is split by the PCLM of the
exposures, as in [`mt_ungroup()`](mt_ungroup.md) (`open_split = "pclm"`,
Rizzi, Gampe and Eilers 2015).

The Whittaker smoother minimises
`sum(w (y - z)^2) + lambda sum((D z)^2)` for log rates `y`, weights `w`
and differences `D` of order `pord` (divided differences over age, so
that groups of different widths are smoothed on the age scale). Its
smoothing parameter is chosen for every fit separately, by generalised
cross-validation within `lambda_range`: for each year (over age) or for
each age (over time) and stratum. Over time, with deaths as weights,
ages whose year-to-year variation is large relative to their Poisson
noise (typically ages with many deaths) are smoothed little or not at
all, and ages with few deaths may be smoothed to a straight line of the
log rates. The P-spline fits choose their smoothing parameters by
`criterion` within `lambda_range`.

Where a method needs one age per cell (the Whittaker smoother over age
and the GAMs), an open last interval such as 90+ is placed at its lower
bound plus half the width of the interval below it (92.5 after 85-89), a
convention of this package. The P-splines instead expand the open
interval into single years up to `omega`, as
[`mt_ungroup()`](mt_ungroup.md) does.

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

Eilers PHC (2003). A perfect smoother. *Analytical Chemistry* 75(14):
3631-3636.

Whittaker ET (1923). On a new method of graduation. *Proceedings of the
Edinburgh Mathematical Society* 41: 63-75.

Wood SN (2017). *Generalized Additive Models: An Introduction with R*
(2nd ed.). CRC Press.

## Examples

``` r
x <- mt_example(years = 2001:2003, deaths = "single")
s <- mt_smooth(x)
head(attr(s, "smooth"))
#>    measure stratum period     method lambda_age lambda_period       ed
#>     <char>   <int> <char>     <char>      <num>         <num>    <num>
#> 1:  deaths       1   2001 pspline_1d   61501.39            NA 6.073523
#> 2:  deaths       1   2002 pspline_1d   29356.92            NA 7.071097
#> 3:  deaths       1   2003 pspline_1d   23158.84            NA 7.425411
#> 4:  deaths       2   2001 pspline_1d   91819.37            NA 5.632183
#> 5:  deaths       2   2002 pspline_1d   20711.89            NA 7.711461
#> 6:  deaths       2   2003 pspline_1d   76056.64            NA 5.855391
#>    max_rel_diff at_bound open_iterations open_converged
#>           <num>   <lgcl>           <int>         <lgcl>
#> 1:    0.1938111    FALSE              NA             NA
#> 2:    0.2021510    FALSE              NA             NA
#> 3:    0.1978243    FALSE              NA             NA
#> 4:    0.2241836    FALSE              NA             NA
#> 5:    0.2471245    FALSE              NA             NA
#> 6:    0.2149863    FALSE              NA             NA
xt <- mt_example(deaths = "single", years = 2000:2015)
st <- mt_smooth(xt, along = "period", method = "whittaker")
mt_smooth_check(xt, st)
#>    stratum    sex  label    mad_log dispersion trend_ratio correlation      edf
#>      <int> <char> <char>      <num>      <num>       <num>       <num>    <num>
#> 1:       1 female female 0.01237280  0.9935176   0.2413490   0.9996135 434.4595
#> 2:       2   male   male 0.01574177  0.9369526   0.1131588   0.9996879 297.3892
```
