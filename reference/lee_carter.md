# Lee-Carter model

`log m[x, t] = a[x] + b[x] k[t]` with the constraints `sum(b) = 1` and
`sum(k) = 0`, fitted per stratum:

## Usage

``` r
lee_carter(
  method = c("poisson", "svd"),
  jump_off = c("fitted", "observed"),
  bootstrap = FALSE,
  n_boot = 100,
  adjust = c("dt", "dxt", "e0", "none"),
  choose_period = FALSE,
  min_years = 20
)
```

## Arguments

- method:

  `"poisson"` or `"svd"`.

- jump_off:

  `"fitted"` (forecasts start from the fitted rates of the last year) or
  `"observed"` (from the observed rates of the last year).

- bootstrap:

  Semi-parametric bootstrap of the parameters.

- n_boot:

  Number of bootstrap fits.

- adjust:

  Re-estimation of `k[t]` for `method = "svd"`: `"dt"`, `"dxt"`, `"e0"`
  or `"none"`.

- choose_period:

  Choose the fitting period (Booth, Maindonald and Smith 2002).

- min_years:

  `choose_period` considers fitting periods of more than `min_years`
  years (`minperiod` of
  [`demography::lca()`](https://pkg.robjhyndman.com/demography/reference/lca.html)).

## Value

A model specification.

## Details

- `method = "poisson"`: Poisson log-bilinear model (Brouhns, Denuit and
  Vermunt 2002), maximum likelihood by Newton updates;

- `method = "svd"`: original method (Lee and Carter 1992): singular
  value decomposition of centred log rates, with `k[t]` re-estimated in
  a second stage (`adjust`, as in
  [`demography::lca()`](https://pkg.robjhyndman.com/demography/reference/lca.html)):
  so that the fitted deaths match the observed total deaths of each year
  (`"dt"`, Lee and Carter 1992), by Poisson maximum likelihood of the
  deaths by age given `a` and `b` (`"dxt"`, Booth, Maindonald and Smith
  2002), so that the fitted life expectancy at the first age matches the
  observed one (`"e0"`, Lee and Miller 2001), or not at all (`"none"`).

`choose_period = TRUE` fits the model to the last years only, starting
in the year that minimises the ratio of the mean deviance of a fit with
a linear `k[t]` to the mean deviance of the Lee-Carter fit (Booth,
Maindonald and Smith 2002): the longest recent period over which the
decline was close to linear
(`demography::lca(chooseperiod = TRUE, breakmethod = "bms")` gives the
same period). Fitted rates of earlier years are missing.

The period index follows a random walk with drift. Draws include the
innovations and the uncertainty of the drift. With `bootstrap = TRUE`
the parameters are re-estimated on `n_boot` Poisson resamples of the
fitted deaths (semi-parametric bootstrap of Brouhns et al. 2005), which
adds the uncertainty of `a` and `b`.

## References

Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality.
*Journal of the American Statistical Association* 87(419): 659-671.

Brouhns N, Denuit M, Vermunt JK (2002). A Poisson log-bilinear
regression approach to the construction of projected lifetables.
*Insurance: Mathematics and Economics* 31(3): 373-393.
[doi:10.1016/S0167-6687(02)00185-3](https://doi.org/10.1016/S0167-6687%2802%2900185-3)

Booth H, Maindonald J, Smith L (2002). Applying Lee-Carter under
conditions of variable mortality decline. *Population Studies* 56(3):
325-336.
[doi:10.1080/00324720215935](https://doi.org/10.1080/00324720215935)

Lee R, Miller T (2001). Evaluating the performance of the Lee-Carter
method for forecasting mortality. *Demography* 38(4): 537-549.
[doi:10.1353/dem.2001.0036](https://doi.org/10.1353/dem.2001.0036)

Brouhns N, Denuit M, Van Keilegom I (2005). Bootstrapping the Poisson
log-bilinear model for mortality forecasting. *Scandinavian Actuarial
Journal* 2005(3): 212-224.

## Examples

``` r
x <- mt_example(deaths = "single")
fit <- mt_fit(lee_carter(), x)
fc <- mt_forecast(fit, h = 10, n_draws = 200, seed = 1)
```
