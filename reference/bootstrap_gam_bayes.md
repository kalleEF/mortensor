# GAM Bayesian posterior propagation prediction intervals

Hybrid approach: calibration is done once on the posterior-mean smooth
(standard `mx_smooth`). Per simulation, a new smooth surface is drawn
from the GAM posterior and all ensemble models are re-fit on that draw,
yielding a weighted-average forecast. Quantile PIs are computed from the
resulting simulation cube.

## Usage

``` r
bootstrap_gam_bayes(
  gam_fits,
  ages,
  years,
  pop_arr,
  groups,
  h,
  n_sim,
  ens_names,
  ew,
  calib_coda,
  calib_direct,
  calib_all,
  width = 5L,
  dispatch_fn = dispatch_forecast
)
```

## Arguments

- gam_fits:

  List of [mgcv::gam](https://rdrr.io/pkg/mgcv/man/gam.html) objects,
  one per group.

- ages:

  Numeric vector of age-group starts.

- years:

  Integer vector of observation years.

- pop_arr:

  3D population array \[n_a x n_t x n_g\].

- groups:

  data.table with group information.

- h:

  Forecast horizon.

- n_sim:

  Number of posterior draws.

- ens_names:

  Character vector of ensemble model names.

- ew:

  Named numeric vector of ensemble weights (sum to 1).

- calib_coda:

  Calibrated Tucker3-CoDa parameters.

- calib_direct:

  Calibrated Tucker3-Direct parameters.

- calib_all:

  Calibrated parameters for all models.

- width:

  Age-group width (default 5).

- dispatch_fn:

  Forecast dispatch function (default
  [dispatch_forecast](https://kalleEF.github.io/mortensor/reference/dispatch_forecast.md)).

## Value

List with `mx_lower`, `mx_upper`, `mx_median` (3D arrays).
