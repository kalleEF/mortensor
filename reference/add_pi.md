# Add prediction intervals to a mort_forecast

Computes bootstrap prediction intervals and attaches them to the
forecast object.

## Usage

``` r
add_pi(
  fc,
  n_sim = 500L,
  method = c("bootstrap", "gam_bayes"),
  cap = 1.5,
  calibration = NULL,
  ...
)
```

## Arguments

- fc:

  A `mort_forecast` object.

- n_sim:

  Number of bootstrap simulations (default 500).

- method:

  One of `"bootstrap"` (standard residual/ARIMA bootstrap) or
  `"gam_bayes"` (GAM posterior propagation). The latter requires that
  GAM fits are stored in the forecast extras.

- cap:

  Maximum PI width on log scale (default 1.5).

- calibration:

  Optional `mort_calibration` object (needed for `"gam_bayes"` method to
  access ensemble weights and GAM fits).

- ...:

  Ignored.

## Value

An updated `mort_forecast` with `mx_lower`, `mx_upper`, `mx_median`
filled in.
