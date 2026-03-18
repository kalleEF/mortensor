# Reconcile multiple cause-specific forecasts to sum to an all-cause forecast

Takes independently produced `mort_forecast` objects (one per cause) and
an all-cause forecast, then applies WLS reconciliation so disease
forecasts sum to the all-cause total. Each reconciled result remains a
valid `mort_forecast`.

## Usage

``` r
reconcile(forecasts, all_cause, recon_cap = 0.5)
```

## Arguments

- forecasts:

  Named list of `mort_forecast` objects (one per disease).

- all_cause:

  A `mort_forecast` object for all-cause mortality.

- recon_cap:

  Soft cap: max \|delta_k\| as fraction of point forecast (default 0.5).

## Value

A named list of reconciled `mort_forecast` objects.
