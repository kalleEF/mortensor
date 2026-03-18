# Bootstrap prediction intervals for Tucker3 forecasts

Simulates future B-mode paths via ARIMA bootstrap, adds shuffled
historical residuals, and reconstructs mortality rates.

## Usage

``` r
bootstrap_forecast(
  fc_obj,
  n_sim,
  h,
  ages,
  n_groups,
  method = "coda",
  width = 5
)
```

## Arguments

- fc_obj:

  Tucker3 forecast object (with `return_full = TRUE`).

- n_sim:

  Number of bootstrap simulations.

- h:

  Forecast horizon.

- ages:

  Integer vector of age-group starts.

- n_groups:

  Number of groups (dim 3).

- method:

  `"coda"` or `"direct"`.

- width:

  Age-group width.

## Value

List with `mx_lower`, `mx_upper`, `mx_median` (3D arrays).
