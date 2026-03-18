# Bootstrap PIs for Lee-Carter / FDM models

Uses the demography package's point forecast plus residual bootstrap.

## Usage

``` r
bootstrap_lc_fdm(
  mx_array,
  pop_array,
  years,
  ages,
  groups_dt,
  h,
  n_sim,
  model_type = "lca",
  adjust = "dt",
  chooseperiod = FALSE,
  order = 6,
  method = "arima",
  mx_raw = NULL
)
```

## Arguments

- mx_array:

  3D array of (smoothed) mortality rates.

- pop_array:

  3D population array.

- years:

  Integer vector.

- ages:

  Integer vector.

- groups_dt:

  data.frame of groups.

- h:

  Forecast horizon.

- n_sim:

  Number of simulations.

- model_type:

  `"lca"` or `"fdm"`.

- adjust:

  LC adjustment type.

- chooseperiod:

  LC choose-period flag.

- order:

  FDM order.

- method:

  Forecast method for time series.

- mx_raw:

  Optional raw (unsmoothed) rates for residual computation.

## Value

4D array \[n_a x h x n_g x n_sim\].
