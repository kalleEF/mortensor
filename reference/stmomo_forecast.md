# Stochastic Mortality Model (StMoMo) forecast

Fits an APC, CBD, RH, or M6 stochastic mortality model per group using
the StMoMo package and produces forecasts.

## Usage

``` r
stmomo_forecast(
  mx_array,
  pop_array,
  years,
  ages,
  groups_dt,
  h,
  model_type = "APC"
)
```

## Arguments

- mx_array:

  3D array `[n_a x n_t x n_g]` of mortality rates.

- pop_array:

  3D array of population exposure.

- years:

  Integer vector of calendar years.

- ages:

  Integer vector of age-group starts.

- groups_dt:

  data.frame with group definitions.

- h:

  Integer forecast horizon.

- model_type:

  One of `"APC"`, `"CBD"`, `"RH"`, `"M6"`.

## Value

A list with `mx_forecast` `[n_a x h x n_g]` (NULL if StMoMo not
available).
