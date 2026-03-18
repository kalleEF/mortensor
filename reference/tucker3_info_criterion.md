# Information criterion for Tucker3 rank selection

Computes BIC or AIC for a Tucker3 decomposition with given ranks. Much
faster than cross-validation and useful for initial rank screening.

## Usage

``` r
tucker3_info_criterion(
  mx_array,
  R_A,
  R_B,
  R_C,
  method = c("coda", "direct"),
  criterion = c("bic", "aic"),
  ages = NULL,
  width = 5
)
```

## Arguments

- mx_array:

  3D array of mortality rates.

- R_A, R_B, R_C:

  Integer Tucker3 ranks.

- method:

  `"coda"` or `"direct"`.

- criterion:

  `"bic"` (default) or `"aic"`.

- ages:

  Integer vector of age-group starts (required for `"coda"`).

- width:

  Age-group width (default 5).

## Value

A list with `criterion` (BIC/AIC value), `fit_pct`, `n_params`, and
`ranks`.
