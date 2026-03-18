# Blend GAM predictions with a parametric tail in log-space

The GAM predicts all target ages (including beyond the open-interval
midpoint). A Kannisto or Gompertz model fitted on GAM predictions at
ages 70+ provides stable old-age behaviour. Logistic weights produce a
C-inf smooth transition: \\w(x) = 1 / (1 + exp(-(x - open\\mid) /
\tau))\\.

## Usage

``` r
blend_parametric_tail(
  mx_gam,
  target_ages,
  open_mid,
  tau = 2,
  method = "kannisto"
)
```

## Arguments

- mx_gam:

  Numeric vector of GAM-predicted mx for all target ages.

- target_ages:

  Integer vector of single-year ages.

- open_mid:

  Centre of the blend (open-interval midpoint, e.g.\\ 95).

- tau:

  Logistic blend width (default 2; larger = wider transition).

- method:

  `"kannisto"` or `"gompertz"`.

## Value

Named list with `mx_blended`, `mx_param`, `blend_weights`.
