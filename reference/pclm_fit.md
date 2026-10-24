# Fit a penalized composite link model

Latent values `gamma = exp(B beta)` on a fine grid are linked to
observed counts by `mu = C (o * gamma)`, where `o` is an offset
(exposures for rates, interval widths for counts). The coefficients are
estimated by penalized Poisson IRLS; the smoothing parameter is chosen
by BIC or AIC on a log scale unless fixed.

## Usage

``` r
pclm_fit(
  y,
  C,
  o,
  B,
  pen,
  lambda = NULL,
  criterion = c("bic", "aic"),
  lambda_range = c(0.01, 1e+08),
  beta = NULL
)
```

## Arguments

- y:

  Observed counts.

- C:

  Composite link matrix (observations x fine cells), dense or sparse.

- o:

  Offset per fine cell.

- B:

  Basis matrix (fine cells x coefficients), or a tensor basis over fine
  ages and periods from `pclm_tensor()`, whose products are computed in
  array form (Currie, Durban and Eilers 2006).

- pen:

  List of penalty matrices (one per smoothing parameter).

- lambda:

  Fixed smoothing parameters, or `NULL` to select them.

- criterion:

  `"bic"` or `"aic"`.

- lambda_range:

  Search range for each smoothing parameter.

- beta:

  Optional starting values.

## Value

A list with `gamma`, `beta`, `mu` (fitted observed counts), `lambda`,
`ed` (effective dimension), `dev`, `ic`, `se_eta`, `converged` and
`at_bound`.

## References

Currie ID, Durban M, Eilers PHC (2006). Generalized linear array models
with applications to multidimensional smoothing. *Journal of the Royal
Statistical Society Series B* 68(2): 259-280.
[doi:10.1111/j.1467-9868.2006.00543.x](https://doi.org/10.1111/j.1467-9868.2006.00543.x)
