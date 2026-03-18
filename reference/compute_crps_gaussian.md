# CRPS for a Gaussian predictive distribution (closed-form)

Operates on the log scale: `obs` and `mu` should be log(mx), `sigma` is
the log-scale standard deviation.

## Usage

``` r
compute_crps_gaussian(obs, mu, sigma)
```

## Arguments

- obs:

  Numeric vector of observed log(mx).

- mu:

  Numeric vector of predictive mean log(mx).

- sigma:

  Numeric vector of predictive SD (log scale).

## Value

Scalar mean CRPS.
