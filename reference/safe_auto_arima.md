# Safe ARIMA fitting with multiple fallbacks

Tries `auto.arima` first, then falls back to `Arima(0,1,0)`, then to a
trivial constant model `Arima(0,0,0)`.

## Usage

``` r
safe_auto_arima(
  ts_r,
  stepwise = TRUE,
  approximation = TRUE,
  max.p = 2,
  max.q = 1,
  max.d = 1
)
```

## Arguments

- ts_r:

  A `ts` object.

- stepwise:

  Logical; stepwise search?

- approximation:

  Logical; use approximations?

- max.p, max.q, max.d:

  Maximum ARIMA orders.

## Value

An ARIMA model object.
