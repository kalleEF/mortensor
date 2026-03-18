# Mean Absolute Scaled Error (MASE) on log(mx)

Compares forecast accuracy to a naive random-walk baseline.

## Usage

``` r
compute_mase(actual_test, predicted_test, train_series)
```

## Arguments

- actual_test:

  Numeric vector of observed test-period rates.

- predicted_test:

  Numeric vector of predicted test-period rates.

- train_series:

  Numeric vector of training-period rates.

## Value

Scalar MASE.
