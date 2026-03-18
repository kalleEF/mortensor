# Inverse CLR transform to death distribution

Converts CLR (centered log-ratio) coordinates back to a dx composition
using the softmax (inverse CLR) transform.

## Usage

``` r
inv_clr_to_dx(clr_vec)
```

## Arguments

- clr_vec:

  Numeric vector of CLR coordinates.

## Value

Numeric vector of dx values (positive, sums to 1).
