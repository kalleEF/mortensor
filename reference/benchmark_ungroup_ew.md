# Ungrouping benchmark on England & Wales (HMD via StMoMo)

Groups observed single-year deaths of men aged 30-100 to 30-34, ...,
85-89 and 90-100, ungroups them with single-year exposures as offset and
compares the result with the observed single-year rates. Errors are
median absolute log errors; the Poisson noise of the observed rates is a
floor (about 0.012 at 30-89 and 0.021 at 90-100 for 1980-2011).

## Usage

``` r
benchmark_ungroup_ew(years = 1980:2011, dims = 1, ...)
```

## Arguments

- years:

  Years used.

- dims:

  1 or 2.

- ...:

  Passed to [`mt_ungroup()`](mt_ungroup.md).

## Value

A list with `male_30_89`, `male_90_100`, `seconds`.
