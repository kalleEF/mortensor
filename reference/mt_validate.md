# Check mortality data for common problems

Runs plausibility checks and returns one row per finding:

## Usage

``` r
mt_validate(data, hierarchy = NULL, tol = 0.005)
```

## Arguments

- data:

  An [mt_data](mt_data.md) object.

- hierarchy:

  Optional [`mt_hierarchy()`](mt_hierarchy.md), or a list with `var`
  (stratum column) and `total` (value of the total).

- tol:

  Relative tolerance for hierarchy sums.

## Value

A `data.table` of class `mt_validation` with columns `check`,
`severity`, `measure`, `stratum`, `n`, `detail`.

## Details

- deaths with zero exposure in the same cell (error);

- central death rates above 1 below age 100 (warning);

- missing values (warning);

- strata without any value of a measure, for example without deaths
  (warning);

- levels that look like totals of other levels of the same variable, for
  example `"Insgesamt"` or `"all"` among the regions (warning);

- zero exposures (info);

- grid cells without deaths or exposures (info);

- hierarchy sums, when `hierarchy` is given
  ([`mt_hierarchy()`](mt_hierarchy.md)): every aggregate stratum in the
  data is compared with the sum of its bottom strata per cell, for
  deaths and exposures (exposures are summed over the variables whose
  exposures add up and taken once over those whose levels share an
  exposure, such as causes of death).

## Examples

``` r
mt_validate(mt_example())
#> <mt_validation> no findings
```
