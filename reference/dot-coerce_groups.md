# Coerce groups to a validated data.frame

Accepts a data.frame, named list, or character/factor vector and returns
a data.frame with one row per group. At most 3 group columns are
supported. All columns are coerced to character.

## Usage

``` r
.coerce_groups(groups, n_g)
```

## Arguments

- groups:

  Groups specification (data.frame, named list, or character vector).

- n_g:

  Expected number of groups (3rd dimension of mx array).

## Value

A data.frame with character columns.
