# Create a mort_data object from pre-built arrays

Convenience wrapper around
[`mort_data()`](https://kalleEF.github.io/mortensor/reference/mort_data.md)
for users who already have a 3D mortality-rate array (e.g. from the
Human Mortality Database, the `demography` package, or any other
source).

## Usage

``` r
read_mortality_from_arrays(
  mx,
  ages,
  years,
  groups,
  pop = NULL,
  cause = NULL,
  age_labels = NULL
)
```

## Arguments

- mx:

  3D array `[n_ages x n_years x n_groups]` of mortality rates.

- ages:

  Integer vector of age-group start values.

- years:

  Integer vector of calendar years.

- groups:

  A data.frame with one row per group (3rd dimension of `mx`). Column
  names define the strata (e.g. `sex`, `region`).

- pop:

  Optional 3D array of population exposure (same dims as `mx`).

- cause:

  Optional cause-of-death label.

- age_labels:

  Optional character vector of age-group labels.

## Value

A
[mort_data](https://kalleEF.github.io/mortensor/reference/mort_data.md)
object.
