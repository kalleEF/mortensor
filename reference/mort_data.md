# Create a mortality data object

Constructs a `mort_data` S3 object that encodes stratified mortality
data as a 3D array with metadata about ages, years, and groups (strata).

## Usage

``` r
mort_data(
  mx,
  ages,
  years,
  groups,
  pop = NULL,
  cause = NULL,
  age_labels = NULL,
  smoothed = FALSE
)
```

## Arguments

- mx:

  3D array `[n_ages x n_years x n_groups]` of mortality rates.

- ages:

  Integer vector of age-group start values (e.g. `c(30, 35, 40, 90)`).

- years:

  Integer vector of calendar years.

- groups:

  Group metadata. Accepted formats:

  - A `data.frame` with one row per group and columns defining strata
    (e.g. `data.frame(sex = c("M","F"), quintile = c(1,2))`). This is
    the recommended format for flexible plotting with
    `facet_by`/`color_by`.

  - A named list with equal-length vectors (converted to data.frame).

  - A character vector of group labels (stored as single column
    `"group"` — faceting by sub-dimensions will not be possible).

  At most 3 group columns are supported.

- pop:

  Optional 3D array of population exposure (same dims as `mx`).

- cause:

  Optional character label for the cause of death.

- age_labels:

  Optional character vector of age-group labels (e.g. `"30-34"`,
  `"35-39"`). Derived from `ages` if omitted.

- smoothed:

  Logical; has the data been smoothed?

## Value

A `mort_data` object (list with class attribute).
