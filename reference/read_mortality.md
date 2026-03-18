# Read mortality data from a long-format table

The primary entry point for loading mortality data into mortensor.
Accepts any long-format `data.frame` or `data.table` containing
mortality rates (or death counts + population) with columns for year,
age, and optionally one or more stratification variables (e.g. sex,
region, deprivation).

## Usage

``` r
read_mortality(
  dt,
  mx_col = NULL,
  ages,
  years = NULL,
  group_cols = NULL,
  pop_col = NULL,
  deaths_col = NULL,
  cause = NULL,
  age_labels = NULL,
  year_col = "year",
  age_col = "age_group"
)
```

## Arguments

- dt:

  A `data.frame` or `data.table` in long format.

- mx_col:

  Character name of the mortality-rate column. Either `mx_col` or both
  `deaths_col` and `pop_col` must be provided.

- ages:

  Integer vector of age-group start values (e.g. `c(30, 35, 40, 90)`),
  in the same order used to sort the age column.

- years:

  Integer vector of calendar years to keep (default: all in `dt`).

- group_cols:

  Character vector of column names defining groups/strata (e.g.
  `c("sex", "region")`). If `NULL`, auto-detected as all columns that
  are not `year_col`, `age_col`, `mx_col`, `deaths_col`, or `pop_col`.

- pop_col:

  Optional character name of a population-exposure column. Required when
  using `deaths_col`.

- deaths_col:

  Optional character name of a death-count column. If provided together
  with `pop_col`, mortality rates are computed as `deaths / population`.

- cause:

  Optional cause-of-death label (e.g. `"IHD"`).

- age_labels:

  Optional character vector of age-group labels (e.g.
  `c("30-34", "35-39")`). If `NULL`, labels are derived from `ages`.

- year_col:

  Character name of the year column (default `"year"`).

- age_col:

  Character name of the age column (default `"age_group"`).

## Value

A
[mort_data](https://kalleEF.github.io/mortensor/reference/mort_data.md)
object.

## Expected input format

A long-format table with one row per year x age x group combination.
Required columns (names configurable):

- Year column:

  Integer calendar year (e.g. 2001, 2002). Specify via `year_col`.

- Age column:

  Age-group identifier, either a label like `"30-34"` or an integer
  start age like `30`. Specify via `age_col`.

- Rate column or deaths + population columns:

  Either a single column of pre-computed mortality rates (specify via
  `mx_col`), or two columns: one for death counts (`deaths_col`) and one
  for population exposure (`pop_col`), in which case
  `mx = deaths / population` is computed internally.

- Group columns (optional):

  Any additional columns that define strata (e.g. `"sex"`, `"region"`).
  Specify via `group_cols`. If `NULL`, `age_col` and `year_col` are
  auto-detected and all remaining non-value columns are treated as group
  columns.

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal example - rates already computed, single stratum
dt <- data.table::data.table(
  year      = rep(2001:2020, each = 13),
  age_group = rep(seq(30, 90, 5), 20),
  mx        = runif(260, 0.001, 0.05)
)
md <- read_mortality(dt, mx_col = "mx", ages = seq(30, 90, 5))

# With deaths + population, two strata
dt2 <- data.table::data.table(
  year   = rep(2001:2020, each = 26),
  age    = rep(rep(seq(30, 90, 5), each = 2), 20),
  sex    = rep(c("Male", "Female"), 260),
  deaths = rpois(520, 50),
  pop    = rpois(520, 10000)
)
md2 <- read_mortality(dt2, deaths_col = "deaths", pop_col = "pop",
                      ages = seq(30, 90, 5), age_col = "age",
                      group_cols = "sex")
} # }
```
