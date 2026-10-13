# Estimate suppressed death counts

Official statistics often suppress small cells for confidentiality (for
example `"."` in tables of the German Federal Statistical Office, `"<3"`
or `"<5"` elsewhere) while publishing totals: deaths of a region over
all causes, deaths of the whole country by cause. `mt_fill_suppressed()`
estimates the missing death counts of the bottom strata from two sources
of information:

## Usage

``` r
mt_fill_suppressed(
  data,
  hierarchy = NULL,
  lower = 0,
  upper = Inf,
  maxit = 10000,
  tol = 1e-06
)
```

## Arguments

- data:

  An [mt_data](mt_data.md) object in which some death cells are missing
  (cells that other strata have; suppressed values are left out when the
  data are read).

- hierarchy:

  An [`mt_hierarchy()`](mt_hierarchy.md) that describes the published
  totals in the data, or `NULL` when there are none.

- lower, upper:

  Bounds of the suppressed counts.

- maxit, tol:

  Iterations and tolerance (in deaths) of the raking; cells are raked
  until no estimate changes by more than `tol`.

## Value

The `mt_data` object with filled deaths; a summary in attribute `"fill"`
(cells filled, cells left missing, inconsistent cells, iterations,
largest relative gap to a published total).

## Details

- **published totals** (`hierarchy`): the missing cells of every
  aggregate with a known total must add up to that total minus its known
  cells. They are estimated by iterative proportional fitting (raking):
  starting values, the expected deaths at the pooled death rate of the
  same age and period, are scaled constraint by constraint until every
  total is matched. With two crossed variables (regions and causes) the
  result is the log-linear independence fit of the missing cells given
  both margins. Missing aggregates are then filled as sums of their
  parts.

- **bounds of the suppressed counts** (`lower`, `upper`): the range that
  the suppression rule implies, for example `lower = 1, upper = 2` when
  counts of 1 and 2 are suppressed, or `upper = 4` for `"<5"`. Estimates
  are kept inside the bounds during raking. Cells that no published
  total covers get the expected deaths at the pooled rate, limited to
  the bounds; this needs a finite `upper`.

Without bounds, cells that no published total covers stay missing (with
a warning). The estimates are plausible values that agree with the
published information; they carry no information beyond it, so forecasts
of strata with many estimated cells are less certain than their
intervals suggest. Filled cells are marked with provenance
`"estimated"`. When the known cells of an aggregate already exceed its
total (inconsistent data), the missing cells of that aggregate are set
to the lower bound and the case is reported.

## Examples

``` r
d <- data.frame(region = c("A", "B", "C", "all"), year = 2020, age = "90+",
                deaths = c(10, NA, NA, 40), exposure = c(100, 200, 200, 500))
x <- mt_data(deaths = d, exposure = d, by = "region")
#> Warning: 2 missing deaths value(s) left out; suppressed counts can be estimated with mt_fill_suppressed().
y <- mt_fill_suppressed(x, mt_hierarchy(region = "all"))
mt_obs(y)
#>    stratum region  label age_lo age_hi period_lo period_hi value provenance
#>      <int> <char> <char>  <num>  <num>     <num>     <num> <num>     <char>
#> 1:       1      A      A     90    Inf      2020      2021    10   observed
#> 2:       2      B      B     90    Inf      2020      2021    15  estimated
#> 3:       3      C      C     90    Inf      2020      2021    15  estimated
#> 4:       4    all    all     90    Inf      2020      2021    40   observed
attr(y, "fill")
#>    filled missing inconsistent iterations max_gap
#>     <int>   <int>        <int>      <int>   <num>
#> 1:      2       0            0          2       0
# counts below 3 suppressed, no totals: expected deaths limited to 1-2
d2 <- data.frame(region = c("A", "B"), year = 2020, age = "90+", deaths = c(30, NA),
                 exposure = c(1000, 50))
z <- mt_fill_suppressed(mt_data(deaths = d2, exposure = d2, by = "region"), lower = 1, upper = 2)
#> Warning: 1 missing deaths value(s) left out; suppressed counts can be estimated with mt_fill_suppressed().
mt_obs(z)
#>    stratum region  label age_lo age_hi period_lo period_hi value provenance
#>      <int> <char> <char>  <num>  <num>     <num>     <num> <num>     <char>
#> 1:       1      A      A     90    Inf      2020      2021  30.0   observed
#> 2:       2      B      B     90    Inf      2020      2021   1.5  estimated
```
