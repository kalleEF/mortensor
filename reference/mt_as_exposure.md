# Exposures from population counts

Approximates person-years at risk from population stocks and stores them
as the exposure table:

## Usage

``` r
mt_as_exposure(data, type)
```

## Arguments

- data:

  An [mt_data](mt_data.md) object with a population table.

- type:

  The kind of population counts, a required choice: `"midyear"`
  (mid-year estimates, used as exposures), `"jan1"` (1 January) or
  `"dec31"` (31 December); see Details.

## Value

The `mt_data` object with an exposure table.

## Details

- `"midyear"`: the mid-year population of a calendar year is used as the
  exposure of that year;

- `"jan1"`: the exposure of year t is the mean of the January 1
  populations of years t and t + 1 in the same age interval. The last
  year has no successor and is dropped.

- `"dec31"`: populations on 31 December of year t (as published by the
  German Federal Statistical Office); the exposure of year t is the mean
  of the populations on 31 December of years t - 1 and t. The first year
  has no predecessor and is dropped.

These are the standard approximations of person-years from population
counts (Preston, Heuveline and Guillot 2001): the population at the
middle of the year, or the mean of the populations at its start and its
end. The HMD uses more refined cohort-based formulas. The approximation
used is recorded in the provenance column.

## References

Preston SH, Heuveline P, Guillot M (2001). *Demography: Measuring and
Modeling Population Processes*. Blackwell.

## Examples

``` r
p <- data.frame(year = 2001:2003, age = 70, population = c(1000, 990, 985))
x <- mt_as_exposure(mt_add_population(mt_data(), p), type = "jan1")
mt_obs(x, "exposure")
#>    stratum  label age_lo age_hi period_lo period_hi value      provenance
#>      <int> <char>  <num>  <num>     <num>     <num> <num>          <char>
#> 1:       1    all     70     71      2001      2002 995.0 population_jan1
#> 2:       1    all     70     71      2002      2003 987.5 population_jan1
```
