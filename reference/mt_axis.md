# Age and period axes

An axis is an ordered set of contiguous half-open intervals
`[lower, upper)`. Only the last interval may be open (`upper = Inf`), as
for an age group such as 90+. Axes describe the resolution of observed
data and of the modelling grid. Period axes may have gaps (for example
the census years 2001, 2011 and 2021), which describe observed data;
models need consecutive years.

## Usage

``` r
mt_axis(lower, upper, labels = NULL, kind = c("age", "period"))

single_ages(from = 0, to = 100, open = TRUE)

age_groups(breaks, open = TRUE, last_upper = NULL)

single_years(from, to)

period_groups(breaks, last_upper)
```

## Arguments

- lower, upper:

  Numeric vectors of interval bounds. `upper` may end with `Inf` for an
  open last interval.

- labels:

  Optional character labels; generated when `NULL`.

- kind:

  `"age"` or `"period"`.

- from, to:

  First and last value (inclusive).

- open:

  Whether the last interval is open (`to+`).

- breaks:

  Lower bounds of the age groups.

- last_upper:

  Upper bound of the last group when `open = FALSE`.

## Value

An object of class `mt_axis`.

## Examples

``` r
single_ages(30, 100)
#> <mt_axis: age, 71 intervals> 30 31 32 33 ... 98 99 100+
age_groups(c(0, 1, seq(5, 90, 5)))
#> <mt_axis: age, 20 intervals> 0 1-4 5-9 10-14 ... 80-84 85-89 90+
single_years(2001, 2019)
#> <mt_axis: period, 19 intervals> 2001 2002 2003 2004 ... 2017 2018 2019
```
