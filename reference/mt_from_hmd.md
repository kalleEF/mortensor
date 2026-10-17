# Human Mortality Database tables

`mt_from_hmd()` converts HMD tables in the layout returned by
[`HMDHFDplus::readHMDweb()`](https://rdrr.io/pkg/HMDHFDplus/man/readHMDweb.html)
(columns `Year`, `Age`, `Female`, `Male`, `Total`, and `OpenInterval`)
into an [mt_data](mt_data.md) object with strata `sex`. Without column
`OpenInterval`, an age written with `+` (`"110+"`) is the open interval.
`mt_read_hmd()` downloads `Deaths_1x1` and `Exposures_1x1` for one
country with HMDHFDplus; an HMD account is needed
(<https://www.mortality.org>).

## Usage

``` r
mt_from_hmd(deaths, exposures, sexes = c("Female", "Male"), years = NULL)

mt_read_hmd(
  country,
  username = Sys.getenv("HMD_USER"),
  password = Sys.getenv("HMD_PASSWORD"),
  years = NULL
)
```

## Arguments

- deaths, exposures:

  Data frames from
  [`HMDHFDplus::readHMDweb()`](https://rdrr.io/pkg/HMDHFDplus/man/readHMDweb.html).

- sexes:

  Columns to use.

- years:

  Optional subset of years.

- country:

  HMD country code, for example `"SWE"`.

- username, password:

  HMD credentials (default: environment variables `HMD_USER` and
  `HMD_PASSWORD`).

## Value

An [mt_data](mt_data.md) object.

## Examples

``` r
d <- data.frame(Year = 2000, Age = 0:2, Female = c(10, 1, 1), Male = c(12, 2, 1),
                Total = c(22, 3, 2), OpenInterval = c(FALSE, FALSE, TRUE))
e <- transform(d, Female = 1000, Male = 1000, Total = 2000)
mt_from_hmd(d, e)
#> <mt_data> 2 strata by sex, omega = 110
#>   deaths     single years 0-1, 2+; periods 2000-2000 (width 1); 6 cells
#>   exposure   single years 0-1, 2+; periods 2000-2000 (width 1); 6 cells
#>   grid       <mt_axis: age, 3 intervals> 0 1 2+
#>              <mt_axis: period, 1 interval> 2000
```
