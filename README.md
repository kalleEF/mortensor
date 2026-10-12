# mortensor

<!-- badges: start -->
[![R-CMD-check](https://github.com/kalleEF/mortensor/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kalleEF/mortensor/actions/workflows/R-CMD-check.yaml)
[![validation](https://github.com/kalleEF/mortensor/actions/workflows/validation.yaml/badge.svg)](https://github.com/kalleEF/mortensor/actions/workflows/validation.yaml)
<!-- badges: end -->

Flexible probabilistic mortality forecasting from data observed at different
granularities: age groups and open age intervals, single or grouped periods,
population strata and causes.

**Status:** version 0.2 is a rewrite in progress (see `dev/progress.md`).
Version 0.1.1 is available from the commit tagged in `dev/progress.md`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("kalleEF/mortensor")
```

Requires R >= 4.3.0. Hard dependencies are data.table and Matrix; models and
data adapters that build on other packages (forecast, bamp, StMoMo,
demography, ggplot2) use them when installed.

## Example

```r
library(mortensor)

# simulated data: deaths in 5-year age groups and 90+, exposures by single year
x <- mt_example()
x
mt_validate(x)

# estimate single-year deaths with the penalized composite link model,
# using the single-year exposures as offset
u <- mt_ungroup(x)
attr(u, "ungroup")

# Human Mortality Database data for England & Wales (via StMoMo)
ew <- as_mt_data(StMoMo::EWMaleData)
```

## Development

See `CLAUDE.md` for conventions and commands, `validation/README.md` for the
validation pipeline, and `dev/implementation-plan.md` for the plan.

## Licence

GPL-3.
