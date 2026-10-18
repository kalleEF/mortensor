# Nordpred age-period-cohort model with a power link

The Nordpred method of Møller et al. (2003) for projections of cancer
incidence and mortality: an age-period-cohort model for 5-year age
groups and 5-year periods, fitted by Poisson maximum likelihood with a
power link, whose linear trend (the drift) is projected and damped in
later periods. The deaths `D` of age group `a` in period `p` with
exposure `E` are Poisson with mean

## Usage

``` r
nordpred_apc(
  link = c("power5", "log"),
  n_periods = 4:6,
  recent = NULL,
  cut_trend = c(0, 0.25, 0.5, 0.75, 0.75),
  est_age = NULL,
  use_age = NULL
)
```

## Arguments

- link:

  `"power5"` (Nordpred's power link) or `"log"` (multiplicative Poisson
  age-period-cohort model).

- n_periods:

  Candidate numbers of 5-year periods in the prediction base, at least 3
  each (nordpred's default `4:6`). A single number fixes the base.

- recent:

  Project the slope of the last 10 years (`TRUE`), the average drift
  (`FALSE`), or choose by the test described above (`NULL`).

- cut_trend:

  Proportions of the drift cut in the future 5-year periods.

- est_age:

  Lower bound of the youngest age group in the regression (nordpred's
  `startestage`). Default: the first age group of the grid.

- use_age:

  Lower bound of the youngest age group projected with the regression
  (nordpred's `startuseage`). Default: `est_age`.

## Value

A model specification.

## Details

`E (A[a] + d p + P[p] + C[c])^5`

(`link = "power5"`), with age effects `A`, drift `d`, non-linear period
effects `P` and cohort effects `C` (cohort `c = p - a` on the scale of
5-year steps), or `E exp(A[a] + d p + P[p] + C[c])` (`link = "log"`, the
multiplicative model). The power 5 levels off the exponential growth of
the multiplicative model (Møller et al. 2003). The model is fitted as a
generalised linear model with a Poisson variance function and the power
link `mu^(1/5)` on the rates, written for this package.

The specification follows Møller et al. (2003) and the R package
nordpred (version 1.1, by Møller and Weedon-Fekjær, GPL-2), whose
documentation and code define the details below. This implementation is
written independently of nordpred's code. Its point projections agree
with those of nordpred on the same data to a relative difference below
`1e-7` (see the unit tests).

- **Identification.** The period effects of the first and the last
  period and the cohort effects of the oldest and the youngest cohort
  are zero, so that the drift carries the linear trend (the
  identification of nordpred).

- **Prediction base.** `n_periods` lists the candidate numbers of recent
  5-year periods to fit. The model is fitted to the largest; while its
  deviance rejects the fit at the 1% level (chi-squared goodness of
  fit), the oldest period is dropped, down to the smallest candidate.

- **Trend.** With `recent = NULL` the projected drift is the slope of
  the last 10 years (the change of the period part of the linear
  predictor between the last two periods) if a quadratic period term is
  significant at the 5% level in a Poisson model with log link (age,
  drift, quadratic term and cohort effects), and the average drift `d`
  otherwise.

- **Cut trend.** In future period `j` the drift is multiplied by
  `1 - cut_trend[j]` (cumulated over the periods): the default reduces
  it by 0%, 25%, 50%, 75% and 75% in the five 5-year periods after the
  last observation.

- **Cohorts.** A future cell keeps the effect of its cohort if that
  cohort is older than the cohort of age group `use_age` in the last
  observed period, and takes the effect of that cohort otherwise.

- **Young ages.** Age groups below `use_age` are projected at the mean
  of the observed rates of the last two 5-year periods; age groups below
  `est_age` do not enter the regression.

**Annual data (package-specific).** mortensor models work on single
calendar years. The deaths and exposures of the training years are
summed to 5-year periods that end with the last training year; up to
four earlier years that do not fill a 5-year period are left out. The
model is fitted and projected on these periods. Annual rates are
obtained by linear interpolation of the linear predictor (the fifth root
of the rate, or its logarithm with `link = "log"`) between the midpoints
of the periods: the fitted value of the last observed period at its
midpoint (two years before the last training year) and the projections
of the future periods at theirs. Where the linear predictor changes
linearly over the periods, its mean over the five years of a future
period equals the Nordpred projection of that period. Fitted rates are
the fitted period rates in each year of the periods fitted, and missing
for earlier years and for age groups below `est_age`. The age groups
must be 5 years wide; the last may be open or wider (it counts as one
5-year step for the cohorts, as nordpred's 85+). Use
[`mt_regroup()`](mt_regroup.md) or [`step_regroup()`](mt_prep.md) for
single years of age.

**Uncertainty (package-specific).** Nordpred gives point projections and
documents no method for prediction intervals. Here each draw takes the
coefficients from their asymptotic normal distribution (estimates and
the inverse Fisher information of the Poisson model) and projects them
as above. Draws describe the uncertainty of the coefficients only, with
no variation of future period effects and no overdispersion, so
intervals are narrow when the data are large; check coverage in a
backtest and widen them with [`recalibrated()`](recalibrated.md) if
needed. With the power link, linear predictors below 0.001 of a draw are
set to 0.001 (a safeguard). Age groups below `use_age` have no
uncertainty. Predictive draws (`mt_forecast(predictive = TRUE)`) add
Poisson variation of the deaths given the rates.

Beyond five future periods (25 years) the last value of `cut_trend` is
used (nordpred projects at most five periods). The candidate number of
periods is chosen with the link of the model (nordpred uses the power
link for this choice also when `link = "log"`).

## References

Møller B, Fekjær H, Hakulinen T, Sigvaldason H, Storm HH, Talbäck M,
Haldorsen T (2003). Prediction of cancer incidence in the Nordic
countries: empirical comparison of different approaches. *Statistics in
Medicine* 22(17): 2751-2766.
[doi:10.1002/sim.1481](https://doi.org/10.1002/sim.1481)

## Examples

``` r
x <- mt_example(deaths = "single", years = 1990:2019)
x5 <- mt_regroup(x, age = seq(30, 85, 5))
fit <- mt_fit(nordpred_apc(), x5)
fc <- mt_forecast(fit, h = 10, n_draws = 100, seed = 1)
```
