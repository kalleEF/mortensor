# Bohk-Ewald and Rau model of age-specific rates of mortality improvement

The probabilistic model of Bohk-Ewald and Rau (2017) forecasts rates of
mortality improvement (ROMIs) by age, whose trends over time differ
between ages, so that the age pattern of mortality decline changes over
the forecast. For each age `x` and year `t` the ROMI is

## Usage

``` r
bohk_rau(
  core = c("exponential", "linear"),
  limits = c(0.005, 0.035),
  smooth = TRUE,
  draws = c("quantiles", "samples"),
  mcmc = c(chains = 5, iterations = 5200, burn_in = 200, thin = 5)
)
```

## Arguments

- core:

  `"exponential"` or `"linear"` core model.

- limits:

  Lower and upper limit of the forecast ROMIs, or `NULL` for none.

- smooth:

  Smooth the death rates with two-dimensional P-splines before computing
  the ROMIs (`FALSE`: the observed rates, with 0.5 deaths in cells
  without deaths, for example when the rates were smoothed in a
  preprocessing step).

- draws:

  `"quantiles"` (equation 20 of the paper) or `"samples"`.

- mcmc:

  Settings of the sampler, a named vector: `chains`, `iterations` per
  chain (burn-in included), `burn_in` and `thin`. Elements left out keep
  their defaults.

## Value

A model specification.

## Details

`rho[x, t] = 1 - m[x, t + 1] / m[x, t]`

(equation 1 of the paper; positive when mortality falls), computed from
death rates smoothed over age and time with two-dimensional P-splines
(Currie, Durban and Eilers 2004). One of two Bayesian core models
describes the ROMIs of every age over the years of the base period:

- `core = "exponential"` (the core model of the paper's applications):
  `rho[x, t] ~ N(exp(b1[x] + b2[x] log t), sigma^2)`,
  `b1[x] ~ N(theta1[x], sigma1^2)`, `b2[x] ~ N(theta2[x], sigma2^2)`,
  with `theta1[x]` and `theta2[x]` the least-squares estimates of
  `log rho[x, t] = theta1[x] + theta2[x] log t` for each age and
  `sigma`, `sigma1`, `sigma2` uniform on (0, 1). Improvements can level
  off; mortality increases cannot be forecast.

- `core = "linear"`: `rho[x, t] ~ N(b1[x] + b2[x] t, sigma^2)` with
  `(b1[x], b2[x])` bivariate normal with means `mu1`, `mu2` and
  covariance from standard deviations `omega1`, `omega2` and correlation
  `r`; uniform priors `sigma`, `omega1`, `omega2` on (0, 1), `mu1`,
  `mu2` on (-0.1, 0.1) and `r` on (-1, 1). Trends of the improvements
  continue linearly.

The posterior and the future ROMIs (each one the curve of a posterior
sample plus normal noise with that sample's `sigma`) are sampled by
Markov chain Monte Carlo. Forecast ROMIs outside `limits` are set to the
nearest limit (the paper uses 0.005 and 0.035: mortality declines by at
least 0.5% and at most 3.5% a year; it excludes implausible forecasts
"by setting thresholds", which this package reads as setting the values
to the thresholds). Death rates follow from the last smoothed rates of
the base period by `m[x, t] = m[x, t - 1] (1 - rho[x, t - 1])` (equation
20). With `draws = "quantiles"` (the paper's method) draw `j` of `n`
applies the quantile `(j - 0.5) / n` of the forecast ROMIs at every age
and year, so that the draws give the quantiles of the death rates as in
the paper. These draws are comonotone across ages and years; their order
is permuted at random, so that the draws of different strata are
independent. With `draws = "samples"` (package-specific) every draw is
the path of one posterior sample with its own noise. Such paths give
much narrower intervals of summary measures such as life expectancy,
because the yearly noise of the improvements averages out over ages and
years.

The defaults follow the paper's application (its Table 3): the
exponential core model, five chains of 5,200 iterations with the first
200 discarded and every fifth kept (5,000 samples), and the limits 0.005
and 0.035. The paper recommends a base period at least as long as the
forecast horizon, and uses single years of age from 0.

**What differs from the paper.** The paper samples with JAGS; here a
sampler written for this package is used (Metropolis steps for the curve
parameters of each age, with proposals scaled by the curvature at the
least-squares estimates and adapted during the burn-in; conjugate draws
of the variances truncated to the prior range; moves that scale `sigma1`
or `sigma2` together with the deviations of the curve parameters from
`theta`; Gibbs and slice sampling steps for the linear core model). The
paper does not state the following, which are choices of this package
(package-specific, reference unknown):

- time `t` counts the years of the base period from 1 (the change from
  the first to the second year) and continues into the forecast;

- ROMIs below the lower limit (0.001 with `limits = NULL` or a lower
  limit of 0 or below) are set to that limit in the least-squares fit of
  `theta` (the paper sets negative values "to a minimum that is greater
  than zero");

- the P-spline surface has cubic B-splines with knots every 5 years of
  age and time and second-order penalties, with both smoothing
  parameters chosen by BIC (the paper uses the package MortalitySmooth);
  when the grid starts with age 0 as its own interval, the rates of age
  0 are smoothed over time only (as [`mt_smooth()`](mt_smooth.md) leaves
  out the first year of life);

- forecast ROMIs are kept below 1 so that rates stay positive.

The model's option to let a population converge to the mortality trend
of reference populations, by weights that shift from the population to
the references over the forecast years, is not implemented: the paper
leaves the weighting scheme open. Each stratum is forecast from its own
data.

Draws describe the expected (smoothed) rates; predictive draws
(`mt_forecast(predictive = TRUE)`) add Poisson variation of the deaths
given the rates. Fitted rates are the smoothed rates (the observed rates
with `smooth = FALSE`). The Gelman-Rubin statistics of the scalar
parameters and the largest one of the curve parameters are reported by
[`mt_parameters()`](mt_parameters.md) (`rhat`); values above 1.1 give a
warning. With the default settings a fit takes about 6 seconds per
stratum for 71 ages and 36 years.

## References

Bohk-Ewald C, Rau R (2017). Probabilistic mortality forecasting with
varying age-specific survival improvements. *Genus* 73(1): 1.
[doi:10.1186/s41118-016-0017-8](https://doi.org/10.1186/s41118-016-0017-8)

Currie ID, Durban M, Eilers PHC (2004). Smoothing and forecasting
mortality rates. *Statistical Modelling* 4(4): 279-298.

## Examples

``` r
x <- mt_subset(mt_example(deaths = "single", years = 1995:2019), sex == "female")
# shorter chains than the defaults of the paper, for a quick example
spec <- bohk_rau(mcmc = c(chains = 2, iterations = 3200, burn_in = 200, thin = 3))
set.seed(1)
fit <- mt_fit(spec, x)
fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1)
life_expectancy(fc)
#>    stratum    sex  label period     mean   q0.025     q0.1     q0.5     q0.9
#>      <int> <char> <char>  <num>    <num>    <num>    <num>    <num>    <num>
#> 1:       1 female female   2020 55.91851 55.91600 55.91685 55.91851 55.92017
#> 2:       1 female female   2021 56.01064 56.00565 56.00732 56.01064 56.01397
#> 3:       1 female female   2022 56.10238 56.09490 56.09740 56.10238 56.10735
#> 4:       1 female female   2023 56.19371 56.18379 56.18708 56.19372 56.20033
#> 5:       1 female female   2024 56.28466 56.27227 56.27638 56.28466 56.29291
#>      q0.975    at
#>       <num> <num>
#> 1: 55.92100    30
#> 2: 56.01563    30
#> 3: 56.10985    30
#> 4: 56.20365    30
#> 5: 56.29705    30
```
