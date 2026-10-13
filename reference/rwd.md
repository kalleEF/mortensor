# Random walk with drift on log death rates

Benchmark model: the log death rate of every age follows a random walk
with drift, `y[t] = y[t-1] + mu + e[t]`. The drift is estimated from the
last `drift_window` years. Draws include drift uncertainty and
innovations:

## Usage

``` r
rwd(
  drift_window = NULL,
  innovations = c("bootstrap", "gaussian"),
  observation_noise = FALSE
)
```

## Arguments

- drift_window:

  Number of most recent years used for drift and innovations (`NULL`:
  all years).

- innovations:

  `"bootstrap"` or `"gaussian"`.

- observation_noise:

  Model the log rates as a random walk plus noise.

## Value

A model specification.

## Details

- `"bootstrap"` (default): innovations are age vectors resampled from
  the observed year-to-year changes, which keeps the correlation between
  ages; the drift of each path is the mean of a bootstrap sample of
  those changes;

- `"gaussian"`: independent normal innovations per age with the
  estimated variance, and normal drift uncertainty.

Zero death counts are replaced by 0.5 before taking logs.

With `observation_noise = TRUE` the log rate is a random walk with drift
observed with noise (a local level model with drift). Year-to-year
changes of observed rates then contain the noise of two years, and a
plain random walk accumulates it at every step, which makes its
intervals for observed rates too wide when deaths are few or rates
fluctuate. The variances are estimated per age from the variance and
first autocovariance of the changes (the noise variance is the Poisson
part, one over the deaths, plus an extra part common to all ages), the
level is Kalman filtered, and the forecast draws are paths of the level;
predictive draws add the observation noise. In backtests the variant had
a lower CRPS than the plain random walk on simulated, HMD and ONS data,
and calibrated intervals on simulated data; on real data over 5 to 10
years its 95% intervals were too narrow (coverage 0.78 to 0.87), because
changes of trend exceed the uncertainty of a constant drift. The plain
random walk, which is wider, covered 0.92 to 0.97.

## Examples

``` r
x <- mt_example(deaths = "single")
fc <- mt_forecast(mt_fit(rwd(drift_window = 15), x), h = 5, n_draws = 100, seed = 1)
```
