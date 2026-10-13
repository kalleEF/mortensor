# Simulated Lee-Carter-type data for calibration checks

The generator of the 0.1 review: `log m = a_x + b_x k_t + shift_g + e`,
with `k_t` a random walk with drift (-0.015, sd 0.02) per stratum, iid
noise e (sd 0.03), 5-year age groups 30-34 ... 90-94 and 10 strata (2
sexes x 5 deprivation quintiles). Deaths are Poisson.

## Usage

``` r
simulate_lc_dgp(
  n_t = 40,
  seed = 1,
  pop = 2e+05,
  sd_k = 0.02,
  drift = -0.015,
  sd_obs = 0.03
)
```

## Arguments

- n_t:

  Number of years.

- seed:

  Random seed.

- pop:

  Population at age 30 per stratum.

## Value

An [mt_data](mt_data.md) object.
