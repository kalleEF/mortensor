# Life Tables and Life Expectancy

## Overview

Life tables convert age-specific mortality rates into interpretable
summary measures — most importantly, life expectancy. **mortensor**
provides:

- [`life_table()`](https://kalleEF.github.io/mortensor/reference/life_table.md):
  Construct a full abridged period life table from a vector of
  age-specific rates
- [`life_expectancy()`](https://kalleEF.github.io/mortensor/reference/life_expectancy.md):
  S3 generic computing life expectancy at the youngest observed age,
  with methods for `mort_data` (observed) and `mort_forecast` (forecast,
  optionally with PIs)
- [`compute_dx()`](https://kalleEF.github.io/mortensor/reference/compute_dx.md)
  /
  [`dx_to_mx()`](https://kalleEF.github.io/mortensor/reference/dx_to_mx.md):
  Convert between mortality rates and the death distribution (used
  internally by Tucker3-CoDa)

## Abridged life tables: `life_table()`

### Basic usage

``` r
library(mortensor)

# A vector of age-specific mortality rates (13 five-year groups, 30–90+)
mx <- c(0.00052, 0.00071, 0.00112, 0.00198, 0.00356,
        0.00621, 0.01102, 0.01987, 0.03542, 0.06011,
        0.09876, 0.15432, 0.24567)
ages <- seq(30, 90, by = 5)

lt <- life_table(mx, ages = ages, width = 5)
print(lt)
```

### Output columns

The returned `data.frame` contains the standard life-table columns:

| Column | Name    | Description                              |
|--------|---------|------------------------------------------|
| `age`  | Age     | Start of age interval                    |
| `n`    | Width   | Interval width (open = `Inf`)            |
| `mx`   | $m_{x}$ | Central death rate                       |
| `qx`   | $q_{x}$ | Probability of dying in interval         |
| `px`   | $p_{x}$ | Probability of surviving ($1 - q_{x}$)   |
| `lx`   | $l_{x}$ | Number surviving to age $x$ (from radix) |
| `dx`   | $d_{x}$ | Deaths in interval                       |
| `Lx`   | $L_{x}$ | Person-years lived in interval           |
| `Tx`   | $T_{x}$ | Total person-years remaining             |
| `ex`   | $e_{x}$ | Life expectancy at age $x$               |

### The $\left. m_{x}\rightarrow q_{x} \right.$ conversion

For closed intervals of width $n$, the standard actuarial conversion is:

$$q_{x} = \frac{n \cdot m_{x}}{1 + \left( n - a_{x} \right) \cdot m_{x}}$$

where $a_{x}$ is the average years lived by those dying in the interval
(defaulting to $n/2$ for adults). For the open interval, $q_{x} = 1$.

### Parameters

``` r
lt <- life_table(
  mx              = mx,
  ages            = ages,
  width           = 5,                  # scalar or vector of widths
  radix           = 100000,             # initial cohort size
  open_age_method = "reciprocal"        # "reciprocal", "kannisto", or "gompertz"
)
```

### Open age interval methods

The open interval (e.g. 90+) requires special treatment since it has
infinite width. Three methods are available:

**Reciprocal** (default): $e_{90 +} = 1/m_{90 +}$. Simple and widely
used.

**Kannisto**: Fits the Kannisto logistic model
$\mu(x) = ae^{bx}/\left( 1 + ae^{bx} \right)$ to the last few age groups
and integrates the survival function to obtain $e_{90 +}$.

**Gompertz**: Fits $\log\mu(x) = \alpha + \beta x$ and integrates.

``` r
# Compare open-age methods
lt_recip   <- life_table(mx, ages, width = 5, open_age_method = "reciprocal")
lt_kann    <- life_table(mx, ages, width = 5, open_age_method = "kannisto")
lt_gomp    <- life_table(mx, ages, width = 5, open_age_method = "gompertz")

cat(sprintf("e30: reciprocal=%.2f, kannisto=%.2f, gompertz=%.2f\n",
            lt_recip$ex[1], lt_kann$ex[1], lt_gomp$ex[1]))
```

The differences are typically small (\< 0.5 years) for well-behaved data
but can matter for cause-specific mortality where $m_{90 +}$ is
uncertain.

## Life expectancy from observed data

The
[`life_expectancy()`](https://kalleEF.github.io/mortensor/reference/life_expectancy.md)
method for `mort_data` computes $e_{x}$ at the youngest observed age for
every year and group:

``` r
data(example_mortality)

le_obs <- life_expectancy(example_mortality, width = 5,
                          open_age_method = "reciprocal")
head(le_obs)
#    year      group     ex
# 1: 2001   Male_Q1  44.23
# 2: 2001   Male_Q2  46.12
# 3: 2001 Female_Q1  48.67
# ...
```

This is useful for:

- Tracking historical trends in life expectancy
- Validating that smoothed data produces sensible life tables
- Comparing life expectancy across population strata (e.g. deprivation
  gaps)

``` r
sm <- smooth_mortality(example_mortality, method = "gam",
                       covid_years = c(2020L, 2021L))

le_smooth <- life_expectancy(sm)
# Compare smoothed vs raw
merge(le_obs, le_smooth, by = c("year", "group"),
      suffixes = c("_raw", "_smooth"))
```

## Life expectancy from forecasts

The
[`life_expectancy()`](https://kalleEF.github.io/mortensor/reference/life_expectancy.md)
method for `mort_forecast` computes $e_{x}$ for forecast years, and
propagates prediction intervals when available:

``` r
cal <- calibrate(sm, holdout_frac = 0.2, ensemble = TRUE)
fc  <- forecast(cal, h = 10)
fc  <- add_pi(fc, n_sim = 500, method = "bootstrap")

le_fc <- life_expectancy(fc, width = 5)
head(le_fc)
#    year      group     ex  ex_lower  ex_upper
# 1: 2021   Male_Q1  47.12     45.89     48.31
# 2: 2021   Male_Q2  48.56     47.21     49.87
# ...
```

### How PIs propagate to life expectancy

When the forecast has a `sims` array (from
[`add_pi()`](https://kalleEF.github.io/mortensor/reference/add_pi.md)),
the method:

1.  For each simulation $b = 1,\ldots,B$: extract the simulated
    mortality surface $m_{x,t,g}^{(b)}$ and compute a full life table →
    $e_{x}^{(b)}$.
2.  Take the 2.5th and 97.5th percentiles across $B$ simulations.

This correctly propagates age-correlated uncertainty into the life
expectancy PI.

### Life expectancy from disaggregated data

After disaggregation to single-year ages, the life table uses width = 1:

``` r
fc_1yr <- disaggregate(fc, target_ages = 30:100,
                       method = "natural_gompertz")
le_1yr <- life_expectancy(fc_1yr)

# Compare 5-year vs 1-year life expectancy
cat(sprintf("e30 from 5-year groups: %.2f\n", le_fc$ex[1]))
cat(sprintf("e30 from single years:  %.2f\n", le_1yr$ex[1]))
```

Disaggregated life tables are generally more accurate because they avoid
the $a_{x}$ assumption (average years lived in interval) that abridged
tables require.

## Death distributions: `compute_dx()` and `dx_to_mx()`

These functions convert between mortality rates and the life-table death
distribution, used internally by Tucker3-CoDa:

### `compute_dx()`

``` r
# Convert mortality rates to a death distribution
dx <- compute_dx(mx, ages = seq(30, 90, 5), width = 5L)
print(dx)
# Numeric vector summing to ~1, representing the share of deaths at each age

sum(dx)  # ≈ 1.0
```

The death distribution $d_{x}$ is:

$$d_{x} = l_{x} \cdot q_{x}/l_{0}$$

where $l_{x}$ and $q_{x}$ come from the standard life table. This
converts rates to a compositional (simplex) representation, which is the
input to the CLR transform in Tucker3-CoDa.

### `dx_to_mx()`

``` r
# Invert: death distribution → mortality rates
mx_recovered <- dx_to_mx(dx, ages = seq(30, 90, 5), width = 5L)

# Should match original (up to numerical precision)
max(abs(log(mx) - log(mx_recovered)))  # < 1e-10
```

The inversion uses Gompertz extrapolation for the open age interval to
recover $m_{90 +}$ from $d_{90 +}$.

### Why use death distributions?

The Tucker3-CoDa method works on death distributions rather than rates
because:

1.  **Compositional**: $d_{x}$ sums to 1, so we can apply log-ratio
    transforms (CLR/ILR) that respect the simplex geometry
2.  **Level-invariant**: The distribution captures the *shape* of
    mortality across ages, not the overall level — forecasting shape
    changes is often more stable than forecasting levels
3.  **Bounded**: Unlike log(mx), CLR-transformed $d_{x}$ avoids issues
    with near-zero rates at young ages

## Building custom life tables

### Variable-width age groups

``` r
# Some data have irregular groupings (e.g. 0, 1-4, 5-9, 10-14, ...)
mx_infant <- c(0.005, 0.001, 0.0003, 0.0002, 0.00015, 0.00035,
               0.00082, 0.0021, 0.0054, 0.014, 0.035, 0.089)
ages_wide <- c(0, 1, 5, 10, 15, 25, 35, 45, 55, 65, 75, 85)
widths    <- c(1, 4, 5, 5, 10, 10, 10, 10, 10, 10, 10, Inf)

lt_wide <- life_table(mx_infant, ages = ages_wide, width = widths)
```

### Single-year life tables

``` r
# After disaggregation, rates are at single-year ages
mx_1yr <- runif(71, 0.0005, 0.3)  # dummy single-year rates for ages 30–100
lt_1yr <- life_table(mx_1yr, ages = 30:100, width = 1)
lt_1yr$ex[1]  # life expectancy at age 30
```

## Complete example: historical and forecast life expectancy

``` r
library(mortensor)
data(example_mortality)

# --- Historical life expectancy ---
le_hist <- life_expectancy(example_mortality)

# --- Forecast pipeline ---
sm  <- smooth_mortality(example_mortality, method = "gam",
                        covid_years = c(2020L, 2021L))
cal <- calibrate(sm, holdout_frac = 0.2, ensemble = TRUE)
fc  <- forecast(cal, h = 10)
fc  <- add_pi(fc, n_sim = 500, method = "bootstrap")

# --- Forecast life expectancy (with PIs) ---
le_fc <- life_expectancy(fc)

# --- Combine historical and forecast ---
le_hist[, type := "observed"]
le_fc[, type := "forecast"]
le_all <- rbind(le_hist, le_fc, fill = TRUE)

# --- Plot (example using ggplot2 directly) ---
library(ggplot2)
ggplot(le_all, aes(x = year, y = ex, colour = group)) +
  geom_line(data = le_all[type == "observed"]) +
  geom_line(data = le_all[type == "forecast"], linetype = "dashed") +
  geom_ribbon(data = le_all[type == "forecast"],
              aes(ymin = ex_lower, ymax = ex_upper, fill = group),
              alpha = 0.15, colour = NA) +
  labs(x = "Year", y = "Life expectancy at age 30",
       title = "Historical and forecast life expectancy") +
  theme_minimal()

# --- Disaggregated life expectancy ---
fc_1yr <- disaggregate(fc, target_ages = 30:100,
                       method = "natural_gompertz")
le_1yr <- life_expectancy(fc_1yr)
print(le_1yr)
```

## Summary

| Function                                                                                | Input                          | Output                                  | Key use                             |
|-----------------------------------------------------------------------------------------|--------------------------------|-----------------------------------------|-------------------------------------|
| [`life_table()`](https://kalleEF.github.io/mortensor/reference/life_table.md)           | mx vector + ages               | data.frame with qx, lx, dx, Lx, Tx, ex  | Full abridged life table            |
| [`life_expectancy()`](https://kalleEF.github.io/mortensor/reference/life_expectancy.md) | `mort_data` or `mort_forecast` | data.table with year, group, ex (+ PIs) | Summary measure across years/groups |
| [`compute_dx()`](https://kalleEF.github.io/mortensor/reference/compute_dx.md)           | mx vector                      | Death distribution (sums to ~1)         | Input to Tucker3-CoDa               |
| [`dx_to_mx()`](https://kalleEF.github.io/mortensor/reference/dx_to_mx.md)               | dx vector                      | Mortality rates                         | Inverse of compute_dx               |

The life-table machinery integrates seamlessly with the forecast
pipeline: `forecast() → add_pi() → disaggregate() → life_expectancy()`
produces group-specific life expectancy projections with prediction
intervals.
