# mortensor

**Cause-Specific Mortality Forecasting with Tensor Decomposition**

`mortensor` is a modular R package for forecasting age-specific mortality rates
across population strata (e.g. sex × deprivation quintile) and causes of
death. It combines Tucker tensor decomposition with compositional data analysis,
ensemble time-series methods, and rigorous uncertainty quantification.

---

## Why mortensor?

Mortality data are inherently multi-way: rates vary by **age**, **period**,
and **population group** simultaneously. Univariate methods (Lee-Carter, ARIMA
per series) ignore this structure. `mortensor` models the full
age × period × group array jointly via Tucker3 decomposition, extracting
low-rank patterns that generalise better out-of-sample.

Key design goals:

- **Stratified by default** — designed for sex × deprivation (or similar)
  cross-classifications, not just aggregate data
- **Ensemble-first** — automatic model comparison and inverse-MASE weighting
  over Tucker3, Lee-Carter, FDM, StMoMo, ARIMA, ETS, and RWD
- **Coherent multi-cause** — WLS reconciliation and Tucker4-ILR ensure
  disease-specific forecasts sum to all-cause
- **Honest uncertainty** — residual bootstrap and GAM-Bayes posterior
  propagation produce well-calibrated 95% prediction intervals
- **Actuarially complete** — disaggregates 5-year groups to single-year ages
  and computes full abridged life tables

---

## Installation

```r
# From GitHub
remotes::install_github("kalleEF/mortensor")

# From a local clone
devtools::install("path/to/mortensor")
```

### Optional backends

```r
install.packages("demography")  # Lee-Carter and FDM models
install.packages("StMoMo")      # Stochastic mortality models
install.packages("nnls")        # NNLS ensemble stacking
install.packages("readxl")      # ONS Excel import
```

---

## Quick start

```r
library(mortensor)

# Load the example dataset (mort_data object)
data(example_mortality)

# 1. Smooth — 2D GAM P-spline, exclude COVID years
sm <- smooth_mortality(example_mortality, method = "gam",
                       covid_years = c(2020L, 2021L))

# 2. Calibrate — compare Tucker3, Lee-Carter, FDM, ETS, ARIMA, RWD;
#    build inverse-MASE ensemble from top 3
cal <- calibrate(sm, holdout_frac = 0.2, ensemble = TRUE)
print(cal)

# 3. Forecast
fc <- forecast(cal, h = 10)

# 4. Add 95% prediction intervals (bootstrap)
fc <- add_pi(fc, n_sim = 500, method = "bootstrap")

# 5. Disaggregate to single-year ages (30–100)
fc_1yr <- disaggregate(fc, target_ages = 30:100,
                       method = "natural_gompertz")

# 6. Life expectancy (with PI bounds)
le <- life_expectancy(fc)
print(le)

# 7. Export to CSV
export(fc, dir = "results/")
```

---

## Pipeline overview

```
mort_data ──► smooth_mortality() ──► calibrate() ──► forecast()
                                                          │
                                                     add_pi()
                                                          │
                                                  disaggregate()
                                                          │
                                              life_expectancy() / export()
```

Every step returns an S3 object with `print()`, `summary()`, `plot()`, and
`export()` methods.

---

## Core features

### Tensor decomposition (Tucker3)

Mortality rates $m_{x,t,g}$ are arranged in a 3D array (age × year × group)
and decomposed as:

$$\mathcal{X} \approx \mathcal{G} \times_1 \mathbf{A} \times_2 \mathbf{B} \times_3 \mathbf{C}$$

The temporal mode $\mathbf{B}$ is the only component that changes over time
— its columns are forecast individually using ARIMA, ETS, or random walk with
drift. Two transform strategies are provided:

- **Tucker3-CoDa**: CLR transform on death *shares* (compositional, bounded)
- **Tucker3-Direct**: log-transform on mortality *rates* (simpler, often
  competitive for all-cause)

Ranks $(R_A, R_B, R_C)$ and the time-series method are selected by
rolling-origin cross-validation or BIC/AIC.

### Model ensemble

`calibrate()` evaluates all available models on a holdout period and builds
an ensemble via **inverse-MASE weighting**:

| Model | Backend |
|-------|---------|
| Tucker3-CoDa | rTensor + compositions |
| Tucker3-Direct | rTensor |
| Lee-Carter | demography |
| FDM | demography |
| StMoMo APC | StMoMo |
| ARIMA per series | forecast |
| ETS (damped) | forecast |
| Random walk + drift | — |

### Prediction intervals

Two methods are available via `add_pi()`:

| Method | Description | Speed |
|--------|-------------|-------|
| `"bootstrap"` | Residual resampling + ARIMA/ETS simulation | Fast (seconds) |
| `"gam_bayes"` | GAM posterior draws → full pipeline re-forecast | Slower (minutes) |

Both produce `mx_lower` / `mx_upper` / `mx_median` arrays and a full
`sims` array for custom quantile extraction.

### Age disaggregation

Three methods disaggregate abridged (5-year) age groups to single years:

| Method | Interpolation | Old-age extrapolation |
|--------|--------------|----------------------|
| `"hyman_kannisto"` | Monotone Hermite spline | Kannisto logistic |
| `"natural_gompertz"` | Natural cubic spline | Gompertz log-linear |
| `"gam"` | 2D GAM (age × year) | Blended parametric tail |

All methods enforce the group-mean constraint (mean of single-year rates
matches the original 5-year grouped rate).

### Multi-cause coherence

For cause-of-death forecasting, two coherence strategies are available:

**WLS reconciliation** (`reconcile()`): Forecast each cause independently,
then adjust via weighted least squares so that diseases sum to all-cause.
Works with any per-cause model.

**Tucker4-ILR** (`run_coherent_forecast()`): Model cause *shares* jointly in
a 4D tensor via ILR-transformed simplex coordinates. Coherence is guaranteed
by construction.

### Scoring rules

`mortensor` implements all major proper scoring rules for forecast evaluation:

| Metric | Function | Type |
|--------|----------|------|
| Combined MASE | `compute_mase_combined()` | Point |
| Cell MASE | `compute_mase_array()` | Point |
| Trend MASE | `compute_mase_group_trend()` | Point |
| RMSE (log) | `compute_rmse_log()` | Point |
| CRPS (empirical) | `compute_crps_empirical()` | Probabilistic |
| CRPS (Gaussian) | `compute_crps_gaussian()` | Probabilistic |
| Interval score | `compute_interval_score_array()` | Probabilistic |
| WIS | `compute_wis()` | Probabilistic |
| Log-score | `compute_log_score_gaussian()` | Probabilistic |

---

## S3 classes

| Class | Constructor | Key methods |
|-------|-------------|-------------|
| `mort_data` | `mort_data()`, `read_mortality()` | `print`, `summary`, `plot`, `[`, `as.data.table`, `life_expectancy` |
| `mort_calibration` | `calibrate()` | `print`, `summary`, `plot`, `forecast`, `export` |
| `mort_forecast` | `forecast()` | `print`, `summary`, `plot`, `coef`, `fitted`, `residuals`, `export`, `life_expectancy` |
| `mort_config` | `mort_config()` | `print` |
| `mort_shared_config` | `mort_shared_config()` | `print` |

---

## Batch pipelines

For end-to-end processing without interactive inspection:

```r
# Single cause
config <- mort_config(cause_label = "All causes", forecast_h = 10,
                      n_sim = 500, disagg_method = "natural_gompertz")
result <- run_forecast(data, config)

# Multi-cause with WLS reconciliation
shared <- mort_shared_config(output_dir = "results/", forecast_h = 10)
result <- run_reconciled_forecast(disease_configs, all_cause_config, shared)

# Multi-cause Tucker4-ILR (coherent by construction)
result <- run_coherent_forecast(disease_configs, all_cause_config, shared)
```

---

## Vignettes

| Vignette | Topic |
|----------|-------|
| [Getting Started](vignettes/getting-started.Rmd) | Installation, data loading, full pipeline walkthrough |
| [Tensor Decomposition](vignettes/tensor-decomposition.Rmd) | Tucker3 CoDa/Direct, Tucker4-ILR, B-mode forecasting |
| [Model Calibration](vignettes/model-calibration-and-selection.Rmd) | CV/BIC rank selection, ensemble building, holdout comparison |
| [Prediction Intervals](vignettes/prediction-intervals.Rmd) | Bootstrap, GAM-Bayes, custom quantiles |
| [Age Disaggregation](vignettes/age-disaggregation.Rmd) | Hyman+Kannisto, Natural+Gompertz, GAM, quality checks |
| [Life Tables](vignettes/life-tables.Rmd) | Abridged life tables, life expectancy, dx/mx conversions |
| [Multi-Cause Forecasting](vignettes/multi-cause-forecasting.Rmd) | WLS reconciliation, Tucker4-ILR, batch pipelines |
| [Scoring & Evaluation](vignettes/scoring-and-evaluation.Rmd) | MASE, CRPS, interval score, WIS, model comparison |
| [Visualization Gallery](vignettes/visualization-gallery.Rmd) | All 18+ plot types with examples |

---

## Dependencies

**Imports** (installed automatically): `data.table`, `rTensor`, `compositions`,
`forecast`, `mgcv`, `ggplot2`, `stats`

**Suggests** (optional): `demography`, `StMoMo`, `MASS`, `nnls`, `readxl`,
`knitr`, `rmarkdown`, `testthat`

---

## Citation

If you use `mortensor` in published work, please cite:

> [Author(s)]. (2024). *mortensor: Cause-Specific Mortality Forecasting with
> Tensor Decomposition*. R package version 0.1.0.

---

## License

MIT © [Author(s)]
