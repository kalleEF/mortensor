# Model Calibration, Selection, and Ensembles

## Overview

Choosing the right model — and its hyperparameters — is critical for
mortality forecasting. **mortensor** provides a unified calibration
framework that:

1.  Searches over Tucker3 ranks $\left( R_{A},R_{B},R_{C} \right)$ and
    time-series methods
2.  Calibrates all competing models (Lee-Carter, FDM, StMoMo, ARIMA,
    ETS, RWD)
3.  Evaluates on a held-out period using proper scoring rules
4.  Builds inverse-MASE or NNLS ensemble weights
5.  Returns a `mort_calibration` object ready for
    [`forecast()`](https://kalleEF.github.io/mortensor/reference/forecast.md)

## Available models

| Model name       | Function                                                                                                | Description                             |
|------------------|---------------------------------------------------------------------------------------------------------|-----------------------------------------|
| `tucker3_coda`   | [`tucker3_coda_forecast()`](https://kalleEF.github.io/mortensor/reference/tucker3_coda_forecast.md)     | Tucker3 on CLR-transformed death shares |
| `tucker3_direct` | [`tucker3_direct_forecast()`](https://kalleEF.github.io/mortensor/reference/tucker3_direct_forecast.md) | Tucker3 on log(mx) directly             |
| `lee_carter`     | [`lca_forecast()`](https://kalleEF.github.io/mortensor/reference/lca_forecast.md)                       | Lee-Carter (via `demography`)           |
| `fdm`            | [`fdm_forecast()`](https://kalleEF.github.io/mortensor/reference/fdm_forecast.md)                       | Functional demographic model            |
| `stmomo`         | [`stmomo_forecast()`](https://kalleEF.github.io/mortensor/reference/stmomo_forecast.md)                 | StMoMo APC/CBD (via `StMoMo`)           |
| `arima`          | [`arima_forecast()`](https://kalleEF.github.io/mortensor/reference/arima_forecast.md)                   | Per-series auto.arima on log(mx)        |
| `ets`            | [`ets_forecast()`](https://kalleEF.github.io/mortensor/reference/ets_forecast.md)                       | Per-series ETS (damped) on log(mx)      |
| `rwd`            | [`rwd_forecast()`](https://kalleEF.github.io/mortensor/reference/rwd_forecast.md)                       | Per-series random walk with drift       |

Models requiring optional packages (`demography`, `StMoMo`) are silently
skipped if the package is not installed.

## The `calibrate()` function

The main entry point for model selection is
[`calibrate()`](https://kalleEF.github.io/mortensor/reference/calibrate.md):

``` r
library(mortensor)
data(example_mortality)

# Smooth first
sm <- smooth_mortality(example_mortality, method = "gam",
                       covid_years = c(2020L, 2021L))

# Calibrate with defaults
cal <- calibrate(sm)
print(cal)
summary(cal)
```

### Key parameters

``` r
cal <- calibrate(
  data           = sm,                  # mort_data (smoothed recommended)
  models         = c("tucker3_coda",    # which models to evaluate
                     "tucker3_direct",
                     "lee_carter", "fdm",
                     "ets", "arima", "rwd"),
  holdout_years  = NULL,                # explicit holdout years (overrides frac)
  holdout_frac   = 0.2,                # fraction of years for holdout
  ensemble       = TRUE,               # build ensemble from top models
  ensemble_top_k = 3L,                 # how many models in ensemble
  config         = mort_config()        # optional configuration object
)
```

### What happens inside `calibrate()`

1.  **Train/test split**: The last `holdout_frac` of years (or
    `holdout_years`) are withheld for evaluation.

2.  **Smoothing**: If `data$smoothed` is `FALSE`,
    [`smooth_mortality()`](https://kalleEF.github.io/mortensor/reference/smooth_mortality.md)
    is called automatically.

3.  **Tucker3 calibration**:
    [`calibrate_tucker3()`](https://kalleEF.github.io/mortensor/reference/calibrate_tucker3.md)
    searches over rank combinations for both CoDa and Direct variants,
    selecting the best
    $\left( R_{A},R_{B},R_{C},\text{ts\_method} \right)$.

4.  **Other model calibration**:
    [`calibrate_all_models()`](https://kalleEF.github.io/mortensor/reference/calibrate_all_models.md)
    evaluates Lee-Carter, FDM, StMoMo, RWD, ARIMA, and ETS via
    expanding-window CV.

5.  **Holdout comparison**:
    [`compare_models()`](https://kalleEF.github.io/mortensor/reference/compare_models.md)
    runs every calibrated model on the training set, forecasts into the
    holdout period, and scores with combined MASE (0.7 cell-level + 0.3
    group-trend).

6.  **Ensemble construction**: The top-$k$ models are combined using
    inverse-MASE weights (or optionally NNLS stacking).

## Tucker3 rank selection in detail

### Cross-validation (recommended)

Rolling-origin CV (`rank_selection = "cv"`) trains on expanding windows
and evaluates on the next $h$ years:

``` r
calib_coda <- calibrate_tucker3(
  mx_array       = sm$mx,
  years          = sm$years,
  ages           = sm$ages,
  method         = "coda",
  rank_selection = "cv",
  h_val          = 3,           # 3-year validation windows
  min_train      = 12           # at least 12 training years
)

# Result contains the optimal configuration
str(calib_coda)
# $R_A       : 3
# $R_B       : 2
# $R_C       : 2
# $ts_method : "ets"
# $cv_score   : 0.847
```

The CV procedure:

1.  For each rank combination, fit Tucker3 on years
    $\left\lbrack 1,T_{0} \right\rbrack$, forecast $h$ steps, score
    against held-out truth using the chosen `cv_metric`.
2.  Roll forward: train on $\left\lbrack 1,T_{0} + 1 \right\rbrack$,
    forecast, score.
3.  Average scores across origins.
4.  Select ranks with lowest average score.
5.  Given best ranks, compare `ts_method` $\in$`{arima, ets, rwd}`.

By default, the scoring metric is cell-level MASE
(`cv_metric = "mase"`). You can also use `"mase_combined"` (0.7 × cell +
0.3 × trend) or `"mase_trend"` (trend only). See the [CV metric
selection](#cv-metric-selection) section below.

### Information criteria (fast)

For quick exploration, BIC or AIC can select ranks without refitting
time-series models:

``` r
calib_bic <- calibrate_tucker3(
  mx_array       = sm$mx,
  years          = sm$years,
  ages           = sm$ages,
  method         = "coda",
  rank_selection = "bic"        # or "aic"
)
```

BIC penalises model complexity more than AIC, favouring lower ranks.
This is fast (seconds vs minutes for CV) but does not account for
forecast performance.

### Direct comparison of selection methods

``` r
# Both methods
calib_cv  <- calibrate_tucker3(sm$mx, sm$years, sm$ages,
                                method = "coda", rank_selection = "cv")
calib_bic <- calibrate_tucker3(sm$mx, sm$years, sm$ages,
                                method = "coda", rank_selection = "bic")

cat(sprintf("CV  ranks: (%d, %d, %d) — MASE: %.3f\n",
            calib_cv$R_A, calib_cv$R_B, calib_cv$R_C, calib_cv$cv_score))
cat(sprintf("BIC ranks: (%d, %d, %d) — BIC: %.1f\n",
            calib_bic$R_A, calib_bic$R_B, calib_bic$R_C, calib_bic$ic_value))
```

## Calibrating non-tensor models

[`calibrate_all_models()`](https://kalleEF.github.io/mortensor/reference/calibrate_all_models.md)
handles the remaining model zoo:

``` r
calib_all <- calibrate_all_models(
  mx_array  = sm$mx,
  pop_array = sm$pop,        # population exposure (optional for some models)
  years     = sm$years,
  ages      = sm$ages,
  groups    = sm$groups,
  h_val     = 4,             # 4-year validation windows
  n_origins = 4,             # 4 rolling origins
  width     = 5              # 5-year age groups
)

# Inspect per-model calibration results
names(calib_all)
# [1] "lc" "fdm" "stmomo" "rwd" "arima" "ets"

# Each entry has model-specific parameters + cv_score
calib_all$lc
# $adjust     : "dt"
# $order      : 6
# $method     : "arima"
# $cv_score    : 0.912

calib_all$ets
# $damped     : TRUE
# $cv_score    : 1.034
```

## CV metric selection

By default, cross-validation selects hyperparameters using cell-level
MASE (`cv_metric = "mase"`). The `cv_metric` parameter lets you choose
an alternative scoring function applied to **all** models during
calibration:

| `cv_metric`       | Function                                                                                                  | Description                   |
|-------------------|-----------------------------------------------------------------------------------------------------------|-------------------------------|
| `"mase"`          | [`compute_mase_array()`](https://kalleEF.github.io/mortensor/reference/compute_mase_array.md)             | Cell-level MASE (default)     |
| `"mase_combined"` | [`compute_mase_combined()`](https://kalleEF.github.io/mortensor/reference/compute_mase_combined.md)       | 0.7 × cell + 0.3 × trend MASE |
| `"mase_trend"`    | [`compute_mase_group_trend()`](https://kalleEF.github.io/mortensor/reference/compute_mase_group_trend.md) | Group-trend MASE only         |

Using `"mase_combined"` aligns calibration with the metric used for
final model ranking in
[`compare_models()`](https://kalleEF.github.io/mortensor/reference/compare_models.md),
which can lead to more consistent selections:

``` r
# Calibrate Tucker3 using combined MASE
calib_combined <- calibrate_tucker3(
  mx_array  = sm$mx,
  years     = sm$years,
  ages      = sm$ages,
  method    = "coda",
  cv_metric = "mase_combined"
)

# The top-level calibrate() passes cv_metric to all models
cal <- calibrate(sm, cv_metric = "mase_combined")
```

## Model comparison and scoring

The holdout evaluation uses
[`compare_models()`](https://kalleEF.github.io/mortensor/reference/compare_models.md)
internally. The comparison data.table is stored in the calibration
object:

``` r
cal <- calibrate(sm, holdout_frac = 0.2, ensemble = TRUE)

# View the comparison table
print(cal$comparison)
#              model      mase mase_trend mase_combined
# 1:    Tucker3-CoDa  0.847     0.721        0.809
# 2:  Tucker3-Direct  0.891     0.756        0.851
# 3:      Lee-Carter  0.923     0.812        0.890
# 4:             FDM  0.967     0.845        0.930
# 5:      ETS-Damped  1.034     0.923        0.991
# 6: ARIMA-perSeries  1.089     0.967        1.052
# 7:             RWD  1.123     1.001        1.087

# The combined MASE is 0.7 * cell + 0.3 * trend
```

The **combined MASE** balances cell-level accuracy (capturing age-group
patterns) with group-trend accuracy (capturing overall decline rates).
This prevents selecting a model that is accurate on average but misses
important trends.

### Visualise model comparison

``` r
plot(cal)  # bar chart of MASE by model
```

## Ensemble construction

### Inverse-MASE weighting (default)

The top-$k$ models are weighted inversely proportional to their MASE:

$$w_{i} = \frac{1/\text{MASE}_{i}}{\sum\limits_{j = 1}^{k}1/\text{MASE}_{j}}$$

``` r
# Inspect ensemble
cat("Ensemble members:", cal$ensemble_names, "\n")
cat("Weights:", round(cal$ensemble_weights, 3), "\n")

# Example output:
# Ensemble members: Tucker3-CoDa Tucker3-Direct Lee-Carter
# Weights: 0.412 0.321 0.267
```

### Controlling ensemble membership

``` r
# Only Tucker3 models in ensemble
cal_t3 <- calibrate(sm,
                    models = c("tucker3_coda", "tucker3_direct"),
                    ensemble = TRUE, ensemble_top_k = 2)

# No ensemble — use single best model
cal_single <- calibrate(sm, ensemble = FALSE)
cat("Best single model:", cal_single$best_model, "\n")

# Larger ensemble
cal_big <- calibrate(sm, ensemble = TRUE, ensemble_top_k = 5)
```

### NNLS stacking (advanced)

If the `nnls` package is installed, non-negative least-squares can
determine weights by minimising holdout forecast error directly:

``` r
# NNLS is used automatically if nnls is installed and ensemble = TRUE
# The method that produces lower holdout MASE (inverse-MASE or NNLS) is kept
```

## Dispatching forecasts

The
[`forecast()`](https://kalleEF.github.io/mortensor/reference/forecast.md)
method on `mort_calibration` uses
[`dispatch_forecast()`](https://kalleEF.github.io/mortensor/reference/dispatch_forecast.md)
to route to the correct model function:

``` r
# forecast() internally calls dispatch_forecast() for each ensemble member,
# then combines using the stored weights
fc <- forecast(cal, h = 10)

# For a single-model calibration, it dispatches directly
fc_single <- forecast(cal_single, h = 10)
cat("Model used:", fc_single$model, "\n")
```

The dispatch table:

| Calibration model | Dispatches to                                                                                           |
|-------------------|---------------------------------------------------------------------------------------------------------|
| `tucker3_coda`    | [`tucker3_coda_forecast()`](https://kalleEF.github.io/mortensor/reference/tucker3_coda_forecast.md)     |
| `tucker3_direct`  | [`tucker3_direct_forecast()`](https://kalleEF.github.io/mortensor/reference/tucker3_direct_forecast.md) |
| `lee_carter`      | [`lca_forecast()`](https://kalleEF.github.io/mortensor/reference/lca_forecast.md)                       |
| `fdm`             | [`fdm_forecast()`](https://kalleEF.github.io/mortensor/reference/fdm_forecast.md)                       |
| `stmomo`          | [`stmomo_forecast()`](https://kalleEF.github.io/mortensor/reference/stmomo_forecast.md)                 |
| `arima`           | [`arima_forecast()`](https://kalleEF.github.io/mortensor/reference/arima_forecast.md)                   |
| `ets`             | [`ets_forecast()`](https://kalleEF.github.io/mortensor/reference/ets_forecast.md)                       |
| `rwd`             | [`rwd_forecast()`](https://kalleEF.github.io/mortensor/reference/rwd_forecast.md)                       |

## Expanding-window CV for any model

You can run expanding-window cross-validation for any model
independently:

``` r
cv_result <- cv_expanding_model(
  mx_array  = sm$mx,
  pop_array = sm$pop,
  years     = sm$years,
  ages      = sm$ages,
  groups    = sm$groups,
  model     = "lee_carter",
  h_val     = 4,
  n_origins = 4,
  width     = 5
)

cat("Lee-Carter CV score:", round(cv_result$cv_score, 3), "\n")
```

## Using `mort_config` for calibration settings

The `mort_config` object centralises pipeline settings:

``` r
cfg <- mort_config(
  forecast_h       = 10L,
  n_sim            = 500L,
  smoothing_method = "gam",
  covid_years      = c(2020L, 2021L, 2022L),
  holdout_train    = 0.8,     # 80% train / 20% holdout
  holdout_test     = 0.2
)

cal <- calibrate(sm, config = cfg)
```

## Complete example: from data to calibrated ensemble

``` r
library(mortensor)
data(example_mortality)

# 1. Smooth
sm <- smooth_mortality(example_mortality, method = "gam",
                       covid_years = c(2020L, 2021L))

# 2. Calibrate all models with 20% holdout, top-3 ensemble
cal <- calibrate(sm, holdout_frac = 0.2,
                 ensemble = TRUE, ensemble_top_k = 3)

# 3. Inspect results
summary(cal)
plot(cal)

# 4. Check ensemble
cat("Best model:", cal$best_model, "\n")
cat("Ensemble:", paste(cal$ensemble_names, collapse = ", "), "\n")
cat("Weights:", paste(round(cal$ensemble_weights, 3), collapse = ", "), "\n")

# 5. Forecast
fc <- forecast(cal, h = 10)
summary(fc)

# 6. Export calibration details
export(cal, dir = "calibration_results/")
```

## Summary

- [`calibrate()`](https://kalleEF.github.io/mortensor/reference/calibrate.md)
  is the one-stop function for model selection and ensemble building
- Tucker3 ranks are selected via CV (accurate) or BIC/AIC (fast)
- All models are compared on a holdout set using combined MASE
- Inverse-MASE weighting produces robust ensembles
- The `mort_calibration` object plugs directly into
  [`forecast()`](https://kalleEF.github.io/mortensor/reference/forecast.md),
  [`add_pi()`](https://kalleEF.github.io/mortensor/reference/add_pi.md),
  and
  [`export()`](https://kalleEF.github.io/mortensor/reference/export.md)
