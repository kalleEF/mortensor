# mortensor production review: handoff notes

Working notes from a first review session (2026-10-12). That session could not
install `rTensor`, `compositions`, `demography` or `StMoMo` because CRAN was
blocked by the environment's network policy, so parts of the review are based
on reading the code only. Continue from here in a session with CRAN access.

## Goal

The owner wants a full evaluation of the package's current state and a
concrete roadmap to turn it into a fully functional, production-grade
mortality forecasting R package. Known pain points from the owner: many
functions do not really work and the vignettes are largely unusable.

Expected deliverable: a written review (current state, verified defects,
methodological issues, prioritised roadmap) plus a short chat summary.

## Environment setup that worked

R is not preinstalled in the cloud container. Ubuntu 24.04 apt provides R 4.3.3
and most dependencies:

```sh
apt-get update
DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  r-base-core r-base-dev pandoc qpdf libcurl4-openssl-dev libssl-dev \
  libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
  libfreetype-dev libpng-dev libtiff-dev libjpeg-dev \
  r-cran-data.table r-cran-forecast r-cran-mgcv r-cran-ggplot2 r-cran-rlang \
  r-cran-readxl r-cran-testthat r-cran-knitr r-cran-rmarkdown r-cran-nnls \
  r-cran-mass r-cran-roxygen2 r-cran-devtools r-cran-rcmdcheck r-cran-lintr \
  r-cran-codetools r-cran-generics
```

Not available via apt; install from CRAN (needs `cloud.r-project.org` allowed):

```r
install.packages(c("rTensor", "compositions", "demography", "StMoMo"),
                 repos = "https://cloud.r-project.org")
```

`packagemanager.posit.co` was also blocked (403). R 4.3.3 is useful for
reproducing the `%||%` defect below; a test on R >= 4.4 is also worthwhile.

## Confirmed by running code (R 4.3.3)

1. **`%||%` is never defined.** It is used in 23 functions (for example
   `R/methods.R:89`, `R/workflow.R:289`, `R/mort_data.R:98`) and only exists in
   base R from 4.4.0. DESCRIPTION declares `Depends: R (>= 4.1.0)`. On R 4.3.3
   `print(example_mortality)` fails with `could not find function "%||%"`.
2. **`forecast()` and `accuracy()` generics are masked.** `R/workflow.R:28` and
   `R/metrics.R:40` define new generics with the same names as
   `forecast::forecast` and `forecast::accuracy` (forecast is in Imports).
   After `library(forecast)`, `forecast(<mort_calibration>)` dispatches to
   `forecast.default` and errors ("non-numeric argument to binary operator").
   Fix: import and re-export the forecast generics and register methods on
   them.
3. Static check: every exported name in NAMESPACE is defined; no duplicate
   top-level definitions.

## Found by reading code (verify by running)

### Model selection and evaluation

- **Holdout leakage.** `calibrate()` smooths all years including the holdout
  (`R/workflow.R:166-172`), then tunes Tucker3 ranks and all other models on
  the full series (`R/workflow.R:196-224`). `cv_rolling_origin()` sets
  `last_origin <- n_yr - h_val` (`R/calibration.R:22`), so its last fold
  evaluates exactly on the holdout years. The holdout comparison is therefore
  optimistic and favours heavily tuned models. Same pattern in
  `run_forecast()` (`R/pipeline.R:111-135`).
- **Ensemble scored on the data it was fitted on.** `compare_models()` picks
  ensemble size and NNLS stacking weights on the holdout and reports the same
  holdout score (`R/comparison.R:451-528`); `calibrate()` then prefers the
  ensemble if it is within 5% of the best single model. NNLS weights are
  estimated on log(mx) but applied to mx.
- **`models` argument of `calibrate()` is ignored by the comparison.**
  `compare_models()` always runs all eight models (`R/comparison.R:344-346`);
  user names (`"lee_carter"`) differ from internal names (`"Lee-Carter"`).
- **Incomplete models can win.** Models with NA cells are dropped per series
  inside MASE and can still be selected as best in `calibrate()`
  (`R/workflow.R:263-267`).
- **Grid size.** Default CV grid is 5 x 7 x 3 ranks x 3 ts methods = 315
  Tucker fits per origin and per variant (`R/calibration.R:161-164`). Time it.
- **Hyperparameter grids partly invalid.** ETS candidates `"AAdN"`, `"ANN"`
  with `damped = TRUE`, and `"MMN"` on negative log rates all error and fall
  back silently to `auto.arima` (`R/calibration.R:394`,
  `R/models_timeseries.R:93`). RWD drift windows `0` and `length(years)` are
  identical (`R/calibration.R:362`).
- **Single-group data breaks Tucker3.** `2:min(4, n_g - 1)` yields descending
  sequences for `n_g <= 2`, and `R_C = 0` for `n_g = 1`.
- **`run_forecast()` defaults disable all CV.** `holdout_train`/`holdout_test`
  default to NULL, so `h_val = 0` (`R/pipeline.R:111`) and every CV score is
  Inf; defaults are used silently.

### Prediction intervals

- **No growth of forecast uncertainty with horizon** for RWD/ARIMA/ETS
  (`R/dispatch.R:256-298`) and Lee-Carter/FDM (`R/bootstrap.R:159-162`):
  i.i.d. residuals are added to the point forecast. Lee-Carter kt uncertainty
  (the main source of LC PI width) is ignored.
- **`add_pi()` on an ensemble** simulates only the first ensemble component
  (`R/workflow.R:330-339`, `447-452`).
- **`ensemble_bootstrap()` averages draws across models** instead of sampling
  from the mixture, which shrinks spread (`R/bootstrap.R:230-251`).
- **Re-centring and capping.** PIs are shifted to the point forecast and
  capped at +/-1.5 on the log scale (`R/workflow.R:508-522`,
  `R/pipeline.R:330-343`); `mx_median` is overwritten with the point forecast.
- **`gam_bayes`** propagates only smoothing uncertainty, not forecast
  innovation uncertainty (`R/bootstrap.R:277-341`).
- Arbitrary fallbacks: `rnorm(h, 0, 0.1)` when simulation fails.

### Silent failures

- Failed model cells are set to 0 inside ensemble averages, biasing the
  ensemble downward (`R/workflow.R:320`, `R/comparison.R:470`,
  `R/bootstrap.R:238`, `R/bootstrap.R:323`).
- Lee-Carter/FDM/StMoMo failures per group return NA with no warning
  (`R/models_demographic.R`, `R/models_stmomo.R`).
- `build_array()` sorts age labels as text and fills cells with 1e-12 when the
  number of ages does not match (`R/array_utils.R:32-43`); the floor applies
  only if the column name contains "mx".

### Demographic methods

- **`life_expectancy()` assumes 5-year groups** (`width = 5` default,
  `R/life_table.R:245`, `291`), so e_x of a disaggregated single-year forecast
  is wrong. Widths should be derived from `diff(ages)`.
- **e_x bounds** come from pointwise mx bounds (`R/life_table.R:320-327`)
  rather than from simulated paths.
- `life_table()` doc says Kannisto is the default open-age method; code
  defaults to reciprocal (`R/life_table.R:90-99`). Kannisto/Gompertz open
  interval integration starts survival at exp(-mu_0) (off by one,
  `R/life_table.R:190-197`); ages use group starts instead of midpoints.
- **CoDa loses the open age group.** `compute_dx()` ignores mx of the open
  interval (`R/life_table.R:25-27`); `dx_to_mx()` replaces it with a Gompertz
  extrapolation from the two previous groups (`R/life_table.R:67-73`). The
  90+ forecast of Tucker3-CoDa is not learned from 90+ data.
- **Disaggregation uses an unweighted group mean**
  (`R/disaggregation.R:80-96`, `128-145`); the correct constraint is
  exposure-weighted, and for the open interval the simple mean over 90-100
  biases old-age rates downward. 90+ open interval and midpoints are
  hardcoded (`R/disaggregation.R:60`, `79`). Consider the penalized composite
  link model (`ungroup` package) as a reference method; verify citations
  before citing.
- StMoMo is fitted with age index 1..n instead of ages
  (`R/models_stmomo.R:22`, `38`); CBD with logit link uses central exposures
  and returns q rather than m. Fitted values for StMoMo always use APC
  (`R/dispatch.R:92`).
- Smoothing is a Gaussian GAM on log(mx), ignoring exposure; a Poisson model
  with log-exposure offset is the usual choice when deaths and population are
  available.

### Multi-cause

- `run_reconciled_forecast()` and `run_coherent_forecast()` call
  `run_forecast(<mort_config>)` (`R/pipeline.R:685`, `910`), which only works
  with the ONS Excel input path; not usable with generic data.
- WLS reconciliation only closes the gap exactly if a cause named "residual"
  exists (`R/reconciliation.R:65`).
- `compute_cause_shares()` adds a remainder part even when the disease list
  already contains "residual" (`default_disease_order()`, `R/config.R:32`),
  creating a near-zero component (1e-6) in the ILR basis.
- `enforce_coherence()` hardcodes `sex == "Male"/"Female"` and `quintile`
  (`R/coherence.R:26-28`); it enforces female <= male, which is a constraint
  rather than demographic coherence in the Li-Lee sense.

### Architecture and package hygiene

- Two parallel implementations: `run_forecast()` (pipeline) and the S3
  workflow (`calibrate` -> `forecast` -> `add_pi` -> `disaggregate`) diverge
  in behaviour (ensemble PIs, coherence, output tables).
- `disaggregate(method = "gam")` stops (`R/workflow.R:566-569`);
  `plot(type = "disagg")` and single-year `export()` rely on
  `extras$mx_1yr_dt`, which only `run_forecast()` sets
  (`R/methods.R:176`, `1230`).
- `cat()` used for progress in ~170 places (not suppressible); plot methods
  call `print(p)` (`R/methods.R:137` etc.); exported `plot_*` helpers write
  PNGs via `ggsave()` as a side effect.
- Tests are thin and never run the main workflow end to end; the
  `build_array` test passes while the function returns an all-NA array
  (`tests/testthat/test-core.R:131-145`).
- No R CMD check CI workflow (only pkgdown).
- Committed files that should not be in the repo:
  `.claude/settings.local.json` (contains a local Windows path) and
  `tests/testthat/Rplots.pdf`.

## Still to do

1. Confirm CRAN access, install the four missing packages, install the
   package from this repo unmodified.
2. `R CMD build` and `R CMD check --as-cran`; record all errors, warnings,
   notes.
3. Run `testthat`.
4. Knit each vignette in `vignettes/`; record which chunks fail, why, and
   runtime.
5. Run the README quick start and the full workflow with all models on
   `example_mortality`; time `calibrate()` with `rank_selection = "cv"` and
   `"bic"`.
6. Empirically confirm the read-only findings above, especially leakage,
   PI coverage (simple backtest), the ensemble NA-to-0 bias, `add_pi()` on
   ensembles, and `life_expectancy()` on disaggregated forecasts.
7. Check `compute_wis()` normalisation (`R/metrics.R:429-448`, divides by the
   sum of weights) against the WIS definition in Bracher et al. (2021),
   verifying the reference before citing it.
8. Write the final review and prioritised roadmap.
