# Progress log

Status per work package of `dev/implementation-plan.md`. Branch for this
work: `claude/vigilant-albattani-2ssw1x` (plays the role of `dev-0.2`).

| WP | Status | Commit | Notes |
|---|---|---|---|
| P0.1, P0.2, P0.5 | done | 41058ba | 0.1.1: generics and `%\|\|%` fixed, GPL-3, baseline in `dev/baseline/` |
| P0.3 | done | this branch | Workflows: R-CMD-check (R 4.3, oldrel, release, devel; macOS, Windows), test-coverage, lint, validation, pkgdown |
| P0.4 | dropped | | 0.1.x vignettes not revived; the clean break replaces them |
| P0.6 | done | this branch | `CLAUDE.md`, this file |
| Q0.1 | done | this branch | `validation/` pipeline with quick and full suites |
| D1.1 | done | this branch | Skeleton, conditions, checks, name-clash test |
| D1.2 | done | this branch | `mt_axis`, `parse_age_labels()` |
| D1.3 | partial | this branch | Strata table; hierarchy sums checked in `mt_validate()`; reconciliation later (P5.4) |
| D1.4, D1.5 | done | this branch | `mt_data()`, `mt_add_*()`, grid, `mt_array()`, `mt_window()`, `mt_granularity()`, `mt_validate()` |
| D1.6 | done | this branch | `mt_regroup()`, `mt_operator()` |
| D1.7 | done | this branch | `mt_as_exposure()` (mid-year, January 1) |
| D1.8 | done | this branch | Native PCLM (1D, 2D), rate-based with stationary weights in the open interval |
| D1.9 | done | this branch | `mt_ungroup()`; E&W benchmark 0.015 / 0.024 (1D), 0.015 / 0.023 (2D) |
| D1.10 | done | this branch | `mt_smooth()` (Poisson P-splines, 1D/2D), `mt_prep()` with `step_ungroup()`, `step_smooth()`, `step_regroup()`, `step_grid()`, `step_shocks()` |
| D1.11 | done | this branch | `life_table()` (interval widths, reciprocal or Gompertz open interval), `life_expectancy()` |
| D1.12 | partial | this branch | `as_mt_data()` for data frames, StMoMoData, demogdata; `mt_example()`; ONS reader not ported |
| M2.1 | done | this branch | `new_mt_spec()`, `mt_fit()`, `mt_forecast()` (incl. predictive draws), per-stratum failures in diagnostics; contract test `expect_mt_model()` |
| M2.2 | done | this branch | `mt_draws`, `mt_summarise()`, `mt_deaths()` |
| M2.3 | done | this branch | `rwd()` (bootstrap or Gaussian), `arima_series()`, `ets_series()` |
| M2.4 | done | this branch | `lee_carter()` Poisson and SVD, jump-off, bootstrap; agrees with StMoMo (6e-5) and demography (1e-6 on fitted log rates) |
| M2.5 | done | this branch | `tucker3()` native HOOI, log and CoDa (open interval via auxiliary series), RW/ARIMA/ETS factors, residual noise; regression vs 0.1.0 not reproduced (different preprocessing; see notes) |
| M2.6 | done | this branch | `mt_backtest()` scores at the observed resolution; leakage test |
| M2.7 | done | this branch | `mt_score()`: CRPS, WIS (Bracher et al. 2021), coverage, PIT, MASE |
| M2.8 | done | this branch | `mt_tune()` (inner rolling origin) and `mt_tuned()` (tuning refitted inside each fit, so nested backtests do not leak) |
| C4.1 | done | this branch | `mt_ensemble()` (equal, numeric, inner-backtest CRPS weights) |
| B3.1-B3.4 | done | this branch | `apc_bayes()` + `mcmc_control()`: fit at the death resolution (`periods_per_agegroup`), explicit convergence check with refits, draws mapped to the grid with shape factors, `apc_effects()`; DIC-based selection not wrapped |
| C4.2 | done | this branch | `life_expectancy()`, `prob_dying()`, `asr()` with `standard_population("esp2013")`, `yll()`, `mt_deaths()`; all path by path |
| C4.3 | done | this branch | `mt_plot_forecast()`, `mt_plot_age()`, `mt_plot_life_expectancy()`, `mt_plot_rates()`, `mt_plot_ungroup()`, `mt_plot_coverage()` (tests check that plots build; no vdiffr snapshots) |
| C4.4 | done | this branch | `mt_export()` (summary or long draws, CSV) |
| C4.5 | done | this branch | Vignettes: getting started (simulated), granularity and ungrouping (HMD E&W, France), forecasting and backtesting (HMD France, E&W), Bayesian APC with bamp (HMD France; precomputed via `dev/precompute-vignettes.R`) |
| C4.7 | partial | this branch | Migration vignette (0.1 to 0.2 function map), pkgdown reference index; tagging a release is the owner's decision |
| D1.12 (HMD) | done | this branch | `mt_from_hmd()`, `mt_read_hmd()` (HMDHFDplus, credentials from `HMD_USER`/`HMD_PASSWORD`) |

| P5.1 | done | this branch | `hyndman_ullah()` (demography; draws via `simulate()`); 95% coverage 0.90 on the benchmark, 0.73 on France |
| P5.2 | done | this branch | `stmomo()` APC, CBD, RH, M6, LC with integer ages, initial exposures for logit links, `simulate()` draws; RH fails to converge at some origins (recorded) |
| P5.3 | done | this branch | `cause_shares()` (Tucker4 on ILR cause shares, any all-cause model), `mt_example_causes()` |
| P5.4 | done | this branch | `mt_reconcile()` draw by draw (MinT with draw covariance, WLS, OLS) |
| P5.5 | done | this branch | `mt_cohere()`, `li_lee()`, `product_ratio()` (Lee-Carter-type version) |
| P5.6 | done | this branch | `mt_ensemble(weights = "stack")` |
| P5.7 | partial | this branch | Whittaker and Gaussian GAM (mgcv) smoothing in `mt_smooth()`; GAM-Bayes smoothing uncertainty not ported |
| P5.8 | partial | this branch | Monotone (Hyman) splitting of deaths; natural-Gompertz and Hyman-Kannisto disaggregation not ported (PCLM is the documented method) |
| P5.9 | partial | this branch | Multi-cause vignette on `mt_example_causes()` (HMD has no causes); release 0.3.0 open |
| C4.6 | partial | this branch | 2D PCLM 87 s to 45 s (standard errors only for the final fit, Cholesky solves); no further work |

## Decisions taken during implementation

- 2026-10-12: Exported verbs use the `mt_` prefix; `mt_obs()` returns an
  observation table (the class name in the plan is not used as a class).
- 2026-10-12: The open age interval is expanded into single years up to
  omega inside the PCLM; for rates the exposure of the open interval is split
  with stationary-population weights from the current fit (three outer
  iterations). This replaces the artificial bin of `ungroup`.
- 2026-10-12: Default smoothing selection BIC on log10 lambda in [-2, 8],
  knots every 2 years, cubic B-splines, second-order differences
  (alternatives tested on E&W: third-order differences and wider knots gave
  the same accuracy within 0.003).
- 2026-10-12: 2D ungrouping takes about 90 s for 71 ages x 32 years; it is
  an option, 1D per year is the default.

- 2026-10-12: Backtests score predictive draws of observed log rates: smooth
  models (Lee-Carter, Tucker3) add extra-Poisson noise estimated from
  in-sample residuals and Poisson variation given the test-year exposures;
  models on observed log rates (random walk, ARIMA, ETS) are used as they
  are.
- 2026-10-12: Calibration on the simulation benchmark (6 replicates, h =
  1-10, 95% intervals): Lee-Carter 0.959 (min 0.924), Tucker3 0.947 (0.923),
  Tucker3-CoDa 0.943 (0.917), random walk 0.987 (over-covers: it treats iid
  noise as persistent). Tucker3 needs `residual_noise = TRUE` (default;
  0.87 without).
- 2026-10-12: England & Wales men 30-100, origins 1996-2001, h = 10:
  Lee-Carter CRPS 0.069 with 95% coverage 0.75; random walk 0.059 and 0.97.
  Lee-Carter fitted from 1981 with observed jump-off: CRPS 0.054, coverage
  0.74. Lee-Carter intervals are too narrow on these data (known limitation:
  only period-index uncertainty).
- 2026-10-12: The 0.1.0 regression baseline (`dev/baseline/`) used GAM
  smoothing and the old data structures; exact reproduction is not a goal.
  The Tucker3 core is tested by exact recovery of a low-rank array instead.

- 2026-10-12: bamp defaults changed to RW2 age, RW1 period and cohort: with
  RW2 period/cohort bamp's own forecasts exploded on French HMD data
  (ages 90-94, 15 years ahead, 90% interval 0.004-0.63 against 0.11-0.17
  with RW1). RW1 forecasts are biased when mortality keeps declining
  (backtest France 1996+10: CRPS 0.067 vs Lee-Carter 0.037, PIT mean 0.23).
- 2026-10-12: Splitting an open exposure interval (85+) into single years
  from its total alone gave errors of 27% at 85-99 and a factor of about 6
  at 100+ (France); `mt_ungroup()` now warns
  (`mt_warning_ungroup_open_counts`).
- 2026-10-12: lintr's `object_name_linter` is disabled (matrix notation such
  as `D`, `E`, `XtWX` is used throughout; lintr 3.1.1 failed with custom
  regexes); line length 120.

- 2026-10-12: `fdm()` renamed `hyndman_ullah()` (clash with `demography::fdm`).
- 2026-10-12: `cause_shares()` with noisy cause-specific rates over-covered
  (95% coverage 0.99) unless the rates are smoothed first
  (`step_smooth()`); with smoothing 0.94-0.97. Residual cause added only when
  the listed causes leave more than 0.1% unexplained (`residual = "auto"`).
- 2026-10-12: France 1986-1995 origins, h = 10 (CRPS, 95% coverage): Li-Lee
  0.058/0.95, product-ratio 0.059/0.91, Lee-Carter 0.061/0.95, StMoMo APC
  0.067/0.95, RH 0.067/0.97 (2 of 4 origins failed), FDM 0.069/0.73, CBD
  0.194/0.67 (CBD is designed for older ages; ages 30+ used here).

## Open questions

See `dev/implementation-plan.md` Section 17.
