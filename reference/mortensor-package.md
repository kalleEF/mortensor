# mortensor: flexible mortality forecasting across data granularities

The package is organised in layers:

## Details

- **Data**: [`mt_data()`](mt_data.md) stores deaths, exposures and
  populations observed at any age and period resolution
  ([`age_groups()`](mt_axis.md), [`single_ages()`](mt_axis.md),
  [`single_years()`](mt_axis.md),
  [`parse_age_labels()`](parse_age_labels.md));
  [`as_mt_data()`](as_mt_data.md) and [`mt_from_hmd()`](mt_from_hmd.md)
  read Human Mortality Database data; [`mt_validate()`](mt_validate.md)
  checks plausibility; [`mt_example()`](mt_example.md) simulates data
  with known truth.

- **Harmonisation**: [`mt_regroup()`](mt_regroup.md) aggregates exactly,
  [`mt_ungroup()`](mt_ungroup.md) estimates single-year data with the
  penalized composite link model, [`mt_smooth()`](mt_smooth.md) smooths
  rates; [`mt_prep()`](mt_prep.md) bundles these steps so that backtests
  refit them on every training window.

- **Models**: [`rwd()`](rwd.md), [`lee_carter()`](lee_carter.md),
  [`p_splines()`](p_splines.md), [`cp_splines()`](cp_splines.md),
  [`arima_series()`](arima_series.md),
  [`ets_series()`](arima_series.md),
  [`hyndman_ullah()`](hyndman_ullah.md),
  [`coherent_fdm()`](coherent_fdm.md), [`stmomo()`](stmomo.md)
  (including the Plat model), [`stanmomo()`](stanmomo.md),
  [`li_lee()`](li_lee.md), [`coda()`](coda.md),
  [`bohk_rau()`](bohk_rau.md), [`nordpred_apc()`](nordpred_apc.md),
  [`apc_bayes()`](apc_bayes.md) (package bamp), the experimental
  [`tucker3()`](tucker3.md), [`coda_multiway()`](coda_multiway.md) and
  [`coda_causes()`](coda_causes.md), and [`ensemble()`](ensemble.md),
  fitted with [`mt_fit()`](mt_fit.md) (optionally tuned with
  [`tuned()`](mt_tune.md)) and forecast with
  [`mt_forecast()`](mt_forecast.md) as simulated paths
  ([mt_draws](mt_draws.md)). The article "Methods and sources" lists the
  published source of every method.

- **Hierarchies**: [`mt_hierarchy()`](mt_hierarchy.md) describes crossed
  or nested strata (for example sex x region x deprivation quintile with
  their totals), [`mt_aggregate_strata()`](mt_aggregate_strata.md) adds
  the aggregates, [`mt_reconcile()`](mt_reconcile.md) and
  [`reconciled()`](reconciled.md) make forecasts add up draw by draw,
  [`mt_order_strata()`](mt_order_strata.md) imposes orderings.

- **Wrappers for any model**: [`ensemble()`](ensemble.md),
  [`reconciled()`](reconciled.md), [`recalibrated()`](recalibrated.md)
  (spread fitted to recent forecast errors), [`anchored()`](anchored.md)
  (jump-off at the observed level), [`bootstrapped()`](bootstrapped.md)
  and [`tuned()`](mt_tune.md).

- **Evaluation**: [`mt_backtest()`](mt_backtest.md) and
  [`mt_score()`](mt_score.md), for rates, life expectancy and
  age-standardised rates.

- **The past**: [`mt_observed()`](mt_observed.md),
  [`mt_contrast()`](mt_contrast.md) and
  [`mt_inequality()`](mt_inequality.md).

- **Derived quantities and plots**: [`life_table()`](life_table.md),
  [`life_expectancy()`](life_expectancy.md),
  [`mt_deaths()`](mt_deaths.md), [`mt_plot_forecast()`](mt_plot.md) and
  the other `mt_plot_*()` functions.

- **Decomposition**: [`mt_decompose()`](mt_decompose.md) splits
  differences in life expectancy and other summaries by age, cause and
  subgroup composition, for observed data, forecasts or both.

## Caution

mortensor is experimental software written largely with AI coding
assistants ("vibe-coded"). It has unit tests that check its methods
against published results and reference implementations, but it has not
been peer reviewed and its interfaces still change. Use it at your own
risk and check results independently before using them in publications,
reports or decisions.

## See also

Useful links:

- <https://github.com/kalleEF/mortensor>

- Report bugs at <https://github.com/kalleEF/mortensor/issues>

## Author

**Maintainer**: Karl Emmert-Fees <karl.emmert-fees@tum.de> \[copyright
holder\]
