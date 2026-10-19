#' mortensor: flexible mortality forecasting across data granularities
#'
#' The package is organised in layers:
#'
#' * **Data**: [mt_data()] stores deaths, exposures and populations observed
#'   at any age and period resolution ([age_groups()], [single_ages()],
#'   [single_years()], [parse_age_labels()]); [as_mt_data()] and
#'   [mt_from_hmd()] read Human Mortality Database data; [mt_validate()]
#'   checks plausibility; [mt_example()] simulates data with known truth.
#' * **Harmonisation**: [mt_regroup()] aggregates exactly, [mt_ungroup()]
#'   estimates single-year data with the penalized composite link model,
#'   [mt_smooth()] smooths rates; [mt_prep()] bundles these steps so that
#'   backtests refit them on every training window.
#' * **Models**: [rwd()], [lee_carter()], [p_splines()], [cp_splines()],
#'   [arima_series()], [ets_series()], [hyndman_ullah()], [coherent_fdm()],
#'   [stmomo()] (including the Plat model), [stanmomo()], [li_lee()],
#'   [coda()], [bohk_rau()], [nordpred_apc()], [apc_bayes()] (package bamp),
#'   the experimental [tucker3()] and [coda_causes()], and [ensemble()],
#'   fitted with [mt_fit()] (optionally tuned with [tuned()]) and forecast
#'   with [mt_forecast()] as simulated paths ([mt_draws]).
#'   [product_ratio()] and [cause_shares()] are deprecated. The article
#'   "Methods and sources" lists the published source of every method.
#' * **Hierarchies**: [mt_hierarchy()] describes crossed or nested strata
#'   (for example sex x region x deprivation quintile with their totals),
#'   [mt_aggregate_strata()] adds the aggregates, [mt_reconcile()] and
#'   [reconciled()] make forecasts add up draw by draw, [mt_order_strata()]
#'   imposes orderings.
#' * **Wrappers for any model**: [ensemble()], [reconciled()],
#'   [recalibrated()] (spread fitted to recent forecast errors), [anchored()]
#'   (jump-off at the observed level), [bootstrapped()] and [tuned()].
#' * **Evaluation**: [mt_backtest()] and [mt_score()], for rates, life
#'   expectancy and age-standardised rates.
#' * **The past**: [mt_observed()], [mt_contrast()] and [mt_inequality()].
#' * **Derived quantities and plots**: [life_table()], [life_expectancy()],
#'   [mt_deaths()], [mt_plot_forecast()] and the other `mt_plot_*()`
#'   functions.
#' * **Decomposition**: [mt_decompose()] splits differences in life
#'   expectancy and other summaries by age, cause and subgroup composition,
#'   for observed data, forecasts or both.
#'
#' @section Caution:
#' mortensor is experimental software written largely with AI coding
#' assistants ("vibe-coded"). It has unit tests that check its methods
#' against published results and reference implementations, but it has not
#' been peer reviewed and its interfaces still change. Use it
#' at your own risk and check results independently before using them in
#' publications, reports or decisions.
#'
#' @keywords internal
#' @import data.table
#' @importFrom stats median quantile sd var setNames
#' @importFrom utils head tail
"_PACKAGE"
## data.table non-standard evaluation symbols
utils::globalVariables(c(
  "e1", "e2", "n_h", "h_max", "cells", "skipped", ".", ".N", ".SD", ".I", ".GRP", ":=", "..cols",
  "age_lo", "age_hi", "period_lo", "period_hi", "value", "measure",
  "stratum", "age", "period", "draw", "rate", "deaths", "exposure",
  "provenance", "grid_age", "grid_period", "horizon", "origin", "model",
  "observed", "q", "metric", "target", "estimate", "lower", "upper",
  "nalo", "nahi", "nplo", "nphi", "area", "pw", "full", "value.d", "value.e",
  "value.next", "k", "ratio", "label", "sex", "d", "e", "g_lo", "g_hi", "N",
  "a1", "a2", "p1", "p2", "a", "p", "i.value", "i.provenance", "lr", "scale", "component", "draw",
  "crps", "wis", "ae", "cov50", "cov80", "cov95", "pit", "at", "status", "rep", "parts",
  ".data", "q0.5", "lo", "hi", "med", "x", "x_hi", "level", "coverage", "v", "region", "quintile", "m", "severity",
  "lev", "w", "year", "dep", "population", "i.population", "cause", "dic", "index", "factor", "kind", "score",
  "group", "mode", "time", "cell", "n_members", "symbol", "raw", "age_group", "cause_label", "idx",
  "group_", "new_", "n_", "first_", "members", "n", "s", ".value", "nxt", "u", "g", "n_lev", "sum_all", "var",
  "age_i", "effect", "q0.05", "q0.95", "residual", "p_positive", "ages", "panel", "age_f", "panel2", "col", "row",
  "band", "ref", "noise", "density",
  "comparison", "xvar", "fvar", "max_residual", "q0.025", "q0.975", "forecast", "fitted", "prep", "V1",
  "unit", "loss", "selected", "se", "method", "col", "setting", "score", "part", "from", "to"
))
