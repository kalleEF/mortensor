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
#' * **Models**: [rwd()], [lee_carter()], [tucker3()], [arima_series()],
#'   [ets_series()], [apc_bayes()] (package bamp) and [mt_ensemble()], fitted
#'   with [mt_fit()] and forecast with [mt_forecast()] as simulated paths
#'   ([mt_draws]).
#' * **Evaluation**: [mt_backtest()] and [mt_score()].
#' * **Derived quantities and plots**: [life_table()], [life_expectancy()],
#'   [mt_deaths()], [mt_plot_forecast()] and the other `mt_plot_*()`
#'   functions.
#'
#' @keywords internal
#' @import data.table
#' @importFrom stats median quantile sd var setNames
#' @importFrom utils head tail
"_PACKAGE"
## data.table non-standard evaluation symbols
utils::globalVariables(c(
  ".", ".N", ".SD", ".I", ".GRP", ":=", "..cols",
  "age_lo", "age_hi", "period_lo", "period_hi", "value", "measure",
  "stratum", "age", "period", "draw", "rate", "deaths", "exposure",
  "provenance", "grid_age", "grid_period", "horizon", "origin", "model",
  "observed", "q", "metric", "target", "estimate", "lower", "upper",
  "nalo", "nahi", "nplo", "nphi", "area", "pw", "full", "value.d", "value.e",
  "value.next", "k", "ratio", "label", "sex", "d", "e", "g_lo", "g_hi", "N",
  "a1", "a2", "p1", "p2", "a", "p", "i.value", "i.provenance", "lr", "scale",
  "crps", "wis", "ae", "cov50", "cov80", "cov95", "pit", "at", "status", "rep",
  ".data", "q0.5", "lo", "hi", "med", "x", "x_hi", "level", "coverage", "v"
))
