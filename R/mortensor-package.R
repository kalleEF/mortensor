#' mortensor: flexible mortality forecasting across data granularities
#'
#' The package is organised in layers:
#'
#' * **Data**: [mt_data()] stores deaths, exposures and populations observed
#'   at any age and period resolution ([age_groups()], [single_ages()],
#'   [single_years()]); [as_mt_data()] converts HMD-based objects from
#'   \pkg{StMoMo} and \pkg{demography}; [mt_validate()] checks plausibility.
#' * **Harmonisation**: [mt_regroup()] aggregates exactly, [mt_ungroup()]
#'   estimates single-year data with the penalized composite link model.
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
  "a1", "a2", "p1", "p2", "a", "p", "i.value", "i.provenance"
))
