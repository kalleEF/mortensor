# Internal package environment for runtime state
.pkg_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  .pkg_env$has_stmomo <- requireNamespace("StMoMo", quietly = TRUE)
  .pkg_env$has_demography <- requireNamespace("demography", quietly = TRUE)
  .pkg_env$has_readxl <- requireNamespace("readxl", quietly = TRUE)
  invisible()
}

# Suppress R CMD check notes for data.table NSE bindings
utils::globalVariables(c(
  ".", ".N", ".SD", ":=",
  "age", "age_group", "age_start", "all_causes", "cause", "complete",
  "deaths_all", "deaths_cause", "deaths_disease_sum",
  "fitted", "hi95", "imd_decile", "lo95", "log_mx",
  "mase", "mase_combined", "mase_rank", "method", "model",
  "mx", "mx_1yr", "mx_1yr_hi95", "mx_1yr_lo95",
  "mx_adjusted", "mx_all", "mx_all_fc", "mx_fitted",
  "mx_hi", "mx_hi95", "mx_hi95_adj",
  "mx_lo", "mx_lo95", "mx_lo95_adj", "mx_lo_s",
  "mx_hi_s",
  "mx_median", "mx_median_adj", "mx_obs", "mx_point", "mx_smooth",
  "mx_sum", "pop_fill", "population", "quintile", "ratio", ".group",
  "reorder", "sd", "sex", "share", "simulate", "smooth",
  "stage", "type", "value", "w", "xmax", "xmin", "y", "year",
  "grp",
  "log_lo", "log_hi", "residual"
))
