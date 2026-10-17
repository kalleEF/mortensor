#' Different models for different groups of strata
#'
#' Fits one model specification to each group of strata, defined by the
#' levels of a stratum variable, and combines their forecasts into one set
#' of draws. A typical use is a trusted all-cause forecast from one model
#' and cause-specific forecasts from a model that pools sparse causes,
#' reconciled so that the causes add up to all causes:
#'
#' `reconciled(stratified("cause", all = lee_carter(), .default = li_lee(common_by = c("sex", "cause"))), h,
#' method = "proportional", fixed = all_cause_labels)`.
#'
#' Joint models (for example [li_lee()]) are joint within their group only;
#' the draws of different groups are independent (inside [reconciled()] the
#' empirical reordering then restores the co-movement of the training years).
#' Each group is fitted and forecast as if it were the whole data, so
#' backtests, [anchored()], [recalibrated()] and [reconciled()] work as for
#' any model.
#'
#' @param var Stratum variable that defines the groups.
#' @param ... Named model specifications: each name is a level of `var`
#'   whose strata the model is fitted to.
#' @param .default Model for the strata of all other levels (`NULL`: every
#'   level must be named).
#' @return A model specification.
#' @examples
#' x <- mt_example(deaths = "single", years = 1990:2019)
#' spec <- stratified("sex", female = lee_carter(), male = rwd())
#' fc <- mt_forecast(mt_fit(spec, x), h = 5, n_draws = 50, seed = 1)
#' fc$model
#' @export
stratified <- function(var, ..., .default = NULL) {
  specs <- list(...)
  if (!is.character(var) || length(var) != 1) mt_abort("`var` must be one stratum variable.", class = "argument")
  if (!length(specs) || is.null(names(specs)) || any(!nzchar(names(specs))) || anyDuplicated(names(specs))) {
    mt_abort("Give the models as named arguments, one per level of `var`.", class = "argument")
  }
  all_specs <- c(specs, if (!is.null(.default)) list(.default))
  if (!all(vapply(all_specs, inherits, TRUE, "mt_spec"))) {
    mt_abort("Every model must be a model specification.", class = "argument")
  }
  parts <- c(sprintf("%s: %s", names(specs), vapply(specs, `[[`, "", "label")),
             if (!is.null(.default)) sprintf("other: %s", .default$label))
  new_mt_spec("stratified", list(var = var, specs = specs, default = .default),
              min_periods = max(vapply(all_specs, `[[`, 0, "min_periods")),
              smooth = all(vapply(all_specs, function(s) isTRUE(s$smooth), TRUE)),
              label = sprintf("By %s (%s)", var, paste(parts, collapse = "; ")))
}

#' @export
mt_fit.mt_spec_stratified <- function(spec, data, prep = NULL, cores = 1, ...) {
  check_class(data, "mt_data", "data")
  p <- spec$params
  st <- data$strata
  if (!p$var %in% setdiff(names(st), c("stratum", "label"))) {
    mt_abort(sprintf("\"%s\" is not a stratum variable.", p$var), class = "argument")
  }
  lev <- as.character(st[[p$var]])
  unknown <- setdiff(names(p$specs), lev)
  if (length(unknown)) {
    mt_abort(sprintf("%s: no strata with %s = %s.", spec$label, p$var, paste(unknown, collapse = ", ")),
             class = "argument")
  }
  groups <- lapply(names(p$specs), function(n) which(lev == n))
  names(groups) <- names(p$specs)
  specs <- p$specs
  rest <- which(!lev %in% names(p$specs))
  if (length(rest)) {
    if (is.null(p$default)) {
      mt_abort(sprintf("No model for %s = %s; name it or give `.default`.", p$var,
                       paste(unique(lev[rest]), collapse = ", ")), class = "argument")
    }
    groups$other <- rest
    specs$other <- p$default
  }
  fits <- lapply(names(groups), function(g) {
    mt_fit(specs[[g]], strata_subset(data, groups[[g]]), prep = prep, cores = cores)
  })
  combine_group_fits(fits, groups, data, spec)
}
