#' Preprocessing pipelines
#'
#' `mt_prep()` bundles preprocessing steps that are applied to the data
#' before a model is fitted. Because the steps are part of the fit, a
#' backtest ([mt_backtest()]) re-runs them on every training window, so that
#' ungrouping or smoothing never sees data after the forecast origin.
#'
#' Steps:
#'
#' * `step_ungroup(...)`: [mt_ungroup()];
#' * `step_smooth(...)`: [mt_smooth()]. The model is fitted to the
#'   smoothed rates; residuals and observation noise refer to the raw
#'   deaths. Smoothing over time (`along = "period"`, or `dims = 2`)
#'   narrows the forecast intervals of models with a period index, so check
#'   coverage in a backtest and consider [recalibrated()];
#' * `step_regroup(age, period)`: [mt_regroup()];
#' * `step_grid(age, period)`: [mt_set_grid()];
#' * `step_window(years)`: keeps the last `years` years of the training
#'   data, a fitting period chosen relative to the forecast origin (for
#'   example after a change of trend; compare windows with [tuned()] or
#'   [mt_backtest()]);
#' * `step_shocks(periods, action)`: periods with exceptional mortality (for
#'   example 2020-2022) are either replaced by log-linear interpolation of
#'   the rates of the neighbouring periods (`"interpolate"`, deaths become
#'   exposure times the interpolated rate) or removed (`"exclude"`, only at
#'   the start or end of the data).
#'
#' @param ... Steps, applied in order.
#' @return An object of class `mt_prep`; `mt_prep_apply()` returns the
#'   processed [mt_data] object.
#' @examples
#' prep <- mt_prep(step_ungroup(), step_shocks(2005, action = "interpolate"))
#' prep
#' x <- mt_prep_apply(prep, mt_example(years = 2001:2008))
#' @export
mt_prep <- function(...) {
  steps <- list(...)
  ok <- vapply(steps, inherits, logical(1), what = "mt_step")
  if (!all(ok)) mt_abort("All arguments of mt_prep() must be steps (step_*()).", class = "prep")
  structure(list(steps = steps), class = "mt_prep")
}

new_step <- function(type, args) {
  structure(list(type = type, args = args), class = c(paste0("mt_step_", type), "mt_step"))
}

#' @rdname mt_prep
#' @export
step_ungroup <- function(...) new_step("ungroup", list(...))

#' @rdname mt_prep
#' @export
step_smooth <- function(...) new_step("smooth", list(...))

#' @rdname mt_prep
#' @param age,period Axes for `step_regroup()` and `step_grid()`.
#' @export
step_regroup <- function(age = NULL, period = NULL) new_step("regroup", list(age = age, period = period))

#' @rdname mt_prep
#' @export
step_grid <- function(age = NULL, period = NULL) new_step("grid", list(age = age, period = period))

#' @rdname mt_prep
#' @param years Number of most recent training years kept by
#'   `step_window()`.
#' @export
step_window <- function(years) {
  check_number(years, "years", lower = 2, integerish = TRUE)
  new_step("window", list(years = years))
}

#' @rdname mt_prep
#' @param periods Calendar years affected by a shock.
#' @param action `"interpolate"` or `"exclude"`.
#' @export
step_shocks <- function(periods, action = c("interpolate", "exclude")) {
  action <- match.arg(action)
  new_step("shocks", list(periods = periods, action = action))
}

#' @rdname mt_prep
#' @param prep An `mt_prep` object or `NULL` (no preprocessing).
#' @param data An [mt_data] object.
#' @export
mt_prep_apply <- function(prep, data) {
  if (is.null(prep)) return(data)
  check_class(prep, "mt_prep", "prep")
  for (st in prep$steps) data <- apply_step(st, data)
  data
}

# The preprocessed data and, when the pipeline smooths, the same data
# without the smoothing steps (raw rates on the modelling grid): models are
# fitted to the smoothed rates, observation noise and residuals come from
# the raw rates
mt_prep_apply_raw <- function(prep, data) {
  if (is.null(prep)) return(list(data = data, raw = NULL))
  check_class(prep, "mt_prep", "prep")
  raw <- NULL
  for (st in prep$steps) {
    if (st$type == "smooth") {
      if (is.null(raw)) raw <- data
      data <- apply_step(st, data)
    } else {
      data <- apply_step(st, data)
      if (!is.null(raw)) raw <- apply_step(st, raw)
    }
  }
  list(data = data, raw = raw)
}

apply_step <- function(step, data) {
  a <- step$args
  switch(step$type,
         ungroup = do.call(mt_ungroup, c(list(data), a)),
         smooth = do.call(mt_smooth, c(list(data), a)),
         regroup = mt_regroup(data, a$age, a$period),
         grid = mt_set_grid(data, a$age, a$period),
         shocks = apply_shocks(data, a$periods, a$action),
         window = {
           last <- max(mt_grid(data)$period$lower)
           mt_window(data, from = last - a$years + 1)
         },
         mt_abort(sprintf("Unknown step '%s'.", step$type), class = "prep"))
}

apply_shocks <- function(data, periods, action) {
  g <- mt_grid(data)
  yrs <- g$period$lower
  hit <- yrs %in% periods
  if (!any(hit)) return(data)
  if (action == "exclude") {
    interior <- hit & cumsum(!hit) > 0 & rev(cumsum(rev(!hit))) > 0
    if (any(interior)) {
      mt_abort("Shock periods inside the data cannot be excluded; use action = \"interpolate\".", class = "prep")
    }
    keep <- yrs[!hit]
    out <- mt_window(data, from = min(keep), to = max(keep))
    attr(out, "excluded_periods") <- yrs[hit]
    return(out)
  }
  D <- mt_array(data, "deaths")
  E <- mt_array(data, "exposure")
  base <- which(!hit)
  if (max(which(hit)) > max(base) || min(which(hit)) < min(base)) {
    mt_abort("Shock periods at the start or end cannot be interpolated; use action = \"exclude\".", class = "prep")
  }
  lr <- log(pmax(D, 0.5) / E)
  for (j in which(hit)) {
    lo <- max(base[base < j])
    hi <- min(base[base > j])
    w <- (j - lo) / (hi - lo)
    D[, j, ] <- E[, j, ] * exp((1 - w) * lr[, lo, ] + w * lr[, hi, ])
  }
  obs <- copy(data$obs$deaths)
  idx <- obs_grid_index(obs, g)
  sel <- !is.na(idx$p1) & hit[pmax(idx$p1, 1L)]
  obs[sel, value := D[cbind(idx$a1[sel], idx$p1[sel], obs$stratum[sel])]]
  obs[sel, provenance := "interpolated"]
  data$obs$deaths <- obs
  data
}

#' @export
print.mt_prep <- function(x, ...) {
  cat("<mt_prep>", length(x$steps), "step(s)\n")
  for (st in x$steps) {
    a <- st$args
    vals <- vapply(a, function(v) paste(format(v), collapse = ","), "")
    desc <- if (length(a)) paste(names(a), vals, sep = "=", collapse = ", ") else ""
    cat(sprintf("  step_%s(%s)\n", st$type, desc))
  }
  invisible(x)
}
