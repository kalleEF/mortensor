#' Lee-Carter forecast via demography package
#'
#' Fits Lee-Carter models per group using [demography::lca()] and forecasts.
#'
#' @param mx_array 3D array `[n_a x n_t x n_g]` of mortality rates.
#' @param pop_array 3D array of population exposure (same dims).
#' @param years Integer vector of calendar years.
#' @param ages Integer vector of age-group starts.
#' @param groups_dt A data.frame with group definitions (must have `sex`,
#'   `quintile` columns, or similar identifiers).
#' @param h Integer forecast horizon.
#' @param return_fitted Logical; return fitted values?
#' @param adjust Adjustment method: `"dt"`, `"dxt"`, `"e0"`, or `"none"`.
#' @param chooseperiod Logical; let lca pick the best fitting period?
#' @return A list with `mx_forecast` `[n_a x h x n_g]` and optionally
#'   `mx_fitted`.
#' @export
lca_forecast <- function(mx_array, pop_array, years, ages, groups_dt, h,
                         return_fitted = FALSE,
                         adjust = "dt", chooseperiod = FALSE) {
  if (!requireNamespace("demography", quietly = TRUE))
    stop("Install the 'demography' package to use lca_forecast().")

  n_a <- dim(mx_array)[1]; n_t <- dim(mx_array)[2]; n_g <- dim(mx_array)[3]
  mx_fc <- array(NA_real_, dim = c(n_a, h, n_g))
  mx_fit <- if (return_fitted) array(NA_real_, dim = c(n_a, n_t, n_g)) else NULL

  for (k in seq_len(n_g)) {
    rate_mat <- pmax(mx_array[, , k], 1e-12)
    pop_mat  <- pmax(pop_array[, , k], 1)
    lbl <- make_group_label(groups_dt, k)
    ddata <- tryCatch(
      demography::demogdata(data = rate_mat, pop = pop_mat,
                            ages = ages, years = years,
                            type = "mortality", label = lbl, name = "total"),
      error = function(e) NULL)
    if (is.null(ddata)) next
    lc_fit <- tryCatch(
      suppressWarnings(demography::lca(ddata, series = "total",
                                       max.age = max(ages),
                                       adjust = adjust,
                                       chooseperiod = chooseperiod)),
      error = function(e) NULL)
    if (is.null(lc_fit)) next
    lc_fc <- tryCatch(forecast::forecast(lc_fit, h = h),
                      error = function(e) NULL)
    if (!is.null(lc_fc) && !is.null(lc_fc$rate[[1]]))
      mx_fc[, , k] <- as.matrix(lc_fc$rate[[1]])
    if (return_fitted && !is.null(lc_fit)) {
      fv <- lc_fit$fitted$y
      if (is.list(fv)) fv <- as.matrix(fv[[1]])
      mx_fit[, , k] <- exp(fv)
    }
  }
  list(mx_forecast = mx_fc, mx_fitted = mx_fit)
}

#' Hyndman-Ullah Functional Demographic Model (FDM) forecast
#'
#' @inheritParams lca_forecast
#' @param order Number of principal components (default 6).
#' @param method Time-series extrapolation: `"arima"` or `"ets"`.
#' @return A list with `mx_forecast` and optionally `mx_fitted`.
#' @export
fdm_forecast <- function(mx_array, pop_array, years, ages, groups_dt, h,
                         return_fitted = FALSE,
                         order = 6, method = "arima") {
  if (!requireNamespace("demography", quietly = TRUE))
    stop("Install the 'demography' package to use fdm_forecast().")

  n_a <- dim(mx_array)[1]; n_t <- dim(mx_array)[2]; n_g <- dim(mx_array)[3]
  mx_fc <- array(NA_real_, dim = c(n_a, h, n_g))
  mx_fit <- if (return_fitted) array(NA_real_, dim = c(n_a, n_t, n_g)) else NULL
  order_fdm <- min(order, n_a - 1)

  for (k in seq_len(n_g)) {
    rate_mat <- pmax(mx_array[, , k], 1e-12)
    pop_mat  <- pmax(pop_array[, , k], 1)
    lbl <- make_group_label(groups_dt, k)
    ddata <- tryCatch(
      demography::demogdata(data = rate_mat, pop = pop_mat,
                            ages = ages, years = years,
                            type = "mortality", label = lbl, name = "total"),
      error = function(e) NULL)
    if (is.null(ddata)) next
    fdm_fit <- tryCatch(
      suppressWarnings(demography::fdm(ddata, series = "total",
                                       order = order_fdm,
                                       max.age = max(ages))),
      error = function(e) NULL)
    if (is.null(fdm_fit)) next
    fdm_fc <- tryCatch(forecast::forecast(fdm_fit, h = h, method = method),
                       error = function(e) NULL)
    if (!is.null(fdm_fc) && !is.null(fdm_fc$rate[[1]]))
      mx_fc[, , k] <- as.matrix(fdm_fc$rate[[1]])
    if (return_fitted && !is.null(fdm_fit)) {
      fv <- fdm_fit$fitted$y
      if (is.list(fv)) fv <- as.matrix(fv[[1]])
      mx_fit[, , k] <- exp(fv)
    }
  }
  list(mx_forecast = mx_fc, mx_fitted = mx_fit)
}


# ---- Internal helper ----

#' @keywords internal
make_group_label <- function(groups_dt, k) {
  parts <- c()
  if ("sex" %in% names(groups_dt))
    parts <- c(parts, tolower(groups_dt$sex[k]))
  if ("quintile" %in% names(groups_dt))
    parts <- c(parts, paste0("q", groups_dt$quintile[k]))
  if ("deprivation" %in% names(groups_dt))
    parts <- c(parts, paste0("d", groups_dt$deprivation[k]))
  if (length(parts) == 0) parts <- paste0("g", k)
  paste(parts, collapse = "_")
}
