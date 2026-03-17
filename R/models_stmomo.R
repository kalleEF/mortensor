#' Stochastic Mortality Model (StMoMo) forecast
#'
#' Fits an APC, CBD, RH, or M6 stochastic mortality model per group using
#' the StMoMo package and produces forecasts.
#'
#' @param mx_array 3D array `[n_a x n_t x n_g]` of mortality rates.
#' @param pop_array 3D array of population exposure.
#' @param years Integer vector of calendar years.
#' @param ages Integer vector of age-group starts.
#' @param groups_dt data.frame with group definitions.
#' @param h Integer forecast horizon.
#' @param model_type One of `"APC"`, `"CBD"`, `"RH"`, `"M6"`.
#' @return A list with `mx_forecast` `[n_a x h x n_g]` (NULL if StMoMo
#'   not available).
#' @export
stmomo_forecast <- function(mx_array, pop_array, years, ages,
                            groups_dt, h, model_type = "APC") {
  if (!requireNamespace("StMoMo", quietly = TRUE))
    stop("Install the 'StMoMo' package to use stmomo_forecast().")

  n_a <- dim(mx_array)[1]; n_g <- dim(mx_array)[3]
  age_idx <- seq_len(n_a)
  mx_fc <- array(NA_real_, dim = c(n_a, h, n_g))
  n_fail <- 0

  for (k in seq_len(n_g)) {
    tryCatch({
      rate_mat <- pmax(mx_array[, , k], 1e-12)
      pop_mat  <- pmax(pop_array[, , k], 1)
      dxt <- round(rate_mat * pop_mat)
      model <- switch(toupper(model_type),
        CBD  = StMoMo::cbd(link = "logit"),
        RH   = StMoMo::rh(link = "log", cohortAgeFun = "1"),
        M6   = StMoMo::m6(link = "log"),
        StMoMo::apc(link = "log")
      )
      sfit <- StMoMo::fit(model, Dxt = dxt, Ext = pop_mat,
                          ages = age_idx, years = years,
                          wxt = 1 * (dxt > 0))
      fc <- forecast::forecast(sfit, h = h)
      mx_fc[, , k] <- fc$rates
    }, error = function(e) { n_fail <<- n_fail + 1; NULL })
  }

  if (all(is.na(mx_fc))) return(list(mx_forecast = NULL))
  list(mx_forecast = mx_fc)
}
