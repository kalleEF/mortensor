#' WLS reconciliation of independent disease forecasts
#'
#' Adjusts independently-forecast disease components so they sum
#' to the all-cause forecast.  Adjustment is proportional to each
#' component's forecast variance (wide PIs absorb more gap).
#'
#' @param mx_all_fc  3D array \[n_a x h x n_g\] of all-cause point forecast.
#' @param mx_indep   Named list of 3D arrays: per-disease point forecasts.
#' @param bs_indep   Named list of bootstrap results (each with
#'   `mx_lower`, `mx_upper`, `mx_median`).
#' @param disease_order Character vector of disease names.
#' @param recon_cap  Soft cap: max |delta_k| as fraction of point forecast.
#' @return Named list with `mx_recon`, `bs_recon`, `delta`, `sigma_sq`.
#' @export
wls_reconcile <- function(mx_all_fc, mx_indep, bs_indep,
                          disease_order, recon_cap = 0.5) {
  n_a <- dim(mx_all_fc)[1]; h <- dim(mx_all_fc)[2]; n_g <- dim(mx_all_fc)[3]

  # Estimate per-cell variance from PI width
  sigma_sq <- list()
  for (nm in disease_order) {
    pi_width <- bs_indep[[nm]]$mx_upper - bs_indep[[nm]]$mx_lower
    s_k <- pmax(pi_width / (2 * 1.96), 1e-10)
    sigma_sq[[nm]] <- s_k^2
  }

  # Gap = all-cause - sum(independent)
  mx_total_indep <- array(0, dim = c(n_a, h, n_g))
  for (nm in disease_order)
    mx_total_indep <- mx_total_indep + mx_indep[[nm]]
  gap <- mx_all_fc - mx_total_indep

  # Sum of variances
  var_total <- array(0, dim = c(n_a, h, n_g))
  for (nm in disease_order)
    var_total <- var_total + sigma_sq[[nm]]
  var_total <- pmax(var_total, 1e-20)

  # WLS adjustment: delta_k = gap * sigma_k^2 / sum(sigma_j^2)
  delta <- list()
  mx_recon <- list()
  bs_recon <- list()

  for (nm in disease_order) {
    d_k <- gap * sigma_sq[[nm]] / var_total

    # Soft cap on adjustment
    cap_abs <- recon_cap * pmax(mx_indep[[nm]], 1e-12)
    d_k <- pmin(pmax(d_k, -cap_abs), cap_abs)

    delta[[nm]] <- d_k
    mx_recon[[nm]] <- pmax(mx_indep[[nm]] + d_k, 1e-12)
    bs_recon[[nm]] <- list(
      mx_lower  = pmax(bs_indep[[nm]]$mx_lower  + d_k, 1e-12),
      mx_upper  = pmax(bs_indep[[nm]]$mx_upper  + d_k, 1e-12),
      mx_median = pmax(bs_indep[[nm]]$mx_median + d_k, 1e-12))
  }

  # Verify coherence; close residual gap if needed
  mx_recon_check <- array(0, dim = c(n_a, h, n_g))
  for (nm in disease_order)
    mx_recon_check <- mx_recon_check + mx_recon[[nm]]
  max_gap <- max(abs(mx_recon_check - mx_all_fc) / pmax(mx_all_fc, 1e-12))

  if (max_gap > 1e-4 && "residual" %in% disease_order) {
    final_gap <- mx_all_fc - mx_recon_check
    mx_recon[["residual"]] <- pmax(mx_recon[["residual"]] + final_gap, 1e-12)
    bs_recon[["residual"]]$mx_lower  <-
      pmax(bs_recon[["residual"]]$mx_lower  + final_gap, 1e-12)
    bs_recon[["residual"]]$mx_upper  <-
      pmax(bs_recon[["residual"]]$mx_upper  + final_gap, 1e-12)
    bs_recon[["residual"]]$mx_median <-
      pmax(bs_recon[["residual"]]$mx_median + final_gap, 1e-12)
  }

  list(mx_recon = mx_recon, bs_recon = bs_recon,
       delta = delta, sigma_sq = sigma_sq)
}
