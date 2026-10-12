#' Smooth death rates on the grid
#'
#' Poisson P-spline smoothing of death rates with the exposures as offset
#' (Currie, Durban and Eilers 2004 for the 2D case): the same engine as
#' [mt_ungroup()] with every grid cell as its own observation. Deaths are
#' replaced by exposure times the smoothed rate; exposures are unchanged.
#' Deaths and exposures must be available on the grid (use [mt_ungroup()]
#' first for grouped data).
#'
#' @param data An [mt_data] object.
#' @param dims 1 (over age, separately per period) or 2 (age-period
#'   surface; P-splines only).
#' @param method `"pspline"` (Poisson P-splines, default) or `"whittaker"`
#'   (Whittaker-Henderson smoothing of log rates over age with weights equal
#'   to the deaths, the approximate inverse variance of a log rate;
#'   smoothing parameter by generalised cross-validation unless `lambda` is
#'   given; an option ported from mortensor 0.1).
#' @inheritParams mt_ungroup
#' @return An `mt_data` object with smoothed deaths; diagnostics as attribute
#'   `"smooth"`.
#' @references
#' Currie ID, Durban M, Eilers PHC (2004). Smoothing and forecasting
#' mortality rates. *Statistical Modelling* 4(4): 279-298.
#' @examples
#' x <- mt_example(years = 2001:2003, deaths = "single")
#' s <- mt_smooth(x)
#' head(attr(s, "smooth"))
#' @export
mt_smooth <- function(data, dims = 1, lambda = NULL, criterion = c("bic", "aic"), kr = 2, deg = 3,
                      pord = 2, lambda_range = c(1e-2, 1e8), method = c("pspline", "whittaker")) {
  check_class(data, "mt_data", "data")
  criterion <- match.arg(criterion)
  method <- match.arg(method)
  check_number(dims, "dims", lower = 1, upper = 2, integerish = TRUE)
  g <- mt_grid(data)
  for (m in c("deaths", "exposure")) {
    if (is.null(data$obs[[m]])) mt_abort(sprintf("Smoothing needs a %s table.", m), class = "missing_measure")
    if (needs_ungroup(data$obs[[m]], g)) {
      mt_abort(sprintf("The %s table is coarser than the grid; use mt_ungroup() first.", m), class = "resolution")
    }
  }
  if (method == "whittaker") {
    res <- whittaker_smooth_obs(data, g, lambda, pord)
    data$obs$deaths <- res$obs
    data$grid <- g
    attr(data, "smooth") <- res$diag
    return(data)
  }
  ctl <- list(lambda = lambda, criterion = criterion, kr = kr, deg = deg, pord = pord,
              lambda_range = lambda_range, preserve_totals = FALSE, keep_observed = FALSE,
              provenance = "smoothed", method_label = "pspline")
  res <- ungroup_rates(data$obs$deaths, data$obs$exposure, g, data$omega, ctl, dims)
  data$obs$deaths <- res$obs
  data$grid <- g
  attr(data, "smooth") <- res$diag
  data
}

whittaker_smooth_obs <- function(data, g, lambda, pord) {
  D <- mt_array(data, "deaths")
  E <- mt_array(data, "exposure")
  na <- dim(D)[1]
  Dm <- diff(diag(na), differences = pord)
  P <- crossprod(Dm)
  out <- D
  diag_rows <- list()
  for (s in seq_len(dim(D)[3])) {
    for (j in seq_len(dim(D)[2])) {
      y <- log(pmax(D[, j, s], 0.5) / E[, j, s])
      w <- pmax(D[, j, s], 0.5)
      fit_l <- function(l) {
        A <- diag(w) + l * P
        z <- solve(A, w * y)
        tr <- sum(diag(solve(A, diag(w))))
        list(z = z, gcv = sum(w * (y - z)^2) / (na - tr)^2 * na)
      }
      l <- lambda %||% 10^stats::optimize(function(ll) fit_l(10^ll)$gcv, c(-2, 8))$minimum
      z <- fit_l(l)$z
      out[, j, s] <- E[, j, s] * exp(z)
      diag_rows[[length(diag_rows) + 1L]] <- data.table(measure = "deaths", stratum = s, period = g$period$labels[j],
                                                        method = "whittaker", lambda_age = l)
    }
  }
  obs <- data.table(stratum = rep(seq_len(dim(D)[3]), each = na * dim(D)[2]),
                    a = rep(seq_len(na), dim(D)[2] * dim(D)[3]),
                    p = rep(rep(seq_len(dim(D)[2]), each = na), dim(D)[3]),
                    value = as.vector(out), provenance = "smoothed")
  list(obs = grid_table(obs, g), diag = rbindlist(diag_rows))
}
