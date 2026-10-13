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
#'   given) or `"gam"` (Gaussian
#'   additive model of log rates on age per period with deaths as weights,
#'   `mgcv::gam()` with REML).
#' @param ages Range of ages that is smoothed, `c(lower, upper)`; other
#'   ages keep their observed deaths and do not enter the fit. By default
#'   all ages except the first year of life when the grid starts with age 0
#'   as its own interval (infant mortality does not follow the age pattern
#'   of the following ages).
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
                      pord = 2, lambda_range = c(1e-2, 1e8), method = c("pspline", "whittaker", "gam"),
                      ages = NULL) {
  check_class(data, "mt_data", "data")
  criterion <- match.arg(criterion)
  method <- match.arg(method)
  check_number(dims, "dims", lower = 1, upper = 2, integerish = TRUE)
  g <- mt_grid(data)
  if (is.null(ages)) ages <- if (g$age$lower[1] == 0 && g$age$upper[1] == 1 && length(g$age) > 3) c(1, Inf) else
    c(-Inf, Inf)
  if (length(ages) != 2) mt_abort("`ages` must be c(lower, upper).", class = "argument")
  inside <- g$age$lower >= ages[1] - 1e-9 & g$age$upper <= ages[2] + 1e-9
  if (!any(inside)) mt_abort("No grid age is inside `ages`.", class = "argument")
  if (!all(inside)) {
    # ages outside the range (by default the first year of life, whose rate
    # does not follow the age pattern of the other ages) are neither smoothed
    # nor used in the fit
    data$grid <- g
    sm <- mt_smooth(mt_subset(data, ages = ages), dims = dims, lambda = lambda, criterion = criterion, kr = kr,
                    deg = deg, pord = pord, lambda_range = lambda_range, method = method, ages = c(-Inf, Inf))
    d0 <- data$obs$deaths
    keep <- !(d0$age_lo >= ages[1] - 1e-9 & d0$age_hi <= ages[2] + 1e-9)
    data$obs$deaths <- rbind(sm$obs$deaths, d0[keep])
    setorderv(data$obs$deaths, c("stratum", "period_lo", "age_lo"))
    attr(data, "smooth") <- attr(sm, "smooth")
    return(data)
  }
  for (m in c("deaths", "exposure")) {
    if (is.null(data$obs[[m]])) mt_abort(sprintf("Smoothing needs a %s table.", m), class = "missing_measure")
    if (needs_ungroup(data$obs[[m]], g)) {
      mt_abort(sprintf("The %s table is coarser than the grid; use mt_ungroup() first.", m), class = "resolution")
    }
  }
  if (method %in% c("whittaker", "gam")) {
    res <- if (method == "whittaker") whittaker_smooth_obs(data, g, lambda, pord) else gam_smooth_obs(data, g)
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

# Matrix of divided differences of order k at the points x, scaled by the
# mean spacing so that equal spacing gives the ordinary differences
divided_differences <- function(x, k) {
  n <- length(x)
  M <- diag(n)
  h <- mean(diff(x))
  for (j in seq_len(k)) {
    m <- nrow(M)
    span <- x[(j + 1):n] - x[1:(n - j)]
    M <- (M[-1, , drop = FALSE] - M[-m, , drop = FALSE]) / (span / (j * h))
  }
  M
}

whittaker_smooth_obs <- function(data, g, lambda, pord) {
  D <- mt_array(data, "deaths")
  E <- mt_array(data, "exposure")
  na <- dim(D)[1]
  # divided differences at the interval midpoints, so that age groups of
  # different widths are smoothed on the age scale
  mid <- axis_midpoints(g$age)
  n <- length(mid)
  if (is_open(g$age) && n > 1) mid[n] <- g$age$lower[n] + (g$age$upper[n - 1] - g$age$lower[n - 1]) / 2
  Dm <- divided_differences(mid, pord)
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

gam_smooth_obs <- function(data, g) {
  need_pkg("mgcv", "for method = \"gam\"")
  D <- mt_array(data, "deaths")
  E <- mt_array(data, "exposure")
  na <- dim(D)[1]
  age <- axis_midpoints(g$age, data$omega)
  out <- D
  diag_rows <- list()
  for (s in seq_len(dim(D)[3])) {
    for (j in seq_len(dim(D)[2])) {
      df <- data.frame(y = log(pmax(D[, j, s], 0.5) / E[, j, s]), w = pmax(D[, j, s], 0.5), age = age)
      wts <- df$w / mean(df$w)
      fit <- mgcv::gam(y ~ s(age, k = min(20, na - 1)), weights = wts, data = df, method = "REML")
      out[, j, s] <- E[, j, s] * exp(stats::fitted(fit))
      diag_rows[[length(diag_rows) + 1L]] <- data.table(measure = "deaths", stratum = s, period = g$period$labels[j],
                                                        method = "gam", edf = sum(fit$edf))
    }
  }
  obs <- data.table(stratum = rep(seq_len(dim(D)[3]), each = na * dim(D)[2]),
                    a = rep(seq_len(na), dim(D)[2] * dim(D)[3]),
                    p = rep(rep(seq_len(dim(D)[2]), each = na), dim(D)[3]),
                    value = as.vector(out), provenance = "smoothed")
  list(obs = grid_table(obs, g), diag = rbindlist(diag_rows))
}
