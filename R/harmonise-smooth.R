#' Smooth death rates on the grid
#'
#' Smooths death rates over age, over time or over both, with the exposures
#' as weights or offset. Deaths are replaced by exposure times the smoothed
#' rate; exposures are unchanged. Cells with missing deaths do not enter the
#' fit and stay missing (only the 2D GAM predicts them). Deaths and
#' exposures must be available on the grid (use [mt_ungroup()] first for
#' grouped data). In a forecasting
#' pipeline, use it as a preprocessing step (`mt_prep(step_smooth(...))`):
#' models are then fitted to the smoothed rates, while their residuals,
#' observed rates, backtest scales and jump-off offsets come from the deaths
#' before smoothing ("smooth input, raw residuals"). Predictive draws add
#' Poisson deaths given the model's rates, as without smoothing.
#'
#' Methods and directions:
#'
#' * `dims = 1, along = "age"` (default): each year separately over age,
#'   with Poisson P-splines and the exposures as offset (`"pspline"`; the
#'   one-dimensional case of Currie, Durban and Eilers 2004), the weighted
#'   Whittaker smoother of log rates with the deaths as weights
#'   (`"whittaker"`; Whittaker 1923; Eilers 2003), or a quasi-Poisson
#'   additive model `s(age)` with the exposures as offset (`"gam"`; Wood
#'   2017), smoothing parameter by REML and overdispersion estimated;
#' * `dims = 1, along = "period"`: each age separately over time, the
#'   weighted Whittaker smoother of log rates with the deaths as weights
#'   (the approximate inverse variances of the log rates);
#' * `dims = 2`: the age-period surface, with a Poisson P-spline surface
#'   (`"pspline"`; Currie, Durban and Eilers 2004) or a quasi-Poisson
#'   additive model with a tensor-product P-spline `te(age, year)` and the
#'   exposures as offset (`"gam"`; Wood 2017), smoothing parameters by REML
#'   and overdispersion estimated, so that year-to-year variation beyond
#'   Poisson noise is not taken for signal.
#'
#' The Whittaker smoother minimises `sum(w (y - z)^2) + lambda sum((D z)^2)`
#' for log rates `y`, weights `w` and differences `D` of order `pord`
#' (divided differences over age, so that groups of different widths are
#' smoothed on the age scale). Its smoothing parameter is chosen for every
#' fit separately, by generalised cross-validation within `lambda_range`:
#' for each year (over age) or for each age (over time) and stratum. Over
#' time, with deaths as weights, ages whose year-to-year variation is large
#' relative to their Poisson noise (typically ages with many deaths) are
#' smoothed little or not at all, and ages with few deaths may be smoothed to
#' a straight line of the log rates. The P-spline fits choose their
#' smoothing parameters by `criterion` within `lambda_range`.
#'
#' Where a method needs one age per cell (the Whittaker smoother over age
#' and the GAMs), an open last interval such as 90+ is placed at its lower
#' bound plus half the width of the interval below it (92.5 after 85-89), a
#' convention of this package. The P-splines instead expand the open
#' interval into single years up to `omega`, as [mt_ungroup()] does.
#'
#' Smoothing over time removes year-to-year fluctuations from the fitted
#' trends and from the starting level of forecasts; the last years are
#' estimated from one side only. It also removes the year-to-year variation
#' from which models such as Lee-Carter estimate the uncertainty of their
#' period index, so their intervals become too narrow unless the spread is
#' re-estimated against the raw observations ([recalibrated()]). Smoothing
#' over age alone changes little for large populations and helps for sparse
#' data such as cause-specific rates of small regions. Whether smoothing
#' improves a forecast is a question for a backtest; the vignette
#' "Smoothing death rates" shows one. Use [mt_smooth_check()] to see how
#' much was removed and [mt_plot_smooth()] to look at it. The 2D P-spline
#' surface is slow for many ages and years (minutes per stratum); the 2D
#' GAM takes about a second.
#'
#' @param data An [mt_data] object.
#' @param dims 1 (over one direction, see `along`) or 2 (age-period
#'   surface; `"pspline"` or `"gam"`).
#' @param method `"pspline"` (Poisson P-splines, default), `"whittaker"` or
#'   `"gam"` (package \pkg{mgcv}); see Details. With `dims = 1` and
#'   `along = "period"` the default is `"whittaker"`.
#' @param along For `dims = 1`: `"age"` or `"period"` (`"whittaker"` only).
#' @param ages Range of ages that is smoothed, `c(lower, upper)`; other
#'   ages keep their observed deaths and do not enter the fit. By default
#'   all ages except the first year of life when the grid starts with age 0
#'   as its own interval (infant mortality does not follow the age pattern
#'   of the following ages).
#' @param k Basis dimensions of the tensor-product smooth for `dims = 2`,
#'   `method = "gam"`: ages and years (default: up to 10 each).
#' @param lambda Fixed smoothing parameter(s) of the P-splines or the
#'   Whittaker smoother; `NULL` selects them.
#' @param lambda_range Search range of the smoothing parameters of the
#'   P-splines and the Whittaker smoother (in both directions).
#' @inheritParams mt_ungroup
#' @return An `mt_data` object with smoothed deaths; diagnostics as attribute
#'   `"smooth"`, one row per fit, with the smoothing parameter, the effective
#'   degrees of freedom (`edf`, or `ed` for P-splines) and, for the Whittaker
#'   smoother, `at_bound` when the selected parameter is at an end of
#'   `lambda_range` before the fit reached its limit there (interpolation at
#'   the lower end, a polynomial of degree `pord - 1` at the upper end); such
#'   fits are also reported with a warning.
#' @references
#' Currie ID, Durban M, Eilers PHC (2004). Smoothing and forecasting
#' mortality rates. *Statistical Modelling* 4(4): 279-298.
#'
#' Eilers PHC (2003). A perfect smoother. *Analytical Chemistry* 75(14):
#' 3631-3636.
#'
#' Whittaker ET (1923). On a new method of graduation. *Proceedings of the
#' Edinburgh Mathematical Society* 41: 63-75.
#'
#' Wood SN (2017). *Generalized Additive Models: An Introduction with R*
#' (2nd ed.). CRC Press.
#' @examples
#' x <- mt_example(years = 2001:2003, deaths = "single")
#' s <- mt_smooth(x)
#' head(attr(s, "smooth"))
#' xt <- mt_example(deaths = "single", years = 2000:2015)
#' st <- mt_smooth(xt, along = "period", method = "whittaker")
#' mt_smooth_check(xt, st)
#' @export
mt_smooth <- function(data, dims = 1, lambda = NULL, criterion = c("bic", "aic"), kr = 2, deg = 3,
                      pord = 2, lambda_range = c(1e-2, 1e8), method = c("pspline", "whittaker", "gam"),
                      ages = NULL, along = c("age", "period"), k = NULL) {
  check_class(data, "mt_data", "data")
  criterion <- match.arg(criterion)
  along <- match.arg(along)
  # smoothing over time in one direction is Whittaker smoothing
  method <- if (missing(method) && dims == 1 && along == "period") "whittaker" else match.arg(method)
  check_number(dims, "dims", lower = 1, upper = 2, integerish = TRUE)
  if (dims == 2 && method == "whittaker") {
    mt_abort("Whittaker smoothing works in one direction; use dims = 1 with `along`, or method \"pspline\" or \"gam\".",
             class = "argument")
  }
  if (dims == 1 && along == "period" && method != "whittaker") {
    mt_abort("Smoothing over time with dims = 1 uses method = \"whittaker\"; for P-splines or GAMs use dims = 2.",
             class = "argument")
  }
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
                    deg = deg, pord = pord, lambda_range = lambda_range, method = method, ages = c(-Inf, Inf),
                    along = along, k = k)
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
    res <- if (method == "whittaker" && along == "period") {
      whittaker_time_obs(data, g, lambda, pord, lambda_range)
    } else if (method == "whittaker") {
      whittaker_smooth_obs(data, g, lambda, pord, lambda_range)
    } else if (dims == 2) {
      gam2d_smooth_obs(data, g, k)
    } else {
      gam_smooth_obs(data, g)
    }
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

# Weighted Whittaker smoother (Whittaker 1923; Eilers 2003) of one series y
# with weights w (0 for cells left out) and penalty matrix P of differences
# of order pord: the solution of (W + lambda P) z = W y, with lambda chosen
# by generalised cross-validation, n RSS / (n - ED)^2 with ED the trace of
# the hat matrix, on a log10 scale within lambda_range unless fixed.
# at_bound: the selected value is at an end of the range although the fit
# has not reached its limit there (ED = n, interpolation, at the lower end;
# ED = pord, a polynomial, at the upper end), so the range constrained it
whittaker_fit <- function(y, w, P, lambda = NULL, lambda_range = c(1e-2, 1e8), pord = 2) {
  n <- sum(w > 0)
  m <- length(y)
  fit_l <- function(l) {
    A <- diag(w, m) + l * P
    z <- solve(A, w * y)
    ed <- sum(diag(solve(A, diag(w, m))))
    list(z = z, ed = ed, gcv = n * sum(w * (y - z)^2) / max(n - ed, 1e-8)^2)
  }
  lr <- log10(lambda_range)
  ll <- if (is.null(lambda)) stats::optimize(function(v) fit_l(10^v)$gcv, lr)$minimum else NA_real_
  f <- fit_l(if (is.null(lambda)) 10^ll else lambda)
  at_bound <- !is.na(ll) && (ll < lr[1] + 0.1 && f$ed < n - 0.1 || ll > lr[2] - 0.1 && f$ed > pord + 0.1)
  list(z = f$z, lambda = if (is.null(lambda)) 10^ll else lambda, edf = f$ed, at_bound = at_bound)
}

whittaker_smooth_obs <- function(data, g, lambda, pord, lambda_range) {
  D <- mt_array(data, "deaths")
  E <- mt_array(data, "exposure")
  # divided differences at the representative ages, so that age groups of
  # different widths are smoothed on the age scale
  Dm <- divided_differences(cell_ages(g$age$lower, g$age$upper), pord)
  P <- crossprod(Dm)
  out <- D
  sm <- array(FALSE, dim(D))
  diag_rows <- list()
  for (s in seq_len(dim(D)[3])) {
    for (j in seq_len(dim(D)[2])) {
      # missing cells get weight 0: the fit uses the observed cells only
      ok <- smooth_cells(D[, j, s], E[, j, s])
      if (sum(ok) <= pord) next
      y <- ifelse(ok, log(pmax(D[, j, s], 0.5) / E[, j, s]), 0)
      w <- ifelse(ok, pmax(D[, j, s], 0.5), 0)
      f <- whittaker_fit(y, w, P, lambda, lambda_range, pord)
      out[ok, j, s] <- E[ok, j, s] * exp(f$z[ok])
      sm[ok, j, s] <- TRUE
      diag_rows[[length(diag_rows) + 1L]] <- data.table(measure = "deaths", stratum = s, period = g$period$labels[j],
                                                        method = "whittaker", lambda_age = f$lambda, edf = f$edf,
                                                        at_bound = f$at_bound)
    }
  }
  warn_whittaker_bound(diag_rows)
  list(obs = smoothed_obs(data, g, out, sm), diag = rbindlist(diag_rows))
}

warn_whittaker_bound <- function(diag_rows) {
  b <- vapply(diag_rows, function(r) isTRUE(r$at_bound), logical(1))
  if (any(b)) {
    mt_warn(sprintf(paste("The smoothing parameter is at the boundary of `lambda_range` in %d of %d Whittaker fits",
                          "before the fit reached its limit (interpolation or a polynomial); widen `lambda_range` or",
                          "fix `lambda`."), sum(b), length(b)), class = "lambda_bound", call = NULL)
  }
}

# Cells that enter a smooth: observed deaths and positive exposures
smooth_cells <- function(d, e) !is.na(d) & is.finite(e) & e > 0

# The observed deaths on the grid with the smoothed values in the cells
# marked in `sm` ([age, period, stratum]); no cells are added for missing
# deaths, and cells that were not smoothed keep their value and provenance
smoothed_obs <- function(data, g, out, sm) {
  obs <- data$obs$deaths
  idx <- obs_grid_index(obs, g, "deaths")
  in_grid <- !is.na(idx$a1)
  obs <- obs[in_grid]
  k <- cbind(idx$a1[in_grid], idx$p1[in_grid], obs$stratum)
  sel <- sm[k]
  v <- out[k]
  obs[sel, `:=`(value = v[sel], provenance = "smoothed")]
  setorderv(obs, c("stratum", "period_lo", "age_lo"))
  obs[]
}

# Quasi-Poisson GAM s(age) per stratum and year, exposures as offset, REML
gam_smooth_obs <- function(data, g) {
  need_pkg("mgcv", "for method = \"gam\"")
  D <- mt_array(data, "deaths")
  E <- mt_array(data, "exposure")
  age <- cell_ages(g$age$lower, g$age$upper)
  out <- D
  sm <- array(FALSE, dim(D))
  diag_rows <- list()
  for (s in seq_len(dim(D)[3])) {
    for (j in seq_len(dim(D)[2])) {
      ok <- smooth_cells(D[, j, s], E[, j, s])
      if (sum(ok) < 4) next
      df <- data.frame(d = D[ok, j, s], e = E[ok, j, s], age = age[ok])
      fit <- mgcv::gam(d ~ s(age, k = min(20, sum(ok) - 1)) + offset(log(e)), family = stats::quasipoisson(),
                       data = df, method = "REML")
      out[ok, j, s] <- stats::fitted(fit)
      sm[ok, j, s] <- TRUE
      diag_rows[[length(diag_rows) + 1L]] <- data.table(measure = "deaths", stratum = s, period = g$period$labels[j],
                                                        method = "gam", edf = sum(fit$edf), dispersion = fit$scale)
    }
  }
  list(obs = smoothed_obs(data, g, out, sm), diag = rbindlist(diag_rows))
}

# Weighted Whittaker smoothing over time, each age and stratum separately
# with its own smoothing parameter (generalised cross-validation); the
# deaths are the weights, as over age
whittaker_time_obs <- function(data, g, lambda, pord, lambda_range) {
  D <- mt_array(data, "deaths")
  E <- mt_array(data, "exposure")
  na <- dim(D)[1]
  np <- dim(D)[2]
  if (np <= pord + 1) mt_abort("Too few periods for smoothing over time.", class = "argument")
  P <- crossprod(diff(diag(np), differences = pord))
  out <- D
  sm <- array(FALSE, dim(D))
  diag_rows <- list()
  for (s in seq_len(dim(D)[3])) {
    for (a in seq_len(na)) {
      # missing cells get weight 0: the fit uses the observed cells only
      ok <- smooth_cells(D[a, , s], E[a, , s])
      if (sum(ok) <= pord) next
      y <- ifelse(ok, log(pmax(D[a, , s], 0.5) / E[a, , s]), 0)
      w <- ifelse(ok, pmax(D[a, , s], 0.5), 0)
      f <- whittaker_fit(y, w, P, lambda, lambda_range, pord)
      out[a, ok, s] <- E[a, ok, s] * exp(f$z[ok])
      sm[a, ok, s] <- TRUE
      diag_rows[[length(diag_rows) + 1L]] <- data.table(measure = "deaths", stratum = s, age = g$age$labels[a],
                                                        method = "whittaker", along = "period",
                                                        lambda_period = f$lambda, edf = f$edf, at_bound = f$at_bound)
    }
  }
  warn_whittaker_bound(diag_rows)
  list(obs = smoothed_obs(data, g, out, sm), diag = rbindlist(diag_rows))
}

# Quasi-Poisson tensor-product P-spline surface te(age, year) per stratum,
# exposures as offset, REML
gam2d_smooth_obs <- function(data, g, k = NULL) {
  need_pkg("mgcv", "for method = \"gam\"")
  D <- mt_array(data, "deaths")
  E <- mt_array(data, "exposure")
  na <- dim(D)[1]
  np <- dim(D)[2]
  if (na < 5 || np < 5) mt_abort("The age-period GAM needs at least 5 ages and 5 periods.", class = "argument")
  k <- k %||% c(min(10, na - 1), min(10, np - 1))
  age <- cell_ages(g$age$lower, g$age$upper)
  grid <- data.frame(age = rep(age, np), year = rep(g$period$lower, each = na))
  out <- D
  diag_rows <- list()
  for (s in seq_len(dim(D)[3])) {
    df <- cbind(grid, d = as.vector(D[, , s]), e = as.vector(E[, , s]))
    use <- !is.na(df$d) & df$e > 0
    fit <- mgcv::gam(d ~ te(age, year, bs = "ps", k = k) + offset(log(e)), family = stats::quasipoisson(),
                     data = df[use, ], method = "REML")
    pr <- exp(stats::predict(fit, newdata = df, type = "link"))
    out[, , s] <- array(pr, c(na, np))
    diag_rows[[length(diag_rows) + 1L]] <- data.table(measure = "deaths", stratum = s, method = "gam", along = "both",
                                                      edf = sum(fit$edf), dispersion = fit$scale)
  }
  obs <- data.table(stratum = rep(seq_len(dim(D)[3]), each = na * np), a = rep(seq_len(na), np * dim(D)[3]),
                    p = rep(rep(seq_len(np), each = na), dim(D)[3]), value = as.vector(out), provenance = "smoothed")
  list(obs = grid_table(obs, g), diag = rbindlist(diag_rows))
}

#' Check and plot the result of smoothing
#'
#' `mt_smooth_check()` compares smoothed with raw rates per stratum:
#'
#' * `mad_log`: median absolute difference of log rates;
#' * `dispersion`: the Pearson estimate of the dispersion of the raw deaths
#'   around the smoothed deaths, the sum of squared Pearson residuals
#'   divided by `n - edf` (McCullagh and Nelder 1989), where `n` is the
#'   number of smoothed cells with raw deaths and `edf` the effective
#'   degrees of freedom of the fits of the stratum (column `edf`, from the
#'   attribute `"smooth"` of [mt_smooth()]; when it is not available, for
#'   example after further preprocessing steps, the sum is divided by `n`
#'   and `edf` is `NA`). About 1 when smoothing removed Poisson noise only,
#'   above 1 when the raw rates vary more (for example year effects), well
#'   below 1 when the smoothing follows the noise;
#' * `trend_ratio`: variance of year-to-year changes of the smoothed log
#'   rates relative to that of the raw log rates (mean over ages); small
#'   values mean that smoothing removed most of the annual variation;
#' * `correlation`: correlation of raw and smoothed log rates.
#'
#' The measures use the smoothed cells only (provenance `"smoothed"`), so
#' ages left out of the smoothing (see `ages` in [mt_smooth()]) do not
#' dilute them.
#'
#' `mt_plot_smooth()` shows raw rates as points and smoothed rates as lines
#' over time, one panel per age.
#'
#' @param data The [mt_data] object before smoothing.
#' @param smoothed The result of [mt_smooth()] (or of [mt_prep()] steps).
#' @param ages Age labels to show (default: four ages across the range).
#' @param strata Stratum labels to show (default: all).
#' @param colour Stratum variable that gives the colours.
#' @return `mt_smooth_check()`: a `data.table` with the stratum columns, the
#'   measures above and `edf`. `mt_plot_smooth()`: a ggplot object.
#' @references
#' McCullagh P, Nelder JA (1989). *Generalized Linear Models* (2nd ed.).
#' Chapman and Hall.
#' @examples
#' x <- mt_example(deaths = "single", years = 2000:2015)
#' s <- mt_smooth(x, along = "period", method = "whittaker")
#' mt_smooth_check(x, s)
#' if (requireNamespace("ggplot2", quietly = TRUE)) mt_plot_smooth(x, s)
#' @export
mt_smooth_check <- function(data, smoothed) {
  check_class(data, "mt_data", "data")
  check_class(smoothed, "mt_data", "smoothed")
  g <- mt_grid(smoothed)
  r <- mt_set_grid(data, age = g$age, period = g$period)
  D <- mt_array(r, "deaths")
  Ds <- mt_array(smoothed, "deaths")
  E <- mt_array(smoothed, "exposure")
  if (!identical(dim(D), dim(Ds))) mt_abort("`data` and `smoothed` must have the same strata.", class = "argument")
  # only the smoothed cells (all cells when no cell is marked as smoothed)
  so <- smoothed$obs$deaths
  idx <- obs_grid_index(so, g, "deaths")
  k <- !is.na(idx$a1) & so$provenance == "smoothed"
  if (any(k)) {
    sm <- array(FALSE, dim(Ds))
    sm[cbind(idx$a1[k], idx$p1[k], so$stratum[k])] <- TRUE
    D[!sm] <- NA
    Ds[!sm] <- NA
  }
  # effective degrees of freedom of the fits, summed per stratum
  edf <- rep(NA_real_, dim(D)[3])
  dg <- attr(smoothed, "smooth")
  if (is.data.frame(dg) && nrow(dg) && "stratum" %in% names(dg)) {
    v <- dg[["edf"]] %||% dg[["ed"]]
    if (!is.null(v)) {
      sv <- tapply(v, dg$stratum, sum)
      edf[as.integer(names(sv))] <- sv
    }
  }
  rows <- lapply(seq_len(dim(D)[3]), function(s) {
    lr <- log(pmax(D[, , s], 0.5) / E[, , s])
    ls <- log(pmax(Ds[, , s], 1e-10) / E[, , s])
    lr <- matrix(lr, dim(D)[1])
    ls <- matrix(ls, dim(D)[1])
    pear <- (D[, , s] - Ds[, , s]) / sqrt(pmax(Ds[, , s], 1e-10))
    tr <- if (ncol(lr) > 2) {
      vr <- apply(lr, 1, function(v) stats::var(diff(v), na.rm = TRUE))
      vs <- apply(ls, 1, function(v) stats::var(diff(v), na.rm = TRUE))
      mean((vs / vr)[is.finite(vs / vr)])
    } else {
      NA_real_
    }
    n <- sum(!is.na(pear))
    df_res <- if (is.finite(edf[s]) && n - edf[s] > 0) n - edf[s] else n
    data.table(stratum = s, mad_log = stats::median(abs(lr - ls), na.rm = TRUE),
               dispersion = sum(pear^2, na.rm = TRUE) / df_res, trend_ratio = tr,
               correlation = suppressWarnings(stats::cor(as.vector(lr), as.vector(ls), use = "complete.obs")),
               edf = if (df_res < n) edf[s] else NA_real_)
  })
  out <- merge(smoothed$strata, rbindlist(rows), by = "stratum", sort = FALSE)
  out[]
}

#' @rdname mt_smooth_check
#' @param along `"period"` (rates over time, one panel per age) or `"age"`
#'   (rates over age, one panel per year in `periods`).
#' @param periods Years shown with `along = "age"` (default: first and last).
#' @export
mt_plot_smooth <- function(data, smoothed, ages = NULL, strata = NULL, colour = NULL, along = c("period", "age"),
                           periods = NULL) {
  need_ggplot()
  check_class(data, "mt_data", "data")
  along <- match.arg(along)
  multi <- is.list(smoothed) && !inherits(smoothed, "mt_data")
  sms <- if (multi) smoothed else list(smoothed)
  if (multi && (is.null(names(sms)) || any(!nzchar(names(sms))))) {
    mt_abort("Give several smoothed data sets as a named list, for example list(age = s1, time = s2).",
             class = "argument")
  }
  for (k in seq_along(sms)) check_class(sms[[k]], "mt_data", "smoothed")
  s1 <- sms[[1]]
  g <- mt_grid(s1)
  ages <- ages %||% g$age$labels[unique(round(seq(1, length(g$age), length.out = 4)))]
  periods <- periods %||% range(g$period$lower)
  st <- s1$strata
  strata <- check_strata_labels(strata %||% if (multi) st$label[1] else st$label, st$label)
  r <- mt_set_grid(data, age = g$age, period = g$period)
  E <- mt_array(s1, "exposure")
  dn <- list(g$age$labels, as.character(g$period$lower), st$label)
  long <- function(a) {
    dimnames(a) <- dn
    d <- array_long(a)[label %in% strata]
    d[, age_lo := g$age$lower[match(age, g$age$labels)]]
    if (along == "period") d[age %in% ages][, age_f := factor(age, levels = ages)] else d[period %in% periods]
  }
  pr <- long(mt_array(r, "deaths") / E)
  ps <- rbindlist(lapply(seq_along(sms), function(k) {
    long(mt_array(sms[[k]], "deaths") / mt_array(sms[[k]], "exposure"))[, method := names(sms)[k] %||% ""]
  }))
  xvar <- if (along == "period") "period" else "age_lo"
  if (multi) {
    ps[, method := factor(method, levels = names(sms))]
    p <- ggplot2::ggplot(ps, ggplot2::aes(x = .data[[xvar]], y = rate, colour = method,
                                          group = interaction(method, label))) +
      ggplot2::geom_point(data = pr, ggplot2::aes(group = label), colour = "grey35", size = 0.8, alpha = 0.6) +
      ggplot2::geom_line(linewidth = 0.6) +
      scale_colour_mortensor(type = "categorical")
    name <- NULL
  } else {
    lay <- plot_layout(st, strata, colour)
    pr <- merge(pr, lay$map, by = "label", sort = FALSE)
    ps <- merge(ps, lay$map, by = "label", sort = FALSE)
    p <- ggplot2::ggplot(ps, ggplot2::aes(x = .data[[xvar]], y = rate, colour = col, group = label)) +
      ggplot2::geom_point(data = pr, size = 0.7, alpha = 0.45) +
      ggplot2::geom_line(linewidth = 0.6) +
      plot_colours(lay)
    name <- lay$name
  }
  several <- length(strata) > 1 && multi
  if (along == "period") {
    p <- p + if (several) {
      ggplot2::facet_grid(ggplot2::vars(age_f), ggplot2::vars(label), scales = "free_y")
    } else if (!multi && plot_layout(st, strata, colour)$has_panel) {
      ggplot2::facet_grid(ggplot2::vars(age_f), ggplot2::vars(panel), scales = "free_y")
    } else {
      ggplot2::facet_wrap(ggplot2::vars(age_f), scales = "free_y")
    }
  } else {
    p <- p + if (several || (!multi && length(strata) > 1)) {
      ggplot2::facet_grid(ggplot2::vars(period), ggplot2::vars(label))
    } else {
      ggplot2::facet_wrap(ggplot2::vars(period))
    }
  }
  p + ggplot2::scale_y_log10() +
    ggplot2::labs(x = if (along == "period") NULL else "Age", y = "Death rate (log scale)", colour = name,
                  title = "Observed and smoothed death rates",
                  subtitle = if (multi && length(strata) == 1) sprintf("%s; points: observed; lines: smoothed", strata)
                  else "Points: observed; lines: smoothed") +
    theme_mortensor(base_size = if (length(ages) > 2) 9 else 11)
}
