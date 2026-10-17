#' Smooth death rates on the grid
#'
#' Smooths death rates over age, over time or over both, with the exposures
#' as weights or offset. Deaths are replaced by exposure times the smoothed
#' rate; exposures are unchanged. Deaths and exposures must be available on
#' the grid (use [mt_ungroup()] first for grouped data). In a forecasting
#' pipeline, use it as a preprocessing step (`mt_prep(step_smooth(...))`):
#' models are then fitted to the smoothed rates, while the observation noise
#' of their predictive draws, residuals and backtest scales come from the
#' rates before smoothing ("smooth input, raw residuals").
#'
#' Methods and directions:
#'
#' * `dims = 1, along = "age"` (default): each year separately over age,
#'   with Poisson P-splines (`"pspline"`), Whittaker smoothing of log rates
#'   with the deaths as weights (`"whittaker"`), or a Gaussian additive
#'   model of log rates with the deaths as weights (`"gam"`);
#' * `dims = 1, along = "period"`: each age separately over time, Whittaker
#'   smoothing of log rates with the deaths as weights and one smoothing
#'   parameter per stratum (generalised cross-validation pooled over ages).
#'   The weights are scaled within each age, so the parameter gives the same
#'   smoothness at every age; ages with many deaths are not left
#'   unsmoothed;
#' * `dims = 2`: the age-period surface, with a Poisson P-spline surface
#'   (`"pspline"`; Currie, Durban and Eilers 2004) or a quasi-Poisson
#'   additive model with a tensor-product P-spline `te(age, year)` and the
#'   exposures as offset (`"gam"`; Wood 2017), smoothing parameters by REML
#'   and overdispersion estimated, so that year-to-year variation beyond
#'   Poisson noise is not taken for signal.
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
#' @inheritParams mt_ungroup
#' @return An `mt_data` object with smoothed deaths; diagnostics as attribute
#'   `"smooth"`.
#' @references
#' Currie ID, Durban M, Eilers PHC (2004). Smoothing and forecasting
#' mortality rates. *Statistical Modelling* 4(4): 279-298.
#'
#' Wood SN (2017). *Generalized Additive Models: An Introduction with R*
#' (2nd ed.). Chapman and Hall/CRC.
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
      whittaker_time_obs(data, g, lambda, pord)
    } else if (method == "whittaker") {
      whittaker_smooth_obs(data, g, lambda, pord)
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
  sm <- array(FALSE, dim(D))
  diag_rows <- list()
  for (s in seq_len(dim(D)[3])) {
    for (j in seq_len(dim(D)[2])) {
      # missing cells get weight 0: the fit uses the observed cells only
      ok <- smooth_cells(D[, j, s], E[, j, s])
      n_ok <- sum(ok)
      if (n_ok <= pord) next
      y <- ifelse(ok, log(pmax(D[, j, s], 0.5) / E[, j, s]), 0)
      w <- ifelse(ok, pmax(D[, j, s], 0.5), 0)
      fit_l <- function(l) {
        A <- diag(w, na) + l * P
        z <- solve(A, w * y)
        tr <- sum(diag(solve(A, diag(w, na))))
        list(z = z, gcv = sum(w * (y - z)^2) / (n_ok - tr)^2 * n_ok)
      }
      l <- lambda %||% 10^stats::optimize(function(ll) fit_l(10^ll)$gcv, c(-2, 8))$minimum
      z <- fit_l(l)$z
      out[ok, j, s] <- E[ok, j, s] * exp(z[ok])
      sm[ok, j, s] <- TRUE
      diag_rows[[length(diag_rows) + 1L]] <- data.table(measure = "deaths", stratum = s, period = g$period$labels[j],
                                                        method = "whittaker", lambda_age = l)
    }
  }
  list(obs = smoothed_obs(data, g, out, sm), diag = rbindlist(diag_rows))
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

gam_smooth_obs <- function(data, g) {
  need_pkg("mgcv", "for method = \"gam\"")
  D <- mt_array(data, "deaths")
  E <- mt_array(data, "exposure")
  na <- dim(D)[1]
  age <- axis_midpoints(g$age, data$omega)
  out <- D
  sm <- array(FALSE, dim(D))
  diag_rows <- list()
  for (s in seq_len(dim(D)[3])) {
    for (j in seq_len(dim(D)[2])) {
      ok <- smooth_cells(D[, j, s], E[, j, s])
      if (sum(ok) < 4) next
      df <- data.frame(y = log(pmax(D[ok, j, s], 0.5) / E[ok, j, s]), w = pmax(D[ok, j, s], 0.5), age = age[ok])
      wts <- df$w / mean(df$w)
      fit <- mgcv::gam(y ~ s(age, k = min(20, sum(ok) - 1)), weights = wts, data = df, method = "REML")
      out[ok, j, s] <- E[ok, j, s] * exp(stats::fitted(fit))
      sm[ok, j, s] <- TRUE
      diag_rows[[length(diag_rows) + 1L]] <- data.table(measure = "deaths", stratum = s, period = g$period$labels[j],
                                                        method = "gam", edf = sum(fit$edf))
    }
  }
  list(obs = smoothed_obs(data, g, out, sm), diag = rbindlist(diag_rows))
}

# Whittaker smoothing over time per age, one smoothing parameter per stratum
# chosen by generalised cross-validation pooled over the ages
whittaker_time_obs <- function(data, g, lambda, pord) {
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
    # matrices [age, period], also with a single age group
    Ds <- matrix(D[, , s], na, np)
    Es <- matrix(E[, , s], na, np)
    # missing cells get weight 0: the fit uses the observed cells only
    O <- smooth_cells(Ds, Es)
    Y <- ifelse(O, log(pmax(Ds, 0.5) / Es), 0)
    W <- ifelse(O, pmax(Ds, 0.5), 0)
    # deaths as weights, scaled within each age: one smoothing parameter then
    # gives the same smoothness at every age (unscaled, ages with many deaths
    # were barely smoothed and ages with few deaths smoothed strongly)
    wm <- rowSums(W) / pmax(rowSums(O), 1)
    W <- W / ifelse(wm > 0, wm, 1)
    ok <- rowSums(O) > pord
    if (!any(ok)) next
    fit_l <- function(l) {
      rss <- 0
      tr <- 0
      Z <- Y
      for (a in which(ok)) {
        A <- diag(W[a, ], np) + l * P
        z <- solve(A, W[a, ] * Y[a, ])
        Z[a, ] <- z
        rss <- rss + sum(W[a, ] * (Y[a, ] - z)^2)
        tr <- tr + sum(diag(solve(A, diag(W[a, ], np))))
      }
      n <- sum(O[ok, ])
      list(Z = Z, gcv = n * rss / (n - tr)^2)
    }
    l <- lambda %||% 10^stats::optimize(function(ll) fit_l(10^ll)$gcv, c(-2, 6))$minimum
    Z <- fit_l(l)$Z
    sel <- O & ok
    Ds[sel] <- (Es * exp(Z))[sel]
    out[, , s] <- Ds
    sm[, , s] <- sel
    diag_rows[[length(diag_rows) + 1L]] <- data.table(measure = "deaths", stratum = s, method = "whittaker",
                                                      along = "period", lambda_period = l)
  }
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
  age <- axis_midpoints(g$age, data$omega)
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
#' * `dispersion`: mean squared Pearson residual of the raw deaths around
#'   the smoothed deaths; about 1 when smoothing removed Poisson noise only,
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
#' @return `mt_smooth_check()`: a `data.table` with the stratum columns and
#'   the measures above. `mt_plot_smooth()`: a ggplot object.
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
    data.table(stratum = s, mad_log = stats::median(abs(lr - ls), na.rm = TRUE),
               dispersion = mean(pear^2, na.rm = TRUE), trend_ratio = tr,
               correlation = suppressWarnings(stats::cor(as.vector(lr), as.vector(ls), use = "complete.obs")))
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
