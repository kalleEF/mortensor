# Spline-based disaggregation of grouped rates (ported from mortensor 0.1,
# disaggregation.R), with group means preserved exposure-weighted.

# Kannisto fit on log rates: m(x) = a e^{bx} / (1 + a e^{bx})
fit_kannisto <- function(mx, ages) {
  log_mx <- log(pmax(mx, 1e-12))
  start_b <- max(0.05, diff(range(log_mx)) / max(diff(range(ages)), 1))
  start_a <- exp(mean(log_mx) - start_b * mean(ages))
  fit <- tryCatch(stats::optim(c(log(start_a), start_b), function(par) {
    a <- exp(par[1])
    b <- par[2]
    mu <- a * exp(b * ages) / (1 + a * exp(b * ages))
    sum((log(mu) - log_mx)^2)
  }, method = "L-BFGS-B", lower = c(-30, 0.01), upper = c(5, 0.5)), error = function(e) NULL)
  if (is.null(fit)) {
    cf <- stats::coef(stats::lm(log_mx ~ ages))
    return(function(x) exp(cf[1] + cf[2] * x))
  }
  a <- exp(fit$par[1])
  b <- fit$par[2]
  function(x) a * exp(b * x) / (1 + a * exp(b * x))
}

# Group of each fine cell by containment (groups sorted by age); NA for cells
# in no group, such as cells in the gap of a missing group or above a closed
# top group
cell_groups <- function(x_lo, x_hi, g_lo, g_hi, omega) {
  grp <- findInterval(x_lo + 1e-9, g_lo)
  grp[grp == 0] <- NA
  grp[!is.na(grp) & pmin(x_hi, omega) > pmin(g_hi, omega)[grp] + 1e-9] <- NA
  grp
}

# One schedule: grouped rates -> rates on fine cells.
#   m_g: group rates; g_lo, g_hi: group bounds (Inf open); x_lo, x_hi: fine cells;
#   e: fine exposures (weights of the group-mean constraint); omega: maximum age
spline_disaggregate <- function(m_g, g_lo, g_hi, x_lo, x_hi, e, omega, method = "hyman_kannisto",
                                exact = TRUE) {
  ng <- length(m_g)
  open <- is.infinite(g_hi[ng])
  hi_f <- pmin(g_hi, omega)
  mids <- (g_lo + hi_f) / 2
  closed <- if (open) seq_len(ng - 1) else seq_len(ng)
  log_m <- log(pmax(m_g, 1e-12))
  xm <- (x_lo + pmin(x_hi, omega)) / 2
  if (method == "hyman_kannisto") {
    old <- closed[g_lo[closed] >= 60]
    if (length(old) < 3) old <- utils::tail(closed, min(4, length(closed)))
    kann <- fit_kannisto(m_g[old], mids[old])
    if (open) {
      virt <- seq(g_lo[ng] + 2.5, omega, by = 5)
      vm <- kann(virt)
      pts <- c(mids[closed], virt)
      lp <- c(log_m[closed], log(vm))
    } else {
      pts <- mids
      lp <- log_m
    }
    if (all(diff(lp) >= -0.05)) {
      for (i in seq_along(lp)[-1]) if (lp[i] <= lp[i - 1]) lp[i] <- lp[i - 1] + 1e-6
      f <- stats::splinefun(pts, lp, method = "hyman")
    } else {
      f <- stats::splinefun(pts, lp, method = "natural")
    }
  } else {
    pts <- mids[closed]
    lp <- log_m[closed]
    nc <- length(closed)
    if (nc >= 2) {
      slope <- (lp[nc] - lp[nc - 1]) / (pts[nc] - pts[nc - 1])
      step <- pts[nc] - pts[nc - 1]
      pts <- c(pts, pts[nc] + step, pts[nc] + 2 * step)
      lp <- c(lp, lp[nc] + slope * step, lp[nc] + 2 * slope * step)
    }
    f <- stats::splinefun(pts, lp, method = "natural")
  }
  m <- exp(f(xm))
  grp <- cell_groups(x_lo, x_hi, g_lo, g_hi, omega)
  w <- pmax(e, 1e-12)
  for (it in 1:15) {
    fac <- vapply(seq_len(ng), function(k) {
      i <- which(grp == k)
      if (!length(i) || m_g[k] <= 0) return(1)
      m_g[k] / (sum(w[i] * m[i]) / sum(w[i]))
    }, 0)
    if (max(abs(fac - 1)) < 1e-8) break
    cf <- stats::splinefun(mids, log(fac), method = "natural")
    m <- m * exp(cf(xm))
  }
  if (exact) {
    for (k in seq_len(ng)) {
      i <- which(grp == k)
      if (length(i) && m_g[k] > 0) m[i] <- m[i] * m_g[k] / (sum(w[i] * m[i]) / sum(w[i]))
    }
  }
  pmax(m, 1e-12)
}

# Deaths as rates by spline methods (1D per stratum and period) or 2D GAM
ungroup_rates_spline <- function(D, E, g, omega, method, tail = "kannisto", tail_width = 2) {
  idx <- obs_grid_index(D, g, "deaths")
  D <- cbind(D, idx)[!is.na(a1)]
  if (any(D$p1 != D$p2)) mt_abort("Spline ungrouping needs deaths for single periods.", class = "ungroup")
  Earr <- obs_to_matrix(E, g)
  out <- list()
  for (s in unique(D$stratum)) {
    Ds <- D[stratum == s]
    if (method == "gam") {
      res <- gam_disaggregate(Ds, matrix(Earr[, , s], dim(Earr)[1]), g, omega, tail, tail_width)
      out[[length(out) + 1L]] <- res
      next
    }
    for (p in sort(unique(Ds$p1))) {
      o <- Ds[p1 == p][order(age_lo)]
      e <- Earr[, p, s]
      Eg <- vapply(seq_len(nrow(o)), function(i) sum(e[o$a1[i]:o$a2[i]]), 0)
      m <- spline_disaggregate(o$value / pmax(Eg, 1e-12), o$age_lo, o$age_hi, g$age$lower, g$age$upper, e, omega,
                               method)
      a <- seq_along(g$age$lower)
      cov <- unlist(Map(seq, o$a1, o$a2))
      out[[length(out) + 1L]] <- data.table(stratum = s, a = a[cov], p = p, value = m[cov] * e[cov],
                                            provenance = method)
    }
  }
  grid_table(rbindlist(out), g)
}

# 2D GAM of grouped log rates (te(age, year)), Kannisto blend above the start
# of the last age group (open or wide), exposure-weighted correction of group means per year
gam_disaggregate <- function(Ds, Es, g, omega, tail = "kannisto", tail_width = 2) {
  need_pkg("mgcv", "for method = \"gam\"")
  periods <- sort(unique(Ds$p1))
  Eg <- vapply(seq_len(nrow(Ds)), function(i) sum(Es[Ds$a1[i]:Ds$a2[i], Ds$p1[i]]), 0)
  hi <- pmin(Ds$age_hi, omega)
  df <- data.frame(age = ifelse(is.infinite(Ds$age_hi), Ds$age_lo + 5, (Ds$age_lo + hi) / 2),
                   year = g$period$lower[Ds$p1], y = log(pmax(Ds$value, 0.5) / pmax(Eg, 1e-12)),
                   w = pmax(Ds$value, 0.5))
  ny <- length(unique(df$year))
  ka <- min(10, length(unique(df$age)) - 1)
  # P-spline bases need at least 4 functions; with few years an age smooth
  # plus year effects
  form <- if (ny >= 5) {
    stats::as.formula(sprintf("y ~ te(age, year, bs = \"ps\", k = c(%d, %d))", ka, min(10, ny - 1)))
  } else if (ny > 1) {
    stats::as.formula(sprintf("y ~ s(age, bs = \"ps\", k = %d) + factor(year)", ka))
  } else {
    stats::as.formula(sprintf("y ~ s(age, bs = \"ps\", k = %d)", ka))
  }
  wts <- df$w / mean(df$w)
  fit <- mgcv::gam(form, data = df, weights = wts, method = "REML")
  xm <- (g$age$lower + pmin(g$age$upper, omega)) / 2
  open_mid <- max(Ds$age_lo) + 2.5
  out <- list()
  for (p in periods) {
    o <- Ds[p1 == p][order(age_lo)]
    m <- exp(stats::predict(fit, newdata = data.frame(age = xm, year = g$period$lower[p])))
    anchor <- which(xm >= 70 & xm <= open_mid)
    if (length(anchor) >= 3) {
      tl <- if (tail == "kannisto") {
        fit_kannisto(m[anchor], xm[anchor])
      } else {
        cf <- stats::coef(stats::lm(log(m[anchor]) ~ xm[anchor]))
        function(x) exp(cf[1] + cf[2] * x)
      }
      wb <- 1 / (1 + exp(-(xm - open_mid) / tail_width))
      m <- exp((1 - wb) * log(m) + wb * log(pmax(tl(xm), 1e-12)))
    }
    e <- Es[, p]
    Eg_p <- vapply(seq_len(nrow(o)), function(i) sum(e[o$a1[i]:o$a2[i]]), 0)
    m_g <- o$value / pmax(Eg_p, 1e-12)
    grp <- cell_groups(g$age$lower, g$age$upper, o$age_lo, o$age_hi, omega)
    for (k in seq_len(nrow(o))) {
      i <- which(grp == k)
      if (length(i) && m_g[k] > 0) m[i] <- m[i] * m_g[k] / (sum(e[i] * m[i]) / sum(e[i]))
    }
    cov <- unlist(Map(seq, o$a1, o$a2))
    out[[length(out) + 1L]] <- data.table(stratum = o$stratum[1], a = cov, p = p, value = m[cov] * e[cov],
                                          provenance = "gam")
  }
  rbindlist(out)
}

#' Disaggregate forecast draws to a finer age grid
#'
#' Maps draws of age groups to finer ages (for example single years), path
#' by path, when a model was forecast at the resolution of the observed
#' deaths. Two approaches:
#'
#' * `"shape"`: age-specific factors relative to the group rate, estimated
#'   from fine-resolution rates (`shape_rates`), exposure-weighted to 1
#'   within each group. [mt_shape_rates()] estimates them from training data
#'   at the finer resolution, constant (mean of the last years) or changing
#'   over the forecast years (an age-period GAM extrapolated in time);
#' * `"hyman_kannisto"` or `"natural_gompertz"`: a spline through the group
#'   rates of each draw (monotone Hyman spline with a Kannisto tail, or a
#'   natural spline with a log-linear tail), corrected so that the
#'   exposure-weighted mean of the fine rates equals each group rate.
#'
#' @param draws An [mt_draws] object on age groups.
#' @param age The finer [mt_axis]; its cells must nest in the groups.
#' @param exposure Exposures on the finer ages, `[age, stratum]` or
#'   `[age, period, stratum]`, used as weights.
#' @param method `"shape"`, `"hyman_kannisto"` or `"natural_gompertz"`.
#' @param shape_rates For `"shape"`: rates on the finer ages, `[age, stratum]`
#'   or `[age, period, stratum]` for the forecast periods (see
#'   [mt_shape_rates()]).
#' @param omega Maximum age.
#' @return An `mt_draws` object on `age`.
#' @examples
#' x <- mt_example(years = 2000:2019, deaths = "single")
#' xg <- mt_regroup(x, age = c(seq(30, 85, 5), 90))
#' fc <- mt_forecast(mt_fit(lee_carter(), xg), h = 3, n_draws = 20, seed = 1)
#' e <- mt_array(x, "exposure")[, 20, ]
#' fs <- mt_disaggregate(fc, mt_grid(x)$age, exposure = e)
#' dim(fs$rates)
#' @export
mt_disaggregate <- function(draws, age, exposure, method = c("hyman_kannisto", "natural_gompertz", "shape"),
                            shape_rates = NULL, omega = 110) {
  check_class(draws, "mt_draws", "draws")
  check_class(age, "mt_axis", "age")
  method <- match.arg(method)
  ga <- draws$age
  d <- dim(draws$rates)
  if (length(dim(exposure)) == 2) exposure <- array(exposure, c(length(age), 1, d[3]))
  if (dim(exposure)[1] != length(age) || dim(exposure)[3] != d[3]) {
    mt_abort("`exposure` must be [age, stratum] or [age, period, stratum] on the finer ages.", class = "argument")
  }
  if (dim(exposure)[2] == 1) exposure <- exposure[, rep(1, d[2]), , drop = FALSE]
  grp <- findInterval(age$lower, ga$lower)
  if (any(grp == 0) || any(pmin(age$upper, omega) > pmin(ga$upper[pmax(grp, 1)], omega) + 1e-9)) {
    mt_abort("The finer ages do not nest in the age groups of the draws.", class = "argument")
  }
  out <- array(NA_real_, c(length(age), d[2], d[3], d[4]))
  if (method == "shape") {
    if (is.null(shape_rates)) mt_abort("`shape_rates` is needed for method = \"shape\".", class = "argument")
    if (length(dim(shape_rates)) == 2) shape_rates <- array(shape_rates, c(nrow(shape_rates), 1, ncol(shape_rates)))
    if (dim(shape_rates)[1] != length(age) || dim(shape_rates)[3] != d[3] ||
        !dim(shape_rates)[2] %in% c(1, d[2])) {
      mt_abort("`shape_rates` must be [age, stratum] or [age, period, stratum] on the finer ages and forecast years.",
               class = "argument")
    }
    for (s in seq_len(d[3])) {
      for (t in seq_len(d[2])) {
        e <- exposure[, t, s]
        r <- shape_rates[, min(t, dim(shape_rates)[2]), s]
        norm <- tapply(e * r, grp, sum) / tapply(e, grp, sum)
        nk <- as.vector(norm[as.character(grp)])
        # groups whose shape rates (or exposures) are all 0 get a flat shape
        f <- ifelse(is.finite(nk) & nk > 0, r / nk, 1)
        out[, t, s, ] <- draws$rates[grp, t, s, ] * f
      }
    }
  } else {
    for (s in seq_len(d[3])) {
      for (t in seq_len(d[2])) {
        for (j in seq_len(d[4])) {
          # missing group rates (for example a stratum whose fit failed) stay
          # missing; the spline goes through the other groups
          rg <- draws$rates[, t, s, j]
          ok <- is.finite(rg)
          if (sum(ok) < 2) next
          fr <- spline_disaggregate(rg[ok], ga$lower[ok], ga$upper[ok], age$lower, age$upper, exposure[, t, s],
                                    omega, method)
          out[, t, s, j] <- ifelse(ok[grp], fr, NA_real_)
        }
      }
    }
  }
  res <- draws
  dimnames(out) <- list(age = age$labels, period = dimnames(draws$rates)$period,
                        stratum = dimnames(draws$rates)$stratum, draw = NULL)
  res$rates <- out
  res$age <- age
  res$observed <- NULL
  res$fitted <- NULL
  res$exposure_last <- exposure[, 1, , drop = FALSE]
  dimnames(res$exposure_last) <- list(age = age$labels, period = dimnames(exposure)[[2]][1],
                                      stratum = dimnames(draws$rates)$stratum)
  res
}

#' Age shapes within groups for disaggregating forecasts
#'
#' Rates on a fine age grid that give [mt_disaggregate()] (`method =
#' "shape"`) the age pattern within each age group of a forecast:
#'
#' * `"last"`: the rates of the last `years` years (deaths and exposures
#'   summed), the same for every forecast year;
#' * `"gam"`: a quasi-Poisson age-period GAM (`te(age, year)`, exposures as
#'   offset, package \pkg{mgcv}) fitted to the last `years` years per
#'   stratum and extrapolated to `periods`, so the age pattern within groups
#'   changes over the forecast years.
#'
#' @param data An [mt_data] object at the fine age resolution (for example
#'   ungrouped training data).
#' @param periods Forecast years.
#' @param method `"last"` or `"gam"`.
#' @param years Number of last observed years used.
#' @return An array `[age, period, stratum]` of rates.
#' @examples
#' x <- mt_example(years = 2000:2019, deaths = "single")
#' sh <- mt_shape_rates(x, periods = 2020:2022)
#' xg <- mt_regroup(x, age = c(seq(30, 85, 5), 90))
#' fc <- mt_forecast(mt_fit(lee_carter(), xg), h = 3, n_draws = 20, seed = 1)
#' fs <- mt_disaggregate(fc, mt_grid(x)$age, exposure = mt_array(x, "exposure")[, 20, ],
#'                       method = "shape", shape_rates = sh)
#' @export
mt_shape_rates <- function(data, periods, method = c("last", "gam"), years = 5) {
  check_class(data, "mt_data", "data")
  method <- match.arg(method)
  check_number(years, "years", lower = 1, integerish = TRUE)
  g <- mt_grid(data)
  D <- mt_array(data, "deaths")
  E <- mt_array(data, "exposure")
  np <- dim(D)[2]
  use <- seq(max(1, np - years + 1), np)
  out <- array(NA_real_, c(dim(D)[1], length(periods), dim(D)[3]),
               dimnames = list(g$age$labels, as.character(periods), data$strata$label))
  if (method == "last") {
    r <- apply(D[, use, , drop = FALSE], c(1, 3), sum) / apply(E[, use, , drop = FALSE], c(1, 3), sum)
    out[] <- r[, rep(seq_len(dim(D)[3]), each = length(periods))]
    return(out)
  }
  need_pkg("mgcv", "for method = \"gam\"")
  if (length(use) < 5) mt_abort("The GAM shape needs at least 5 years.", class = "argument")
  age <- axis_midpoints(g$age, data$omega)
  for (s in seq_len(dim(D)[3])) {
    df <- data.frame(age = rep(age, length(use)), year = rep(g$period$lower[use], each = length(age)),
                     d = as.vector(D[, use, s]), e = as.vector(E[, use, s]))
    df <- df[!is.na(df$d) & df$e > 0, ]
    fit <- mgcv::gam(d ~ te(age, year, bs = "ps", k = c(min(10, length(age) - 1), min(5, length(use) - 1))) +
                       offset(log(e)), family = stats::quasipoisson(), data = df, method = "REML")
    nd <- data.frame(age = rep(age, length(periods)), year = rep(periods, each = length(age)), e = 1)
    out[, , s] <- exp(stats::predict(fit, newdata = nd, type = "link"))
  }
  out
}

#' Check and plot disaggregated forecasts
#'
#' `mt_disaggregate_check()` compares fine-age draws with the grouped draws
#' they come from: the largest relative difference between the
#' exposure-weighted mean of the fine rates and the group rate (0 when the
#' group totals are kept) and the largest absolute second difference of
#' the median log rates over age (large values mean kinks or oscillations).
#' `mt_plot_disaggregation()` shows, for one forecast year, the group medians
#' as steps and the fine-age medians as lines (`type = "profile"`), or the
#' absolute second differences of the fine log rates by age as a heat map
#' (`type = "smoothness"`).
#'
#' @param grouped The [mt_draws] on age groups.
#' @param fine The result of [mt_disaggregate()] for `grouped`.
#' @param exposure The exposures used by [mt_disaggregate()].
#' @param period Forecast year shown (default: the first).
#' @param strata Stratum labels shown (default: all).
#' @param type `"profile"` or `"smoothness"`.
#' @return `mt_disaggregate_check()`: a `data.table` per stratum and period
#'   with `max_group_error` and `max_second_difference`.
#'   `mt_plot_disaggregation()`: a ggplot object.
#' @examples
#' x <- mt_example(years = 2000:2019, deaths = "single")
#' xg <- mt_regroup(x, age = c(seq(30, 85, 5), 90))
#' fc <- mt_forecast(mt_fit(lee_carter(), xg), h = 3, n_draws = 20, seed = 1)
#' e <- mt_array(x, "exposure")[, 20, ]
#' fs <- mt_disaggregate(fc, mt_grid(x)$age, exposure = e)
#' mt_disaggregate_check(fc, fs, e)
#' if (requireNamespace("ggplot2", quietly = TRUE)) mt_plot_disaggregation(fc, fs)
#' @export
mt_disaggregate_check <- function(grouped, fine, exposure) {
  check_class(grouped, "mt_draws", "grouped")
  check_class(fine, "mt_draws", "fine")
  d <- dim(fine$rates)
  if (length(dim(exposure)) == 2) exposure <- array(exposure, c(d[1], 1, d[3]))
  grp <- findInterval(fine$age$lower, grouped$age$lower)
  rows <- list()
  for (s in seq_len(d[3])) {
    for (t in seq_len(d[2])) {
      e <- exposure[, min(t, dim(exposure)[2]), s]
      Fm <- matrix(fine$rates[, t, s, ], d[1])
      wm <- rowsum(Fm * e, grp) / as.vector(rowsum(e, grp))
      G <- matrix(grouped$rates[sort(unique(grp)), t, s, ], length(unique(grp)))
      med <- log(apply(Fm, 1, stats::median))
      d2 <- if (d[1] > 2) max(abs(diff(med, differences = 2))) else NA_real_
      rows[[length(rows) + 1L]] <- data.table(stratum = s, period = fine$periods[t],
                                              max_group_error = max(abs(wm / G - 1)), max_second_difference = d2)
    }
  }
  out <- merge(fine$strata, rbindlist(rows), by = "stratum", sort = FALSE)
  out[, stratum := NULL]
  out[]
}

#' @rdname mt_disaggregate_check
#' @export
mt_plot_disaggregation <- function(grouped, fine, period = NULL, strata = NULL, type = c("profile", "smoothness")) {
  need_ggplot()
  check_class(grouped, "mt_draws", "grouped")
  check_class(fine, "mt_draws", "fine")
  type <- match.arg(type)
  period <- period %||% fine$periods[1]
  t <- match(period, fine$periods)
  if (is.na(t)) mt_abort("`period` is not a forecast year.", class = "argument")
  strata <- check_strata_labels(strata %||% fine$strata$label, fine$strata$label)
  si <- match(strata, fine$strata$label)
  fm <- apply(fine$rates[, t, si, , drop = FALSE], c(1, 3), stats::median)
  fa <- fine$age
  if (type == "smoothness") {
    sd2 <- apply(log(fm), 2, function(v) c(NA, abs(diff(v, differences = 2)), NA))
    tab <- data.table(age = rep(fa$lower, length(si)), label = rep(strata, each = length(fa)), value = as.vector(sd2))
    return(ggplot2::ggplot(tab, ggplot2::aes(x = age, y = label, fill = pmin(value, 1))) +
             ggplot2::geom_tile() +
             ggplot2::scale_fill_viridis_c(option = "inferno", limits = c(0, 1), name = "|second difference|") +
             ggplot2::labs(x = "Age", y = NULL, title = "Smoothness of the disaggregated rates",
                           subtitle = sprintf("Median log rates, %s; values near 0 are smooth", period)) +
             theme_mortensor(legend = "right"))
  }
  gm <- apply(grouped$rates[, t, si, , drop = FALSE], c(1, 3), stats::median)
  ga <- grouped$age
  hi <- pmin(ga$upper, max(fa$upper[is.finite(fa$upper)], ga$lower + 5))
  gt <- data.table(label = rep(strata, each = length(ga)), lo = rep(ga$lower, length(si)),
                   hi = rep(hi, length(si)), rate = as.vector(gm))
  ft <- data.table(label = rep(strata, each = length(fa)), age = rep(fa$lower + 0.5, length(si)), rate = as.vector(fm))
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = gt, ggplot2::aes(x = lo, xend = hi, y = rate, yend = rate), colour = "grey55",
                          linewidth = 0.8) +
    ggplot2::geom_line(data = ft, ggplot2::aes(x = age, y = rate), colour = "#D7191C", linewidth = 0.6) +
    ggplot2::facet_wrap(ggplot2::vars(label)) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Age", y = "Death rate (log scale)", title = "Disaggregated forecast",
                  subtitle = sprintf("%s, medians; grey: age groups; red: finer ages", period)) +
    theme_mortensor()
}
