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
  grp <- findInterval(x_lo, g_lo)
  grp[grp == 0] <- NA
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
ungroup_rates_spline <- function(D, E, g, omega, method) {
  idx <- obs_grid_index(D, g, "deaths")
  D <- cbind(D, idx)[!is.na(a1)]
  if (any(D$p1 != D$p2)) mt_abort("Spline ungrouping needs deaths for single periods.", class = "ungroup")
  Earr <- obs_to_matrix(E, g)
  out <- list()
  for (s in unique(D$stratum)) {
    Ds <- D[stratum == s]
    if (method == "gam") {
      res <- gam_disaggregate(Ds, Earr[, , s], g, omega)
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
gam_disaggregate <- function(Ds, Es, g, omega) {
  need_pkg("mgcv", "for method = \"gam\"")
  periods <- sort(unique(Ds$p1))
  Eg <- vapply(seq_len(nrow(Ds)), function(i) sum(Es[Ds$a1[i]:Ds$a2[i], Ds$p1[i]]), 0)
  hi <- pmin(Ds$age_hi, omega)
  df <- data.frame(age = ifelse(is.infinite(Ds$age_hi), Ds$age_lo + 5, (Ds$age_lo + hi) / 2),
                   year = g$period$lower[Ds$p1], y = log(pmax(Ds$value, 0.5) / pmax(Eg, 1e-12)),
                   w = pmax(Ds$value, 0.5))
  form <- stats::as.formula(sprintf("y ~ te(age, year, bs = \"ps\", k = c(%d, %d))",
                                    min(10, length(unique(df$age)) - 1), min(10, length(unique(df$year)) - 1)))
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
      kann <- fit_kannisto(m[anchor], xm[anchor])
      wb <- 1 / (1 + exp(-(xm - open_mid) / 2))
      m <- exp((1 - wb) * log(m) + wb * log(pmax(kann(xm), 1e-12)))
    }
    e <- Es[, p]
    Eg_p <- vapply(seq_len(nrow(o)), function(i) sum(e[o$a1[i]:o$a2[i]]), 0)
    m_g <- o$value / pmax(Eg_p, 1e-12)
    grp <- findInterval(g$age$lower, o$age_lo)
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
#'   from fine-resolution rates (`shape_rates`, for example ungrouped training
#'   data of the last years), exposure-weighted to 1 within each group;
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
#' @param shape_rates For `"shape"`: rates on the finer ages, `[age, stratum]`.
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
    for (s in seq_len(d[3])) {
      for (t in seq_len(d[2])) {
        e <- exposure[, t, s]
        r <- shape_rates[, s]
        norm <- tapply(e * r, grp, sum) / tapply(e, grp, sum)
        f <- as.vector(r / norm[as.character(grp)])
        out[, t, s, ] <- draws$rates[grp, t, s, ] * f
      }
    }
  } else {
    for (s in seq_len(d[3])) {
      for (t in seq_len(d[2])) {
        for (j in seq_len(d[4])) {
          out[, t, s, j] <- spline_disaggregate(draws$rates[, t, s, j], ga$lower, ga$upper, age$lower, age$upper,
                                                exposure[, t, s], omega, method)
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
  res$exposure_last <- exposure[, 1, , drop = FALSE]
  res
}
