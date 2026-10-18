#' Lee-Carter model
#'
#' `log m[x, t] = a[x] + b[x] k[t]` with the constraints `sum(b) = 1` and
#' `sum(k) = 0`, fitted per stratum:
#'
#' * `method = "poisson"`: Poisson log-bilinear model (Brouhns, Denuit and
#'   Vermunt 2002), maximum likelihood by Newton updates;
#' * `method = "svd"`: original method (Lee and Carter 1992): singular value
#'   decomposition of centred log rates, with `k[t]` re-estimated in a
#'   second stage (`adjust`, as in `demography::lca()`): so that the fitted
#'   deaths match the observed total deaths of each year (`"dt"`, Lee and
#'   Carter 1992), by Poisson maximum likelihood of the deaths by age given
#'   `a` and `b` (`"dxt"`, Booth, Maindonald and Smith 2002), so that the
#'   fitted life expectancy at the first age matches the observed one
#'   (`"e0"`, Lee and Miller 2001), or not at all (`"none"`).
#'
#' `choose_period = TRUE` fits the model to the last years only, starting in
#' the year that minimises the ratio of the mean deviance of a fit with a
#' linear `k[t]` to the mean deviance of the Lee-Carter fit (Booth,
#' Maindonald and Smith 2002): the longest recent period over which the
#' decline was close to linear (`demography::lca(chooseperiod = TRUE,
#' breakmethod = "bms")` gives the same period). Fitted rates of earlier
#' years are missing.
#'
#' Compared with `demography::lca()`, `k[t]` is centred again after the
#' second stage (a constant shift of `k` and `a`, with identical fitted
#' rates), and life tables use `a_x = n/2` at age 0 (as [life_table()]), so
#' `adjust = "e0"` matches a slightly different life expectancy than
#' demography's life table where age 0 is in the data.
#'
#' The period index follows a random walk with drift. Draws include the
#' innovations and the uncertainty of the drift. With `bootstrap = TRUE` the
#' parameters are re-estimated on `n_boot` Poisson resamples of the fitted
#' deaths of the fitting period (semi-parametric bootstrap with fitted deaths
#' of Renshaw and Haberman 2008), which adds the uncertainty of `a` and `b`.
#' The bootstrap of Brouhns, Denuit and Van Keilegom (2005) resamples the
#' observed deaths instead: `bootstrapped(lee_carter(), type = "data")`
#' ([bootstrapped()]).
#'
#' Age groups without deaths in the fitting period get a constant rate (0.5
#' deaths over the period) and no trend; with `adjust = "e0"` they enter the
#' life tables at that rate.
#'
#' **Safeguards of this package for sparse counts.** Where the Poisson
#' likelihood has no finite maximum (`method = "poisson"` with very few
#' deaths), the fit falls back to two package-specific estimators (reference
#' unknown), each with a warning (class `mt_warning_model_fit`; within
#' [mt_fit()] it is reported with the stratum): first to the Poisson fit with
#' a ridge penalty that shrinks `b` towards equal shares `1 / n_ages`, and
#' if that also fails to the SVD estimates (`method = "svd"` with `adjust =
#' "dt"`). Their forecasts are not those of the Poisson Lee-Carter model;
#' wider age groups ([mt_regroup()]) or pooling ([li_lee()]) avoid them.
#'
#' Draws describe the expected rates. Predictive draws
#' (`mt_forecast(predictive = TRUE)`) add Poisson variation of the deaths
#' given the rates and the exposures, the observation model of the Poisson
#' Lee-Carter model (Brouhns, Denuit and Vermunt 2002); no other observation
#' noise is added. If intervals are too narrow on real data, that is a
#' property of the model: see [recalibrated()], [ensemble()], or
#' overdispersed models such as `stanmomo("lc", family = "nb")`
#' ([stanmomo()]; negative binomial Lee-Carter model of Delwarde, Denuit and
#' Partrat 2007).
#'
#' `rotate = TRUE` rotates the age pattern of decline as life expectancy
#' rises (Li, Lee and Gerland 2013). In long-lived populations the decline
#' of mortality slows down at young ages and speeds up at old ages, which a
#' fixed `b[x]` does not capture over long horizons. For each path and year
#' `t`, the life expectancy at birth `e0` of the unrotated forecast gives the
#' weight `w = (0.5 (1 + sin(pi / 2 (2 u - 1))))^power` with
#' `u = (e0 - e0_start) / (e0_end - e0_start)`; `w` is 0 below `e0_start` and
#' 1 from `e0_end`. The age pattern of that year is `B[x, t] = (1 - w) b[x] +
#' w bu[x]`. The ultimate pattern `bu` has the mean of `b` over ages 15-64 at
#' all ages below 65; from age 65 it keeps the shape of `b`, scaled so that
#' the first age group from 65 equals that mean; it sums to 1 like `b`. The
#' index is then solved again, path by path and year by year, so that the
#' rotated rates have the same life expectancy at birth as the unrotated
#' ones: rotation changes the age pattern of the forecast and keeps its life
#' expectancy at birth. Decline thus moves from young ages to the ages where
#' `bu` exceeds `b`. Where `b` at 65 is already above its mean over ages
#' 15-64 (as for French women), it moves mostly to ages 15-64. Rates change
#' from the jump-off year `T`, `log m[x, t] = log m[x, T] + B[x, t] (k*[t] -
#' k[T])`, so the jump-off rates stay as they are. Fitted rates, and paths
#' with `e0` below `e0_start`, are unchanged. In MortCast (Sevcikova et al.
#' 2016) the path of life expectancy comes from a separate model
#' (bayesLife); here it is that of the Lee-Carter draws. Rotation needs ages
#' from 0 (life expectancy at birth) with age groups starting at 15 and 65.
#' With a closed last age group, a warning notes that life expectancy is
#' then truncated at the highest age. Rotated forecasts compute life tables
#' for every path and year and are much slower (about 20 times for single
#' years of age).
#'
#' @param method `"poisson"` or `"svd"`.
#' @param adjust Re-estimation of `k[t]` for `method = "svd"`: `"dt"`,
#'   `"dxt"`, `"e0"` or `"none"`.
#' @param choose_period Choose the fitting period (Booth, Maindonald and
#'   Smith 2002).
#' @param min_years `choose_period` considers fitting periods of more than
#'   `min_years` years (`minperiod` of `demography::lca()`).
#' @param jump_off `"fitted"` (forecasts start from the fitted rates of the
#'   last year) or `"observed"` (from the observed rates of the last year;
#'   Lee and Miller 2001).
#' @param bootstrap Semi-parametric bootstrap of the parameters.
#' @param n_boot Number of bootstrap fits.
#' @param rotate Rotate the age pattern of decline `b` towards an ultimate
#'   pattern as life expectancy at birth rises (Li, Lee and Gerland 2013).
#' @param rotation Settings of the rotation, a named vector: life expectancy
#'   at birth at which the rotation starts (`e0_start`) and is complete
#'   (`e0_end`), and the exponent of the weights (`power`). The defaults are
#'   those of Li, Lee and Gerland (2013) and MortCast; elements left out keep
#'   their defaults.
#' @return A model specification.
#' @references
#' Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality.
#' *Journal of the American Statistical Association* 87(419): 659-671.
#'
#' Brouhns N, Denuit M, Vermunt JK (2002). A Poisson log-bilinear regression
#' approach to the construction of projected lifetables. *Insurance:
#' Mathematics and Economics* 31(3): 373-393.
#' \doi{10.1016/S0167-6687(02)00185-3}
#'
#' Booth H, Maindonald J, Smith L (2002). Applying Lee-Carter under
#' conditions of variable mortality decline. *Population Studies* 56(3):
#' 325-336. \doi{10.1080/00324720215935}
#'
#' Lee R, Miller T (2001). Evaluating the performance of the Lee-Carter
#' method for forecasting mortality. *Demography* 38(4): 537-549.
#' \doi{10.1353/dem.2001.0036}
#'
#' Delwarde A, Denuit M, Partrat C (2007). Negative binomial version of the
#' Lee-Carter model for mortality forecasting. *Applied Stochastic Models in
#' Business and Industry* 23(5): 385-401.
#'
#' Brouhns N, Denuit M, Van Keilegom I (2005). Bootstrapping the Poisson
#' log-bilinear model for mortality forecasting. *Scandinavian Actuarial
#' Journal* 2005(3): 212-224.
#'
#' Renshaw AE, Haberman S (2008). On simulation-based approaches to risk
#' measurement in mortality with specific reference to Poisson Lee-Carter
#' modelling. *Insurance: Mathematics and Economics* 42(2): 797-816.
#'
#' Li N, Lee R, Gerland P (2013). Extending the Lee-Carter method to model
#' the rotation of age patterns of mortality decline for long-term
#' projections. *Demography* 50(6): 2037-2051.
#' \doi{10.1007/s13524-013-0232-2}
#'
#' Sevcikova H, Li N, Kantorova V, Gerland P, Raftery AE (2016). Age-specific
#' mortality and fertility rates for probabilistic population projections.
#' In Schoen R (ed.) *Dynamic Demographic Analysis*. Springer, 285-310.
#' \doi{10.1007/978-3-319-26603-9_15}
#' @examples
#' x <- mt_example(deaths = "single")
#' fit <- mt_fit(lee_carter(), x)
#' fc <- mt_forecast(fit, h = 10, n_draws = 200, seed = 1)
#' \donttest{
#' # rotation needs ages from 0: French women (HMD, via demography)
#' if (requireNamespace("demography", quietly = TRUE)) {
#'   fr <- as_mt_data(demography::fr.mort, series = "female", years = 1970:2006)
#'   fr <- mt_regroup(fr, age = c(0, 1, seq(5, 100, 5)))
#'   fc_rot <- mt_forecast(mt_fit(lee_carter(rotate = TRUE), fr), h = 50, n_draws = 100, seed = 1)
#'   life_expectancy(fc_rot)
#' }
#' }
#' @export
lee_carter <- function(method = c("poisson", "svd"), jump_off = c("fitted", "observed"), bootstrap = FALSE,
                       n_boot = 100, adjust = c("dt", "dxt", "e0", "none"), choose_period = FALSE,
                       min_years = 20, rotate = FALSE, rotation = c(e0_start = 80, e0_end = 102, power = 0.5)) {
  method <- match.arg(method)
  jump_off <- match.arg(jump_off)
  adjust <- match.arg(adjust)
  check_flag(choose_period, "choose_period")
  check_number(min_years, "min_years", lower = 5, integerish = TRUE)
  check_flag(rotate, "rotate")
  rotation <- rotation_settings(rotation)
  new_mt_spec("lee_carter", list(method = method, jump_off = jump_off, bootstrap = bootstrap, n_boot = n_boot,
                                 adjust = adjust, choose_period = choose_period, min_years = min_years,
                                 rotate = rotate, rotation = rotation),
              min_periods = if (choose_period) min_years else 5,
              label = paste0("Lee-Carter", if (method == "svd") sprintf(" (SVD, %s)", adjust) else "",
                             if (rotate) rotation_label(rotation) else "",
                             if (choose_period) ", chosen period" else "",
                             if (jump_off == "observed") ", observed jump-off" else "",
                             if (bootstrap) ", bootstrap" else ""))
}

# Settings of the rotation: defaults completed by the elements given
rotation_settings <- function(rotation) {
  out <- c(e0_start = 80, e0_end = 102, power = 0.5)
  if (is.list(rotation)) rotation <- unlist(rotation)
  nm <- names(rotation)
  if (!is.numeric(rotation) || (length(rotation) && (is.null(nm) || !all(nm %in% names(out)) || anyDuplicated(nm)))) {
    mt_abort("`rotation` must be a named numeric vector with elements `e0_start`, `e0_end` and `power`.",
             class = "argument", call = NULL)
  }
  out[nm] <- rotation
  if (anyNA(out) || any(!is.finite(out)) || out[["e0_end"]] <= out[["e0_start"]] || out[["power"]] <= 0) {
    mt_abort("`rotation` needs finite values with `e0_start` < `e0_end` and `power` > 0.", class = "argument",
             call = NULL)
  }
  out
}

rotation_label <- function(rotation) {
  default <- c(e0_start = 80, e0_end = 102, power = 0.5)
  if (isTRUE(all.equal(rotation, default))) return(", rotated")
  sprintf(", rotated (e0 %g-%g, power %g)", rotation[["e0_start"]], rotation[["e0_end"]], rotation[["power"]])
}

#' @export
fit_model.mt_spec_lee_carter <- function(spec, D, E, ages, periods, ...) {
  p <- spec$params
  if (isTRUE(p$rotate)) rotation_check(ages)
  zero <- rowSums(D) == 0
  if (sum(!zero) < 2) mt_abort("Fewer than two age groups with deaths.", class = "model_fit")
  fallback <- character()
  # Poisson estimates; where the Poisson likelihood has no finite maximum
  # (very sparse counts) the package's safeguards: a ridge penalty on b, then
  # SVD estimates (`fallback` records which). `fixed`: constant rates of the
  # age groups left out of D (NA for the rows of D), for life tables over all
  # ages
  est_one <- function(D, E, fixed = NULL) {
    if (p$method == "svd") return(lc_svd(D, E, adjust = p$adjust %||% "dt", ages = ages, fixed = fixed))
    tryCatch(lc_poisson(D, E), mt_error_model_fit = function(e) {
      tryCatch({
        out <- lc_poisson(D, E, ridge = 1)
        fallback <<- c(fallback, "ridge")
        out
      }, mt_error_model_fit = function(e) {
        fallback <<- c(fallback, "svd")
        lc_svd(D, E)
      })
    })
  }
  fit_lc <- function(D, E) {
    if (!any(zero)) return(est_one(D, E))
    # age groups without any death: constant rate of 0.5 deaths over the
    # training years, no trend (the Poisson estimate would be minus infinity)
    fixed <- rep(NA_real_, nrow(D))
    fixed[zero] <- 0.5 / rowSums(E[zero, , drop = FALSE])
    e <- est_one(D[!zero, , drop = FALSE], E[!zero, , drop = FALSE], fixed)
    a <- b <- numeric(nrow(D))
    a[!zero] <- e$a
    b[!zero] <- e$b
    a[zero] <- log(fixed[zero])
    e$a <- a
    e$b <- b
    e
  }
  first <- 1L
  if (isTRUE(p$choose_period)) {
    first <- bms_start(D, E, est_one, p$min_years)
    D <- D[, first:ncol(D), drop = FALSE]
    E <- E[, first:ncol(E), drop = FALSE]
    zero <- rowSums(D) == 0
  }
  # fallbacks of the final fit only (not of the candidate periods)
  fallback <- character()
  est <- fit_lc(D, E)
  if (any(zero)) {
    mt_warn(sprintf(paste("No deaths in age group(s) %s in the training years: constant rate (0.5 deaths in all",
                          "years), no trend."), paste(ages$labels[zero], collapse = ", ")), class = "sparse")
  }
  if ("ridge" %in% fallback) {
    mt_warn(paste("The Poisson likelihood of the Lee-Carter model has no finite maximum (sparse counts); fitted",
                  "with a ridge penalty that shrinks the age pattern of change b towards equal shares (a",
                  "safeguard of this package, not the Poisson Lee-Carter model)."), class = "model_fit")
  }
  if ("svd" %in% fallback) {
    mt_warn(paste("The Poisson likelihood of the Lee-Carter model has no finite maximum (sparse counts), also with",
                  "a ridge penalty; the SVD estimates (method = \"svd\", adjust = \"dt\") are used instead (a",
                  "safeguard of this package)."), class = "model_fit")
  }
  est$last_log_rate <- log_rates(D[, ncol(D)], E[, ncol(E)])
  est$fitted <- exp(est$a + outer(est$b, est$k))
  est$jump_off <- p$jump_off
  if (isTRUE(p$bootstrap)) {
    # resamples of the (chosen) fitting period
    fitted <- E * est$fitted
    est$boot <- lapply(seq_len(p$n_boot), function(i) {
      Db <- matrix(stats::rpois(length(fitted), fitted), nrow(fitted))
      tryCatch(fit_lc(Db, E), error = function(e) NULL)
    })
    est$boot <- Filter(Negate(is.null), est$boot)
  }
  if (isTRUE(p$rotate)) {
    est$rotation <- c(list(age = ages), as.list(p$rotation))
    est$b_ultimate <- lc_ultimate_b(est$b, ages, zero)
    if (is.null(est$b_ultimate)) {
      mt_abort(paste("Rotation needs a positive age pattern of decline b at ages 15-64 (mean) and in the first age",
                     "group from 65."), class = "model_fit")
    }
    # bootstrap fits without a valid ultimate pattern are dropped
    if (length(est$boot)) {
      est$boot <- Filter(Negate(is.null), lapply(est$boot, function(s) {
        s$b_ultimate <- lc_ultimate_b(s$b, ages, zero)
        if (is.null(s$b_ultimate)) NULL else s
      }))
    }
  }
  if (first > 1L) {
    # fitted rates only for the chosen period
    est$fitted <- cbind(matrix(NA_real_, nrow(D), first - 1L), est$fitted)
    est$first_year <- periods[first]
  }
  est
}

#' @export
simulate_model.mt_spec_lee_carter <- function(spec, fit, h, n_draws, ...) {
  sets <- if (length(fit$boot)) fit$boot else list(fit)
  which_set <- rep_len(seq_along(sets), n_draws)
  out <- array(NA_real_, c(length(fit$a), h, n_draws))
  for (i in seq_along(sets)) {
    js <- which(which_set == i)
    if (!length(js)) next
    s <- sets[[i]]
    kp <- rw_paths(s$k[length(s$k)], s$drift, matrix(s$dk - s$drift, nrow = 1), h, length(js), "gaussian",
                   sigma = s$sigma_k, exp_out = FALSE)
    for (m in seq_along(js)) {
      kk <- kp[1, , m]
      base <- if (fit$jump_off == "observed") fit$last_log_rate - s$b * s$k[length(s$k)] else s$a
      out[, , js[m]] <- exp(base + outer(s$b, kk))
    }
    if (!is.null(fit$rotation)) out[, , js] <- lc_rotate(out[, , js, drop = FALSE], matrix(kp[1, , ], h), s, fit)
  }
  out
}

# Rotation of the age pattern of mortality decline (Li, Lee and Gerland
# 2013). rotation_weights() and lc_ultimate_b() are adapted from MortCast
# 2.8-0 (rotate.leecarter() and ultimate.bx() in R/LC.R), GPL (>= 2), by Hana
# Sevcikova, Nan Li and Patrick Gerland. Changes: the ultimate pattern is
# defined on any age axis with groups starting at 0, 15 and 65 (the mean
# over ages 15-64 is weighted by the widths of the age groups, which gives
# MortCast's mean for 5-year groups and single years), age groups without
# deaths keep no trend, and the sum is that of b. In MortCast the index is
# solved by bisection (to 0.01 years) for an external path of life
# expectancy; here rotate_rates() solves it for the life expectancy of each
# unrotated Lee-Carter path, and rates change from the jump-off year.

# Ages from 0 (life expectancy at birth) with age groups starting at 15 and
# 65 (the ultimate pattern)
rotation_check <- function(age) {
  lo <- age$lower
  if (lo[1] != 0 || !all(c(15, 65) %in% lo)) {
    mt_abort(sprintf(paste("lee_carter(rotate = TRUE) needs ages from 0 with age groups starting at 15 and 65",
                           "(life expectancy at birth, ultimate pattern of decline); the grid has ages %s to %s.",
                           "Extend or regroup the ages with mt_regroup()."), age$labels[1], age$labels[length(age)]),
             class = "model_data", call = NULL)
  }
  if (is.finite(age$upper[length(age)])) {
    mt_warn(sprintf(paste("The last age group (%s) is closed: the rotation uses life expectancy up to age %s, below",
                          "life expectancy at birth."), age$labels[length(age)], format(age$upper[length(age)])),
            class = "model_data", call = NULL)
  }
  invisible(age)
}

# Weights of the ultimate pattern for life expectancies e0 (MortCast's
# rotate.leecarter(): 0 below e0_start, 1 from e0_end)
rotation_weights <- function(e0, e0_start = 80, e0_end = 102, power = 0.5) {
  u <- (e0 - e0_start) / (e0_end - e0_start)
  w <- (0.5 * (1 + sin(pi / 2 * (2 * u - 1))))^power
  w[!is.na(e0) & e0 < e0_start] <- 0
  w[!is.na(e0) & e0 >= e0_end] <- 1
  w
}

# Ultimate age pattern of decline (MortCast's ultimate.bx()): the mean of b
# over ages 15-64 at ages below 65; from 65 the shape of b, scaled so that
# the first age group from 65 equals that mean; the same sum as b. Age groups
# without deaths (`fixed`) keep b = 0. NULL where the mean or b at 65 is not
# positive.
lc_ultimate_b <- function(b, age, fixed = NULL) {
  fixed <- fixed %||% rep(FALSE, length(b))
  lo <- age$lower
  n <- age$upper - lo
  young <- lo >= 15 & age$upper <= 65 & !fixed
  i65 <- which(lo == 65)
  if (!any(young) || length(i65) != 1L) return(NULL)
  m <- sum(b[young] * n[young]) / sum(n[young])
  if (!isTRUE(m > 0) || !isTRUE(b[i65] > 0)) return(NULL)
  bu <- ifelse(lo >= 65, b * (m / b[i65]), m)
  bu[fixed] <- 0
  bu / sum(bu) * sum(b)
}

# Rotated rates for the unrotated draws R [age, h, draw] of parameter set s,
# with index paths K [h, draw]. Rates change from the jump-off year: log m =
# J + B (k* - k[T]), J the jump-off log rates (fitted or observed)
lc_rotate <- function(R, K, s, fit) {
  kj <- s$k[length(s$k)]
  J <- if (fit$jump_off == "observed") fit$last_log_rate else s$a + s$b * kj
  d <- dim(R)
  array(rotate_rates(matrix(R, d[1]), J, s$b, s$b_ultimate, as.vector(K) - kj, fit$rotation), d)
}

# Rotation of unrotated rates M [age, path-year] = exp(J + b kap0): per
# column, the weight of the ultimate pattern bu from the life expectancy at
# birth of M, and the index solved so that exp(J + B kap) keeps that life
# expectancy. Columns below e0_start are returned as they are. `rot` holds
# the age axis and the settings (fit$rotation).
rotate_rates <- function(M, J, b, bu, kap0, rot) {
  e0 <- life_expectancy(M, rot$age)
  w <- rotation_weights(e0, rot$e0_start, rot$e0_end, rot$power)
  act <- which(!is.na(w) & w > 0)
  if (!length(act)) return(M)
  B <- outer(b, 1 - w[act]) + outer(bu, w[act])
  kap <- lc_solve_e0(J, B, e0[act], kap0[act], rot$age)
  ok <- !is.na(kap)
  if (!all(ok)) {
    mt_warn(sprintf(paste("Rotation: no index matches the life expectancy of %d of %d path-years; their rates are",
                          "not rotated."), sum(!ok), length(ok)), class = "model", call = NULL)
  }
  M[, act[ok]] <- exp(J + B[, ok, drop = FALSE] * rep(kap[ok], each = nrow(M)))
  M
}

# Index kap (one per column of B) such that the life expectancy at birth of
# exp(J + B[, j] kap[j]) equals e0[j]: secant steps for all columns at once
# from `start`, uniroot() for the columns that do not converge; NA where no
# root is found
lc_solve_e0 <- function(J, B, e0, start, age, tol = 1e-10, maxit = 30) {
  na <- length(J)
  f <- function(kap, cols) {
    life_expectancy(exp(J + B[, cols, drop = FALSE] * rep(kap, each = na)), age) - e0[cols]
  }
  out <- rep(NA_real_, length(start))
  x0 <- start
  f0 <- f(x0, seq_along(x0))
  conv <- is.finite(f0) & abs(f0) < tol
  out[conv] <- x0[conv]
  x1 <- x0 - 0.01 * (1 + abs(x0))
  act <- which(!conv & is.finite(f0))
  for (it in seq_len(maxit)) {
    if (!length(act)) break
    xa <- x1[act]
    f1 <- f(xa, act)
    conv <- is.finite(f1) & abs(f1) < tol
    out[act[conv]] <- xa[conv]
    step <- f1 * (xa - x0[act]) / (f1 - f0[act])
    x0[act] <- xa
    f0[act] <- f1
    x1[act] <- xa - step
    act <- act[!conv & is.finite(step)]
  }
  for (j in which(is.na(out))) {
    out[j] <- tryCatch(stats::uniroot(function(kap) f(kap, j), start[j] + c(-1, 1) * (1 + abs(start[j])),
                                      extendInt = "yes", tol = 1e-12)$root,
                       error = function(e) NA_real_)
  }
  out
}

lc_finish <- function(a, b, k) {
  s <- sum(b)
  kbar <- mean(k)
  a <- a + b * kbar
  k <- (k - kbar) * s
  b <- b / s
  dk <- diff(k)
  drift <- mean(dk)
  list(a = a, b = b, k = k, dk = dk, drift = drift, sigma_k = sqrt(sum((dk - drift)^2) / max(length(dk) - 1, 1)))
}

# `fixed`: for adjust = "e0" with rows of `ages` left out of D, their constant
# rates (NA for the rows of D), so that life tables cover the whole age range
lc_svd <- function(D, E, adjust = "dt", ages = NULL, fixed = NULL) {
  y <- log_rates(D, E)
  a <- rowMeans(y)
  s <- svd(y - a, nu = 1, nv = 1)
  b <- s$u[, 1]
  k <- s$d[1] * s$v[, 1]
  if (sum(b) < 0) {
    b <- -b
    k <- -k
  }
  # second stage: re-estimate k (Lee and Carter: total deaths; Booth et al.:
  # deaths by age; Lee and Miller: life expectancy)
  if (adjust == "e0" && is.null(ages)) adjust <- "dt"
  all_ages <- function(m) {
    if (is.null(fixed)) return(array(m, c(nrow(m), ncol(m), 1)))
    out <- matrix(fixed, length(fixed), ncol(m))
    out[is.na(fixed), ] <- m
    array(out, c(nrow(out), ncol(out), 1))
  }
  if (adjust == "e0") e_obs <- life_expectancy(all_ages(pmax(D, 0.5) / E), ages)[, 1]
  for (t in seq_len(ncol(D))) {
    if (adjust == "none") break
    if (adjust == "dxt") {
      kt <- k[t]
      for (it in 1:50) {
        mu <- E[, t] * exp(a + b * kt)
        step <- sum((D[, t] - mu) * b) / sum(mu * b^2)
        kt <- kt + step
        if (abs(step) < 1e-10) break
      }
      if (is.finite(kt)) k[t] <- kt
      next
    }
    f <- if (adjust == "dt") {
      function(kt) sum(E[, t] * exp(a + b * kt)) - sum(D[, t])
    } else {
      function(kt) life_expectancy(all_ages(matrix(exp(a + b * kt))), ages)[1, 1] - e_obs[t]
    }
    rng <- k[t] + c(-1, 1) * (10 + 10 * abs(k[t]))
    # without a root in reach, total deaths keep the SVD index (with `b` of
    # both signs the nearest root of total deaths can lie far off and makes
    # the index erratic); life expectancy searches for the nearest root
    k[t] <- tryCatch(stats::uniroot(f, rng, extendInt = "yes", tol = 1e-10)$root, error = function(e) {
      if (adjust == "e0") nearest_root(f, k[t], 0.1 * (1 + abs(k[t]))) else k[t]
    })
  }
  out <- lc_finish(a, b, k)
  out$method <- "svd"
  out
}

# Root of f nearest to x0, for f that is not monotone (life expectancy as a
# function of k when b has both signs, with sparse counts): steps of growing
# size in both directions until the sign changes; x0 if there is none
nearest_root <- function(f, x0, step) {
  f0 <- f(x0)
  if (!is.finite(f0) || f0 == 0) return(x0)
  prev <- c(x0, x0)
  for (d in step * 2^(-4:12)) {
    for (i in 1:2) {
      x <- x0 + c(-1, 1)[i] * d
      fx <- f(x)
      if (is.finite(fx) && sign(fx) != sign(f0)) {
        return(tryCatch(stats::uniroot(f, sort(c(prev[i], x)), tol = 1e-10)$root, error = function(e) x0))
      }
      prev[i] <- x
    }
  }
  x0
}

lc_poisson <- function(D, E, init = NULL, maxit = 500, tol = 1e-8, ridge = 0) {
  if (is.null(init)) init <- lc_svd(D, E)
  a <- init$a
  b <- if (ridge > 0) rep(1 / nrow(D), nrow(D)) else init$b
  k <- init$k
  b0 <- 1 / nrow(D)
  ll_old <- -Inf
  for (it in seq_len(maxit)) {
    Dh <- E * exp(a + outer(b, k))
    a <- a + rowSums(D - Dh) / rowSums(Dh)
    Dh <- E * exp(a + outer(b, k))
    k <- k + colSums((D - Dh) * b) / colSums(Dh * b^2)
    a <- a + b * mean(k)
    k <- k - mean(k)
    Dh <- E * exp(a + outer(b, k))
    kk <- matrix(k, nrow(D), ncol(D), byrow = TRUE)
    # a ridge penalty shrinks b towards equal shares where the counts are
    # too sparse for a finite maximum likelihood estimate
    b <- b + (rowSums((D - Dh) * kk) - ridge * (b - b0)) / (rowSums(Dh * kk^2) + ridge)
    Dh <- E * exp(a + outer(b, k))
    ll <- sum(D * log(Dh) - Dh)
    if (!is.finite(ll)) mt_abort("Lee-Carter fit diverged.", class = "model_fit")
    if (abs(ll - ll_old) < tol * abs(ll)) break
    ll_old <- ll
  }
  out <- lc_finish(a, b, k)
  out$method <- "poisson"
  out$loglik <- ll
  out$iterations <- it
  out
}

# First year of the fitting period with the smallest ratio of the mean
# deviance of a linear-k fit to the mean deviance of the Lee-Carter fit
# (Booth, Maindonald and Smith 2002)
bms_start <- function(D, E, est_one, min_years) {
  nt <- ncol(D)
  na <- nrow(D)
  if (nt <= min_years + 1) return(1L)
  starts <- seq_len(nt - min_years)
  ratio <- vapply(starts, function(s) {
    cols <- s:nt
    Ds <- D[, cols, drop = FALSE]
    Es <- E[, cols, drop = FALSE]
    est <- tryCatch(est_one(Ds, Es), error = function(e) NULL)
    if (is.null(est)) return(Inf)
    n <- length(cols)
    # linear k through its mean with the mean annual change as slope
    klin <- mean(est$k) + mean(diff(est$k)) * (seq_len(n) - (n + 1) / 2)
    dev_base <- poisson_dev(Ds, Es * exp(est$a + outer(est$b, est$k)))
    dev_lin <- poisson_dev(Ds, Es * exp(est$a + outer(est$b, klin)))
    (dev_lin / (na * (n - 2))) / (dev_base / ((na - 1) * (n - 2)))
  }, 0)
  starts[which.min(ratio)]
}
