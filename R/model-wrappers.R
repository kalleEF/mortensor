#' Functional demographic model of Hyndman and Ullah (demography)
#'
#' Wrapper around `demography::fdm()` (Hyndman and Ullah 2007): the rates of
#' each year are smoothed over age with `demography::smooth.demogdata()`
#' (`smooth = TRUE`, by default monotonic regression splines above age 65,
#' weighted by the exposures), the log rates are decomposed into a mean
#' function and principal components, and the scores are forecast with
#' ARIMA or ETS models. Draws are simulated with `demography`'s `simulate()`
#' method for functional forecasts (uncertainty of the coefficients and
#' resampled model errors). Predictive draws (`mt_forecast(predictive =
#' TRUE)`) add Poisson variation of the deaths given the rates and the
#' exposures, as for other models of expected rates, and no other
#' observation noise. If intervals are too narrow on real data, that is a
#' property of the model (see [recalibrated()], [ensemble()], or
#' overdispersed models such as `stanmomo(family = "nb")`).
#'
#' Without smoothing, the model errors that `simulate()` resamples contain
#' the observation noise, which predictive draws then add a second time (as
#' Poisson variation): keep `smooth = TRUE` unless the rates were smoothed
#' in a preprocessing step. `smooth.demogdata()` works on single years of age; for age groups
#' the smoothed values at the observed ages are used (a weaker smoothing).
#' The smoothing settings (`smooth_args`, such as `k`, the number of knots)
#' and the number of components can be chosen by forecast accuracy with
#' [tuned()]. When the default monotonic spline cannot be fitted (for
#' example for age ranges that start in middle age), the rates are smoothed
#' without the monotonicity constraint and the fit diagnostics say so.
#'
#' @param order Number of principal components.
#' @param method Forecast method for the coefficients (`"arima"`, `"ets"`,
#'   `"rwdrift"`, ...), passed to `forecast()`.
#' @param smooth Smooth the rates of each year over age before the
#'   decomposition.
#' @param smooth_args List of arguments of `demography::smooth.demogdata()`
#'   (for example `list(k = 20)` or `list(method = "spline")`).
#' @param ... Further arguments of `demography::fdm()`.
#' @return A model specification.
#' @references
#' Hyndman RJ, Ullah MS (2007). Robust forecasting of mortality and fertility
#' rates: a functional data approach. *Computational Statistics & Data
#' Analysis* 51(10): 4942-4956.
#' @examples
#' if (requireNamespace("demography", quietly = TRUE)) {
#'   x <- mt_example(deaths = "single", years = 1995:2019)
#'   fc <- mt_forecast(mt_fit(hyndman_ullah(order = 2), x), h = 5, n_draws = 50, seed = 1)
#' }
#' @export
hyndman_ullah <- function(order = 6, method = "arima", smooth = TRUE, smooth_args = list(), ...) {
  need_pkg("demography", "for hyndman_ullah()")
  check_flag(smooth, "smooth")
  if (!is.list(smooth_args)) mt_abort("`smooth_args` must be a list.", class = "argument")
  new_mt_spec("fdm", list(order = order, method = method, smooth = smooth, smooth_args = smooth_args,
                          args = list(...)), min_periods = 8,
              label = paste0("Functional demographic model", if (!smooth) " (unsmoothed)"))
}

#' @export
fit_model.mt_spec_fdm <- function(spec, D, E, ages, periods, ...) {
  p <- spec$params
  rates <- pmax(D, 0.5) / E
  dimnames(rates) <- dimnames(E) <- list(ages$lower, periods)
  dd <- demography::demogdata(data = rates, pop = E, ages = ages$lower, years = periods, type = "mortality",
                              label = "mortensor", name = "s")
  if (!isFALSE(p$smooth)) dd <- fdm_smooth(dd, ages$lower, p$smooth_args)
  fit <- suppressWarnings(do.call(demography::fdm, c(list(dd, series = "s", order = min(p$order, nrow(D) - 1)),
                                                     p$args)))
  fv <- fit$fitted$y
  if (is.list(fv)) fv <- fv[[1]]
  list(fit = fit, method = p$method, fitted = exp(as.matrix(fv)))
}

# Smooth the rates of each year over age with demography::smooth.demogdata()
# and keep the observed ages (the smoother returns single years of age)
fdm_smooth <- function(dd, ages, args) {
  try_smooth <- function(a) {
    tryCatch(suppressWarnings(do.call(demography::smooth.demogdata, c(list(dd), a))), error = function(e) e)
  }
  sm <- try_smooth(args)
  if (inherits(sm, "error") && is.null(args$b) && is.null(args$method)) {
    # the monotonic spline above age b can fail to start; smooth without
    # the constraint (a penalised regression spline)
    first <- conditionMessage(sm)
    sm <- try_smooth(c(args, list(b = 1e3)))
    if (!inherits(sm, "error")) {
      warning(sprintf("smooth.demogdata() failed (%s); smoothed without the monotonicity constraint.", first),
              call. = FALSE)
    }
  }
  if (inherits(sm, "error")) {
    mt_abort(paste("smooth.demogdata() failed:", conditionMessage(sm), "Try other `smooth_args` or `smooth = FALSE`."),
             class = "model", call = NULL)
  }
  keep <- match(ages, sm$age)
  if (anyNA(keep)) {
    mt_abort("smooth.demogdata() did not return every observed age.", class = "model", call = NULL)
  }
  r <- as.matrix(sm$rate[[1]])[keep, , drop = FALSE]
  r[!is.finite(r) | r <= 0] <- as.matrix(dd$rate[[1]])[!is.finite(r) | r <= 0]
  dimnames(r) <- dimnames(dd$rate[[1]])
  dd$rate[[1]] <- r
  dd
}

#' @export
simulate_model.mt_spec_fdm <- function(spec, fit, h, n_draws, ...) {
  fc <- suppressWarnings(forecast::forecast(fit$fit, h = h, method = fit$method))
  sims <- stats::simulate(fc, nsim = n_draws)
  array(as.numeric(sims), c(dim(sims)[1], h, n_draws))
}

#' Generalised age-period-cohort models (StMoMo)
#'
#' Wrapper around \pkg{StMoMo} (Villegas, Kaishev and Millossovich 2018):
#' age-period-cohort (`"apc"`), Cairns-Blake-Dowd (`"cbd"`, logit link),
#' Renshaw-Haberman (`"rh"`), M6, M7 and M8 (`"m6"`, `"m7"`, `"m8"`, logit
#' link), the Plat model (`"plat"`, logit link) and Lee-Carter (`"lc"`). `"rh"` is StMoMo's Renshaw-Haberman model
#' with `cohortAgeFun = "1"` (the default of StMoMo's `rh()`): `log m(x,t) =
#' a(x) + b(x) k(t) + g(t - x)`, a cohort effect without age modulation, not
#' the original 2006 specification of the Renshaw-Haberman model, whose
#' cohort effect is modulated by a second age function. Models are fitted per stratum with real
#' ages (the lower bounds of single-year intervals; models with a cohort
#' effect need single years of age), deaths and central exposures; for the
#' logit-link models initial exposures are approximated as central exposures
#' plus half the deaths, and forecast death probabilities are converted to
#' central rates (`m = -log(1 - q)`). Cohorts observed in fewer than `clip`
#' cells get zero weight; their fitted rates use the effect of the nearest
#' estimated cohort. Draws come from `simulate()` of the fitted model (period
#' indexes simulated as a multivariate random walk with drift, the cohort
#' index as an ARIMA process). Predictive draws (`mt_forecast(predictive =
#' TRUE)`) add Poisson variation of the deaths given the rates and the
#' exposures, and no other observation noise.
#'
#' M6, M7 and M8 are the cohort extensions of the CBD model of Cairns et al.
#' (2009), with `x` the age, `xbar` the mean age of the data and `s2` the
#' mean of `(x - xbar)^2`:
#'
#' * M6: `logit q(x,t) = k1(t) + (x - xbar) k2(t) + g(t - x)`.
#' * M7: M6 plus a quadratic age term `((x - xbar)^2 - s2) k3(t)`; the cohort
#'   effect has no linear or quadratic trend (constraints of StMoMo's
#'   `m7()`).
#' * M8: `logit q(x,t) = k1(t) + (x - xbar) k2(t) + (xc - x) g(t - x)`: the
#'   cohort effect fades with age and vanishes at age `xc`.
#'
#' The CBD family is meant for older ages (StMoMo's examples use ages 55 to
#' 89). For M7 the cohort index is forecast with ARIMA(2,0,0) with a non-zero
#' mean, as in the StMoMo vignette, since its constraints remove any trend;
#' the other cohort models use ARIMA(1,1,0) with drift, except the Plat
#' model.
#'
#' **Plat model.** Plat (2009) combines the CBD model with features of the
#' Lee-Carter model. `"plat"` is the reduced Plat model for older ages as the
#' StMoMo vignette defines it (Villegas, Kaishev and Millossovich 2018,
#' Section 4): `logit q(x,t) = a(x) + k1(t) + (xbar - x) k2(t) + g(t - x)`,
#' built with `StMoMo::StMoMo()` from the vignette's age functions (`1` and
#' `xbar - x`, `xbar` the mean age of the data) and its constraint function.
#' The constraints centre both period indexes and remove the level, the
#' linear and the quadratic trend from the cohort effect (regression of the
#' cohort effect on the year of birth and its square, as in Haberman and
#' Renshaw 2011). The period indexes are forecast with a multivariate random
#' walk with drift and the cohort index with ARIMA(2,0,0) with a non-zero
#' mean (the vignette's Table 3). The vignette uses the logit link with
#' initial exposures for all its models; Plat (2009) uses the log link with
#' Poisson deaths, and for all ages adds a third period term `(xbar - x)+
#' k3(t)`, which he suggests dropping when only older ages are of interest.
#' The full model with that term is not offered. In the vignette's 50-year
#' forecasts for England and Wales men aged 55-89, the rates of this model
#' increase, which the vignette calls implausible: the central trend is
#' linked to the estimated cohort effect.
#'
#' @param model `"apc"`, `"cbd"`, `"rh"`, `"m6"`, `"m7"`, `"m8"`, `"plat"` or
#'   `"lc"`.
#' @param clip Cohorts with fewer observations get zero weight.
#' @param gc_order ARIMA order for the cohort index. Default: `c(2, 0, 0)`
#'   for M7 and the Plat model, `c(1, 1, 0)` for the other cohort models.
#' @param xc Age at which the cohort effect of M8 vanishes (`model = "m8"`
#'   only). Cairns et al. (2009) estimate it; StMoMo's `m8()` needs a fixed
#'   value, and its example uses the highest age of the data (89 for ages
#'   55 to 89). Default (`NULL`): the highest age of the data of each fit.
#'   Forecasts can depend strongly on it; compare values in a backtest
#'   ([mt_compare()]) or choose one with [tuned()].
#' @return A model specification.
#' @references
#' Cairns AJG, Blake D, Dowd K, Coughlan GD, Epstein D, Ong A, Balevich I
#' (2009). A quantitative comparison of stochastic mortality models using
#' data from England and Wales and the United States. *North American
#' Actuarial Journal* 13(1): 1-35.
#'
#' Haberman S, Renshaw A (2011). A comparative study of parametric mortality
#' projection models. *Insurance: Mathematics and Economics* 48(1): 35-55.
#'
#' Plat R (2009). On stochastic mortality modeling. *Insurance: Mathematics
#' and Economics* 45(3): 393-404. \doi{10.1016/j.insmatheco.2009.08.006}
#'
#' Villegas AM, Kaishev VK, Millossovich P (2018). StMoMo: an R package for
#' stochastic mortality modeling. *Journal of Statistical Software* 84(3):
#' 1-38. \doi{10.18637/jss.v084.i03}
#' @examples
#' if (requireNamespace("StMoMo", quietly = TRUE)) {
#'   x <- mt_example(deaths = "single", years = 1995:2019)
#'   fc <- mt_forecast(mt_fit(stmomo("cbd"), x), h = 5, n_draws = 50, seed = 1)
#'   old <- mt_regroup(mt_set_grid(x, age = single_ages(55, 89, open = FALSE)),
#'                     age = single_ages(55, 89, open = FALSE))
#'   fc7 <- mt_forecast(mt_fit(stmomo("m7"), old), h = 5, n_draws = 50, seed = 1)
#'   fcp <- mt_forecast(mt_fit(stmomo("plat"), old), h = 5, n_draws = 50, seed = 1)
#' }
#' @export
stmomo <- function(model = c("apc", "cbd", "rh", "m6", "m7", "m8", "plat", "lc"), clip = 3, gc_order = NULL,
                   xc = NULL) {
  model <- match.arg(model)
  need_pkg("StMoMo", "for stmomo()")
  if (!is.null(xc)) {
    if (model != "m8") mt_abort("`xc` is used by model \"m8\" only.", class = "argument")
    check_number(xc, "xc")
  }
  if (!is.null(gc_order) && (!is.numeric(gc_order) || length(gc_order) != 3 || anyNA(gc_order) || any(gc_order < 0))) {
    mt_abort("`gc_order` must be an ARIMA order c(p, d, q).", class = "argument")
  }
  gc_order <- gc_order %||% if (model %in% c("m7", "plat")) c(2, 0, 0) else c(1, 1, 0)
  params <- c(list(model = model, clip = clip, gc_order = gc_order), if (model == "m8") list(xc = xc))
  new_mt_spec("stmomo", params, min_periods = 8,
              label = sprintf("StMoMo %s%s", toupper(model), if (is.null(xc)) "" else sprintf(" (xc = %g)", xc)))
}

#' @export
fit_model.mt_spec_stmomo <- function(spec, D, E, ages, periods, ...) {
  p <- spec$params
  x <- ages$lower
  if (p$model %in% c("apc", "rh", "m6", "m7", "m8", "plat") && any(utils::head(ages$upper - ages$lower, -1) != 1)) {
    mt_abort("StMoMo cohort models need single years of age.", class = "model_data")
  }
  logit <- p$model %in% c("cbd", "m6", "m7", "m8", "plat")
  mod <- switch(p$model,
                apc = StMoMo::apc(link = "log"),
                cbd = StMoMo::cbd(link = "logit"),
                rh = StMoMo::rh(link = "log", cohortAgeFun = "1"),
                m6 = StMoMo::m6(link = "logit"),
                m7 = StMoMo::m7(link = "logit"),
                m8 = StMoMo::m8(link = "logit", xc = p$xc %||% max(x)),
                plat = stmomo_plat(),
                lc = StMoMo::lc(link = "log"))
  Dx <- round(D)
  Ex <- if (logit) E + D / 2 else E
  dimnames(Dx) <- dimnames(Ex) <- list(x, periods)
  wxt <- StMoMo::genWeightMat(ages = x, years = periods, clip = p$clip)
  utils::capture.output(fit <- suppressWarnings(StMoMo::fit(mod, Dxt = Dx, Ext = Ex, ages = x, years = periods,
                                                            wxt = wxt, verbose = FALSE)))
  if (!isTRUE(fit$conv)) mt_abort("StMoMo fit did not converge.", class = "model_fit")
  fv <- stmomo_fitted(fit)
  list(fit = fit, logit = logit, gc_order = p$gc_order, fitted = if (logit) -log(1 - pmin(fv, 0.999)) else fv)
}

# The reduced Plat model of the StMoMo vignette. stmomo_plat() and
# stmomo_plat_constraints() are taken from the code of the vignette of StMoMo
# 0.4.1 (Villegas, Kaishev and Millossovich 2018, Section 4, functions `f2`
# and `constPlat`), GPL (>= 2), by Andres M. Villegas, Pietro Millossovich and
# Vladimir Kaishev. Changes: names and layout only.
stmomo_plat <- function() {
  StMoMo::StMoMo(link = "logit", staticAgeFun = TRUE, periodAgeFun = c("1", stmomo_plat_f2), cohortAgeFun = "1",
                 constFun = stmomo_plat_constraints)
}

stmomo_plat_f2 <- function(x, ages) mean(ages) - x

# Constraints: the cohort effect without level, linear and quadratic trend
# (transformation 13 of the vignette with constants from a regression on the
# cohort and its square), then both period indexes centred (transformation 14)
stmomo_plat_constraints <- function(ax, bx, kt, b0x, gc, wxt, ages) {
  n_years <- dim(wxt)[2]
  x <- ages
  t <- 1:n_years
  c <- (1 - utils::tail(ages, 1)):(n_years - ages[1])
  xbar <- mean(x)
  phi <- stats::coef(stats::lm(gc ~ 1 + c + I(c^2), na.action = stats::na.omit))
  gc <- gc - phi[1] - phi[2] * c - phi[3] * c^2
  kt[2, ] <- kt[2, ] + 2 * phi[3] * t
  kt[1, ] <- kt[1, ] + phi[2] * t + phi[3] * (t^2 - 2 * xbar * t)
  ax <- ax + phi[1] - phi[2] * x + phi[3] * x^2
  ci <- rowMeans(kt, na.rm = TRUE)
  ax <- ax + ci[1] + ci[2] * (xbar - x)
  kt[1, ] <- kt[1, ] - ci[1]
  kt[2, ] <- kt[2, ] - ci[2]
  list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
}

# Fitted rates; cohorts with zero weight (`clip`) have no estimated effect,
# so their cells take the effect of the nearest estimated cohort
stmomo_fitted <- function(fit) {
  gc <- fit$gc
  if (is.null(gc) || !anyNA(gc) || all(is.na(gc))) return(stats::fitted(fit, type = "rates"))
  ok <- which(!is.na(gc))
  gc[] <- gc[ok[vapply(seq_along(gc), function(i) which.min(abs(ok - i)), integer(1))]]
  stats::predict(fit, years = fit$years, kt = fit$kt, gc = gc, type = "rates")
}

#' @export
simulate_model.mt_spec_stmomo <- function(spec, fit, h, n_draws, ...) {
  # StMoMo's simulate() needs h >= 2; the first h years are kept
  utils::capture.output(sim <- suppressWarnings(stats::simulate(fit$fit, nsim = n_draws, h = max(h, 2),
                                                                gc.order = fit$gc_order)))
  r <- sim$rates[, seq_len(h), , drop = FALSE]
  if (fit$logit) r <- -log(1 - pmin(r, 0.999))
  array(as.numeric(r), c(dim(r)[1], h, n_draws))
}
