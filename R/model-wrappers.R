#' Functional demographic model of Hyndman and Ullah (demography)
#'
#' Wrapper around `demography::fdm()` (Hyndman and Ullah 2007): log death
#' rates are smoothed over age and decomposed into principal components whose
#' scores are forecast with ARIMA or ETS models. Draws are simulated with
#' `demography`'s `simulate()` method for functional forecasts (uncertainty
#' of the coefficients and model error); predictive draws add the
#' extra-Poisson and Poisson variation as for other models of smooth rates.
#'
#' @param order Number of principal components.
#' @param method Forecast method for the coefficients (`"arima"`, `"ets"`,
#'   `"rwdrift"`, ...), passed to `forecast()`.
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
hyndman_ullah <- function(order = 6, method = "arima", ...) {
  need_pkg("demography", "for hyndman_ullah()")
  new_mt_spec("fdm", list(order = order, method = method, args = list(...)), min_periods = 8,
              label = "Functional demographic model")
}

#' @export
fit_model.mt_spec_fdm <- function(spec, D, E, ages, periods, ...) {
  p <- spec$params
  rates <- pmax(D, 0.5) / E
  dimnames(rates) <- dimnames(E) <- list(ages$lower, periods)
  dd <- demography::demogdata(data = rates, pop = E, ages = ages$lower, years = periods, type = "mortality",
                              label = "mortensor", name = "s")
  fit <- suppressWarnings(do.call(demography::fdm, c(list(dd, series = "s", order = min(p$order, nrow(D) - 1)),
                                                     p$args)))
  fv <- fit$fitted$y
  if (is.list(fv)) fv <- fv[[1]]
  list(fit = fit, method = p$method, fitted = exp(as.matrix(fv)))
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
#' Renshaw-Haberman (`"rh"`), M6 (`"m6"`, logit link) and Lee-Carter
#' (`"lc"`). Models are fitted per stratum with real ages (the lower bounds
#' of single-year intervals; models with a cohort effect need single years of
#' age), deaths and central exposures; for the logit-link models
#' initial exposures are approximated as central exposures plus half the
#' deaths, and forecast death probabilities are converted to central rates
#' (`m = -log(1 - q)`). Cohorts observed in fewer than `clip` cells get zero
#' weight. Draws come from `simulate()` of the fitted model (period and
#' cohort indexes simulated as multivariate random walks with drift and
#' ARIMA, respectively).
#'
#' @param model `"apc"`, `"cbd"`, `"rh"`, `"m6"` or `"lc"`.
#' @param clip Cohorts with fewer observations get zero weight.
#' @param gc_order ARIMA order for the cohort index.
#' @return A model specification.
#' @references
#' Villegas AM, Kaishev VK, Millossovich P (2018). StMoMo: an R package for
#' stochastic mortality modeling. *Journal of Statistical Software* 84(3):
#' 1-38. \doi{10.18637/jss.v084.i03}
#' @examples
#' if (requireNamespace("StMoMo", quietly = TRUE)) {
#'   x <- mt_example(deaths = "single", years = 1995:2019)
#'   fc <- mt_forecast(mt_fit(stmomo("cbd"), x), h = 5, n_draws = 50, seed = 1)
#' }
#' @export
stmomo <- function(model = c("apc", "cbd", "rh", "m6", "lc"), clip = 3, gc_order = c(1, 1, 0)) {
  model <- match.arg(model)
  need_pkg("StMoMo", "for stmomo()")
  new_mt_spec("stmomo", list(model = model, clip = clip, gc_order = gc_order), min_periods = 8,
              label = sprintf("StMoMo %s", toupper(model)))
}

#' @export
fit_model.mt_spec_stmomo <- function(spec, D, E, ages, periods, ...) {
  p <- spec$params
  x <- ages$lower
  if (p$model %in% c("apc", "rh", "m6") && any(utils::head(ages$upper - ages$lower, -1) != 1)) {
    mt_abort("StMoMo cohort models need single years of age.", class = "model_data")
  }
  logit <- p$model %in% c("cbd", "m6")
  mod <- switch(p$model,
                apc = StMoMo::apc(link = "log"),
                cbd = StMoMo::cbd(link = "logit"),
                rh = StMoMo::rh(link = "log", cohortAgeFun = "1"),
                m6 = StMoMo::m6(link = "logit"),
                lc = StMoMo::lc(link = "log"))
  Dx <- round(D)
  Ex <- if (logit) E + D / 2 else E
  dimnames(Dx) <- dimnames(Ex) <- list(x, periods)
  wxt <- StMoMo::genWeightMat(ages = x, years = periods, clip = p$clip)
  utils::capture.output(fit <- suppressWarnings(StMoMo::fit(mod, Dxt = Dx, Ext = Ex, ages = x, years = periods,
                                                            wxt = wxt, verbose = FALSE)))
  if (!isTRUE(fit$conv)) mt_abort("StMoMo fit did not converge.", class = "model_fit")
  fv <- stats::fitted(fit, type = "rates")
  list(fit = fit, logit = logit, gc_order = p$gc_order, fitted = if (logit) -log(1 - pmin(fv, 0.999)) else fv)
}

#' @export
simulate_model.mt_spec_stmomo <- function(spec, fit, h, n_draws, ...) {
  utils::capture.output(sim <- suppressWarnings(stats::simulate(fit$fit, nsim = n_draws, h = h,
                                                                gc.order = fit$gc_order)))
  r <- sim$rates
  if (fit$logit) r <- -log(1 - pmin(r, 0.999))
  array(as.numeric(r), c(dim(r)[1], h, n_draws))
}
