# Product-ratio functional demographic model.
# Adapted from demography 2.0.1 (coherentfdm() and forecast.fdmpr()),
# GPL (>= 3), by Rob J Hyndman, with contributions from Heather Booth,
# Leonie Tickle and John Maindonald. Changes: the product is the geometric
# mean of the bottom strata of each group (aggregates of a hierarchy get a
# ratio but do not enter the product), the rates of each stratum are
# smoothed over age with smooth.demogdata() before the decomposition, and
# forecasts are simulated paths (product and ratios simulated independently,
# rescaled path by path so that the ratios of each group multiply to one).
# The principal components are weighted (weight = TRUE, beta = 0.05 of
# ftsa's ftsm()), as in Hyndman, Booth and Yasmeen (2013); demography does
# not weight them by default.

#' Coherent functional demographic model (product-ratio, demography)
#'
#' The product-ratio method of Hyndman, Booth and Yasmeen (2013), with the
#' functional data models of the demography package: for a group of strata
#' (for example women and men, or the deprivation quintiles of one sex)
#'
#' * the **product**, the geometric mean of the strata's rates, is modelled
#'   with a functional demographic model ([hyndman_ullah()]) whose scores are
#'   forecast with ARIMA models, so it carries the common trend;
#' * the **ratio** of each stratum to the product is modelled with a
#'   functional demographic model whose scores are forecast with stationary
#'   ARFIMA models, so the differences between the strata stay bounded and
#'   revert towards their recent level.
#'
#' Forecasts of the strata are products of the two. This is the method of
#' `demography::coherentfdm()` and `forecast()` for its result; [li_lee()]
#' is a related model with a common Lee-Carter factor.
#'
#' Draws simulate the product and each ratio (uncertainty of the scores and
#' resampled model errors, as `simulate()` for functional forecasts in
#' demography) and multiply them; the ratios of each path are rescaled so
#' that they multiply to one, as `forecast()` in demography does for the
#' point forecasts.
#'
#' **Relation to the paper.** As in Hyndman, Booth and Yasmeen (2013,
#' Sections 2.1-2.4): the rates of each stratum are smoothed over age first
#' (weighted penalised regression splines, monotone above age 65, through
#' `demography::smooth.demogdata()`); the product is the geometric mean of
#' the smoothed rates of the group (their Eq. 3), and every stratum has a
#' ratio model (their Eq. 4; they note that J - 1 ratio models would
#' suffice); the scores of the product are forecast with ARIMA models
#' chosen by the algorithm of Hyndman and Khandakar (2008), those of the
#' ratios with stationary ARFIMA models (orders by the same algorithm,
#' maximum likelihood); there is no jump-off adjustment; `order1 = order2 =
#' 6` are the numbers of components of their applications. The principal
#' components are weighted, with weight `beta (1 - beta)^i` for the year `i`
#' years before the last and `beta = 0.05` as in their applications
#' (Section 2.2 and Discussion); `demography::coherentfdm()` does not
#' weight them by default (`weighted = FALSE`). The departures:
#'
#' * the fractional differencing parameter of the ratio models lies in
#'   `drange = c(0, 0.5)` (the default of demography); the paper names
#'   `-0.5 < d < 0.5` as the stationary range;
#' * prediction intervals come from simulated paths with resampled model
#'   errors; the paper gives analytic variances for the rates (their Eq. 8)
#'   and simulates for life expectancy. Predictive draws add Poisson deaths
#'   in place of the observational variance of their Eq. 8;
#' * the ratios of each simulated path are rescaled to multiply to one, and
#'   with a hierarchy (see below) the aggregates get a ratio but do not
#'   enter the product. Both are package-specific (reference unknown);
#'   demography rescales the point forecasts in the same way, and the paper
#'   states the constraint for the observed ratios.
#' Predictive draws (`mt_forecast(predictive = TRUE)`) add Poisson variation
#' of the deaths given the rates and the exposures, and no other observation
#' noise. If intervals are too narrow on real data, that is a property of the
#' model (see [recalibrated()], [ensemble()], or
#' overdispersed models such as `stanmomo(family = "nb")`).
#'
#' With `common_by` (for example `"sex"`), each level of that stratum
#' variable forms its own group with its own product (for example the
#' quintiles of each sex). If the data carry a hierarchy (see
#' [mt_aggregate_strata()]), the product is computed from the bottom strata;
#' aggregates get their own ratio to it. Every group needs at least two
#' strata ([hyndman_ullah()] models a single stratum).
#'
#' @param order1,order2 Number of principal components of the product and
#'   of the ratios.
#' @param common_by Stratum variable(s) defining groups with their own
#'   product.
#' @param smooth Smooth the rates of each stratum and year over age before
#'   the decomposition (see [hyndman_ullah()]).
#' @param smooth_args List of arguments of `demography::smooth.demogdata()`.
#' @param K Number of most recent years used to forecast the ratios.
#' @param drange Range of the fractional differencing parameter of the
#'   ARFIMA models of the ratios (`c(0, 0.5)`: stationary).
#' @param weighted,beta Weighted principal components of the product and of
#'   the ratios, with geometrically decaying weights `beta (1 - beta)^i` of
#'   the years (Hyndman, Booth and Yasmeen 2013: `beta = 0.05`).
#'   `weighted = FALSE` gives the unweighted decomposition, the default of
#'   `demography::coherentfdm()`.
#' @return A model specification.
#' @references
#' Hyndman RJ, Booth H, Yasmeen F (2013). Coherent mortality forecasting: the
#' product-ratio method with functional time series models. *Demography*
#' 50(1): 261-283. \doi{10.1007/s13524-012-0145-5}
#'
#' Hyndman RJ, Khandakar Y (2008). Automatic time series forecasting: the
#' forecast package for R. *Journal of Statistical Software* 27(3): 1-22.
#' \doi{10.18637/jss.v027.i03}
#' @examples
#' if (requireNamespace("demography", quietly = TRUE)) {
#'   x <- mt_example(deaths = "single", years = 1990:2019)
#'   fit <- mt_fit(coherent_fdm(order1 = 2, order2 = 1), x)
#'   fc <- mt_forecast(fit, h = 10, n_draws = 50, seed = 1)
#' }
#' @export
coherent_fdm <- function(order1 = 6, order2 = 6, common_by = NULL, smooth = TRUE, smooth_args = list(), K = 100,
                         drange = c(0, 0.5), weighted = TRUE, beta = 0.05) {
  need_pkg("demography", "for coherent_fdm()")
  check_flag(smooth, "smooth")
  check_flag(weighted, "weighted")
  check_number(beta, "beta", lower = 1e-6, upper = 1)
  if (!is.list(smooth_args)) mt_abort("`smooth_args` must be a list.", class = "argument")
  label <- if (is.null(common_by)) "Coherent fdm" else
    sprintf("Coherent fdm (product by %s)", paste(common_by, collapse = ", "))
  new_mt_spec("coherent_fdm", list(order1 = order1, order2 = order2, common_by = common_by, smooth = smooth,
                                   smooth_args = smooth_args, K = K, drange = drange, weighted = weighted,
                                   beta = beta),
              joint = TRUE, min_periods = 8, label = label)
}

#' @export
fit_model.mt_spec_coherent_fdm <- function(spec, D, E, ages, periods, strata = NULL, hierarchy = NULL, ...) {
  p <- spec$params
  ns <- dim(D)[3]
  group <- rep(1L, ns)
  if (!is.null(p$common_by)) {
    if (is.null(strata) || !all(p$common_by %in% names(strata))) {
      mt_abort("`common_by` must name stratum variables.", class = "argument")
    }
    key <- do.call(paste, c(lapply(p$common_by, function(v) strata[[v]]), sep = "."))
    group <- match(key, unique(key))
  }
  pool <- rep(TRUE, ns)
  if (!is.null(hierarchy) && !is.null(strata)) pool <- hierarchy_structure(strata, hierarchy)$is_bottom
  pool_size <- vapply(seq_len(max(group)), function(g) {
    if (any(group == g & pool)) sum(group == g & pool) else sum(group == g)
  }, 1L)
  if (any(pool_size < 2)) {
    mt_abort(paste("coherent_fdm() needs at least two strata in every group: the product is their geometric mean,",
                   "and the ratio of a single stratum to it is always one. Use hyndman_ullah() for a single",
                   "stratum."), class = "model_data")
  }
  rates <- lapply(seq_len(ns), function(s) {
    dd <- fdm_demogdata(pmax(D[, , s], 0.5) / E[, , s], E[, , s], ages, periods, "s")
    if (isTRUE(p$smooth)) dd <- fdm_smooth(dd, ages$lower, p$smooth_args)
    as.matrix(dd$rate[[1]])
  })
  order1 <- min(p$order1, dim(D)[1] - 1)
  order2 <- min(p$order2, dim(D)[1] - 1)
  # weighted principal components (Hyndman, Booth and Yasmeen 2013); specs
  # saved before the argument existed were unweighted
  weighted <- isTRUE(p$weighted)
  groups <- lapply(seq_len(max(group)), function(g) {
    use <- which(group == g & pool)
    if (!length(use)) use <- which(group == g)
    prod <- exp(Reduce(`+`, lapply(rates[use], log)) / length(use))
    Eg <- Reduce(`+`, lapply(use, function(s) E[, , s]))
    fit_p <- suppressWarnings(demography::fdm(fdm_demogdata(prod, Eg, ages, periods, "product"),
                                              series = "product", order = order1, weight = weighted,
                                              beta = p$beta))
    list(fit = fit_p, pool = use)
  })
  ratios <- lapply(seq_len(ns), function(s) {
    g <- groups[[group[s]]]
    prod_obs <- exp(Reduce(`+`, lapply(rates[g$pool], log)) / length(g$pool))
    r <- rates[[s]] / prod_obs
    r[r > 1e9] <- NA
    suppressWarnings(demography::fdm(fdm_demogdata(r, E[, , s], ages, periods, "ratio"), series = "ratio",
                                     order = order2, weight = weighted, beta = p$beta))
  })
  fitted <- array(NA_real_, dim(D))
  for (s in seq_len(ns)) {
    fitted[, , s] <- exp(as.matrix(groups[[group[s]]]$fit$fitted$y) + as.matrix(ratios[[s]]$fitted$y))
  }
  list(groups = groups, group = group, ratios = ratios, K = p$K, drange = p$drange, fitted = fitted)
}

# Evaluate quietly: the model searches of the forecast package print the
# errors of rejected candidate models through try(), to the connection in
# option try.outFile (knitr sets it to the document), so it is set to
# stderr and captured with the messages
quiet_try <- function(expr) {
  out <- NULL
  op <- options(try.outFile = stderr())
  on.exit(options(op), add = TRUE)
  utils::capture.output(out <- suppressWarnings(expr), type = "message")
  out
}

# A demogdata object of one series
fdm_demogdata <- function(rates, pop, ages, periods, name) {
  dimnames(rates) <- dimnames(pop) <- list(ages$lower, periods)
  demography::demogdata(data = rates, pop = pop, ages = ages$lower, years = periods, type = "mortality",
                        label = "mortensor", name = name)
}

#' @export
simulate_model.mt_spec_coherent_fdm <- function(spec, fit, h, n_draws, ...) {
  ns <- length(fit$ratios)
  sim_p <- lapply(fit$groups, function(g) {
    fc <- quiet_try(forecast::forecast(g$fit, h = h, method = "arima"))
    stats::simulate(fc, nsim = n_draws)
  })
  sim_r <- lapply(fit$ratios, function(r) {
    ny <- length(r$year)
    r$weights <- 0 * r$weights + 1
    if (fit$K < ny) r$weights[seq_len(ny - fit$K)] <- 0
    fc <- quiet_try(forecast::forecast(r, h = h, method = "arfima", estim = "mle", drange = fit$drange))
    stats::simulate(fracdiff_sim_models(fc), nsim = n_draws)
  })
  na <- dim(sim_p[[1]])[1]
  out <- array(NA_real_, c(na, h, ns, n_draws))
  for (g in seq_along(fit$groups)) {
    pool <- fit$groups[[g]]$pool
    # rescale the ratios of each path so that those of the group multiply to one
    adj <- exp(Reduce(`+`, lapply(sim_r[pool], function(a) log(array(a, c(na, h, n_draws))))) / length(pool))
    for (s in which(fit$group == g)) {
      out[, , s, ] <- array(sim_p[[g]], c(na, h, n_draws)) * array(sim_r[[s]], c(na, h, n_draws)) / adj
    }
  }
  out
}

# Simulation of the ARFIMA models of the ratio coefficients. Adapted from
# forecast 8.21.1 (simulate.fracdiff() and unfracdiff()), GPL-3, by Rob J
# Hyndman and co-authors. Change: demography's simulate() draws each path
# coefficient by coefficient, and simulate.fracdiff() rebuilds the same
# fixed-coefficient ARMA model of the differenced series at every call;
# here it is built once per coefficient. The draws are unchanged.
fracdiff_sim_models <- function(fc) {
  for (j in seq_along(fc$coeff)) {
    cf <- fc$coeff[[j]]
    nested <- inherits(cf$model, "forecast")
    mod <- if (nested) cf$model$model else cf$model
    if (!inherits(mod, "fracdiff")) next
    x <- forecast::getResponse(mod)
    if (is.null(x) || anyNA(x)) next
    meanx <- mean(x)
    xx <- x - meanx
    y <- fracdiff::diffseries(xx, d = mod$d)
    if (!is.null(stats::tsp(x))) stats::tsp(y) <- stats::tsp(x)
    arma <- forecast::Arima(y, order = c(length(mod$ar), 0, length(mod$ma)), include.mean = FALSE,
                            fixed = c(mod$ar, -mod$ma))
    sim <- structure(list(arma = arma, xx = xx, meanx = meanx, d = mod$d), class = "mt_fracdiff_sim")
    if (nested) fc$coeff[[j]]$model$model <- sim else fc$coeff[[j]]$model <- sim
  }
  fc
}

#' @exportS3Method stats::simulate
#' @noRd
simulate.mt_fracdiff_sim <- function(object, nsim, seed = NULL, future = TRUE, bootstrap = FALSE, innov = NULL,
                                     ...) {
  ysim <- stats::simulate(object$arma, nsim, seed, future = future, bootstrap = bootstrap, innov = innov)
  unfracdiff(object$xx, ysim, length(object$xx), nsim, object$d) + object$meanx
}

# Inverse of fractional differencing of the simulated values y (h steps)
# given the past values x (forecast 8.21.1, as above)
unfracdiff <- function(x, y, n, h, d) {
  bin_c <- (-1)^(0:(n + h)) * choose(d, 0:(n + h))
  xnew <- lhs <- numeric(h)
  rhs <- cumsum(y)
  bs <- cumsum(bin_c[1:h])
  b <- bin_c[(1:n) + 1]
  xnew[1] <- rhs[1] <- y[1] - sum(b * rev(x))
  if (h > 1) {
    for (k in 2:h) {
      b <- b + bin_c[(1:n) + k]
      rhs[k] <- rhs[k] - sum(b * rev(x))
      lhs[k] <- sum(rev(xnew[1:(k - 1)]) * bs[2:k])
      xnew[k] <- rhs[k] - lhs[k]
    }
  }
  tspx <- stats::tsp(x)
  if (is.null(tspx)) tspx <- c(1, length(x), 1)
  stats::ts(xnew, frequency = tspx[3], start = tspx[2] + 1 / tspx[3])
}
