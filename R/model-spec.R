#' Model specifications
#'
#' Every model is created by a constructor (`rwd()`, `lee_carter()`,
#' `tucker3()`, `arima_series()`, `ets_series()`, `apc_bayes()`,
#' `mt_ensemble()`) that returns an `mt_spec` object. Specifications are
#' fitted with [mt_fit()] and forecast with [mt_forecast()].
#'
#' `new_mt_spec()` is the constructor for developers who add models: a model
#' class `mt_spec_<name>` needs methods for the internal generics
#' `fit_model()` and `simulate_model()`. Per-stratum models receive matrices
#' `[age, period]` of deaths and exposures; joint models (`joint = TRUE`)
#' receive arrays `[age, period, stratum]`. `simulate_model()` returns rates
#' `[age, horizon, draw]` (per stratum) or `[age, horizon, stratum, draw]`
#' (joint).
#'
#' @param model Model name (used in the class).
#' @param params List of parameters.
#' @param joint Whether the model is fitted jointly across strata.
#' @param min_periods Minimum number of periods in the training data.
#' @param label Display name.
#' @param smooth Whether the draws describe smooth underlying rates (`TRUE`,
#'   for example Lee-Carter; such fits return in-sample `fitted` rates, from
#'   which extra-Poisson noise is estimated) or observed rates (`FALSE`, for
#'   example a random walk on observed log rates).
#' @return An object of class `c("mt_spec_<model>", "mt_spec")`.
#' @keywords internal
#' @export
new_mt_spec <- function(model, params = list(), joint = FALSE, min_periods = 5, label = model, smooth = TRUE) {
  structure(list(model = model, params = params, joint = joint, min_periods = min_periods, label = label,
                 smooth = smooth),
            class = c(paste0("mt_spec_", model), "mt_spec"))
}

#' @export
print.mt_spec <- function(x, ...) {
  p <- x$params
  show <- function(v) {
    if (is.null(v)) "NULL" else if (is.list(v)) "<list>" else paste(format(v), collapse = ",")
  }
  desc <- if (length(p)) paste(names(p), vapply(p, show, ""), sep = " = ", collapse = ", ") else ""
  cat(sprintf("<mt_spec> %s%s\n", x$label, if (nzchar(desc)) paste0(": ", desc) else ""))
  invisible(x)
}

fit_model <- function(spec, D, E, ages, periods, ...) UseMethod("fit_model")
simulate_model <- function(spec, fit, h, n_draws, ...) UseMethod("simulate_model")

#' @export
fit_model.default <- function(spec, D, E, ages, periods, ...) {
  mt_abort(sprintf("No fit_model() method for model '%s'.", spec$model), class = "model")
}

#' @export
simulate_model.default <- function(spec, fit, h, n_draws, ...) {
  mt_abort(sprintf("No simulate_model() method for model '%s'.", spec$model), class = "model")
}

#' Fit a model to mortality data
#'
#' Applies the preprocessing pipeline, takes deaths and exposures on the grid
#' and fits the model per stratum (or jointly). A failure in one stratum is
#' recorded in the diagnostics and does not stop the other strata.
#'
#' @param spec A model specification (`mt_spec`).
#' @param data An [mt_data] object.
#' @param prep An [mt_prep] pipeline or `NULL`.
#' @param cores Number of cores for per-stratum fits (forked processes; 1 on
#'   Windows).
#' @param ... Passed to methods.
#' @return An object of class `mt_fit`.
#' @examples
#' x <- mt_example(deaths = "single")
#' fit <- mt_fit(lee_carter(), x)
#' fit
#' @export
mt_fit <- function(spec, data, prep = NULL, cores = 1, ...) UseMethod("mt_fit")

#' @export
mt_fit.mt_spec <- function(spec, data, prep = NULL, cores = 1, ...) {
  check_class(data, "mt_data", "data")
  d <- mt_prep_apply(prep, data)
  md <- model_inputs(d, spec)
  st <- d$strata
  run <- function(expr_fun) {
    t0 <- proc.time()[["elapsed"]]
    w <- new.env(parent = emptyenv())
    w$msg <- character()
    res <- withCallingHandlers(tryCatch(expr_fun(), error = function(e) e),
                               warning = function(cnd) {
                                 w$msg <- c(w$msg, conditionMessage(cnd))
                                 invokeRestart("muffleWarning")
                               })
    secs <- proc.time()[["elapsed"]] - t0
    ok <- !inherits(res, "error")
    list(fit = if (ok) res else NULL, status = if (ok) "ok" else "failed",
         message = if (ok) paste(unique(w$msg), collapse = " | ") else conditionMessage(res), seconds = secs)
  }
  if (spec$joint) {
    r <- run(function() fit_model(spec, md$D, md$E, md$age, md$periods, strata = st, hierarchy = d$hierarchy))
    fits <- r$fit
    diagnostics <- data.table(stratum = st$label, status = r$status, message = r$message, seconds = r$seconds)
  } else {
    one <- function(s) run(function() fit_model(spec, md$D[, , s], md$E[, , s], md$age, md$periods))
    res <- par_lapply(seq_len(nrow(st)), one, cores)
    fits <- lapply(res, `[[`, "fit")
    diagnostics <- data.table(stratum = st$label, status = vapply(res, `[[`, "", "status"),
                              message = vapply(res, `[[`, "", "message"),
                              seconds = vapply(res, `[[`, 0, "seconds"))
  }
  if (all(diagnostics$status == "failed")) {
    mt_abort(sprintf("The model failed in every stratum: %s", diagnostics$message[1]), class = "model_fit")
  }
  structure(list(spec = spec, prep = prep, fits = fits, diagnostics = diagnostics, age = md$age,
                 periods = md$periods, strata = st, omega = d$omega,
                 observed = md$D / md$E, train_D = md$D, train_E = md$E,
                 exposure_last = md$E[, length(md$periods), , drop = FALSE],
                 noise_sd = extra_noise_sd(spec, fits, md$D, md$E)),
            class = "mt_fit")
}

# Extra-Poisson noise on the log scale per age and stratum, from in-sample fitted rates
extra_noise_sd <- function(spec, fits, D, E) {
  out <- matrix(0, dim(D)[1], dim(D)[3])
  if (!isTRUE(spec$smooth)) return(out)
  fitted <- if (spec$joint) fits$fitted else NULL
  for (s in seq_len(dim(D)[3])) {
    f <- if (spec$joint) fitted[, , s] else fits[[s]]$fitted
    if (is.null(f)) next
    r <- log_rates(D[, , s], E[, , s]) - log(f)
    r[!is.finite(r)] <- NA
    pois <- rowMeans(1 / pmax(f * E[, , s], 0.5), na.rm = TRUE)
    n <- rowSums(!is.na(r))
    v <- rowSums((r - rowMeans(r, na.rm = TRUE))^2, na.rm = TRUE) / pmax(n - 1, 1) - pois
    v[n < 3 | !is.finite(v)] <- 0
    out[, s] <- sqrt(pmax(v, 0))
  }
  out
}

model_inputs <- function(d, spec) {
  g <- mt_grid(d)
  for (m in c("deaths", "exposure")) {
    if (is.null(d$obs[[m]])) mt_abort(sprintf("Models need a %s table.", m), class = "model_data")
  }
  D <- tryCatch(mt_array(d, "deaths"), mt_error_resolution = function(e) {
    mt_abort("Deaths are coarser than the grid; add step_ungroup() to the preprocessing.", class = "model_data",
             call = NULL)
  })
  E <- tryCatch(mt_array(d, "exposure"), mt_error_resolution = function(e) {
    mt_abort("Exposures are coarser than the grid; add step_ungroup() to the preprocessing.", class = "model_data",
             call = NULL)
  })
  if (anyNA(D) || anyNA(E)) {
    mt_abort(paste("Deaths and exposures must be complete on the grid (no missing cells); suppressed cells can be",
                   "estimated from published totals with mt_fill_suppressed()."), class = "model_data")
  }
  if (any(E <= 0)) mt_abort("Exposures must be positive on the grid.", class = "model_data")
  if (any(g$period$upper - g$period$lower != 1)) {
    mt_abort("Models need single calendar years as periods.", class = "model_data")
  }
  if (length(g$period) < spec$min_periods) {
    mt_abort(sprintf("%s needs at least %d periods.", spec$label, spec$min_periods), class = "model_data")
  }
  list(D = D, E = E, age = g$age, periods = g$period$lower)
}

par_lapply <- function(x, f, cores = 1) {
  if (cores > 1 && .Platform$OS.type == "unix") {
    parallel::mclapply(x, f, mc.cores = cores, mc.set.seed = TRUE)
  } else {
    lapply(x, f)
  }
}

#' @export
print.mt_fit <- function(x, ...) {
  cat(sprintf("<mt_fit> %s; %d strata; ages %s; periods %s-%s\n", x$spec$label, nrow(x$strata),
              paste(x$age$labels[c(1, length(x$age))], collapse = " to "), min(x$periods), max(x$periods)))
  if (!is.null(x$prep)) cat(sprintf("  preprocessing: %d step(s)\n", length(x$prep$steps)))
  bad <- x$diagnostics[status != "ok"]
  if (nrow(bad)) cat(sprintf("  failed: %s\n", paste(bad$stratum, collapse = ", ")))
  invisible(x)
}

#' Forecast draws from a fitted model
#'
#' Simulates future death rates on the grid. The result holds simulated paths
#' (draws), so that uncertainty propagates correctly to derived quantities
#' such as life expectancy.
#'
#' @param fit An `mt_fit` object.
#' @param h Forecast horizon in years.
#' @param to Alternatively, the last year to forecast.
#' @param n_draws Number of simulated paths.
#' @param seed Random seed (`NULL` leaves the random number generator as
#'   is).
#' @param predictive If `TRUE`, draws describe observable rates: for
#'   models of smooth rates, extra-Poisson noise (estimated from in-sample
#'   residuals) and Poisson variation of deaths given `exposure` are added,
#'   and zero counts are replaced by 0.5 as in the training data. Used by
#'   [mt_backtest()].
#' @param exposure Exposures for predictive draws: array `[age, period,
#'   stratum]` or matrix `[age, stratum]` (default: last training year).
#' @param ... Passed to methods.
#' @return An [mt_draws] object.
#' @examples
#' x <- mt_example(deaths = "single")
#' fc <- mt_forecast(mt_fit(rwd(), x), h = 5, n_draws = 100, seed = 1)
#' fc
#' @export
mt_forecast <- function(fit, h = 10, to = NULL, n_draws = 500, seed = NULL, predictive = FALSE,
                        exposure = NULL, ...) {
  UseMethod("mt_forecast")
}

#' @export
mt_forecast.mt_fit <- function(fit, h = 10, to = NULL, n_draws = 500, seed = NULL, predictive = FALSE,
                               exposure = NULL, ...) {
  last <- max(fit$periods)
  if (!is.null(to)) h <- to - last
  check_number(h, "h", lower = 1, integerish = TRUE)
  check_number(n_draws, "n_draws", lower = 1, integerish = TRUE)
  spec <- fit$spec
  na <- length(fit$age)
  ns <- nrow(fit$strata)
  arr <- array(NA_real_, c(na, h, ns, n_draws))
  with_seed(seed, {
    if (spec$joint) {
      if (!is.null(fit$fits)) arr[] <- simulate_model(spec, fit$fits, h, n_draws)
    } else {
      for (s in seq_len(ns)) {
        if (!is.null(fit$fits[[s]])) arr[, , s, ] <- simulate_model(spec, fit$fits[[s]], h, n_draws)
      }
    }
    if (predictive) arr <- predictive_rates(arr, fit, exposure, spec)
  })
  out <- new_mt_draws(arr, fit$age, last + seq_len(h), fit$strata, model = spec$label, origin = last, fit = fit,
                      seed = seed)
  out$predictive <- predictive
  out
}

# Observed-rate draws: extra-Poisson noise (smooth models) and Poisson deaths given exposures
predictive_rates <- function(arr, fit, exposure, spec) {
  d <- dim(arr)
  E <- exposure %||% fit$exposure_last
  if (length(dim(E)) == 2) {
    if (!all(dim(E) == d[c(1, 3)])) mt_abort("`exposure` must be [age, stratum] or [age, period, stratum].",
                                             class = "argument")
    E <- array(E, c(d[1], 1, d[3]))
  }
  if (dim(E)[2] == 1) E <- E[, rep(1, d[2]), , drop = FALSE]
  if (!all(dim(E) == d[1:3])) mt_abort("`exposure` must be [age, period, stratum] matching the forecast.",
                                       class = "argument")
  if (isTRUE(spec$smooth)) {
    sd <- array(fit$noise_sd[, rep(seq_len(d[3]), each = d[2])], d[1:3])
    arr <- arr * exp(stats::rnorm(length(arr), 0, as.vector(sd)))
    mu <- arr * as.vector(E)
    ok <- !is.na(mu)
    deaths <- array(NA_real_, d)
    deaths[ok] <- stats::rpois(sum(ok), mu[ok])
    arr <- pmax(deaths, 0.5) / as.vector(E)
  }
  arr
}
