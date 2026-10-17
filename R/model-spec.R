#' Model specifications
#'
#' Every model is created by a constructor (`rwd()`, `lee_carter()`,
#' `tucker3()`, `arima_series()`, `ets_series()`, `apc_bayes()`,
#' `ensemble()`) that returns an `mt_spec` object. Specifications are
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
#'   Windows). Each stratum is fitted with its own random number stream,
#'   started from the current generator, so fits that use random numbers
#'   (such as `lee_carter(bootstrap = TRUE)`) are reproduced by `set.seed()`
#'   with any number of cores.
#' @param ... Passed to methods.
#' @return An object of class `mt_fit`.
#' @details Smooth inside the fit (`prep = mt_prep(step_smooth(...))`), not
#'   before it: data returned by [mt_smooth()] carry no raw deaths, so
#'   residuals and observation noise would refer to the smoothed rates and
#'   the intervals would be too narrow. Such data give a warning.
#' @examples
#' x <- mt_example(deaths = "single")
#' fit <- mt_fit(lee_carter(), x)
#' fit
#' @export
mt_fit <- function(spec, data, prep = NULL, cores = 1, ...) {
  if (!isTRUE(.mt_state$fitting) && inherits(data, "mt_data") && !is.null(attr(data, "smooth"))) {
    mt_warn(paste("`data` were smoothed with mt_smooth(): residuals and observation noise then refer to the",
                  "smoothed rates and the forecast intervals are too narrow. Pass the raw data and smooth inside",
                  "the fit: prep = mt_prep(step_smooth(...))."), class = "smoothed_data")
  }
  old <- .mt_state$fitting
  .mt_state$fitting <- TRUE
  on.exit(.mt_state$fitting <- old)
  UseMethod("mt_fit")
}

# State shared by nested calls (a warning about the data is given once per fit)
.mt_state <- new.env(parent = emptyenv())

#' @export
mt_fit.mt_spec <- function(spec, data, prep = NULL, cores = 1, ...) {
  check_class(data, "mt_data", "data")
  pr <- mt_prep_apply_raw(prep, data)
  d <- pr$data
  md <- model_inputs(d, spec)
  raw_D <- raw_deaths(pr$raw, md)
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
    # matrices [age, period] also with a single age group
    slice <- function(A, s) array(A[, , s], dim(A)[1:2], dimnames(A)[1:2])
    one <- function(s) run(function() fit_model(spec, slice(md$D, s), slice(md$E, s), md$age, md$periods))
    # one random number stream per stratum, so that results do not depend on `cores`
    res <- par_lapply(seq_len(nrow(st)), one, cores, streams = TRUE)
    fits <- lapply(res, `[[`, "fit")
    diagnostics <- data.table(stratum = st$label, status = vapply(res, `[[`, "", "status"),
                              message = vapply(res, `[[`, "", "message"),
                              seconds = vapply(res, `[[`, 0, "seconds"))
  }
  if (all(diagnostics$status == "failed")) {
    mt_abort(sprintf("The model failed in every stratum: %s", diagnostics$message[1]), class = "model_fit")
  }
  report_fit_problems(diagnostics, spec$label)
  # observation noise, observed rates and residuals refer to the raw rates
  # (before smoothing); the model is fitted to the preprocessed rates. A model
  # of observed rates (rwd(), arima_series()) fitted to smoothed rates
  # forecasts smoothed rates: its predictive draws add the noise of the raw
  # rates around the smoothed ones
  D_obs <- raw_D %||% md$D
  structure(list(spec = spec, prep = prep, fits = fits, diagnostics = diagnostics, age = md$age,
                 periods = md$periods, strata = st, omega = d$omega,
                 observed = D_obs / md$E, train_D = md$D, raw_D = raw_D, train_E = md$E,
                 exposure_last = md$E[, length(md$periods), , drop = FALSE], level_totals = data$level_totals,
                 noise_sd = extra_noise_sd(spec, fits, D_obs, md$E, smoothed = if (!is.null(raw_D)) md$D / md$E),
                 smoothed_input = !is.null(raw_D)),
            class = "mt_fit")
}

# Raw deaths [age, period, stratum] on the modelling grid, when the
# preprocessing smooths (NULL otherwise)
raw_deaths <- function(raw, md) {
  if (is.null(raw)) return(NULL)
  D <- tryCatch(mt_array(raw, "deaths"), error = function(e) NULL)
  if (is.null(D) || !identical(dim(D), dim(md$D)) || anyNA(D)) return(NULL)
  dimnames(D) <- dimnames(md$D)
  D
}

# One warning per fit for failed strata and for warnings inside strata
report_fit_problems <- function(diagnostics, label) {
  bad <- diagnostics[status == "failed"]
  if (nrow(bad)) {
    mt_warn(sprintf(paste("%s failed in %d of %d strata (%s: %s); their forecasts are missing. Sparse strata often",
                          "fit better with wider age groups (mt_regroup()), pooling (li_lee(common_by = )) or",
                          "smoothing (step_smooth())."), label, nrow(bad), nrow(diagnostics),
                    paste(utils::head(bad$stratum, 3), collapse = ", "), bad$message[1]), class = "model_fit")
  }
  msg <- diagnostics[status == "ok" & nzchar(message)]
  if (nrow(msg)) {
    mt_warn(sprintf("%s: warnings in %d stratum/strata, for example %s: %s (see the diagnostics of the fit).", label,
                    nrow(msg), msg$stratum[1], msg$message[1]), class = "model_fit")
  }
}

# Extra-Poisson noise on the log scale per age and stratum, from in-sample
# fitted rates; for models of observed rates fitted to smoothed rates
# (`smoothed`, [age, period, stratum]), around the smoothed rates
extra_noise_sd <- function(spec, fits, D, E, smoothed = NULL) {
  out <- matrix(0, dim(D)[1], dim(D)[3])
  smooth <- isTRUE(spec$smooth)
  if (!smooth && is.null(smoothed)) return(out)
  fitted <- if (smooth && spec$joint) fits$fitted else NULL
  na <- dim(D)[1]
  for (s in seq_len(dim(D)[3])) {
    f <- if (!smooth) smoothed[, , s] else if (spec$joint) fitted[, , s] else fits[[s]]$fitted
    if (is.null(f)) next
    # matrices [age, period] also with a single age group
    f <- matrix(f, na)
    Es <- matrix(E[, , s], na)
    r <- log_rates(matrix(D[, , s], na), Es) - log(f)
    r[!is.finite(r)] <- NA
    pois <- rowMeans(1 / pmax(f * Es, 0.5), na.rm = TRUE)
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
  coarse <- function(m) {
    nm <- c(deaths = "deaths", exposure = "exposures")
    other <- setdiff(c("deaths", "exposure"), m)
    mt_abort(sprintf(paste("The %s are in coarser age (or period) groups than the grid, which follows the %s.",
                           "Either model the groups, by aggregating all tables to groups they share with",
                           "mt_regroup(x, age = \"common\"), or estimate single years with mt_ungroup() or",
                           "step_ungroup()."), nm[[m]], nm[[other]]),
             class = "model_data", call = NULL)
  }
  D <- tryCatch(mt_array(d, "deaths"), mt_error_resolution = function(e) coarse("deaths"))
  E <- tryCatch(mt_array(d, "exposure"), mt_error_resolution = function(e) coarse("exposure"))
  if (anyNA(D) || anyNA(E)) mt_abort(describe_missing(D, E, g, d$strata), class = "model_data", call = NULL)
  if (any(E <= 0)) {
    i <- which(E <= 0, arr.ind = TRUE)[1, ]
    mt_abort(sprintf(paste("Exposures must be positive on the grid; %d cell(s) are zero or negative, for example %s,",
                           "age %s, %s. At the highest ages, close the age range with an open interval, for example",
                           "mt_regroup(x, age = 0:100)."), sum(E <= 0), d$strata$label[i[3]], g$age$labels[i[1]],
                     g$period$labels[i[2]]),
             class = "model_data", call = NULL)
  }
  r <- D / E
  high <- which(r > 1 & is.finite(g$age$upper)[slice.index(r, 1)])
  if (length(high)) {
    i <- arrayInd(high[1], dim(r))
    mt_warn(sprintf(paste("Death rates above 1 in %d cell(s), for example %s, age %s, %s (%s deaths, exposure %s);",
                          "check the data with mt_validate()."), length(high), d$strata$label[i[3]],
                    g$age$labels[i[1]], g$period$labels[i[2]], format(D[high[1]]), format(E[high[1]])),
            class = "data_quality")
  }
  if (any(g$period$upper - g$period$lower != 1)) {
    mt_abort("Models need single calendar years as periods.", class = "model_data")
  }
  np <- length(g$period)
  if (!isFALSE(spec$consecutive) && np > 1 && any(abs(g$period$lower[-1] - g$period$upper[-np]) > 1e-9)) {
    mt_abort(sprintf(paste("Models need consecutive years; the data have gaps (periods %s). Describe and decompose",
                           "such data with life_expectancy(), mt_decompose() and the other functions for observed",
                           "data, or keep a consecutive span with mt_subset(periods = )."),
                     paste(utils::head(g$period$labels, 6), collapse = ", ")), class = "model_data")
  }
  if (length(g$period) < spec$min_periods) {
    mt_abort(sprintf("%s needs at least %d periods.", spec$label, spec$min_periods), class = "model_data")
  }
  list(D = D, E = E, age = g$age, periods = g$period$lower)
}

# Which cells are missing on the grid, and the likely reason
describe_missing <- function(D, E, g, strata) {
  msg <- "Deaths and exposures must be complete on the grid."
  for (m in c("deaths", "exposures")) {
    A <- if (m == "deaths") D else E
    if (!anyNA(A)) next
    na <- is.na(A)
    per_period <- apply(na, 2, mean)
    per_stratum <- apply(na, 3, mean)
    full_p <- g$period$labels[per_period == 1]
    full_s <- strata$label[per_stratum == 1]
    part <- sum(na) - sum(na[, per_period == 1, ]) - sum(na[, per_period < 1, per_stratum == 1])
    if (length(full_p)) {
      msg <- c(msg, sprintf(paste("No %s in %s (all strata): restrict the years with mt_subset(periods = ), for",
                                  "example when 31 December populations cover a year without deaths."), m,
                            paste(utils::head(full_p, 5), collapse = ", ")))
    }
    if (length(full_s)) {
      msg <- c(msg, sprintf(paste("No %s for stratum/strata %s: check the coding of the stratum variables in the",
                                  "tables (see the warnings of mt_data()) or remove the strata with mt_subset()."), m,
                            paste(utils::head(full_s, 5), collapse = ", ")))
    }
    if (part > 0) {
      i <- which(na & (per_period < 1)[slice.index(na, 2)] & (per_stratum < 1)[slice.index(na, 3)],
                 arr.ind = TRUE)[1, ]
      msg <- c(msg, sprintf(paste("%d %s cell(s) are missing, for example %s, age %s, %s. Suppressed deaths are",
                                  "estimated with mt_fill_suppressed(); single missing years can be excluded with",
                                  "mt_subset(periods = )."), part, m, strata$label[i[3]], g$age$labels[i[1]],
                            g$period$labels[i[2]]))
    }
  }
  paste(msg, collapse = " ")
}

# lapply(), in forked processes with cores > 1. With `streams = TRUE`, task i
# runs with the i-th of independent L'Ecuyer-CMRG random number streams
# started from one draw of the caller's generator, so results follow
# set.seed() and do not depend on `cores`; the caller's generator continues
# after that one draw.
par_lapply <- function(x, f, cores = 1, streams = FALSE) {
  if (streams) {
    seeds <- rng_streams(length(x))
    g <- f
    f <- function(i) with_stream(seeds[[i]], g(x[[i]]))
    x <- seq_along(x)
  }
  if (cores > 1 && .Platform$OS.type == "unix") {
    parallel::mclapply(x, f, mc.cores = cores, mc.set.seed = TRUE)
  } else {
    lapply(x, f)
  }
}

# n independent L'Ecuyer-CMRG streams (as parallel::clusterSetRNGStream())
# from one draw of the current generator, which is otherwise left as it was
rng_streams <- function(n) {
  seed <- sample.int(.Machine$integer.max, 1L)
  old <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  s <- get(".Random.seed", envir = globalenv())
  out <- vector("list", n)
  for (i in seq_len(n)) {
    out[[i]] <- s
    s <- parallel::nextRNGStream(s)
  }
  out
}

# Evaluate `expr` with the random number state `seed` (from rng_streams());
# the caller's state is restored afterwards
with_stream <- function(seed, expr) {
  env <- globalenv()
  old <- if (exists(".Random.seed", envir = env, inherits = FALSE)) get(".Random.seed", envir = env) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = env) else assign(".Random.seed", old, envir = env))
  assign(".Random.seed", seed, envir = env)
  expr
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
#'   models of smooth rates, and for any model fitted to smoothed rates
#'   (`step_smooth()`), extra-Poisson noise (estimated from in-sample
#'   residuals of the raw rates) and Poisson variation of deaths given
#'   `exposure` are added;
#'   cells can have zero deaths (rate 0). Used by [mt_backtest()], which
#'   sums the deaths to the observed cells before taking logarithms.
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
  if (!is.null(to)) {
    if (!missing(h)) mt_warn("Both `h` and `to` are given; `to` is used.", class = "argument")
    h <- to - last
  }
  check_number(h, "h", lower = 1, integerish = TRUE)
  check_number(n_draws, "n_draws", lower = 1, integerish = TRUE)
  spec <- fit$spec
  na <- length(fit$age)
  ns <- nrow(fit$strata)
  arr <- array(NA_real_, c(na, h, ns, n_draws),
               dimnames = list(age = fit$age$labels, period = as.character(last + seq_len(h)),
                               stratum = fit$strata$label, draw = NULL))
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
  check_draws_plausible(arr, fit$age, fit$strata, spec$label)
  out <- new_mt_draws(arr, fit$age, last + seq_len(h), fit$strata, model = spec$label, origin = last, fit = fit,
                      seed = seed)
  out$predictive <- predictive
  out
}

# Warn when forecast death rates exceed 1 at any age, the open interval
# included (most of the population would die within a year): the model is
# then unstable, usually because of sparse counts or a steep trend
check_draws_plausible <- function(arr, age, strata, label) {
  bad <- !is.na(arr) & arr > 1
  if (!any(bad)) return(invisible())
  share <- mean(apply(bad, 4, any))
  by_s <- apply(bad, 3, any)
  ages <- which(apply(bad, 1, any))
  mt_warn(sprintf(paste("%s: forecast death rates above 1 in %.1f%% of the draws (strata %s; ages %s).",
                        "The model is unstable here, often because of sparse counts or a steep trend extrapolated",
                        "far ahead; consider wider age groups",
                        "(mt_regroup()), pooling (li_lee(common_by = )), smoothing (step_smooth()) or an ensemble."),
                  label, 100 * share, paste(utils::head(strata$label[by_s], 3), collapse = ", "),
                  paste(utils::head(age$labels[ages], 4), collapse = ", ")), class = "implausible")
  invisible()
}

# Observed-rate draws: extra-Poisson noise (models of smooth rates, or fitted to smoothed rates) and
# Poisson deaths given exposures.
# Cells may have zero deaths (rate 0); scores that take logarithms floor the
# deaths after summing them to the observed cells, not per grid cell
predictive_rates <- function(arr, fit, exposure, spec) {
  d <- dim(arr)
  E <- expand_exposure(exposure %||% fit$exposure_last, d, dimnames(arr))
  if (isTRUE(spec$smooth) || isTRUE(fit$smoothed_input)) {
    sd <- array(fit$noise_sd[, rep(seq_len(d[3]), each = d[2])], d[1:3])
    arr <- arr * exp(stats::rnorm(length(arr), 0, as.vector(sd)))
    mu <- arr * as.vector(E)
    ok <- !is.na(mu)
    deaths <- array(NA_real_, d)
    deaths[ok] <- stats::rpois(sum(ok), mu[ok])
    arr <- deaths / pmax(as.vector(E), 1e-12)
  }
  arr
}
