#' Forecast draws
#'
#' An `mt_draws` object holds simulated death rates as an array
#' `[age, period, stratum, draw]` on the modelling grid, the observed rates
#' of the training data, and metadata (model, forecast origin, seed,
#' diagnostics). Summaries are computed from the paths:
#' `mt_summarise()` gives the mean and quantiles per cell,
#' [life_expectancy()] gives life expectancy path by path.
#'
#' @name mt_draws
#' @param x,draws An `mt_draws` object.
#' @param probs Quantile levels.
#' @return `mt_summarise()` and `as.data.table()` return a `data.table` with
#'   the stratum columns, `age_lo`, `age`, `period`, `mean` and one column per
#'   quantile (`q0.025`, ...).
#' @examples
#' x <- mt_example(deaths = "single")
#' fc <- mt_forecast(mt_fit(rwd(), x), h = 3, n_draws = 50, seed = 1)
#' head(mt_summarise(fc))
NULL

new_mt_draws <- function(rates, age, periods, strata, model, origin, fit = NULL, seed = NULL) {
  dimnames(rates) <- list(age = age$labels, period = as.character(periods), stratum = strata$label,
                          draw = NULL)
  structure(list(rates = rates, age = age, periods = periods, strata = strata, model = model, origin = origin,
                 seed = seed, observed = fit$observed, observed_periods = fit$periods,
                 exposure_last = fit$exposure_last, diagnostics = fit$diagnostics),
            class = "mt_draws")
}

#' @export
print.mt_draws <- function(x, ...) {
  d <- dim(x$rates)
  cat(sprintf("<mt_draws> %s; origin %s; %d ages x %d periods (%s-%s) x %d strata x %d draws\n", x$model,
              x$origin, d[1], d[2], min(x$periods), max(x$periods), d[3], d[4]))
  invisible(x)
}

#' @rdname mt_draws
#' @export
mt_summarise <- function(draws, probs = c(0.025, 0.1, 0.5, 0.9, 0.975)) {
  check_class(draws, "mt_draws", "draws")
  d <- dim(draws$rates)
  M <- matrix(draws$rates, ncol = d[4])
  ok <- rowSums(is.na(M)) == 0
  Q <- matrix(NA_real_, nrow(M), length(probs))
  if (any(ok)) Q[ok, ] <- row_quantiles(M[ok, , drop = FALSE], probs)
  out <- data.table(stratum = rep(seq_len(d[3]), each = d[1] * d[2]),
                    age_lo = rep(draws$age$lower, d[2] * d[3]), age = rep(draws$age$labels, d[2] * d[3]),
                    period = rep(rep(draws$periods, each = d[1]), d[3]), mean = rowMeans(M))
  qn <- paste0("q", format(probs, trim = TRUE, scientific = FALSE))
  for (k in seq_along(probs)) set(out, j = qn[k], value = Q[, k])
  out <- merge(draws$strata, out, by = "stratum", sort = FALSE)
  setorderv(out, c("stratum", "period", "age_lo"))
  out[]
}

#' @rdname mt_draws
#' @param keep.rownames,... Ignored.
#' @method as.data.table mt_draws
#' @export
as.data.table.mt_draws <- function(x, keep.rownames = FALSE, ...) mt_summarise(x)

# Quantiles of each row (type 7), vectorised by sorting
row_quantiles <- function(M, probs) {
  n <- ncol(M)
  S <- t(apply(M, 1, sort))
  if (n == 1L) S <- t(S)
  h <- (n - 1) * probs + 1
  lo <- floor(h)
  hi <- pmin(ceiling(h), n)
  w <- h - lo
  S[, lo, drop = FALSE] * rep(1 - w, each = nrow(M)) + S[, hi, drop = FALSE] * rep(w, each = nrow(M))
}

life_expectancy_draws <- function(x, at = NULL, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  e <- life_expectancy(x$rates, x$age, at = at, ...)
  d <- dim(e)
  M <- matrix(e, ncol = d[3])
  out <- data.table(stratum = rep(seq_len(d[2]), each = d[1]), period = rep(x$periods, d[2]),
                    mean = rowMeans(M))
  ok <- rowSums(is.na(M)) == 0
  Q <- matrix(NA_real_, nrow(M), length(probs))
  if (any(ok)) Q[ok, ] <- row_quantiles(M[ok, , drop = FALSE], probs)
  qn <- paste0("q", format(probs, trim = TRUE, scientific = FALSE))
  for (k in seq_along(probs)) set(out, j = qn[k], value = Q[, k])
  out <- merge(x$strata, out, by = "stratum", sort = FALSE)
  out[, at := at %||% x$age$lower[1]]
  setorderv(out, c("stratum", "period"))
  out[]
}

#' Expected deaths from draws and future exposures
#'
#' Multiplies rate draws by exposures and optionally adds Poisson noise, so
#' that the result describes observable death counts.
#'
#' @param draws An [mt_draws] object.
#' @param exposure Future exposures: an array `[age, period, stratum]` or a
#'   matrix `[age, stratum]` used for every period. Default: the exposures of
#'   the last training year.
#' @param noise Add Poisson noise.
#' @param seed Random seed.
#' @return An array `[age, period, stratum, draw]` of deaths.
#' @examples
#' x <- mt_example(deaths = "single")
#' fc <- mt_forecast(mt_fit(rwd(), x), h = 2, n_draws = 20, seed = 1)
#' dim(mt_deaths(fc))
#' @export
mt_deaths <- function(draws, exposure = NULL, noise = TRUE, seed = NULL) {
  check_class(draws, "mt_draws", "draws")
  d <- dim(draws$rates)
  E <- exposure %||% draws$exposure_last
  if (length(dim(E)) == 2) E <- array(E, c(d[1], 1, d[3]))
  if (length(dim(E)) != 3 || dim(E)[1] != d[1] || dim(E)[3] != d[3]) {
    mt_abort("`exposure` does not match the draws.", class = "argument")
  }
  if (dim(E)[2] == 1) E <- E[, rep(1, d[2]), , drop = FALSE]
  if (dim(E)[2] != d[2]) mt_abort("`exposure` does not match the forecast periods.", class = "argument")
  mu <- draws$rates * as.vector(E)
  if (!noise) return(mu)
  with_seed(seed, array(stats::rpois(length(mu), mu), dim(mu), dimnames(mu)))
}
