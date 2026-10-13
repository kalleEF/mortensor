#' Rolling-origin backtest
#'
#' For each forecast origin, the data are cut at the origin with
#' [mt_window()], the preprocessing pipeline and the model are fitted on the
#' training window only, and predictive draws for the next `h` years are
#' compared with the observations. Nothing after the origin enters the fit,
#' including ungrouping and smoothing (they are part of `prep`).
#'
#' Evaluation happens at the resolution of the observed deaths: grid draws of
#' deaths (Poisson given the exposures of the test years, which are treated
#' as known) are summed to each observed cell, so models on single years of
#' age can be scored against data in 5-year groups. Exposures of the test
#' years are ungrouped if needed (no deaths of the test years are used for
#' that).
#'
#' @param spec A model specification.
#' @param data An [mt_data] object.
#' @param origins Last training years. Default: every year from the
#'   `n_origins`-th last possible origin, such that horizon 1 is observed.
#' @param h Forecast horizon.
#' @param prep Preprocessing pipeline ([mt_prep()]).
#' @param n_draws Number of draws per origin.
#' @param n_origins Number of origins when `origins` is `NULL`.
#' @param seed Random seed.
#' @param cores Cores used over origins (forked processes; 1 on Windows).
#' @return An object of class `mt_backtest`: a `data.table` with one row per
#'   observed cell, origin and horizon: `model`, `origin`, `horizon`,
#'   stratum columns, `age_lo`, `age_hi`, `period`, `observed` (log rate),
#'   `median`, `mean`, quantiles, and the scores `crps`, `wis`, `ae`
#'   (absolute error of the median), `cov50`, `cov80`, `cov95`, `pit` and the
#'   MASE scale `scale` (mean absolute year-to-year change in the training
#'   data). Summarise with [mt_score()].
#' @examples
#' x <- mt_example(deaths = "single", years = 1990:2012)
#' bt <- mt_backtest(rwd(), x, origins = c(2005, 2007), h = 5, n_draws = 100)
#' mt_score(bt)
#' @export
mt_backtest <- function(spec, data, origins = NULL, h = 10, prep = NULL, n_draws = 200, n_origins = 5,
                        seed = 1, cores = 1) {
  check_class(data, "mt_data", "data")
  last <- max(mt_grid(data)$period$lower)
  origins <- origins %||% seq(last - h - n_origins + 2, last - 1)
  origins <- origins[origins < last]
  if (!length(origins)) mt_abort("No origin leaves test data.", class = "backtest")
  one <- function(o) {
    res <- tryCatch(backtest_origin(spec, data, o, h, prep, n_draws, seed + o), error = function(e) e)
    if (inherits(res, "error")) {
      mt_warn(sprintf("Origin %s failed: %s", o, conditionMessage(res)), class = "backtest_origin")
      return(NULL)
    }
    res
  }
  out <- rbindlist(par_lapply(origins, one, cores), fill = TRUE)
  if (!nrow(out)) mt_abort("All origins failed.", class = "backtest")
  setattr(out, "class", c("mt_backtest", class(out)))
  out
}

backtest_origin <- function(spec, data, origin, h, prep, n_draws, seed) {
  train <- mt_window(data, to = origin)
  fit <- mt_fit(spec, train, prep = prep)
  last_fit <- max(fit$periods)
  to <- min(origin + h, max(mt_grid(data)$period$lower))
  test <- mt_window(data, from = last_fit + 1, to = to)
  test <- mt_set_grid(test, age = fit$age)
  E_grid <- test_exposures(test, fit)
  fc <- mt_forecast(fit, to = to, n_draws = n_draws, seed = seed, predictive = TRUE, exposure = E_grid)
  D_draws <- fc$rates * as.vector(E_grid)
  obs <- test$obs$deaths
  obs <- obs[period_lo > origin]
  gt <- mt_grid(test)
  idx <- obs_grid_index(obs, gt, "deaths", strict = TRUE)
  ok <- !is.na(idx$a1)
  obs <- obs[ok]
  idx <- idx[ok]
  na <- length(gt$age)
  np <- length(gt$period)
  probs <- wis_probs()
  rows <- list()
  for (s in unique(obs$stratum)) {
    if (anyNA(D_draws[, , s, ])) next
    sel <- obs$stratum == s
    C <- build_operator(idx[sel], na, np)
    Ds <- matrix(D_draws[, , s, ], ncol = n_draws)
    Ecell <- as.numeric(C %*% as.vector(E_grid[, , s]))
    Dcell <- as.matrix(C %*% Ds)
    X <- log(pmax(Dcell, 0.5) / Ecell)
    o <- obs[sel]
    y <- log(pmax(o$value, 0.5) / Ecell)
    sc <- score_samples(X, y, probs)
    scale <- mase_scale(fit, s, o)
    rows[[length(rows) + 1L]] <- cbind(
      data.table(model = spec$label, origin = origin, horizon = o$period_lo - origin, stratum = s,
                 age_lo = o$age_lo, age_hi = o$age_hi, period = o$period_lo, observed = y, scale = scale),
      sc)
  }
  if (!length(rows)) mt_abort("The model failed in every stratum.", class = "backtest")
  out <- rbindlist(rows)
  out <- merge(data$strata, out, by = "stratum", sort = FALSE)
  out[]
}

# Exposures of the forecast periods on the grid (ungrouped if needed; deaths not used)
test_exposures <- function(test, fit) {
  t2 <- test
  t2$obs$deaths <- NULL
  g <- mt_grid(t2)
  if (needs_ungroup(t2$obs$exposure, g)) t2 <- mt_ungroup(t2)
  E <- mt_array(t2, "exposure")
  if (anyNA(E)) mt_abort("Exposures of the test periods are incomplete.", class = "backtest")
  E
}

# Mean absolute year-to-year change of the observed log rates in the training
# data (after preprocessing), aggregated to each scored age cell
mase_scale <- function(fit, s, o) {
  cells <- unique(data.table(age_lo = o$age_lo, age_hi = o$age_hi))
  a1 <- match_bound(cells$age_lo, fit$age$lower)
  a2 <- match_bound(cells$age_hi, fit$age$upper)
  D <- fit$train_D[, , s]
  E <- fit$train_E[, , s]
  cells[, scale := vapply(seq_len(.N), function(i) {
    if (is.na(a1[i]) || is.na(a2[i])) return(NA_real_)
    r <- a1[i]:a2[i]
    lr <- log(pmax(colSums(D[r, , drop = FALSE]), 0.5) / colSums(E[r, , drop = FALSE]))
    mean(abs(diff(lr)))
  }, 0)]
  cells$scale[match(paste(o$age_lo, o$age_hi), paste(cells$age_lo, cells$age_hi))]
}

wis_probs <- function() {
  alpha <- c(0.02, 0.05, seq(0.1, 0.9, 0.1))
  sort(unique(c(alpha / 2, 0.5, 1 - alpha / 2)))
}

#' Proper scores from forecast samples
#'
#' Sample-based CRPS, weighted interval score (Bracher et al. 2021, with the
#' 23 quantiles of the COVID-19 forecast hubs: 11 central intervals plus the
#' median), interval coverage, PIT and the absolute error of the median.
#'
#' @param X Matrix of samples (one row per target).
#' @param y Observed values.
#' @param probs Quantile levels (must contain 0.5 and symmetric pairs).
#' @return A `data.table` with `median`, `mean`, `crps`, `wis`, `ae`,
#'   `cov50`, `cov80`, `cov95`, `pit`.
#' @references
#' Bracher J, Ray EL, Gneiting T, Reich NG (2021). Evaluating epidemic
#' forecasts in an interval format. *PLoS Computational Biology* 17(2):
#' e1008618. \doi{10.1371/journal.pcbi.1008618}
#'
#' Gneiting T, Raftery AE (2007). Strictly proper scoring rules, prediction,
#' and estimation. *Journal of the American Statistical Association*
#' 102(477): 359-378. \doi{10.1198/016214506000001437}
#' @keywords internal
score_samples <- function(X, y, probs = wis_probs()) {
  n <- ncol(X)
  S <- t(apply(X, 1, sort))
  if (n == 1L) S <- t(S)
  Q <- row_quantiles(X, probs)
  colnames(Q) <- paste0("q", probs)
  w <- (2 * seq_len(n) - n - 1)
  crps <- rowMeans(abs(X - y)) - as.numeric(S %*% w) / n^2
  med <- Q[, which.min(abs(probs - 0.5))]
  wis <- wis_from_quantiles(Q, y, probs)
  inside <- function(lo, hi) as.numeric(y >= Q[, lo] & y <= Q[, hi])
  pq <- function(p) which.min(abs(probs - p))
  less <- rowSums(X < y)
  eq <- rowSums(X == y)
  pit <- (less + stats::runif(length(y)) * eq) / n
  out <- data.table(median = med, mean = rowMeans(X), crps = crps, wis = wis, ae = abs(med - y),
                    cov50 = inside(pq(0.25), pq(0.75)), cov80 = inside(pq(0.1), pq(0.9)),
                    cov95 = inside(pq(0.025), pq(0.975)), pit = pit)
  out[]
}

wis_from_quantiles <- function(Q, y, probs) {
  lower <- probs[probs < 0.5]
  K <- length(lower)
  med <- Q[, which.min(abs(probs - 0.5))]
  total <- 0.5 * abs(y - med)
  for (p in lower) {
    alpha <- 2 * p
    l <- Q[, which.min(abs(probs - p))]
    u <- Q[, which.min(abs(probs - (1 - p)))]
    is <- (u - l) + (2 / alpha) * (l - y) * (y < l) + (2 / alpha) * (y - u) * (y > u)
    total <- total + (alpha / 2) * is
  }
  total / (K + 0.5)
}

#' Summarise backtest scores
#'
#' @param backtest An [mt_backtest()] result (or several combined with
#'   `rbind()`).
#' @param by Grouping columns, for example `c("model", "horizon")`.
#' @return A `data.table` with the mean CRPS, WIS and absolute error, MASE
#'   (mean of absolute error over the scale), empirical coverage of the 50%,
#'   80% and 95% intervals, mean PIT and the number of scored cells.
#' @examples
#' x <- mt_example(deaths = "single", years = 1990:2012)
#' bt <- mt_backtest(rwd(), x, origins = 2007, h = 5, n_draws = 100)
#' mt_score(bt, by = c("model", "horizon"))
#' @export
mt_score <- function(backtest, by = "model") {
  if (!is.data.frame(backtest)) mt_abort("`backtest` must be a backtest table.", class = "argument")
  bt <- as.data.table(backtest)
  bt[, list(crps = mean(crps), wis = mean(wis), mae = mean(ae), mase = mean(ae / scale, na.rm = TRUE),
            cov50 = mean(cov50), cov80 = mean(cov80), cov95 = mean(cov95), pit_mean = mean(pit), n = .N),
     by = by]
}

#' @export
print.mt_backtest <- function(x, ...) {
  if (!all(c("model", "origin", "crps", "cov95") %in% names(x))) {
    # a derived table (for example x[, .N, by = horizon]) keeps the class
    print(unclass_backtest(x), ...)
    return(invisible(x))
  }
  cat(sprintf("<mt_backtest> %s; origins %s; %d scored cells\n", paste(unique(x$model), collapse = ", "),
              paste(unique(x$origin), collapse = ", "), nrow(x)))
  print(mt_score(x))
  invisible(x)
}

unclass_backtest <- function(x) {
  y <- copy(x)
  setattr(y, "class", setdiff(class(y), "mt_backtest"))
  y
}
