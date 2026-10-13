#' Tune model parameters on training data
#'
#' `mt_tune()` evaluates every combination of parameters in `grid` with an
#' inner rolling-origin backtest on `data` and returns the scores and the best
#' specification. `tuned()` wraps the same procedure into a model
#' specification: each time it is fitted (for example at every origin of an
#' outer [mt_backtest()]), the tuning uses only the training data it
#' receives, so the selection never sees the years being forecast.
#'
#' @param spec_fun A function that returns a model specification from named
#'   arguments (for example `function(drift_window) rwd(drift_window)`).
#' @param grid A data frame of parameter combinations (columns match the
#'   arguments of `spec_fun`).
#' @param data An [mt_data] object.
#' @param prep Preprocessing pipeline.
#' @param inner_origins Number of inner origins (the last years that leave
#'   `inner_h` years of test data).
#' @param inner_h Horizon of the inner backtest.
#' @param metric Column of [mt_score()] to minimise (`"crps"`, `"wis"`,
#'   `"mae"`, `"mase"`).
#' @param n_draws Draws per inner forecast.
#' @param strategy `"full"` evaluates every row of `grid`; `"two_stage"`
#'   first evaluates a coarse grid (every second value of each column,
#'   including the first and the last) and then the rows of `grid` within one
#'   step of the best coarse combination in every column. For large grids,
#'   such as Tucker ranks, this needs far fewer fits.
#' @return `mt_tune()`: a list with `scores` (one row per combination) and
#'   `best` (the chosen specification). `tuned()`: a model specification.
#' @examples
#' x <- mt_example(deaths = "single", years = 1990:2012)
#' tn <- mt_tune(function(drift_window) rwd(drift_window = drift_window),
#'               data.frame(drift_window = c(5, 10, 20)), x, inner_origins = 2, inner_h = 3)
#' tn$scores
#' @export
mt_tune <- function(spec_fun, grid, data, prep = NULL, inner_origins = 3, inner_h = 5, metric = "crps",
                    n_draws = 100, strategy = c("full", "two_stage")) {
  if (!is.function(spec_fun)) mt_abort("`spec_fun` must be a function.", class = "argument")
  strategy <- match.arg(strategy)
  grid <- as.data.frame(grid)
  last <- max(mt_grid(data)$period$lower)
  origins <- seq(last - inner_origins - inner_h + 1, last - inner_h)
  score_rows <- function(rows) {
    rbindlist(lapply(rows, function(i) {
      args <- as.list(grid[i, , drop = FALSE])
      sp <- do.call(spec_fun, args)
      bt <- suppressWarnings(tryCatch(mt_backtest(sp, data, origins = origins, h = inner_h, prep = prep,
                                                  n_draws = n_draws), error = function(e) NULL))
      val <- if (is.null(bt)) NA_real_ else mt_score(bt)[[metric]]
      cbind(data.table(row = i), as.data.table(grid[i, , drop = FALSE]), value = val)
    }))
  }
  if (strategy == "full") {
    scores <- score_rows(seq_len(nrow(grid)))
  } else {
    pos <- lapply(grid, function(v) match(v, sort(unique(v))))
    coarse <- Reduce(`&`, lapply(pos, function(k) {
      n <- max(k)
      k %% 2 == 1 | k == n
    }))
    scores <- score_rows(which(coarse))
    if (all(is.na(scores$value))) mt_abort("Every parameter combination failed.", class = "tune")
    b <- scores$row[which.min(scores$value)]
    near <- Reduce(`&`, lapply(pos, function(k) abs(k - k[b]) <= 1))
    scores <- rbind(scores, score_rows(setdiff(which(near), scores$row)))
  }
  setnames(scores, "value", metric)
  if (all(is.na(scores[[metric]]))) mt_abort("Every parameter combination failed.", class = "tune")
  best_i <- scores$row[which.min(scores[[metric]])]
  scores[, row := NULL]
  list(scores = scores[], best = do.call(spec_fun, as.list(grid[best_i, , drop = FALSE])), origins = origins)
}

#' @rdname mt_tune
#' @param label Display name.
#' @export
tuned <- function(spec_fun, grid, inner_origins = 3, inner_h = 5, metric = "crps", n_draws = 100,
                     label = "Tuned model", strategy = c("full", "two_stage")) {
  strategy <- match.arg(strategy)
  new_mt_spec("tuned", list(spec_fun = spec_fun, grid = as.data.frame(grid), inner_origins = inner_origins,
                            inner_h = inner_h, metric = metric, n_draws = n_draws, strategy = strategy),
              min_periods = inner_origins + inner_h + 5, label = label)
}

#' @export
mt_fit.mt_spec_tuned <- function(spec, data, prep = NULL, cores = 1, ...) {
  p <- spec$params
  tn <- mt_tune(p$spec_fun, p$grid, data, prep = prep, inner_origins = p$inner_origins, inner_h = p$inner_h,
                metric = p$metric, n_draws = p$n_draws, strategy = p$strategy %||% "full")
  fit <- mt_fit(tn$best, data, prep = prep, cores = cores)
  fit$tuning <- tn$scores
  fit$spec$label <- sprintf("%s: %s", spec$label, tn$best$label)
  fit
}
