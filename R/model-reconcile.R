#' Reconcile forecast draws over a hierarchy
#'
#' Makes independently forecast strata add up, draw by draw. The hierarchy is
#' one stratum variable with a total level and its parts (for example
#' `cause = all ~ cancer + cvd + other`), within each combination of the
#' other stratum variables. Causes of death share the exposures of their
#' population, so their rates add up and reconciliation works on rates; for
#' parts with their own populations (for example regions) supply `exposure`
#' and reconciliation works on deaths.
#'
#' For each cell (age, year, draw) the base forecasts `y` (total and parts)
#' are projected onto the coherent subspace, `S (S' W^-1 S)^-1 S' W^-1 y`,
#' with summing matrix `S`:
#'
#' * `"ols"`: `W = I`;
#' * `"wls"`: `W` diagonal with the variance of each series' draws (the
#'   0.1 approach: uncertain series absorb more of the gap);
#' * `"mint"`: `W` the covariance of the series' draws, shrunk towards its
#'   diagonal with intensity `shrink`. This follows the minimum-trace idea of
#'   Wickramasuriya et al. (2019), with the covariance of the forecast draws
#'   in place of in-sample residuals.
#'
#' Reconciled values are bounded below by a small positive number.
#'
#' @param draws An [mt_draws] object containing the total and the parts.
#' @param var Stratum variable of the hierarchy.
#' @param total Level of `var` that is the total.
#' @param parts Levels that are the parts (default: all other levels).
#' @param method `"mint"`, `"wls"` or `"ols"`.
#' @param shrink Shrinkage intensity for `"mint"`.
#' @param exposure Optional exposures `[age, period, stratum]` to reconcile
#'   deaths instead of rates.
#' @return The reconciled `mt_draws`, with the largest relative gap before
#'   reconciliation in `$reconciliation`.
#' @references
#' Wickramasuriya SL, Athanasopoulos G, Hyndman RJ (2019). Optimal forecast
#' reconciliation for hierarchical and grouped time series through trace
#' minimization. *Journal of the American Statistical Association*
#' 114(526): 804-819. \doi{10.1080/01621459.2018.1448825}
#' @examples
#' d <- mt_example(deaths = "single", years = 2000:2019)
#' # a toy hierarchy: total = female + male on the death scale
#' fc <- mt_forecast(mt_fit(rwd(), d), h = 3, n_draws = 50, seed = 1)
#' @export
mt_reconcile <- function(draws, var, total, parts = NULL, method = c("mint", "wls", "ols"), shrink = 0.5,
                         exposure = NULL) {
  check_class(draws, "mt_draws", "draws")
  method <- match.arg(method)
  st <- draws$strata
  if (!var %in% names(st)) mt_abort("`var` is not a stratum variable.", class = "argument")
  parts <- parts %||% setdiff(unique(st[[var]]), total)
  others <- setdiff(names(st), c("stratum", "label", var))
  key <- if (length(others)) do.call(paste, c(lapply(others, function(o) st[[o]]), sep = "\r")) else rep("", nrow(st))
  X <- draws$rates
  if (!is.null(exposure)) X <- X * as.vector(exposure)
  gaps <- numeric(0)
  for (k in unique(key)) {
    ids <- st$stratum[key == k]
    tot <- ids[st[[var]][ids] == total]
    prt <- ids[st[[var]][ids] %in% parts]
    if (length(tot) != 1 || length(prt) < 1) next
    nb <- length(prt)
    S <- rbind(rep(1, nb), diag(nb))
    d <- dim(X)
    Y <- cbind(as.vector(X[, , tot, ]), sapply(prt, function(p) as.vector(X[, , p, ])))
    gaps <- c(gaps, max(abs(Y[, 1] - rowSums(Y[, -1, drop = FALSE])) / pmax(Y[, 1], 1e-12)))
    cell <- rep(seq_len(d[1] * d[2]), d[4])
    R <- matrix(NA_real_, nrow(Y), nb + 1)
    for (c0 in unique(cell)) {
      rows <- which(cell == c0)
      Yc <- Y[rows, , drop = FALSE]
      W <- switch(method,
                  ols = diag(nb + 1),
                  wls = diag(pmax(apply(Yc, 2, stats::var), 1e-20), nb + 1),
                  mint = {
                    V <- stats::cov(Yc)
                    V <- (1 - shrink) * V + shrink * diag(diag(V), nb + 1)
                    V + diag(1e-20, nb + 1)
                  })
      Wi <- tryCatch(solve(W), error = function(e) MASS_ginv(W))
      G <- solve(t(S) %*% Wi %*% S, t(S) %*% Wi)
      R[rows, ] <- Yc %*% t(S %*% G)
    }
    R <- pmax(R, 1e-12)
    X[, , tot, ] <- R[, 1]
    for (j in seq_len(nb)) X[, , prt[j], ] <- R[, j + 1]
  }
  if (!length(gaps)) mt_abort("No total with parts found for the hierarchy.", class = "argument")
  if (!is.null(exposure)) X <- X / as.vector(exposure)
  draws$rates <- X
  draws$reconciliation <- list(method = method, max_relative_gap_before = max(gaps))
  draws
}
