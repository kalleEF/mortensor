#' Li-Lee coherent multi-population model
#'
#' Joint model for several strata (for example women and men, or regions)
#' whose forecasts do not diverge in the long run (Li and Lee 2005):
#' `log m[x,t,i] = a[x,i] + B[x] K[t] + b[x,i] k[t,i]`. The common factor
#' `B K` is a Lee-Carter model of the pooled rates (sum of deaths over sum of
#' exposures, SVD with `K` refitted to the total deaths of each year; a year
#' where the search finds no such `K` keeps the SVD index, as with `adjust =
#' "dt"` in [lee_carter()]); the stratum factors `b k` are the first
#' `ratio_order` components of an SVD of the residuals of each stratum.
#'
#' With `common_by` (for example `"sex"`), `li_lee()` fits one common factor
#' per level of that stratum variable (women and men each with their own
#' trend), and the other strata (for example regions and deprivation
#' quintiles) deviate from their group's common factor. If the data carry a
#' hierarchy (see [mt_aggregate_strata()]), the common factors are estimated
#' from the bottom strata only, so aggregates are not counted twice.
#'
#' The common period index follows a random walk with drift (with drift
#' uncertainty) and the stratum-specific indexes `k[t,i]` follow AR(1)
#' processes, so the differences between strata stay bounded. Two choices
#' of this package differ from Li and Lee (2005), who model `k[t,i]` as an
#' AR(1) process with a constant:
#'
#' * the AR(1) processes have no constant (intercept): the coefficient is
#'   the least-squares estimate through the origin, and the stratum indexes
#'   revert to zero, their mean over the fitting period;
#' * the AR(1) coefficient is bounded to `|phi| <= 0.99`, a package-specific
#'   safeguard against unit roots and explosive deviations (reference
#'   unknown). A warning (class `mt_warning_li_lee`) names the strata where
#'   the bound binds; their deviations then revert very slowly, and the
#'   assumption of coherence may not suit them.
#'
#' Draws include the innovations of all indexes and the drift uncertainty of
#' the common index; `a`, `B`, `b` and the AR(1) coefficients are fixed.
#' Predictive draws (`mt_forecast(predictive = TRUE)`) add Poisson variation
#' of the deaths given the rates and the exposures, and no other observation
#' noise. If intervals are too narrow on real data, that is a property of
#' the model (see [recalibrated()], [ensemble()], or
#' overdispersed models such as `stanmomo(family = "nb")`).
#'
#' @param ratio_order Number of stratum-specific factors. It is capped by
#'   the number of age groups and the number of years minus 2, with a
#'   warning (class `mt_warning_argument`) that names the requested and the
#'   fitted number; [mt_parameters()] reports the fitted number
#'   (`ratio_order`).
#' @param jump_off `"fitted"` (forecasts start from the fitted rates of the
#'   last year) or `"observed"` (forecasts are multiplied by the ratio of
#'   observed to fitted rates in the last year, per age and stratum; Lee and
#'   Miller 2001). With few factors per stratum the fitted rates of the last
#'   year can be far from the observed ones; the observed jump-off then
#'   avoids a step between the last observed and the first forecast year.
#' @param common_by Stratum variable(s) defining groups with their own common
#'   factor.
#' @return A model specification.
#' @references
#' Li N, Lee R (2005). Coherent mortality forecasts for a group of
#' populations: an extension of the Lee-Carter method. *Demography* 42(3):
#' 575-594. \doi{10.1353/dem.2005.0021}
#'
#' Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality.
#' *Journal of the American Statistical Association* 87(419): 659-671.
#'
#' Lee R, Miller T (2001). Evaluating the performance of the Lee-Carter
#' method for forecasting mortality. *Demography* 38(4): 537-549.
#' @seealso [coherent_fdm()] for the product-ratio method.
#' @examples
#' x <- mt_example(deaths = "single")
#' fc <- mt_forecast(mt_fit(li_lee(), x), h = 10, n_draws = 100, seed = 1)
#' @export
li_lee <- function(ratio_order = 1, common_by = NULL, jump_off = c("fitted", "observed")) {
  jump_off <- match.arg(jump_off)
  new_mt_spec("multipop", list(ratio_order = ratio_order, common_by = common_by, jump_off = jump_off), joint = TRUE,
              min_periods = 8,
              label = paste0(if (is.null(common_by)) "Li-Lee" else sprintf("Li-Lee (common by %s)",
                                                                           paste(common_by, collapse = ", ")),
                             if (jump_off == "observed") ", observed jump-off" else ""))
}

#' @export
fit_model.mt_spec_multipop <- function(spec, D, E, ages, periods, strata = NULL, hierarchy = NULL, ...) {
  p <- spec$params
  ns <- dim(D)[3]
  y <- log(pmax(D, 0.5) / E)
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
  commons <- lapply(seq_len(max(group)), function(g) {
    use <- group == g & pool
    if (!any(use)) use <- group == g
    lc_svd(apply(D[, , use, drop = FALSE], 1:2, sum), apply(E[, , use, drop = FALSE], 1:2, sum))
  })
  base_of <- function(cm) outer(cm$b, cm$k)
  # stratum factors: at most one per age group and the number of years minus 2
  r <- min(p$ratio_order, nrow(D), ncol(D) - 2)
  if (r != p$ratio_order) {
    mt_warn(sprintf(paste("`ratio_order = %s` exceeds what the data allow (the number of age groups, %d, and the",
                          "number of years minus 2, %d); %d stratum factor(s) are fitted."), format(p$ratio_order),
                    nrow(D), ncol(D) - 2, r), class = "argument")
  }
  strata_fits <- lapply(seq_len(ns), function(s) {
    res <- y[, , s] - base_of(commons[[group[s]]])
    a <- rowMeans(res)
    sv <- svd(res - a, nu = r, nv = r)
    b <- sv$u[, seq_len(r), drop = FALSE]
    k <- sweep(sv$v[, seq_len(r), drop = FALSE], 2, sv$d[seq_len(r)], "*")
    ar <- lapply(seq_len(r), function(j) ar1_fit(k[, j]))
    list(a = a, b = b, k = k, ar = ar)
  })
  bound <- vapply(strata_fits, function(st) any(vapply(st$ar, `[[`, logical(1), "bounded")), logical(1))
  if (any(bound)) {
    who <- if (!is.null(strata$label)) strata$label[bound] else which(bound)
    mt_warn(sprintf(paste("Li-Lee: the AR(1) coefficient of the stratum index is bounded at |phi| = 0.99 in %d",
                          "stratum/strata (%s): their deviations from the common trend are close to a unit root",
                          "and revert very slowly; coherence may not suit them."), sum(bound),
                    paste(utils::head(who, 5), collapse = ", ")), class = "li_lee")
  }
  fitted <- array(NA_real_, dim(D))
  for (s in seq_len(ns)) {
    st <- strata_fits[[s]]
    fitted[, , s] <- exp(base_of(commons[[group[s]]]) + st$a + st$b %*% t(st$k))
  }
  nt <- dim(D)[2]
  jump <- if (identical(p$jump_off, "observed")) y[, nt, , drop = TRUE] - log(fitted[, nt, , drop = TRUE]) else NULL
  if (!is.null(jump)) jump <- matrix(jump, dim(D)[1], ns)
  list(commons = commons, group = group, strata = strata_fits, fitted = fitted, jump = jump, ratio_order = r)
}

# AR(1) without constant, least squares through the origin; the coefficient
# is bounded to |phi| <= 0.99 (a package safeguard; `bounded` reports it)
ar1_fit <- function(k) {
  n <- length(k)
  phi_ls <- sum(k[-1] * k[-n]) / sum(k[-n]^2)
  phi <- max(min(phi_ls, 0.99), -0.99)
  e <- k[-1] - phi * k[-n]
  list(phi = phi, sigma = sqrt(sum(e^2) / max(n - 2, 1)), last = k[n], bounded = isTRUE(abs(phi_ls) > 0.99))
}

#' @export
simulate_model.mt_spec_multipop <- function(spec, fit, h, n_draws, ...) {
  na <- length(fit$commons[[1]]$a)
  ns <- length(fit$strata)
  kps <- lapply(fit$commons, function(cm) {
    rw_paths(cm$k[length(cm$k)], cm$drift, matrix(cm$dk - cm$drift, nrow = 1), h, n_draws, "gaussian",
             sigma = cm$sigma_k, exp_out = FALSE)
  })
  out <- array(NA_real_, c(na, h, ns, n_draws))
  for (j in seq_len(n_draws)) {
    for (s in seq_len(ns)) {
      g <- fit$group[s]
      cm <- fit$commons[[g]]
      base <- outer(cm$b, kps[[g]][1, , j])
      st <- fit$strata[[s]]
      kk <- vapply(st$ar, function(a) {
        v <- numeric(h)
        prev <- a$last
        for (t in seq_len(h)) {
          prev <- a$phi * prev + stats::rnorm(1, 0, a$sigma)
          v[t] <- prev
        }
        v
      }, numeric(h))
      lp <- base + st$a + st$b %*% t(matrix(kk, h))
      if (!is.null(fit$jump)) lp <- lp + fit$jump[, s]
      out[, , s, j] <- exp(lp)
    }
  }
  out
}

#' Enforce an ordering between strata in forecast draws
#'
#' For each draw, age and year where the rate of the `lower` stratum exceeds
#' that of the `upper` stratum (for example female above male mortality),
#' both rates are replaced so that their mean is kept and their ratio equals
#' the mean observed ratio of the last `years` training years (bounded to
#' 0.02-0.98), applied path by path
#' and reported as the share of adjusted cells.
#'
#' This is a post-hoc adjustment of this package without a published source
#' (reference unknown): it changes the draws and leaves the fitted model as
#' it is, so the adjusted draws no longer follow the model's distribution. Models
#' that keep strata together by construction are [li_lee()],
#' [coherent_fdm()] and `coda(coherent = TRUE)` ([coda()]); [reconciled()]
#' makes totals and their parts coherent.
#'
#' @param draws An [mt_draws] object.
#' @param var Stratum variable that distinguishes the pair (for example
#'   `"sex"`).
#' @param lower,upper Values of `var` that should satisfy lower <= upper.
#' @param years Years used for the observed ratio.
#' @return The adjusted `mt_draws`, with the share of adjusted cells in
#'   `$coherence`.
#' @examples
#' x <- mt_example(deaths = "single")
#' fc <- mt_forecast(mt_fit(rwd(), x), h = 10, n_draws = 50, seed = 1)
#' fc2 <- mt_order_strata(fc, "sex", lower = "female", upper = "male")
#' fc2$coherence
#' @export
mt_order_strata <- function(draws, var, lower, upper, years = 5) {
  check_class(draws, "mt_draws", "draws")
  st <- draws$strata
  if (!var %in% names(st)) mt_abort("`var` is not a stratum variable.", class = "argument")
  others <- setdiff(names(st), c("stratum", "label", var))
  lo <- st[get(var) == lower]
  up <- st[get(var) == upper]
  pairs <- if (length(others)) merge(lo, up, by = others, suffixes = c(".lo", ".up")) else
    data.table(stratum.lo = lo$stratum, stratum.up = up$stratum)
  if (!nrow(pairs)) mt_abort("No pairs of strata found.", class = "argument")
  obs <- draws$observed
  nt <- dim(obs)[2]
  use <- max(1, nt - years + 1):nt
  adj <- numeric(0)
  for (i in seq_len(nrow(pairs))) {
    a <- pairs$stratum.lo[i]
    b <- pairs$stratum.up[i]
    r <- rowMeans(obs[, use, a, drop = FALSE]) / rowMeans(obs[, use, b, drop = FALSE])
    r <- pmin(pmax(r, 0.02), 0.98)
    L <- draws$rates[, , a, , drop = FALSE]
    U <- draws$rates[, , b, , drop = FALSE]
    bad <- L > U
    rr <- array(r, dim(L))
    avg <- (L + U) / 2
    L[bad] <- (avg * 2 * rr / (1 + rr))[bad]
    U[bad] <- (avg * 2 / (1 + rr))[bad]
    draws$rates[, , a, ] <- L
    draws$rates[, , b, ] <- U
    adj <- c(adj, mean(bad))
  }
  draws$coherence <- data.table(lower = st$label[pairs$stratum.lo], upper = st$label[pairs$stratum.up],
                                share_adjusted = adj)
  draws
}
