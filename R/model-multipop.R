#' Coherent multi-population models
#'
#' Joint models for several strata (for example women and men, or regions)
#' whose forecasts do not diverge in the long run:
#'
#' * `li_lee()` (Li and Lee 2005): `log m[x,t,i] = a[x,i] + B[x] K[t] +
#'   b[x,i] k[t,i]`. The common factor `B K` is a Lee-Carter model of the
#'   pooled rates (sum of deaths over sum of exposures); the stratum factors
#'   `b k` are fitted to the residuals.
#' * `product_ratio()` (Hyndman, Booth and Yasmeen 2013): the geometric mean
#'   of the stratum rates (product) is modelled with a Lee-Carter-type
#'   factor, and the ratios of each stratum to the product with factors whose
#'   scores are stationary. This is a Lee-Carter-type version of the method;
#'   the original uses functional data models.
#'
#' With `common_by` (for example `"sex"`), `li_lee()` fits one common factor
#' per level of that stratum variable (women and men each with their own
#' trend), and the other strata (for example regions and deprivation
#' quintiles) deviate from their group's common factor. If the data carry a
#' hierarchy (see [mt_aggregate_strata()]), the common factors are estimated
#' from the bottom strata only, so aggregates are not counted twice.
#'
#' In both models the common period index follows a random walk with drift
#' (with drift uncertainty) and the stratum-specific scores follow AR(1)
#' processes without drift, so the differences between strata stay bounded.
#' Draws include the innovations of all indexes; predictive draws add
#' extra-Poisson and Poisson variation.
#'
#' @param ratio_order Number of stratum-specific (ratio) factors.
#' @param jump_off `"fitted"` (forecasts start from the fitted rates of the
#'   last year) or `"observed"` (forecasts are multiplied by the ratio of
#'   observed to fitted rates in the last year, per age and stratum; Lee and
#'   Miller 2001). With few factors per stratum the fitted rates of the last
#'   year can be far from the observed ones; the observed jump-off then
#'   avoids a step between the last observed and the first forecast year.
#' @param common_by Stratum variable(s) defining groups with their own common
#'   factor (`li_lee()` only).
#' @return A model specification.
#' @references
#' Li N, Lee R (2005). Coherent mortality forecasts for a group of
#' populations: an extension of the Lee-Carter method. *Demography* 42(3):
#' 575-594. \doi{10.1353/dem.2005.0021}
#'
#' Lee R, Miller T (2001). Evaluating the performance of the Lee-Carter
#' method for forecasting mortality. *Demography* 38(4): 537-549.
#'
#' Hyndman RJ, Booth H, Yasmeen F (2013). Coherent mortality forecasting: the
#' product-ratio method with functional time series models. *Demography*
#' 50(1): 261-283.
#' @examples
#' x <- mt_example(deaths = "single")
#' fc <- mt_forecast(mt_fit(li_lee(), x), h = 10, n_draws = 100, seed = 1)
#' @export
li_lee <- function(ratio_order = 1, common_by = NULL, jump_off = c("fitted", "observed")) {
  jump_off <- match.arg(jump_off)
  new_mt_spec("multipop", list(type = "li_lee", ratio_order = ratio_order, common_by = common_by,
                               jump_off = jump_off), joint = TRUE,
              min_periods = 8,
              label = paste0(if (is.null(common_by)) "Li-Lee" else sprintf("Li-Lee (common by %s)",
                                                                           paste(common_by, collapse = ", ")),
                             if (jump_off == "observed") ", observed jump-off" else ""))
}

#' @rdname li_lee
#' @export
product_ratio <- function(ratio_order = 1, jump_off = c("fitted", "observed")) {
  jump_off <- match.arg(jump_off)
  new_mt_spec("multipop", list(type = "product_ratio", ratio_order = ratio_order, jump_off = jump_off), joint = TRUE,
              min_periods = 8,
              label = paste0("Product-ratio", if (jump_off == "observed") ", observed jump-off" else ""))
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
    if (p$type == "li_lee") {
      lc_svd(apply(D[, , use, drop = FALSE], 1:2, sum), apply(E[, , use, drop = FALSE], 1:2, sum))
    } else {
      lc_from_log(apply(y[, , use, drop = FALSE], 1:2, mean))
    }
  })
  base_of <- function(cm) if (p$type == "li_lee") outer(cm$b, cm$k) else cm$a + outer(cm$b, cm$k)
  r <- min(p$ratio_order, ncol(D) - 2)
  strata_fits <- lapply(seq_len(ns), function(s) {
    res <- y[, , s] - base_of(commons[[group[s]]])
    a <- rowMeans(res)
    sv <- svd(res - a, nu = r, nv = r)
    b <- sv$u[, seq_len(r), drop = FALSE]
    k <- sweep(sv$v[, seq_len(r), drop = FALSE], 2, sv$d[seq_len(r)], "*")
    ar <- lapply(seq_len(r), function(j) ar1_fit(k[, j]))
    list(a = a, b = b, k = k, ar = ar)
  })
  fitted <- array(NA_real_, dim(D))
  for (s in seq_len(ns)) {
    st <- strata_fits[[s]]
    fitted[, , s] <- exp(base_of(commons[[group[s]]]) + st$a + st$b %*% t(st$k))
  }
  nt <- dim(D)[2]
  jump <- if (identical(p$jump_off, "observed")) y[, nt, , drop = TRUE] - log(fitted[, nt, , drop = TRUE]) else NULL
  if (!is.null(jump)) jump <- matrix(jump, dim(D)[1], ns)
  list(commons = commons, group = group, base_type = p$type, strata = strata_fits, fitted = fitted, jump = jump)
}

# Lee-Carter by SVD of a matrix of log rates (no second-stage refit of k)
lc_from_log <- function(y) {
  a <- rowMeans(y)
  s <- svd(y - a, nu = 1, nv = 1)
  b <- s$u[, 1]
  k <- s$d[1] * s$v[, 1]
  if (sum(b) < 0) {
    b <- -b
    k <- -k
  }
  lc_finish(a, b, k)
}

ar1_fit <- function(k) {
  n <- length(k)
  phi <- sum(k[-1] * k[-n]) / sum(k[-n]^2)
  phi <- max(min(phi, 0.99), -0.99)
  e <- k[-1] - phi * k[-n]
  list(phi = phi, sigma = sqrt(sum(e^2) / max(n - 2, 1)), last = k[n])
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
      if (fit$base_type == "product_ratio") base <- base + cm$a
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
