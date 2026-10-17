#' Observed rates as draws
#'
#' Turns observed deaths and exposures into an [mt_draws] object whose paths
#' are Poisson draws of the death counts given the exposures. All functions
#' for draws then also describe the past, with intervals for the observed
#' rates: [life_expectancy()], [asr()], [prob_dying()], [mt_contrast()],
#' [mt_inequality()]. The derived quantities accept an [mt_data] object
#' directly and call `mt_observed()` with its defaults.
#'
#' @param data An [mt_data] object with deaths and exposures on the grid.
#' @param n_draws Number of draws.
#' @param seed Random seed (default 1, so that results are reproducible).
#' @return An `mt_draws` object with the observed periods; draw 1 holds the
#'   observed rates.
#' @examples
#' x <- mt_example(deaths = "single", years = 2015:2019)
#' life_expectancy(x, at = 65)
#' @export
mt_observed <- function(data, n_draws = 500, seed = 1) {
  check_class(data, "mt_data", "data")
  md <- model_inputs(data, list(min_periods = 1, label = "Observed rates", consecutive = FALSE))
  D <- md$D
  E <- md$E
  arr <- with_seed(seed, array(stats::rpois(length(D) * n_draws, rep(as.vector(D), n_draws)),
                               c(dim(D), n_draws)))
  arr[, , , 1] <- D
  arr <- arr / as.vector(E)
  fit <- list(observed = D / E, periods = md$periods, train_D = D, train_E = E,
              exposure_last = E[, length(md$periods), , drop = FALSE])
  out <- new_mt_draws(arr, md$age, md$periods, data$strata, model = "Observed rates", origin = NA, fit = fit,
                      seed = seed)
  out$exposure <- E
  out
}

# Draws from draws or from data
as_draws <- function(x) {
  if (inherits(x, "mt_data")) return(mt_observed(x))
  check_class(x, "mt_draws", "draws")
  x
}

#' Differences and ratios between strata, path by path
#'
#' Compares a derived quantity (life expectancy, age-standardised rate,
#' probability of dying, or death rate of one age group) between the levels
#' of a stratum variable and a reference level, path by path, so that the
#' interval of the difference includes the joint uncertainty of both
#' strata. Strata are matched on all other stratum variables (for example
#' the deprivation gap within each sex).
#'
#' @param x An [mt_draws] object (forecast) or an [mt_data] object
#'   (observed, see [mt_observed()]).
#' @param var Stratum variable to compare across (for example `"quintile"`).
#' @param ref Reference level of `var`.
#' @param quantity `"life_expectancy"`, `"asr"`, `"prob_dying"` or
#'   `"rate"`.
#' @param type `"difference"` (level minus reference) or `"ratio"`.
#' @param probs Quantile levels.
#' @param ... Arguments of the quantity: `at` for life expectancy,
#'   `standard` and `per` for `asr()`, `from` and `to` for `prob_dying()`,
#'   `age` (an age label) for `"rate"`.
#' @return A `data.table` with the other stratum variables, `level`, `ref`,
#'   `period`, `mean`, quantiles and `p_positive`, the share of paths with
#'   a positive difference (or a ratio above 1).
#' @examples
#' x <- mt_example(deaths = "single", years = 2001:2015)
#' fc <- mt_forecast(mt_fit(lee_carter(), x), h = 5, n_draws = 100, seed = 1)
#' mt_contrast(fc, var = "sex", ref = "female", at = 30)
#' @export
mt_contrast <- function(x, var, ref, quantity = c("life_expectancy", "asr", "prob_dying", "rate"),
                        type = c("difference", "ratio"), probs = c(0.025, 0.5, 0.975), ...) {
  quantity <- match.arg(quantity)
  type <- match.arg(type)
  x <- as_draws(x)
  st <- x$strata
  if (!var %in% names(st)) {
    mt_abort(sprintf("`var` must be a stratum variable (%s).", paste(setdiff(names(st), c("stratum", "label")),
                                                                     collapse = ", ")), class = "argument")
  }
  if (!ref %in% st[[var]]) {
    mt_abort(sprintf("`ref` must be a level of `%s` (%s).", var, paste(unique(st[[var]]), collapse = ", ")),
             class = "argument")
  }
  q <- quantity_paths(x, quantity, ...)
  others <- setdiff(names(st), c("stratum", "label", var))
  key <- if (length(others)) do.call(paste, c(lapply(others, function(o) st[[o]]), sep = "\r")) else rep("", nrow(st))
  rows <- list()
  for (i in which(st[[var]] != ref)) {
    j <- which(st[[var]] == ref & key == key[i])
    if (length(j) != 1) next
    v <- if (type == "difference") q[, i, , drop = FALSE] - q[, j, , drop = FALSE] else
      q[, i, , drop = FALSE] / q[, j, , drop = FALSE]
    M <- matrix(v, nrow = dim(q)[1])
    Q <- row_quantiles(M, probs)
    r <- data.table(level = st[[var]][i], ref = ref, period = x$periods, mean = rowMeans(M))
    for (k in seq_along(probs)) set(r, j = qnames(probs)[k], value = Q[, k])
    r[, p_positive := rowMeans(M > if (type == "difference") 0 else 1)]
    if (length(others)) r <- cbind(st[i, others, with = FALSE][rep(1, nrow(r))], r)
    rows[[length(rows) + 1L]] <- r
  }
  if (!length(rows)) mt_abort("No strata to compare with the reference.", class = "argument")
  out <- rbindlist(rows)
  out[, `:=`(quantity = quantity, type = type)]
  out[]
}

# Paths [period, stratum, draw] of a derived quantity
quantity_paths <- function(x, quantity, ...) {
  a <- list(...)
  switch(quantity,
         life_expectancy = {
           e <- life_expectancy(x$rates, x$age, at = a$at)
           e
         },
         asr = asr(x, standard = a$standard %||% "esp2013", per = a$per %||% 1e5, paths = TRUE),
         prob_dying = {
           if (is.null(a$from) || is.null(a$to)) mt_abort("Give `from` and `to` for prob_dying.", class = "argument")
           prob_dying(x, from = a$from, to = a$to, paths = TRUE)
         },
         rate = {
           i <- match(a$age %||% NA, x$age$labels)
           if (is.na(i)) {
             mt_abort(sprintf("Give `age`, one of the age labels (%s).", paste(x$age$labels, collapse = ", ")),
                      class = "argument")
           }
           array(x$rates[i, , , ], dim(x$rates)[-1])
         })
}

#' Slope and relative index of inequality, path by path
#'
#' Summarises a gradient over the ordered levels of a stratum variable (for
#' example deprivation quintiles) by a regression of the quantity on the
#' relative rank of each level in the population, weighted by population
#' size: the slope index of inequality (SII) is the difference between the
#' fitted values at the top (rank 1) and the bottom (rank 0) of the
#' population, the relative index of inequality (RII) their ratio. The
#' relative rank of a level is the midpoint of its cumulative population
#' share, with the levels ordered as in `order`. Computed path by path, with
#' the population shares of the exposures of each year for observed rates
#' and of the last observed year for forecasts.
#'
#' @param x An [mt_draws] or [mt_data] object.
#' @param var Ordered stratum variable (for example `"quintile"`).
#' @param order Levels of `var` from rank 0 to rank 1 (default: their order
#'   in the data, without the total and the groups of a hierarchy stored in
#'   `x` by [mt_aggregate_strata()], [reconciled()] or
#'   `mt_group_levels(total = )`). For deprivation, list the least deprived
#'   level first so that a positive SII means higher values in more deprived
#'   groups. Levels not listed, such as the total of a hierarchy, are left
#'   out.
#' @inheritParams mt_contrast
#' @return A `data.table` with the other stratum variables, `index` (`"sii"`
#'   or `"rii"`), `period`, `mean` and quantiles.
#' @references
#' Mackenbach JP, Kunst AE (1997). Measuring the magnitude of
#' socio-economic inequalities in health: an overview of available measures
#' illustrated with two examples from Europe. *Social Science & Medicine*
#' 44(6): 757-771.
#' @examples
#' x <- mt_example_regions(years = 2010:2019, regions = 2)
#' mt_inequality(x, "quintile", order = paste0("Q", 5:1), quantity = "asr")
#' @export
mt_inequality <- function(x, var, order = NULL, quantity = c("asr", "life_expectancy", "prob_dying", "rate"),
                          probs = c(0.025, 0.5, 0.975), ...) {
  quantity <- match.arg(quantity)
  h <- object_hierarchy(x)
  x <- as_draws(x)
  st <- x$strata
  if (!var %in% names(st)) mt_abort("`var` must be a stratum variable.", class = "argument")
  if (is.null(order)) {
    # totals and groups of a hierarchy are not levels of the gradient
    hv <- if (!is.null(h)) h$vars[[var]]
    order <- setdiff(unique(st[[var]]), c(hv$total, names(hv$groups)))
  }
  if (anyDuplicated(order) || !all(order %in% st[[var]]) || length(order) < 2) {
    mt_abort(sprintf("`order` must list at least two levels of `%s`, each once.", var), class = "argument")
  }
  q <- quantity_paths(x, quantity, ...)
  # population weights per year: the exposures of each year for observed
  # rates, those of the last observed year for forecasts
  np <- length(x$periods)
  pop <- if (!is.null(x$exposure) && dim(x$exposure)[2] == np) {
    apply(x$exposure, c(2, 3), sum)
  } else {
    E <- x$exposure_last
    matrix(colSums(matrix(E, nrow = dim(E)[1])), np, nrow(st), byrow = TRUE)
  }
  others <- setdiff(names(st), c("stratum", "label", var))
  key <- if (length(others)) do.call(paste, c(lapply(others, function(o) st[[o]]), sep = "\r")) else rep("", nrow(st))
  rows <- list()
  for (k in unique(key)) {
    sel <- which(key == k & st[[var]] %in% order)
    idx <- sel[order(match(st[[var]][sel], order))]
    if (length(idx) < 2) next
    sii <- rii <- matrix(NA_real_, dim(q)[1], dim(q)[3])
    for (t in seq_len(dim(q)[1])) {
      w <- pop[t, idx] / sum(pop[t, idx])
      rank <- cumsum(w) - w / 2
      X <- cbind(1, rank)
      H <- solve(crossprod(X, X * w), t(X * w))   # weighted least squares
      b <- H %*% matrix(q[t, idx, ], length(idx))
      sii[t, ] <- b[2, ]
      rii[t, ] <- (b[1, ] + b[2, ]) / b[1, ]
    }
    for (ix in c("sii", "rii")) {
      M <- if (ix == "sii") sii else rii
      Q <- row_quantiles(M, probs)
      r <- data.table(index = ix, period = x$periods, mean = rowMeans(M))
      for (j in seq_along(probs)) set(r, j = qnames(probs)[j], value = Q[, j])
      if (length(others)) r <- cbind(st[idx[1], others, with = FALSE][rep(1, nrow(r))], r)
      rows[[length(rows) + 1L]] <- r
    }
  }
  out <- rbindlist(rows)
  out[, quantity := quantity]
  out[]
}
