# Indices of lifespan variation.
# Adapted from LifeIneq 00.04.09 (2024-03-07; ineq_edag(), ineq_var(),
# ineq_sd(), ineq_cov(), ineq_H(), ineq_gini() and ineq_iqr() in
# R/indices.R; https://github.com/alysonvanraalte/LifeIneq), GPL-3, by
# Alyson van Raalte, Tim Riffe and Christian Dudel. Changes: computed for
# many life tables at once (the columns of a matrix) at one conditioning
# age; the inputs come from lt_columns(); the variance uses the distance of
# each age at death from the conditioning age (LifeIneq's ineq_var() uses the
# distance of the first ages from the first age, which differs for abridged
# tables); the Gini coefficient uses the sorted-values form of the double sum
# (identical result, linear instead of quadratic in the number of ages); the
# quartiles interpolate the survival curve linearly within closed intervals
# and exponentially in the open interval (LifeIneq: a monotone spline through
# the survival curve after a random perturbation of tied values).

#' Lifespan variation from rates or forecast draws
#'
#' Indices of the variation of the age at death in the period life table,
#' conditional on survival to age `at`: how unequal lifespans are, as
#' opposed to how long they are on average (life expectancy). For forecast
#' draws they are computed path by path, like [life_expectancy()].
#'
#' The ages at death of the life table are taken as `x + a_x` for the deaths
#' `d_x` of each interval (the discrete formulas of LifeIneq). With `l`
#' the survivors and `e` the life expectancy at age `at`, `y_x = x + a_x`
#' and `p_x = d_x / l(at)` for the intervals from `at`:
#'
#' * `"edagger"`: life disparity, the average remaining life expectancy lost
#'   at death (Vaupel 1986), `sum(d_x * ebar_x) / l(at)`, where `ebar_x`
#'   interpolates life expectancy linearly to the average age at death in
#'   the interval, `e_x + a_x / n_x * (e_{x+n} - e_x)` (the open interval
#'   uses its `e_x`);
#' * `"var"`, `"sd"`: variance and standard deviation of the age at death,
#'   `sum(p_x * (y_x - at - e)^2)`;
#' * `"cv"`: coefficient of variation, the standard deviation divided by the
#'   mean age at death `at + e` (or by `e` with `distribution =
#'   "remaining_life"`);
#' * `"entropy"`: the life-table entropy of Keyfitz (1977) and Leser (1955),
#'   `H = e-dagger / e`, the elasticity of life expectancy to a
#'   proportional change of the rates at all ages (Goldman and Lord 1986);
#' * `"gini"`: Gini coefficient of the ages at death, the mean absolute
#'   difference between two ages at death divided by twice the mean age at
#'   death `at + e` (Hanada 1983; Shkolnikov, Andreev and Begun 2003), or by
#'   twice the mean remaining life `e` with `distribution =
#'   "remaining_life"`;
#' * `"iqr"`: interquartile range, the distance between the ages by which a
#'   quarter and three quarters of the survivors at age `at` have died.
#'   The survival curve is interpolated linearly within closed intervals
#'   (deaths spread evenly, as with the default `a_x = n / 2`) and
#'   exponentially in the open interval (with the mean `a_x` of that
#'   interval).
#'
#' The absolute measures (`edagger`, `sd`, `iqr`) are in years, `var` in
#' years squared; `cv`, `entropy` and `gini` have no unit.
#'
#' The discrete formulas treat everyone dying in an interval as dying at age
#' `x + a_x`, so `var`, `sd`, `cv` and `gini` leave out the spread of the
#' ages at death within each interval (about `n^2 / 12` years squared for
#' deaths spread evenly over `n` years). Compared with ages at death spread
#' evenly within closed intervals and exponentially in the open one, for
#' French life tables (HMD) the standard deviation is lower by less than
#' 0.01 years with single years of age to 110+, by up to 0.14 years with
#' 5-year groups to 100+, by 0.3 to 0.6 years with an open interval from 90
#' and by 1 to 2 years with one from 85 (Gini coefficient: up to 0.002 lower
#' with 5-year groups to 100+, 0.016 with 85+). `edagger`, `entropy`
#' and `iqr` are not affected by the open interval (remaining life
#' expectancy does not change within an interval of constant hazard). Use
#' single years of age up to 100 or more where possible (for example from
#' [mt_ungroup()]).
#'
#' The life table needs an open last age interval (otherwise the deaths do
#' not add up to the survivors at age `at`). The open interval uses the
#' life-table closure (`open`, see [life_table()]); its average time lived
#' `a_x` is its life expectancy.
#'
#' @param x A matrix or array of rates with age in the first dimension, an
#'   [mt_draws] object, or an [mt_data] object (observed rates with Poisson
#'   intervals, see [mt_observed()]).
#' @param measure One or more of `"edagger"`, `"sd"`, `"var"`, `"cv"`,
#'   `"entropy"`, `"gini"`, `"iqr"` (one only for arrays and `paths =
#'   TRUE`).
#' @param age An [mt_axis] (taken from `x` for draws).
#' @param at Age on which the measures are conditioned (start of an
#'   interval; default the first age).
#' @param distribution Denominator of the relative measures `"cv"` and
#'   `"gini"`: the mean age at death (`"age_at_death"`) or the mean
#'   remaining life (`"remaining_life"`) at age `at`.
#' @param probs Quantiles reported for draws.
#' @param paths For draws, return the values of every path, an array
#'   `[period, stratum, draw]`.
#' @param ... Passed to the life-table calculation (`ax`, `open`, `k`), see
#'   [life_table()].
#' @return For arrays, an array without the age dimension. For draws and
#'   data, a `data.table` with the stratum columns, `period`, `mean`,
#'   quantiles, `measure` and `at`, computed path by path.
#' @references
#' Goldman N, Lord G (1986). A new look at entropy and the life table.
#' *Demography* 23(2): 275-282.
#'
#' Hanada K (1983). A formula of Gini's concentration ratio and its
#' application to life tables. *Journal of the Japan Statistical Society,
#' Japanese Issue* 13(2): 95-98.
#'
#' Keyfitz N (1977). Mortality comparisons; the male-female ratio. In
#' *Applied Mathematical Demography*, 54-76. John Wiley and Sons.
#'
#' Leser CEV (1955). Variations in mortality and life expectation.
#' *Population Studies* 9(1): 67-71.
#'
#' Shkolnikov VM, Andreev EM, Begun AZ (2003). Gini coefficient as a life
#' table function: computation from discrete data, decomposition of
#' differences and empirical examples. *Demographic Research* 8: 305-358.
#'
#' van Raalte A, Riffe T, Dudel C (2024). LifeIneq: Calculate indices of
#' lifespan variation. R package version 00.04.09.
#' <https://github.com/alysonvanraalte/LifeIneq>
#'
#' Vaupel JW (1986). How change in age-specific mortality affects life
#' expectancy. *Population Studies* 40(1): 147-157.
#' @examples
#' m <- 5e-5 * exp(0.095 * (0:100))
#' lifespan_variation(m, single_ages(0, 100), measure = "edagger")
#' x <- mt_example(deaths = "single", years = 2015:2019)
#' lifespan_variation(x, measure = c("edagger", "gini"), at = 65)
#' fc <- mt_forecast(mt_fit(lee_carter(), x), h = 5, n_draws = 100, seed = 1)
#' lifespan_variation(fc, measure = "sd", at = 30)
#' @export
lifespan_variation <- function(x, measure = "edagger", age = NULL, at = NULL,
                               distribution = c("age_at_death", "remaining_life"),
                               probs = c(0.025, 0.1, 0.5, 0.9, 0.975), paths = FALSE, ...) {
  distribution <- match.arg(distribution)
  measure <- check_lifespan_measure(measure)
  if (inherits(x, "mt_data")) x <- mt_observed(x)
  draws <- inherits(x, "mt_draws")
  if (length(measure) > 1 && (!draws || paths)) {
    mt_abort("Give one `measure` for an array of rates or with `paths = TRUE`.", class = "argument")
  }
  if (draws) {
    age <- x$age
    rates <- x$rates
  } else {
    check_class(age, "mt_axis", "age")
    rates <- x
  }
  d <- dim(rates) %||% length(rates)
  M <- matrix(rates, nrow = d[1])
  if (nrow(M) != length(age)) mt_abort("`x` needs one row per age interval of `age`.", class = "lifetable")
  i <- if (is.null(at)) 1L else match(at, age$lower)
  if (length(i) != 1 || is.na(i)) mt_abort("`at` must be the start of an age interval.", class = "lifetable")
  v <- lifespan_columns(M, age, i, measure, distribution, ...)
  if (!draws) return(if (length(d) > 1) array(v[[1]], d[-1], dimnames = dimnames(x)[-1]) else v[[1]])
  if (paths) {
    return(array(v[[1]], d[-1], dimnames = list(period = x$periods, stratum = x$strata$label, draw = NULL)))
  }
  out <- rbindlist(lapply(measure, function(ms) {
    r <- summarise_paths(array(v[[ms]], d[-1]), x, probs)
    r[, `:=`(measure = ms, at = age$lower[i])]
  }))
  out[]
}

lifespan_measures <- c("edagger", "sd", "var", "cv", "entropy", "gini", "iqr")

check_lifespan_measure <- function(measure) {
  if (!is.character(measure) || !length(measure) || !all(measure %in% lifespan_measures)) {
    mt_abort(sprintf("`measure` must be one or more of %s.", paste0("\"", lifespan_measures, "\"", collapse = ", ")),
             class = "argument")
  }
  unique(measure)
}

# Indices of lifespan variation (a list with one vector per measure) for the
# life tables of the columns of M (rates [age, schedule]), conditional on
# survival to the start of age interval i
lifespan_columns <- function(M, age, i, measure, distribution = "age_at_death", ...) {
  if (!is_open(age)) {
    mt_abort("Lifespan variation needs an open last age interval, so that every death is in the life table.",
             class = "lifetable")
  }
  na <- nrow(M)
  r <- lt_columns(M, age, ...)
  rows <- i:na
  k <- length(rows)
  x <- age$lower[rows]
  w <- (age$upper - age$lower)[rows]
  l <- r$lx[rows, , drop = FALSE]
  d <- r$dx[rows, , drop = FALSE]
  a <- r$ax[rows, , drop = FALSE]
  e <- r$ex[rows, , drop = FALSE]
  # open interval: time lived by those who die in it (1 / m for the reciprocal
  # closure, the life expectancy of the extension for "gompertz"/"kannisto")
  a[k, ] <- ifelse(d[k, ] > 0, r$Lx[na, ] / d[k, ], a[k, ])
  l0 <- l[1, ]
  e0 <- e[1, ]
  na_out <- rep(NA_real_, ncol(M))
  out <- stats::setNames(rep(list(na_out), length(measure)), measure)
  ok <- colSums(!is.finite(l) | !is.finite(d) | !is.finite(a) | !is.finite(e)) == 0 & l0 > 0
  if (!any(ok)) return(out)
  l <- l[, ok, drop = FALSE]
  d <- d[, ok, drop = FALSE]
  a <- a[, ok, drop = FALSE]
  e <- e[, ok, drop = FALSE]
  l0 <- l0[ok]
  e0 <- e0[ok]
  p <- d / rep(l0, each = k)
  y <- x + a
  mu <- x[1] + e0
  centre <- if (distribution == "age_at_death") mu else e0
  if (any(c("var", "sd", "cv") %in% measure)) variance <- colSums(p * (y - rep(mu, each = k))^2)
  if (any(c("edagger", "entropy") %in% measure)) {
    e_next <- rbind(e[-1, , drop = FALSE], e[k, , drop = FALSE])
    ebar <- e + a / w * (e_next - e)
    edagger <- colSums(d * ebar) / l0
  }
  for (ms in measure) {
    out[[ms]][ok] <- switch(ms,
                            edagger = edagger,
                            entropy = edagger / e0,
                            var = variance,
                            sd = sqrt(variance),
                            cv = sqrt(variance) / centre,
                            gini = gini_columns(p, y) / centre,
                            iqr = lifespan_quantile(l, x, w, a, 0.25) - lifespan_quantile(l, x, w, a, 0.75))
  }
  out
}

# Half the mean absolute difference between two ages at death,
# sum_{j<k} p_j p_k |y_j - y_k|, per column (sorted form; the double sum where
# the ages at death of a column are not increasing)
gini_columns <- function(p, y) {
  k <- nrow(p)
  if (k == 1) return(rep(0, ncol(p)))
  cum <- matrix(apply(p, 2, cumsum), nrow = k)
  cum_prev <- rbind(0, cum[-k, , drop = FALSE])
  s <- colSums(p * y * (cum_prev + cum - 1))
  unsorted <- which(colSums(diff(y) < 0) > 0)
  for (j in unsorted) s[j] <- sum(abs(outer(y[, j], y[, j], "-")) * outer(p[, j], p[, j])) / 2
  s
}

# Age by which a share 1 - surv of the survivors at the first age has died:
# linear interpolation of the survival curve in closed intervals, exponential
# (mean time lived a) in the open last interval
lifespan_quantile <- function(l, x, w, a, surv) {
  k <- nrow(l)
  s <- l / rep(l[1, ], each = k)
  j <- pmax(colSums(s >= surv), 1L)
  cols <- seq_len(ncol(l))
  sj <- s[cbind(j, cols)]
  out <- numeric(ncol(l))
  closed <- j < k
  if (any(closed)) {
    jc <- j[closed]
    s_next <- s[cbind(jc + 1L, cols[closed])]
    out[closed] <- x[jc] + w[jc] * (sj[closed] - surv) / (sj[closed] - s_next)
  }
  if (any(!closed)) out[!closed] <- x[k] + a[cbind(k, cols[!closed])] * log(sj[!closed] / surv)
  out
}
