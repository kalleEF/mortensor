#' Decompose differences in mortality by age, cause and subgroup
#'
#' Splits the difference in a mortality summary between two sets of rates
#' into contributions of ages, causes of death and population subgroups. It
#' is a post-processing step for observed data, forecasts, or both:
#'
#' * with `var` and `ref`, each level of a stratum variable is compared with
#'   a reference level (for example the deprivation gap Q5 - Q1 in each sex),
#'   in each of `periods`;
#' * with `y`, `y` is compared with `x` for matching strata and periods (two
#'   forecasts, or observed data and a forecast);
#' * otherwise the change from the first of `periods` to each later one is
#'   decomposed for every stratum. For forecast draws, `periods` may mix
#'   observed years (the training data held in the draws) and forecast
#'   years.
#'
#' Forecasts are decomposed draw by draw, so every contribution has an
#' interval, and sums over ages, causes or subgroups ([summary()]) are formed
#' path by path. Observed data give the decomposition of the observed rates
#' (`n_draws = 1`) or, with `n_draws > 1`, Poisson draws of the deaths as in
#' [mt_observed()]. Totals and contributions are means over draws;
#' `summary()` adds quantiles.
#'
#' @section Observed and forecast years:
#' Forecast draws hold the observed rates of their training years. A change
#' from an observed year to a forecast year then contains the step from the
#' observed to the model's fitted rates of that year (the jump-off, see
#' [mt_jump_off()]) and the sampling noise of the observed year, besides
#' the forecast change itself. The totals report the jump-off share
#' (`jump_off`). `observed = "fitted"` uses the fitted rates of observed
#' years, so that the decomposition describes the model's change alone;
#' `observed = "poisson"` uses Poisson draws of the observed deaths, so that
#' their sampling uncertainty enters the intervals.
#'
#' Two objects (`y`) are compared draw by draw. Draws of different models
#' are independent, so by default (`pairing = "independent"`) the draws of
#' `y` are shuffled before pairing; draws of models with the same seed would
#' otherwise share their random numbers and give intervals that are too
#' narrow. Use `pairing = "draw"` when `y` is derived from the draws of `x`
#' (for example a reconciliation of the same draws).
#'
#' @section Methods:
#' * `"arriaga"` (life expectancy): the symmetric Arriaga decomposition, the
#'   average of the decompositions from 1 to 2 and (with the sign reversed)
#'   from 2 to 1, computed on the package's life table (`...` passes `ax`
#'   and `open` to it). Causes get the age-specific contribution in
#'   proportion to their rate differences (an implied sensitivity, Riffe et
#'   al. 2024, eqs. 12-13 of the working paper); where the all-cause
#'   difference is close to zero, the sensitivity is a numerical derivative
#'   at the average rates. Exact.
#' * `"horiuchi"` (any quantity): the line-integral method of Horiuchi,
#'   Wilmoth and Pletcher (2008) with `n_steps` steps; exact up to a small
#'   discretisation error, reported in the totals.
#' * `"stepwise"` (any quantity): stepwise replacement (Andreev, Shkolnikov
#'   and Begun 2002) from the youngest to the oldest age, averaged over
#'   both directions of the comparison; with causes, all causes of an age
#'   are replaced before the next age (DemoDecomp and LEdecomp replace cause
#'   by cause, which changes the split between causes slightly). Exact;
#'   results depend on the order. Arriaga's decomposition from 1 to 2 is
#'   stepwise replacement from the oldest to the youngest age, so for life
#'   expectancy without causes the two symmetric versions are identical.
#' * `"kitagawa"` (crude death rate): the rate effect (age-specific rates at
#'   the average age structure) and the composition effect (age structure
#'   at the average rates) of Kitagawa (1955). Exact.
#'
#' Draws whose rates give no valid life table (for example rates far above
#' any observed value, see the warning of [mt_forecast()]) are left out with
#' a warning, and a warning reports contributions that do not add up to
#' the difference.
#'
#' @section Subgroups:
#' With `groups`, the populations being compared consist of the levels of a
#' stratum variable, and the decomposition separates the effect of their
#' composition. `composition` must be chosen, because the two options
#' describe different quantities:
#'
#' * `"radix"` (Riffe, Tursun-zade and Trias Llimós 2024, "Arriaga meets
#'   Kitagawa"): life expectancy (or the probability of dying) of the
#'   population is the mean of the subgroups' values, weighted by their
#'   population shares at the starting age (`at`), as when subgroup life
#'   tables are blended at the radix ("current conditions" in the terms of
#'   Vaupel 2002, which the authors use). A Kitagawa decomposition of that
#'   mean gives one composition effect and, for each subgroup, its
#'   difference times its average share, which is split by age (and cause)
#'   with `method`. The authors state that this holds when membership is
#'   fixed from the starting age on and later shares change only through
#'   mortality: education from age 30 or 35, sex, country of birth. It does
#'   not hold for areas such as deprivation quintiles or regions, between
#'   which people move and which are reclassified.
#' * `"pooled"`: the quantity of the pooled rates, which weight the
#'   subgroup rates by their population shares at each age (the rates of the
#'   aggregate strata, as published for whole populations; "current
#'   rates"). The rates of every subgroup, age and cause and the shares at
#'   every age are the parameters of a `"horiuchi"` (default) or
#'   `"stepwise"` decomposition; the composition effects are reported by
#'   age, summed over subgroups, whose separate composition terms are not
#'   identified. With many subgroups, causes and draws this takes minutes;
#'   lower `n_steps` or the number of draws for a first look.
#'
#' Shares come from the exposures: of each observed year, and for forecast
#' years from `exposure` (default: the last observed year, so composition
#' does not change in the forecast unless exposures, for example a
#' population projection, are given).
#'
#' @section Totals among causes and groups:
#' Levels of `cause` and `groups` that are totals of the others (for example
#' the cause `"all"` or the region `"Germany"`) must not be counted as
#' parts. They are left out when the hierarchy declares them: `hierarchy`,
#' or what the object records, the total of [mt_group_levels()] (`total`,
#' also in fits and draws of such data) and the hierarchy of
#' [mt_aggregate_strata()] or [reconciled()]. A level that equals the sum
#' of the others in the observed deaths without being declared stops the
#' decomposition.
#'
#' @param x An [mt_data] object (observed) or an [mt_draws] object
#'   (forecast).
#' @param y Optional second [mt_data] or [mt_draws] object with the same
#'   age grid and strata.
#' @param var,ref Stratum variable and its reference level for comparisons
#'   between strata.
#' @param periods Years to compare (see Details). Defaults: the last year
#'   for comparisons between strata or objects; for changes, the first and
#'   last observed year of data, or the last observed and the last forecast
#'   year of draws.
#' @param quantity `"life_expectancy"`, `"prob_dying"`, `"asr"`,
#'   `"crude_rate"`, or a function `f(m, age)` that takes a matrix of rates
#'   `[age, n]` (ages of the grid, `age` the [mt_axis]) and returns `n`
#'   values.
#' @param method `"arriaga"`, `"horiuchi"`, `"stepwise"` or `"kitagawa"`
#'   (default: `"arriaga"` for life expectancy, `"kitagawa"` for the crude
#'   rate, `"horiuchi"` otherwise and for pooled composition).
#' @param cause Stratum variable holding causes of death, whose rates add
#'   up to all causes.
#' @param groups Stratum variable holding population subgroups.
#' @param hierarchy An [mt_hierarchy()] that declares total and group levels
#'   of `cause` and `groups` (default: the hierarchy stored in `x` or `y`).
#' @param composition `"radix"` or `"pooled"`; required with `groups` (see
#'   Subgroups).
#' @param observed How observed years held in forecast draws enter:
#'   `"data"` (observed rates), `"fitted"` (the model's fitted rates) or
#'   `"poisson"` (Poisson draws of the observed deaths). See Observed and
#'   forecast years.
#' @param pairing `"independent"` or `"draw"`: how the draws of `x` and `y`
#'   are paired.
#' @param at Starting age (life expectancy at `at`; lower bound of
#'   `prob_dying`); default the first grid age.
#' @param to Upper age for `"prob_dying"`.
#' @param standard,per Standard population and multiplier for `"asr"`
#'   (see [asr()]); `per` also scales `"crude_rate"`.
#' @param exposure Exposures for forecast years `[age, stratum]` or
#'   `[age, period, stratum]` (see [mt_deaths()]).
#' @param n_steps Number of steps of `"horiuchi"`.
#' @param n_draws For [mt_data]: 1 (the observed rates) or the number of
#'   Poisson draws.
#' @param probs Quantile levels of the summaries.
#' @param seed Seed of Poisson draws and of the shuffling of `y`.
#' @param ... Passed to the life table (`ax`, `open`, `k`, see
#'   [life_table()]) for life expectancy.
#' @return An `mt_decomposition` object. `summary()` returns contributions
#'   summed (path by path) over the variables not in `by`, with mean and
#'   quantiles; `summary(by = NULL)` returns one row per comparison (its
#'   number in `comparison`) with the values in the two populations, their
#'   difference (mean and quantiles), the jump-off share where it applies,
#'   the sum of the contributions and the largest difference between the
#'   two over draws (`max_residual`). `summary(change = c(p1, p2))`
#'   decomposes the change between two periods of comparisons between
#'   strata or objects (for example of a gap), path by path. `print()`
#'   shows the totals.
#' @references
#' Arriaga EE (1984). Measuring and explaining the change in life
#' expectancies. *Demography* 21(1): 83-96. \doi{10.2307/2061029}
#'
#' Kitagawa EM (1955). Components of a difference between two rates.
#' *Journal of the American Statistical Association* 50(272): 1168-1194.
#' \doi{10.1080/01621459.1955.10501299}
#'
#' Andreev EM, Shkolnikov VM, Begun AZ (2002). Algorithm for decomposition
#' of differences between aggregate demographic measures and its application
#' to life expectancies, healthy life expectancies, parity-progression
#' ratios and total fertility rates. *Demographic Research* 7: 499-522.
#'
#' Horiuchi S, Wilmoth JR, Pletcher SD (2008). A decomposition method based
#' on a model of continuous change. *Demography* 45(4): 785-801.
#' \doi{10.1353/dem.0.0033}
#'
#' Riffe T, Tursun-zade R, Trias Llimós S (2024). Arriaga meets Kitagawa.
#' Life expectancy decomposition with population subgroups. MPIDR Working
#' Paper WP-2024-029. \doi{10.4054/MPIDR-WP-2024-029}. Published in
#' *Geographical Analysis* (2026), \doi{10.1111/gean.70033}.
#'
#' Vaupel JW (2002). Life expectancy at current rates vs. current
#' conditions: a reflexion stimulated by Bongaarts and Feeney's "How long
#' do we live?". *Demographic Research* 7: 365-378.
#' @examples
#' x <- mt_example(deaths = "single", years = 2000:2015)
#' # change 2000-2015 in life expectancy at 30, by age, per sex
#' d <- mt_decompose(x, periods = c(2000, 2015), at = 30)
#' d
#' summary(d, by = "age")[1:5]
#' # the sex gap in forecast years, with intervals, and its change
#' fc <- mt_forecast(mt_fit(lee_carter(), x), h = 10, n_draws = 100, seed = 1)
#' g <- mt_decompose(fc, var = "sex", ref = "male", periods = c(2016, 2025), at = 30)
#' summary(g, by = NULL)
#' summary(g, by = "age", change = c(2016, 2025))[1:5]
#' # observed 2015 against the forecast of 2025, with the jump-off share
#' mt_decompose(fc, periods = c(2015, 2025), at = 30)
#' @export
mt_decompose <- function(x, y = NULL, var = NULL, ref = NULL, periods = NULL,
                         quantity = c("life_expectancy", "prob_dying", "asr", "crude_rate"), method = NULL,
                         cause = NULL, groups = NULL, hierarchy = NULL, composition = NULL,
                         observed = c("data", "fitted", "poisson"), pairing = c("independent", "draw"),
                         at = NULL, to = NULL, standard = "esp2013", per = 1e5, exposure = NULL, n_steps = 20,
                         n_draws = 1, probs = c(0.025, 0.5, 0.975), seed = 1, ...) {
  observed <- match.arg(observed)
  pairing <- match.arg(pairing)
  custom <- is.function(quantity)
  qname <- if (custom) "custom" else match.arg(quantity)
  grouped <- !is.null(groups) && qname != "crude_rate"
  if (grouped && is.null(composition)) {
    mt_abort(paste("With `groups`, choose `composition`: \"pooled\" for the quantity of the pooled rates (areas,",
                   "regions, deprivation groups), or \"radix\" for subgroups blended at the starting age",
                   "(membership fixed from `at` on, as education or sex; Riffe et al. 2024)."), class = "argument")
  }
  composition <- if (grouped) match.arg(composition, c("radix", "pooled"))
  pooled <- identical(composition, "pooled")
  method <- method %||% switch(qname, life_expectancy = if (pooled) "horiuchi" else "arriaga",
                               crude_rate = "kitagawa", "horiuchi")
  check_choice(method, c("arriaga", "horiuchi", "stepwise", "kitagawa"), "method")
  if (method == "arriaga" && qname != "life_expectancy") {
    mt_abort("`method = \"arriaga\"` decomposes life expectancy; use \"horiuchi\" or \"stepwise\".",
             class = "argument")
  }
  if ((method == "kitagawa") != (qname == "crude_rate")) {
    mt_abort("The crude death rate is decomposed with `method = \"kitagawa\"`, which applies to it only.",
             class = "argument")
  }
  if (grouped) {
    if (composition == "radix" && !qname %in% c("life_expectancy", "prob_dying")) {
      mt_abort("`composition = \"radix\"` applies to life expectancy and prob_dying; use \"pooled\".",
               class = "argument")
    }
    if (pooled && method == "arriaga") {
      mt_abort("`composition = \"pooled\"` needs `method = \"horiuchi\"` or \"stepwise\".", class = "argument")
    }
  }
  check_number(n_steps, "n_steps", lower = 1, integerish = TRUE)
  check_number(n_draws, "n_draws", lower = 1, integerish = TRUE)
  if (!is.null(var) && !is.null(y)) mt_abort("Give either `var` or `y`.", class = "argument")

  sx <- decomp_source(x, n_draws, seed, exposure, "x")
  sy <- if (!is.null(y)) decomp_source(y, n_draws, seed + 1, exposure, "y")
  if (!is.null(sy) && !identical(sy$age$labels, sx$age$labels)) {
    mt_abort("`x` and `y` must have the same age grid.", class = "argument")
  }
  if (!is.null(sy) && pairing == "independent" && dim(sy$rates)[4] > 1) {
    sy$rates <- sy$rates[, , , with_seed(seed + 2, sample.int(dim(sy$rates)[4])), drop = FALSE]
  }
  q <- decomp_quantity(qname, quantity, sx$age, at, to, standard, per, ...)
  hierarchy <- hierarchy %||% object_hierarchy(x) %||% if (!is.null(y)) object_hierarchy(y)

  px <- decomp_populations(sx, var, cause, groups, hierarchy, qname)
  py <- if (!is.null(sy)) decomp_populations(sy, var, cause, groups, hierarchy, qname)
  if (!is.null(py) && (!identical(px$cause_levels, py$cause_levels) ||
                       !identical(px$group_levels, py$group_levels))) {
    mt_abort("`x` and `y` must have the same causes and groups.", class = "argument")
  }
  comps <- decomp_comparisons(sx, sy, px, py, var, ref, periods)

  cache <- list()
  get_period <- function(src, tag, p, how) {
    key <- paste(tag, p, how)
    if (is.null(cache[[key]])) cache[[key]] <<- decomp_period(src, p, how, seed)
    cache[[key]]
  }
  side <- function(a, idx) {
    nc <- nrow(idx)
    ng <- ncol(idx)
    R <- a$R[, as.vector(idx), , drop = FALSE]
    list(R = array(R, c(dim(R)[1], nc, ng, dim(R)[3])), E = a$E[, idx[1, ], drop = FALSE], observed = a$observed,
         forecast = a$forecast)
  }
  tables <- list()
  paths <- list()
  totals <- list()
  frozen <- FALSE
  for (k in seq_along(comps$list)) {
    cmp <- comps$list[[k]]
    a1 <- get_period(cmp$src1, cmp$tag1, cmp$p1, observed)
    a2 <- get_period(cmp$src2, cmp$tag2, cmp$p2, observed)
    s1 <- side(a1, cmp$idx1)
    s2 <- side(a2, cmp$idx2)
    n <- max(dim(s1$R)[4], dim(s2$R)[4])
    s1$R <- recycle_draws(s1$R, n)
    s2$R <- recycle_draws(s2$R, n)
    res <- decomp_core(s1, s2, q, method, composition, grouped, n_steps)
    # share of the difference due to observed instead of fitted rates
    jump <- NA_real_
    if (observed == "data" && (isTRUE(s1$observed) || isTRUE(s2$observed))) {
      jump <- 0
      for (sd in 1:2) {
        a <- if (sd == 1) a1 else a2
        if (!isTRUE(a$observed)) next
        src <- if (sd == 1) cmp$src1 else cmp$src2
        if (is.null(src$fitted)) {
          jump <- NA_real_
          break
        }
        f <- side(get_period(src, if (sd == 1) cmp$tag1 else cmp$tag2, if (sd == 1) cmp$p1 else cmp$p2, "fitted"),
                  if (sd == 1) cmp$idx1 else cmp$idx2)
        o <- if (sd == 1) s1 else s2
        o$R <- o$R[, , , 1, drop = FALSE]
        step <- decomp_value(f, q, method, composition, grouped) - decomp_value(o, q, method, composition, grouped)
        jump <- jump + if (sd == 1) step else -step
      }
    }
    if (grouped && is.null(exposure) && (isTRUE(s1$forecast) || isTRUE(s2$forecast))) frozen <- TRUE
    tab <- data.table(comparison = k, component = res$rows$component,
                      group = if (!grouped && is.null(groups)) NA_character_ else
                        px$group_levels[res$rows$group_i],
                      age = sx$age$labels[q$rows][res$rows$age_i],
                      cause = if (is.null(cause)) NA_character_ else px$cause_levels[res$rows$cause_i])
    tables[[k]] <- tab
    paths[[k]] <- res$paths
    totals[[k]] <- list(v1 = res$v1, v2 = res$v2, jump = jump)
  }
  n_max <- max(vapply(paths, ncol, integer(1)))
  P <- do.call(rbind, lapply(paths, recycle_cols, n_max))
  from <- do.call(rbind, lapply(totals, function(t) recycle_cols(matrix(t$v1, 1), n_max)))
  to <- do.call(rbind, lapply(totals, function(t) recycle_cols(matrix(t$v2, 1), n_max)))
  table <- rbindlist(tables)
  # draws without a valid life table are left out
  bad <- apply(!is.finite(P), 2, any) | apply(!is.finite(from) | !is.finite(to), 2, any)
  if (all(bad)) mt_abort("No draw gives a valid decomposition; check the rates for implausible values.",
                         class = "decompose")
  if (any(bad)) {
    mt_warn(sprintf(paste("%d of %d draws left out: their rates give no valid life table (rates far above",
                          "plausible values; see the warning of mt_forecast())."), sum(bad), length(bad)),
            class = "decompose")
    P <- P[, !bad, drop = FALSE]
    from <- from[, !bad, drop = FALSE]
    to <- to[, !bad, drop = FALSE]
  }
  dif <- to - from
  resid <- apply(abs(dif - rowsum(P, table$comparison, reorder = TRUE)), 1, max) / pmax(1, abs(rowMeans(dif)))
  tol <- if (method == "horiuchi") 1e-3 else 1e-6
  if (any(resid > tol)) {
    mt_warn(sprintf(paste("The contributions differ from the difference by up to %.3g (relative) in comparison",
                          "%d.%s"), max(resid), which.max(resid),
                    if (method == "horiuchi") " Increase `n_steps`." else
                      " Check the draws for implausible rates (see mt_forecast())."), class = "decompose")
  }
  structure(list(table = table, paths = P, meta = comps$meta, from = from, to = to,
                 jump_off = vapply(totals, function(t) t$jump, numeric(1)),
                 quantity = qname, method = method, composition = composition, frozen = frozen,
                 at = q$at, probs = probs, age_labels = sx$age$labels[q$rows],
                 cause_levels = px$cause_levels, group_levels = px$group_levels),
            class = "mt_decomposition")
}

# Hierarchy stored in data or draws: from mt_aggregate_strata() or
# reconciled() (`$hierarchy`) and from mt_group_levels(total = )
# (`$level_totals`)
object_hierarchy <- function(obj) {
  vars <- if (!is.null(obj$hierarchy)) as_hierarchy(obj$hierarchy)$vars else list()
  for (v in names(obj$level_totals)) if (is.null(vars[[v]])) vars[[v]] <- obj$level_totals[[v]]
  if (!length(vars)) return(NULL)
  add <- vapply(vars, `[[`, "", "additive")
  structure(list(vars = vars, additive = if (all(add == "rates")) "rates" else "deaths"), class = "mt_hierarchy")
}

# Rates and exposures of a source: observed data (Poisson draws) or draws
decomp_source <- function(obj, n_draws, seed, exposure, name) {
  if (inherits(obj, "mt_data")) {
    d <- mt_observed(obj, n_draws = n_draws, seed = seed)
    return(list(age = d$age, strata = d$strata, periods = d$periods, rates = d$rates, exposure = d$exposure,
                obs_periods = NULL, train_D = d$train_D, train_E = d$train_E, label = "observed", data = TRUE))
  }
  check_class(obj, "mt_draws", name)
  E <- if (!is.null(exposure)) {
    expand_exposure(exposure, dim(obj$rates), dimnames(obj$rates))
  } else if (!is.null(obj$exposure) && identical(dim(obj$exposure), dim(obj$rates)[1:3])) {
    obj$exposure
  } else {
    expand_exposure(obj$exposure_last, dim(obj$rates), dimnames(obj$rates))
  }
  list(age = obj$age, strata = obj$strata, periods = obj$periods, rates = obj$rates, exposure = E,
       obs_periods = obj$observed_periods, observed = obj$observed, obs_exposure = obj$train_E,
       fitted = obj$fitted, train_D = obj$train_D, train_E = obj$train_E, label = obj$model,
       data = identical(obj$model, "Observed rates"))
}

# Rates [age, stratum, draw] and exposures [age, stratum] of one period; for
# observed years of draws the observed rates, the fitted rates or Poisson
# draws of the observed deaths
decomp_period <- function(src, p, how = "data", seed = 1) {
  d <- dim(src$rates)
  k <- match(p, src$periods)
  if (!is.na(k)) {
    return(list(R = array(src$rates[, k, , ], c(d[1], d[3], d[4])), E = matrix(src$exposure[, k, ], d[1]),
                observed = isTRUE(src$data), forecast = !isTRUE(src$data)))
  }
  k <- match(p, src$obs_periods)
  if (!is.na(k) && !is.null(src$observed)) {
    E <- matrix(src$obs_exposure[, k, ], d[1])
    R <- switch(how,
                data = array(src$observed[, k, ], c(d[1], d[3], 1)),
                fitted = {
                  if (is.null(src$fitted)) {
                    mt_abort(sprintf("`observed = \"fitted\"`: the draws of %s hold no fitted rates.", src$label),
                             class = "argument")
                  }
                  array(src$fitted[, k, ], c(d[1], d[3], 1))
                },
                poisson = {
                  D <- src$train_D[, k, ]
                  dd <- with_seed(seed, stats::rpois(length(D) * d[4], rep(as.vector(D), d[4])))
                  array(dd, c(d[1], d[3], d[4])) / as.vector(pmax(E, 1e-12))
                })
    return(list(R = R, E = E, observed = TRUE, forecast = FALSE))
  }
  mt_abort(sprintf("Year %s is not in the data or draws.", p), class = "argument")
}

recycle_draws <- function(R, n) {
  if (dim(R)[4] == n) return(R)
  if (dim(R)[4] != 1) mt_abort("The two sides have different numbers of draws.", class = "argument")
  R[, , , rep(1L, n), drop = FALSE]
}

recycle_cols <- function(M, n) if (ncol(M) == n) M else M[, rep(1L, n), drop = FALSE]

# Quantity: function of rates [ages used, n] and the ages it uses
decomp_quantity <- function(qname, quantity, age, at, to, standard, per, ...) {
  na <- length(age)
  i0 <- if (is.null(at)) 1L else match(at, age$lower)
  if (is.na(i0)) mt_abort("`at` must be the start of an age interval.", class = "argument")
  sub_age <- function(rows) list(lower = age$lower[rows], upper = age$upper[rows])
  lt_args <- list(...)
  switch(qname,
         life_expectancy = {
           rows <- seq(i0, na)
           ag <- sub_age(rows)
           if (!is.null(lt_args$ax) && length(lt_args$ax) == na) lt_args$ax <- lt_args$ax[rows]
           if (!identical(lt_args$open %||% "reciprocal", "reciprocal")) ag <- age_subset(age, rows)
           lt <- function(M) do.call(lt_columns, c(list(M, ag), lt_args))
           list(rows = rows, at = age$lower[i0], fun = function(M) lt(M)$ex[1, ], lt = lt)
         },
         prob_dying = {
           j <- if (is.null(to)) NA else match(to, age$lower)
           if (is.na(j) || j <= i0) {
             mt_abort("`to` must be the start of an age interval above `at` (below the open interval).",
                      class = "argument")
           }
           rows <- seq(i0, j - 1)
           nn <- age$upper[rows] - age$lower[rows]
           list(rows = rows, at = age$lower[i0], fun = function(M) {
             Q <- lt_qx(M, matrix(nn, nrow(M), ncol(M)), matrix(nn / 2, nrow(M), ncol(M)))
             1 - apply(1 - Q, 2, prod)
           })
         },
         asr = {
           if (is.character(standard)) standard <- standard_population(standard)
           w <- standard_weights(age, standard)
           rows <- which(w > 0)
           ws <- w[rows] / sum(w) * per
           list(rows = rows, at = NA, fun = function(M) colSums(M * ws))
         },
         crude_rate = list(rows = seq_len(na), at = NA, per = per),
         custom = {
           rows <- seq_len(na)
           list(rows = rows, at = NA, fun = function(M) {
             v <- quantity(M, age)
             if (length(v) != ncol(M)) mt_abort("`quantity` must return one value per column.", class = "argument")
             v
           })
         })
}

# An mt_axis-like list of selected intervals (for the Gompertz and Kannisto
# closures of the life table, which need the axis)
age_subset <- function(age, rows) {
  out <- age
  for (nm in intersect(names(age), c("lower", "upper", "labels"))) out[[nm]] <- age[[nm]][rows]
  out
}

# Populations: strata grouped by everything except causes and groups
decomp_populations <- function(src, var, cause, groups, hierarchy, qname = "life_expectancy") {
  st <- src$strata
  vars <- setdiff(names(st), c("stratum", "label"))
  given <- c(var, cause, groups)
  bad <- setdiff(given, vars)
  if (length(bad)) {
    mt_abort(sprintf("Not a stratum variable: %s (stratum variables: %s).", paste(bad, collapse = ", "),
                     paste(vars, collapse = ", ")), class = "argument")
  }
  if (anyDuplicated(given)) mt_abort("`var`, `cause` and `groups` must differ.", class = "argument")
  h <- if (!is.null(hierarchy)) as_hierarchy(hierarchy)
  aggregate_levels <- function(v) {
    hv <- if (!is.null(h)) h$vars[[v]]
    if (is.null(hv)) character() else c(hv$total, names(hv$groups))
  }
  keep <- rep(TRUE, nrow(st))
  for (v in c(cause, groups)) keep <- keep & !st[[v]] %in% aggregate_levels(v)
  # totals that are not declared would be counted twice
  if (!is.null(src$train_D) && length(c(cause, groups))) {
    tl <- total_levels_array(src$train_D, st[keep], c(cause, groups))
    if (nrow(tl)) {
      mt_abort(sprintf(paste("Level \"%s\" of `%s` equals the sum of the other levels and would be counted as a",
                             "part. Declare it, for example hierarchy = mt_hierarchy(%s = \"%s\"%s), or remove",
                             "it with mt_subset()."), tl$level[1], tl$var[1], tl$var[1], tl$level[1],
                       if (identical(tl$var[1], cause)) ", additive = \"rates\"" else ""), class = "hierarchy")
    }
  }
  # life expectancy of single causes is not meaningful
  if (qname %in% c("life_expectancy", "prob_dying") && !is.null(src$train_D)) {
    for (v in setdiff(vars, c(var, cause, groups))) {
      if (nrow(total_levels_array(src$train_D, st, v)) && shares_exposure(src$train_E, st, v)) {
        mt_warn(sprintf(paste("The levels of `%s` share one exposure and one of them is their sum, as causes of",
                              "death do: life expectancy of a single level is not meaningful. Give cause = \"%s\"",
                              "to decompose by cause."), v, v), class = "decompose")
      }
    }
  }
  s <- st[keep]
  cl <- if (!is.null(cause)) unique(s[[cause]]) else NA_character_
  gl <- if (!is.null(groups)) unique(s[[groups]]) else NA_character_
  pop_vars <- setdiff(vars, c(cause, groups))
  key <- if (length(pop_vars)) do.call(paste, c(lapply(pop_vars, function(v) s[[v]]), sep = "\r")) else
    rep("", nrow(s))
  pops <- unique(cbind(s[, pop_vars, with = FALSE], pop_key = key))
  idx <- lapply(pops$pop_key, function(k) {
    r <- which(key == k)
    m <- matrix(NA_integer_, length(cl), length(gl))
    ci <- if (is.null(cause)) rep(1L, length(r)) else match(s[[cause]][r], cl)
    gi <- if (is.null(groups)) rep(1L, length(r)) else match(s[[groups]][r], gl)
    m[cbind(ci, gi)] <- s$stratum[r]
    if (anyNA(m)) {
      mt_abort(sprintf("Stratum %s lacks some causes or groups.", gsub("\r", ".", k)), class = "argument")
    }
    m
  })
  list(pops = pops, idx = idx, pop_vars = pop_vars, cause_levels = if (!is.null(cause)) cl,
       group_levels = if (!is.null(groups)) gl)
}

# Levels of `vars` whose deaths equal the sum of the other levels (as in
# mt_validate()); D [age, period, stratum] with strata `st` (rows of
# st$stratum index the third dimension)
total_levels_array <- function(D, st, vars) {
  out <- list()
  for (v in vars) {
    levs <- unique(st[[v]])
    if (length(levs) < 3) next
    others <- setdiff(names(st), c("stratum", "label", v))
    key <- if (length(others)) do.call(paste, c(lapply(others, function(o) st[[o]]), sep = "\r")) else
      rep("", nrow(st))
    for (l in levs) {
      hit <- 0
      n <- 0
      for (k in unique(key[st[[v]] == l])) {
        i <- st$stratum[key == k & st[[v]] == l]
        j <- st$stratum[key == k & st[[v]] != l]
        if (length(i) != 1 || length(j) != length(levs) - 1) next
        a <- D[, , i]
        b <- apply(D[, , j, drop = FALSE], c(1, 2), sum)
        ok <- is.finite(a) & is.finite(b) & a > 0
        hit <- hit + sum(abs(a[ok] - b[ok]) <= pmax(0.005 * a[ok], 0.5))
        n <- n + sum(ok)
      }
      if (n >= 5 && hit / n >= 0.9) out[[length(out) + 1L]] <- data.table(var = v, level = l)
    }
  }
  if (length(out)) rbindlist(out) else data.table(var = character(), level = character())
}

# Do the levels of `v` share one exposure (as causes of death do)?
shares_exposure <- function(E, st, v) {
  if (is.null(E)) return(FALSE)
  others <- setdiff(names(st), c("stratum", "label", v))
  key <- if (length(others)) do.call(paste, c(lapply(others, function(o) st[[o]]), sep = "\r")) else
    rep("", nrow(st))
  all(vapply(unique(key), function(k) {
    i <- st$stratum[key == k]
    if (length(i) < 2) return(TRUE)
    ref <- E[, , i[1]]
    all(vapply(i[-1], function(j) isTRUE(all.equal(E[, , j], ref, tolerance = 1e-8)), logical(1)))
  }, logical(1)))
}

decomp_comparisons <- function(sx, sy, px, py, var, ref, periods) {
  last_any <- function(s) max(c(s$periods, s$obs_periods))
  acc <- new.env(parent = emptyenv())
  acc$out <- list()
  acc$meta <- list()
  add <- function(m, src1, tag1, p1, idx1, src2, tag2, p2, idx2) {
    acc$out[[length(acc$out) + 1L]] <- list(src1 = src1, tag1 = tag1, p1 = p1, idx1 = idx1, src2 = src2,
                                            tag2 = tag2, p2 = p2, idx2 = idx2)
    acc$meta[[length(acc$meta) + 1L]] <- m
  }
  pv <- px$pop_vars
  if (!is.null(var)) {
    if (is.null(ref) || !ref %in% px$pops[[var]]) {
      mt_abort(sprintf("`ref` must be a level of `%s`.", var), class = "argument")
    }
    periods <- periods %||% max(sx$periods)
    others <- setdiff(pv, var)
    okey <- if (length(others)) do.call(paste, c(lapply(others, function(o) px$pops[[o]]), sep = "\r")) else
      rep("", nrow(px$pops))
    for (i in which(px$pops[[var]] != ref)) {
      j <- which(px$pops[[var]] == ref & okey == okey[i])
      if (length(j) != 1) next
      for (p in periods) {
        m <- data.table(from = ref, to = px$pops[[var]][i], period = p)
        if (length(others)) m <- cbind(px$pops[i, others, with = FALSE], m)
        add(m, sx, "x", p, px$idx[[j]], sx, "x", p, px$idx[[i]])
      }
    }
  } else if (!is.null(sy)) {
    periods <- periods %||% min(last_any(sx), last_any(sy))
    for (i in seq_len(nrow(px$pops))) {
      j <- match(px$pops$pop_key[i], py$pops$pop_key)
      if (is.na(j)) next
      for (p in periods) {
        m <- data.table(from = sx$label, to = sy$label, period = p)
        if (length(pv)) m <- cbind(px$pops[i, pv, with = FALSE], m)
        add(m, sx, "x", p, px$idx[[i]], sy, "y", p, py$idx[[j]])
      }
    }
  } else {
    periods <- periods %||% if (length(sx$obs_periods)) c(max(sx$obs_periods), max(sx$periods)) else
      range(sx$periods)
    if (length(periods) < 2) mt_abort("Give at least two `periods` to decompose a change.", class = "argument")
    for (i in seq_len(nrow(px$pops))) {
      for (p in periods[-1]) {
        m <- data.table(from = as.character(periods[1]), to = as.character(p))
        if (length(pv)) m <- cbind(px$pops[i, pv, with = FALSE], m)
        add(m, sx, "x", periods[1], px$idx[[i]], sx, "x", p, px$idx[[i]])
      }
    }
  }
  if (!length(acc$out)) mt_abort("No comparisons: check `var`, `ref` and the strata of `y`.", class = "argument")
  meta <- rbindlist(acc$meta)
  meta[, comparison := seq_len(.N)]
  list(list = acc$out, meta = meta)
}

# One comparison: s1, s2 with R [age, cause, group, draw] and E [age, group]
decomp_core <- function(s1, s2, q, method, composition, grouped, n_steps) {
  rows <- q$rows
  R1 <- s1$R[rows, , , , drop = FALSE]
  R2 <- s2$R[rows, , , , drop = FALSE]
  d <- dim(R1)
  nr <- d[1]
  nc <- d[2]
  ng <- d[3]
  n <- d[4]
  allc <- function(R) {
    M <- R[, 1, , , drop = FALSE]
    if (nc > 1) for (c in 2:nc) M <- M + R[, c, , , drop = FALSE]
    array(M, c(nr, ng, n))
  }
  M1 <- allc(R1)
  M2 <- allc(R2)
  label_rows <- function(component, age_i, cause_i, group_i) {
    data.table(component = component, age_i = age_i, cause_i = cause_i, group_i = group_i)
  }
  v1 <- decomp_value(s1, q, method, composition, grouped)
  v2 <- decomp_value(s2, q, method, composition, grouped)
  if (method == "kitagawa") {
    c1 <- s1$E / sum(s1$E)
    c2 <- s2$E / sum(s2$E)
    cb <- (c1 + c2) / 2
    dc <- c2 - c1
    per <- q$per
    rate <- array(0, c(nr, nc, ng, n))
    compo <- array(0, c(nr, nc, n))
    for (g in seq_len(ng)) for (c in seq_len(nc)) {
      rate[, c, g, ] <- (R2[, c, g, ] - R1[, c, g, ]) * cb[, g] * per
      compo[, c, ] <- compo[, c, ] + (R1[, c, g, ] + R2[, c, g, ]) / 2 * dc[, g] * per
    }
    tab <- rbind(label_rows("rate", rep(seq_len(nr), nc * ng), rep(rep(seq_len(nc), each = nr), ng),
                            rep(seq_len(ng), each = nr * nc)),
                 label_rows("composition", rep(seq_len(nr), nc), rep(seq_len(nc), each = nr), NA_integer_))
    return(list(rows = tab, paths = rbind(matrix(rate, ncol = n), matrix(compo, ncol = n)), v1 = v1, v2 = v2))
  }
  # decomposition of one homogeneous population: [age, cause, draw]
  by_age <- function(Ra, Rb, Ma, Mb) {
    if (method == "arriaga") {
      dl <- (arriaga_dir(Ma, Mb, q$lt) - arriaga_dir(Mb, Ma, q$lt)) / 2
      if (nc == 1) return(array(dl, c(nr, 1, n)))
      return(cause_split(dl, Ra, Rb, Ma, Mb, q$fun))
    }
    f <- function(P) {
      M <- P[seq_len(nr), , drop = FALSE]
      if (nc > 1) for (c in 2:nc) M <- M + P[(c - 1) * nr + seq_len(nr), , drop = FALSE]
      q$fun(M)
    }
    P1 <- matrix(Ra, nr * nc)
    P2 <- matrix(Rb, nr * nc)
    ord <- order(rep(seq_len(nr), nc))
    out <- if (method == "horiuchi") horiuchi_core(f, P1, P2, n_steps) else stepwise_core(f, P1, P2, ord)
    array(out, c(nr, nc, n))
  }
  if (!grouped) {
    cont <- by_age(array(R1[, , 1, ], c(nr, nc, n)), array(R2[, , 1, ], c(nr, nc, n)),
                   matrix(M1[, 1, ], nr), matrix(M2[, 1, ], nr))
    tab <- label_rows("rate", rep(seq_len(nr), nc), rep(seq_len(nc), each = nr), NA_integer_)
    return(list(rows = tab, paths = matrix(cont, ncol = n), v1 = v1, v2 = v2))
  }
  if (composition == "radix") {
    w1 <- s1$E[rows[1], ] / sum(s1$E[rows[1], ])
    w2 <- s2$E[rows[1], ] / sum(s2$E[rows[1], ])
    paths <- list()
    q1 <- q2 <- matrix(0, ng, n)
    for (g in seq_len(ng)) {
      Ma <- matrix(M1[, g, ], nr)
      Mb <- matrix(M2[, g, ], nr)
      q1[g, ] <- q$fun(Ma)
      q2[g, ] <- q$fun(Mb)
      cont <- by_age(array(R1[, , g, ], c(nr, nc, n)), array(R2[, , g, ], c(nr, nc, n)), Ma, Mb)
      paths[[g]] <- matrix(cont, ncol = n) * (w1[g] + w2[g]) / 2
    }
    kt <- kitagawa_mean(w1, w2, q1, q2)
    tab <- rbind(label_rows("rate", rep(seq_len(nr), nc * ng), rep(rep(seq_len(nc), each = nr), ng),
                            rep(seq_len(ng), each = nr * nc)),
                 label_rows("composition", NA_integer_, NA_integer_, NA_integer_))
    return(list(rows = tab, paths = rbind(do.call(rbind, paths), matrix(kt$composition, 1)), v1 = v1, v2 = v2))
  }
  # pooled: rates of every group and the shares by age are the parameters
  sh1 <- s1$E[rows, , drop = FALSE] / rowSums(s1$E[rows, , drop = FALSE])
  sh2 <- s2$E[rows, , drop = FALSE] / rowSums(s2$E[rows, , drop = FALSE])
  nrate <- nr * nc * ng
  f <- function(P) {
    k <- ncol(P)
    M <- matrix(0, nr, k)
    for (g in seq_len(ng)) {
      Mg <- matrix(0, nr, k)
      for (c in seq_len(nc)) Mg <- Mg + P[((g - 1) * nc + c - 1) * nr + seq_len(nr), , drop = FALSE]
      M <- M + Mg * P[nrate + (g - 1) * nr + seq_len(nr), , drop = FALSE]
    }
    q$fun(M)
  }
  P1 <- rbind(matrix(R1, nrate), matrix(sh1, nr * ng, n))
  P2 <- rbind(matrix(R2, nrate), matrix(sh2, nr * ng, n))
  age_of <- c(rep(seq_len(nr), nc * ng), rep(seq_len(nr), ng))
  out <- if (method == "horiuchi") horiuchi_core(f, P1, P2, n_steps) else stepwise_core(f, P1, P2, order(age_of))
  compo <- matrix(0, nr, n)
  for (g in seq_len(ng)) compo <- compo + out[nrate + (g - 1) * nr + seq_len(nr), , drop = FALSE]
  tab <- rbind(label_rows("rate", rep(seq_len(nr), nc * ng), rep(rep(seq_len(nc), each = nr), ng),
                          rep(seq_len(ng), each = nr * nc)),
               label_rows("composition", seq_len(nr), NA_integer_, NA_integer_))
  list(rows = tab, paths = rbind(out[seq_len(nrate), , drop = FALSE], compo), v1 = v1, v2 = v2)
}

# Value of the quantity for one side (R [age, cause, group, draw], E [age,
# group]): the population of the comparison, with subgroups blended at the
# radix or pooled
decomp_value <- function(s, q, method, composition, grouped) {
  rows <- q$rows
  R <- s$R[rows, , , , drop = FALSE]
  d <- dim(R)
  M <- R[, 1, , , drop = FALSE]
  if (d[2] > 1) for (c in 2:d[2]) M <- M + R[, c, , , drop = FALSE]
  M <- array(M, c(d[1], d[3], d[4]))
  if (method == "kitagawa") {
    cc <- s$E / sum(s$E)
    return(colSums(matrix(M * as.vector(cc), d[1] * d[3])) * q$per)
  }
  if (!grouped) return(q$fun(matrix(M[, 1, ], d[1])))
  if (identical(composition, "radix")) {
    w <- s$E[rows[1], ] / sum(s$E[rows[1], ])
    return(colSums(do.call(rbind, lapply(seq_len(d[3]), function(g) q$fun(matrix(M[, g, ], d[1])))) * w))
  }
  sh <- s$E[rows, , drop = FALSE] / rowSums(s$E[rows, , drop = FALSE])
  Mp <- matrix(0, d[1], d[4])
  for (g in seq_len(d[3])) Mp <- Mp + matrix(M[, g, ], d[1]) * sh[, g]
  q$fun(Mp)
}

# Kitagawa decomposition of a difference of weighted means sum(w * v):
# composition sum(dw * mean(v)) and value effects dv * mean(w) per element
kitagawa_mean <- function(w1, w2, v1, v2) {
  v1 <- as.matrix(v1)
  v2 <- as.matrix(v2)
  list(composition = colSums((w2 - w1) * (v1 + v2) / 2), value = (v2 - v1) * (w1 + w2) / 2)
}

# Arriaga's decomposition from M1 to M2 (columns are schedules), ages from the
# start of the life table; sums to e(M2) - e(M1) at the first age
arriaga_dir <- function(M1, M2, lt) {
  a <- lt(M1)
  b <- lt(M2)
  na <- nrow(M1)
  out <- matrix(0, na, ncol(M1))
  if (na > 1) {
    i <- seq_len(na - 1)
    out[i, ] <- a$lx[i, ] * (b$Lx[i, ] / b$lx[i, ] - a$Lx[i, ] / a$lx[i, ]) +
      b$Tx[i + 1, ] * (a$lx[i, ] / b$lx[i, ] - a$lx[i + 1, ] / b$lx[i + 1, ])
  }
  out[na, ] <- a$lx[na, ] * (b$Tx[na, ] / b$lx[na, ] - a$Tx[na, ] / a$lx[na, ])
  out
}

# Age contributions split by cause with the implied sensitivity
# dl / (Mb - Ma); a numerical derivative where the all-cause rates hardly
# differ (contributions that are not finite stay so and the draw is dropped)
cause_split <- function(dl, Ra, Rb, Ma, Mb, fun) {
  dM <- Mb - Ma
  s <- dl / dM
  bad <- is.finite(dl) & (abs(dM) <= 1e-8 * (Ma + Mb) | !is.finite(s))
  if (any(bad)) {
    Mbar <- (Ma + Mb) / 2
    for (i in which(rowSums(bad) > 0)) {
      cols <- which(bad[i, ])
      h <- pmax(1e-6 * Mbar[i, cols], 1e-12)
      Mp <- Mm <- Mbar[, cols, drop = FALSE]
      Mp[i, ] <- Mp[i, ] + h
      Mm[i, ] <- Mm[i, ] - h
      s[i, cols] <- (fun(Mp) - fun(Mm)) / (2 * h)
    }
  }
  out <- Rb - Ra
  for (c in seq_len(dim(out)[2])) out[, c, ] <- out[, c, ] * s
  out
}

# Horiuchi, Wilmoth and Pletcher (2008): f maps parameters [np, k] to k
# values; P1, P2 are [np, n]; the derivative along the straight path is
# taken by central differences at the midpoints of n_steps intervals
horiuchi_core <- function(f, P1, P2, n_steps) {
  np <- nrow(P1)
  n <- ncol(P1)
  D <- (P2 - P1) / n_steps
  out <- matrix(0, np, n)
  chunk <- max(1L, floor(4e6 / (2 * np * n)))
  blocks <- split(seq_len(np), ceiling(seq_len(np) / chunk))
  for (k in seq_len(n_steps)) {
    Pm <- P1 + D * (k - 0.5)
    for (ch in blocks) {
      m <- length(ch)
      big <- Pm[, rep(seq_len(n), 2 * m), drop = FALSE]
      for (jj in seq_len(m)) {
        j <- ch[jj]
        cp <- (2 * jj - 2) * n + seq_len(n)
        big[j, cp] <- Pm[j, ] + D[j, ] / 2
        big[j, cp + n] <- Pm[j, ] - D[j, ] / 2
      }
      v <- f(big)
      for (jj in seq_len(m)) {
        cp <- (2 * jj - 2) * n + seq_len(n)
        out[ch[jj], ] <- out[ch[jj], ] + v[cp] - v[cp + n]
      }
    }
  }
  out
}

# Stepwise replacement (Andreev, Shkolnikov and Begun 2002) in the order
# `ord`, averaged over both directions
stepwise_core <- function(f, P1, P2, ord, symmetric = TRUE) {
  one <- function(A, B) {
    out <- matrix(0, nrow(A), ncol(A))
    cur <- A
    v0 <- f(cur)
    for (j in ord) {
      cur[j, ] <- B[j, ]
      v1 <- f(cur)
      out[j, ] <- v1 - v0
      v0 <- v1
    }
    out
  }
  fw <- one(P1, P2)
  if (!symmetric) return(fw)
  (fw - one(P2, P1)) / 2
}

#' @rdname mt_decompose
#' @param object An `mt_decomposition` object.
#' @param by Variables of the contributions to keep (`"component"`,
#'   `"group"`, `"age"`, `"cause"`); the others are summed path by path.
#'   `NULL` gives the totals.
#' @param change Two periods: the change of each comparison between them
#'   (for comparisons between strata or objects in several periods), formed
#'   path by path.
#' @export
summary.mt_decomposition <- function(object, by = c("component", "group", "age", "cause"), probs = object$probs,
                                     change = NULL, ...) {
  if (!missing(by) && !is.null(by)) {
    absent <- setdiff(by, c("component", "group", "age", "cause"))
    absent <- c(absent, intersect(by, c("group", "age", "cause"))[vapply(intersect(by, c("group", "age", "cause")),
                                                                         function(b) all(is.na(object$table[[b]])),
                                                                         logical(1))])
    if (length(absent)) {
      mt_abort(sprintf("The decomposition has no `%s`.", paste(absent, collapse = "`, `")), class = "argument")
    }
  }
  if (!is.null(change)) object <- decomp_change(object, change)
  if (is.null(by)) return(decomp_totals(object, probs)[])
  out <- decomp_summary(object, by, probs)
  out[, comparison := NULL]
  out[]
}

# Change between two periods of the same comparisons, path by path, as a
# decomposition object
decomp_change <- function(object, change) {
  meta <- object$meta
  if (!"period" %in% names(meta) || length(change) != 2) {
    mt_abort("`change` needs two periods of comparisons between strata or objects.", class = "argument")
  }
  keys <- setdiff(names(meta), c("comparison", "period"))
  kk <- do.call(paste, c(lapply(keys, function(k) meta[[k]]), sep = "\r"))
  c1 <- meta$comparison[meta$period == change[1]]
  c2 <- meta$comparison[meta$period == change[2]]
  if (!length(c1) || !length(c2)) mt_abort("Both `change` periods must be in the decomposition.", class = "argument")
  pairs <- data.table(c2 = c2, c1 = c1[match(kk[c2], kk[c1])])[!is.na(c1)]
  tab <- object$table
  paths <- list()
  tables <- list()
  for (i in seq_len(nrow(pairs))) {
    r1 <- which(tab$comparison == pairs$c1[i])
    r2 <- which(tab$comparison == pairs$c2[i])
    paths[[i]] <- object$paths[r2, , drop = FALSE] - object$paths[r1, , drop = FALSE]
    tables[[i]] <- copy(tab[r2])[, comparison := i]
  }
  m <- meta[pairs$c2]
  m$period <- rep(paste(change[1], "to", change[2]), nrow(m))
  m$comparison <- seq_len(nrow(m))
  d1 <- object$to[pairs$c1, , drop = FALSE] - object$from[pairs$c1, , drop = FALSE]
  d2 <- object$to[pairs$c2, , drop = FALSE] - object$from[pairs$c2, , drop = FALSE]
  object$table <- rbindlist(tables)
  object$paths <- do.call(rbind, paths)
  object$meta <- m
  object$from <- d1
  object$to <- d2
  object$jump_off <- rep(NA_real_, nrow(m))
  object
}

# Contributions summed path by path over the variables not in `by`
decomp_summary <- function(object, by, probs) {
  tab <- object$table
  by <- intersect(c("component", "group", "age", "cause"), by)
  by <- by[vapply(by, function(b) !all(is.na(tab[[b]])), logical(1))]
  keys <- tab[, c("comparison", by), with = FALSE]
  grp <- do.call(paste, c(as.list(keys), sep = "\r"))
  grp <- factor(grp, levels = unique(grp))
  S <- rowsum(object$paths, grp, reorder = FALSE)
  out <- keys[!duplicated(grp)]
  Q <- row_quantiles(S, probs)
  out[, mean := rowMeans(S)]
  for (k in seq_along(probs)) set(out, j = qnames(probs)[k], value = Q[, k])
  merge(object$meta, out, by = "comparison", sort = FALSE)
}

decomp_totals <- function(object, probs = object$probs) {
  tot <- rowsum(object$paths, object$table$comparison, reorder = TRUE)
  dif <- object$to - object$from
  Q <- row_quantiles(dif, probs)
  out <- cbind(object$meta, from_value = rowMeans(object$from), to_value = rowMeans(object$to),
               difference = rowMeans(dif))
  for (k in seq_along(probs)) set(out, j = paste0("difference_", qnames(probs)[k]), value = Q[, k])
  if (any(!is.na(object$jump_off))) set(out, j = "jump_off", value = object$jump_off)
  out[, `:=`(sum_contributions = rowMeans(tot), max_residual = apply(abs(dif - tot), 1, max))]
  setcolorder(out, "comparison")
  out
}

#' @export
print.mt_decomposition <- function(x, ...) {
  what <- if (x$quantity == "life_expectancy") sprintf("life expectancy at %s", x$at) else x$quantity
  cat(sprintf("<mt_decomposition> %s, method %s%s; %d comparison(s), %s\n", what, x$method,
              if (!is.null(x$composition)) sprintf(", %s composition", x$composition) else "",
              nrow(x$meta), if (ncol(x$paths) > 1) sprintf("means over %d draws", ncol(x$paths)) else "1 draw"))
  if (identical(x$composition, "radix")) {
    cat(sprintf(paste0("Values: means of the subgroups' values weighted by their shares at age %s (subgroups ",
                       "blended at the radix),\nnot the values of the pooled rates.\n"), x$at))
  } else if (identical(x$composition, "pooled")) {
    cat("Values: of the pooled rates (subgroups weighted by their shares at each age).\n")
  }
  tot <- summary(x, by = NULL)
  num <- setdiff(names(tot)[vapply(tot, is.numeric, logical(1))], c("comparison", "period", "max_residual"))
  if (ncol(x$paths) == 1) tot[, (grep("^difference_q", names(tot), value = TRUE)) := NULL]
  num <- intersect(num, names(tot))
  tot[, (num) := lapply(.SD, round, 3), .SDcols = num]
  tot[, max_residual := signif(max_residual, 2)]
  print(tot)
  if (any(!is.na(x$jump_off))) {
    cat(paste("jump_off: part of the difference due to observed instead of fitted rates in observed years of the",
              "draws\n(see mt_jump_off()); observed = \"fitted\" decomposes the model's change alone.\n"))
  }
  if (isTRUE(x$frozen)) {
    cat(paste("Shares in forecast years are those of the last observed year; give `exposure` (a population",
              "projection)\nfor composition effects in the forecast.\n"))
  }
  invisible(x)
}

#' @export
as.data.table.mt_decomposition <- function(x, keep.rownames = FALSE, ...) summary(x)

#' Plot a decomposition
#'
#' Bars of the mean contributions by age, cause or subgroup, stacked by a
#' second variable, with the 95% interval of the sum per bar when the
#' decomposition has more than one draw. The composition effect of
#' decompositions with subgroups is shown as its own bar.
#'
#' @param decomposition An [mt_decompose()] result.
#' @param by Variable on the axis: `"age"`, `"cause"` or `"group"`.
#' @param fill Variable that stacks the bars (`"cause"`, `"group"`), or
#'   `NULL`.
#' @param comparisons Comparisons to show: numbers (column `comparison` of
#'   `summary(decomposition, by = NULL)`) or a named list of values of its
#'   columns, for example `list(to = "Q5")`; default: all, as panels.
#' @param change Two periods: plot the change between them (see
#'   [summary.mt_decomposition()]).
#' @return A ggplot object.
#' @examples
#' x <- mt_example(deaths = "single", years = 2000:2015)
#' mt_plot_decomposition(mt_decompose(x, periods = c(2000, 2015), at = 30))
#' @export
mt_plot_decomposition <- function(decomposition, by = c("age", "cause", "group"), fill = NULL,
                                  comparisons = NULL, change = NULL) {
  need_ggplot()
  check_class(decomposition, "mt_decomposition", "decomposition")
  by <- match.arg(by)
  if (!is.null(change)) decomposition <- decomp_change(decomposition, change)
  tab <- decomposition$table
  if (all(is.na(tab[[by]]))) mt_abort(sprintf("The decomposition has no `%s`.", by), class = "argument")
  if (!is.null(fill) && (identical(fill, by) || !fill %in% c("cause", "group") || all(is.na(tab[[fill]])))) {
    mt_abort("`fill` must be another variable of the decomposition (\"cause\" or \"group\").", class = "argument")
  }
  meta <- copy(decomposition$meta)
  keep <- select_comparisons(meta, comparisons)
  strata_cols <- setdiff(names(meta), c("comparison", "from", "to", "period"))
  vs <- paste(meta$to, "vs", meta$from)
  if ("period" %in% names(meta)) vs <- paste0(vs, ", ", meta$period)
  meta[, panel := if (length(strata_cols)) {
    paste0(do.call(paste, c(lapply(.SD, as.character), sep = ", ")), ": ", vs)
  } else {
    vs
  }, .SDcols = strata_cols]
  prep <- function(d) {
    d <- d[comparison %in% keep]
    d[, panel := meta$panel[match(comparison, meta$comparison)]]
    d[component == "composition", (by) := "composition"]
    lev <- switch(by, age = decomposition$age_labels, cause = decomposition$cause_levels,
                  group = decomposition$group_levels)
    d[, xvar := factor(get(by), levels = c(lev, if (any(d[[by]] == "composition")) "composition"))]
    d
  }
  s <- prep(decomp_summary(decomposition, c("component", by, fill), c(0.025, 0.975)))
  tot <- prep(decomp_summary(decomposition, c("component", by), c(0.025, 0.975)))
  p <- ggplot2::ggplot(s, ggplot2::aes(x = xvar, y = mean))
  if (!is.null(fill)) {
    s[component == "composition", (fill) := "composition"]
    lv <- unique(s[[fill]])
    s[, fvar := factor(get(fill), levels = c(setdiff(lv, "composition"), intersect(lv, "composition")))]
    p <- p + ggplot2::geom_col(ggplot2::aes(fill = fvar), width = 0.8) + scale_fill_mortensor(type = "categorical")
  } else {
    p <- p + ggplot2::geom_col(fill = mt_palette(1, type = "categorical"), width = 0.8)
  }
  if (ncol(decomposition$paths) > 1) {
    p <- p + ggplot2::geom_errorbar(data = tot, ggplot2::aes(ymin = q0.025, ymax = q0.975), width = 0.3,
                                    colour = "grey30")
  }
  unit <- if (decomposition$quantity == "life_expectancy") " (years)" else ""
  p <- p + ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
    ggplot2::facet_wrap(ggplot2::vars(panel), labeller = ggplot2::label_wrap_gen(width = 40)) +
    ggplot2::labs(x = NULL, y = paste0("Contribution", unit), fill = NULL,
                  subtitle = sprintf("Difference in %s", decomp_label(decomposition)),
                  caption = sprintf("Method: %s%s", decomposition$method,
                                    if (ncol(decomposition$paths) > 1) "; bars: means; lines: 95% intervals"
                                    else "")) +
    theme_mortensor()
  p <- p + if (by == "age") {
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  } else {
    ggplot2::coord_flip()
  }
  p
}

# Comparisons selected by number or by a named list of meta values
select_comparisons <- function(meta, comparisons) {
  if (is.null(comparisons)) return(meta$comparison)
  if (is.list(comparisons)) {
    bad <- setdiff(names(comparisons), names(meta))
    if (is.null(names(comparisons)) || length(bad)) {
      mt_abort(sprintf("`comparisons` must name columns of the totals (%s).",
                       paste(setdiff(names(meta), "comparison"), collapse = ", ")), class = "argument")
    }
    ok <- Reduce(`&`, lapply(names(comparisons), function(v) meta[[v]] %in% comparisons[[v]]))
    keep <- meta$comparison[ok]
  } else {
    keep <- intersect(comparisons, meta$comparison)
  }
  if (!length(keep)) mt_abort("No comparison matches `comparisons`.", class = "argument")
  keep
}

decomp_label <- function(d) {
  switch(d$quantity, life_expectancy = sprintf("life expectancy at %s", d$at),
         prob_dying = "the probability of dying", asr = "the age-standardised rate",
         crude_rate = "the crude death rate", "the quantity")
}
