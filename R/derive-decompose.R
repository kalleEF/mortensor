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
#' [mt_observed()]. Observed years inside forecast draws are the same in
#' every draw.
#'
#' @section Methods:
#' * `"arriaga"` (life expectancy): the symmetric Arriaga decomposition, the
#'   average of the decompositions from 1 to 2 and (with the sign reversed)
#'   from 2 to 1, computed on the package's life table. Causes get the
#'   age-specific contribution in proportion to their rate differences
#'   (an implied sensitivity, Riffe et al. 2024, eqs. 12-13); where the
#'   all-cause difference is close to zero, the sensitivity is a numerical
#'   derivative at the average rates. Exact.
#' * `"horiuchi"` (any quantity): the line-integral method of Horiuchi,
#'   Wilmoth and Pletcher (2008) with `n_steps` steps; exact up to a small
#'   discretisation error, reported in the totals.
#' * `"stepwise"` (any quantity): stepwise replacement (Andreev, Shkolnikov
#'   and Begun 2002) from the youngest to the oldest age, averaged over
#'   both directions of the comparison. Exact; results depend on the order.
#'   Arriaga's decomposition from 1 to 2 is stepwise replacement from the
#'   oldest to the youngest age, so for life expectancy without causes the
#'   two symmetric versions give identical results.
#' * `"kitagawa"` (crude death rate): the rate effect (age-specific rates at
#'   the average age structure) and the composition effect (age structure
#'   at the average rates) of Kitagawa (1955). Exact.
#'
#' @section Subgroups:
#' With `groups`, the populations being compared consist of the levels of a
#' stratum variable (for example deprivation quintiles or education
#' groups), and the decomposition separates the effect of their
#' composition:
#'
#' * `composition = "radix"` (Riffe, Tursun-zade and Trias Llimós 2024,
#'   "Arriaga meets Kitagawa"): life expectancy (or the probability of
#'   dying) of the population is the mean of the subgroups' values, weighted
#'   by their population shares at the starting age (`at`), as when
#'   subgroup life tables are blended at the radix. A Kitagawa
#'   decomposition of that mean gives one composition effect and, for each
#'   subgroup, its difference times its average share, which is split by
#'   age (and cause) with `method`.
#' * `composition = "pooled"`: the quantity of the pooled rates, which
#'   weight the subgroup rates by their population shares at each age (the
#'   rates of the aggregate strata; Vaupel 2002, "current rates"). The rates
#'   of every subgroup, age and cause and the shares at every age are the
#'   parameters of a `"horiuchi"` or `"stepwise"` decomposition; the
#'   composition effects are reported by age, summed over subgroups, whose
#'   separate composition terms are not identified.
#'
#' Shares come from the exposures: of each observed year, and for forecast
#' years from `exposure` (default: the last observed year, so composition
#' does not change in the forecast unless exposures are given).
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
#' @param hierarchy Optional [mt_hierarchy()]: total and group levels of
#'   `cause` and `groups` are left out (they are sums of the others).
#' @param composition `"radix"` or `"pooled"` (see Subgroups).
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
#' @param seed Seed of the Poisson draws.
#' @return An `mt_decomposition` object. `summary()` returns contributions
#'   summed (path by path) over the variables not in `by`, with mean and
#'   quantiles; `summary(by = NULL)` returns one row per comparison (its
#'   number in `comparison`) with the values in the two populations, their
#'   difference (mean and quantiles), the sum of the contributions and the
#'   largest difference between the two over draws (`max_residual`).
#'   `print()` shows the totals.
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
#' Paper WP-2024-029. \doi{10.4054/MPIDR-WP-2024-029}
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
#' # the sex gap in a forecast year, with intervals
#' fc <- mt_forecast(mt_fit(lee_carter(), x), h = 10, n_draws = 100, seed = 1)
#' summary(mt_decompose(fc, var = "sex", ref = "male", at = 30), by = NULL)
#' # observed 2015 against the forecast of 2025
#' mt_decompose(fc, periods = c(2015, 2025), at = 30)
#' @export
mt_decompose <- function(x, y = NULL, var = NULL, ref = NULL, periods = NULL,
                         quantity = c("life_expectancy", "prob_dying", "asr", "crude_rate"), method = NULL,
                         cause = NULL, groups = NULL, hierarchy = NULL, composition = c("radix", "pooled"),
                         at = NULL, to = NULL, standard = "esp2013", per = 1e5, exposure = NULL, n_steps = 20,
                         n_draws = 1, probs = c(0.025, 0.5, 0.975), seed = 1) {
  composition <- match.arg(composition)
  custom <- is.function(quantity)
  qname <- if (custom) "custom" else match.arg(quantity)
  pooled <- !is.null(groups) && composition == "pooled"
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
  if (!is.null(groups) && qname != "crude_rate") {
    if (composition == "radix" && !qname %in% c("life_expectancy", "prob_dying")) {
      mt_abort("`composition = \"radix\"` applies to life expectancy and prob_dying; use \"pooled\".",
               class = "argument")
    }
    if (composition == "pooled" && method == "arriaga") {
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
  q <- decomp_quantity(qname, quantity, sx$age, at, to, standard, per)

  px <- decomp_populations(sx$strata, var, cause, groups, hierarchy)
  py <- if (!is.null(sy)) decomp_populations(sy$strata, var, cause, groups, hierarchy)
  if (!is.null(py) && (!identical(px$cause_levels, py$cause_levels) ||
                       !identical(px$group_levels, py$group_levels))) {
    mt_abort("`x` and `y` must have the same causes and groups.", class = "argument")
  }
  comps <- decomp_comparisons(sx, sy, px, py, var, ref, periods)

  cache <- list()
  side <- function(src, tag, p, idx) {
    key <- paste(tag, p)
    if (is.null(cache[[key]])) cache[[key]] <<- decomp_period(src, p)
    a <- cache[[key]]
    nc <- nrow(idx)
    ng <- ncol(idx)
    R <- a$R[, as.vector(idx), , drop = FALSE]
    list(R = array(R, c(dim(R)[1], nc, ng, dim(R)[3])), E = a$E[, idx[1, ], drop = FALSE])
  }
  tables <- list()
  paths <- list()
  totals <- list()
  for (k in seq_along(comps$list)) {
    cmp <- comps$list[[k]]
    s1 <- side(cmp$src1, cmp$tag1, cmp$p1, cmp$idx1)
    s2 <- side(cmp$src2, cmp$tag2, cmp$p2, cmp$idx2)
    n <- max(dim(s1$R)[4], dim(s2$R)[4])
    s1$R <- recycle_draws(s1$R, n)
    s2$R <- recycle_draws(s2$R, n)
    res <- decomp_core(s1, s2, q, method, composition, grouped = !is.null(groups), n_steps)
    tab <- data.table(comparison = k, component = res$rows$component,
                      group = if (is.null(groups)) NA_character_ else px$group_levels[res$rows$group_i],
                      age = sx$age$labels[q$rows][res$rows$age_i],
                      cause = if (is.null(cause)) NA_character_ else px$cause_levels[res$rows$cause_i])
    tables[[k]] <- tab
    paths[[k]] <- res$paths
    totals[[k]] <- list(v1 = res$v1, v2 = res$v2)
  }
  n_max <- max(vapply(paths, ncol, integer(1)))
  P <- do.call(rbind, lapply(paths, recycle_cols, n_max))
  structure(list(table = rbindlist(tables), paths = P, meta = comps$meta,
                 from = do.call(rbind, lapply(totals, function(t) recycle_cols(matrix(t$v1, 1), n_max))),
                 to = do.call(rbind, lapply(totals, function(t) recycle_cols(matrix(t$v2, 1), n_max))),
                 quantity = qname, method = method, composition = if (!is.null(groups)) composition,
                 at = q$at, probs = probs, age_labels = sx$age$labels[q$rows],
                 cause_levels = px$cause_levels, group_levels = px$group_levels),
            class = "mt_decomposition")
}

# Rates and exposures of a source: observed data (Poisson draws) or draws
decomp_source <- function(obj, n_draws, seed, exposure, name) {
  if (inherits(obj, "mt_data")) {
    d <- mt_observed(obj, n_draws = n_draws, seed = seed)
    return(list(age = d$age, strata = d$strata, periods = d$periods, rates = d$rates, exposure = d$exposure,
                obs_periods = NULL, label = "observed"))
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
       label = obj$model)
}

# Rates [age, stratum, draw] and exposures [age, stratum] of one period
decomp_period <- function(src, p) {
  d <- dim(src$rates)
  k <- match(p, src$periods)
  if (!is.na(k)) {
    return(list(R = array(src$rates[, k, , ], c(d[1], d[3], d[4])), E = matrix(src$exposure[, k, ], d[1])))
  }
  k <- match(p, src$obs_periods)
  if (!is.na(k) && !is.null(src$observed)) {
    return(list(R = array(src$observed[, k, ], c(d[1], d[3], 1)), E = matrix(src$obs_exposure[, k, ], d[1])))
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
decomp_quantity <- function(qname, quantity, age, at, to, standard, per) {
  na <- length(age)
  i0 <- if (is.null(at)) 1L else match(at, age$lower)
  if (is.na(i0)) mt_abort("`at` must be the start of an age interval.", class = "argument")
  sub_age <- function(rows) list(lower = age$lower[rows], upper = age$upper[rows])
  switch(qname,
         life_expectancy = {
           rows <- seq(i0, na)
           ag <- sub_age(rows)
           list(rows = rows, at = age$lower[i0], fun = function(M) lt_columns(M, ag)$ex[1, ], le = TRUE, ag = ag)
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
             Q <- pmin(nn * M / (1 + nn / 2 * M), 1)
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

# Populations: strata grouped by everything except causes and groups
decomp_populations <- function(st, var, cause, groups, hierarchy) {
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

decomp_comparisons <- function(sx, sy, px, py, var, ref, periods) {
  last_any <- function(s) max(c(s$periods, s$obs_periods))
  out <- list()
  meta <- list()
  add <- function(m, src1, tag1, p1, idx1, src2, tag2, p2, idx2) {
    out[[length(out) + 1L]] <<- list(src1 = src1, tag1 = tag1, p1 = p1, idx1 = idx1, src2 = src2, tag2 = tag2,
                                     p2 = p2, idx2 = idx2)
    meta[[length(meta) + 1L]] <<- m
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
  if (!length(out)) mt_abort("No comparisons: check `var`, `ref` and the strata of `y`.", class = "argument")
  meta <- rbindlist(meta)
  meta[, comparison := seq_len(.N)]
  list(list = out, meta = meta)
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
    v1 <- colSums(matrix(M1 * as.vector(c1), nr * ng)) * per
    v2 <- colSums(matrix(M2 * as.vector(c2), nr * ng)) * per
    tab <- rbind(label_rows("rate", rep(seq_len(nr), nc * ng), rep(rep(seq_len(nc), each = nr), ng),
                            rep(seq_len(ng), each = nr * nc)),
                 label_rows("composition", rep(seq_len(nr), nc), rep(seq_len(nc), each = nr), NA_integer_))
    return(list(rows = tab, paths = rbind(matrix(rate, ncol = n), matrix(compo, ncol = n)), v1 = v1, v2 = v2))
  }
  # decomposition of one homogeneous population: [age, cause, draw]
  by_age <- function(Ra, Rb, Ma, Mb) {
    if (method == "arriaga") {
      dl <- (arriaga_dir(Ma, Mb, q$ag) - arriaga_dir(Mb, Ma, q$ag)) / 2
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
    return(list(rows = tab, paths = matrix(cont, ncol = n), v1 = q$fun(matrix(M1[, 1, ], nr)),
                v2 = q$fun(matrix(M2[, 1, ], nr))))
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
    return(list(rows = tab, paths = rbind(do.call(rbind, paths), matrix(kt$composition, 1)),
                v1 = colSums(q1 * w1), v2 = colSums(q2 * w2)))
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
  list(rows = tab, paths = rbind(out[seq_len(nrate), , drop = FALSE], compo), v1 = f(P1), v2 = f(P2))
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
arriaga_dir <- function(M1, M2, ag) {
  a <- lt_columns(M1, ag)
  b <- lt_columns(M2, ag)
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
# dl / (Mb - Ma); a numerical derivative where the all-cause rates hardly differ
cause_split <- function(dl, Ra, Rb, Ma, Mb, fun) {
  dM <- Mb - Ma
  s <- dl / dM
  bad <- !is.finite(s) | abs(dM) <= 1e-8 * (Ma + Mb)
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
#' @param ... Unused.
#' @export
summary.mt_decomposition <- function(object, by = c("component", "group", "age", "cause"), probs = object$probs,
                                     ...) {
  if (is.null(by)) return(decomp_totals(object, probs)[])
  out <- decomp_summary(object, by, probs)
  out[, comparison := NULL]
  out[]
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
  out[, `:=`(sum_contributions = rowMeans(tot), max_residual = apply(abs(dif - tot), 1, max))]
  setcolorder(out, "comparison")
  out
}

#' @export
print.mt_decomposition <- function(x, ...) {
  what <- if (x$quantity == "life_expectancy") sprintf("life expectancy at %s", x$at) else x$quantity
  cat(sprintf("<mt_decomposition> %s, method %s%s; %d comparison(s), %d draw(s)\n", what, x$method,
              if (!is.null(x$composition)) sprintf(", %s composition", x$composition) else "",
              nrow(x$meta), ncol(x$paths)))
  tot <- summary(x, by = NULL)
  num <- setdiff(names(tot)[vapply(tot, is.numeric, logical(1))], c("comparison", "period", "max_residual"))
  if (ncol(x$paths) == 1) tot[, (grep("^difference_q", names(tot), value = TRUE)) := NULL]
  num <- intersect(num, names(tot))
  tot[, (num) := lapply(.SD, round, 3), .SDcols = num]
  tot[, max_residual := signif(max_residual, 2)]
  print(tot)
  invisible(x)
}

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
#' @param comparisons Numbers of the comparisons to show (rows of
#'   `summary(decomposition, by = NULL)`; default: all, as panels).
#' @return A ggplot object.
#' @examples
#' x <- mt_example(deaths = "single", years = 2000:2015)
#' mt_plot_decomposition(mt_decompose(x, periods = c(2000, 2015), at = 30))
#' @export
mt_plot_decomposition <- function(decomposition, by = c("age", "cause", "group"), fill = NULL,
                                  comparisons = NULL) {
  need_ggplot()
  check_class(decomposition, "mt_decomposition", "decomposition")
  by <- match.arg(by)
  tab <- decomposition$table
  if (all(is.na(tab[[by]]))) mt_abort(sprintf("The decomposition has no `%s`.", by), class = "argument")
  if (!is.null(fill) && (identical(fill, by) || !fill %in% c("cause", "group") || all(is.na(tab[[fill]])))) {
    mt_abort("`fill` must be another variable of the decomposition (\"cause\" or \"group\").", class = "argument")
  }
  meta <- copy(decomposition$meta)
  keep <- comparisons %||% meta$comparison
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
  p <- p + ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
    ggplot2::facet_wrap(ggplot2::vars(panel)) +
    ggplot2::labs(x = NULL, y = sprintf("Contribution to the difference in %s", decomp_label(decomposition)),
                  fill = NULL,
                  subtitle = sprintf("Method: %s%s", decomposition$method,
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

decomp_label <- function(d) {
  switch(d$quantity, life_expectancy = sprintf("life expectancy at %s (years)", d$at),
         prob_dying = "the probability of dying", asr = "the age-standardised rate",
         crude_rate = "the crude death rate", "the quantity")
}
