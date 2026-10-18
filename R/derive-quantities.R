#' Derived quantities from forecast draws or observed data
#'
#' All quantities are computed path by path, so their intervals reflect the
#' joint uncertainty across ages and years. For observed data (an [mt_data]
#' object) the paths are Poisson draws of the observed death counts (see
#' [mt_observed()]), so the intervals are those of the observed rates.
#'
#' * `prob_dying()`: probability of dying between ages `from` and `to` (for
#'   example the SDG indicator 3.4.1 uses 30 to 70), `1 - l(to) / l(from)`
#'   from the period life table of [life_table()] (Preston, Heuveline and
#'   Guillot 2001);
#' * `asr()`: age-standardised death rate per `per` persons (default
#'   100,000), by direct standardisation with the weights of a standard
#'   population (see Age-standardised rates);
#' * `yll()`: years of life lost, the deaths of each age interval times the
#'   remaining life expectancy at the average age at death in the interval
#'   (see Years of life lost).
#'
#' @section Age-standardised rates:
#' By default (`age_groups = "standard"`) the rates are first combined to
#' the age groups of the standard population, as for the European Standard
#' Population (Eurostat 2013): the rate of a standard age group is the sum
#' of the deaths divided by the sum of the exposures of the grid cells in
#' it, that is the mean of their rates weighted by their exposures. The ASR
#' is the mean of these rates weighted by the standard population. The
#' exposures are those of each observed year for observed data and, for
#' forecasts, `exposure` (default: the exposures of the last training year,
#' as in [mt_deaths()]). Rates of a group without exposures (all zero, or
#' draws that hold none) are averaged with equal weights.
#'
#' Grid cells wider than the standard groups (for example 0-4 for the groups
#' 0 and 1-4 of the ESP 2013, or an open interval 90+ for the groups 90-94
#' and 95+) take the sum of the weights of the standard groups they cover:
#' rates and weights are combined over the smallest age groups that are
#' unions both of grid cells and of standard groups. Only the ages of the
#' grid enter, and the weights are renormalised to them (columns `ages`,
#' `standard` and `per` state this; an ASR for ages 50+ is per 100,000 of
#' the standard population aged 50+); a standard group that the grid covers
#' only in part enters with the covered share of its weight.
#'
#' `age_groups = "grid"` instead shares the weight of each standard group
#' evenly among the grid cells within it (the weight of the open standard
#' group equally among the grid cells it covers) and weights the rates of
#' the grid cells with these shares. This is specific to this package
#' (reference unknown). Where the grid is finer than the standard groups,
#' every cell of a group gets the same weight, although the older cells of
#' a group have higher death rates and, at older ages, fewer persons, so the
#' ASR is higher than with the standard's age groups: by 1.8% to 3.6% for
#' French rates on single years of age to 100+ (HMD, 1950-2006, ESP 2013),
#' and by more when the grid extends to higher ages, over which the weight
#' of the open standard group 95+ is then spread evenly (up to 11% with
#' single years to 110+). On a grid of the standard's age groups (or
#' coarser) both give the same value.
#'
#' @section Years of life lost:
#' `yll()` multiplies the deaths of each age interval by the remaining life
#' expectancy at the average age at death in the interval, `x + a_x`,
#' interpolated linearly between the values at the start and at the end of
#' the interval, `e_x + a_x / n_x * (e_(x+n) - e_x)` (as for e-dagger in
#' [lifespan_variation()]); the last interval uses its `e_x`. With the life
#' expectancy of the forecast's own period life table (default) these are
#' the period expected years of life lost of Murray (1994). With
#' `reference_ex`, the remaining life expectancy of a standard life table at
#' the start of each grid interval, interpolated in the same way with the
#' `a_x` of the forecast's life table, they are the standard expected years
#' of life lost of the same paper. Neither discounting nor age weighting is
#' applied. The deaths are rates times `exposure` (default: the exposures of
#' the last training year), without Poisson noise.
#'
#' @param draws An [mt_draws] object; `prob_dying()` and `asr()` also take an
#'   [mt_data] object for observed rates (with Poisson intervals, see
#'   [mt_observed()]).
#' @param from,to Ages (starts of grid intervals; `to` may be the end of the
#'   last closed interval).
#' @param standard The name of a standard population (see
#'   [standard_population()]) or a data frame with `age_lo`, `age_hi`,
#'   `weight`.
#' @param exposure Exposures (see [mt_deaths()]): for `yll()` the future
#'   exposures that give the deaths; for `asr()` the exposures that weight
#'   the rates within each standard age group (default: those of each
#'   observed year for observed data, of the last training year for
#'   forecasts).
#' @param per Rates per `per` persons (age-standardised rates).
#' @param age_groups `"standard"` (default: direct standardisation on the
#'   age groups of the standard population) or `"grid"` (the standard
#'   weights shared among the grid cells, specific to this package). See
#'   Age-standardised rates.
#' @param reference_ex Reference remaining life expectancy at the start of
#'   each grid interval, or `NULL` for the forecast's own.
#' @param probs Quantile levels.
#' @param paths Return the values of every path, an array `[period, stratum,
#'   draw]`, instead of the summary (for example to compute ratios between
#'   strata path by path).
#' @return A `data.table` with stratum columns, `period`, `mean` and
#'   quantiles, or with `paths = TRUE` an array.
#' @references
#' Eurostat (2013). Revision of the European Standard Population: report of
#' Eurostat's task force. Methodologies and Working papers. Publications
#' Office of the European Union, Luxembourg. \doi{10.2785/11470}
#'
#' Murray CJL (1994). Quantifying the burden of disease: the technical
#' basis for disability-adjusted life years. *Bulletin of the World Health
#' Organization* 72(3): 429-445.
#'
#' Preston SH, Heuveline P, Guillot M (2001). *Demography: Measuring and
#' Modeling Population Processes*. Blackwell.
#' @examples
#' x <- mt_example(deaths = "single")
#' fc <- mt_forecast(mt_fit(lee_carter(), x), h = 5, n_draws = 100, seed = 1)
#' prob_dying(fc, from = 30, to = 70)
#' asr(fc)
#' yll(fc)
#' @name derived
NULL

summarise_paths <- function(arr, draws, probs, paths = FALSE) {
  d <- dim(arr)
  if (paths) {
    dimnames(arr) <- list(period = draws$periods, stratum = draws$strata$label, draw = NULL)
    return(arr)
  }
  M <- matrix(arr, ncol = d[3])
  out <- data.table(stratum = rep(seq_len(d[2]), each = d[1]), period = rep(draws$periods, d[2]),
                    mean = rowMeans(M))
  ok <- rowSums(is.na(M)) == 0
  Q <- matrix(NA_real_, nrow(M), length(probs))
  if (any(ok)) Q[ok, ] <- row_quantiles(M[ok, , drop = FALSE], probs)
  qn <- qnames(probs)
  for (k in seq_along(probs)) set(out, j = qn[k], value = Q[, k])
  out <- merge(draws$strata, out, by = "stratum", sort = FALSE)
  setorderv(out, c("stratum", "period"))
  out[]
}

#' @rdname derived
#' @export
prob_dying <- function(draws, from, to, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), paths = FALSE) {
  draws <- as_draws(draws)
  a <- draws$age
  i <- match(from, a$lower)
  j <- match(to, a$lower)
  if (is.na(j) && to == a$upper[length(a)]) j <- length(a) + 1L
  if (is.na(i) || is.na(j) || j <= i) mt_abort("`from` and `to` must be interval bounds of the grid with from < to.",
                                               class = "argument")
  d <- dim(draws$rates)
  M <- matrix(draws$rates, nrow = d[1])
  lt <- lt_columns(M, a)
  l_from <- lt$lx[i, ]
  l_to <- if (j > d[1]) lt$lx[d[1], ] - lt$dx[d[1], ] else lt$lx[j, ]
  q <- array(1 - l_to / l_from, d[-1])
  if (paths) return(summarise_paths(q, draws, probs, paths = TRUE))
  out <- summarise_paths(q, draws, probs)
  out[, `:=`(from = from, to = to)]
  out[]
}

#' @rdname derived
#' @export
asr <- function(draws, standard = "esp2013", probs = c(0.025, 0.1, 0.5, 0.9, 0.975), paths = FALSE, per = 1e5,
                exposure = NULL, age_groups = c("standard", "grid")) {
  age_groups <- match.arg(age_groups)
  draws <- as_draws(draws)
  if (is.character(standard)) standard <- standard_population(standard)
  d <- dim(draws$rates)
  M <- matrix(draws$rates, nrow = d[1])
  E <- if (age_groups == "standard") draws_exposure(draws, exposure)
  v <- array(asr_columns(M, draws$age, standard, E, per, age_groups), d[-1])
  out <- summarise_paths(v, draws, probs, paths = paths)
  if (paths) return(out)
  a <- draws$age
  out[, `:=`(ages = sprintf("%s-%s", a$labels[1], a$labels[length(a)]), standard = attr(standard, "name") %||% "user",
             per = per)]
  out[]
}

# Exposures [age, period, stratum] of draws: `exposure` if given, those of
# the observed years of mt_observed(), else those of the last training year
# (the draws' own, matched by position if unnamed); NULL when the draws hold
# none that fit
draws_exposure <- function(draws, exposure = NULL) {
  d <- dim(draws$rates)
  if (!is.null(exposure)) return(expand_exposure(exposure, d, dimnames(draws$rates)))
  if (!is.null(draws$exposure) && identical(dim(draws$exposure), d[1:3])) return(draws$exposure)
  if (is.null(draws$exposure_last)) return(NULL)
  tryCatch(suppressWarnings(expand_exposure(draws$exposure_last, d, dimnames(draws$rates))),
           mt_error_argument = function(e) NULL)
}

# Age-standardised rates of the columns of M (rates [age, n]) with the
# exposures E (a vector or matrix [age, k] recycled over the columns of M,
# n a multiple of k; NULL for equal weights): the rates combined to the age
# groups of asr_groups() (sum of rate times exposure over the sum of the
# exposures; the plain mean where the exposures are zero), weighted by the
# standard population
asr_columns <- function(M, age, standard, E = NULL, per = 1e5, age_groups = "standard") {
  g <- asr_groups(age, standard, age_groups)
  rows <- which(!is.na(g$group))
  M <- M[rows, , drop = FALSE]
  grp <- g$group[rows]
  E <- if (is.null(E)) matrix(1, length(rows), 1) else matrix(E, length(age))[rows, , drop = FALSE]
  den <- rowsum(E, grp, reorder = TRUE)
  R <- rowsum(M * as.vector(E), grp, reorder = TRUE) / as.vector(den)
  zero <- !(den > 0)
  if (any(zero)) {
    plain <- rowsum(M, grp, reorder = TRUE) / tabulate(grp)
    z <- rep_len(as.vector(zero), length(R))
    R[z] <- plain[z]
  }
  colSums(R * g$weight) / sum(g$weight) * per
}

# Age groups of an age-standardised rate: `group` (the group of each grid
# cell, NA outside the standard) and per group its standard `weight` and its
# bounds `lower` and `upper`. "standard": the smallest groups that are
# unions both of grid cells and of standard groups (the standard's age
# groups where the grid is finer); a standard group covered in part by the
# grid counts with the covered share of its weight. "grid": every grid cell
# with the weights of standard_weights()
asr_groups <- function(age, standard, age_groups = "standard") {
  standard <- as.data.table(standard)
  lo <- age$lower
  hi <- age$upper
  if (age_groups == "grid") {
    w <- standard_weights(age, standard)
    grp <- cumsum(w > 0)
    grp[w <= 0] <- NA_integer_
    return(list(group = grp, weight = w[w > 0], lower = lo[w > 0], upper = hi[w > 0]))
  }
  slo <- standard$age_lo
  shi <- standard$age_hi
  ov <- outer(seq_along(lo), seq_along(slo), function(i, k) pmin(hi[i], shi[k]) - pmax(lo[i], slo[k]))
  hit <- ov > 1e-9
  if (!any(hit)) mt_abort("The standard population does not cover the grid ages.", class = "argument")
  # neighbouring cells share a group when a standard group spans their bound
  grp <- rep(NA_integer_, length(lo))
  cur <- 0L
  for (i in seq_along(lo)) {
    if (!any(hit[i, ])) next
    if (i == 1L || is.na(grp[i - 1L]) || !any(hit[i, ] & hit[i - 1L, ])) cur <- cur + 1L
    grp[i] <- cur
  }
  covered <- colSums(ifelse(hit, ov, 0))
  share <- ifelse(is.finite(shi), pmin(covered / (shi - slo), 1), 1)
  share[colSums(hit) == 0] <- 0
  sg <- apply(hit, 2, function(h) if (any(h)) grp[which(h)[1]] else NA_integer_)
  weight <- vapply(seq_len(cur), function(k) sum((standard$weight * share)[which(sg == k)]), 0)
  list(group = grp, weight = weight, lower = vapply(seq_len(cur), function(k) min(lo[grp == k], na.rm = TRUE), 0),
       upper = vapply(seq_len(cur), function(k) max(hi[grp == k], na.rm = TRUE), 0))
}

# Weight of each grid cell: the standard weights of the overlapping standard
# groups, shared in proportion to the overlap (an open standard group is
# shared equally by the grid cells it covers; an open grid cell collects all
# standard groups above its lower bound)
standard_weights <- function(age, standard) {
  standard <- as.data.table(standard)
  w <- numeric(length(age))
  for (k in seq_len(nrow(standard))) {
    lo <- standard$age_lo[k]
    hi <- standard$age_hi[k]
    ov_lo <- pmax(age$lower, lo)
    ov_hi <- pmin(age$upper, hi)
    inside <- ov_hi > ov_lo
    if (!any(inside)) next
    width <- ov_hi[inside] - ov_lo[inside]
    if (is.finite(hi)) {
      w[inside] <- w[inside] + standard$weight[k] * width / (hi - lo)
    } else {
      width[!is.finite(width)] <- 1
      w[inside] <- w[inside] + standard$weight[k] * width / sum(width)
    }
  }
  if (sum(w) == 0) mt_abort("The standard population does not cover the grid ages.", class = "argument")
  w
}

#' Standard populations
#'
#' * `"esp2013"`: the revised European Standard Population 2013 (Eurostat
#'   2013, Revision of the European Standard Population, report
#'   KS-RA-13-028), in 5-year groups with 0, 1-4 and 95+, total 100,000;
#' * `"esp1976"`: the European Standard Population of 1976 (0, 1-4, 5-year
#'   groups, 85+), total 100,000;
#' * `"who2000"`: the WHO World Standard Population, world average
#'   2000-2025 (Ahmad et al. 2001), 5-year groups to 100+; the published
#'   weights sum to 100,035 and are used as published (rates are divided by
#'   the sum of the weights).
#'
#' A national standard is used by passing a data frame with `age_lo`,
#' `age_hi` (`Inf` for the open group) and `weight` to the `standard`
#' argument of [asr()].
#'
#' @param name Name of the standard.
#' @return A `data.table` with `age_lo`, `age_hi`, `weight`; the name is in
#'   attribute `"name"`.
#' @references
#' Ahmad OB, Boschi-Pinto C, Lopez AD, Murray CJL, Lozano R, Inoue M (2001).
#' Age standardization of rates: a new WHO standard. GPE Discussion Paper
#' Series No. 31. World Health Organization, Geneva.
#' @examples
#' standard_population()
#' standard_population("who2000")
#' @export
standard_population <- function(name = c("esp2013", "esp1976", "who2000")) {
  name <- match.arg(name)
  out <- switch(name,
    esp2013 = {
      lo <- c(0, 1, seq(5, 95, 5))
      data.table(age_lo = lo, age_hi = c(lo[-1], Inf),
                 weight = c(1000, 4000, 5500, 5500, 5500, 6000, 6000, 6500, 7000, 7000, 7000, 7000, 6500, 6000,
                            5500, 5000, 4000, 2500, 1500, 800, 200))
    },
    esp1976 = {
      lo <- c(0, 1, seq(5, 85, 5))
      data.table(age_lo = lo, age_hi = c(lo[-1], Inf),
                 weight = c(1600, 6400, rep(7000, 10), 6000, 5000, 4000, 3000, 2000, 1000, 1000))
    },
    who2000 = {
      lo <- seq(0, 100, 5)
      data.table(age_lo = lo, age_hi = c(lo[-1], Inf),
                 weight = c(8860, 8690, 8600, 8470, 8220, 7930, 7610, 7150, 6590, 6040, 5370, 4550, 3720, 2960,
                            2210, 1520, 910, 440, 150, 40, 5))
    })
  setattr(out, "name", name)
  out
}

#' @rdname derived
#' @export
yll <- function(draws, exposure = NULL, reference_ex = NULL, probs = c(0.025, 0.1, 0.5, 0.9, 0.975)) {
  check_class(draws, "mt_draws", "draws")
  deaths <- mt_deaths(draws, exposure = exposure, noise = FALSE)
  d <- dim(draws$rates)
  M <- matrix(draws$rates, nrow = d[1])
  lt <- lt_columns(M, draws$age)
  if (is.null(reference_ex)) {
    ex <- lt$ex
  } else {
    if (length(reference_ex) != d[1]) mt_abort("`reference_ex` needs one value per grid age.", class = "argument")
    ex <- matrix(reference_ex, d[1], prod(d[-1]))
  }
  v <- array(colSums(matrix(deaths, nrow = d[1]) * ex_at_death(ex, lt$ax, draws$age)), d[-1])
  summarise_paths(v, draws, probs)
}

# Remaining life expectancy at the average age at death x + a_x of each
# interval, interpolated linearly between e_x and e_(x+n) (ex and ax
# [age, n]); the last interval keeps its e_x
ex_at_death <- function(ex, ax, age) {
  na <- nrow(ex)
  if (na == 1L) return(ex)
  n <- age$upper - age$lower
  i <- seq_len(na - 1L)
  out <- ex
  e0 <- ex[i, , drop = FALSE]
  out[i, ] <- e0 + ax[i, , drop = FALSE] / n[i] * (ex[i + 1L, , drop = FALSE] - e0)
  out
}

#' Export forecasts and result tables
#'
#' Writes a tidy CSV of forecast draws (their summary, `what = "summary"`,
#' or all draws in long format, `what = "draws"`) or of any result table,
#' such as those of [life_expectancy()], [asr()], [mt_contrast()],
#' [mt_inequality()] or [mt_score()]. The file name of a forecast records
#' the model, the grid and the forecast origin unless `file` is given.
#'
#' @param x An [mt_draws] object, an [mt_decompose()] result (its summary is
#'   written) or a data frame.
#' @param dir Output directory.
#' @param what For draws: `"summary"` or `"draws"`.
#' @param file Optional file name (required for tables other than draws
#'   unless the default `"mortensor_table.csv"` is wanted).
#' @param probs Quantiles for the summary of draws.
#' @param digits Optional number of significant digits of numeric columns.
#' @return The path of the written file, invisibly.
#' @examples
#' x <- mt_example(deaths = "single", years = 2005:2019)
#' fc <- mt_forecast(mt_fit(rwd(), x), h = 2, n_draws = 10, seed = 1)
#' f <- mt_export(fc, tempdir())
#' head(read.csv(f), 3)
#' f2 <- mt_export(life_expectancy(fc, at = 65), tempdir(), file = "e65.csv", digits = 4)
#' @export
mt_export <- function(x, dir, what = c("summary", "draws"), file = NULL,
                      probs = c(0.025, 0.1, 0.5, 0.9, 0.975), digits = NULL) {
  what <- match.arg(what)
  if (!dir.exists(dir)) mt_abort("`dir` does not exist.", class = "argument")
  if (inherits(x, "mt_decomposition")) {
    out <- summary(x, probs = probs)
    path <- file.path(dir, file %||% sprintf("decomposition_%s.csv", x$quantity))
  } else if (is.data.frame(x)) {
    out <- as.data.table(x)
    path <- file.path(dir, file %||% "mortensor_table.csv")
  } else {
    check_class(x, "mt_draws", "x")
    draws <- x
    a <- draws$age
    grid_label <- sprintf("ages%s-%s", a$labels[1], gsub("\\+", "plus", a$labels[length(a)]))
    file <- file %||% sprintf("%s_%s_origin%s_%s.csv", gsub("[^A-Za-z0-9]+", "-", tolower(draws$model)),
                              grid_label, draws$origin, what)
    path <- file.path(dir, file)
    out <- if (what == "summary") {
      mt_summarise(draws, probs)
    } else {
      d <- dim(draws$rates)
      long <- data.table(stratum = rep(rep(seq_len(d[3]), each = d[1] * d[2]), d[4]),
                         age_lo = rep(a$lower, d[2] * d[3] * d[4]), age = rep(a$labels, d[2] * d[3] * d[4]),
                         period = rep(rep(draws$periods, each = d[1]), d[3] * d[4]),
                         draw = rep(seq_len(d[4]), each = d[1] * d[2] * d[3]), rate = as.vector(draws$rates))
      merge(draws$strata, long, by = "stratum", sort = FALSE)
    }
  }
  if (!is.null(digits)) {
    num <- names(out)[vapply(out, is.double, logical(1))]
    if (length(num)) out[, (num) := lapply(.SD, signif, digits = digits), .SDcols = num]
  }
  fwrite(out, path)
  invisible(path)
}
