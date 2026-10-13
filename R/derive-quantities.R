#' Derived quantities from forecast draws or observed data
#'
#' All quantities are computed path by path, so their intervals reflect the
#' joint uncertainty across ages and years. For observed data (an [mt_data]
#' object) the paths are Poisson draws of the observed death counts (see
#' [mt_observed()]), so the intervals are those of the observed rates.
#'
#' * `prob_dying()`: probability of dying between ages `from` and `to` (for
#'   example the SDG indicator 3.4.1 uses 30 to 70), `1 - l(to) / l(from)`
#'   from the period life table;
#' * `asr()`: age-standardised death rate per `per` persons (default
#'   100,000) with weights of a standard population, restricted to the ages
#'   of the grid and renormalised (columns `ages`, `standard` and `per`
#'   state this; an ASR for ages 50+ is per 100,000 of the standard
#'   population aged 50+);
#' * `yll()`: years of life lost, deaths times a reference remaining life
#'   expectancy at each age (by default the forecast's own period life
#'   expectancy at the interval midpoint, approximated by the value at the
#'   start of the interval).
#'
#' @param draws An [mt_draws] object; `prob_dying()` and `asr()` also take an
#'   [mt_data] object for observed rates (with Poisson intervals, see
#'   [mt_observed()]).
#' @param from,to Ages (starts of grid intervals; `to` may be the end of the
#'   last closed interval).
#' @param standard The name of a standard population (see
#'   [standard_population()]) or a data frame with `age_lo`, `age_hi`,
#'   `weight`.
#' @param exposure Future exposures for `yll()` (see [mt_deaths()]).
#' @param per Rates per `per` persons (age-standardised rates).
#' @param reference_ex Reference remaining life expectancy per grid age, or
#'   `NULL` for the forecast's own.
#' @param probs Quantile levels.
#' @param paths Return the values of every path, an array `[period, stratum,
#'   draw]`, instead of the summary (for example to compute ratios between
#'   strata path by path).
#' @return A `data.table` with stratum columns, `period`, `mean` and
#'   quantiles, or with `paths = TRUE` an array.
#' @examples
#' x <- mt_example(deaths = "single")
#' fc <- mt_forecast(mt_fit(lee_carter(), x), h = 5, n_draws = 100, seed = 1)
#' prob_dying(fc, from = 30, to = 70)
#' asr(fc)
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
asr <- function(draws, standard = "esp2013", probs = c(0.025, 0.1, 0.5, 0.9, 0.975), paths = FALSE, per = 1e5) {
  draws <- as_draws(draws)
  if (is.character(standard)) standard <- standard_population(standard)
  w <- standard_weights(draws$age, standard)
  d <- dim(draws$rates)
  M <- matrix(draws$rates, nrow = d[1])
  v <- array(colSums(M * w) / sum(w) * per, d[-1])
  out <- summarise_paths(v, draws, probs, paths = paths)
  if (paths) return(out)
  a <- draws$age
  out[, `:=`(ages = sprintf("%s-%s", a$labels[1], a$labels[length(a)]), standard = attr(standard, "name") %||% "user",
             per = per)]
  out[]
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
  if (is.null(reference_ex)) {
    M <- matrix(draws$rates, nrow = d[1])
    ex <- lt_columns(M, draws$age)$ex
  } else {
    if (length(reference_ex) != d[1]) mt_abort("`reference_ex` needs one value per grid age.", class = "argument")
    ex <- matrix(reference_ex, d[1], prod(d[-1]))
  }
  v <- array(colSums(matrix(deaths, nrow = d[1]) * ex), d[-1])
  summarise_paths(v, draws, probs)
}

#' Export forecasts and result tables
#'
#' Writes a tidy CSV of forecast draws (their summary, `what = "summary"`,
#' or all draws in long format, `what = "draws"`) or of any result table,
#' such as those of [life_expectancy()], [asr()], [mt_contrast()],
#' [mt_inequality()] or [mt_score()]. The file name of a forecast records
#' the model, the grid and the forecast origin unless `file` is given.
#'
#' @param x An [mt_draws] object or a data frame.
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
  if (is.data.frame(x)) {
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
