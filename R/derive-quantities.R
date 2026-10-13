#' Derived quantities from forecast draws
#'
#' All quantities are computed path by path, so their intervals reflect the
#' joint uncertainty across ages and years.
#'
#' * `prob_dying()`: probability of dying between ages `from` and `to` (for
#'   example the SDG indicator 3.4.1 uses 30 to 70), `1 - l(to) / l(from)`
#'   from the period life table;
#' * `asr()`: age-standardised death rate per 100,000 with weights of a
#'   standard population, restricted to the ages of the grid and renormalised;
#' * `yll()`: years of life lost, deaths times a reference remaining life
#'   expectancy at each age (by default the forecast's own period life
#'   expectancy at the interval midpoint, approximated by the value at the
#'   start of the interval).
#'
#' @param draws An [mt_draws] object.
#' @param from,to Ages (starts of grid intervals; `to` may be the end of the
#'   last closed interval).
#' @param standard A standard population from [standard_population()] or a
#'   data frame with `age_lo`, `age_hi`, `weight`.
#' @param exposure Future exposures for `yll()` (see [mt_deaths()]).
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
  check_class(draws, "mt_draws", "draws")
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
asr <- function(draws, standard = standard_population("esp2013"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975),
                paths = FALSE) {
  check_class(draws, "mt_draws", "draws")
  w <- standard_weights(draws$age, standard)
  d <- dim(draws$rates)
  M <- matrix(draws$rates, nrow = d[1])
  v <- array(colSums(M * w) / sum(w) * 1e5, d[-1])
  summarise_paths(v, draws, probs, paths = paths)
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
#' `"esp2013"`: the revised European Standard Population 2013 (Eurostat 2013,
#' Revision of the European Standard Population, report KS-RA-13-028), in
#' 5-year groups with 90-94 and 95+, total 100,000.
#'
#' @param name Name of the standard.
#' @return A `data.table` with `age_lo`, `age_hi`, `weight`.
#' @examples
#' standard_population()
#' @export
standard_population <- function(name = "esp2013") {
  check_choice(name, "esp2013", "name")
  lo <- c(0, 1, seq(5, 95, 5))
  data.table(age_lo = lo, age_hi = c(lo[-1], Inf),
             weight = c(1000, 4000, 5500, 5500, 5500, 6000, 6000, 6500, 7000, 7000, 7000, 7000, 6500, 6000,
                        5500, 5000, 4000, 2500, 1500, 800, 200))
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

#' Export forecasts
#'
#' Writes a tidy CSV of the summary (`what = "summary"`) or of all draws in
#' long format (`what = "draws"`). The file name records the model, the grid
#' and the forecast origin unless `file` is given.
#'
#' @param draws An [mt_draws] object.
#' @param dir Output directory.
#' @param what `"summary"` or `"draws"`.
#' @param file Optional file name.
#' @param probs Quantiles for the summary.
#' @return The path of the written file, invisibly.
#' @examples
#' x <- mt_example(deaths = "single", years = 2005:2019)
#' fc <- mt_forecast(mt_fit(rwd(), x), h = 2, n_draws = 10, seed = 1)
#' f <- mt_export(fc, tempdir())
#' head(read.csv(f), 3)
#' @export
mt_export <- function(draws, dir, what = c("summary", "draws"), file = NULL,
                      probs = c(0.025, 0.1, 0.5, 0.9, 0.975)) {
  check_class(draws, "mt_draws", "draws")
  what <- match.arg(what)
  if (!dir.exists(dir)) mt_abort("`dir` does not exist.", class = "argument")
  a <- draws$age
  grid_label <- sprintf("ages%s-%s", a$labels[1], gsub("\\+", "plus", a$labels[length(a)]))
  file <- file %||% sprintf("%s_%s_origin%s_%s.csv", gsub("[^A-Za-z0-9]+", "-", tolower(draws$model)), grid_label,
                            draws$origin, what)
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
  fwrite(out, path)
  invisible(path)
}
