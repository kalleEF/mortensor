#' Age and period axes
#'
#' An axis is an ordered set of contiguous half-open intervals
#' `[lower, upper)`. Only the last interval may be open (`upper = Inf`), as
#' for an age group such as 90+. Axes describe the resolution of observed data
#' and of the modelling grid.
#'
#' @param lower,upper Numeric vectors of interval bounds. `upper` may end with
#'   `Inf` for an open last interval.
#' @param labels Optional character labels; generated when `NULL`.
#' @param kind `"age"` or `"period"`.
#' @return An object of class `mt_axis`.
#' @examples
#' single_ages(30, 100)
#' age_groups(c(0, 1, seq(5, 90, 5)))
#' single_years(2001, 2019)
#' @export
mt_axis <- function(lower, upper, labels = NULL, kind = c("age", "period")) {
  kind <- match.arg(kind)
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  if (length(lower) == 0L || length(lower) != length(upper)) {
    mt_abort("`lower` and `upper` must be non-empty and of equal length.", class = "axis")
  }
  if (anyNA(lower) || anyNA(upper)) mt_abort("Axis bounds must not be missing.", class = "axis")
  if (any(upper <= lower)) mt_abort("Every interval needs `upper > lower`.", class = "axis")
  n <- length(lower)
  if (n > 1L && any(abs(upper[-n] - lower[-1]) > 1e-9)) {
    mt_abort("Axis intervals must be contiguous and increasing (upper[i] == lower[i + 1]).",
             class = "axis")
  }
  if (any(is.infinite(upper[-n]))) mt_abort("Only the last interval may be open.", class = "axis")
  if (is.null(labels)) labels <- axis_default_labels(lower, upper, kind)
  if (length(labels) != n) mt_abort("`labels` must have one entry per interval.", class = "axis")
  structure(list(lower = lower, upper = upper, labels = as.character(labels), kind = kind),
            class = "mt_axis")
}

axis_default_labels <- function(lower, upper, kind) {
  fmt <- function(x) format(x, trim = TRUE, scientific = FALSE)
  w <- upper - lower
  ifelse(is.infinite(upper), paste0(fmt(lower), "+"),
         ifelse(abs(w - 1) < 1e-9, fmt(lower), paste0(fmt(lower), "-", fmt(upper - 1))))
}

#' @rdname mt_axis
#' @param from,to First and last value (inclusive).
#' @param open Whether the last interval is open (`to+`).
#' @export
single_ages <- function(from = 0, to = 100, open = TRUE) {
  check_number(from, "from", lower = 0, integerish = TRUE)
  check_number(to, "to", lower = from, integerish = TRUE)
  lower <- seq(from, to)
  upper <- c(lower[-1], if (open) Inf else to + 1)
  mt_axis(lower, upper, kind = "age")
}

#' @rdname mt_axis
#' @param breaks Lower bounds of the age groups.
#' @param last_upper Upper bound of the last group when `open = FALSE`.
#' @export
age_groups <- function(breaks, open = TRUE, last_upper = NULL) {
  breaks <- sort(as.numeric(breaks))
  if (!open && is.null(last_upper)) mt_abort("`last_upper` is needed when `open = FALSE`.", class = "axis")
  mt_axis(breaks, c(breaks[-1], if (open) Inf else last_upper), kind = "age")
}

#' @rdname mt_axis
#' @export
single_years <- function(from, to) {
  check_number(from, "from", integerish = TRUE)
  check_number(to, "to", lower = from, integerish = TRUE)
  lower <- seq(from, to)
  mt_axis(lower, lower + 1, kind = "period")
}

#' @rdname mt_axis
#' @export
period_groups <- function(breaks, last_upper) {
  breaks <- sort(as.numeric(breaks))
  mt_axis(breaks, c(breaks[-1], last_upper), kind = "period")
}

#' @export
length.mt_axis <- function(x) length(x$lower)

#' @export
format.mt_axis <- function(x, ...) {
  n <- length(x)
  shown <- if (n > 8) c(x$labels[1:4], "...", x$labels[(n - 2):n]) else x$labels
  sprintf("<mt_axis: %s, %d interval%s> %s", x$kind, n, if (n == 1) "" else "s",
          paste(shown, collapse = " "))
}

#' @export
print.mt_axis <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @export
as.data.frame.mt_axis <- function(x, ...) {
  data.frame(lower = x$lower, upper = x$upper, label = x$labels, stringsAsFactors = FALSE)
}

#' Interval widths and midpoints of an axis
#'
#' @param axis An [mt_axis].
#' @param omega Assumed upper bound used in place of `Inf` for an open last
#'   interval (default 110 for ages).
#' @return Numeric vector.
#' @export
axis_widths <- function(axis, omega = 110) {
  check_class(axis, "mt_axis", "axis")
  upper <- axis$upper
  if (is.infinite(upper[length(upper)])) {
    if (omega <= axis$lower[length(upper)]) mt_abort("`omega` must exceed the start of the open interval.", class = "axis")
    upper[length(upper)] <- omega
  }
  upper - axis$lower
}

#' @rdname axis_widths
#' @export
axis_midpoints <- function(axis, omega = 110) axis$lower + axis_widths(axis, omega) / 2

#' @rdname axis_widths
#' @export
is_open <- function(axis) is.infinite(axis$upper[length(axis$upper)])

#' Map intervals onto an axis
#'
#' For each interval `[lo, hi)`, returns the index of the axis cell that
#' contains it (`mode = "within"`) or that equals it (`mode = "exact"`), and
#' `NA` otherwise.
#'
#' @param axis An [mt_axis].
#' @param lo,hi Interval bounds.
#' @param mode `"within"` or `"exact"`.
#' @return Integer vector.
#' @keywords internal
axis_locate <- function(axis, lo, hi, mode = c("within", "exact")) {
  mode <- match.arg(mode)
  i <- findInterval(lo, axis$lower)
  i[i == 0L] <- NA_integer_
  ok <- !is.na(i)
  if (mode == "exact") {
    ok[ok] <- abs(axis$lower[i[ok]] - lo[ok]) < 1e-9 & abs(axis$upper[i[ok]] - hi[ok]) < 1e-9 |
      (is.infinite(axis$upper[i[ok]]) & is.infinite(hi[ok]) & abs(axis$lower[i[ok]] - lo[ok]) < 1e-9)
  } else {
    ok[ok] <- hi[ok] <= axis$upper[i[ok]] + 1e-9
  }
  i[!ok] <- NA_integer_
  i
}

#' Parse age-group labels
#'
#' Understands the label styles of common sources, for example `"<1"`,
#' `"01-14"`, `"30-34"`, `"90+"`, `"85 and over"`, `"under 1"`, `"0"`,
#' `"110+"`.
#'
#' @param x Character vector of labels.
#' @return A `data.table` with `label`, `age_lo`, `age_hi` (`Inf` for open
#'   intervals), in the order of `x`.
#' @examples
#' parse_age_labels(c("<1", "01-14", "15-19", "90+"))
#' @export
parse_age_labels <- function(x) {
  s <- tolower(trimws(as.character(x)))
  s <- gsub("\u2013|\u2014", "-", s)
  lo <- hi <- rep(NA_real_, length(s))
  num <- function(z) as.numeric(z)
  m <- regmatches(s, regexec("^(?:<|under)\\s*(\\d+)$", s, perl = TRUE))
  for (k in which(lengths(m) == 2)) { lo[k] <- 0; hi[k] <- num(m[[k]][2]) }
  m <- regmatches(s, regexec("^(\\d+)\\s*-\\s*(\\d+)$", s, perl = TRUE))
  for (k in which(lengths(m) == 3)) { lo[k] <- num(m[[k]][2]); hi[k] <- num(m[[k]][3]) + 1 }
  m <- regmatches(s, regexec("^(\\d+)\\s*(?:\\+|and over|and older|plus|or older)$", s, perl = TRUE))
  for (k in which(lengths(m) == 2)) { lo[k] <- num(m[[k]][2]); hi[k] <- Inf }
  m <- regmatches(s, regexec("^(\\d+)$", s, perl = TRUE))
  for (k in which(lengths(m) == 2)) { lo[k] <- num(m[[k]][2]); hi[k] <- lo[k] + 1 }
  bad <- is.na(lo)
  if (any(bad)) {
    mt_abort(sprintf("Cannot parse age label(s): %s.", paste0("\"", unique(x[bad]), "\"", collapse = ", ")),
             class = "axis")
  }
  data.table(label = as.character(x), age_lo = lo, age_hi = hi)
}

#' Build an axis from interval bounds found in data
#' @keywords internal
axis_from_bounds <- function(lo, hi, kind) {
  b <- unique(data.table(lo = lo, hi = hi))
  setorder(b, lo)
  if (anyDuplicated(b$lo)) {
    mt_abort(sprintf("Overlapping %s intervals with the same lower bound.", kind), class = "axis")
  }
  mt_axis(b$lo, b$hi, kind = kind)
}
