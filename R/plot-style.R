#' Plot style of mortensor
#'
#' The plots of mortensor share one layout: a minimal theme, observed values
#' as small transparent points, fitted values as solid lines, forecast
#' medians as dashed lines with light interval ribbons, and a dotted vertical
#' line between the last observed and the first forecast year. Strata with
#' ordered levels (such as deprivation quintiles `Q1` to `Q5`) are coloured
#' on a red-yellow-green scale, other strata (causes, regions, models) with
#' distinct colours.
#'
#' * `theme_mortensor()`: the theme ([ggplot2::theme_minimal()] with the
#'   legend at the bottom);
#' * `scale_colour_mortensor()`, `scale_fill_mortensor()`: the colour
#'   scales, for your own plots of mortensor results;
#' * `mt_palette()`: the colours themselves.
#'
#' The ordered scale runs from red (first level) to green (last level), with
#' a dark yellow in the middle that stays visible on white; `reverse = TRUE`
#' turns it around. Which end is red is a decision about the data: with the
#' Index of Multiple Deprivation, `Q1` is the most deprived quintile.
#'
#' @param base_size Base font size.
#' @param legend Legend position.
#' @param n Number of colours.
#' @param type `"ordered"` (red to green) or `"categorical"`.
#' @param reverse Reverse the order of the colours.
#' @param ... Passed to [ggplot2::discrete_scale()].
#' @return `theme_mortensor()`: a ggplot2 theme; the scales: ggplot2 scales;
#'   `mt_palette()`: a character vector of colours.
#' @examples
#' mt_palette(5)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   x <- mt_example(deaths = "single")
#'   fc <- mt_forecast(mt_fit(lee_carter(), x), h = 10, n_draws = 100, seed = 1)
#'   e <- life_expectancy(fc, at = 65)
#'   ggplot2::ggplot(e, ggplot2::aes(period, q0.5, colour = sex)) +
#'     ggplot2::geom_line() +
#'     scale_colour_mortensor(type = "categorical") +
#'     theme_mortensor()
#' }
#' @name mt_plot_style
NULL

#' @rdname mt_plot_style
#' @export
theme_mortensor <- function(base_size = 11, legend = "bottom") {
  need_ggplot()
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(legend.position = legend, strip.text = ggplot2::element_text(size = ggplot2::rel(0.9)),
                   plot.title.position = "plot")
}

#' @rdname mt_plot_style
#' @export
mt_palette <- function(n, type = c("ordered", "categorical"), reverse = FALSE) {
  type <- match.arg(type)
  check_number(n, "n", lower = 1, integerish = TRUE)
  anchors <- if (type == "ordered") {
    c("#D7191C", "#FDAE61", "#D9B526", "#A6D96A", "#1A9641")
  } else {
    c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
  }
  cols <- if (type == "categorical" && n <= length(anchors)) anchors[seq_len(n)] else colour_ramp(anchors, n)
  if (reverse) rev(cols) else cols
}

#' @rdname mt_plot_style
#' @export
scale_colour_mortensor <- function(..., type = c("ordered", "categorical"), reverse = FALSE) {
  need_ggplot()
  type <- match.arg(type)
  ggplot2::discrete_scale("colour", palette = function(n) mt_palette(n, type, reverse), ...)
}

#' @rdname mt_plot_style
#' @export
scale_fill_mortensor <- function(..., type = c("ordered", "categorical"), reverse = FALSE) {
  need_ggplot()
  type <- match.arg(type)
  ggplot2::discrete_scale("fill", palette = function(n) mt_palette(n, type, reverse), ...)
}

# Linear interpolation between hex colours (without grDevices)
colour_ramp <- function(cols, n) {
  rgb <- vapply(cols, function(h) strtoi(substring(h, c(2, 4, 6), c(3, 5, 7)), 16L), integer(3))
  if (n == 1) return(cols[ceiling(length(cols) / 2)])
  pos <- seq(1, ncol(rgb), length.out = n)
  vapply(pos, function(p) {
    i <- floor(p)
    j <- min(i + 1, ncol(rgb))
    v <- round(rgb[, i] * (1 - (p - i)) + rgb[, j] * (p - i))
    sprintf("#%02X%02X%02X", v[1], v[2], v[3])
  }, "")
}
