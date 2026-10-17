# Plot gaps between strata over time

Draws the result of [`mt_contrast()`](mt_contrast.md) over time:
forecast medians (dashed) with intervals and, if given, the observed
gaps (points). Several forecasts (models, scenarios) can be compared as
a named list, one colour each; otherwise the levels of the contrast get
the colours. Other stratum variables (for example sex) give the panels.

## Usage

``` r
mt_plot_contrast(contrast = NULL, observed = NULL, levels = NULL)
```

## Arguments

- contrast:

  The result of [`mt_contrast()`](mt_contrast.md) for a forecast, or a
  named list of such results.

- observed:

  The result of [`mt_contrast()`](mt_contrast.md) for observed data, or
  `NULL`. Without `contrast`, the observed gaps are drawn with their
  Poisson intervals.

- levels:

  Levels of the compared variable to show (default all).

## Value

A ggplot object.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  x <- mt_example(deaths = "single", years = 1995:2012)
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 10, n_draws = 100, seed = 1)
  mt_plot_contrast(mt_contrast(fc, var = "sex", ref = "male", at = 30),
                   observed = mt_contrast(x, var = "sex", ref = "male", at = 30))
}
```
