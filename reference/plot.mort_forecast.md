# Plot a mort_forecast object

Plot a mort_forecast object

## Usage

``` r
# S3 method for class 'mort_forecast'
plot(
  x,
  type = c("trends", "trends_by_age", "age", "pi", "residuals", "smoothing",
    "components", "model_comparison", "disagg"),
  group = NULL,
  facet_by = NULL,
  color_by = NULL,
  ...
)
```

## Arguments

- x:

  A `mort_forecast` object.

- type:

  One of `"trends"` (age-averaged trends with PIs), `"trends_by_age"`
  (trends faceted by age group), `"age"` (age profiles for forecast
  years), `"pi"` (PI fan chart), `"residuals"` (residual diagnostics),
  `"smoothing"` (observed vs smoothed), `"components"` (Tucker3 B-mode
  temporal components), `"model_comparison"` (bar chart from
  calibration), or `"disagg"` (single-year disaggregated profiles).

- group:

  Integer index to show a single group, or `NULL` (default) for all
  groups.

- facet_by:

  Character vector of column names from the groups data.frame to use for
  faceting. If `NULL` (default), auto-detected from group structure:
  when groups have 2+ columns, the first columns are used for faceting.

- color_by:

  Character column name from the groups data.frame to map to the colour
  aesthetic. If `NULL` (default), auto-detected: the last column of the
  groups data.frame is used.

- ...:

  Additional arguments passed to ggplot.

## Value

A `ggplot` object (invisibly).
