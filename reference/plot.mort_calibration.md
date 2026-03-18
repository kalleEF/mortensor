# Plot a mort_calibration object

Plot a mort_calibration object

## Usage

``` r
# S3 method for class 'mort_calibration'
plot(x, type = c("comparison", "ranks"), ...)
```

## Arguments

- x:

  A `mort_calibration` object.

- type:

  One of `"comparison"` (bar chart of MASE by model) or `"ranks"`
  (Tucker3 rank selection heatmap, if available).

- ...:

  Ignored.

## Value

A `ggplot` object (invisibly).
