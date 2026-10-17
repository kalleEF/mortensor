# Plot a decomposition

Bars of the mean contributions by age, cause or subgroup, stacked by a
second variable, with the 95% interval of the sum per bar when the
decomposition has more than one draw. The composition effect of
decompositions with subgroups is shown as its own bar.

## Usage

``` r
mt_plot_decomposition(
  decomposition,
  by = c("age", "cause", "group"),
  fill = NULL,
  comparisons = NULL,
  change = NULL
)
```

## Arguments

- decomposition:

  An [`mt_decompose()`](mt_decompose.md) result.

- by:

  Variable on the axis: `"age"`, `"cause"` or `"group"`.

- fill:

  Variable that stacks the bars (`"cause"`, `"group"`), or `NULL`.

- comparisons:

  Comparisons to show: numbers (column `comparison` of
  `summary(decomposition, by = NULL)`) or a named list of values of its
  columns, for example `list(to = "Q5")`; default: all, as panels.

- change:

  Two periods: plot the change between them (see
  [`summary.mt_decomposition()`](mt_decompose.md)).

## Value

A ggplot object.

## Examples

``` r
x <- mt_example(deaths = "single", years = 2000:2015)
mt_plot_decomposition(mt_decompose(x, periods = c(2000, 2015), at = 30))
```
