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
  comparisons = NULL
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

  Numbers of the comparisons to show (rows of
  `summary(decomposition, by = NULL)`; default: all, as panels).

## Value

A ggplot object.

## Examples

``` r
x <- mt_example(deaths = "single", years = 2000:2015)
mt_plot_decomposition(mt_decompose(x, periods = c(2000, 2015), at = 30))
```
