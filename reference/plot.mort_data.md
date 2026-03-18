# Plot a mort_data object

Plot a mort_data object

## Usage

``` r
# S3 method for class 'mort_data'
plot(
  x,
  type = c("heatmap", "age_profile", "trend"),
  group = 1L,
  years = NULL,
  ages = NULL,
  ...
)
```

## Arguments

- x:

  A `mort_data` object.

- type:

  One of `"heatmap"` (Lexis surface), `"age_profile"` (log-mx vs age for
  selected years), or `"trend"` (log-mx vs year for selected ages).

- group:

  Integer index of the group to plot (default 1).

- years:

  Integer vector of years to highlight (for `"age_profile"`). Defaults
  to 4 evenly-spaced years.

- ages:

  Integer vector of age starts to highlight (for `"trend"`). Defaults to
  every other age group.

- ...:

  Ignored.

## Value

A `ggplot` object (invisibly).
