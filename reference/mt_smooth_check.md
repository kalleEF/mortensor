# Check and plot the result of smoothing

`mt_smooth_check()` compares smoothed with raw rates per stratum:

## Usage

``` r
mt_smooth_check(data, smoothed)

mt_plot_smooth(data, smoothed, ages = NULL, strata = NULL, colour = NULL)
```

## Arguments

- data:

  The [mt_data](mt_data.md) object before smoothing.

- smoothed:

  The result of [`mt_smooth()`](mt_smooth.md) (or of
  [`mt_prep()`](mt_prep.md) steps).

- ages:

  Age labels to show (default: four ages across the range).

- strata:

  Stratum labels to show (default: all).

- colour:

  Stratum variable that gives the colours.

## Value

`mt_smooth_check()`: a `data.table` with the stratum columns and the
measures above. `mt_plot_smooth()`: a ggplot object.

## Details

- `mad_log`: median absolute difference of log rates;

- `dispersion`: mean squared Pearson residual of the raw deaths around
  the smoothed deaths; about 1 when smoothing removed Poisson noise
  only, above 1 when the raw rates vary more (for example year effects),
  well below 1 when the smoothing follows the noise;

- `trend_ratio`: variance of year-to-year changes of the smoothed log
  rates relative to that of the raw log rates (mean over ages); small
  values mean that smoothing removed most of the annual variation;

- `correlation`: correlation of raw and smoothed log rates.

`mt_plot_smooth()` shows raw rates as points and smoothed rates as lines
over time, one panel per age.

## Examples

``` r
x <- mt_example(deaths = "single", years = 2000:2015)
s <- mt_smooth(x, along = "period", method = "whittaker")
mt_smooth_check(x, s)
#>    stratum    sex  label    mad_log dispersion  trend_ratio correlation
#>      <int> <char> <char>      <num>      <num>        <num>       <num>
#> 1:       1 female female 0.02210234  0.9265103 0.0003089649   0.9994547
#> 2:       2   male   male 0.02110521  0.8958134 0.0005426953   0.9996224
if (requireNamespace("ggplot2", quietly = TRUE)) mt_plot_smooth(x, s)
```
