# MCMC settings for apc_bayes()

Presets: `"fast"` (5,000 iterations, burn-in 1,000, thinning 5; suitable
for backtests), `"final"` (bamp's automatic choice of chain length), and
`"test"` (1,000 iterations, for examples and unit tests only).

## Usage

``` r
mcmc_control(
  preset = c("fast", "final", "test"),
  iterations = NULL,
  burn_in = NULL,
  thin = NULL,
  tuning = NULL,
  parallel = TRUE
)
```

## Arguments

- preset:

  `"fast"`, `"final"` or `"test"`.

- iterations, burn_in, thin:

  Override the preset.

- tuning:

  Tuning iterations.

- parallel:

  Run the four chains in parallel. bamp runs them in a socket cluster
  without seeded random number streams, so results differ slightly
  between runs even when a seed is set.

## Value

A list.

## Examples

``` r
mcmc_control("fast")
#> $options
#> $options$number_of_iterations
#> [1] 5000
#> 
#> $options$burn_in
#> [1] 1000
#> 
#> $options$step
#> [1] 5
#> 
#> $options$tuning
#> [1] 500
#> 
#> 
#> $parallel
#> [1] TRUE
#> 
#> $preset
#> [1] "fast"
#> 
```
