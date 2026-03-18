# Run the coherent multi-disease forecast pipeline (Tucker4-ILR)

Jointly models cause-of-death shares via a 4D Tucker decomposition on
ILR-transformed shares, ensuring compositions sum to 1 by construction.

## Usage

``` r
run_coherent_forecast(disease_configs, all_cause_config, shared)
```

## Arguments

- disease_configs:

  Named list of `mort_config` objects (one per cause).

- all_cause_config:

  `mort_config` for all-causes.

- shared:

  Named list of shared parameters.

## Value

Invisibly, a list with all coherent results.
