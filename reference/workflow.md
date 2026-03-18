# Step-by-step mortality forecasting workflow

These functions provide a modular, interactive workflow for mortality
forecasting as an alternative to the monolithic
[`run_forecast()`](https://kalleEF.github.io/mortensor/reference/run_forecast.md):

## Details

1.  [`read_mortality()`](https://kalleEF.github.io/mortensor/reference/read_mortality.md)
    or
    [`read_mortality_from_arrays()`](https://kalleEF.github.io/mortensor/reference/read_mortality_from_arrays.md)
    – load data

2.  [`smooth_mortality()`](https://kalleEF.github.io/mortensor/reference/smooth_mortality.md)
    – smooth the data (already exists)

3.  [`calibrate()`](https://kalleEF.github.io/mortensor/reference/calibrate.md)
    – calibrate and compare models

4.  [`forecast.mort_calibration()`](https://kalleEF.github.io/mortensor/reference/forecast.mort_calibration.md)
    – produce forecasts

5.  [`add_pi()`](https://kalleEF.github.io/mortensor/reference/add_pi.md)
    – add prediction intervals

6.  [`disaggregate()`](https://kalleEF.github.io/mortensor/reference/disaggregate.md)
    – disaggregate to single-year ages
