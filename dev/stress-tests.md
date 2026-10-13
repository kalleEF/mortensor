# Stress tests with simulated defective data

`dev/stress-tests.R` simulates deaths and populations for three regions
(200,000, 100,000 and 5,000 inhabitants) by sex, 5-year age groups and
2001-2020, then introduces one realistic defect per scenario and runs the
usual pipeline (`mt_data()`, `mt_as_exposure()`, `mt_fit(lee_carter())`,
`mt_forecast()`, `life_expectancy()`). It records whether mortensor stops
with an error, warns or continues silently, and the message. The
scenarios that remain are part of the unit tests
(`tests/testthat/test-robustness.R`).

## Outcomes before and after the changes of this round

| Scenario | Before | After |
|---|---|---|
| Rows "All ages" in age column | error "Cannot parse age label(s)" | rows left out with a warning |
| Level "Persons" next to female and male | silent (counted twice in aggregates) | warning in `mt_validate()` and `mt_aggregate_strata()` |
| Sex "female"/"male" in deaths, "F"/"M" in populations | generic error at fit ("must be complete") | warning in `mt_data()` naming the levels; error at fit names the strata without deaths or exposures |
| Trailing spaces, upper/lower case in labels | generic error | spaces trimmed; case differences named in a warning |
| Missing sex in some rows | generic error | rows left out with a warning |
| Year and age as factors | correct | correct |
| Financial years "2001/02" | error "Cannot parse age label(s)" | read as years starting in 2001 |
| Age labels "1 to 4", "85 years and over", "90 & over", "Age 85+", "1 bis 4 Jahre" | 8 of 14 variants failed | all parsed |
| Single ages with the last age meaning 90+ | silent | warning suggesting `age_open` |
| Thousands separators "12,345" | values lost (warning), fit error | read correctly |
| Decimal commas "1234,5" | values lost, generic error | warning suggests `decimal = ","` |
| Duplicated block of rows | error without cause | error names the unused columns and repeated rows |
| Missing year (all strata) | generic error | error names the year and suggests `mt_subset(periods = )` |
| Year missing in one stratum | generic error | error names the cells |
| Populations missing for a region | generic error | error names the strata; `mt_validate()` warns |
| Deaths to 90+, populations to 85+ | message with wrong grammar | message suggests `mt_regroup(age = "common")` |
| Age groups changing over time (85+ then 85-89, 90+) | regroup to the deaths' groups failed | `mt_regroup(age = "common")` gives 85+ |
| Subtotal 85+ next to 85-89 and 90+ | "duplicate cells" | "overlapping age groups ... subtotals" |
| Zero population with deaths | "Exposures must be positive" | error locates the cells |
| Deaths larger than population | silent | warning locates the cells |
| Negative death count | "Negative deaths values." | count and location |
| Counts "<5" suppressed, no totals | `mt_fill_suppressed()` filled the cells with the populations, silently | error without totals or bounds; with `lower`, `upper` the cells get bounded expected deaths |
| Region of 5,000 people (53% zero cells), Lee-Carter | Poisson fit diverged, stratum silently missing; non-finite forecasts | zero-death age groups get a constant rate, sparse rows a ridge penalty; one warning |
| Same region, random walk | rates up to 1e3 in 64% of draws, silently | warning about implausible rates; `rwd(observation_noise = TRUE)` pools the level variance and stays plausible |
| 400 districts x 2 sexes x 18 ages x 20 years | 4.2 s read, 2.4 s grouping | 2.9 s read, 0.4 s grouping |

## Remaining limits

* Lee-Carter and the plain random walk remain unsuitable for very small
  populations; the warnings point to pooling (`li_lee(common_by = )`),
  wider age groups and smoothing.
* A numeric top age is only recognised as open when it holds clearly more
  than the age below it; otherwise `age_open` must be given.
* "1,234" with a decimal comma is ambiguous and read as 1234 unless
  `decimal = ","` is set.
