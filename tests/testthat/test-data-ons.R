test_that("wide tables are read with one value column per level", {
  d <- data.frame(Year = 2020, Region = rep(c("A", "B"), each = 2), Age = rep(c("85-89", "90+"), 2),
                  Female = 1:4, Male = c("9", ".", "11", "12"))
  p <- data.frame(Year = 2020, Region = rep(c("A", "B"), each = 2), Age = rep(c("85-89", "90+"), 2),
                  Population = 100)
  expect_warning(x <- mt_data(deaths = d, population = p, by = c("Region", "sex"), age = "Age", period = "Year",
                              value = list(deaths = c(female = "Female", male = "Male"), population = "Population")),
                 class = "mt_warning_missing")
  st <- mt_strata(x)
  expect_identical(nrow(st), 4L)
  expect_setequal(st$sex, c("female", "male"))
  expect_equal(sum(mt_obs(x)$value), sum(1:4) + sum(c(9, 11:12)))
  expect_identical(nrow(mt_obs(x, "population")), 8L)
  expect_error(mt_data(deaths = d, by = c("Region", "sex", "cause"), age = "Age", period = "Year",
                       value = list(deaths = c("Female", "Male"))), class = "mt_error_obs")
  expect_error(mt_data(deaths = d, by = c("Region", "sex"), age = "Age", period = "Year",
                       value = list(deaths = c("Female", "Other"))), class = "mt_error_obs")
})

test_that("the ONS data read into consistent strata with the general functions", {
  expect_identical(dim(ons_data), c(7600L, 36L))
  x <- mt_data(deaths = ons_data, population = ons_data, by = c("sex", "imd_decile"), age = "age_group",
               value = c(deaths = "all_causes", population = "population"))
  x <- mt_as_exposure(x)
  expect_equal(sum(mt_obs(x, "deaths")$value), sum(ons_data$all_causes))
  expect_equal(sum(mt_obs(x, "exposure")$value), sum(ons_data$population))
  expect_identical(mt_grid(x)$age$labels[1:3], c("0", "1-4", "5-9"))
  q <- mt_group_levels(x, "imd_decile", function(d) paste0("Q", ceiling(as.numeric(d) / 2)), name = "imd_quintile")
  expect_identical(q$by, c("sex", "imd_quintile"))
  expect_identical(unique(mt_strata(q)$imd_quintile), paste0("Q", 1:5))
  expect_equal(sum(mt_obs(q, "exposure")$value), sum(ons_data$population))
  causes <- names(ons_data)[-(1:6)]
  y <- mt_data(deaths = ons_data, population = ons_data, by = c("sex", "imd_decile", "cause"), age = "age_group",
               value = list(deaths = causes, population = "population"))
  y <- mt_subset(mt_as_exposure(y), ages = c(50, Inf), periods = c(2015, 2019))
  y <- mt_group_levels(y, "cause", ons_cause_groups, other = "rest", total = "all_causes", additive = "rates")
  expect_identical(mt_grid(y)$age$labels[1], "50-54")
  d <- mt_obs(y)
  expect_equal(sum(d[cause == "all_causes"]$value), sum(d[cause != "all_causes"]$value))
  expect_true(all(d[cause == "rest"]$value >= 0))
  v <- mt_validate(y, hierarchy = mt_hierarchy(cause = "all_causes", additive = "rates"))
  expect_true(all(v$severity == "info"))
  e <- mt_obs(y, "exposure")
  expect_equal(e[cause == "rest"]$value, e[cause == "all_causes"]$value)
  expect_warning(mt_group_levels(y, "cause", list(big = c("ihd", "rest", "lung_cancer", "stroke", "dementia",
                                                          "other_cancers")), other = "small",
                                 total = "all_causes", additive = "rates"), NA)
})

test_that("grouped ages are the default grid and regroup to a table's intervals", {
  d <- data.frame(year = 2020, age = c("0-4", "5-9", "10+"), deaths = 1:3)
  p <- data.frame(year = 2020, age = c(0:14, "15+"), exposure = 100)
  p$age[16] <- "15+"
  x <- mt_data(deaths = d, exposure = p)
  expect_identical(length(mt_grid(x)$age), 16L)
  expect_error(mt_fit(rwd(), x), class = "mt_error_model_data")
  y <- mt_regroup(x, age = "deaths")
  expect_identical(mt_grid(y)$age$labels, c("0-4", "5-9", "10+"))
  expect_equal(mt_obs(y, "exposure")$value, c(500, 500, 600))
})
