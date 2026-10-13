# Decomposition of differences by age, cause and subgroup

ons_small <- function(by_cause = FALSE) {
  if (by_cause) {
    causes <- names(mortensor::ons_data)[-(1:7)]
    x <- mt_data(deaths = mortensor::ons_data, population = mortensor::ons_data, by = c("sex", "imd_decile", "cause"),
                 age = "age_group", value = list(deaths = c(all = "all_causes", causes), population = "population"))
    x <- mt_as_exposure(x, type = "midyear")
    x <- mt_subset(x, ages = c(50, Inf), periods = c(2017, 2019))
    x <- mt_group_levels(x, "imd_decile", function(d) paste0("Q", ceiling(as.numeric(d) / 2)), name = "quintile")
    return(mt_group_levels(x, "cause", mortensor::ons_cause_groups, other = "rest", total = "all", additive = "rates"))
  }
  x <- mt_data(deaths = mortensor::ons_data, population = mortensor::ons_data, by = c("sex", "imd_decile"),
               age = "age_group", value = c(deaths = "all_causes", population = "population"))
  x <- mt_as_exposure(mt_subset(x, periods = c(2001, 2019)), type = "midyear")
  mt_group_levels(x, "imd_decile", function(d) paste0("Q", ceiling(as.numeric(d) / 2)), name = "quintile")
}

test_that("decompositions of life expectancy add up and agree with their definitions", {
  x <- mt_example(deaths = "single", years = 2000:2015)
  for (m in c("arriaga", "stepwise", "horiuchi")) {
    d <- mt_decompose(x, periods = c(2000, 2015), at = 30, method = m)
    tot <- summary(d, by = NULL)
    e <- life_expectancy(mt_array(x, "deaths") / mt_array(x, "exposure"), mt_grid(x)$age, at = 30)
    expect_equal(tot$difference, unname(e["2015", ] - e["2000", ]), tolerance = 1e-10)
    expect_lt(max(tot$max_residual), if (m == "horiuchi") 1e-4 else 1e-10)
  }
  # forward Arriaga is stepwise replacement from the oldest age
  m <- mt_array(x, "deaths") / mt_array(x, "exposure")
  age <- mt_grid(x)$age
  rows <- match(30, age$lower):length(age)
  ag <- list(lower = age$lower[rows], upper = age$upper[rows])
  m1 <- matrix(m[rows, "2000", 1])
  m2 <- matrix(m[rows, "2015", 1])
  f <- function(M) mortensor:::lt_columns(M, ag)$ex[1, ]
  sw <- mortensor:::stepwise_core(f, m1, m2, rev(seq_along(rows)), symmetric = FALSE)
  expect_equal(mortensor:::arriaga_dir(m1, m2, ag), sw, tolerance = 1e-12)
  # ages below `at` are not part of the decomposition
  d <- mt_decompose(x, periods = c(2000, 2015), at = 30)
  expect_setequal(unique(d$table$age), age$labels[rows])
})

test_that("decompositions match DemoDecomp and LEdecomp", {
  skip_if_not_installed("DemoDecomp")
  x <- mt_example(deaths = "single", years = 2000:2015)
  m <- mt_array(x, "deaths") / mt_array(x, "exposure")
  age <- mt_grid(x)$age
  rows <- match(30, age$lower):length(age)
  ag <- list(lower = age$lower[rows], upper = age$upper[rows])
  e30 <- function(v) mortensor:::lt_columns(matrix(v, ncol = 1), ag)$ex[1, 1]
  m1 <- m[rows, "2000", "female"]
  m2 <- m[rows, "2015", "female"]
  get <- function(method) {
    d <- mt_decompose(x, periods = c(2000, 2015), at = 30, method = method)
    summary(d, by = "age")[sex == "female"]$mean
  }
  expect_equal(get("horiuchi"), DemoDecomp::horiuchi(e30, m1, m2, N = 20), tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(get("stepwise"), DemoDecomp::stepwise_replacement(e30, m1, m2, symmetrical = TRUE, direction = "up"),
               tolerance = 1e-10, ignore_attr = TRUE)
  skip_if_not_installed("LEdecomp")
  expect_equal(get("arriaga"), LEdecomp::arriaga_sym(m1, m2, age = age$lower[rows], nx = rep(1, length(rows)),
                                                     sex1 = "f", sex2 = "f", closeout = TRUE),
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("the Kitagawa step reproduces Riffe et al. (2024), Table 1", {
  # rounded inputs of the table: shares and e(35) of men (1) and women (2)
  k <- mortensor:::kitagawa_mean(c(0.10, 0.52, 0.38), c(0.09, 0.41, 0.50), c(42.47, 44.52, 46.92),
                                 c(48.21, 50.29, 51.30))
  expect_equal(as.vector(k$value), c(0.54, 2.69, 1.92), tolerance = 0.02)
  expect_equal(k$composition, 0.23, tolerance = 0.05)
  expect_equal(sum(k$value) + k$composition, sum(c(0.09, 0.41, 0.50) * c(48.21, 50.29, 51.30)) -
                 sum(c(0.10, 0.52, 0.38) * c(42.47, 44.52, 46.92)), tolerance = 1e-12)
})

test_that("causes add up to the all-cause decomposition", {
  x <- ons_small(by_cause = TRUE)
  h <- mt_hierarchy(cause = "all", additive = "rates")
  dc <- mt_decompose(x, var = "quintile", ref = "Q1", periods = 2019, cause = "cause", hierarchy = h)
  expect_lt(max(summary(dc, by = NULL)$max_residual), 1e-10)
  xa <- mt_subset(x, strata = cause == "all")
  da <- mt_decompose(xa, var = "quintile", ref = "Q1", periods = 2019)
  sc <- summary(dc, by = "age")
  sa <- summary(da, by = "age")
  expect_equal(sc$mean, sa$mean[match(paste(sc$sex, sc$to, sc$age), paste(sa$sex, sa$to, sa$age))],
               tolerance = 1e-8)
  expect_setequal(dc$cause_levels, setdiff(unique(x$strata$cause), "all"))
  dh <- mt_decompose(x, var = "quintile", ref = "Q1", periods = 2019, cause = "cause", hierarchy = h,
                     method = "horiuchi")
  expect_lt(max(summary(dh, by = NULL)$max_residual), 1e-3)
})

test_that("subgroups: radix (Arriaga-Kitagawa) and pooled composition", {
  x <- ons_small()
  dr <- mt_decompose(x, periods = c(2001, 2019), groups = "quintile")
  tot <- summary(dr, by = NULL)
  expect_lt(max(tot$max_residual), 1e-10)
  # total = radix-weighted mean of the quintiles' life expectancy
  E <- mt_array(x, "exposure")
  e <- life_expectancy(mt_array(x, "deaths") / E, mt_grid(x)$age)
  st <- x$strata
  w <- E[1, "2001", st$label[st$sex == "female"]]
  expect_equal(tot[sex == "female"]$from_value, sum(w / sum(w) * e["2001", st$label[st$sex == "female"]]),
               tolerance = 1e-10)
  expect_equal(nrow(summary(dr, by = "component")[component == "composition"]), 2)
  # pooled: total = life expectancy of the summed deaths and exposures
  dp <- mt_decompose(x, periods = c(2001, 2019), groups = "quintile", composition = "pooled", n_steps = 10)
  pooled <- mt_group_levels(x, "quintile", list(all = paste0("Q", 1:5)))
  ep <- life_expectancy(mt_array(pooled, "deaths") / mt_array(pooled, "exposure"), mt_grid(x)$age)
  expect_equal(summary(dp, by = NULL)[sex == "female"]$from_value, unname(ep["2001", "female.all"]),
               tolerance = 1e-8)
  expect_lt(max(summary(dp, by = NULL)$max_residual), 1e-3)
  expect_error(mt_decompose(x, groups = "quintile", composition = "pooled", method = "arriaga"),
               class = "mt_error_argument")
  expect_error(mt_decompose(x, groups = "quintile", quantity = "asr"), class = "mt_error_argument")
})

test_that("other quantities, Kitagawa and custom functions", {
  x <- ons_small()
  dk <- mt_decompose(x, periods = c(2001, 2019), quantity = "crude_rate")
  expect_lt(max(summary(dk, by = NULL)$max_residual), 1e-8)
  expect_setequal(unique(dk$table$component), c("rate", "composition"))
  da <- mt_decompose(x, periods = c(2001, 2019), quantity = "asr", method = "stepwise")
  a <- asr(mt_observed(x, n_draws = 1))
  expect_equal(summary(da, by = NULL)[sex == "female" & quintile == "Q1"]$from_value,
               a[label == "female.Q1" & period == 2001]$mean, tolerance = 1e-10)
  dp <- mt_decompose(x, periods = c(2001, 2019), quantity = "prob_dying", at = 30, to = 75)
  expect_true(all(dp$table$age %in% c("30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69",
                                      "70-74")))
  cu <- mt_decompose(x, periods = c(2001, 2019), quantity = function(m, age) colSums(m), method = "stepwise")
  expect_lt(max(summary(cu, by = NULL)$max_residual), 1e-12)
  expect_error(mt_decompose(x, quantity = "asr", method = "arriaga"), class = "mt_error_argument")
  expect_error(mt_decompose(x, quantity = "crude_rate", method = "horiuchi"), class = "mt_error_argument")
  expect_error(mt_decompose(x, var = "region", ref = "a"), class = "mt_error_argument")
  expect_error(mt_decompose(x, var = "quintile", ref = "Q9"), class = "mt_error_argument")
})

test_that("forecasts are decomposed draw by draw, alone or with observed years", {
  x <- mt_example(deaths = "single", years = 2000:2015)
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 10, n_draws = 60, seed = 1)
  d <- mt_decompose(fc, var = "sex", ref = "male", at = 30)
  expect_equal(ncol(d$paths), 60)
  s <- summary(d, by = NULL)
  expect_lt(s$max_residual, 1e-10)
  e <- life_expectancy(fc, at = 30, paths = TRUE)
  expect_equal(s$difference, mean(e["2025", "female", ] - e["2025", "male", ]), tolerance = 1e-10)
  sa <- summary(d, by = "age")
  expect_true(all(sa$q0.025 <= sa$q0.975) && any(sa$q0.025 < sa$q0.975))
  # observed 2015 against the forecast of 2025
  d2 <- mt_decompose(fc, periods = c(2015, 2025), at = 30)
  e15 <- life_expectancy(fc$observed, fc$age, at = 30)
  expect_equal(summary(d2, by = NULL)$from_value, unname(e15["2015", ]), tolerance = 1e-10)
  # two forecasts compared
  fr <- mt_forecast(mt_fit(rwd(), x), h = 10, n_draws = 60, seed = 2)
  d3 <- mt_decompose(fc, fr, periods = 2020, at = 30)
  expect_equal(unique(d3$meta$to), fr$model)
  # observed data with Poisson draws
  d4 <- mt_decompose(x, periods = c(2000, 2015), at = 30, n_draws = 50)
  expect_equal(ncol(d4$paths), 50)
  skip_if_not_installed("ggplot2")
  expect_s3_class(mt_plot_decomposition(d), "ggplot")
  xc <- ons_small(by_cause = TRUE)
  dc <- mt_decompose(xc, var = "quintile", ref = "Q1", periods = 2019, cause = "cause",
                     hierarchy = mt_hierarchy(cause = "all", additive = "rates"))
  expect_s3_class(mt_plot_decomposition(dc, by = "age", fill = "cause", comparisons = 4), "ggplot")
  expect_s3_class(mt_plot_decomposition(dc, by = "cause", comparisons = 4), "ggplot")
  dr <- mt_decompose(ons_small(), periods = c(2001, 2019), groups = "quintile")
  expect_s3_class(mt_plot_decomposition(dr, by = "group"), "ggplot")
  expect_error(mt_plot_decomposition(d, by = "cause"), class = "mt_error_argument")
})
