long_fixture <- function() {
  d <- expand.grid(sex = c("f", "m"), year = 2001:2003, age = c("30-34", "35-39", "40+"),
                   stringsAsFactors = FALSE)
  d$deaths <- seq_len(nrow(d))
  d
}

test_that("mt_data ingests labels, strata and computes a default grid", {
  x <- mt_data(deaths = long_fixture(), by = "sex")
  expect_s3_class(x, "mt_data")
  expect_equal(mt_strata(x)$label, c("f", "m"))
  g <- mt_grid(x)
  expect_equal(g$age$lower[1], 30)
  expect_true(is_open(g$age))
  expect_equal(tail(g$age$labels, 1), "40+")
  expect_length(g$period, 3)
  o <- mt_obs(x, "deaths")
  expect_equal(nrow(o), 18)
  expect_equal(sum(o$value), sum(long_fixture()$deaths))
})

test_that("malformed inputs give classed errors", {
  d <- long_fixture()
  expect_error(mt_data(deaths = rbind(d, d[1, ]), by = "sex"), class = "mt_error_obs")
  d2 <- d
  d2$age[d2$age == "35-39"] <- "33-39"
  expect_error(mt_data(deaths = d2, by = "sex"), class = "mt_error_obs")
  d3 <- d
  d3$deaths[1] <- -1
  expect_error(mt_data(deaths = d3, by = "sex"), class = "mt_error_obs")
  expect_error(mt_data(deaths = d, by = "region"), class = "mt_error_obs")
  d4 <- d
  d4$deaths[2] <- 1.5
  expect_warning(mt_data(deaths = d4, by = "sex"), class = "mt_warning_noninteger")
  d5 <- d
  d5$deaths[3] <- NA
  expect_warning(mt_data(deaths = d5, by = "sex"), class = "mt_warning_missing")
})

test_that("tables of different resolution share the grid; coarse cells are detected", {
  x <- mt_example(years = 2001:2003)
  g <- mt_grid(x)
  expect_equal(length(g$age), 71)
  expect_error(mt_array(x, "deaths"), class = "mt_error_resolution")
  a <- mt_array(x, "deaths", coarse = "na")
  expect_true(all(is.na(a)))
  e <- mt_array(x, "exposure")
  expect_false(anyNA(e))
  expect_equal(dimnames(e)$stratum, c("female", "male"))
})

test_that("grid must be compatible with the observations", {
  x <- mt_example(years = 2001:2002)
  expect_error(mt_set_grid(x, age = age_groups(seq(30, 100, 10))), class = "mt_error_grid")
  y <- mt_set_grid(x, age = single_ages(30, 105))
  expect_equal(length(mt_grid(y)$age), 76)
})

test_that("mt_regroup preserves totals and rejects straddling cells", {
  x <- mt_example(years = 2001:2004)
  r <- mt_regroup(x, age = c(30, 50, 70, 90), period = c(2001, 2003))
  for (m in c("deaths", "exposure")) {
    expect_equal(sum(mt_obs(r, m)$value), sum(mt_obs(x, m)$value), tolerance = 1e-10)
  }
  expect_equal(mt_grid(r)$age$labels, c("30-49", "50-69", "70-89", "90+"))
  expect_equal(mt_grid(r)$period$labels, c("2001-2002", "2003-2004"))
  expect_error(mt_regroup(x, age = c(30, 52, 90)), class = "mt_error_regroup")
})

test_that("regrouping drops partly covered cells and handles an open interval that changes over time", {
  d <- data.frame(year = rep(2001:2002, c(3, 4)), age_lo = c(80, 85, 90, 80, 85, 90, 95),
                  age_hi = c(85, 90, Inf, 85, 90, 95, Inf), deaths = 1:7)
  x <- mt_data(deaths = d, age = c("age_lo", "age_hi"))
  expect_equal(tail(mt_grid(x)$age$labels, 1), "95+")
  r <- mt_regroup(x, age = c(80, 90))
  o <- mt_obs(r, "deaths")
  expect_equal(o[period_lo == 2001 & age_lo == 90]$value, 3)
  expect_equal(o[period_lo == 2002 & age_lo == 90]$value, 13)
  expect_warning(mt_regroup(x, age = c(80, 85), period = c(2000, 2002)), class = "mt_warning_regroup_partial")
})

test_that("operator rows sum the grid cells inside each observed cell", {
  x <- mt_example(years = 2001:2002)
  C <- mt_operator(x, "deaths", 1)
  e <- mt_array(x, "exposure")[, , 1]
  agg <- as.numeric(C %*% as.vector(e))
  expect_equal(nrow(C), 26)
  expect_equal(sum(agg), sum(e))
  expect_equal(range(Matrix::rowSums(C)), c(5, 11))
})

test_that("mt_window keeps only cells inside the window", {
  x <- mt_example(years = 2001:2010)
  w <- mt_window(x, from = 2003, to = 2006)
  expect_equal(range(mt_obs(w, "deaths")$period_lo), c(2003, 2006))
  expect_equal(mt_grid(w)$period$labels, as.character(2003:2006))
})

test_that("mt_validate flags impossible values and doubled totals", {
  d <- data.frame(cause = rep(c("all", "a", "b"), each = 2), year = 2001, age = c(70, 71),
                  deaths = c(20, 30, 4, 5, 6, 10), exposure = c(100, 0, 100, 100, 100, 100))
  x <- mt_data(deaths = d, exposure = d, by = "cause")
  v <- mt_validate(x, hierarchy = list(var = "cause", total = "all"))
  expect_s3_class(v, "mt_validation")
  expect_true("deaths without exposure" %in% v$check)
  expect_true(any(grepl("twice the sum of parts", v$detail)))
  expect_output(print(v), "finding")
})

test_that("exposures from January 1 populations average consecutive years", {
  p <- data.frame(year = 2001:2003, age = 70, population = c(1000, 990, 980))
  x <- mt_as_exposure(mt_add_population(mt_data(), p), type = "jan1")
  e <- mt_obs(x, "exposure")
  expect_equal(e$value, c(995, 985))
  expect_equal(unique(e$provenance), "population_jan1")
})

test_that("adapters convert HMD-based objects", {
  skip_if_not_installed("StMoMo")
  ew <- as_mt_data(StMoMo::EWMaleData)
  expect_equal(length(mt_grid(ew)$age), 101)
  expect_equal(sum(mt_obs(ew, "deaths")$value), sum(StMoMo::EWMaleData$Dxt))
  skip_if_not_installed("demography")
  fr <- suppressWarnings(as_mt_data(demography::fr.mort, years = 2000:2001))
  expect_equal(tail(mt_grid(fr)$age$labels, 1), "110+")
  expect_equal(mt_strata(fr)$sex, c("female", "male"))
})

test_that("printing and conversion work", {
  x <- mt_example(years = 2001:2002)
  expect_output(print(x), "mt_data")
  dt <- data.table::as.data.table(x)
  expect_setequal(unique(dt$measure), c("deaths", "exposure"))
  expect_true(all(c("sex", "label", "age_lo", "value") %in% names(dt)))
})
