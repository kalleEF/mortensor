test_that("ONS Excel tables are read with deciles combined to quintiles", {
  skip_if_not_installed("readxl")
  f <- test_path("fixtures", "ons-synthetic.xlsx")
  raw <- as.data.frame(readxl::read_excel(f, sheet = "1", skip = 2))
  expect_warning(x <- mt_read_ons(f, f), class = "mt_warning_ons_fill")
  st <- mt_strata(x)
  expect_identical(sort(unique(st$quintile)), paste0("Q", 1:5))
  expect_identical(sort(unique(st$sex)), c("female", "male"))
  expect_equal(sum(mt_obs(x, "deaths")$value), sum(raw[["All causes"]]))
  e <- mt_obs(x, "exposure")
  expect_equal(e[period_lo == 2021]$value, e[period_lo == 2020]$value)
  expect_identical(mt_grid(x)$age$labels, c("85-89", "90+"))
  expect_identical(dim(mt_array(x, "deaths")), c(2L, 3L, 10L))
})

test_that("cause levels sum the listed columns and the rest closes the total", {
  skip_if_not_installed("readxl")
  f <- test_path("fixtures", "ons-synthetic.xlsx")
  x <- suppressWarnings(mt_read_ons(f, f, causes = list(circ = "Stroke", dem = "dementia", other = "rest"),
                                    deprivation = "decile", years = 2019:2020))
  st <- mt_strata(x)
  expect_setequal(unique(st$cause), c("all", "circ", "dem", "other"))
  expect_setequal(unique(st$decile), paste0("D", 1:10))
  d <- mt_obs(x, "deaths")
  tot <- d[, list(all = sum(value[cause == "all"]), parts = sum(value[cause != "all"])),
           by = list(sex, decile, age_lo, period_lo)]
  expect_equal(tot$all, tot$parts)
  expect_error(mt_read_ons(f, f, causes = list(a = "nonexistent")), class = "mt_error_ons")
  expect_error(mt_read_ons(f, f, fill_population = FALSE), class = "mt_error_ons")
})

test_that("data frames in the ONS layout are accepted, with regions", {
  d <- data.frame(Year = 2020, Region = rep(c("A", "B"), each = 8), Sex = rep(c("Males", "Females"), each = 4),
                  Age = rep(c("85-89", "90+"), 8), IMD = rep(rep(1:2, each = 2), 4), `All causes` = 1:16,
                  check.names = FALSE)
  p <- cbind(d[, 1:5], Population = 100)
  x <- mt_read_ons(d, p, region = "Region")
  expect_identical(x$by, c("sex", "region", "quintile"))
  expect_equal(nrow(mt_strata(x)), 4)
  expect_equal(sum(mt_obs(x, "deaths")$value), sum(1:16))
  expect_error(mt_read_ons(d[, -5], p), class = "mt_error_ons")
  x <- mt_read_ons(d, p, region = "Region", deprivation = "decile")
  expect_equal(nrow(mt_strata(x)), 8)
  names(d)[5] <- names(p)[5] <- "IMD quintile"
  expect_error(mt_read_ons(d, p, deprivation = "decile"), class = "mt_error_ons")
})
