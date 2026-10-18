test_that("relabelling a level to an existing level combines both", {
  d <- data.frame(sex = c("F", "female", "male"), year = 2001, age = "90+", deaths = c(5, 6, 7),
                  exposure = c(50, 60, 70))
  x <- mt_data(deaths = d, exposure = d, by = "sex")
  y <- mt_relabel(x, "sex", c(F = "female"))
  expect_identical(mt_strata(y)$sex, c("female", "male"))
  expect_equal(mt_obs(y)[sex == "female"]$value, 11)
  expect_equal(mt_obs(y, "exposure")[sex == "female"]$value, 110)
  expect_equal(sum(mt_obs(y)$value), 18)
  # the same with mt_group_levels() and an unlisted level named like a group
  z <- mt_group_levels(x, "sex", list(female = "F"))
  expect_equal(mt_obs(z)[sex == "female"]$value, 11)
  # a group named like the total level would be overwritten
  expect_error(mt_group_levels(x, "sex", list(male = "F"), total = "male"), class = "mt_error_argument")
})

test_that("numeric total, other and keep levels are used as labels", {
  d <- data.frame(code = 1:4, year = 2020, age = "90+", deaths = c(60, 10, 20, 30), exposure = 1000)
  x <- mt_data(deaths = d, exposure = d, by = "code")
  y <- mt_group_levels(x, "code", list(g23 = c(2, 3)), other = "rest", total = 1, additive = "rates")
  y_chr <- mt_group_levels(x, "code", list(g23 = c("2", "3")), other = "rest", total = "1", additive = "rates")
  expect_identical(mt_obs(y), mt_obs(y_chr))
  expect_equal(mt_obs(y)[match(c("g23", "1", "rest"), label)]$value, c(30, 60, 30))
  expect_identical(y$level_totals$code$total, "1")
  k <- mt_group_levels(x, "code", list(g23 = c(2, 3)), other = "rest", keep = 4)
  expect_setequal(mt_strata(k)$code, c("g23", "4", "rest"))
  d0 <- data.frame(code = 0:3, year = 2020, age = "90+", deaths = c(60, 10, 20, 30), exposure = 1000)
  x0 <- mt_data(deaths = d0, exposure = d0, by = "code")
  y0 <- mt_group_levels(x0, "code", list(g12 = c(1, 2)), other = "rest", total = 0, additive = "rates")
  expect_equal(mt_obs(y0)[match(c("g12", "0", "rest"), label)]$value, c(30, 60, 30))
})

test_that("GBE values with several thousands separators are read in full", {
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f))
  writeLines(c("\"Sterbefaelle, Gliederungsmerkmale: Jahre, Region, Alter, Geschlecht\"",
               ";;;\"2021\";\"2021\";\"2022\";\"2022\"",
               ";;;\"Alle Todesursachen\";\"A00-B99\";\"Alle Todesursachen\";\"A00-B99\"",
               paste0("\"Deutschland\";\"Alle Altersgruppen\";\"Beide Geschlechter\";\"1.023.687\";\"12.345\";",
                      "\"1.066.341\";\".\";"),
               "\"Bremen\";\"Alle Altersgruppen\";\"Beide Geschlechter\";\"7.432\";\"-\";\"7.500,5\";\"X\";",
               "\"Quelle: Statistisches Bundesamt\""), f, useBytes = TRUE)
  g <- mt_read_gbe(f)
  expect_setequal(unique(g$region), c("Deutschland", "Bremen"))
  v <- function(r, y, c) g[region == r & year == y & cause == c]$value
  expect_equal(v("Deutschland", 2021, "Alle Todesursachen"), 1023687)
  expect_equal(v("Deutschland", 2022, "Alle Todesursachen"), 1066341)
  expect_equal(v("Deutschland", 2021, "A00-B99"), 12345)
  expect_equal(v("Bremen", 2021, "Alle Todesursachen"), 7432)
  expect_equal(v("Bremen", 2022, "Alle Todesursachen"), 7500.5)
  expect_equal(v("Bremen", 2021, "A00-B99"), 0)
  expect_true(is.na(v("Deutschland", 2022, "A00-B99")))
})

test_that("cells reported at another resolution are not taken as suppressed", {
  d <- rbind(data.frame(region = "A", year = 2020, age = c("85-89", "90+"), deaths = c(10, 20), exposure = c(100, 50)),
             data.frame(region = "B", year = 2020, age = "85+", deaths = 25, exposure = 120))
  x <- mt_data(deaths = d, exposure = d, by = "region")
  y <- mt_fill_suppressed(x, lower = 1, upper = 2)
  expect_identical(mt_obs(y), mt_obs(x))
  expect_identical(attr(y, "fill")$filled, 0L)
  # a suppressed cell next to a stratum at another resolution is still filled
  d2 <- rbind(d, data.frame(region = "C", year = 2020, age = c("85-89", "90+"), deaths = c(30, NA),
                            exposure = c(300, 40)))
  x2 <- suppressWarnings(mt_data(deaths = d2, exposure = d2, by = "region"))
  y2 <- mt_fill_suppressed(x2, lower = 1, upper = 2)
  est <- mt_obs(y2)[provenance == "estimated"]
  expect_identical(est$label, "C")
  expect_identical(est$age_lo, 90)
  expect_true(est$value >= 1 && est$value <= 2)
  expect_identical(nrow(mt_obs(y2)[label == "B"]), 1L)
  # with published totals: a member at another resolution leaves that cell unconstrained
  d3 <- rbind(d2, data.frame(region = "all", year = 2020, age = c("85-89", "90+"), deaths = c(50, 60),
                             exposure = c(500, 200)))
  x3 <- suppressWarnings(mt_data(deaths = d3, exposure = d3, by = "region"))
  y3 <- mt_fill_suppressed(x3, mt_hierarchy(region = "all"), lower = 1, upper = 2)
  o3 <- mt_obs(y3)
  expect_identical(nrow(o3[label == "B"]), 1L)
  expect_true(o3[label == "C" & age_lo == 90]$value <= 2)
})

test_that("missing numeric ages and periods are unknown", {
  d <- data.frame(year = 2001, age = c(80, 81, NA), deaths = c(10, 10, 4))
  expect_warning(x <- mt_data(deaths = d), class = "mt_warning_unknown_age")
  expect_equal(mt_obs(x)$value, c(12, 12))
  d2 <- data.frame(year = 2001, age = c("80", "81", NA), deaths = c(10, 10, 4))
  expect_warning(x2 <- mt_data(deaths = d2), class = "mt_warning_unknown_age")
  expect_equal(mt_obs(x2)$value, c(12, 12))
  d3 <- data.frame(year = c(2001, NA), age = 80, deaths = c(10, 3))
  expect_warning(x3 <- mt_data(deaths = d3), class = "mt_warning_missing")
  expect_equal(mt_obs(x3)$value, 10)
  d4 <- data.frame(year = c(2001, 2001, 2001), lo = c(80, 85, NA), hi = c(85, Inf, NA), deaths = c(10, 10, 4))
  expect_warning(x4 <- mt_data(deaths = d4, age = c("lo", "hi")), class = "mt_warning_unknown_age")
  expect_equal(mt_obs(x4)$value, c(12, 12))
})

test_that("two-column bounds given as factors are read as their labels", {
  d <- data.frame(year = 2020, age_lo = factor(c(80, 85, 90)), age_hi = factor(c(85, 90, Inf)), deaths = c(10, 20, 30))
  x <- mt_data(deaths = d, age = c("age_lo", "age_hi"))
  expect_equal(mt_obs(x)$age_lo, c(80, 85, 90))
  expect_equal(mt_obs(x)$age_hi, c(85, 90, Inf))
  d$age_hi <- c("85", "90", "open")
  expect_error(mt_data(deaths = d, age = c("age_lo", "age_hi")), class = "mt_error_obs")
})

test_that("mt_from_hmd() reads a '+' age as open without column OpenInterval", {
  d <- data.frame(Year = 2000, Age = c("108", "109", "110+"), Female = c(10, 5, 4), Male = c(12, 2, 1),
                  Total = c(22, 7, 5))
  e <- transform(d, Female = 100, Male = 100, Total = 200)
  x <- mt_from_hmd(d, e)
  o <- mt_obs(x)[label == "female"]
  expect_equal(o$age_hi, c(109, 110, Inf))
})

test_that("mt_example_regions() works for a single year", {
  x <- mt_example_regions(years = 2019, regions = 2)
  expect_s3_class(x, "mt_data")
  expect_identical(unique(mt_obs(x)$period_lo), 2019)
  expect_identical(dim(attr(x, "truth")), c(11L, 1L, 20L))
})


test_that("suppressed cells without a published total get the truncated Poisson mean", {
  d <- data.frame(region = c("A", "B", "C"), year = 2020, age = "90+", deaths = c(30, NA, NA),
                  exposure = c(1000, 50, 100))
  x <- suppressWarnings(mt_data(deaths = d, exposure = d, by = "region"))
  y <- mt_fill_suppressed(x, lower = 1, upper = 2)
  est <- mt_obs(y)[provenance == "estimated"][order(label)]
  # expected deaths 1.5 and 3 at the pooled rate 0.03: E[Y | 1 <= Y <= 2] = (1 + lambda) / (1 + lambda / 2)
  lam <- c(1.5, 3)
  expect_equal(est$value, (1 + lam) / (1 + lam / 2), tolerance = 1e-10)
  # wider bounds: the mean over 0..4
  z <- mt_fill_suppressed(x, lower = 0, upper = 4)
  k <- 0:4
  tm <- vapply(lam, function(l) sum(k * stats::dpois(k, l)) / sum(stats::dpois(k, l)), 0)
  expect_equal(mt_obs(z)[provenance == "estimated"][order(label)]$value, tm, tolerance = 1e-10)
  # limits: no deaths expected, far above the bounds, bounds without an integer, wide bounds
  expect_equal(mortensor:::truncated_poisson_mean(c(0, 1e6), 1, 2), c(1, 2), tolerance = 1e-5)
  expect_equal(mortensor:::truncated_poisson_mean(0.5, 0.2, 0.8), 0.5)
  expect_equal(mortensor:::truncated_poisson_mean(5, 0, 1e6), 5, tolerance = 1e-8)
})
