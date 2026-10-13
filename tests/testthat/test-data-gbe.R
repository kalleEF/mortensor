test_that("GBE exports are parsed with carried row labels and statistical symbols", {
  g <- mt_read_gbe(test_path("fixtures", "gbe-synthetic.csv"))
  expect_setequal(unique(g$region), c("Deutschland", "Bremen"))
  expect_identical(sort(unique(g$year)), c(2019L, 2020L))
  expect_equal(nrow(g), 2 * 2 * 2 * 2 * 3)
  expect_true(all(g[symbol == "-"]$value == 0))
  expect_true(all(is.na(g[symbol %in% c(".", "X")]$value)))
  expect_equal(g[region == "Bremen" & year == 2020 & sex == "Weiblich" & grepl("^Alle", cause) &
                   grepl("^85", age)]$value, 4)
})

test_that("suppressed cells are raked to the published totals", {
  d <- data.frame(region = c("A", "B", "C", "all"), year = 2020, age = "90+", deaths = c(10, NA, NA, 40),
                  exposure = c(100, 300, 100, 500))
  x <- suppressWarnings(mt_data(deaths = d, exposure = d, by = "region"))
  expect_warning(mt_validate(x), NA)
  y <- mt_fill_suppressed(x, mt_hierarchy(region = "all"))
  o <- mt_obs(y)
  expect_equal(o[label == "B"]$value, 22.5)
  expect_equal(o[label == "C"]$value, 7.5)
  expect_identical(o[label == "B"]$provenance, "estimated")
  # two crossed variables: the margins identify the cells here
  tr <- round(outer(c(A = 50, B = 30, C = 20), c(c1 = 0.5, c2 = 0.3, c3 = 0.2)) * 10)
  g <- expand.grid(region = c("A", "B", "C", "all"), cause = c("c1", "c2", "c3", "all"), stringsAsFactors = FALSE)
  full <- rbind(cbind(tr, all = rowSums(tr)), all = c(colSums(tr), sum(tr)))
  g$deaths <- full[cbind(match(g$region, rownames(full)), match(g$cause, colnames(full)))]
  g$exposure <- 1000
  g$year <- 2020
  g$age <- "90+"
  obs <- g
  obs$deaths[obs$region == "B" & obs$cause %in% c("c2", "c3")] <- NA
  obs$deaths[obs$region == "C" & obs$cause == "c3"] <- NA
  x2 <- suppressWarnings(mt_data(deaths = obs, exposure = obs, by = c("region", "cause")))
  y2 <- mt_fill_suppressed(x2, mt_hierarchy(region = "all", cause = "all"))
  o2 <- mt_obs(y2)[provenance == "estimated"]
  expect_equal(o2[label == "B.c2"]$value, 90, tolerance = 1e-6)
  expect_equal(o2[label == "C.c3"]$value, 40, tolerance = 1e-6)
  expect_lt(attr(y2, "fill")$max_gap, 1e-6)
  expect_error(mt_fit(rwd(), x), class = "mt_error_model_data")
})

test_that("31 December populations give exposures and the German data read consistently", {
  p <- data.frame(year = 2018:2020, age = "90+", population = c(100, 110, 130))
  x <- mt_as_exposure(mt_add_population(mt_data(), p), type = "dec31")
  e <- mt_obs(x, "exposure")
  expect_equal(e$period_lo, c(2019, 2020))
  expect_equal(e$value, c(105, 120))
  y <- mt_read_destatis(de_deaths, de_population, years = 2017:2019)
  expect_identical(nrow(mt_strata(y)), 32L)
  expect_identical(mt_grid(y)$age$labels[c(1, 2, 18)], c("0", "1-14", "90+"))
  dd <- as.data.frame(de_deaths)
  pub <- sum(dd$deaths[dd$region == "Deutschland" & dd$cause == "all" & dd$year %in% 2017:2019])
  expect_equal(sum(mt_obs(y)$value), pub, tolerance = 1e-8)
  z <- mt_read_destatis(de_deaths, de_population, causes = de_cause_groups, years = 2018:2019, min_age = 80)
  expect_setequal(unique(mt_strata(z)$cause), c("all", names(de_cause_groups), "rest"))
  expect_true(all(mt_validate(z, hierarchy = mt_hierarchy(cause = "all", additive = "rates"))$severity == "info"))
})
