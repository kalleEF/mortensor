genesis_frame <- function(time, vars, value) {
  out <- data.frame(time = time)
  for (i in seq_along(vars)) {
    out[[sprintf("%d_variable_code", i)]] <- vars[[i]]$code
    out[[sprintf("%d_variable_attribute_code", i)]] <- vars[[i]]$attr_code %||% ""
    out[[sprintf("%d_variable_attribute_label", i)]] <- vars[[i]]$label
  }
  out$value <- value
  out
}

test_that("GENESIS flat tables give deaths and exposures from 31 December populations", {
  umlaut <- paste0("m", intToUtf8(0xe4), "nnlich")
  ages_d <- c("unter 1 Jahr", "1 bis unter 85 Jahre", "85 Jahre und mehr")
  d <- genesis_frame(rep(2020:2021, each = 6),
                     list(list(code = "DINSG", label = "Deutschland"),
                          list(code = "GES", label = rep(rep(c(umlaut, "weiblich"), each = 3), 2)),
                          list(code = "ALT027", label = rep(ages_d, 4))),
                     value = c("10", "-", "300", "8", "5", "400", "11", "3", ".", "9", "4", "410"))
  pop_ages <- c("unter 1 Jahr", paste0(1:84, "-J", intToUtf8(0xe4), "hrige"), "85 Jahre und mehr")
  grid <- expand.grid(age = pop_ages, sex = c(umlaut, "weiblich"), land = c("A", "B"), time = 2019:2021,
                      stringsAsFactors = FALSE)
  p <- genesis_frame(paste0(grid$time, "-12-31"),
                     list(list(code = "DLAND", label = grid$land), list(code = "GES", label = grid$sex),
                          list(code = "ALT103", label = grid$age)),
                     value = 100 + (grid$time - 2019) * 10)
  expect_warning(x <- mt_read_genesis(d, p), class = "mt_warning_genesis")
  expect_identical(x$by, "sex")
  expect_identical(mt_grid(x)$age$labels, c("0", "1-84", "85+"))
  e <- mt_obs(x, "exposure")
  # two Laender, population 100 at the end of 2019 and 110 at the end of 2020: exposure 2020 = 2 x 105
  expect_equal(e[period_lo == 2020 & age_lo == 0 & sex == "male"]$value, 210)
  expect_equal(e[period_lo == 2020 & age_lo == 1 & sex == "male"]$value, 84 * 210)
  dd <- mt_obs(x, "deaths")
  expect_equal(dd[period_lo == 2020 & age_lo == 1 & sex == "male"]$value, 0)
  expect_equal(nrow(dd), 11)
  expect_error(mt_read_genesis(d[, 1:3], p), class = "mt_error_genesis")
})

test_that("GENESIS cause variables become a cause stratum with a closing rest", {
  ages_d <- c("unter 1 Jahr", "1 bis unter 85 Jahre", "85 Jahre und mehr")
  d <- genesis_frame(rep(2020, 6),
                     list(list(code = "GES", label = "weiblich"), list(code = "ALT027", label = rep(ages_d, 2)),
                          list(code = "TODUR1", attr_code = rep(c("", "TDU-07"), each = 3),
                               label = rep(c("Insgesamt", "Krankheiten des Kreislaufsystems"), each = 3))),
                     value = c(10, 50, 300, 0, 10, 120))
  p <- genesis_frame(rep(c("2019-12-31", "2020-12-31"), each = 3),
                     list(list(code = "GES", label = "weiblich"), list(code = "ALT103", label = rep(ages_d, 2))),
                     value = 1000)
  x <- mt_read_genesis(d, p, causes = TRUE)
  expect_setequal(mt_strata(x)$cause, c("all", "TDU-07", "rest"))
  expect_identical(attr(x, "cause_labels")$cause, "TDU-07")
  o <- mt_obs(x)
  expect_equal(o[cause == "rest"]$value, c(10, 40, 180))
  y <- mt_read_genesis(d, p)
  expect_false("cause" %in% names(mt_strata(y)))
})
