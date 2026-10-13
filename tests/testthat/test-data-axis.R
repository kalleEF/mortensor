test_that("axis constructors produce contiguous intervals", {
  a <- single_ages(30, 100)
  expect_s3_class(a, "mt_axis")
  expect_length(a, 71)
  expect_true(is_open(a))
  expect_equal(a$labels[c(1, 71)], c("30", "100+"))
  g <- age_groups(c(0, 1, seq(5, 90, 5)))
  expect_equal(g$labels[1:3], c("0", "1-4", "5-9"))
  expect_equal(tail(g$labels, 1), "90+")
  y <- single_years(2001, 2005)
  expect_false(is_open(y))
  expect_equal(axis_widths(y), rep(1, 5))
})

test_that("axis validation rejects gaps, overlaps and inner open intervals", {
  expect_error(mt_axis(c(0, 5), c(4, 10)), class = "mt_error_axis")
  expect_error(mt_axis(c(0, 5), c(Inf, 10)), class = "mt_error_axis")
  expect_error(mt_axis(c(5, 0), c(10, 5)), class = "mt_error_axis")
})

test_that("widths and midpoints use omega for the open interval", {
  g <- age_groups(c(80, 85, 90))
  expect_equal(axis_widths(g, omega = 110), c(5, 5, 20))
  expect_equal(axis_midpoints(g, omega = 110), c(82.5, 87.5, 100))
  expect_error(axis_widths(g, omega = 90), class = "mt_error_axis")
})

test_that("parse_age_labels handles common label styles", {
  p <- parse_age_labels(c("<1", "01-14", "15-19", "85 and over", "90+", "0", "1-4", "110+", "under 1"))
  expect_equal(p$age_lo, c(0, 1, 15, 85, 90, 0, 1, 110, 0))
  expect_equal(p$age_hi, c(1, 15, 20, Inf, Inf, 1, 5, Inf, 1))
  expect_error(parse_age_labels("thirty"), class = "mt_error_axis")
})

test_that("axis_locate maps intervals by containment and exact match", {
  g <- age_groups(c(30, 35, 40, 90))
  expect_equal(axis_locate(g, c(30, 33, 36, 95), c(31, 34, 40, Inf)), c(1L, 1L, 2L, 4L))
  expect_true(is.na(axis_locate(g, 33, 37)))
  expect_equal(axis_locate(g, c(30, 30), c(35, 34), mode = "exact"), c(1L, NA))
  expect_equal(axis_locate(g, 90, Inf, mode = "exact"), 4L)
})

test_that("German GENESIS age labels are parsed in any locale", {
  x <- c(paste0("68-J", intToUtf8(0xe4), "hrige"), "unter 1 Jahr", "1 Jahr bis unter 15 Jahre", "85 Jahre und mehr")
  a <- parse_age_labels(x)
  expect_equal(a$age_lo, c(68, 0, 1, 85))
  expect_equal(a$age_hi, c(69, 1, 15, Inf))
})
