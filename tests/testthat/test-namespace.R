test_that("exports do not clash with common packages", {
  ours <- getNamespaceExports("mortensor")
  pkgs <- c("base", "stats", "utils", "graphics", "methods", "data.table",
            "dplyr", "forecast", "generics", "demography", "StMoMo", "ggplot2")
  clashes <- character()
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) next
    hit <- intersect(ours, getNamespaceExports(p))
    if (length(hit)) clashes <- c(clashes, paste0(p, "::", hit))
  }
  expect_identical(clashes, character())
})

test_that("conditions carry mortensor classes", {
  err <- tryCatch(mortensor:::mt_abort("boom", class = "x"), error = function(e) e)
  expect_s3_class(err, "mt_error_x")
  expect_s3_class(err, "mt_error")
})

test_that("datasets of suggested packages are found", {
  skip_if_not_installed("demography")
  skip_if_not_installed("StMoMo")
  expect_s3_class(mortensor:::pkg_data("fr.mort", "demography"), "demogdata")
  expect_s3_class(mortensor:::pkg_data("EWMaleData", "StMoMo"), "StMoMoData")
  expect_error(mortensor:::pkg_data("no_such_data", "StMoMo"), class = "mt_error")
})
