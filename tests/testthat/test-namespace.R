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
