# Contract test reused for every model: shapes, reproducibility, positivity
expect_mt_model <- function(spec, data = NULL, h = 3, n_draws = 20) {
  data <- data %||% mt_example(deaths = "single", years = 1995:2012)
  fit <- mt_fit(spec, data)
  testthat::expect_s3_class(fit, "mt_fit")
  testthat::expect_true(all(fit$diagnostics$status == "ok"))
  fc1 <- mt_forecast(fit, h = h, n_draws = n_draws, seed = 42)
  fc2 <- mt_forecast(fit, h = h, n_draws = n_draws, seed = 42)
  testthat::expect_s3_class(fc1, "mt_draws")
  g <- mt_grid(data)
  testthat::expect_equal(dim(fc1$rates), c(length(g$age), h, nrow(mt_strata(data)), n_draws))
  testthat::expect_identical(fc1$rates, fc2$rates)
  testthat::expect_true(all(is.finite(fc1$rates) & fc1$rates > 0))
  testthat::expect_equal(fc1$periods, max(g$period$lower) + seq_len(h))
  s <- mt_summarise(fc1)
  testthat::expect_true(all(c("q0.025", "q0.5", "q0.975") %in% names(s)))
  testthat::expect_true(all(s$q0.025 <= s$q0.5 & s$q0.5 <= s$q0.975))
  invisible(list(fit = fit, fc = fc1))
}

`%||%` <- function(x, y) if (is.null(x)) y else x
