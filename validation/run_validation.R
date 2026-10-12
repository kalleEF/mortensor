#!/usr/bin/env Rscript
# mortensor validation pipeline
#
# Usage:
#   Rscript validation/run_validation.R [--suite quick|full] [--out DIR]
#                                       [--record] [--load-all] [--selftest]
#
# --suite     quick (default; every push) or full (weekly, before releases)
# --out       output directory for results.csv and report.md
# --record    append the results to validation/history.csv (commit the file)
# --load-all  use pkgload::load_all() on the working tree instead of the
#             installed package
# --selftest  add a check that breaches its threshold (pipeline must fail)
#
# Exit status 1 when a metric breaches its threshold or a check errors.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(suite = "quick", out = "validation-output", record = FALSE, load_all = FALSE, selftest = FALSE)
i <- 1
while (i <= length(args)) {
  a <- args[i]
  if (a == "--suite") { opt$suite <- args[i + 1]; i <- i + 1 }
  else if (a == "--out") { opt$out <- args[i + 1]; i <- i + 1 }
  else if (a == "--record") opt$record <- TRUE
  else if (a == "--load-all") opt$load_all <- TRUE
  else if (a == "--selftest") opt$selftest <- TRUE
  else stop("Unknown argument: ", a)
  i <- i + 1
}
stopifnot(opt$suite %in% c("quick", "full"))

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
vdir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else "validation"
root <- dirname(vdir)

if (opt$load_all) {
  pkgload::load_all(root, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(mortensor))
}
source(file.path(vdir, "checks.R"), local = TRUE)
if (opt$selftest) {
  checks$selftest <- list(suite = "quick", description = "Deliberately failing check",
                          run = function() data.frame(metric = "always_breaches", value = 1))
}

suites <- if (opt$suite == "full") c("quick", "full") else "quick"
todo <- Filter(function(ch) ch$suite %in% suites, checks)
thresholds <- utils::read.csv(file.path(vdir, "thresholds.csv"), stringsAsFactors = FALSE,
                              comment.char = "#", strip.white = TRUE)
if (opt$selftest) {
  thresholds <- rbind(thresholds, data.frame(check = "selftest", metric = "always_breaches",
                                             comparator = "<=", threshold = 0, note = "selftest"))
}

rows <- list()
for (nm in names(todo)) {
  cat(sprintf("[%s] %s ... ", format(Sys.time(), "%H:%M:%S"), nm))
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(todo[[nm]]$run(), error = function(e) e)
  secs <- proc.time()[["elapsed"]] - t0
  if (inherits(res, "error")) {
    cat("ERROR:", conditionMessage(res), "\n")
    rows[[nm]] <- data.frame(check = nm, metric = "error", value = NA_real_, seconds = secs,
                             message = conditionMessage(res))
  } else {
    cat(sprintf("%.1fs\n", secs))
    rows[[nm]] <- data.frame(check = nm, metric = res$metric, value = as.numeric(res$value),
                             seconds = secs, message = "")
  }
}
res <- do.call(rbind, rows)
res <- merge(res, thresholds, by = c("check", "metric"), all.x = TRUE, sort = FALSE)
cmp <- function(v, op, t) {
  ifelse(is.na(op), NA, ifelse(op == "<=", v <= t, ifelse(op == ">=", v >= t, NA)))
}
res$pass <- cmp(res$value, res$comparator, res$threshold)
res$status <- ifelse(res$metric == "error", "ERROR", ifelse(is.na(res$pass), "info", ifelse(res$pass, "pass", "FAIL")))
res <- res[order(match(res$check, names(todo))), ]

sha <- tryCatch(system2("git", c("-C", root, "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
if (!length(sha)) sha <- "unknown"
meta <- sprintf("suite %s, commit %s, %s, mortensor %s, %s", opt$suite, sha[1], R.version.string,
                as.character(utils::packageVersion("mortensor")), format(Sys.time(), "%Y-%m-%d %H:%M"))

dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(res, file.path(opt$out, "results.csv"), row.names = FALSE)

fmt <- function(x) ifelse(is.na(x), "", formatC(x, digits = 4, format = "g"))
lines <- c("# mortensor validation report", "", meta, "",
           sprintf("**%d pass, %d fail, %d error, %d info**", sum(res$status == "pass"), sum(res$status == "FAIL"),
                   sum(res$status == "ERROR"), sum(res$status == "info")), "",
           "| Check | Metric | Value | Threshold | Status | Seconds |", "|---|---|---|---|---|---|",
           sprintf("| %s | %s | %s | %s %s | %s | %.1f |", res$check, res$metric, fmt(res$value),
                   ifelse(is.na(res$comparator), "", res$comparator), fmt(res$threshold), res$status, res$seconds))
errs <- res[res$status == "ERROR", ]
if (nrow(errs)) lines <- c(lines, "", "## Errors", "", sprintf("- %s: %s", errs$check, errs$message))
descr <- vapply(todo, function(ch) ch$description, "")
lines <- c(lines, "", "## Checks", "", sprintf("- **%s**: %s", names(descr), descr))
writeLines(lines, file.path(opt$out, "report.md"))
cat(paste(lines, collapse = "\n"), "\n")

if (opt$record) {
  hist_file <- file.path(vdir, "history.csv")
  h <- data.frame(date = format(Sys.Date()), commit = sha[1], suite = opt$suite,
                  r_version = paste(R.version$major, R.version$minor, sep = "."),
                  check = res$check, metric = res$metric, value = res$value, status = res$status)
  utils::write.table(h, hist_file, sep = ",", row.names = FALSE, col.names = !file.exists(hist_file),
                     append = file.exists(hist_file))
}

if (any(res$status %in% c("FAIL", "ERROR"))) {
  cat("\nValidation FAILED\n")
  quit(status = 1)
}
cat("\nValidation passed\n")
