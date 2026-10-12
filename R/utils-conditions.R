# Classed conditions without extra dependencies -------------------------------

mt_abort <- function(message, class = NULL, ..., call = sys.call(-1)) {
  cond <- errorCondition(paste(message, collapse = "\n"),
                         class = c(if (!is.null(class)) paste0("mt_error_", class),
                                   "mt_error"),
                         call = call, ...)
  stop(cond)
}

mt_warn <- function(message, class = NULL, ..., call = sys.call(-1)) {
  cond <- warningCondition(paste(message, collapse = "\n"),
                           class = c(if (!is.null(class)) paste0("mt_warning_", class),
                                     "mt_warning"),
                           call = call, ...)
  warning(cond)
}

mt_inform <- function(...) {
  if (isTRUE(getOption("mortensor.verbose", FALSE))) message(...)
  invisible(NULL)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

need_pkg <- function(pkg, why) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    mt_abort(sprintf("Package '%s' is needed %s. Install it with install.packages(\"%s\").",
                     pkg, why, pkg), class = "missing_package", call = NULL)
  }
  invisible(TRUE)
}

# Load a dataset shipped with a suggested package
pkg_data <- function(name, pkg) {
  need_pkg(pkg, sprintf("for the dataset %s", name))
  env <- new.env(parent = emptyenv())
  utils::data(list = name, package = pkg, envir = env)
  env[[name]]
}
