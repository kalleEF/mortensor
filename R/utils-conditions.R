# Classed conditions without extra dependencies -------------------------------

mt_abort <- function(message, class = NULL, ..., call = sys.call(-1)) {
  cond <- errorCondition(paste(message, collapse = "\n"),
                         class = c(if (!is.null(class)) paste0("mt_error_", class),
                                   "mt_error"),
                         call = tidy_call(call), ...)
  stop(cond)
}

mt_warn <- function(message, class = NULL, ..., call = sys.call(-1)) {
  cond <- warningCondition(paste(message, collapse = "\n"),
                           class = c(if (!is.null(class)) paste0("mt_warning_", class),
                                     "mt_warning"),
                           call = tidy_call(call), ...)
  warning(cond)
}

# The call shown with a condition: calls made through do.call() carry the
# function and the data themselves, which would be printed in full; such
# calls are shortened to the function name, or dropped
tidy_call <- function(call) {
  if (!is.call(call)) return(call)
  if (!is.name(call[[1]])) return(NULL)
  if (sum(nchar(deparse(call, nlines = 20L))) > 300) return(as.call(list(call[[1]])))
  call
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
  out <- tryCatch(getExportedValue(pkg, name), error = function(e) NULL)
  if (is.null(out)) {
    env <- new.env(parent = emptyenv())
    suppressWarnings(utils::data(list = name, package = pkg, envir = env))
    out <- env[[name]]
  }
  if (is.null(out)) mt_abort(sprintf("Dataset %s not found in package %s.", name, pkg), class = "missing_package")
  out
}
