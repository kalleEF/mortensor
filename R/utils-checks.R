# Argument checks ---------------------------------------------------------------

check_number <- function(x, name, lower = -Inf, upper = Inf, integerish = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && x >= lower && x <= upper
  if (ok && integerish) ok <- abs(x - round(x)) < 1e-8
  if (!ok) {
    mt_abort(sprintf("`%s` must be a single %s between %s and %s.", name,
                     if (integerish) "whole number" else "number", lower, upper),
             class = "argument", call = NULL)
  }
  invisible(x)
}

check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    mt_abort(sprintf("`%s` must be TRUE or FALSE.", name), class = "argument")
  }
  invisible(x)
}

check_choice <- function(x, choices, name) {
  if (!is.character(x) || length(x) != 1L || !x %in% choices) {
    mt_abort(sprintf("`%s` must be one of %s.", name,
                     paste0("\"", choices, "\"", collapse = ", ")),
             class = "argument", call = NULL)
  }
  invisible(x)
}

check_class <- function(x, cls, name) {
  if (!inherits(x, cls)) {
    mt_abort(sprintf("`%s` must be an object of class <%s>.", name, cls),
             class = "argument", call = NULL)
  }
  invisible(x)
}
