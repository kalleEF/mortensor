#' Objects re-exported from other packages
#'
#' These generics are re-exported so that `forecast()` and `accuracy()` work
#' on mortensor objects whether or not the \pkg{forecast} package is attached.
#'
#' @importFrom generics forecast accuracy
#' @name reexports
#' @keywords internal
#' @export forecast
#' @export accuracy
NULL

#' @importFrom rlang %||%
NULL
