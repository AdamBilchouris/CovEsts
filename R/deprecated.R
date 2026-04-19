#' Deprecated functions in package \pkg{CovEsts}.
#'
#' These functions are provided for compatibility with older verisons of the \pkg{CovEsts} package.
#' These functions may be removed entirely.
#' @rdname CovEsts-deprecated
#' @name CovEsts-deprecated
#' @keywords internal
#'
#' @param x,name,params Parameters passed to the updated functions.
#'
#' @section \code{kernel}:
#' For \code{kernel()}, use \code{\link{kernel_ec}} instead.
#' @section \code{kernel_symm}:
#' For \code{kernel_symm()}, use \code{\link{kernel_symm_ec}} instead.
#' @section \code{window}:
#' For \code{window()}, use \code{\link{window_ec}} instead.
#' @section \code{window_symm}:
#' For \code{window_symm()}, use \code{\link{window_symm_ec}} instead.
#'
NULL

#' @export
#' @rdname CovEsts-deprecated
kernel <- function(x, name, params=c(1)) {
  .Deprecated("kernel_ec", "CovEsts")
  kernel_ec(x = x, name = name, params = params)
}

#' @export
#' @rdname CovEsts-deprecated
kernel_symm <- function(x, name, params=c(1)) {
  .Deprecated("kernel_symm_ec", "CovEsts")
  kernel_symm_ec(x = x, name = name, params = params)
}

#' @export
#' @rdname CovEsts-deprecated
window <- function(x, name, params=c(1)) {
  .Deprecated("window_ec", "CovEsts")
  window_ec(x = x, name = name, params = params)
}

#' @export
#' @rdname CovEsts-deprecated
window_symm <- function(x, name, params=c(1)) {
  .Deprecated("window_symm_ec", "CovEsts")
  window_symm_ec(x = x, name = name, params = params)
}
