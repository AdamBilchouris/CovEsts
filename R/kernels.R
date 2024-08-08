# == Isotropic kernels ==
# These are taken from the following paper:
# Genton, Marc. (2001). Classes of Kernels for Machine Learning: A Statistics Perspective.
# Journal of Machine Learning Research. 2. 299-312. 10.1162/15324430260185646.

#' Isotropic Gaussian Kernel.
#'
#' This computes the isotropic Gaussian kernel, which is valid for \eqn{\mathbb{R}^{d}.}
#' \deqn{a(\tau) = \exp(-\tau^{2} / \theta).}
#'
#' @references Genton, Marc. (2001). Classes of Kernels for Machine Learning: A Statistics Perspective. Journal of Machine Learning Research. 2. 299-312. 10.1162/15324430260185646.
#'
#' @param x A vector of values of at least length 1.
#' @param theta Scale parameter, \eqn{\theta > 0.}
#'
#' @return A vector of values.
#' @export
#'
#' @examples
#' x <- 1
#' kernel_gaussian(x, 1)
#' x <- c(1, 2, 3)
#' kernel_gaussian(x, 0.2)
kernel_gaussian <- function(x, theta) {
  stopifnot(theta > 0, length(x) >= 1, all(x >= 0))
  return(exp(-(x^2) / theta))
}

#' Isotropic Exponential Kernel.
#'
#' This computes the isotropic Gaussian kernel, which is valid for \eqn{\mathbb{R}^{d}.}
#' \deqn{a(\tau) = \exp(-\tau / \theta).}
#'
#' @references
#' Genton, Marc. (2001). Classes of Kernels for Machine Learning: A Statistics Perspective. Journal of Machine Learning Research. 2. 299-312. 10.1162/15324430260185646.
#'
#' @param x A vector of values of at least length 1.
#' @param theta Scale parameter, \eqn{\theta > 0.}
#'
#' @return A vector of values.
#' @export
#'
#' @examples
#' x <- 1
#' kernel_exponential(x, 1)
#' x <- c(1, 2, 3)
#' kernel_exponential(x, 0.2)
kernel_exponential <- function(x, theta) {
  stopifnot(theta > 0, length(x) >= 1, all(x >= 0))
  return(exp(-(x) / theta))
}

