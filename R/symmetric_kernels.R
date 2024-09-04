# == Symmetric kernels ==
# These are to be used with Hall's estimator.
# Unlike the isotropic case, these are symmetric probability distributions.

#' Symmetric Gaussian Kernel.
#'
#' This computes the Gaussian kernel. It has the properties of a probability distribution, unlike `kernel_gaussian(x, theta)`.
#' \deqn{a(\tau;\theta) = \sqrt{\pi \theta} \exp(-\tau^{2} / \theta).}
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
#' kernel_symm_gaussian(x, 1)
#' x <- c(-3, -2, -1, 0, 1, 2, 3)
#' kernel_symm_gaussian(x, 0.2)
kernel_symm_gaussian <- function(x, theta, ...) {
  stopifnot(theta > 0, length(x) >= 1, !any(is.na(x)))
  returnVals <- exp(-(x^2) / theta)
  return(sqrt(pi * theta)^(-1) * returnVals)
}

#' Symmetric Wave (Cardinal Sine) Kernel.
#'
#' This computes the wave (cardinal sine) kernel. It has the properties of a probability distribution, unlike `kernel_wave(x, theta)`.
#' \deqn{a(\tau;\theta) = \begin{cases}
#' (\sqrt{a^{-2}} / \pi) \frac{\theta}{\tau} \sin\left( \frac{\tau}{\theta} \right), & \tau \neq = 0 \\
#' 1, & \tau = 0
#' \end{cases} .
#' }
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
#' kernel_symm_wave(x, 1)
#' x <- c(-3, -2, -1, 0, 1, 2, 3)
#' kernel_symm_wave(x, 0.2)
kernel_symm_wave <- function(x, theta, ...) {
  stopifnot(theta > 0, length(x) >= 1, !any(is.na(x)))
  returnVals <- sapply(x, function(t) ifelse(t == 0, 1, ifelse(t == Inf, 0, (theta / t) * sin(t / theta))))
  return((sqrt(theta^(-2)) / pi) * returnVals)
}

#' Symmetric Rational Quadratic Kernel.
#'
#' This computes the rational quadratic kernel. It has the properties of a probability distribution, unlike `kernel_rational_quadratic(x, theta)`.
#' \deqn{a(\tau;\theta) = (\pi \sqrt{\theta})^{-1} (1 - \frac{\tau^{2}}{\tau^{2} + \theta}).}
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
#' kernel_symm_rational_quadratic(x, 1)
#' x <- c(-3, -2, -1, 0, 1, 2, 3)
#' kernel_symm_rational_quadratic(x, 0.2)
kernel_symm_rational_quadratic <- function(x, theta, ...) {
  stopifnot(theta > 0, length(x) >= 1, !any(is.na(x)))
  returnVals <- sapply(x, function(t) ifelse(t == Inf, 0, 1 - (t^2 / (t^2 + theta))))
  return(((pi * sqrt(theta))^(-1)) * returnVals)
}

#' Symmetric Bessel Kernel.
#'
#' This computes the Bessel kernel. It has the properties of a probability distribution, unlike `kernel_bessel_j(x, theta)`.
#' \deqn{a(\tau; \nu) = (\Gamma(\frac{1}{2} + \nu)/(2 \sqrt{\pi} \theta \Gamma(1 + \nu))) ( 2^{\nu} \Gamma(\nu + 1) J_{\nu}(\tau) \tau^{-\nu}) ,}
#' where \eqn{J_{\nu}(\cdot)} is the Bessel function of the first kind.
#'
#' @references
#' Table 4.2 of Hristopulos, D. T. (2020). Random Fields for Spatial Data Modeling:  A Primer for Scientists and Engineers. In Advances in Geographic Information Science. Springer Netherlands. https://doi.org/10.1007/978-94-024-1918-4
#'
#' @param x A vector of values of at least length 1.
#' @param theta Scale parameter, \eqn{\theta > 0.}
#' @param nu Parameter defining strength of periodicity \eqn{\nu \geq \frac{\text{dim}}{2} - 1.}
#' @param dim The dimension of the isotropic kernel.
#'
#' @return A vector of values.
#' @export
#'
#' @examples
#' x <- 1
#' kernel_symm_bessel_j(x, 1, -0.5, 1)
#' x <- c(-3, -2, -1, 0, 1, 2, 3)
#' kernel_symm_bessel_j(x, 1, 2, 2)
kernel_symm_bessel_j <- function(x, theta, nu, dim) {
  stopifnot(theta > 0, nu >= (dim/2) - 1, length(x) >= 1, !any(is.na(x)))
  returnVal <- sapply(x, function(t) ifelse(t == 0, 1, ifelse(t == Inf, 0, (2^nu) * gamma(nu + 1) * (besselJ(t / theta, nu) / ((t / theta)^nu)))))
  return((gamma((1/2) + nu) / (2*sqrt(pi) * theta * gamma(1 + nu))) * returnVal)
}
