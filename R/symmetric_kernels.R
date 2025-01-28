# == Symmetric kernels ==
# These are to be used with Hall's estimator.
# Unlike the isotropic case, these are symmetric probability distributions.

#' 1D isotropic symmetric kernels.
#'
#' These have the properties of symmetric probability   distributions.
#'
#' \strong{Symmetric Gaussian Kernel}.
#' \deqn{a(\tau;\theta) = \sqrt{\pi \theta} \exp(-\tau^{2} / \theta), \theta > 0}.
#' The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.
#'
#' \strong{Symmetric Wave Kernel}.
#' \deqn{a(\tau;\theta) = \begin{cases}
#' (\sqrt{a^{-2}} / \pi) \frac{\theta}{\tau} \sin\left( \frac{\tau}{\theta} \right), & \tau \neq = 0 \\
#' 1, & \tau = 0
#' \end{cases} ,
#' \theta > 0.}
#' The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}
#'
#' \strong{Symmetric Rational Quadratic Kernel}.
#' \deqn{a(\tau;\theta) = (\pi \sqrt{\theta})^{-1} (1 - \frac{\tau^{2}}{\tau^{2} + \theta}), \theta > 0.}
#' The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}
#'
#' \strong{Symmetric Besesel Kernel}.
#' \deqn{a(\tau; \nu) = (\Gamma(\frac{1}{2} + \nu)/(2 \sqrt{\pi} \theta \Gamma(1 + \nu))) ( 2^{\nu} \Gamma(\nu + 1) J_{\nu}(\tau) \tau^{-\nu}), \theta > 0, \nu \geq \frac{d}{2} - 1.}
#' where \eqn{J_{\nu}(\cdot)} is the Bessel function of the first kind and \eqn{d} is the dimension.
#' The \code{params} argument is of the form \code{c(}\eqn{\theta, \nu, d}\code{)}.
#'
#' @param x A vector of values of at least length 1.
#' @param name The name of the kernel. Options are: "gaussian", "wave", "rational_quadratic", "bessel_j".
#' @param params A vector of parameters for the kernel. See the documentation below for the position of the parameters. All kernels will have a scale parameter as the first value in the vector.
#'
#' @return A vector of values.
#' @export
#'
#' @examples
#' x <- c(-3, -2, -1, 0, 1, 2, 3)
#' kernel_symm(x, "gaussian", c(0.2))
#' kernel_symm(x, "wave", c(0.2))
#' kernel_symm(x, "rational_quadratic", c(0.2))
#' kernel_symm(x, "bessel_j", c(0.2, 2, 2))
kernel_symm <- function(x, name, params=c(1)) {
  stopifnot(length(x) >= 1, !any(is.na(x)), is.vector(params), length(params) > 0, is.numeric(params[1]), params[1] > 0)
  stopifnot(name %in% c("gaussian", "wave", "rational_quadratic", "bessel_j"))

  if(length(params) == 1) {
    if(name == "gaussian") {
      returnVals <- exp(-(x^2) / params[1])
      return(sqrt(pi * params[1])^(-1) * returnVals)
    }
    else if(name == "wave") {
      returnVals <- sapply(x, function(t) ifelse(t == 0, 1, ifelse(t == Inf, 0, (params[1] / t) * sin(t / params[1]))))
      return((sqrt(params[1]^(-2)) / pi) * returnVals)
    }
    else if(name == "rational_quadratic") {
      returnVals <- sapply(x, function(t) ifelse(t == Inf, 0, 1 - (t^2 / (t^2 + params[1]))))
      return(((pi * sqrt(params[1]))^(-1)) * returnVals)
    }
    stop(paste0("Unknown kernel: ", name))
  }

  else if(length(params) == 3) {
    if(name == "bessel_j") {
      returnVal <- sapply(x, function(t) ifelse(t == 0, 1, ifelse(t == Inf, 0, (2^params[2]) * gamma(params[2] + 1) * (besselJ(t / params[1], params[2]) / ((t / params[1])^params[2])))))
      return((gamma((1/2) + params[2]) / (2*sqrt(pi) * params[1] * gamma(1 + params[2]))) * returnVal)
    }
    stop(paste0("Unknown kernel: ", name))
  }

  stop(paste0("length(params) is 1 or 3."))
}
