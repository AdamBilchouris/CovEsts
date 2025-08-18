# == Isotropic kernels ==
# These are taken from the following paper and book:
# Genton, Marc. (2001). Classes of Kernels for Machine Learning: A Statistics Perspective.
# Journal of Machine Learning Research. 2. 299-312. 10.1162/15324430260185646.
#
# Hristopulos, D. T. (2020). Random Fields for Spatial Data Modeling:  A Primer for Scientists and Engineers.
# In Advances in Geographic Information Science. Springer Netherlands. https://doi.org/10.1007/978-94-024-1918-4

#' 1D Isotropic Kernels.
#'
#' This function computes one of the isotropic kernels listed below.
#' Unlike [kernel_symm], these kernels are only defined when \eqn{x \geq 0.} They are used as
#' kernel multipliers in estimators [corrected_est] and [kernel_est].
#'
#' \strong{Gaussian Kernel}.
#' The isotropic Gaussian kernel, which is positive-definite for \eqn{{R}^{d}, d \in N,} is defined as
#' \deqn{a(x;\theta) = \exp(-x^{2} / \theta).}
#' The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.
#'
#' \strong{Exponential Kernel}.
#' The isotropic exponential kernel, which is positive-definite for \eqn{{R}^{d}, d \in N,} is defined as
#' \deqn{a(x;\theta) = \exp(-x / \theta).}
#' The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.
#'
#' \strong{Isotropic Wave (Cardinal Sine) Kernel}.
#' The isotropic wave (cardinal sine) kernel, which is positive-definite for \eqn{{R}^{d}, d \leq 3,} is given by
#' \deqn{a(x;\theta) = \left\{ \begin{array}{ll}
#' \frac{\theta}{x} \sin\left( \frac{x}{\theta} \right), & x \neq 0 \\
#' 1, & x = 0
#' \end{array} . \right.
#' }
#' The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.
#'
#' \strong{Isotropic Rational Quadratic Kernel}.
#' The isotropic rational quadratic kernel, which is positive-definite for \eqn{{R}^{d}, d \in N,} is defined as
#' \deqn{a(x;\theta) = 1 - \frac{x^{2}}{x^{2} + \theta}.}
#' The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.
#'
#' \strong{Isotropic Spherical Kernel}.
#' The isotropic spherical kernel, which is positive-definite for \eqn{{R}^{3}, d \leq 3,} is given by
#' \deqn{a(x;\theta) = \left\{
#' \begin{array}{ll}
#' 1 - \frac{3}{2}\frac{x}{\theta} + \frac{1}{2}\left( \frac{x}{\theta} \right)^{3}, & x < \theta \\
#' 0, & \mbox{otherwise}
#' \end{array} . \right.
#' }
#' The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.
#'
#' \strong{Isotropic Circular Kernel}.
#' The isotropic circular kernel, which is positive-definite for \eqn{{R}^{d}, d \leq 2,} is given by
#' \deqn{a(x;\theta) = \left\{
#' \begin{array}{ll}
#' \frac{2}{\pi}\arccos\left( \frac{x}{\theta} \right) - \frac{2}{\pi}\frac{x}{\theta} \sqrt{ 1 - \left( \frac{x}{\theta} \right)^{2} }, & x < \theta \\
#' 0, & \mbox{otherwise}
#' \end{array} . \right.
#' }
#' The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.
#'
#' \strong{Isotropic Matérn Kernel}.
#' The isotropic Matérn kernel, which is positive-definite for \eqn{{R}^{d}, d \in N,} and when \eqn{\nu > 0,} is defined as
#' \deqn{a(x; \theta, \nu) = \left(\sqrt{2\nu} \frac{x}{\theta} \right)^{\nu} \left(2^{\nu - 1} \Gamma(\nu) \right)^{-1} K_{\nu}\left( \sqrt{2\nu} \frac{x}{\theta} \right) ,}
#' where \eqn{K_{\nu}(\cdot)} is the modified Bessel function of the second kind.
#' The \code{params} argument is of the form \code{c(}\eqn{\theta, \nu}\code{)}.
#'
#' \strong{Isotropic Bessel Kernel}.
#' The isotropic Bessel kernel, which is positive-definite for \eqn{{R}^{d}, d \in N,} and when \eqn{\nu \geq \frac{d}{2} - 1,} is given by
#' \deqn{a(x; \theta, \nu) = 2^{\nu} \Gamma(\nu + 1) J_{\nu}(x / \theta) (x / \theta)^{-\nu} ,}
#' where \eqn{J_{\nu}(\cdot)} is the Bessel function of the first kind.
#' The \code{params} argument is of the form \code{c(}\eqn{\theta, \nu, d} \code{)}.
#'
#' \strong{Isotropic Cauchy Kernel}.
#' The isotropic Cauchy kernel, which is positive-definite for \eqn{{R}^{d}, d \in N,} and when \eqn{0 < \alpha \leq 2} and \eqn{\beta \geq 0,} is defined by
#' \deqn{a(x ; \theta, \alpha, \beta) = (1 + (x / \theta)^{\alpha})^{-(\beta / \alpha)} .}
#' The \code{params} argument is of the form \code{c(}\eqn{\theta, \alpha, \beta} \code{)}.
#'
#' @references
#' Genton, M. (2001). Classes of Kernels for Machine Learning: A Statistics Perspective. Journal of Machine Learning Research. 2, 299-312. 10.1162/15324430260185646.
#'
#' Table 4.2 in Hristopulos, D. T. (2020). Random Fields for Spatial Data Modeling:  A Primer for Scientists and Engineers. Springer. https://doi.org/10.1007/978-94-024-1918-4
#'
#' @param x A vector or matrix of arguments of at least length 1 for which the kernel is computed at.
#' @param name The name of the kernel. Options are: gaussian, exponential, wave, rational_quadratic, spherical, circular, bessel_j, matern, and cauchy.
#' @param params A vector of parameters for the kernel. See the documentation below for the position of the parameters. All kernels have a scale parameter as the first value in the vector.
#'
#' @return A vector or matrix of kernel values.
#' @export
#'
#' @examples
#' x <- c(0.2, 0.4, 0.6)
#' theta <- 0.9
#' kernel(x, "gaussian", c(theta))
#' kernel(x, "exponential", c(theta))
#' kernel(x, "wave", c(theta))
#' kernel(x, "rational_quadratic", c(theta))
#' kernel(x, "spherical", c(theta))
#' kernel(x, "circular", c(theta))
#' nu <- 1
#' kernel(x, "matern", c(theta, nu))
#' dim <- 1
#' kernel(x, "bessel_j", c(theta, nu, dim))
#' alpha <- 1
#' beta <- 2
#' kernel(x, "cauchy", c(theta, alpha, beta))
#' curve(kernel(x, "gaussian", c(theta)), from = 0, to = 5)
#' curve(kernel(x, "exponential", c(theta)), from = 0, to = 5)
#' curve(kernel(x, "wave", c(theta)), from = 0, to = 5)
#' curve(kernel(x, "rational_quadratic", c(theta)), from = 0, to = 5)
#' curve(kernel(x, "spherical", c(theta)), from = 0, to = 5)
#' curve(kernel(x, "circular", c(theta)), from = 0, to = 5)
#' curve(kernel(x, "matern", c(theta, nu)), from = 0, to = 5)
#' curve(kernel(x, "bessel_j", c(theta, nu, dim)), from = 0, to = 5)
#' curve(kernel(x, "cauchy", c(theta, alpha, beta)), from = 0, to = 5)
kernel <- function(x, name, params=c(1)) {
  stopifnot(length(x) >= 1, all(x >= 0), length(params) > 0, is.numeric(params[1]), params[1] > 0)
  stopifnot(name %in% c("gaussian", "exponential", "wave", "rational_quadratic",
                        "spherical", "circular", "bessel_j", "matern", "cauchy"))

  # Gaussian, Exponential, Wave, Rational Quadratic, Spherical, Circular
  if(length(params) == 1) {
    if(name == "gaussian") {
      return(exp(-(x^2) / params[1]))
    }
    else if(name == "exponential") {
      return(exp(-x / params[1]))
    }
    else if(name == "wave") {
      returnVal <- x
      returnVal_indices <- which(returnVal != 0 & returnVal != Inf)
      returnVal[returnVal == 0] <- 1
      returnVal[returnVal == Inf] <- 0
      returnVal[returnVal_indices] <- (params[1] / returnVal[returnVal_indices]) * sin(returnVal[returnVal_indices] / params[1])
      return(returnVal)
    }
    else if(name == "rational_quadratic") {
      returnVal <- x
      returnVal_indices <- which(returnVal != Inf)
      returnVal[returnVal == Inf] <- 0
      returnVal[returnVal_indices] <- 1 - (returnVal[returnVal_indices]^2 / (returnVal[returnVal_indices]^2 + params[1]))
      return(returnVal)
    }
    else if(name == "spherical") {
      returnVal <- x
      returnVal_indices <- which(returnVal < params[1])
      returnVal[returnVal_indices] <- 1 - ((3/2) * (returnVal[returnVal_indices] / params[1])) + ((1/2) * (returnVal[returnVal_indices] / params[1])^3)
      returnVal[-returnVal_indices] <- 0
      return(returnVal)
    }
    else if(name == "circular") {
      returnVal <- x
      returnVal_indices <- which(returnVal < params[1])
      returnVal[returnVal_indices] <- ((2 / pi) * acos(returnVal[returnVal_indices]  / params[1])) - ((2 / pi) * (returnVal[returnVal_indices]  / params[1]) * sqrt(1 - (returnVal[returnVal_indices]  / params[1])^2))
      returnVal[-returnVal_indices] <- 0
      return(returnVal)
    }
    stop(paste0("Unknown kernel: ", name))
  }
  # Matern
  else if(length(params) == 2) {
    if(name == "matern") {
      stopifnot(params[2] > 0)
      returnVal <- x
      returnVal_indices <- which(returnVal != 0 & returnVal != Inf)
      returnVal[returnVal == 0] <- 1
      returnVal[returnVal == Inf] <- 0
      returnVal[returnVal_indices] <- (((sqrt(2 * params[2]) * returnVal[returnVal_indices] / params[1])^params[2]) / (2^(params[2] - 1) * gamma(params[2]))) * besselK(sqrt(2 * params[2]) * returnVal[returnVal_indices] / params[1], params[2])
      return(returnVal)
    }
    stop(paste0("Unknown kernel: ", name))
  }
  else if(length(params) == 3) {
    if(name == "bessel_j") {
      stopifnot(params[3] >= 1, params[2] >= (params[3] / 2) - 1)
      returnVal <- x
      returnVal_indices <- which(returnVal != 0 & returnVal != Inf)
      returnVal[returnVal == 0] <- 1
      returnVal[returnVal == Inf] <- 0
      returnVal[returnVal_indices] <- (2^params[2]) * gamma(params[2] + 1) * (besselJ(returnVal[returnVal_indices] / params[1], params[2]) / ((returnVal[returnVal_indices] / params[1])^params[2]))
      return(returnVal)
    }
    else if(name == "cauchy") {
      stopifnot(params[2] > 0, params[2] <= 2, params[3] >= 0)
      return((1 + (x / params[1])^params[2])^(-(params[3] / params[2])))
    }
    stop(paste0("Unknown kernel: ", name))
  }
  stop(paste0("Unknown kernel: ", name))
}

# == Symmetric kernels ==
# These are to be used with Hall's estimator.
# Unlike the isotropic case, these are symmetric probability distributions.

#' 1D Isotropic Symmetric Kernels.
#'
#' These functions computes values of kernels that have the properties of symmetric probability distributions.
#' For a kernel \eqn{a(x)}, the standardised version is \eqn{a(x) / \int_{-\infty}^{\infty} a(x) dx}, so that the integral is 1.
#' The symmetric kernels are different to [kernel] and are used in the functions [adjusted_est] and [truncated_est].
#'
#' \strong{Symmetric Gaussian Kernel}.
#' The symmetric Gaussian kernel is defined as
#' \deqn{a(x;\theta) = \sqrt{\pi \theta} \exp(-x^{2} / \theta), \theta > 0.}
#' The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.
#'
#' \strong{Symmetric Wave Kernel}.
#' The wave (cardinal sine) kernel is given by
#' \deqn{a(x;\theta) = \left\{ \begin{array}{ll}
#' (\sqrt{\theta^{2}} \pi)^{-1} \frac{\theta}{x} \sin\left( \frac{x}{\theta} \right), & x \neq 0 \\
#' 1, & x = 0
#' \end{array},\right.}
#' where \eqn{\theta > 0.}
#' The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}
#'
#' \strong{Symmetric Rational Quadratic Kernel}.
#' The symmetric rational quadratic kernel is given by
#' \deqn{a(x;\theta) = (\pi \sqrt{\theta})^{-1} \left(1 - \frac{x^{2}}{x^{2} + \theta}\right), \theta > 0.}
#' The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}
#'
#' \strong{Symmetric Besesel Kernel}.
#' The symmetric Bessel kernel, which is valid when \eqn{\nu \geq \frac{d}{2} - 1,} is given by
#' \deqn{a(x; \theta, \nu) = \left(\Gamma\left(\frac{1}{2} + \nu\right)/(2 \sqrt{\pi} \theta \Gamma(1 + \nu))\right) ( 2^{\nu} \Gamma(\nu + 1) J_{\nu}(x / \theta) (x / \theta)^{-\nu}), \,\theta > 0, \nu \geq \frac{d}{2} - 1,}
#' where \eqn{J_{\nu}(\cdot)} is the Bessel function of the first kind and \eqn{d} is the dimension.
#' The \code{params} argument is of the form \code{c(}\eqn{\theta, \nu, d}\code{)}.
#'
#' @param x A vector or matrix of arguments of at least length 1 for which the kernel is computed at. Each value can be negative as well as positive.
#' @param name The name of the kernel. Options are: gaussian, wave, rational_quadratic, bessel_j.
#' @param params A vector of parameters for the kernel. See the documentation below for the position of the parameters. All kernels will have a scale parameter as the first value in the vector.
#'
#' @return A vector or matrix of values.
#' @export
#'
#' @examples
#' x <- c(-2, -1, 0, 1, 2)
#' theta <- 1
#' kernel_symm(x, "gaussian", c(theta))
#' kernel_symm(x, "wave", c(theta))
#' kernel_symm(x, "rational_quadratic", c(theta))
#' dim <- 1
#' nu <- 1
#' kernel_symm(x, "bessel_j", c(theta, nu, dim))
#' curve(kernel_symm(x, "gaussian", c(theta)), from = -5, to = 5)
#' curve(kernel_symm(x, "wave", c(theta)), from = -5, to = 5)
#' curve(kernel_symm(x, "rational_quadratic", c(theta)), from = -5, to = 5)
#' curve(kernel_symm(x, "bessel_j", c(theta, nu, dim)), from = -5, to = 5)
kernel_symm <- function(x, name, params=c(1)) {
  stopifnot(length(x) >= 1, !any(is.na(x)), length(params) > 0, is.numeric(params[1]), params[1] > 0)
  stopifnot(name %in% c("gaussian", "wave", "rational_quadratic", "bessel_j"))

  if(length(params) == 1) {
    if(name == "gaussian") {
      returnVals <- exp(-(x^2) / params[1])
      return(sqrt(pi * params[1])^(-1) * returnVals)
    }
    else if(name == "wave") {
      returnVal <- x
      returnVal_indices <- which(returnVal != 0 & returnVal != Inf)
      returnVal[returnVal == 0] <- 1
      returnVal[returnVal == Inf] <- 0
      returnVal[returnVal_indices] <- (params[1] / returnVal[returnVal_indices]) * sin(returnVal[returnVal_indices] / params[1])
      return((sqrt(params[1]^(2)) * pi)^(-1) * returnVal)
    }
    else if(name == "rational_quadratic") {
      returnVal <- x
      returnVal_indices <- which(returnVal != Inf)
      returnVal[returnVal == Inf] <- 0
      returnVal[returnVal_indices] <- 1 - (returnVal[returnVal_indices]^2 / (returnVal[returnVal_indices]^2 + params[1]))
      return(((pi * sqrt(params[1]))^(-1)) * returnVal)
    }
    stop(paste0("Unknown kernel: ", name))
  }

  else if(length(params) == 3) {
    if(name == "bessel_j") {
      returnVal <- x
      returnVal_indices <- which(returnVal != 0 & returnVal != Inf)
      returnVal[returnVal == 0] <- 1 /  ( (2^params[2]) * gamma(params[2] + 1) )
      returnVal[returnVal == Inf] <- 0
      returnVal[returnVal_indices] <- (besselJ(abs(returnVal)[returnVal_indices] / params[1], params[2]) / ((abs(returnVal[returnVal_indices]) / params[1])^params[2]))
      return((gamma((1/2 + params[2])) / (2*sqrt(pi) * params[1] * gamma(1 + params[2]))) * (2^params[2]) * gamma(params[2] + 1) * returnVal)
    }
    stop(paste0("Unknown kernel: ", name))
  }

  stop(paste0("length(params) is 1 or 3."))
}
