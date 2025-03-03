# == Isotropic kernels ==
# These are taken from the following paper and book:
# Genton, Marc. (2001). Classes of Kernels for Machine Learning: A Statistics Perspective.
# Journal of Machine Learning Research. 2. 299-312. 10.1162/15324430260185646.
#
# Hristopulos, D. T. (2020). Random Fields for Spatial Data Modeling:  A Primer for Scientists and Engineers.
# In Advances in Geographic Information Science. Springer Netherlands. https://doi.org/10.1007/978-94-024-1918-4

#' Isotropic 1D kernels.
#'
#' This computes one of the isotropic kernels listed below.
#'
#' \strong{Gaussian Kernel}.
#' This computes the isotropic Gaussian kernel, which is valid for \eqn{{R}^{d}.}
#' \deqn{a(\tau;\theta) = \exp(-\tau^{2} / \theta).}
#' The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.
#'
#' \strong{Exponential Kernel}.
#' This computes the isotropic exponential kernel, which is valid for \eqn{{R}^{d}.}
#' \deqn{a(\tau;\theta) = \exp(-\tau / \theta).}
#' The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.
#'
#' \strong{Isotropic Wave (Cardinal Sine) Kernel}.
#' This computes the isotropic wave (cardinal sine) kernel, which is valid for \eqn{{R}^{3}.}
#' \deqn{a(\tau;\theta) = \left\{ \begin{array}{ll}
#' \frac{\theta}{\tau} \sin\left( \frac{\tau}{\theta} \right), & \tau \neq 0 \\
#' 1, & \tau = 0
#' \end{array} . \right.
#' }
#' The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.
#'
#' \strong{Isotropic Rational Quadratic Kernel}.
#' This computes the isotropic rational quadratic kernel, which is valid for \eqn{{R}^{d}.}
#' \deqn{a(\tau;\theta) = 1 - \frac{\tau^{2}}{\tau^{2} + \theta}.}
#' The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.
#'
#' \strong{Isotropic Spherical Kernel}.
#' This computes the isotropic spherical kernel, which is valid for \eqn{{R}^{3}.}
#' \deqn{a(\tau;\theta) = \left\{
#' \begin{array}{ll}
#' 1 - \frac{3}{2}\frac{\tau}{\theta} + \frac{1}{2}\left( \frac{\tau}{\theta} \right)^{3}, & \tau < \theta \\
#' 0, & \mbox{otherwise}
#' \end{array} . \right.
#' }
#' The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.
#'
#' \strong{Isotropic Circular Kernel}.
#' This computes the isotropic circular kernel, which is valid for \eqn{{R}^{2}.}
#' \deqn{a(\tau;\theta) = \left\{
#' \begin{array}{ll}
#' \frac{2}{\pi}\arccos\left( \frac{\tau}{\theta} \right) - \frac{2}{\pi}\frac{\tau}{\theta} \sqrt{ 1 - \left( \frac{\tau}{\theta} \right)^{2} }, & \tau < \theta \\
#' 0, & \mbox{otherwise}
#' \end{array} . \right.
#' }
#' The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.
#'
#' \strong{Isotropic Matérn Kernel}.
#' This computes the isotropic Matérn kernel, which is valid for \eqn{{R}^{d}} for \eqn{\nu > 0.}
#' \deqn{a(\tau; \theta, \nu) = \left(\sqrt{2\nu} \frac{\tau}{\theta} \right)^{\nu} \left(2^{\nu - 1} \Gamma(\nu) \right)^{-1} K_{\nu}\left( \sqrt{2\nu} \frac{\tau}{\theta} \right) ,}
#' where \eqn{K_{\nu}(\cdot)} is the modified Bessel function of the second kind.
#' The \code{params} argument is of the form \code{c(}\eqn{\theta, \nu}\code{)}.
#'
#' \strong{Isotropic Bessel Kernel}.
#' This computes the isotropic Bessel kernel, which is valid for \eqn{{R}^{d}} for \eqn{\nu \geq \frac{d}{2} - 1.}
#' \deqn{a(\tau; \theta, \nu) = 2^{\nu} \Gamma(\nu + 1) J_{\nu}(\tau / \theta) (\tau / \theta)^{-\nu} ,}
#' where \eqn{J_{\nu}(\cdot)} is the Bessel function of the first kind.
#' The \code{params} argument is of the form \code{c(}\eqn{\theta, \nu, d} \code{)}.
#'
#' \strong{Isotropic Cauchy Kernel}.
#' This computes the isotropic Cauchy kernel, which is valid for \eqn{{R}^{d}} for \eqn{0 < \alpha \leq 2} and \eqn{\beta \geq 0.}
#' \deqn{a(\tau ; \theta, \alpha, \beta) = (1 + (\tau / \theta)^{\alpha})^{-(\beta / \alpha)} .}
#' The \code{params} argument is of the form \code{c(}\eqn{\theta, \alpha, \beta} \code{)}.
#'
#' @references
#' Genton, Marc. (2001). Classes of Kernels for Machine Learning: A Statistics Perspective. Journal of Machine Learning Research. 2. 299-312. 10.1162/15324430260185646.
#' Table 4.2 of Hristopulos, D. T. (2020). Random Fields for Spatial Data Modeling:  A Primer for Scientists and Engineers. In Advances in Geographic Information Science. Springer Netherlands. https://doi.org/10.1007/978-94-024-1918-4
#'
#' @param x A vector or matrix of values of at least length 1
#' @param name The name of the kernel. Options are: "gaussian", "exponential", "wave", "rational_quadratic", "spherical", "circular", "bessel_j", "matern", and "cauchy".
#' @param params A vector of parameters for the kernel. See the documentation below for the position of the parameters. All kernels will have a scale parameter as the first value in the vector.
#'
#' @return A vector of values.
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

      # return(sapply(x, function(t) ifelse(t == 0, 1, ifelse(t == Inf, 0, (params[1] / t) * sin(t / params[1])))))
    }
    else if(name == "rational_quadratic") {
      returnVal <- x
      returnVal_indices <- which(returnVal != Inf)
      returnVal[returnVal == Inf] <- 0
      returnVal[returnVal_indices] <- 1 - (returnVal[returnVal_indices]^2 / (returnVal[returnVal_indices]^2 + params[1]))
      return(returnVal)

      # return(sapply(x, function(t) ifelse(t == Inf, 0, 1 - (t^2 / (t^2 + params[1])))))
    }
    else if(name == "spherical") {
      returnVal <- x
      returnVal_indices <- which(returnVal < params[1])
      returnVal[returnVal_indices] <- 1 - ((3/2) * (returnVal[returnVal_indices] / params[1])) + ((1/2) * (returnVal[returnVal_indices] / params[1])^3)
      returnVal[-returnVal_indices] <- 0
      return(returnVal)

      # return(sapply(x, function(t) ifelse(t < params[1], 1 - ((3/2) * (t / params[1])) + ((1/2) * (t / params[1])^3), 0)))
    }
    else if(name == "circular") {
      returnVal <- x
      returnVal_indices <- which(returnVal < params[1])
      returnVal[returnVal_indices] <- ((2 / pi) * acos(returnVal[returnVal_indices]  / params[1])) - ((2 / pi) * (returnVal[returnVal_indices]  / params[1]) * sqrt(1 - (returnVal[returnVal_indices]  / params[1])^2))
      returnVal[-returnVal_indices] <- 0
      return(returnVal)

      # return(sapply(x, function(t) ifelse(t < params[1], ((2 / pi) * acos(t / params[1])) - ((2 / pi) * (t / params[1]) * sqrt(1 - (t / params[1])^2)), 0)))
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

      # return(sapply(x, function(t) ifelse(t == 0, 1, ifelse(t == Inf, 0, (((sqrt(2 * params[2]) * t / params[1])^params[2]) / (2^(params[2] - 1) * gamma(params[2]))) * besselK(sqrt(2 * params[2]) * t / params[1], params[2])))))
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

      # return(sapply(x, function(t) ifelse(t == 0, 1, ifelse(t == Inf, 0, (2^params[2]) * gamma(params[2] + 1) * (besselJ(t / params[1], params[2]) / ((t / params[1])^params[2]))))))
    }
    else if(name == "cauchy") {
      stopifnot(params[2] > 0, params[2] <= 2, params[3] >= 0)
      return((1 + (x / params[1])^params[2])^(-(params[3] / params[2])))
    }
    stop(paste0("Unknown kernel: ", name))
  }
  stop(paste0("Unknown kernel: ", name))
}
