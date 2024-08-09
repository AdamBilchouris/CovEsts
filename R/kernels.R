# == Isotropic kernels ==
# These are taken from the following paper and book:
# Genton, Marc. (2001). Classes of Kernels for Machine Learning: A Statistics Perspective.
# Journal of Machine Learning Research. 2. 299-312. 10.1162/15324430260185646.
#
# Hristopulos, D. T. (2020). Random Fields for Spatial Data Modeling:  A Primer for Scientists and Engineers.
# In Advances in Geographic Information Science. Springer Netherlands. https://doi.org/10.1007/978-94-024-1918-4

#' Isotropic Gaussian Kernel.
#'
#' This computes the isotropic Gaussian kernel, which is valid for \eqn{\mathbb{R}^{d}.}
#' \deqn{a(\tau;\theta) = \exp(-\tau^{2} / \theta).}
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
#' \deqn{a(\tau;\theta) = \exp(-\tau / \theta).}
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

#' Isotropic Wave (Cardinal Sine) Kernel.
#'
#' This computes the isotropic wave (cardinal sine) kernel, which is valid for \eqn{\mathbb{R}^{3}.}
#' \deqn{a(\tau;\theta) = \begin{cases}
#' \frac{\theta}{\tau} \sin\left( \frac{\tau}{\theta} \right), & \tau \neq = 0 \\
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
#' kernel_wave(x, 1)
#' x <- c(1, 2, 3)
#' kernel_wave(x, 0.2)
kernel_wave <- function(x, theta) {
  stopifnot(theta > 0, length(x) >= 1, all(x >= 0))
  return(sapply(x, function(t) ifelse(t == 0, 1, (theta / t) * sin(t / theta))))
}

#' Isotropic Rational Quadratic Kernel.
#'
#' This computes the isotropic rational quadratic kernel, which is valid for \eqn{\mathbb{R}^{d}.}
#' \deqn{a(\tau;\theta) = 1 - \frac{\tau^{2}}{\tau^{2} + \theta}.}
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
#' kernel_rational_quadratic(x, 1)
#' x <- c(1, 2, 3)
#' kernel_rational_quadratic(x, 0.2)
kernel_rational_quadratic <- function(x, theta) {
  stopifnot(theta > 0, length(x) >= 1, all(x >= 0))
  return(1 - (x^2 / (x^2 + theta)))
}

#' Isotropic Spherical Kernel.
#'
#' This computes the isotropic spherical kernel, which is valid for \eqn{\mathbb{R}^{3}.}
#' \deqn{a(\tau;\theta) =
#' \begin{cases}
#' 1 - \frac{3}{2}\frac{\tau}{\theta} + \frac{1}{2}\left( \frac{\tau}{\theta} \right)^{3}, & \tau < \theta \\
#' 0 & \text{otherwise}
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
#' kernel_spherical(x, 2)
#' x <- c(1, 2, 3)
#' kernel_spherical(x, 2)
kernel_spherical <- function(x, theta) {
  stopifnot(theta > 0, length(x) >= 1, all(x >= 0))
  return(sapply(x, function(t) ifelse(t < theta, 1 - ((3/2) * (t / theta)) + ((1/2) * (t / theta)^3), 0)))
}

#' Isotropic Circular Kernel.
#'
#' This computes the isotropic circular kernel, which is valid for \eqn{\mathbb{R}^{2}.}
#' \deqn{a(\tau;\theta) =
#' \begin{cases}
#' \frac{2}{\pi}\arccos\left( \frac{\tau}{\theta} \right) - \frac{2}{\pi}\frac{\tau}{\theta} \sqrt{ 1 - \left( \frac{\tau}{\theta} \right)^{2} }, & \tau < \theta \\
#' 0 & \text{otherwise}
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
#' kernel_circular(x, 2)
#' x <- c(1, 2, 3)
#' kernel_circular(x, 2)
kernel_circular <- function(x, theta) {
  stopifnot(theta > 0, length(x) >= 1, all(x >= 0))
  return(sapply(x, function(t) ifelse(t < theta, ((2 / pi) * acos(t / theta)) - ((2 / pi) * (t / theta) * sqrt(1 - (t / theta)^2)), 0)))
}

#' Isotropic Bessel Kernel.
#'
#' This computes the isotropic Bessel kernel, which is valid for \eqn{\mathbb{R}^{d}} for \eqn{\nu \geq \frac{d}{2} - 1.}
#' \deqn{a(\tau; \nu) = 2^{\nu} \Gamma(\nu + 1) J_{\nu}(\tau) \tau^{-\nu} ,}
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
#' kernel_bessel_j(x, 1, -0.5, 1)
#' x <- c(1, 2, 3)
#' kernel_bessel_j(x, 1, 2, 2)
kernel_bessel_j <- function(x, theta, nu, dim) {
  stopifnot(theta > 0, nu >= (dim/2) - 1, length(x) >= 1, all(x >= 0))
  return(sapply(x, function(t) ifelse(t == 0, 1, (2^nu) * gamma(nu + 1) * (besselJ(t / theta, nu) / ((t / theta)^nu)))))
}

#' Isotropic Matérn Kernel.
#'
#' This computes the isotropic Matérn kernel, which is valid for \eqn{\mathbb{R}^{d}} for \eqn{\nu > 0.}
#' \deqn{a(\tau; \nu) = \left(\sqrt{2\nu} \frac{\tau}{\theta} \right)^{\nu} \left(2^{\nu - 1} \Gamma(\nu) \right)^{-1} K_{\nu}\left( \sqrt{2\nu} \frac{\tau}{\theta} \right) ,}
#' where \eqn{K_{\nu}(\cdot)} is the modified Bessel function of the second kind.
#'
#' @references
#' Table 4.2 of Hristopulos, D. T. (2020). Random Fields for Spatial Data Modeling:  A Primer for Scientists and Engineers. In Advances in Geographic Information Science. Springer Netherlands. https://doi.org/10.1007/978-94-024-1918-4
#' Genton, Marc. (2001). Classes of Kernels for Machine Learning: A Statistics Perspective. Journal of Machine Learning Research. 2. 299-312. 10.1162/15324430260185646.
#'
#' @param x A vector of values of at least length 1.
#' @param theta Scale parameter, \eqn{\theta > 0.}
#' @param nu A parameter \eqn{\nu > 0.}
#'
#' @return A vector of values.
#' @export
#'
#' @examples
#' x <- 1
#' kernel_matern(x, 2, 1)
#' x <- c(1, 2, 3)
#' kernel_matern(x, 0.5, 2)
kernel_matern <- function(x, theta, nu) {
  stopifnot(theta > 0, nu > 0, length(x) >= 1, all(x >= 0))
  return(sapply(x, function(t) ifelse(t == 0, 1, (((sqrt(2 * nu) * t / theta)^nu) / (2^(nu - 1) * gamma(nu))) * besselK(sqrt(2 * nu) * t / theta, nu))))
}

#' Isotropic Cauchy Kernel.
#'
#' This computes the isotropic Cauchy kernel, which is valid for \eqn{\mathbb{R}^{d}} for \eqn{0 < \alpha \leq 2} and \eqn{\beta \geq 0.}
#' \deqn{a(\tau ; \theta, \alpha, \beta) = (1 + (\tau / \theta)^{\alpha})^{-(\beta / \alpha)} .}
#'
#' @references
#' Table 4.2 of Hristopulos, D. T. (2020). Random Fields for Spatial Data Modeling:  A Primer for Scientists and Engineers. In Advances in Geographic Information Science. Springer Netherlands. https://doi.org/10.1007/978-94-024-1918-4
#'
#' @param x A vector of values of at least length 1.
#' @param theta Scale parameter, \eqn{\theta > 0.}
#' @param alpha A parameter for the kernel, \eqn{0 < a \leq 2.}
#' @param beta A parameter for the kernel, \eqn{b > 0.}
#'
#' @return A vector of values.
#' @export
#'
#' @examples
#' x <- 1
#' kernel_cauchy(x, 2, 1, 1)
#' x <- c(1, 2, 3)
#' kernel_cauchy(x, 0, 2, 3)
kernel_cauchy <- function(x, theta, alpha, beta) {
  stopifnot(theta > 0, alpha > 0, alpha <= 2, beta >= 0, length(x) >= 1, all(x >= 0))
  return((1 + (x / theta)^alpha)^(-(beta / alpha)))
}
