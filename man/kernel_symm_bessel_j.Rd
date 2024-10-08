% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symmetric_kernels.R
\name{kernel_symm_bessel_j}
\alias{kernel_symm_bessel_j}
\title{Symmetric Bessel Kernel.}
\usage{
kernel_symm_bessel_j(x, theta, nu, dim)
}
\arguments{
\item{x}{A vector of values of at least length 1.}

\item{theta}{Scale parameter, \eqn{\theta > 0.}}

\item{nu}{Parameter defining strength of periodicity \eqn{\nu \geq \frac{\text{dim}}{2} - 1.}}

\item{dim}{The dimension of the isotropic kernel.}
}
\value{
A vector of values.
}
\description{
This computes the Bessel kernel. It has the properties of a probability distribution, unlike \code{kernel_bessel_j(x, theta)}.
\deqn{a(\tau; \nu) = (\Gamma(\frac{1}{2} + \nu)/(2 \sqrt{\pi} \theta \Gamma(1 + \nu))) ( 2^{\nu} \Gamma(\nu + 1) J_{\nu}(\tau) \tau^{-\nu}) ,}
where \eqn{J_{\nu}(\cdot)} is the Bessel function of the first kind.
}
\examples{
x <- 1
kernel_symm_bessel_j(x, 1, -0.5, 1)
x <- c(-3, -2, -1, 0, 1, 2, 3)
kernel_symm_bessel_j(x, 1, 2, 2)
}
\references{
Table 4.2 of Hristopulos, D. T. (2020). Random Fields for Spatial Data Modeling:  A Primer for Scientists and Engineers. In Advances in Geographic Information Science. Springer Netherlands. https://doi.org/10.1007/978-94-024-1918-4
}
