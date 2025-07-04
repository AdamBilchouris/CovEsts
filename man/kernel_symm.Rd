% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symmetric_kernels.R
\name{kernel_symm}
\alias{kernel_symm}
\title{1D Isotropic Symmetric Kernels.}
\usage{
kernel_symm(x, name, params = c(1))
}
\arguments{
\item{x}{A vector or matrix of arguments of at least length 1 for which the kernel is computed at. Each value can be negative as well as positive.}

\item{name}{The name of the kernel. Options are: gaussian, wave, rational_quadratic, bessel_j.}

\item{params}{A vector of parameters for the kernel. See the documentation below for the position of the parameters. All kernels will have a scale parameter as the first value in the vector.}
}
\value{
A vector or matrix of values.
}
\description{
These functions computes values of kernels that have the properties of symmetric probability distributions.
For a kernel \eqn{a(x)}, the standardised version is \eqn{a(x) / \int_{-\infty}^{\infty} a(x) dx}, so that the integral is 1.
The symmetric kernels are different to \link{kernel} and are used in the functions \link{compute_adjusted_est} and \link{compute_truncated_est}.
}
\details{
\strong{Symmetric Gaussian Kernel}.
The symmetric Gaussian kernel is defined as
\deqn{a(x;\theta) = \sqrt{\pi \theta} \exp(-x^{2} / \theta), \theta > 0.}
The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.

\strong{Symmetric Wave Kernel}.
The wave (cardinal sine) kernel is given by
\deqn{a(x;\theta) = \left\{ \begin{array}{ll}
(\sqrt{\theta^{2}} \pi)^{-1} \frac{\theta}{x} \sin\left( \frac{x}{\theta} \right), & x \neq 0 \\
1, & x = 0
\end{array} ,
\theta > 0 . \right.}
The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}

\strong{Symmetric Rational Quadratic Kernel}.
The symmetric rational quadratic kernel is given by
\deqn{a(x;\theta) = (\pi \sqrt{\theta})^{-1} (1 - \frac{x^{2}}{x^{2} + \theta}), \theta > 0.}
The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}

\strong{Symmetric Besesel Kernel}.
The symmetric Bessel kernel, which is valid when \eqn{\nu \geq \frac{d}{2} - 1,} is given by
\deqn{a(x; \theta, \nu) = (\Gamma(\frac{1}{2} + \nu)/(2 \sqrt{\pi} \theta \Gamma(1 + \nu))) ( 2^{\nu} \Gamma(\nu + 1) J_{\nu}(x / \theta) (x / \theta)^{-\nu}), \theta > 0, \nu \geq \frac{d}{2} - 1.}
where \eqn{J_{\nu}(\cdot)} is the Bessel function of the first kind and \eqn{d} is the dimension.
The \code{params} argument is of the form \code{c(}\eqn{\theta, \nu, d}\code{)}.
}
\examples{
x <- c(-2, -1, 0, 1, 2)
theta <- 1
kernel_symm(x, "gaussian", c(theta))
kernel_symm(x, "wave", c(theta))
kernel_symm(x, "rational_quadratic", c(theta))
dim <- 1
nu <- 1
kernel_symm(x, "bessel_j", c(theta, nu, dim))
curve(kernel_symm(x, "gaussian", c(theta)), from = -5, to = 5)
curve(kernel_symm(x, "wave", c(theta)), from = -5, to = 5)
curve(kernel_symm(x, "rational_quadratic", c(theta)), from = -5, to = 5)
curve(kernel_symm(x, "bessel_j", c(theta, nu, dim)), from = -5, to = 5)
}
