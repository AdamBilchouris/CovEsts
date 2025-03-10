% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernels.R
\name{kernel}
\alias{kernel}
\title{Isotropic 1D kernels.}
\usage{
kernel(x, name, params = c(1))
}
\arguments{
\item{x}{A vector or matrix of values of at least length 1}

\item{name}{The name of the kernel. Options are: "gaussian", "exponential", "wave", "rational_quadratic", "spherical", "circular", "bessel_j", "matern", and "cauchy".}

\item{params}{A vector of parameters for the kernel. See the documentation below for the position of the parameters. All kernels will have a scale parameter as the first value in the vector.}
}
\value{
A vector of values.
}
\description{
This computes one of the isotropic kernels listed below.
}
\details{
\strong{Gaussian Kernel}.
This computes the isotropic Gaussian kernel, which is valid for \eqn{{R}^{d}.}
\deqn{a(\tau;\theta) = \exp(-\tau^{2} / \theta).}
The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.

\strong{Exponential Kernel}.
This computes the isotropic exponential kernel, which is valid for \eqn{{R}^{d}.}
\deqn{a(\tau;\theta) = \exp(-\tau / \theta).}
The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.

\strong{Isotropic Wave (Cardinal Sine) Kernel}.
This computes the isotropic wave (cardinal sine) kernel, which is valid for \eqn{{R}^{3}.}
\deqn{a(\tau;\theta) = \left\{ \begin{array}{ll}
\frac{\theta}{\tau} \sin\left( \frac{\tau}{\theta} \right), & \tau \neq 0 \\
1, & \tau = 0
\end{array} . \right.
}
The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.

\strong{Isotropic Rational Quadratic Kernel}.
This computes the isotropic rational quadratic kernel, which is valid for \eqn{{R}^{d}.}
\deqn{a(\tau;\theta) = 1 - \frac{\tau^{2}}{\tau^{2} + \theta}.}
The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.

\strong{Isotropic Spherical Kernel}.
This computes the isotropic spherical kernel, which is valid for \eqn{{R}^{3}.}
\deqn{a(\tau;\theta) = \left\{
\begin{array}{ll}
1 - \frac{3}{2}\frac{\tau}{\theta} + \frac{1}{2}\left( \frac{\tau}{\theta} \right)^{3}, & \tau < \theta \\
0, & \mbox{otherwise}
\end{array} . \right.
}
The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.

\strong{Isotropic Circular Kernel}.
This computes the isotropic circular kernel, which is valid for \eqn{{R}^{2}.}
\deqn{a(\tau;\theta) = \left\{
\begin{array}{ll}
\frac{2}{\pi}\arccos\left( \frac{\tau}{\theta} \right) - \frac{2}{\pi}\frac{\tau}{\theta} \sqrt{ 1 - \left( \frac{\tau}{\theta} \right)^{2} }, & \tau < \theta \\
0, & \mbox{otherwise}
\end{array} . \right.
}
The \code{params} argument is of the form \code{c(}\eqn{\theta}\code{)}.

\strong{Isotropic Matérn Kernel}.
This computes the isotropic Matérn kernel, which is valid for \eqn{{R}^{d}} for \eqn{\nu > 0.}
\deqn{a(\tau; \theta, \nu) = \left(\sqrt{2\nu} \frac{\tau}{\theta} \right)^{\nu} \left(2^{\nu - 1} \Gamma(\nu) \right)^{-1} K_{\nu}\left( \sqrt{2\nu} \frac{\tau}{\theta} \right) ,}
where \eqn{K_{\nu}(\cdot)} is the modified Bessel function of the second kind.
The \code{params} argument is of the form \code{c(}\eqn{\theta, \nu}\code{)}.

\strong{Isotropic Bessel Kernel}.
This computes the isotropic Bessel kernel, which is valid for \eqn{{R}^{d}} for \eqn{\nu \geq \frac{d}{2} - 1.}
\deqn{a(\tau; \theta, \nu) = 2^{\nu} \Gamma(\nu + 1) J_{\nu}(\tau / \theta) (\tau / \theta)^{-\nu} ,}
where \eqn{J_{\nu}(\cdot)} is the Bessel function of the first kind.
The \code{params} argument is of the form \code{c(}\eqn{\theta, \nu, d} \code{)}.

\strong{Isotropic Cauchy Kernel}.
This computes the isotropic Cauchy kernel, which is valid for \eqn{{R}^{d}} for \eqn{0 < \alpha \leq 2} and \eqn{\beta \geq 0.}
\deqn{a(\tau ; \theta, \alpha, \beta) = (1 + (\tau / \theta)^{\alpha})^{-(\beta / \alpha)} .}
The \code{params} argument is of the form \code{c(}\eqn{\theta, \alpha, \beta} \code{)}.
}
\examples{
x <- c(0.2, 0.4, 0.6)
theta <- 0.9
kernel(x, "gaussian", c(theta))
kernel(x, "exponential", c(theta))
kernel(x, "wave", c(theta))
kernel(x, "rational_quadratic", c(theta))
kernel(x, "spherical", c(theta))
kernel(x, "circular", c(theta))
nu <- 1
kernel(x, "matern", c(theta, nu))
dim <- 1
kernel(x, "bessel_j", c(theta, nu, dim))
alpha <- 1
beta <- 2
kernel(x, "cauchy", c(theta, alpha, beta))
}
\references{
Genton, Marc. (2001). Classes of Kernels for Machine Learning: A Statistics Perspective. Journal of Machine Learning Research. 2. 299-312. 10.1162/15324430260185646.
Table 4.2 of Hristopulos, D. T. (2020). Random Fields for Spatial Data Modeling:  A Primer for Scientists and Engineers. In Advances in Geographic Information Science. Springer Netherlands. https://doi.org/10.1007/978-94-024-1918-4
}
