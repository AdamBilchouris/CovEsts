% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corrected_standard_estimator.R
\name{compute_corrected_standard_est}
\alias{compute_corrected_standard_est}
\title{Kernel Correction of the Standard Estimator.}
\usage{
compute_corrected_standard_est(
  X,
  kernel_name,
  kernel_params = c(),
  N_T = 0.1 * length(X),
  meanX = mean(X),
  maxLag = length(X) - 1,
  pd = TRUE,
  custom_kernel = FALSE,
  type = "autocovariance"
)
}
\arguments{
\item{X}{A vector representing observed values of the time series.}

\item{kernel_name}{The name of the \link{kernel} function to be used. Possible values are:
gaussian, exponential, wave, rational_quadratic, spherical, circular, bessel_j, matern, cauchy.}

\item{kernel_params}{A vector of parameters of the kernel function. See \link{kernel} for parameters.
In the case of gaussian, wave, rational_quadratic, spherical and circular, \code{N_T} takes the place of \eqn{\theta}.
For kernels that require parameters other than \eqn{\theta}, such as the Matern kernel, those parameters are passed.}

\item{N_T}{The range at which the kernel function vanishes at. Recommended to be \eqn{0.1 N} when considering all lags. This parameter may be large for a lag small estimation lag.}

\item{meanX}{The average value of \code{X}. Defaults to \code{mean(X)}.}

\item{maxLag}{An optional parameter that determines the maximum lag to compute the estimated autocovariance function at. Defaults to \code{length(X) - 1}.}

\item{pd}{Whether a positive-definite estimate should be used. Defaults to \code{TRUE}.}

\item{custom_kernel}{If a custom kernel is to be used or not. Defaults to \code{FALSE}. See examples.}

\item{type}{Compute either the 'autocovariance' or 'autocorrelation'. Defaults to 'autocovariance'.}
}
\value{
A vector whose values are the kernel corrected autocovariance estimates.
}
\description{
This function computes the standard autocovariance estimator and applies kernel correction to it.
It considers a kernel \eqn{a(\cdot)} which decays or vanishes to zero (depending on the type of kernel) where \eqn{a(0) = 1.}
The rate or value at which the kernel vanishes is \eqn{N_{T}}, which is recommended to be of order \eqn{0.1 N}, where \eqn{N} is the length of the observation window, however, one may need to play with this value.
\deqn{\widehat{C}_{T}^{(a)}(h) = \widehat{C}(h) a_{T}(h), } where \eqn{a_{T}(h) := a(h / N_{T}).}
}
\details{
The aim of this estimator is gradually bring the estimated values to zero through the use of a kernel multiplier. This can be useful when estimating an
autocovariance function that is short-range dependent as estimators can have large fluctuations as the lag increases, or to deal with the wave artefacts for large lags.
This estimator can be positive-definite depending on whether the choice of \eqn{\widehat{C}(\cdot)} and \eqn{a} are chosen to be positive-definite or not.
}
\examples{
X <- c(1, 2, 3)
compute_corrected_standard_est(X, "gaussian")

X <- rnorm(1000)
Y <- c(X[1], X[2])
for(i in 3:length(X)) { Y[i] <- X[i] - 0.3*X[i - 1] - 0.6*X[i - 2] }
plot(Y)
plot(compute_corrected_standard_est(Y, "bessel_j",
     kernel_params=c(0, 1), N_T=0.2*length(Y)))

# Custom kernel
my_kernel <- function(x, theta, params) {
  stopifnot(theta > 0, length(x) >= 1, all(x >= 0))
  return(sapply(x, function(t) ifelse(t == 0, 1,
         ifelse(t == Inf, 0,
         (sin((t^params[1]) / theta) / ((t^params[1]) / theta)) * cos((t^params[2]) / theta)))))
}
plot(compute_corrected_standard_est(Y,
     "my_kernel", kernel_params=c(2, 0.25), custom_kernel = TRUE))
}
\references{
Yaglom, AM (1987). Correlation Theory of Stationary and Related Random Functions. Volume I: Basic Results. Springer New York. 10.1007/978-1-4612-4628-2.
}
