% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tapered_estimator.R
\name{compute_tapered_est}
\alias{compute_tapered_est}
\title{Compute the tapered autocovariance function over a set of lags.}
\usage{
compute_tapered_est(
  X,
  maxLag,
  rho,
  window_name,
  window_params = c(1),
  custom_window = FALSE,
  type = "autocovariance",
  meanX = mean(X)
)
}
\arguments{
\item{X}{A vector representing observed values of the process.}

\item{maxLag}{The maximum lag to compute the autocovariance function at.}

\item{rho}{A scale parameter in \eqn{(0, 1].}}

\item{window_name}{The name of the \link{window} function to be used. Possible values are:
"tukey", "triangular", "power_sine", "blackman_window", "hann_poisson", "welch". Alternatively, a custom window function can be provided, see the example.}

\item{window_params}{A vector of parameters of the window function.}

\item{custom_window}{If a custom window is to be used or not. Defaults to \code{FALSE}.}

\item{type}{Compute either the 'autocovariance' or 'autocorrelation'. Defaults to 'autocovariance'.}

\item{meanX}{The average value of \code{X}. Defaults to \code{mean(X)}.}
}
\value{
A vector whose values are the tapered autocovariance estimates.
}
\description{
This function computes the tapered autocovariance over a set of lags,
\deqn{\widehat{C}_{N}^{a} (h) = (H_{2, n}(0))^{-1} \sum_{t=1}^{n-h}  (X(t_{i}) - \bar{X} ) ( X(t_{i} + h) - \bar{X} ) [ a((t_{i} - 1/2) / n; \rho) a((t_{i} + h - 1/2) / n; \rho)  ] ,}
where \eqn{a(\cdot)} is a window function, \eqn{\rho \in (0, 1]} is a scale parameter.
For each lag, the tapered autocovariance is computed using the function \link{tapered_cov_single}.
}
\details{
This function computes the tapered autocovariance over a set of lags,
\deqn{\widehat{C}_{N}^{a} (h) = (H_{2, n}(0))^{-1} \sum_{t=1}^{n-h}  (X(t_{i}) - \bar{X} ) ( X(t_{i} + h) - \bar{X} ) [ a((t_{i} - 1/2) / n; \rho) a((t_{i} + h - 1/2) / n; \rho)  ] ,}
where \eqn{a(\cdot)} is a window function, \eqn{\rho \in (0, 1]}  is a scale parameter.
This estimator considers the edge effect during estimation, assigning a lower weight to values closer to the edges and higher weights for observations closer to the middle.
This estimator is positive-definite and asymptotically unbiased.

Internally, this function calls \link{tapered_cov_single} for each lag \eqn{h}.

The set of taper values (\eqn{a((t_{i} - 1/2) / n; \rho)}) are precomputed in order to save time, as they do not vary.
As \eqn{a((t_{i} + h - 1/2) / n; \rho)} is varying with lag \eqn{h}, the set of taper values is truncated to get the desired result.
}
\examples{
X <- c(1, 2, 3)
compute_tapered_est(X, 2, 0.5, "tukey")
}
\references{
Dahlhaus R., Künsch, H. (1987). Edge Effects and Efficient Parameter Estimation for Stationary Random Fields. In Biometrika (Vol. 74 Issue 4, pp. 877-882). 10.1093/biomet/74.4.877
}
