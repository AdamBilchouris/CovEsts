% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standard_est.R
\name{standard_est_single}
\alias{standard_est_single}
\title{Autocovariance at a single lag.}
\usage{
standard_est_single(X, tau, N = length(X), meanX = mean(X), pd = TRUE)
}
\arguments{
\item{X}{A vector of observations representing a time series.}

\item{tau}{The lag at which to compute the autocovariance function.}

\item{N}{The length of the vector X.}

\item{meanX}{The average value of the observations.}

\item{pd}{Whether a positive definite estimate should be used.}
}
\value{
The autocovariance at lag tau for the time series X.
}
\description{
This computes the following two variants of the autocovariance function depending on
the parameter \code{pd}.
}
\details{
For \code{pd = TRUE}:
\deqn{
\widehat{C}(\tau) = \frac{1}{N} \sum_{i=1}^{N-\tau} ( X(t_{i}) - \bar{X} ) ( X(t_{i} + \tau) - \bar{X} ) .
}

For \code{pd = FALSE}:
\deqn{
\widehat{C}(\tau) = \frac{1}{N - \tau} \sum_{i=1}^{N-\tau} ( X(t_{i}) - \bar{X} ) ( X(t_{i} + \tau) - \bar{X} ) .
}
}
\examples{
X <- c(1, 2, 3)
standard_est_single(X, 1, length(X), mean(X), pd=FALSE)
}
