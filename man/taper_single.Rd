% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tapered_estimator.R
\name{taper_single}
\alias{taper_single}
\title{Compute Taper}
\usage{
taper_single(x, rho, window_name, window_params = c(1), custom_window = FALSE)
}
\arguments{
\item{x}{A number between 0 and 1 (inclusive).}

\item{rho}{A scale parameter in \eqn{(0, 1].}}

\item{window_name}{The name of the \link{window} function to be used. Possible values are:
tukey, triangular, power_sine, blackman_window, hann_poisson, welch. Alternatively, a custom window function can be provided, see the example.}

\item{window_params}{A vector of parameters of the window function.}

\item{custom_window}{If a custom window is to be used or not. Defaults to \code{FALSE}.}
}
\value{
A single value being the taper function value at x.
}
\description{
This helper function computes the taper function for a given window.
\deqn{a(x; \rho) = \left\{
\begin{array}{ll}
   w(2x/\rho)     & 0 \leq x < \frac{1}{2} \rho, \\
   1              & \frac{1}{2}\rho \leq x \leq \frac{1}{2} \\
   a(1 - x; \rho) & \frac{1}{2} < x \leq 1
\end{array} ,
\right. }
where \eqn{w(\cdot)} is a continuous increasing function with \eqn{w(0)=0, w(1)=1,}
\eqn{\rho \in (0, 1],} and \eqn{x \in [0, 1].} The possible window choices are found in \link{window}.
}
\examples{
x <- 0.4
taper_single(x, 0.5, "tukey")
my_taper <- function(x, ...) {
  return(x)
}
taper_single(x, 0.5, "my_taper", custom_window = TRUE)
}
