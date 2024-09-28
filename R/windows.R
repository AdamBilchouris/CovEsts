# The functions in this file are intended to be used with the tapered estimator.
# Despite being called windows, they are 'half' of the window.
# At 0, the window has a value of 0, at 1 it has a value of 1.

#' Tukey Window.
#'
#' This function computes the Tukey window.
#' \deqn{w(x) = \frac{1}{2} \cos(\pi x) , x \in [0, 1].}
#'
#' @param x A vector of values at which the window is computed at.
#' @param ... Extra unusued arguments.
#'
#' @return A vector of values.
#' @export
#'
#' @examples
#' X <- c(0.1, 0.5, 0.75)
#' tukey_window(X)
tukey_window <- function(x, ...) {
  stopifnot(length(x) >= 1, all(x >= 0), all(x <= 1))
  return((1/2) * (1- cos(pi * x)))
}

#' Triangular Window.
#'
#' This computes the triangular window.
#' \deqn{w(x) = x, x \in [0, 1].}
#'
#' @param x A vector of values at which the window is computed at.
#' @param ... Extra unusued arguments.
#'
#' @return A vector of values.
#' @export
#'
#' @examples
#' X <- c(0.1, 0.5, 0.75)
#' triangular_window(X)
triangular_window <- function(x, ...) {
  stopifnot(length(x) >= 1, all(x >= 0), all(x <= 1))
  return(x)
}

#' Sine Window.
#'
#' This computes the sine window.
#' \deqn{w(x) = \sin(\pi x / 2), x \in [0, 1].}
#'
#' @param x A vector of values at which the window is computed at.
#' @param ... Extra unusued arguments.
#'
#' @return A vector of values.
#' @export
#'
#' @examples
#' X <- c(0.1, 0.5, 0.75)
#' sine_window(X)
sine_window <- function(x, ...) {
  stopifnot(length(x) >= 1, all(x >= 0), all(x <= 1))
  return(sin(pi * x / 2))
}

#' Power Sine Window.
#'
#' This computes the power sine window.
#' \deqn{w(x) = \sin^{a}(\pi x / 2), x \in [0, 1], a > 0.}
#'
#' @param x A vector of values at which the window is computed at.
#' @param a A positive number indicating the power.
#' @param ... Extra unusued arguments.
#'
#' @return A vector of values.
#' @export
#'
#' @examples
#' X <- c(0.1, 0.5, 0.75)
#' power_sine_window(X, 3.2)
power_sine_window <- function(x, a, ...) {
  stopifnot(length(x) >= 1, a > 0, all(x >= 0), all(x <= 1))
  return((sin(pi * x / 2)^(a)))
}

#' Blackman Window.
#'
#' This computes the Blackman window.
#' \deqn{w(x) = ( (1 - \alpha) / 2) - (1/2) \cos(\pi x) + (\alpha / 2) \cos(2 \pi x), x \in [0, 1], \alpha \in \mathbb{R}.}
#' As \eqn{\alpha} increases, the function
#'
#' @param x A vector of values at which the window is computed at.
#' @param alpha A real number which controls the height of the spike. A negative alpha creates a positive peak and a positive alpha creates a negative peak.
#' @param ... Extra unusued arguments.
#'
#' @return A vector of values.
#' @export
#'
#' @examples
#' X <- c(0.1, 0.5, 0.75)
#' blackman_window(X, 1.34)
blackman_window <- function(x, alpha, ...) {
  stopifnot(length(x) >= 1, all(x >= 0), all(x <= 1))
  return(((1 - alpha) / 2)  - (1 / 2) * cos(pi*x) + (alpha / 2) * cos(2*pi*x))
}

#' Hann-Poisson Window.
#'
#' This computes the Hann-Poisson window.
#' \deqn{w(x) = (1/2) (1 - \cos(\pi x)) \exp( - (\alpha \left|1 - x \right|) ) , x \in [0, 1], \alpha \in \mathbb{R} .}
#'
#' @param x A vector of values at which the window is computed at.
#' @param alpha A real number which controls the height of the spike. A negative alpha creates a large positive peak, try to avoid negative values above \eqn{-5} unless otherwise needed.
#' @param ... Extra unusued arguments.
#'
#' @return A vector of values.
#' @export
#'
#' @examples
#' X <- c(0.1, 0.5, 0.75)
#' hann_poisson_window(X, 0.2)
hann_poisson_window <- function(x, alpha, ...) {
  stopifnot(length(x) >= 1, all(x >= 0), all(x <= 1))
  return((1 / 2) * (1 - cos(pi * x)) * exp(- (alpha * abs(1 - x))))
}

#' Welch Window.
#'
#' This computes the Welch  window.
#' \deqn{w(x) = 1 - (x - 1)^2 , x \in [0, 1] . }
#'
#' @param x A vector of values at which the window is computed at.
#' @param ... Extra unusued arguments.
#'
#' @return A vector of values.
#' @export
#'
#' @examples
#' X <- c(0.1, 0.5, 0.75)
#' welch_window(X)
welch_window <- function(x, ...) {
  stopifnot(length(x) >= 1, all(x >= 0), all(x <= 1))
  return(1 - (x-1)^2)
}

