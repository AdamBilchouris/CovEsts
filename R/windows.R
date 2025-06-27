# The function in this file are intended to be used with the tapered estimator.
# Despite being called windows, they are 'half' of the window.
# At 0, the window has a value of 0, at 1 it has a value of 1.

#' 1D Window Functions.
#'
#' A window function in this context is a continuous nondecreasing function such that at 0 it is 0, and at 1, it is 1.
#' This computes one of the window functions listed below.
#'
#' \strong{Tukey Window}.
#' The Tukey window is defined as
#' \deqn{w(x) = \frac{1}{2} - \frac{1}{2} \cos(\pi x) , x \in [0, 1].}
#' The \code{params} argument is empty, see the example.
#'
#' \strong{Triangular Window}.
#' The triangular window is given by
#' \deqn{w(x) = x, x \in [0, 1].}
#' The \code{params} argument is empty, see the example.
#'
#' \strong{Sine Window}.
#' The sine window is given by
#' \deqn{w(x) = \sin\left(\pi x / 2 \right), x \in [0, 1].}
#' The \code{params} argument is empty, see the example.
#'
#' \strong{Power Sine Window}.
#' The power sine window is given by
#' \deqn{w(x; a) = \sin^{a}(\pi x / 2), x \in [0, 1], a > 0.}
#' The \code{params} argument is of the form \code{c(}\eqn{a}\code{)}
#'
#' \strong{Blackman Window}.
#' The Blackman window is defined as
#' \deqn{w(x; a) = ( (1 - a) / 2) - \frac{1}{2} \cos(\pi x) + \frac{a}{2} \cos(2 \pi x), x \in [0, 1], a \in {R} .}
#' The \code{params} argument is of the form \code{c(}\eqn{a}\code{)}
#' It is recommended that \eqn{a \in [-0.25, 0.25]} to ensure that the window is nondecreasing on \eqn{[0, 1].}
#'
#' \strong{Hann-Poisson Window}.
#' The Hann-Poisson window is defined as
#' \deqn{w(x; a) = \frac{1}{2} (1 - \cos(\pi x)) \exp( - (a \left|1 - x \right|) ) , x \in [0, 1], a \in {R} .}
#' The \code{params} argument is of the form \code{c(}\eqn{a}\code{)}
#'
#' \strong{Welch Window}.
#' The Welch window is given by
#' \deqn{w(x) = 1 - (x - 1)^2 , x \in [0, 1] . }
#' The \code{params} argument is empty, see the example.
#'
#' @param x A vector or matrix of arguments of at least length 1. Each value must be between 0 and 1, inclusive.
#' @param name The name of the window. Options are: tukey, triangular, sine, power_sine, blackman, hann_poisson, welch.
#' @param params A vector of parameters for the windows. See the documentation below for the position of the parameters.
#'
#' @return A vector or matrix of values.
#' @export
#'
#' @examples
#' x <- c(0.2, 0.4, 0.6)
#' window(x, "tukey")
#' window(x, "triangular")
#' window(x, "sine")
#' window(x, "power_sine", c(0.7))
#' window(x, "blackman", c(0.16))
#' window(x, "hann_poisson", c(0.7))
#' window(x, "welch")
#' curve(window(x, "tukey"), from = 0, to = 1)
#' curve(window(x, "triangular"), from = 0, to = 1)
#' curve(window(x, "sine"), from = 0, to = 1)
#' curve(window(x, "power_sine", c(0.7)), from = 0, to = 1)
#' curve(window(x, "blackman", c(0.16)), from = 0, to = 1)
#' curve(window(x, "hann_poisson", c(0.7)), from = 0, to = 1)
#' curve(window(x, "welch"), from = 0, to = 1)
window <- function(x, name, params=c(1)) {
  stopifnot(length(x) >= 1, all(x >= 0), all(x <= 1), length(params) > 0, is.numeric(params[1]))
  stopifnot(name %in% c("tukey", "triangular", "sine", "power_sine",
                        "blackman", "hann_poisson", "welch"))

  # params is lenght 1 by default.
  if(length(params) == 1) {
    if(name == "tukey") {
      return((1/2) * (1- cos(pi * x)))
    }
    else if(name == "triangular") {
      return(x)
    }
    else if(name == "sine") {
      return(sin(pi * x / 2))
    }
    else if(name == "power_sine") {
      stopifnot(params[1] > 0)
      return((sin(pi * x / 2)^(params[1])))
    }
    else if(name == "blackman") {
      stopifnot(params[1] >= -0.25, params[1] <= 0.25)
      return(((1 - params[1]) / 2)  - (1 / 2) * cos(pi*x) + (params[1] / 2) * cos(2*pi*x))
    }
    else if(name == "hann_poisson") {
      stopifnot(params[1] > 0)
      return((1 / 2) * (1 - cos(pi * x)) * exp(- (params[1] * abs(1 - x))))
    }
    else if(name == "welch") {
      return(1 - (x-1)^2)
    }
    stop(paste0("Unknown kernel: ", name))
  }
  stop(paste0("length(params) is not 1."))
}
