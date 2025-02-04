# The function in this file are intended to be used with the tapered estimator.
# Despite being called windows, they are 'half' of the window.
# At 0, the window has a value of 0, at 1 it has a value of 1.

#' 1D window functions.
#'
#' A window function in this context is a continuous nondecreasing function such that at 0 it is 0, and at 1, it is 1.
#' This computes one of the window functions listed below.
#'
#' \strong{Tukey Window}.
#' \deqn{w(x) = \frac{1}{2} \cos(\pi x) , x \in [0, 1].}
#' The \code{params} argument is empty, see the example.
#'
#' \strong{Triangular Window}.
#' \deqn{w(x) = x, x \in [0, 1].}
#' The \code{params} argument is empty, see the example.
#'
#' \strong{Sine Window}.
#' \deqn{w(x) = x, x \in [0, 1].}
#' The \code{params} argument is empty, see the example.
#'
#' \strong{Power Sine Window}.
#' \deqn{w(x) = \sin^{a}(\pi x / 2), x \in [0, 1], a > 0.}
#' The \code{params} argument is of the form \code{c(}\eqn{a}\code{)}
#'
#' \strong{Blackman Window}.
#' \deqn{w(x) = ( (1 - \alpha) / 2) - (1/2) \cos(\pi x) + (\alpha / 2) \cos(2 \pi x), x \in [0, 1], \alpha \in {R} .}
#' The \code{params} argument is of the form \code{c(}\eqn{\alpha}\code{)}
#'
#' \strong{Hann-Poisson Window}.
#' \deqn{w(x) = (1/2) (1 - \cos(\pi x)) \exp( - (\alpha \left|1 - x \right|) ) , x \in [0, 1], \alpha \in {R} .}
#' The \code{params} argument is of the form \code{c(}\eqn{\alpha}\code{)}
#'
#' \strong{Welch Window}.
#' \deqn{w(x) = 1 - (x - 1)^2 , x \in [0, 1] . }
#' The \code{params} argument is empty, see the example.
#'
#' @param x A vector of values of at least length 1.
#' @param name The name of the window. Options are: "tukey", "triangular", "sine", "power_sine", "blackman", "hann_poisson", "welch".
#' @param params A vector of parameters for the windows. See the documentation below for the position of the parameters.
#'
#' @return A vector of values.
#' @export
#'
#' @examples
#' x <- c(0.2, 0.4, 0.6)
#' window(x, "tukey")
#' window(x, "triangular")
#' window(x, "sine")
#' window(x, "power_sine", c(0.7))
#' window(x, "blackman", c(-0.7))
#' window(x, "hann_poisson", c(-0.7))
#' window(x, "welch")
window <- function(x, name, params=c(1)) {
  stopifnot(length(x) >= 1, all(x >= 0), all(x <= 1), is.vector(params), length(params) > 0, is.numeric(params[1]))
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
      return(((1 - params[1]) / 2)  - (1 / 2) * cos(pi*x) + (params[1] / 2) * cos(2*pi*x))
    }
    else if(name == "hann_poisson") {
      return((1 / 2) * (1 - cos(pi * x)) * exp(- (params[1] * abs(1 - x))))
    }
    else if(name == "welch") {
      return(1 - (x-1)^2)
    }
    stop(paste0("Unknown kernel: ", name))
  }
  stop(paste0("length(params) is not 1."))
}
