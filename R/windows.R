# The function in this file are intended to be used with the tapered estimator.
# Despite being called windows, they are 'half' of the window.
# At 0, the window has a value of 0, at 1 it has a value of 1.

.tukey_window <- function(x) {
  return((1/2) * (1 - cos(pi * x)))
}

.triangular_window <- function(x) {
  return(x)
}

.sine_window <- function(x) {
  return(sin(pi * x / 2))
}

.power_sine_window <- function(x, params = 1) {
  return((sin(pi * x  / 2))^params[1])
}

.blackman_window <- function(x, params = 0.16) {
  return(((1 - params[1]) / 2)  - (1 / 2) * cos(pi * x) + (params[1] / 2) * cos(2 * pi * x))
}

.hann_poisson_window <- function(x, params = 1) {
  (1 / 2) * (1 - cos(pi * x)) * exp(- (params[1] * abs(1 - x)))
}

.welch_window <- function(x) {
  return(1 - (x - 1)^2)
}

.validate_window_params <- function(name = c("tukey", "triangular", "sine", "power_sine", "blackman", "hann_poisson", "welch"), params) {
  name <- match.arg(name)

  switch(name,
         tukey = TRUE,
         triangular = TRUE,
         sine = TRUE,
         power_sine = ifelse(params[1] > 0, TRUE, FALSE),
         blackman = ifelse(abs(params[1]) <= 0.25, TRUE, FALSE),
         hann_poisson = ifelse(params[1] > 0, TRUE, FALSE),
         welch = TRUE
  )
}

#' 1D Window Functions.
#'
#' A window function in this context is a continuous nondecreasing function such that at 0 it is 0, and at 1, it is 1.
#' Note that `window()` is **deprecated**, please use `window_ec()` instead.
#' This computes one of the window functions listed below.
#'
#' \strong{Tukey Window}.
#' The Tukey window is defined as
#' \deqn{w(x) = \frac{1}{2} - \frac{1}{2} \cos(\pi x) , x \in [0, 1].}
#' The \code{params} argument is empty.
#'
#' \strong{Triangular Window}.
#' The triangular window is given by
#' \deqn{w(x) = x, x \in [0, 1].}
#' The \code{params} argument is empty.
#'
#' \strong{Sine Window}.
#' The sine window is given by
#' \deqn{w(x) = \sin\left(\pi x / 2 \right), x \in [0, 1].}
#' The \code{params} argument is empty.
#'
#' \strong{Power Sine Window}.
#' The power sine window is given by
#' \deqn{w(x; a) = \sin^{a}(\pi x / 2), x \in [0, 1], a > 0.}
#' The \code{params} argument is of the form \code{c(}\eqn{a}\code{)}.
#'
#' \strong{Blackman Window}.
#' The Blackman window is defined as
#' \deqn{w(x; a) = ( (1 - a) / 2) - \frac{1}{2} \cos(\pi x) + \frac{a}{2} \cos(2 \pi x), x \in [0, 1], a \in [-0.25, 0.25] .}
#' The \code{params} argument is of the form \code{c(}\eqn{a}\code{)}.
#' It is recommended that \eqn{a \in [-0.25, 0.25]} to ensure that the window is nondecreasi#'ng on \eqn{[0, 1].}
#'
#' \strong{Hann-Poisson Window}.
#' The Hann-Poisson window is defined as
#' \deqn{w(x; a) = \frac{1}{2} (1 - \cos(\pi x)) \exp( - (a \left|1 - x \right|) ) , x \in [0, 1], a > 0.}
#' The \code{params} argument is of the form \code{c(}\eqn{a}\code{)}.
#'
#' \strong{Welch Window}.
#' The Welch window is given by
#' \deqn{w(x) = 1 - (x - 1)^2 , x \in [0, 1] . }
#' The \code{params} argument is empty.
#'
#' See the function call examples below.
#'
#' @param x A vector or matrix of arguments of at least length 1. Each value must be between 0 and 1, inclusive.
#' @param name The name of the window. Options are: tukey, triangular, sine, power_sine, blackman, hann_poisson, welch.
#' @param params A vector of parameters for the windows. See the documentation below for the position of the parameters.
#'
#' @return A vector or matrix of values.
#' @export
#'
#' @rdname window_ec
#'
#' @examples
#' x <- c(0.2, 0.4, 0.6)
#' window_ec(x, "tukey")
#' window_ec(x, "triangular")
#' window_ec(x, "sine")
#' window_ec(x, "power_sine", c(0.7))
#' window_ec(x, "blackman", c(0.16))
#' window_ec(x, "hann_poisson", c(0.7))
#' window_ec(x, "welch")
#' curve(window_ec(x, "tukey"), from = 0, to = 1)
#' curve(window_ec(x, "triangular"), from = 0, to = 1)
#' curve(window_ec(x, "sine"), from = 0, to = 1)
#' curve(window_ec(x, "power_sine", c(0.7)), from = 0, to = 1)
#' curve(window_ec(x, "blackman", c(0.16)), from = 0, to = 1)
#' curve(window_ec(x, "hann_poisson", c(0.7)), from = 0, to = 1)
#' curve(window_ec(x, "welch"), from = 0, to = 1)
window_ec <- function(x, name = c("tukey", "triangular", "sine", "power_sine", "blackman", "hann_poisson", "welch"), params=c(1)) {
  stopifnot(length(x) >= 1, is.numeric(x), all(x >= 0), all(x <= 1), length(params) > 0, is.numeric(params[1]))

  name <- match.arg(name)

  if(!.validate_window_params(name, params)) {
    stop("Parameter choice `(", paste0(params, collapse = ", "), ")` is not appropriate for `", name, "` window.")
  }

  returnVal <- switch(name,
         tukey = .tukey_window(x),
         triangular = .triangular_window(x),
         sine = .sine_window(x),
         power_sine = .power_sine_window(x, params),
         blackman = .blackman_window(x, params),
         hann_poisson = .hann_poisson_window(x, params),
         welch = .welch_window(x)
  )

  return(returnVal)
}

#' 1D Symmetric Window Functions.
#'
#' A symmetric window function in this context are traditional window functions, unlike the [window] functions.
#' Note that `window_symm()` is **deprecated**, please use `window_symm_ec()` instead.
#' This computes one of the symmetric window functions listed below, all of which are defined for \eqn{x \in [-1, 1],} and are 0 otherwise.
#'
#' \strong{Tukey Window}.
#' The Tukey window is defined as
#' \deqn{w(x) = \frac{1}{2} + \frac{1}{2} \cos(\pi |x|) , x \in [-1, 1].}
#' The \code{params} argument is empty, see the example.
#'
#' \strong{Triangular Window}.
#' The triangular window is given by
#' \deqn{w(x) = 1 - |x|, x \in [-1, 1].}
#' The \code{params} argument is empty, see the example.
#'
#' \strong{Sine Window}.
#' The sine window is given by
#' \deqn{w(x) = 1 - \sin\left(\pi |x| / 2 \right), x \in [-1, 1].}
#' The \code{params} argument is empty, see the example.
#'
#' \strong{Power Sine Window}.
#' The power sine window is given by
#' \deqn{w(x; a) = 1 - \sin^{a}(\pi |x| / 2), x \in [-1, 1], a > 0.}
#' The \code{params} argument is of the form \code{c(}\eqn{a}\code{)}
#'
#' \strong{Blackman Window}.
#' The Blackman window is defined as
#' \deqn{w(x; a) = 1 +  ( (a - 1) / 2) + \frac{1}{2} \cos(\pi |x|) - \frac{a}{2} \cos(2 \pi |x|), x \in [-1, 1], a \in [-0.25, 0.25] .}
#' The \code{params} argument is of the form \code{c(}\eqn{a}\code{)}.
#' The standard value of \eqn{a} for this window is \eqn{0.16.}
#'
#' \strong{Hann-Poisson Window}.
#' The Hann-Poisson window is defined as
#' \deqn{w(x; a) = 1 - \frac{1}{2} (1 - \cos(\pi |x|)) \exp( - (a \left|1 - |x| \right|) ) , x \in [-1, 1], a > 0 .}
#' The \code{params} argument is of the form \code{c(}\eqn{a}\code{)}
#'
#' \strong{Welch Window}.
#' The Welch window is given by
#' \deqn{w(x) = (|x| - 1)^2 , x \in [-1, 1] . }
#' The \code{params} argument is empty, see the example.
#'
#' @param x A vector or matrix of arguments of at least length 1. Each value must be between 0 and 1, inclusive.
#' @param name The name of the window. Options are: tukey, triangular, sine, power_sine, blackman, hann_poisson, welch.
#' @param params A vector of parameters for the windows. See the documentation below for the position of the parameters.
#'
#' @return A vector or matrix of values.
#' @export
#'
#' @rdname window_symm_ec
#'
#' @examples
#' x <- c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)
#' window_symm_ec(x, "tukey")
#' window_symm_ec(x, "triangular")
#' window_symm_ec(x, "sine")
#' window_symm_ec(x, "power_sine", c(0.7))
#' window_symm_ec(x, "blackman", c(0.16))
#' window_symm_ec(x, "hann_poisson", c(0.7))
#' window_symm_ec(x, "welch")
#' curve(window_symm_ec(x, "tukey"), from = -1, to = 1)
#' curve(window_symm_ec(x, "triangular"), from = -1, to = 1)
#' curve(window_symm_ec(x, "sine"), from = -1, to = 1)
#' curve(window_symm_ec(x, "power_sine", c(0.7)), from = -1, to = 1)
#' curve(window_symm_ec(x, "blackman", c(0.16)), from = -1, to = 1)
#' curve(window_symm_ec(x, "hann_poisson", c(0.7)), from = -1, to = 1)
#' curve(window_symm_ec(x, "welch"), from = -1, to = 1)
window_symm_ec <- function(x, name = c("tukey", "triangular", "sine", "power_sine", "blackman", "hann_poisson", "welch"), params=c(1)) {
  stopifnot(length(x) >= 1, is.numeric(x), !any(is.na(x)), length(params) > 0, is.numeric(params[1]))

  name <- match.arg(name)

  if(!.validate_window_params(name, params)) {
    stop("Parameter choice `(", paste0(params, collapse = ", "), ")` is not appropriate for `", name, "` symmetric window.")
  }

  returnVal <- switch(name,
         tukey = 1 - .tukey_window(abs(x)),
         triangular = 1 - .triangular_window(abs(x)),
         sine = 1 - .sine_window(abs(x)),
         power_sine = 1 - .power_sine_window(abs(x), params),
         blackman = 1 - .blackman_window(abs(x), params),
         hann_poisson = 1 - .hann_poisson_window(abs(x), params),
         welch = 1 - .welch_window(abs(x))
  )
  returnVal[abs(x) > 1] <- 0

  return(returnVal)
}
