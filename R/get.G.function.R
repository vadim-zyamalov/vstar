#' @title
#' Get transition function and its derivative
#' @order 1
#'
#' @description
#' Auxiliary functions returning transition function and its derivative.
#'
#' @param g.function a transition function to be used
#' * "L": Logistic transition function
#' \deqn{\frac{1}{1 + \exp (-\gamma (s_t - c))}}
#' * "E": Exponential transition function
#' \deqn{1 - \exp(-\gamma (s_t - c)^2)}
#' * fn: Custom function receiving transition variable, gammas and threshold
#' values in exactly this order.
#'
#' If any custom function is provided to the `get.G.function` then
#' it's returned as-is.
#'
#' @keywords internal
get.G.function <- function(g.function = "L") {
  if (is.character(g.function)) {
    if (g.function == "L") {
      func <- function(s, g, thr, sd = 1) {
        1 / (1 + exp(-g * (s - thr) / sd))
      }
    } else if (g.function == "E") {
      func <- function(s, g, thr, sd = 1) {
        1 - exp(-g * ((s - thr) / sd)^2)
      }
    } else {
      stop(paste("G: Unknown function \"", g.function, "\"!", sep = ""))
    }
  } else if (is.function(g.function)) {
    func <- g.function
  } else {
    stop("G: Unknown transition function!")
  }

  return(func)
}

#' @rdname get.G.function
#' @order 2
#' @keywords internal
get.G.derivative <- function(g.function = "L") {
  G.func <- get.G.function(g.function)
  if (is.character(g.function)) {
    if (g.function == "L") {
      g.func <- function(s, g, thr, sd = 1) {
        ((s - thr) / sd) * G.func(s, g, thr) * (1 - G.func(s, g, thr))
      }
      c.func <- function(s, g, thr, sd = 1) {
        -g * G.func(s, g, thr) * (1 - G.func(s, g, thr))
      }
    } else if (g.function == "E") {
      g.func <- function(s, g, thr, sd = 1) {
        ((s - thr) / sd)^2 * (1 - G.func(s, g, thr))
      }
      c.func <- function(s, g, thr, sd = 1) {
        -g * ((s - thr) / sd) * (1 - G.func(s, g, thr))
      }
    } else {
      stop(paste("G: Unknown function \"", g.function, "\"!", sep = ""))
    }
  } else if (is.function(g.function)) {
    stop("G: Provide array of derivative functions!")
  } else {
    stop("G: Unknown transition function!")
  }

  return(list(g = g.func, c = c.func))
}
