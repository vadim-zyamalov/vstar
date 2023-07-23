#' @title
#' Get grid values
#'
#' @description
#' An auxiliary function to get two sequences of values of \eqn{\gamma}
#' and \eqn{c} in grid points.
#'
#' @param g.limits a vector of lower and upper bounds of \eqn{\gamma}.
#' @param c.limits a vector of lower and upper bounds of \eqn{c}.
#' @param points a number of points in sequences.
#'
#' @return
#' A list containing sequences of values of \eqn{\gamma} and \eqn{c}.
#'
#' @keywords internal
get.grid <- function(g.limits, c.limits, points) {
  return(list(
    g = seq(
      from = min(g.limits),
      to = max(g.limits),
      length.out = points
    ),
    c = seq(
      from = min(c.limits),
      to = max(c.limits),
      length.out = points
    )
  ))
}
