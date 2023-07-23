#' @title
#' Get a unity vector of a specific legth
#'
#' @param n a number of elements.
#'
#' @return
#' A column of `n` unities.
#'
#' @keywords internal
unity <- function(n) {
  return(matrix(rep(1, n), nrow = n, ncol = 1))
}
