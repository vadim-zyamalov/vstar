#' @title
#' sample values from a given set
#'
#' @param x a data to sample from.
#' @param n a number of values to sample.
#'
#' @return
#' A vector of sampled values.
#'
#' @keywords internal
get.bootstrap <- function(x, n) {
  N <- length(x)
  idx <- sample(1:N, n, replace = TRUE)
  return(x[idx])
}
