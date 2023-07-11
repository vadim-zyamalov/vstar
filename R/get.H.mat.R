#' @title
#' Get the sequence of \eqn{H_{ij}} matrices
#'
#' @description
#' An auxiliary function returning a sequence of \eqn{H_i} matrices
#' used in tests from (Teräsvirta and Yang, 2014).
#'
#' The length of the sequence is exactly the number of elements in the
#' matrix of estimated coefficients \eqn{\hat{B}}.
#' Each matrix \eqn{H_{ij}} is a zero matrix with unity at the position
#' \eqn{(i, j)}. As we don't care about the ordering of these matrices
#' we order them in the order corresponding with R's inner representation
#' of matrices (i.e. up-to-down left-to-right).
#'
#' @param B.mat a matrix \eqn{\hat{B}}.
#'
#' @return
#' A list containing the sequence of \eqn{H_{ij}} matrices.
#'
#' @references
#' T. Teräsvirta and Y. Yang,
#' “Linearity and misspecification tests for vector smooth transition
#' regression models,”
#' Department of Economics and Business Economics, Aarhus University,
#' 2014–04, Feb. 2014.
#'
#' @keywords internal
get.H.mat <- function(B.mat) {
    result <- list()
    for (i in seq_len(length(B.mat))) {
        elems <- rep(0, length(B.mat))
        elems[i] <- 1
        result[[i]] <- matrix(elems, ncol = ncol(B.mat))
    }
    return(result)
}
