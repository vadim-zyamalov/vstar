#' @title
#' Modfied `lag` function
#'
#' @description
#' Inspired by the function `lagn` from GAUSS.
#' Adds the value of `na` to the beginning or the end of the result
#' to keep the initial size of the array.
#'
#' @param x a vector or matrix with the values to be lagged
#' If a vector is provided it will be converted to a column matrix.
#' @param i an integer indicating which lag to use.
#' @param na a value to fill missing points.
#'
#' @return a vector (!) with the computed lags.
#'
#' @keywords internal
lagn <- function(x,
                 i,
                 na = NA) {
    if (!is.matrix(x)) x <- as.matrix(x)
    n.obs <- nrow(x)
    n.var <- ncol(x)
    if (i > 0) {
        return(
            rbind(
                matrix(data = na, nrow = i, ncol = n.var),
                x[1:(n.obs - i), , drop = FALSE]
            )
        )
    } else {
        return(
            rbind(
                x[(1 + abs(i)):n.obs, , drop = FALSE],
                matrix(data = na, nrow = abs(i), ncol = n.var)
            )
        )
    }
}
