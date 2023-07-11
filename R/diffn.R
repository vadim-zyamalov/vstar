#' @title
#' Modfied `diff` function
#'
#' @description
#' Inspired by the function `diffn` from GAUSS.
#' Adds the value of `na` to the beginning of the result
#' to keep the initial size of the array.
#'
#' @param x a vector or matrix with the values to be defferenced.
#' If a vector is provided it will be converted to a column matrix.
#' @param lag an integer indicating which lag to use.
#' @param differences an integer indicating the order of the difference.
#' @param na a value to fill missing points.
#' @param ... further arguments to be passed to or from methods.
#'
#' @return a vector (!) with the computed differences.
#'
#' @keywords internal
diffn <- function(x,
                  lag = 1,
                  differences = 1,
                  na = NA,
                  ...) {
    if (!is.matrix(x)) x <- as.matrix(x)
    n.obs <- nrow(x)
    n.var <- ncol(x)
    tmp.diff <- diff(x, lag = lag, differences = differences, ...)
    return(
        rbind(
            matrix(
                data = na,
                nrow = n.obs - nrow(tmp.diff),
                ncol = n.var
            ),
            tmp.diff
        )
    )
}
