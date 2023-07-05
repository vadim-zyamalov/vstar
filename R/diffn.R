diffn <- function(x,
                  lag = 1,
                  differences = 1,
                  na = NA) {
    if (!is.matrix(x)) x <- as.matrix(x)
    n.obs <- nrow(x)
    n.var <- ncol(x)
    tmp.diff <- diff(x, lag = lag, differences = differences)
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
