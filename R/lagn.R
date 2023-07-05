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
