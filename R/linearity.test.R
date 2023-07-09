#' @importFrom MASS ginv
#' @importFrom matrixcalc matrix.trace
#' @importFrom matrixcalc vec
#'
#' @export
linearity.test <- function(model, J = 1) {
    k  <- model$dim$k
    m  <- model$dim$m
    N  <- model$dim$N
    Nx <- ncol(model$data$X)

    Xn <- NULL
    for (j in 0:J) {
        Xn <- cbind(Xn,
                    (model$data$S^j %x% t(unity(Nx))) * model$data$X)
    }

    vX0  <- model$data$X %x% diag(k)
    vXn <- Xn %x% diag(k)

    vE0 <- vec(t(model$data$Y)) -
        vX0 %*% ginv(t(vX0) %*% vX0) %*% t(vX0) %*% vec(t(model$data$Y))
    E0 <- t(matrix(vE0, ncol = N))
    SSR0 <- t(E0) %*% E0

    vEn <- vE0 -
        vXn %*% ginv(t(vXn) %*% vXn) %*% t(vXn) %*% vE0
    En <- t(matrix(vEn, ncol = N))
    SSRn <- t(En) %*% En

    LM.stat <- N * (k - matrix.trace(ginv(SSR0) %*% SSRn))

    result <- list(LM.stat = LM.stat,
                   crit.10 = qchisq(.90, df = J * k * Nx),
                   crit.5  = qchisq(.95, df = J * k * Nx),
                   crit.1  = qchisq(.99, df = J * k * Nx),
                   p.value = 1 - pchisq(LM.stat, df = J * k * Nx))

    return(result)
}
