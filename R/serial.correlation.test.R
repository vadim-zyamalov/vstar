#' @importFrom MASS ginv
#' @importFrom matrixcalc matrix.trace
#' @importFrom matrixcalc vec
#'
#' @export
serial.correlation.test <- function(model, J = 1, ortogonalize = FALSE) {
    k  <- model$dim$k
    N  <- model$dim$N

    Z <- NULL
    for (j in 1:J) {
        Z <- cbind(Z, lagn(model$residuals, j))
    }
    Z <- na.omit(Z)
    K <- get.K.mat(model, J = J)
    KZ <- cbind(K, Z)

    vK <- K %x% diag(k)
    vKZ <- KZ %x% diag(k)

    U <- model$residuals[(J + 1):N, , drop = FALSE]

    if (ortogonalize) {
        vE0 <- vec(t(U)) -
            vK %*% ginv(t(vK) %*% vK) %*% t(vK) %*% vec(t(U))
        E0 <- t(matrix(vE0, ncol = N - J))
        SSR0 <- t(E0) %*% E0

        vEn <- vE0 -
            vKZ %*% ginv(t(vKZ) %*% vKZ) %*% t(vKZ) %*% vE0
        En <- t(matrix(vEn, ncol = N - J))
        SSRn <- t(En) %*% En
    } else {
        SSR0 <- t(U) %*% U

        vEn <- vec(t(U)) -
            vKZ %*% ginv(t(vKZ) %*% vKZ) %*% t(vKZ) %*% vec(t(U))
        En <- t(matrix(vEn, ncol = N - J))
        SSRn <- t(En) %*% En
    }

    LM.stat <- (N - J) * (k - matrix.trace(ginv(SSR0) %*% SSRn))

    result <- list(LM.stat = LM.stat,
                   crit.10 = qchisq(.90, df = J * k^2),
                   crit.5  = qchisq(.95, df = J * k^2),
                   crit.1  = qchisq(.99, df = J * k^2),
                   p.value = 1 - pchisq(LM.stat, df = J * k^2))

    return(result)
}
