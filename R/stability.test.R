#' @importFrom MASS ginv
#' @importFrom matrixcalc matrix.trace
#' @importFrom matrixcalc vec
#'
#' @export
stability.test <- function(model) {
    k  <- model$dim$k
    m <- model$dim$m
    N  <- model$dim$N
    Nx <- ncol(model$data$X)

    G.func <- get.G.function(model$g.function)
    gg <- unity(N) %x% t(model$g)
    cc <- unity(N) %x% t(model$c)
    G.mat <- G.func(model$data$S, gg, cc)

    Z <- NULL
    for (i in 1:N) {
        PSI <- t(t(c(1, G.mat[i, ]))) %x% diag(k)
        x <- t(model$data$X[i, , drop = FALSE])
        tau <- j / N

        Z <- rbind(Z,
                   t(vec(PSI %x% (tau * x))))
    }

    K <- get.K.mat(model, J = 0)
    KZ <- cbind(K, Z)

    vK <- K %x% diag(k)
    vKZ <- KZ %x% diag(k)

    vE <- vec(t(model$residuals))

    vE0 <- vE -
        vK %*% ginv(t(vK) %*% vK) %*% t(vK) %*% vE
    E0 <- t(matrix(vE0, ncol = N))
    SSR0 <- t(E0) %*% E0

    vEn <- vE0 -
        vKZ %*% ginv(t(vKZ) %*% vKZ) %*% t(vKZ) %*% vE0
    En <- t(matrix(vEn, ncol = N))
    SSRn <- t(En) %*% En

    LM.stat <- N * (k - matrix.trace(ginv(SSR0) %*% SSRn))

    result <- list(LM.stat = LM.stat,
                   crit.10 = qchisq(.90, df = m * k^2 * Nx),
                   crit.5  = qchisq(.95, df = m * k^2 * Nx),
                   crit.1  = qchisq(.99, df = m * k^2 * Nx),
                   p.value = 1 - pchisq(LM.stat, df = m * k^2 * Nx))

    return(result)
}
