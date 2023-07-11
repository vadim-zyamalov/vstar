#' @importFrom MASS ginv
#' @importFrom matrixcalc matrix.trace
#' @importFrom matrixcalc vec
#'
#' @export
nonlinearity.test <- function(model, J = 1, stat.type = "all") {
    if (!"vstar" %in% class(model)) {
        stop("Wrong `model`: an object with estimated VSTAR model is needed!")
    }

    k  <- model$dim$k
    m  <- model$dim$m
    N  <- model$dim$N
    Nx <- ncol(model$data$X)

    Zn <- NULL
    for (j in 1:J) {
        Zn <- cbind(Zn,
                    (model$data$S^j %x% t(unity(Nx))) * model$data$X)
    }

    K <- get.K.mat(model, J = 0)
    KZ <- cbind(K, Zn)

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
                   crit.10 = qchisq(.90, df = k * N * Nx),
                   crit.5  = qchisq(.95, df = k * N * Nx),
                   crit.1  = qchisq(.99, df = k * N * Nx),
                   p.value = 1 - pchisq(LM.stat, df = k * N * Nx))

    return(result)
}
