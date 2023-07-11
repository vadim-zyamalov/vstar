#' @importFrom MASS ginv
#' @importFrom matrixcalc matrix.trace
#' @importFrom matrixcalc vec
#'
#' @export
linearity.test <- function(dataset, J = 1, stat.type = "all") {
    if (!any(c("vstar", "vstar.data") %in% class(dataset))) {
        stop("Wrong `model`: an object with prepared data or estimated model is needed!")
    }

    k  <- dataset$dim$k
    N  <- dataset$dim$N
    Nx <- ncol(dataset$data$X)

    Zn <- NULL
    for (j in 1:J) {
        Zn <- cbind(Zn,
                    (dataset$data$S^j %x% t(unity(Nx))) * dataset$data$X)
    }

    XZ <- cbind(dataset$data$X, Zn)

    vX  <- dataset$data$X %x% diag(k)
    vXZ <- XZ %x% diag(k)

    vE0 <- vec(t(dataset$data$Y)) -
        vX %*% ginv(t(vX) %*% vX) %*% t(vX) %*% vec(t(dataset$data$Y))
    E0 <- t(matrix(vE0, ncol = N))
    SSR0 <- t(E0) %*% E0

    vEn <- vE0 -
        vXZ %*% ginv(t(vXZ) %*% vXZ) %*% t(vXZ) %*% vE0
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
