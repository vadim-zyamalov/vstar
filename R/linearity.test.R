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

    result <- NULL

    if (stat.type %in% c("all", "LM")) {
        LM.stat <- N * (k - matrix.trace(ginv(SSR0) %*% SSRn))
        .df <- J * k * Nx

        result <- rbind(
            result,
            "LM" = c(stat    = LM.stat,
                     crit.10 = qchisq(.90, df = .df),
                     crit.5  = qchisq(.95, df = .df),
                     crit.1  = qchisq(.99, df = .df),
                     p.value = 1 - pchisq(LM.stat, df = .df))
        )
    }

    if (stat.type %in% c("all", "r-LM")) {
        LM.stat <- N * (k - matrix.trace(ginv(SSR0) %*% SSRn))

        .K <- (J + 1) * k * Nx
        .G <- J * k * Nx

        F.stat <- LM.stat * (k * N - .K) / (.G * k * N)

        .df1 <- .G
        .df2 <- k * N - .K

        result <- rbind(
            result,
            "r-LM" = c(stat    = F.stat,
                       crit.10 = qf(.90, df1 = .df1, df2 = .df2),
                       crit.5  = qf(.95, df1 = .df1, df2 = .df2),
                       crit.1  = qf(.99, df1 = .df1, df2 = .df2),
                       p.value = 1 - pf(F.stat, df1 = .df1, df2 = .df2))
        )
    }

    if (stat.type %in% c("all", "Rao")) {
        .s <- sqrt((ncol(Zn)^2 * k^2 - 4) / (k^2 + ncol(Zn)^2 - 5))
        .N <- N - ncol(dataset$data$X) - .5 * (k + ncol(Zn) + 1)

        .df1 <- ncol(Zn) * k
        .df2 <- .N * .s - .5 * ncol(Zn) * k + 1

        F.stat <- ((det(SSR0) / det(SSRn))^(1/.s) - 1) * (.df2 / .df1)

        result <- rbind(
            result,
            "Rao" = c(stat    = F.stat,
                      crit.10 = qf(.90, df1 = .df1, df2 = .df2),
                      crit.5  = qf(.95, df1 = .df1, df2 = .df2),
                      crit.1  = qf(.99, df1 = .df1, df2 = .df2),
                      p.value = 1 - pf(F.stat, df1 = .df1, df2 = .df2))
        )
    }

    if (stat.type %in% c("all", "Bartlett-Wilks")) {
        l.stat <- (.5 * (k + ncol(Zn) + 1) + ncol(dataset$data$X) - N) *
            log(det(SSRn) / det(SSR0))
        .df <- ncol(Zn) * k

        result <- rbind(
            result,
            "Bartlett-Wilks" = c(stat    = l.stat,
                                 crit.10 = qchisq(.90, df = .df),
                                 crit.5  = qchisq(.95, df = .df),
                                 crit.1  = qchisq(.99, df = .df),
                                 p.value = 1 - pchisq(l.stat, df = .df))
        )
    }

    return(result)
}
