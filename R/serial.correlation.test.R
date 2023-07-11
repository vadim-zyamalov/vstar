#' @importFrom MASS ginv
#' @importFrom matrixcalc matrix.trace
#' @importFrom matrixcalc vec
#'
#' @export
serial.correlation.test <- function(model, J = 1,
                                    stat.type = "all", ortogonalize = TRUE) {
    if (!"vstar" %in% class(model)) {
        stop("Wrong `model`: an object with estimated VSTAR model is needed!")
    }

    k  <- model$dim$k
    m  <- model$dim$m
    N  <- model$dim$N
    Nx <- ncol(model$data$X)

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

    result <- NULL

    if (stat.type %in% c("all", "LM")) {
        LM.stat <- N * (k - matrix.trace(ginv(SSR0) %*% SSRn))
        .df <- J * k^2

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

        .K <- m * k * Nx + 2 * (m - 1) + J * k^2
        .G <- J * k^2

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
        .s <- sqrt((ncol(Z)^2 * k^2 - 4) / (k^2 + ncol(Z)^2 - 5))
        .N <- N - ncol(K) - .5 * (k + ncol(Z) + 1)

        .df1 <- ncol(Z) * k
        .df2 <- .N * .s - .5 * ncol(Z) * k + 1

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
        l.stat <- (.5 * (k + ncol(Z) + 1) + ncol(K) - N) *
            log(det(SSRn) / det(SSR0))
        .df <- ncol(Z) * k

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
