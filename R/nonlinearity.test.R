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

    result <- NULL

    if (stat.type %in% c("all", "LM")) {
        LM.stat <- N * (k - matrix.trace(ginv(SSR0) %*% SSRn))
        .df <- k * J * Nx

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

        .K <- k * (m + J) * Nx
        .G <- k * J * Nx

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
        .N <- N - ncol(K) - .5 * (k + ncol(Zn) + 1)

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
        l.stat <- (.5 * (k + ncol(Zn) + 1) + ncol(K) - N) *
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
