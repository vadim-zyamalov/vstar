#' @importFrom MASS ginv
#' @importFrom matrixcalc vec
get.estimates <- function(est, model, g.function) {
    k <- model$dim$k
    m <- model$dim$m
    N <- model$dim$N

    G.func <- get.G.function(g.function)

    vec.Yf <- est$M %*% est$vec.B
    vec.E  <- vec(t(model$data$Y)) - vec.Yf

    B  <- matrix(est$vec.B, ncol = k * m)
    colnames(B) <- paste("R",
                         rep(1:m, each = k),
                         ".",
                         colnames(model$data$Y), sep = "")
    rownames(B) <- colnames(model$data$X)

    Yf <- t(matrix(vec.Yf, ncol = N, byrow = FALSE))
    colnames(Yf) <- paste("fit.", colnames(model$data$Y), sep = "")
    rownames(Yf) <- rownames(model$data$Y)

    E <- t(matrix(vec.E, ncol = N, byrow = FALSE))
    colnames(E) <- paste("resid.", colnames(model$data$Y), sep = "")
    rownames(E) <- rownames(model$data$Y)

    OMEGA <- (1 / N) * t(E) %*% E

    gg <- unity(N) %x% t(est$g)
    cc <- unity(N) %x% t(est$c)

    G.mat <- G.func(model$data$S, gg, cc)

    GX <- cbind(model$data$X,
                (t(unity(ncol(model$data$X))) %x% G.mat) * model$data$X)
    VAR <- matrix(diag(OMEGA %x% MASS::ginv(t(GX) %*% GX)), ncol = k * m)
    colnames(VAR) <- paste("R",
                           rep(1:m, each = k),
                           ".",
                           colnames(model$data$Y), sep = "")
    rownames(VAR) <- colnames(model$data$X)

    return(list(coef = B,
                sd = sqrt(VAR),
                t.stat = B / sqrt(VAR),
                fitted.values = Yf,
                residuals = E,
                cov = OMEGA,
                g.function = g.function))
}
