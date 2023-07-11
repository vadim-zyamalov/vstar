#' @importFrom MASS ginv
#' @importFrom matrixcalc vec
get.estimates <- function(est, dataset, g.function) {
    if (!any(c("vstar", "vstar.data") %in% class(dataset))) {
        stop("Wrong `model`: an object with prepared data or estimated model is needed!")
    }

    k <- dataset$dim$k
    m <- est$m
    N <- dataset$dim$N
    Nx <- ncol(dataset$data$X)

    G.func <- get.G.function(g.function)

    vec.Yf <- est$M %*% est$vec.B
    vec.E  <- vec(t(dataset$data$Y)) - vec.Yf

    B  <- matrix(est$vec.B, ncol = k * m)
    colnames(B) <- paste("R",
                         rep(1:m, each = k),
                         ".",
                         colnames(dataset$data$Y), sep = "")
    rownames(B) <- colnames(dataset$data$X)

    Yf <- t(matrix(vec.Yf, ncol = N, byrow = FALSE))
    colnames(Yf) <- paste("fit.", colnames(dataset$data$Y), sep = "")
    rownames(Yf) <- rownames(dataset$data$Y)

    E <- t(matrix(vec.E, ncol = N, byrow = FALSE))
    colnames(E) <- paste("resid.", colnames(dataset$data$Y), sep = "")
    rownames(E) <- rownames(dataset$data$Y)

    OMEGA <- (1 / N) * t(E) %*% E

    sm <- dataset$data$S %x% t(unity(m - 1))
    gm <- unity(N) %x% t(est$g)
    cm <- unity(N) %x% t(est$c)

    G.mat <- G.func(sm, gm, cm)
    G.mat <- cbind(1, G.mat)

    GX <- (t(unity(Nx)) %x% G.mat) * (dataset$data$X %x% t(unity(m)))
    VAR <- matrix(diag(ginv(t(GX) %*% GX) %x% OMEGA),
                  ncol = k * m,
                  byrow = TRUE)
    colnames(VAR) <- paste("R",
                           rep(1:m, each = k),
                           ".",
                           colnames(dataset$data$Y), sep = "")
    rownames(VAR) <- colnames(dataset$data$X)

    return(list(coef = B,
                sd = sqrt(VAR),
                t.stat = B / sqrt(VAR),
                fitted.values = Yf,
                residuals = E,
                cov = OMEGA,
                g.function = g.function))
}
