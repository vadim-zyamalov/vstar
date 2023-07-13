#' @title
#' Get a set of estimated values
#'
#' @description
#' An auxiliary function returning a set of estimated values of a model.
#' Values are devectorized if needed.
#'
#' @param est a list containing the intermediate results of model estimation.
#' This list emerges as a result of inner loops in functions for estimating
#' VSTAR models (see \link{vstar.grid} and \link{vstar.nls}).
#' @param dataset an object of S3-classes `vstar.data` or `vstar`.
#' @param G.func a transition function to be used.
#' NB: it's a function, not its destription!
#'
#' @return
#' A list containing
#' * A matrix of estimated coefficients.
#' * A matrix of estimated coefficients standard devations.
#' * A matrix of estimated coefficients t-statistics.
#' * A matrix of fitted values.
#' * A matrix of residuals.
#' * An estimate of a covariance matrix.
#' * The value of `g.function`.
#'
#' @keywords internal
get.estimates <- function(est, dataset, G.func) {
    if (!any(c("vstar", "vstar.data") %in% class(dataset))) {
        stop("Wrong `model`: an object with prepared data or estimated model is needed!")
    }

    k <- dataset$dim$k
    m <- est$m
    N <- dataset$dim$N
    Nx <- ncol(dataset$data$X)

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
                cov = OMEGA))
}
