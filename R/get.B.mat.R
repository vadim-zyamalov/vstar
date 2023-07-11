#' @importFrom MASS ginv
#' @importFrom matrixcalc vec
get.B.mat <- function(dataset, m, g, thr, g.function = "L") {
    if (is.null(g.function) && !("g.function" %in% names(dataset))) {
        stop("B: No transition function is provided!")
    }

    k  <- dataset$dim$k
    N  <- dataset$dim$N
    Nx <- ncol(dataset$data$X)

    G.func <- get.G.function(g.function)

    sm <- dataset$data$S %x% t(unity(m - 1))
    gm <- unity(N) %x% t(g)
    cm <- unity(N) %x% t(thr)

    G.mat <- G.func(sm, gm, cm)

    M <- (cbind(unity(N), G.mat) %x% diag(k) %x% t(unity(Nx))) *
            (t(unity(k * m)) %x% dataset$data$X %x% unity(k))

    vec.B <- ginv(t(M) %*% M) %*% t(M) %*% vec(t(dataset$data$Y))

    return(list(vec.B = vec.B,
                M = M))
}
