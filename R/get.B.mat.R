#' @importFrom MASS ginv
#' @importFrom matrixcalc vec
get.B.mat <- function(model, g, thr, g.function = NULL) {
    if (is.null(g.function) && !("g.function" %in% names(model))) {
        stop("B: No transition function is provided!")
    }

    k  <- model$dim$k
    m  <- model$dim$m
    N  <- model$dim$N
    Nx <- ncol(model$data$X)

    G.func <- get.G.function(g.function)

    gg <- unity(N) %x% t(g)
    cc <- unity(N) %x% t(thr)

    G.mat <- G.func(model$data$S, gg, cc)

    M <- (cbind(unity(N), G.mat) %x% diag(k) %x% t(unity(Nx))) *
            (t(unity(k * m)) %x% model$data$X %x% unity(k))

    vec.B <- ginv(t(M) %*% M) %*% t(M) %*% vec(t(model$data$Y))

    return(list(vec.B = vec.B,
                M = M))
}
