#' @importFrom matrixcalc vec
get.K.mat <- function(model, J = 0) {
    if (!"coef" %in% names(model)) {
        stop("K: Please, estimate the model first!")
    }

    k  <- model$dim$k
    m  <- model$dim$m
    N  <- model$dim$N

    H.mat <- get.H.mat(model$coef)

    G.func <- get.G.function(model$g.function)
    G.deriv <- get.G.derivative(model$g.function)

    gg <- unity(N) %x% t(model$g)
    cc <- unity(N) %x% t(model$c)

    G.mat <- G.func(model$data$S, gg, cc)
    G.g.mat <- G.deriv[[1]](model$data$S, gg, cc)
    G.c.mat <- G.deriv[[2]](model$data$S, gg, cc)

    K.mat <- NULL

    for (i in (J + 1):N) {
        K.t <- NULL
        PSI <- t(t(c(1, G.mat[i, ]))) %x% diag(k)
        x <- t(model$data$X[i, , drop = FALSE])

        for (j in seq_len(length(model$coef))) {
            K.t <- cbind(K.t,
                         t(PSI) %*% t(H.mat[[j]]) %*% x)
        }

        for (j in 1:(m - 1)) {
            elems <- rep(0, m - 1)
            elems[j] <- G.g.mat[i, j]
            G.g <- t(c(0, elems)) %x% diag(k)

            elems[j] <- G.c.mat[i, j]
            G.c <- t(c(0, elems)) %x% diag(k)

            K.t <- cbind(K.t,
                         G.g %*% t(model$coef) %*% x,
                         G.c %*% t(model$coef) %*% x)
        }

        K.mat <- rbind(K.mat, t(vec(K.t)))
    }

    return(K.mat)
}
