#' @export
vstar.irf <- function(model,
                      shock = 1,
                      n = 10,
                      history = 1:model$dim$N,
                      iter = 1000,
                      R = 100) {
    k  <- model$dim$k
    m  <- model$dim$m
    N  <- max(history)
    Nx <- ncol(model$data$X)
    Nd <- 0

    Ns <- 0
    boot.s <- NULL
    if (!is.null(model$params$season)) {
        Ns <- model$params$season
        names.s <- paste("(s", 2:Ns, ")", sep = "")
        s.chunk <- if (all(model$data$X[N, names.s] == 0)) {
            c(1, model$data$X[N, names.s])
        } else {
            c(0, model$data$X[N, names.s])
        }
        boot.s <- s.chunk
        for (j in 1:n) {
            s.chunk <- c(s.chunk[Ns], s.chunk[1:(Ns - 1)])
            boot.s <- rbind(boot.s, s.chunk)
        }
        boot.s <- boot.s[, -1]
        Nd <- Nd + Ns - 1
    }

    boot.c <- NULL
    if (model$params$const) {
        boot.c <- 1
        Nd <- Nd + 1
    }

    boot.t <- NULL
    if (model$params$trend) {
        boot.t <- N:(N + n)
        Nd <- Nd + 1
    }

    Nc <- 0
    if (!is.null(model$params$coint.beta)) {
        Nc <- nrow(model$params$coint.beta)
        names.z <- paste("(z", j, ")", sep = "")
    }

    G.func <- model$func$g.function

    boot.E <- model$residuals[history, , drop = FALSE]
    OMEGA <- (1 / length(history)) * t(boot.E) %*% boot.E
    P <- t(chol(OMEGA))
    boot.Ksi <- t(MASS::ginv(P) %*% t(boot.E))

    pb <- progress_bar$new(
        format = ":percent [:bar] :elapsed | ETA: :eta",
        total = iter,
        width = 60
    )
    progress <- function(n) pb$tick()

    results <- foreach(i = 1:iter) %do% {
        progress(i)

        for (r in 1:R) {
            boot.z <- NULL
            if (!is.null(model$params$coint.beta)) {
                for (j in 1:Nc) {
                    boot.z <- cbind(
                        boot.z,
                        get.bootstrap(model$data$X[history, names.z], n))
                }
            }

            boot.exo <- NULL
            if (!is.null(model$params$exog)) for (j in model$params$exog) {
                boot.exo <- cbind(
                    boot.exo,
                    get.bootstrap(model$data$X[history, j], n))
            }

            boot.trans <- c(
                model$data$S[N],
                get.bootstrap(model$data$S[history], n))
            boot.G <- G.func(t(vstar:::unity(m - 1)) %x% boot.trans,
                             model$g,
                             model$c)

            tmp.Ksi <- matrix(
                boot.Ksi[sample(prod(dim(boot.Ksi)), size = k, replace = TRUE)],
                ncol = k)
            tmp.E <- t(P %*% t(tmp.Ksi))
            tmp.shock <- tmp.E[shock]

            tmp.eps <- matrix(rnorm(k * (n + 1)), ncol = k)
            tmp.eps.p <- tmp.eps
            tmp.eps.p[1, shock] <- abs(tmp.shock)
            tmp.eps.m <- tmp.eps
            tmp.eps.m[1, shock] <- -abs(tmp.shock)

            boot.yr.p <- NULL
            boot.yr.m <- NULL
            boot.y    <- NULL
            for (j in 1:(n + 1)) {
                if (j == 1) {
                    boot.x <- model$data$X[N - 1, (1 + Nc):(Nx - Nd)]
                    boot.x <- c(model$data$Y[N - 1, ], boot.x[-(1:k)])
                    boot.xr.p <- boot.x
                    boot.xr.m <- boot.x
                } else {
                    boot.xr.p <- c(boot.yr.p[j - 1, ], boot.xr.p[-(1:k)])
                    boot.xr.m <- c(boot.yr.m[j - 1, ], boot.xr.m[-(1:k)])
                    boot.x    <- c(boot.y[j - 1, ], boot.x[-(1:k)])
                }

                tmp.xr.p <- matrix(
                    c(if (!is.null(boot.z)) boot.z[j, ] else NULL,
                      boot.xr.p,
                      if (!is.null(boot.c)) 1 else NULL,
                      if (!is.null(boot.t)) boot.t[j] else NULL,
                      if (!is.null(boot.s)) boot.s[j, ] else NULL),
                    ncol = Nx)

                tmp.xr.m <- matrix(
                    c(if (!is.null(boot.z)) boot.z[j, ] else NULL,
                      boot.xr.m,
                      if (!is.null(boot.c)) 1 else NULL,
                      if (!is.null(boot.t)) boot.t[j] else NULL,
                      if (!is.null(boot.s)) boot.s[j, ] else NULL),
                    ncol = Nx)

                tmp.x <- matrix(
                    c(if (!is.null(boot.z)) boot.z[j, ] else NULL,
                      boot.x,
                      if (!is.null(boot.c)) 1 else NULL,
                      if (!is.null(boot.t)) boot.t[j] else NULL,
                      if (!is.null(boot.s)) boot.s[j, ] else NULL),
                    ncol = Nx)

                PSI <- t(t(c(1, boot.G[j, ]))) %x% diag(k)

                tmp.yr.p  <- t(PSI) %*% t(model$coef) %*% t(tmp.xr.p) + 
                    tmp.eps.p[j, ]
                boot.yr.p <- rbind(boot.yr.p, t(tmp.yr.p))

                tmp.yr.m  <- t(PSI) %*% t(model$coef) %*% t(tmp.xr.m) + 
                    tmp.eps.m[j, ]
                boot.yr.m <- rbind(boot.yr.m, t(tmp.yr.m))

                tmp.y  <- t(PSI) %*% t(model$coef) %*% t(tmp.x) + 
                    tmp.eps[j, ]
                boot.y <- rbind(boot.y, t(tmp.y))
            }

            if (r == 1) {
                irf.p <- boot.yr.p - boot.y
                irf.m <- boot.yr.m - boot.y
            } else {
                irf.p <- irf.p + (boot.yr.p - boot.y)
                irf.m <- irf.m + (boot.yr.m - boot.y)
            }
        }

        irf.p <- irf.p / R
        irf.m <- irf.m / R

        list(p = irf.p, m = irf.m)
    }

    irf.p <- array(NA, dim = c(n + 1, k, iter))
    irf.m <- array(NA, dim = c(n + 1, k, iter))

    for (ii in 1:iter) {
        irf.p[, , ii] <- results[[ii]]$p
        irf.m[, , ii] <- results[[ii]]$m
    }

    irf.mean.p <- matrix(NA, ncol = k, nrow = n + 1)
    colnames(irf.mean.p) <- colnames(model$data$Y)

    irf.sd.p   <- matrix(NA, ncol = k, nrow = n + 1)
    colnames(irf.sd.p) <- colnames(model$data$Y)

    irf.mean.m <- matrix(NA, ncol = k, nrow = n + 1)
    colnames(irf.mean.m) <- colnames(model$data$Y)

    irf.sd.m   <- matrix(NA, ncol = k, nrow = n + 1)
    colnames(irf.sd.m) <- colnames(model$data$Y)

    for (i in 1:(n + 1)) for (j in 1:k) {
        irf.mean.p[i, j] <- mean(irf.p[i, j, ])
        irf.mean.m[i, j] <- mean(irf.m[i, j, ])

        irf.sd.p[i, j] <- sqrt(var(irf.p[i, j, ]))
        irf.sd.m[i, j] <- sqrt(var(irf.m[i, j, ]))
    }

    result <- list(pos = irf.mean.p,
                   pos.sd = irf.sd.p,
                   neg = irf.mean.m,
                   neg.sd = irf.sd.m,
                   shock = colnames(model$data$Y)[shock],
                   n = n)

    class(result) <- "vstar.irf"

    return(result)
}

#' @export
plot.vstar.irf <- function(res, ...) {
    k <- ncol(res$pos)
    n <- res$n

    names.y <- colnames(res$pos)

    graphs <- list()

    for (i in 1:k) {
        plot.data <- cbind(res$pos[, i], res$neg.sd[, i],
                           res$neg[, i], res$pos.sd[, i])
        colnames(plot.data) <- c("pos", "pos.sd", "neg", "neg.sd")
        rownames(plot.data) <- 0:10

        graphs[[paste("g", i, sep = "")]] <-
            ggplot(data = as.data.frame(plot.data),
                            mapping = aes(x = 0:n)) +
            geom_line(aes(y = pos, colour = "Positive")) +
            geom_line(aes(y = pos + 2 * pos.sd),
                               linetype = "dashed", colour = "red") +
            geom_line(aes(y = pos - 2 * pos.sd),
                               linetype = "dashed", colour = "red") +
            geom_line(aes(y = neg, colour = "Negative")) +
            geom_line(aes(y = neg + 2 * neg.sd),
                               linetype = "dashed", colour = "blue") +
            geom_line(aes(y = neg - 2 * neg.sd),
                               linetype = "dashed", colour = "blue") +
            scale_colour_manual(name = "Shock",
                                         values = c("blue", "red")) +
            scale_x_continuous(breaks = 0:n) +
            labs(title = paste("Response of", names.y[i],
                                        "to", res$shock),
                          x = "Period",
                          y = paste("Response to", res$shock))
    }

    plot_grid(plotlist = graphs, align = "v", ncol = 1)
    #return(graphs)
}

get.bootstrap <- function(x, n) {
    N <- length(x)
    idx <- sample(1:N, n, replace = TRUE)
    return(x[idx])
}
