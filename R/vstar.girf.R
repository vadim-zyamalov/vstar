#' @title
#' Generalized Impulse Response Functions
#'
#' @description
#' A function for constructing generalized impulse response functions using
#' Koop, et al. (1996) approach.
#'
#' @details
#' Exogenous variables are bootstrapped. \eqn{n} values of every exogenous
#' variable is sampled from \eqn{history} values. Determined components
#' are restored to continue the values from \eqn{history}.
#'
#' New values of a transition variable depend on its nature.
#' - If it's a lag of exogenous variable (maybe zero lag) then new values are
#' taken from old and new values of the corresponding exogenous variable.
#' - If it's a lag of endogenous variable then new values are taken from
#' the existing ones, and after they are over the calculated values of the
#' corresponding variable are used.
#' - If it's a custom variable then it's bootstrapped.
#'
#' The error term is bootstrapped using the \eqn{history} values or the model
#' residuals.
#'
#' @param model an object of S3-classes `vstar`.
#' @param shack a number of an equation, in which a shock appears.
#' @param n a number of periods to calculate IRF.
#' @param history a range of observations, which forms the history before the
#' impulse range.
#' @param iter a number of steps to get the distribution of GIRF.
#' @param R a number of repetitions in each step.
#' @param cores a number of cores to be used.
#' If greater than 1 then a grid search will be done in parallel.
#' The number of cores has the upper bound of processor cores number minus 1.
#'
#' @return
#' An object of S3 class `vstar.irf` containing
#' - matrices of means of responses to positive and negative shocks.
#' - matrices of corresponding standard deviations.
#' - the number of the shock.
#' - the length of the response range.
#'
#' @references
#' G. Koop, M. H. Pesaran, and S. M. Potter,
#' “Impulse response analysis in nonlinear multivariate models,”
#' Journal of Econometrics, vol. 74, no. 1, pp. 119–147, 1996.
#'
#' @export
vstar.girf <- function(model,
                       shock = 1,
                       n = 10,
                       history = 1:model$dim$N,
                       iter = 1000,
                       R = 100,
                       cores = 1) {
    k  <- model$dim$k
    p  <- model$dim$p.fixed
    m  <- model$dim$m
    N  <- max(history)
    Nx <- ncol(model$data$X)
    Nd <- 0

    if (length(model$params$trans) == 1 &&
        typeof(model$params$trans) == "character") {
        trans <- model$params$trans
        trans.lag <- 0
    } else if (typeof(model$params$trans) == "list") {
        trans <- model$params$trans[[1]]
        trans.lag <- model$params$trans[[2]]
    } else {
        trans <- NULL
        trans.lag <- NULL
    }

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

    cluster <- makeCluster(min(cores, detectCores() - 1), type = "SOCK")
    registerDoSNOW(cluster)

    results <- foreach(i = 1:iter,
                       .packages = c("vstar"),
                       .options.snow = list(progress = progress)) %dopar% {
        iter.Ksi <- matrix(
            boot.Ksi[sample(prod(dim(boot.Ksi)), size = k, replace = TRUE)],
            ncol = k)
        iter.E <- t(P %*% t(iter.Ksi))
        iter.shock <- iter.E[shock]

        iter.z <- NULL
        if (Nc > 0) for (j in 1:Nc) {
            iter.z <- cbind(
                iter.z,
                get.bootstrap(model$data$X[history, names.z], n))
        }

        iter.exo <- NULL
        if (!is.null(model$params$exog)) for (j in model$params$exog) {
            iter.exo <- cbind(
                iter.exo,
                get.bootstrap(model$data$X[history, j], n))
            colnames(iter.exo) <- model$params$exog
        }

        endo.trans <- FALSE
        diff.trans <- FALSE
        if (!is.null(trans) &&
            !is.null(model$params$exog) &&
            trans %in% model$params$exog) {
                if (trans.lag == 0) {
                    iter.trans <- iter.exo[, trans]
                } else {
                    range1 <- (N - trans.lag):N
                    range2 <- 1:(n - trans.lag)
                    iter.trans   <- model$data$X[history, trans][range1]
                    iter.trans   <- c(iter.trans, iter.exo[, trans][range2])
                    iter.trans.p <- iter.trans
                    iter.trans.m <- iter.trans
                }
        } else if (!is.null(trans) &&
                   trans %in% model$params$endo) {
            endo.trans <- TRUE
            iter.trans   <- model$data$S[N]
            iter.trans.p <- iter.trans
            iter.trans.m <- iter.trans

            if (!trans %in% colnames(model$data$Y)) {
                diff.trans <- TRUE
            }
        } else {
            iter.trans   <- c(model$data$S[N],
                              get.bootstrap(model$data$S[history], n))
            iter.trans.p <- iter.trans
            iter.trans.m <- iter.trans
        }

        iter.irf.p <- NULL
        iter.irf.m <- NULL

        for (r in 1:R) {
            tmp.eps <- matrix(boot.Ksi[sample(prod(dim(boot.Ksi)),
                                              size = k * (n + 1),
                                              replace = TRUE)],
                              ncol = k)
            tmp.eps <- t(P %*% t(tmp.eps))
            #tmp.eps <- matrix(rnorm(k * (n + 1)), ncol = k)
            tmp.eps.p <- tmp.eps
            tmp.eps.p[1, shock] <- abs(iter.shock)
            tmp.eps.m <- tmp.eps
            tmp.eps.m[1, shock] <- -abs(iter.shock)

            step.yr.p <- NULL
            step.yr.m <- NULL
            step.yr   <- NULL

            for (j in 1:(n + 1)) {
                range <- (k * (p - 1) + 1):(k * p)
                if (j == 1) {
                    step.xr   <- model$data$X[N - 1, (1 + Nc):(Nx - Nd)]
                    step.xr   <- c(model$data$Y[N - 1, ], step.xr[-range])
                    step.xr.p <- step.xr
                    step.xr.m <- step.xr
                } else {
                    step.xr.p <- c(step.yr.p[j - 1, ], step.xr.p[-range])
                    step.xr.m <- c(step.yr.m[j - 1, ], step.xr.m[-range])
                    step.xr   <- c(step.yr[j - 1, ],   step.xr[-range])
                }

                boot.G   <- G.func(t(unity(m - 1)) %x% iter.trans[j],
                                   model$g,
                                   model$c)
                boot.G.p <- G.func(t(unity(m - 1)) %x% iter.trans.p[j],
                                   model$g,
                                   model$c)
                boot.G.m <- G.func(t(unity(m - 1)) %x% iter.trans.m[j],
                                   model$g,
                                   model$c)

                step.xr.p <- matrix(
                    c(if (!is.null(iter.z)) iter.z[j, ] else NULL,
                      step.xr.p,
                      if (!is.null(boot.c)) 1 else NULL,
                      if (!is.null(boot.t)) boot.t[j] else NULL,
                      if (!is.null(boot.s)) boot.s[j, ] else NULL,
                      if (!is.null(iter.exo)) iter.exo[j, ] else NULL),
                    ncol = Nx)

                step.xr.m <- matrix(
                    c(if (!is.null(iter.z)) iter.z[j, ] else NULL,
                      step.xr.m,
                      if (!is.null(boot.c)) 1 else NULL,
                      if (!is.null(boot.t)) boot.t[j] else NULL,
                      if (!is.null(boot.s)) boot.s[j, ] else NULL,
                      if (!is.null(iter.exo)) iter.exo[j, ] else NULL),
                    ncol = Nx)

                step.xr <- matrix(
                    c(if (!is.null(iter.z)) iter.z[j, ] else NULL,
                      step.xr,
                      if (!is.null(boot.c)) 1 else NULL,
                      if (!is.null(boot.t)) boot.t[j] else NULL,
                      if (!is.null(boot.s)) boot.s[j, ] else NULL,
                      if (!is.null(iter.exo)) iter.exo[j, ] else NULL),
                    ncol = Nx)

                PSI   <- t(t(c(1, boot.G)))   %x% diag(k)
                PSI.p <- t(t(c(1, boot.G.p))) %x% diag(k)
                PSI.m <- t(t(c(1, boot.G.m))) %x% diag(k)

                tmp.yr.p  <- t(PSI.p) %*% t(model$coef) %*% t(step.xr.p) +
                    tmp.eps.p[j, ]
                step.yr.p <- rbind(step.yr.p, t(tmp.yr.p))

                tmp.yr.m  <- t(PSI.m) %*% t(model$coef) %*% t(step.xr.m) +
                    tmp.eps.m[j, ]
                step.yr.m <- rbind(step.yr.m, t(tmp.yr.m))

                tmp.yr  <- t(PSI) %*% t(model$coef) %*% t(step.xr) +
                    tmp.eps[j, ]
                step.yr <- rbind(step.yr, t(tmp.yr))

                if (endo.trans) {
                    if (trans.lag == 0) {
                        iter.trans.p <- c(
                            iter.trans.p, tmp.yr.p[model$params$endo == trans]
                        )
                        iter.trans.m <- c(
                            iter.trans.m, tmp.yr.m[model$params$endo == trans]
                        )
                        iter.trans <- c(
                            iter.trans, tmp.yr[model$params$endo == trans]
                        )
                    } else {
                        if (j - trans.lag <= 0) {
                            iter.trans.p <- c(
                                iter.trans.p,
                                model$data$Y[history, model$params$endo == trans
                                             ][N - trans.lag + j])
                            iter.trans.m <- c(
                                iter.trans.m,
                                model$data$Y[history, model$params$endo == trans
                                             ][N - trans.lag + j])
                            iter.trans <- c(
                                iter.trans,
                                model$data$Y[history, model$params$endo == trans
                                             ][N - trans.lag + j])
                        } else {
                            iter.trans.p <- c(
                                iter.trans.p,
                                step.yr.p[
                                j - trans.lag, model$params$endo == trans]
                            )
                            iter.trans.m <- c(
                                iter.trans.m,
                                step.yr.m[
                                j - trans.lag, model$params$endo == trans]
                            )
                            iter.trans <- c(
                                iter.trans,
                                step.yr[
                                j - trans.lag, model$params$endo == trans]
                            )
                        }
                    }

                    if (diff.trans) {
                        iter.trans.p[j + 1] <-
                            iter.trans.p[j + 1] + iter.trans.p[j]
                        iter.trans.m[j + 1] <-
                            iter.trans.m[j + 1] + iter.trans.m[j]
                        iter.trans[j + 1] <-
                            iter.trans[j + 1] + iter.trans[j]
                    }
                }
            }

            iter.irf.p <- cbind(iter.irf.p, vec(step.yr.p - step.yr))
            iter.irf.m <- cbind(iter.irf.m, vec(step.yr.m - step.yr))
        }

        list(p = rowMeans(iter.irf.p),
             var.p = apply(iter.irf.p, 1, var) / R,
             m = rowMeans(iter.irf.m),
             var.m = apply(iter.irf.m, 1, var) / R)
    }

    stopCluster(cluster)

    irf.p     <- NULL
    irf.m     <- NULL
    irf.var.p <- NULL
    irf.var.m <- NULL

    for (ii in 1:iter) {
        irf.p     <- cbind(irf.p, results[[ii]]$p)
        irf.var.p <- cbind(irf.var.p, results[[ii]]$var.p)
        irf.m     <- cbind(irf.m, results[[ii]]$m)
        irf.var.m <- cbind(irf.var.m, results[[ii]]$var.m)
    }

    res.irf.p <- matrix(rowMeans(irf.p), ncol = k)
    colnames(res.irf.p) <- colnames(model$data$Y)

    res.irf.m <- matrix(rowMeans(irf.m), ncol = k)
    colnames(res.irf.m) <- colnames(model$data$Y)

    res.irf.var.p <- matrix(rowMeans(irf.var.p), ncol = k)
    colnames(res.irf.var.p) <- colnames(model$data$Y)

    res.irf.var.m <- matrix(rowMeans(irf.var.m), ncol = k)
    colnames(res.irf.var.p) <- colnames(model$data$Y)

    result <- list(pos = res.irf.p,
                   pos.sd = sqrt(res.irf.var.p),
                   neg = res.irf.m,
                   neg.sd = sqrt(res.irf.var.m),
                   shock = colnames(model$data$Y)[shock],
                   n = n)

    class(result) <- "vstar.irf"

    return(result)
}
