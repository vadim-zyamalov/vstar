#' @importFrom stats nls
#' @importFrom stats nls.control
#' @importFrom matrixcalc vec
#' @importFrom optimx optimx
#'
#' @export
vstar.nls <- function(model,
                      tol = 1e-6,
                      iter = 10000,
                      algorithm = "default",
                      control = nls.control()) {
    m <- model$dim$m
    
    vec.y <- vec(t(model$data$Y))
    G.func <- get.G.function(model$g.function)

    yf <- function(g, thr) { # nolint
        BM <- get.B.mat(model, g, thr, G.func)
        BM$M %*% BM$vec.B
    }

    Q <- function(gc) {
        g <- gc[1:(m - 1)]
        thr <- gc[m:length(gc)]
        BM <- get.B.mat(model, m, g, thr, G.func)
        vec.E <- vec(t(model$data$Y)) - BM$M %*% BM$vec.B
        drop(t(vec.E) %*% vec.E)
    }

    iter.g <- model$estimates$g
    iter.c <- model$estimates$c
    iter.B <- model$estimates$vec.B # nolint

    converged <- FALSE

    for (step in 1:iter) {
        cat(step, "\n")
        #iter.nls <- nls(vec.y ~ yf(gx, cx),
        #                start = list(gx = iter.g, cx = iter.c),
        #                algorithm = algorithm,
        #                control = control)

        #result <- coef(iter.nls)
        iter.nls <- optim(c(iter.g, iter.c), Q)
        result <- iter.nls$par

        dd <- result - c(iter.g, iter.c)

        #iter.g <- result[startsWith(names(result), "gx")]
        iter.g <- result[1:(m - 1)]
        #iter.c <- result[startsWith(names(result), "cx")]
        iter.c <- result[m:length(result)]

        iter.BM <- get.B.mat(model, m, iter.g, iter.c, G.func)
        iter.B <- iter.BM$vec.B

        if (sqrt(dd %*% dd) < tol) {
            converged <- TRUE
            cat("Tolerance level", tol, "achieved in", step, "steps.\n")
            break
        }
    }

    if (!converged) {
        cat("Reached limit of", iter, "steps.\n")
        cat("Possible absence of convergence!\n")
    }

    #nls.summary <- summary(iter.nls)
    #g.sd <- nls.summary$coefficients[, 2][names(iter.g)]
    #c.sd <- nls.summary$coefficients[, 2][names(iter.c)]

    nls.result <- list(vec.B = iter.BM$vec.B,
                       vec.E = vec.y - iter.BM$M %*% iter.BM$vec.B,
                       M = iter.BM$M,
                       SSR = drop(t(model$estimates$vec.E) %*%
                                model$estimates$vec.E),
                       m = m,
                       g = iter.g,
                       #g.sd = g.sd,
                       c = iter.c)
                       #c.sd = c.sd)

    final.est <- get.estimates(nls.result, model, model$g.function)

    result <- list(coef = final.est$coef,
                   sd = final.est$sd,
                   t.stat = final.est$t.stat,
                   g = iter.g,
                   #g.sd = g.sd,
                   c = iter.c,
                   #c.sd = c.sd,
                   fitted.values = final.est$fitted.values,
                   residuals = final.est$residuals,
                   cov = final.est$cov,
                   g.function = final.est$g.function,
                   estimates = nls.result,
                   params = model$params,
                   dim = model$dim,
                   data = model$data)

    class(result) <- "vstar"

    return(result)
}
