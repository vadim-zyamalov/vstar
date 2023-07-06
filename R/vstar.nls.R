#' @importFrom stats nls
#' @importFrom stats nls.control
#' @importFrom matrixcalc vec
#'
#' @export
vstar.nls <- function(model,
                      tol = 1e-6,
                      iter = 10000,
                      control = nls.control()) {
    vec.y <- vec(t(model$data$Y))
    G.func <- get.G.function(model$g.function)

    yf <- function(g, thr, vec.B) { # nolint
        BM <- get.B.mat(model, g, thr, G.func)
        BM$M %*% BM$vec.B
    }

    iter.g <- model$estimates$g
    iter.c <- model$estimates$c
    iter.B <- model$estimates$vec.B # nolint

    converged <- FALSE

    for (step in 1:iter) {
        resnls <- nls(vec.y ~ yf(gx, cx, iter.B),
                      start = list(gx = iter.g, cx = iter.c),
                      control = control)

        result <- coef(resnls)
        dd <- result - c(iter.g, iter.c)

        iter.g <- result[startsWith(names(result), "gx")]
        iter.c <- result[startsWith(names(result), "cx")]
        iter.BM <- get.B.mat(model, iter.g, iter.c, G.func)
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

    nls.result <- list(vec.B = iter.BM$vec.B,
                       vec.E = vec.y - iter.BM$M %*% iter.BM$vec.B,
                       M = iter.BM$M,
                       SSR = drop(t(model$estimates$vec.E) %*%
                                model$estimates$vec.E),
                       g = iter.g,
                       c = iter.c)

    final.est <- get.estimates(nls.result, model, model$g.function)

    result <- list(coef = final.est$coef,
                   sd = final.est$sd,
                   t.stat = final.est$t.stat,
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
