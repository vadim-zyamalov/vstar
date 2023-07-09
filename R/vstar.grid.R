#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom arrangements icombinations
#' @importFrom arrangements ncombinations
#' @importFrom arrangements ipermutations
#' @importFrom arrangements npermutations
#' @importFrom MASS ginv
#' @importFrom matrixcalc vec
#'
#' @export
vstar.grid <- function(model,
                       g.function = "L",
                       gamma.limits = c(.1, 100),
                       points = 200,
                       trim = .15,
                       gap = .1) {
    m <- model$dim$m
    G.func <- get.G.function(g.function)

    if (trim < 1) {
        min.c.ind <- trunc(model$dim$N * trim)
    } else {
        min.c.ind <- trim
    }
    max.c.ind <- model$dim$N - min.c.ind

    min.c <- model$data$S[order(model$data$S)][min.c.ind]
    max.c <- model$data$S[order(model$data$S)][max.c.ind]

    if (gap < 1) {
        gap <- trunc(model$dim$N * gap)
    }

    grid.data <- get.grid(gamma.limits,
                          c(min.c, max.c),
                          points)

    progress.bar <- txtProgressBar(max = ncombinations(n = points,
                                                       k = m - 1),
                                   style = 3)
    progress <- function(n) setTxtProgressBar(progress.bar, n)

    cores <- detectCores()
    cluster <- makeCluster(max(cores - 1, 1), type = "SOCK")
    registerDoSNOW(cluster)

    c.iter <- icombinations(n = points, k = m - 1)

    grid.result <- foreach(cc = c.iter,
                           .combine = function(A, B) {
                               if (A$SSR < B$SSR) A else B
                           },
                           .packages = c("arrangements", "foreach", "vstar"),
                           .options.snow = list(progress = progress)) %dopar% {
        cc <- drop(cc)

        if (length(cc) > 1 && !all(diff(cc) > gap)) {
            list(SSR = Inf)
        }

        c.vals <- grid.data$c[cc]
        g.iter <- ipermutations(n = points, k = m - 1)

        foreach(gg = g.iter,
                     .combine = function(A, B) {
                        if (A$SSR < B$SSR) A else B
                     }) %do% {
            gg <- drop(gg)

            g.vals <- grid.data$g[gg]

            BM <- get.B.mat(model, g.vals, c.vals, G.func)

            vec.E  <- vec(t(model$data$Y)) - BM$M %*% BM$vec.B

            SSR <- drop(t(vec.E) %*% vec.E)

            list(vec.B = BM$vec.B,
                 vec.E = vec.E,
                 M = BM$M,
                 SSR = SSR,
                 g = g.vals,
                 c = c.vals)
        }
    }

    stopCluster(cluster)

    final.est <- get.estimates(grid.result, model, g.function)

    result <- list(coef = final.est$coef,
                   sd = final.est$sd,
                   t.stat = final.est$t.stat,
                   g = grid.result$g,
                   c = grid.result$c,
                   fitted.values = final.est$fitted.values,
                   residuals = final.est$residuals,
                   cov = final.est$cov,
                   g.function = final.est$g.function,
                   estimates = grid.result,
                   params = model$params,
                   dim = model$dim,
                   data = model$data)

    class(result) <- "vstar"

    return(result)
}
