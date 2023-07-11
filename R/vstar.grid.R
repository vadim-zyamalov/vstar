#' @title
#' Estimate VSTAR using grid search
#'
#' @description
#' This function estimates the parameters of VSTAR model using
#' the provided `vstar.data` object with the data prepared.
#'
#' @details
#'
#' As stated by (Teräsvirta and Yang, 2014), estimating VSTAR directly by NLS
#' is quite hard due to the big number of parameters.
#' They proposed to estimate the model by grid search.
#' The main idea is that with \eqn{\gamma} and \eqn{c} fixed the system
#' becomes linear, and it can be estimated by OLS.
#' So one can go through a grid of possible combinations of \eqn{\gamma}
#' and \eqn{c} values searching for the one with the minimum value of SSR.
#'
#' It should be noted that the number of steps increases substantially as
#' the number of regimes \eqn{m} grows. Let's set the number of steps along one
#' directon of the grid to 100 and denote it as \eqn{p}.
#' If we have only two regimes then the number of steps will be
#' \deqn{p^2 = 10,000}
#'
#' But if we have three regimes, then the number of steps will be
#' \deqn{p^2 \times \frac{p(p+1)}{2} = \frac{1}{2}(p^4 + p^3) = 50,500,000}
#' The second factor in the first expression has this form due to the
#' restriction \eqn{c_1 < c_2}, and actually it's an area of a "dicrete"
#' triangle.
#'
#' If we have four regimes, then the number of steps will be
#' \deqn{p^3 \times \frac{1}{6}p(p+1)(p+2) = 171,700,000,000}
#'
#' It's quite easy to show, that for an arbitrary number of regimes the number
#' of step will be
#' \eqn{p^{m-1} \times \frac{1}{m!} p(p+1)(p+2)\dots(p+m-1)}
#'
#' So if you want to estimate a model with a big number of regimes you should be
#' sure you know what you do!
#'
#' @param dataset an object of S3-classes `vstar.data` or `vstar`.
#' @param m a number of regimes. It should be not less than 2.
#' @param g.function a transition function to be used
#' * "L": Logistic transition function
#' \deqn{\frac{1}{1 + \exp (-\gamma (s_t - c))}}
#' * "E": Exponential transition function
#' \deqn{1 - \exp(-\gamma (s_t - c)^2)}
#' * fn: Custom function receiving transition variable, gammas and threshold
#' values in exactly this order.
#' @param g.derivative a vector of transition functions derivatives:
#' the first is by \eqn{\gamma}, and the second is by \eqn{c}.
#' If you use "L" or "E" as `g.function` then you shouldn't bother with this
#' parameter.
#' But if you provide a custom derivative, then you must provide
#' custom derivatives as well.
#' Although they are not needed by the estimation itself,
#' they're stored for future testing.
#' If you're not going to make diagnostic testing just pass NULL.
#' @param gamma.limits a vector of lower and upper bounds of \eqn{\gamma}.
#' @param points a number of points in grid sequences.
#' @param trim a number (if greater than 1) or a share of the transition
#' variable values (in ascending order) to be ignored during grid construction.
#' Affects only the \eqn{c} part of the grid.
#' @param gap a number (if greater than 1) or a share of minimum distance
#' between two values of \eqn{c} in grid.
#' Used only if \eqn{m > 2}.
#'
#' @return
#' An object of S3-class `vstar` containing
#' * everything from the input object `dataset`. To the sublist `dim` added is
#' the number of regimes \eqn{m}.
#' * a sublist `estimates` containing the result of the loop over the grid.
#' * a sublist `func` containing a transition function and its derivatives.
#' NB: these are function bodies!
#' * a matrix of coefficients.
#' * a matrix of standard deviations of coefficients.
#' * a matrix of t-statistics of coefficients.
#' * a vector of estimated values of \eqn{\gamma}.
#' * a vector of estimated values of \eqn{c}.
#' * a matrix of fitted values.
#' * a matrix of residuals.
#' * covariance matrix of residuals.
#' * a value of `g.function` parameter.
#'
#' @references
#' T. Teräsvirta and Y. Yang,
#' “Specification, estimation and evaluation of
#' vector smooth transition autoregressive models with applications,”
#' Department of Economics and Business Economics, Aarhus University,
#' 2014–08, Mar. 2014.
#'
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
                       m = 2,
                       g.function = "L",
                       gamma.limits = c(.1, 100),
                       points = 200,
                       trim = .15,
                       gap = .1) {

    p <- model$dim$p
    p.fixed <- model$dim$p.fixed
    k <- model$dim$k
    N <- model$dim$N

    if (!"vstar.data" %in% class(model)) {
        stop("Wrong `model`: an object with prepared data is needed!")
    }
    if (m < 2) {
        stop("`m` too low: at least 2 regimes is needed!")
    }

    G.func <- get.G.function(g.function)

    if (trim < 1) {
        min.c.ind <- trunc(N * trim)
    } else {
        min.c.ind <- trim
    }
    max.c.ind <- N - min.c.ind

    min.c <- model$data$S[order(model$data$S)][min.c.ind]
    max.c <- model$data$S[order(model$data$S)][max.c.ind]

    if (gap < 1) {
        gap <- trunc(N * gap)
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

            BM <- get.B.mat(model, m, g.vals, c.vals, G.func)

            vec.E  <- vec(t(model$data$Y)) - BM$M %*% BM$vec.B

            SSR <- drop(t(vec.E) %*% vec.E)

            list(vec.B = BM$vec.B,
                 vec.E = vec.E,
                 M = BM$M,
                 SSR = SSR,
                 m = m,
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
                   dim = list(p = p,
                              p.fixed = p.fixed,
                              m = m,
                              k = k,
                              N = N),
                   data = model$data)

    class(result) <- "vstar"

    return(result)
}
