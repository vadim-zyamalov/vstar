#' @title
#' Improve VSTAR estimates using NLS
#'
#' @description
#' This function estimates the parameters of VSTAR model by NLS using
#' the provided `vstar` object as the source of data and starting values.
#'
#' @details
#' As stated by (Teräsvirta and Yang, 2014), estimating VSTAR directly by NLS
#' is quite hard due to the big number of parameters.
#' They proposed to estimate the model by grid search.
#' The results of this grid search can be used as the starting values for NLS.
#'
#' The procedure goes as follows:
#' * first, \eqn{\gamma} and \eqn{c} are reestimated using \link[stats]{optim}
#' with "Nelder-Mead" algorithm.
#' * then, with these new values \eqn{B} is reestimated.
#' * Repeat until convergence is reached.
#'
#' @param model an object of S3-classes `vstar`.
#' @param tol a level of tolerance. If the vector of \eqn{\gamma} and \eqn{c}
#' changes less than this value the process stops.
#' @param iter the maximum number of iterations.
#' @param verbose should intermediate steps be shown?
#'
#' @return
#' An object of S3-class `vstar`. See \link{vstar.grid} for details.
#'
#' @references
#' T. Teräsvirta and Y. Yang,
#' “Specification, estimation and evaluation of
#' vector smooth transition autoregressive models with applications,”
#' Department of Economics and Business Economics, Aarhus University,
#' 2014–08, Mar. 2014.
#'
#' @export
vstar.nls <- function(model,
                      tol = 1e-6,
                      iter = 100,
                      verbose = FALSE) {
  if (!"vstar" %in% class(model)) {
    stop("Wrong `model`: an object with estimated VSTAR model is needed!")
  }

  m <- model$dim$m

  vec.y <- vec(t(model$data$Y))
  G.func <- model$func$g.function

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
    if (verbose) {
      cat("Step ", step, ": ", sep = "")
    }

    iter.nls <- optim(c(iter.g, iter.c), Q)
    result <- iter.nls$par

    dd <- result - c(iter.g, iter.c)
    delta <- sqrt(dd %*% dd)

    if (verbose) {
      cat("delta = ", delta, "\n", sep = "")
    }

    iter.g <- result[1:(m - 1)]
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

  nls.result <- list(
    vec.B = iter.BM$vec.B,
    vec.E = vec.y - iter.BM$M %*% iter.BM$vec.B,
    M = iter.BM$M,
    SSR = drop(t(model$estimates$vec.E) %*%
      model$estimates$vec.E),
    m = m,
    g = iter.g,
    c = iter.c
  )

  final.est <- get.estimates(nls.result, model, model$g.function)

  result <- list(
    coef = final.est$coef,
    sd = final.est$sd,
    t.stat = final.est$t.stat,
    g = iter.g,
    c = iter.c,
    fitted.values = final.est$fitted.values,
    residuals = final.est$residuals,
    cov = final.est$cov,
    g.function = model$g.function,
    func = model$func,
    estimates = nls.result,
    params = model$params,
    dim = model$dim,
    data = model$data
  )

  class(result) <- "vstar"

  return(result)
}
