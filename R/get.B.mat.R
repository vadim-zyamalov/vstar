#' @title
#' Get coefficients estimates
#'
#' @description
#' An auxiliary function returning vectorized matrix of estimated coefficients
#' \eqn{\operatorname{vec}\hat{B}} and the matrix of modified regressors
#' \eqn{M}.
#'
#' The code in this function is separated from the other code because all
#' lines of code are densely tied and enter the codebase as the
#' single chunk.
#'
#' In this function implemented is the joint OLS estimation procedure from
#' (Teräsvirta and Yang, 2014).
#'
#' @param dataset an object of S3-classes `vstar.data` or `vstar`.
#' @param m a number of regimes in the model.
#' @param g a vector of smoothness parameter values.
#' @param thr a vector of threshold values.
#' @param g.function a transition function to be used
#' * "L": Logistic transition function
#' \deqn{\frac{1}{1 + \exp (-\gamma (s_t - c))}}
#' * "E": Exponential transition function
#' \deqn{1 - \exp(-\gamma (s_t - c)^2)}
#' * fn: Custom function receiving transition variable, gammas and threshold
#' values in exactly this order.
#'
#' @return
#' A list containing vectorized (!) matrix of estimated coefficients
#' \eqn{\operatorname{vec}\hat{B}} and the matrix of modified regressors
#' \eqn{M}.
#'
#' @references
#' T. Teräsvirta and Y. Yang,
#' “Specification, estimation and evaluation of
#' vector smooth transition autoregressive models with applications,”
#' Department of Economics and Business Economics, Aarhus University,
#' 2014–08, Mar. 2014.
#'
#' @keywords internal
get.B.mat <- function(dataset, m, g, thr, g.function = "L") {
  k <- dataset$dim$k
  N <- dataset$dim$N
  Nx <- ncol(dataset$data$X)

  G.func <- get.G.function(g.function)

  sm <- dataset$data$S %x% t(unity(m - 1))
  gm <- unity(N) %x% t(g)
  cm <- unity(N) %x% t(thr)

  G.mat <- G.func(sm, gm, cm)

  M <- (cbind(unity(N), G.mat) %x% diag(k) %x% t(unity(Nx))) *
    (t(unity(k * m)) %x% dataset$data$X %x% unity(k))

  vec.B <- ginv(t(M) %*% M) %*% t(M) %*% vec(t(dataset$data$Y))

  return(list(
    vec.B = vec.B,
    M = M
  ))
}
