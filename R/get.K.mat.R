#' @title
#' Get matrix of the models first derivatives.
#'
#' @description
#' An auxiliary function for estimating the matrix \eqn{K} of the first
#' derivatives of the model by its parameters \eqn{\{B, \Omega, \Gamma, C\}}.
#'
#' This matrix is used in tests from (Teräsvirta and Yang, 2014).
#'
#' @param model an object of S3-class `vstar`.
#' @param J a number of first observations to be dropped.
#'
#' @return
#' A matrix of first derivatives by \eqn{\{B, \Gamma, C\}}.
#' \eqn{\Omega} is excluded as all derivatives by it are zero.
#'
#' @references
#' T. Teräsvirta and Y. Yang,
#' “Linearity and misspecification tests for vector smooth transition
#' regression models,”
#' Department of Economics and Business Economics, Aarhus University,
#' 2014–04, Feb. 2014.
#'
#' @keywords internal
get.K.mat <- function(model, J = 0) {
  if (!"vstar" %in% class(model)) {
    stop("Wrong `model`: an object with estimated VSTAR model is needed!")
  }
  if (!"coef" %in% names(model)) {
    stop("K: Please, estimate the model first!")
  }

  k <- model$dim$k
  m <- model$dim$m
  N <- model$dim$N

  H.mat <- get.H.mat(model$coef)

  G.func <- model$func$g.function
  G.deriv <- model$func$g.derivative

  sm <- model$data$S %x% t(unity(m - 1))
  gm <- unity(N) %x% t(model$g)
  cm <- unity(N) %x% t(model$c)

  G.mat <- G.func(sm, gm, cm)
  G.g.mat <- G.deriv[[1]](sm, gm, cm)
  G.c.mat <- G.deriv[[2]](sm, gm, cm)

  K.mat <- NULL

  for (i in (J + 1):N) {
    K.t <- NULL
    PSI <- t(t(c(1, G.mat[i, ]))) %x% diag(k)
    x <- t(model$data$X[i, , drop = FALSE])

    for (j in seq_len(length(model$coef))) {
      K.t <- cbind(
        K.t,
        t(PSI) %*% t(H.mat[[j]]) %*% x
      )
    }

    for (j in 1:(m - 1)) {
      elems <- rep(0, m - 1)
      elems[j] <- G.g.mat[i, j]
      G.g <- t(c(0, elems)) %x% diag(k)

      elems[j] <- G.c.mat[i, j]
      G.c <- t(c(0, elems)) %x% diag(k)

      K.t <- cbind(
        K.t,
        G.g %*% t(model$coef) %*% x,
        G.c %*% t(model$coef) %*% x
      )
    }

    K.mat <- rbind(K.mat, t(vec(K.t)))
  }

  return(K.mat)
}
