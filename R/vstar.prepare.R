#' @title
#' Data preparation procedure
#'
#' @description
#' This function prepares data to be used in the estimation and testing
#' procedures.
#'
#' @details
#' All variables names should be either data frames or matrix column names.
#'
#' @param endo a vector of endogenous variables names.
#' @param exog a vector of exogenous variables names.
#' @param trans a transition variable. It can be either
#' * a name of transition variable.
#' * a list of 2 values, name of a variable and its lag, respectively.
#' * a custom variable as a vector or a column matrix.
#' @param const should an intercept be included.
#' @param trend should a trend be included.
#' @param season should season dummies be included.
#' This parameter should be eigher NULL of interger number of seasons.
#' @param p a number of lags of autoregressive processes.
#' @param coint.beta an optional matrix of cointegrating vectors.
#' These should be pre-estimated.
#' This follows the idea from (Rothman, et al., 1999).
#' Cointegrating vectors should be located by rows.
#' Coefficients of constant term and trend (if any) should be the in the last
#' columns!
#' @param coint.const should an intercept be included into cointegration
#' relations.
#' @param coint.trend should a trend be included into cointegration relations.
#' @param na.action an action that should be applied to NA values.
#' @param dataset a data frame or a matrix with the data.
#'
#' @return
#' An object of S3-class `vstar.data` containing three sublists.
#'
#' The first one contains all input parameters except `p`, `na.action`, and
#' `dataset`.
#'
#' The second one contains dimensions of the model: the final number of
#' observations, number of lags, and number of endogenous variables.
#'
#' The last one conatins matrices of LHS \eqn{Y}, RHS \eqn{X},
#' trnsition variable \eqn{S}, and time \eqn{t}. The latter is not used at the
#' moment.
#'
#' @references
#' P. Rothman, D. J. C. van Dijk, and P. H. B. F. Franses,
#' “A multivariate STAR analysis of the relationship between money and output,”
#' Erasmus University Rotterdam, Erasmus School of Economics (ESE),
#' Econometric Institute, EI 9945-/A, Nov. 1999.
#'
#' T. Teräsvirta and Y. Yang,
#' “Specification, estimation and evaluation of
#' vector smooth transition autoregressive models with applications,”
#' Department of Economics and Business Economics, Aarhus University,
#' 2014–08, Mar. 2014.
#'
#' @export
vstar.prepare <- function(endo,
                          exog = NULL,
                          trans,
                          const = FALSE,
                          trend = FALSE,
                          season = NULL,
                          p = 1,
                          coint.beta = NULL,
                          coint.const = FALSE,
                          coint.trend = FALSE,
                          na.action = na.omit,
                          dataset) {
    inner.p <- p
    if (!is.null(coint.beta)) inner.p <- inner.p - 1

    if (is.data.frame(dataset)) {
        inner.y <- as.matrix(dataset[endo])
    } else {
        inner.y <- dataset[, endo, drop = FALSE]
    }
    names.y   <- endo

    inner.exo <- as.matrix(dataset[exog])

    inner.t   <- seq_len(nrow(inner.y))

    if (typeof(trans) == "character" && length(trans) == 1) {
        inner.s <- as.matrix(dataset[trans])
    } else if (typeof(trans) == "list") {
        if (is.data.frame(dataset)) {
            inner.s <- lagn(dataset[trans[[1]]], trans[[2]])
        } else {
            inner.s <- lagn(dataset[, trans[[1]], drop = FALSE], trans[[2]])
        }
    } else {
        inner.s <- as.matrix(trans)
    }

    inner.x <- NULL
    names.x <- NULL

    if (!is.null(coint.beta)) {
        tmp.y <- inner.y

        if (coint.const) {
            tmp.y <- cbind(tmp.y, 1)
        }
        if (coint.trend) {
            tmp.y <- cbind(tmp.y, seq_len(nrow(tmp.y)))
        }

        inner.y <- diffn(inner.y)
        names.y <- paste("d.", endo, sep = "")

        inner.x <- lagn(tmp.y %*% t(coint.beta), 1)
        names.x <- c(names.x,
                     paste("(z", seq_len(nrow(coint.beta)), ")", sep = ""))
    }

    if (inner.p > 0) for (l in 1:inner.p) {
        inner.x <- cbind(inner.x, lagn(inner.y, l))
        names.x <- c(names.x, paste("l", l, ".", names.y, sep = ""))
    }

    if (const) {
        inner.x <- cbind(inner.x, 1)
        names.x <- c(names.x, "(C)")
    }
    if (trend) {
        inner.x <- cbind(inner.x, seq_len(nrow(inner.x)))
        names.x <- c(names.x, "(t)")
    }
    if (!is.null(season)) {
        tmp.seas <- unity(ceiling(nrow(inner.x) / season)) %x% diag(season)
        tmp.seas <- tmp.seas[seq_len(nrow(inner.x)), -1, drop = FALSE]

        inner.x <- cbind(inner.x, tmp.seas)
        names.x <- c(names.x, paste("(s", 2:season, ")", sep = ""))
    }

    if (!is.null(exog)) {
        inner.x <- cbind(inner.x, inner.exo)
        names.x <- c(names.x, exog)
    }

    # Delete NA's
    if (!is.null(na.action)) {
        tmp.data <- cbind(inner.y, inner.x, inner.s, inner.t)
        tmp.data <- na.action(tmp.data)

        inner.t <- tmp.data[, ncol(tmp.data)]
        tmp.data <- tmp.data[, -ncol(tmp.data), drop = FALSE]

        inner.s <- tmp.data[, ncol(tmp.data), drop = FALSE]
        tmp.data <- tmp.data[, -ncol(tmp.data), drop = FALSE]

        inner.y <- tmp.data[, seq_along(endo), drop = FALSE]
        inner.x <- tmp.data[, -seq_along(endo), drop = FALSE]
    }

     colnames(inner.y) <- names.y
     colnames(inner.x) <- names.x
     if (is.character(trans)) {
        colnames(inner.s) <- trans
     } else if (typeof(trans) == "list") {
         colnames(inner.s) <- trans[[1]]
     }

    dimensions <- list(p = p,
                       p.fixed = inner.p,
                       k = length(endo),
                       N = nrow(inner.y))

    params <- list(endo = endo,
                   exog = exog,
                   trans = trans,
                   const = const,
                   trend = trend,
                   season = season,
                   coint.beta = coint.beta,
                   coint.const = coint.const,
                   coint.trend = coint.trend)

    data <- list(Y = inner.y,
                 X = inner.x,
                 S = inner.s,
                 t = inner.t)

    result <- list(params = params,
                   dim = dimensions,
                   data = data)

    class(result) <- "vstar.data"

    return(result)
}
