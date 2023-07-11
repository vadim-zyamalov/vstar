#' @title
#' Printing the results of the model estimation
#'
#' @param model an object of S3-class `vstar`.
#'
#' @keywords internal
#'
#' @export
summary.vstar <- function(model, ...) {
    if (!"vstar" %in% class(model)) {
        stop("Wrong `model`: an object with estimated VSTAR model is needed!")
    }
    if (!"coef" %in% names(model)) {
        stop("Estimate the model first!")
    }

    sum.func <- ""
    if (is.character(model$g.function)) {
        sum.func <- paste(model$g.function, "-", sep = "")
    }

    cat("Results of ", sum.func, "VSTAR model estimation\n\n", sep = "")
    cat("Number of regimes:     ", model$dim$m, "\n")
    cat("Number of lags:        ", model$dim$p.fixed, "\n")
    cat("Number of observations:", model$dim$N, "\n")
    cat("Endogenous variables:\n", model$params$endo, "\n")
    if (!is.null(model$params$exog)) {
        cat("Exogenous variables:\n", model$params$exog, "\n")
    }
    if (length(model$params$trans) == 1) {
        cat("Transition variable:   ", model$params$trans, "\n")
    } else {
        cat("Transition variable:    custom\n")
    }

    cat("\nAdditional model components:\n")
    cat("Constant included:              ", model$params$const, "\n")
    cat("Trend included:                 ", model$params$trend, "\n")
    if (!is.null(model$params$season)) {
        cat("Number of seasons:              ", model$params$season, "\n")
    }
    if (!is.null(model$params$coint.beta)) {
        cat("Number of cointegration vectors:",
            nrow(model$params$coint.beta), "\n")
        cat("Constant in cointegration:      ",
            nrow(model$params$coint.const), "\n")
        cat("Trend in cointegration:         ",
            nrow(model$params$coint.trend), "\n")
    }

    sum.func <- "custom"
    if (is.character(model$g.function)) {
        if (model$g.function == "L") sum.func <- "logistic"
        else if (model$g.function == "L") sum.func <- "exponent"
    }
    cat("\nTransition function:", sum.func, "\n")
    cat("gamma values:\n")
    print(model$g)
    cat("threshold values:\n")
    print(model$c)

    cat("\nMatrix of coefficients:\n")
    print(model$coef)

    cat("\nMatrix of coefficients SD:\n")
    print(model$sd)

    cat("\nResiduals covariance matrix:\n")
    print(model$cov)
}

#' @keywords internal
#' @export
print.vstar <- function(model, ...) {
    summary(model)
}
