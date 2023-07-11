#' @title
#' Extract estimated parameters from VSTAR model
#'
#' @param model an object of S3-class `vstar`.
#'
#' @return
#' List containing:
#' * matrix \eqn{B}.
#' * vector of \eqn{\gamma} values.
#' * vector of threshold values.
#'
#' @keywords internal
#' @export
coef.vstar <- function(model, ...) {
    if (!"vstar" %in% class(model)) {
        stop("Wrong `model`: an object with estimated VSTAR model is needed!")
    }
    if (!"coef" %in% names(model)) {
        stop("Estimate the model first!")
    }

    return(list(coef = model$coef,
                gamma = model$estimates$g,
                c = model$estimates$c))
}
