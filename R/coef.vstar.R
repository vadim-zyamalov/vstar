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
