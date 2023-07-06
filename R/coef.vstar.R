#' @keywords internal
#' @export
coef.vstar <- function(model, ...) {
    if (!"coef" %in% names(model)) {
        stop("Estimate the model first!")
    }
    return(list(coef = model$coef,
                gamma = model$estimates$g,
                c = model$estimates$c))
}
