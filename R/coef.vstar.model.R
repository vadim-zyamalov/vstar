#' @keywords internal
#' @export
coef.vstar.model <- function(object, ...) {
    return(list(coef = object$coef,
                gamma = object$estimates$g,
                c = object$estimates$c))
}
