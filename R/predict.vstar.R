#' @title
#' Model Predictions
#'
#' @description
#' The function returns the fitted values of the estimated model.
#' If `new.data` is not NULL then the fitted values for the new data are
#' calculated.
#'
#' @param model an object of S3-classes `vstar`.
#' @param new.data if not NULL, use this new data as a source of variables for
#' prediction.
#' @param new.trans if not NULL, use this new transition variable for
#' prediction. Use this if the model was estimated with custom variable.
#' @param na.action an action that should be applied to NA values.
#'
#' @return
#' A matrix of fitted.values.
#'
#' @export
predict.vstar <- function(model,
                          new.data = NULL,
                          new.trans = NULL,
                          na.action = na.omit, ...) {
  if (!"vstar" %in% class(model)) {
    stop("Wrong `model`: an object with estimated VSTAR model is needed!")
  }

  if (is.null(new.data)) {
    return(model$fitted.values)
  }

  if (!"vstar.data" %in% class(new.data)) {
    new.data <- vstar.prepare(
      endo = model$params$endo,
      exog = model$params$exog,
      trans = new.trans,
      const = model$params$const,
      trend = model$params$trend,
      season = model$params$season,
      p = model$dim$p,
      coint.beta = model$params$coint.beta,
      coint.const = model$params$coint.const,
      coint.trend = model$params$coint.trend,
      na.action = na.action,
      dataset = new.data
    )
  }

  BM.mat <- get.B.mat(
    new.data,
    model$dim$m,
    model$g,
    model$c,
    model$g.function
  )

  result <- t(matrix(BM.mat$M %*% vec(model$coef), ncol = new.data$dim$N))

  return(result)
}
