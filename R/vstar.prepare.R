#' @export
vstar.prepare <- function(endo,
                          exog = NULL,
                          trans,
                          const = FALSE,
                          trend = FALSE,
                          season = NULL,
                          p = 1,
                          m = 1,
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
        names.x <- c(names.x, paste("z", seq_len(nrow(coint.beta)), sep = ""))
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
        names.x <- c(names.x, paste("s", 2:season, sep = ""))
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
     }

    dimensions <- list(p = p,
                       p.fixed = inner.p,
                       m = m,
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

    class(result) <- "vstar"

    return(result)
}
