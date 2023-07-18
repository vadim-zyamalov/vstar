#' @title
#' Plotting Generalized Impulse Responses
#'
#' @param res an object of S3-class `vstar.irf`.
#'
#' @export
plot.vstar.irf <- function(res, ...) {
    k <- ncol(res$pos)
    n <- res$n

    names.y <- colnames(res$pos)

    graphs <- list()

    for (i in 1:k) {
        plot.data <- cbind(
            res$pos[, i], res$neg.sd[, i],
            res$neg[, i], res$pos.sd[, i]
        )
        colnames(plot.data) <- c("pos", "pos.sd", "neg", "neg.sd")
        rownames(plot.data) <- 0:n

        graphs[[paste("g", i, sep = "")]] <-
            ggplot(
                data = as.data.frame(plot.data),
                mapping = aes(x = 0:n)
            ) +
            geom_line(aes(y = pos, colour = "Positive")) +
            geom_line(aes(y = pos + 2 * pos.sd),
                linetype = "dashed", colour = "red"
            ) +
            geom_line(aes(y = pos - 2 * pos.sd),
                linetype = "dashed", colour = "red"
            ) +
            geom_line(aes(y = neg, colour = "Negative")) +
            geom_line(aes(y = neg + 2 * neg.sd),
                linetype = "dashed", colour = "blue"
            ) +
            geom_line(aes(y = neg - 2 * neg.sd),
                linetype = "dashed", colour = "blue"
            ) +
            scale_colour_manual(
                name = "Shock",
                values = c("blue", "red")
            ) +
            scale_x_continuous(breaks = 0:n) +
            labs(
                title = paste(
                    "Response of", names.y[i],
                    "to", res$shock
                ),
                x = "Period",
                y = paste("Response to", res$shock)
            )
    }

    plot_grid(
        plotlist = graphs,
        align = "v",
        ncol = floor(sqrt(k))
    )
}
