#' @title
#' Plotting Generalized Impulse Responses
#'
#' @param res an object of S3-class `vstar.irf`.
#'
#' @export
plot.vstar.irf <- function(res, ...) {
  k <- length(res$pos)
  n <- res$n

  names.y <- names(res$pos)

  graphs <- list()

  for (i in 1:k) {
    plot.data <- cbind(res$pos[[i]], res$neg[[i]])
    colnames(plot.data) <- c(
      paste("pos", c("lowest", "low", "mid", "high", "highest"), sep = "."),
      paste("neg", c("lowest", "low", "mid", "high", "highest"), sep = ".")
    )
    rownames(plot.data) <- 0:n

    graphs[[paste("g", i, sep = "")]] <-
      ggplot(
        data = as.data.frame(plot.data),
        mapping = aes(x = 0:n)
      ) +
      geom_line(aes(y = pos.mid, colour = "Positive")) +
      geom_ribbon(aes(ymin = pos.low, ymax = pos.high),
        fill = "red", alpha = .2
      ) +
      geom_line(aes(y = pos.lowest), linetype = "dashed", colour = "red") +
      geom_line(aes(y = pos.highest), linetype = "dashed", colour = "red") +
      geom_line(aes(y = neg.mid, colour = "Negative")) +
      geom_ribbon(aes(ymin = neg.low, ymax = neg.high),
                  fill = "blue", alpha = .2
      ) +
      geom_line(aes(y = neg.lowest), linetype = "dashed", colour = "blue") +
      geom_line(aes(y = neg.highest), linetype = "dashed", colour = "blue") +
      scale_colour_manual(name = "Shock", values = c("blue", "red")) +
      scale_x_continuous(breaks = 0:n) +
      labs(title = paste("Response of", names.y[i], "to", res$shock),
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
