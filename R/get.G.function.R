get.G.function <- function(g.function = "L") {
    if (is.character(g.function)) {
        if (g.function == "L") {
            func <- function(s, g, thr, sd = 1) {
                1 / (1 + exp(-g * (s - thr) / sd))
            }
        } else if (g.function == "E") {
            func <- function(s, g, thr, sd = 1) {
                1 - exp(-g * ((s - thr) / sd)^2)
            }
        } else {
            stop(paste("G: Unknown function \"", g.function, "\"!", sep = ""))
        }
    } else if (is.function(g.function)) {
        func <- g.function
    } else {
        stop("G: Unknown transition function!")
    }

    return(func)
}
