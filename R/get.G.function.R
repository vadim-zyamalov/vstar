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

get.G.derivative <- function(g.function = "L") {
    G.func <- get.G.function(g.function)
    if (is.character(g.function)) {
        if (g.function == "L") {
            g.func <- function(s, g, thr, sd = 1) {
                ((s - thr) / sd) * G.func(s, g, thr) * (1 - G.func(s, g, thr))
            }
            c.func <- function(s, g, thr, sd = 1) {
                -g * G.func(s, g, thr) * (1 - G.func(s, g, thr))
            }
        } else if (g.function == "E") {
            g.func <- function(s, g, thr, sd = 1) {
                ((s - thr) / sd)^2 * (1 - G.func(s, g, thr))
            }
            c.func <- function(s, g, thr, sd = 1) {
                -g * ((s - thr) / sd) * (1 - G.func(s, g, thr))
            }
        } else {
            stop(paste("G: Unknown function \"", g.function, "\"!", sep = ""))
        }
    } else if (is.function(g.function)) {
        stop("G: Provide array of derivative functions!")
    } else {
        stop("G: Unknown transition function!")
    }

    return(list(g = g.func, c = c.func))
}
