get.grid <- function(g.limits, c.limits, points) {
    return(list(g = seq(from = min(g.limits),
                        to = max(g.limits),
                        length.out = points),
                c = seq(from = min(c.limits),
                        to = max(c.limits),
                        length.out = points)))
}
