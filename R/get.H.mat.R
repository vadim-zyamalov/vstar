get.H.mat <- function(B.mat) {
    result <- list()
    for (i in seq_len(length(B.mat))) {
        elems <- rep(0, length(B.mat))
        elems[i] <- 1
        result[[i]] <- matrix(elems, ncol = ncol(B.mat))
    }
    return(result)
}
