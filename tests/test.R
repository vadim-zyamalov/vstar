library(vstar)

## -----------------------------------------------------------------------------
N <- 1000

## -----------------------------------------------------------------------------
s <- rnorm(N)

thr_2reg <- median(s)
thr_3reg <- quantile(s, prob = c(.33, .67))

print(thr_2reg)
print(thr_3reg)

## -----------------------------------------------------------------------------
G.func <- vstar:::get.G.function("L")

sm <- s %x% t(vstar:::unity(1))
gm <- vstar:::unity(N) %x% t(4)
cm <- vstar:::unity(N) %x% t(thr_2reg)

G.mat_2reg <- G.func(sm, gm, cm)

sm <- s %x% t(vstar:::unity(2))
gm <- vstar:::unity(N) %x% t(c(4, 4))
cm <- vstar:::unity(N) %x% t(thr_3reg)

G.mat_3reg <- G.func(sm, gm, cm)

## -----------------------------------------------------------------------------
library(ggplot2)

ggplot(mapping = aes(x = s)) +
    geom_line(mapping = aes(y = G.mat_2reg)) +
    labs(x = "Transition variable", y = "Transition functions")

ggplot(mapping = aes(x = s)) +
    geom_line(mapping = aes(y = G.mat_3reg[,1], colour = "1")) +
    geom_line(mapping = aes(y = G.mat_3reg[,2], colour = "2")) +
    labs(x = "Transition variable", y = "Transition functions") +
    scale_colour_manual(name = "G", values = c("blue", "red"))

## -----------------------------------------------------------------------------
A1 <- .5 * diag(2)
A2 <- matrix(c(.8, .1, .1, .8), ncol = 2)
A3 <- -.5 * diag(2)

B_2reg <- cbind(t(A1), t(A2))
B_2reg_alt <- cbind(t(A1), t(A3))
B_3reg <- cbind(t(A1), t(A2), t(A3))

## -----------------------------------------------------------------------------
y.lin <- matrix(0, ncol = 2, nrow = 1001)
y.lin[1, ] <- rnorm(2)

for (k in 1:N) {
    y.lag <- y.lin[k, , drop = FALSE]
    y.lin[k + 1, ] <- t(A1 %*% t(y.lag) + rnorm(2))
}

y.lin <- y.lin[52:(N + 1), 1:2]
colnames(y.lin) <- c("y1", "y2")

plot.ts(y.lin)

## -----------------------------------------------------------------------------
y.2reg <- matrix(0, ncol = 2, nrow = 1001)
y.2reg[1, ] <- rnorm(2)

for (k in 1:N) {
    y.lag <- y.2reg[k, , drop = FALSE]
    PSI <- t(t(c(1, G.mat_2reg[k, ]))) %x% diag(2)
    y.2reg[k + 1, ] <- t(t(PSI) %*% t(B_2reg) %*% t(y.lag) + rnorm(2))
}
y.2reg <- y.2reg[52:(N + 1), 1:2, drop = FALSE]
colnames(y.2reg) <- c("y1", "y2")

plot.ts(y.2reg)

## -----------------------------------------------------------------------------
y.3reg <- matrix(0, ncol = 2, nrow = 1001)
y.3reg[1, ] <- rnorm(2)

for (k in 1:N) {
    y.lag <- y.3reg[k, , drop = FALSE]
    PSI <- t(t(c(1, G.mat_3reg[k, ]))) %x% diag(2)
    y.3reg[k + 1, ] <- t(t(PSI) %*% t(B_3reg) %*% t(y.lag) + rnorm(2))
}
y.3reg <- y.3reg[52:(N + 1), 1:2, drop = FALSE]
colnames(y.3reg) <- c("y1", "y2")

plot.ts(y.3reg)

## -----------------------------------------------------------------------------
y.acorr <- matrix(0, ncol = 2, nrow = 1001)
y.acorr[1, ] <- rnorm(2)

e <- matrix(rnorm(2002), ncol = 2)

for (k in 1:N) {
    y.lag <- y.acorr[k, , drop = FALSE]
    e[k + 1, ] <- e[k + 1, ] + .8 * e[k, ]
    PSI <- t(t(c(1, G.mat_2reg[k, ]))) %x% diag(2)
    y.acorr[k + 1, ] <- t(t(PSI) %*% t(B_2reg) %*% t(y.lag) +
                              t(e[k + 1, , drop = FALSE]))
}
y.acorr <- y.acorr[52:(N + 1), 1:2, drop = FALSE]
colnames(y.acorr) <- c("y1", "y2")

plot.ts(y.acorr)

## -----------------------------------------------------------------------------
y.chcoef <- matrix(0, ncol = 2, nrow = 1001)
y.chcoef[1, ] <- rnorm(2)

for (k in 1:N) {
    y.lag <- y.chcoef[k, , drop = FALSE]
    PSI <- t(t(c(1, G.mat_2reg[k, ]))) %x% diag(2)
    if (k <= 525) {
        y.chcoef[k + 1, ] <- t(t(PSI) %*% t(B_2reg) %*% t(y.lag) +
                                   rnorm(2))
    } else {
        y.chcoef[k + 1, ] <- t(t(PSI) %*% t(B_2reg_alt) %*% t(y.lag) +
                                   rnorm(2))
    }
}
y.chcoef <- y.chcoef[52:(N + 1), 1:2, drop = FALSE]
colnames(y.chcoef) <- c("y1", "y2")

plot.ts(y.chcoef)

## -----------------------------------------------------------------------------
s <- s[51:N]

## -----------------------------------------------------------------------------
model.lin <- vstar.prepare(endo = c("y1", "y2"),
                           p = 1,
                           trans = s,
                           dataset = y.lin)
model.2reg <- vstar.prepare(endo = c("y1", "y2"),
                            p = 1,
                            trans = s,
                            dataset = y.2reg)
model.3reg <- vstar.prepare(endo = c("y1", "y2"),
                            p = 1,
                            trans = s,
                            dataset = y.3reg)
model.acorr <- vstar.prepare(endo = c("y1", "y2"),
                             p = 1,
                             trans = s,
                             dataset = y.acorr)
model.chcoef <- vstar.prepare(endo = c("y1", "y2"),
                              p = 1,
                              trans = s,
                              dataset = y.chcoef)

## -----------------------------------------------------------------------------
linearity.test(model.lin, J = 1)

## -----------------------------------------------------------------------------
round(linearity.test(model.2reg, J = 1), 4)
round(linearity.test(model.3reg, J = 1), 4)
round(linearity.test(model.acorr, J = 1), 4)
round(linearity.test(model.chcoef, J = 1), 4)

## ---- results = FALSE---------------------------------------------------------
result.2reg <- vstar.grid(dataset = model.2reg,
                          m = 2,
                          gamma.limits = c(3, 5),
                          points = 200,
                          cores = 10)
result.3reg <- vstar.grid(dataset = model.3reg,
                          m = 2,
                          gamma.limits = c(3, 5),
                          points = 200,
                          cores = 10)
result.3reg.ok <- vstar.grid(dataset = model.3reg,
                             m = 3,
                             gamma.limits = c(3, 5),
                             points = 20,
                             cores = 10)
result.acorr <- vstar.grid(dataset = model.acorr,
                           m = 2,
                           gamma.limits = c(3, 5),
                           points = 200,
                           cores = 10)
result.chcoef <- vstar.grid(dataset = model.chcoef,
                            m = 2,
                            gamma.limits = c(3, 5),
                            points = 200,
                            cores = 10)

## -----------------------------------------------------------------------------
result.nls.2reg   <- vstar.nls(result.2reg, tol = 1e-6, verbose = TRUE)
result.nls.3reg   <- vstar.nls(result.3reg, tol = 1e-6, verbose = TRUE)
result.nls.acorr  <- vstar.nls(result.acorr, tol = 1e-6, verbose = TRUE)
result.nls.chcoef <- vstar.nls(result.chcoef, tol = 1e-6, verbose = TRUE)

## -----------------------------------------------------------------------------
summary(result.nls.2reg)
summary(result.nls.3reg)
summary(result.nls.acorr)
summary(result.nls.chcoef)

## -----------------------------------------------------------------------------
round(nonlinearity.test(result.nls.2reg, J = 3), 4)
round(nonlinearity.test(result.nls.3reg, J = 3), 4)
round(nonlinearity.test(result.nls.acorr, J = 3), 4)
round(nonlinearity.test(result.nls.chcoef, J = 3), 4)

## -----------------------------------------------------------------------------
round(serial.correlation.test(result.nls.2reg, J = 1, ortogonalize = TRUE), 4)
round(serial.correlation.test(result.nls.3reg, J = 1, ortogonalize = TRUE), 4)
round(serial.correlation.test(result.nls.acorr, J = 1, ortogonalize = TRUE), 4)
round(serial.correlation.test(result.nls.chcoef, J = 1, ortogonalize = TRUE), 4)

## -----------------------------------------------------------------------------
round(stability.test(result.nls.2reg), 4)
round(stability.test(result.nls.3reg), 4)
round(stability.test(result.nls.acorr), 4)
round(stability.test(result.nls.chcoef), 4)
