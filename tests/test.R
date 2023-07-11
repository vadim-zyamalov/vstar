library(vstar)

N <- 1000

x <- matrix(rnorm(2), nrow = 1)
x <- cbind(x, 1)

s <- rnorm(N)
thr <- median(s)
s <- t(t(c(NA, s))) %x% t(vstar:::unity(1))

gg <- vstar:::unity(N + 1) %x% t(.5)
cc <- vstar:::unity(N + 1) %x% t(thr)

G.func <- vstar:::get.G.function("L")

G.mat <- G.func(s, gg, cc)
G.mat <- G.mat[-1, ,drop = FALSE]

A1 <- matrix(c(.5, -.3, .2, .1), ncol = 2)
A1 <- cbind(A1, t(t(c(.1, .2))))
A2 <- matrix(c(.2, .3, -.8, .8), ncol = 2)
A2 <- cbind(A2, t(t(c(.2, .1))))
B <- cbind(t(A1), t(A2))

# --- Только в рамках теста ---
for (k in 1:N) {
    PSI <- t(t(c(1, G.mat[k, ]))) %x% diag(2)
    x <- rbind(
        x,
        cbind(t(t(PSI) %*% t(B) %*% t(x[nrow(x), , drop = FALSE]) + .01 * t(t(rnorm(2)))), 1)
    )
}
y <- x[52:(N + 1), 1:2, drop = FALSE]
colnames(y) <- c("y1", "y2")
s <- s[52:(N + 1)]

y2 <- y[1:50, , drop = FALSE]
s2 <- s[1:50]

model <- vstar.prepare(endo = c("y1", "y2"), const = TRUE, p = 1, trans = s, dataset = y)

res <- vstar.grid(dataset = model, m = 2, gamma.limits = c(0.1, 1), points = 10)

res2 <- vstar.nls(res, tol = 1e-6)

yf2 <- vstar:::predict.vstar(res2, y2, s2)

