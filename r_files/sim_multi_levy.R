require("ghyp")


n <- 5
p <- qr.Q(qr(matrix(rnorm(n^2), n)))
# NIG
nig <- NIG(chi = 2, psi = 2, mu = 0, sigma = diag(rep(1, length(mu))),
           gamma = rep(0, length(mu)), alpha.bar = NULL, data = NULL)
