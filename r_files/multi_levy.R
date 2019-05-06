require("ghyp")

set.seed(42)
n <- 30
R <- matrix(runif(n^2), ncol=n) 
RtR <- R %*% t(R) 
Q <- cov2cor(RtR) 

mu <- rep(0, length(Q[1,]))
nig <- NIG(chi = 2, psi = 2, mu = mu, sigma = Q,
            gamma = rep(0, length(mu)))
ghyp::rghyp(n=50, object = nig)

Q
rghyp(n = 10000, object = nig)


std.t <- student.t(nu = 3.5, mu = mu, sigma = Q,
            gamma = rep(0, length(mu)), data = NULL)
rghyp(n = 10000, object = std.t)


vg <- VG(lambda = 1, mu = mu, sigma = Q,
          gamma = rep(0, length(mu)), data = NULL)
rghyp(n = 10000, object = vg)

# Compound Poisson
# Conditioned on the number jumps: we KNOW the number of jumps:
# - ordered in time (uniforms);
# - Normally distributed jump-size;
# - Make first moments comparable;
# - Fit to the standard data, parametric bootstrap;
