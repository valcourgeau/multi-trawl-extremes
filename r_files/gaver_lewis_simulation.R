# Due to Walker (2000)

rlgprocess <- function(alpha, beta, rho, kappa, timesteps, n){
  if(rho >= 1 | rho <= 0) stop("Rho should be between 0 and 1.")
  
  latent_lambda <- rgamma(n = timesteps*n,
                     shape = alpha,
                     rate = 1.0)
  latent_lambda <- rpois(n = timesteps*n,
                    lambda = latent_lambda*(1-rho)/rho)
  latent_lambda <- rgamma(n = timesteps*n, 
                     shape = latent_lambda,
                     rate = beta/rho)
  latent_lambda <- rho*rgamma(n = timesteps*n,
                          shape = alpha,
                          rate = beta) + latent_lambda
  results <- rep(0, timesteps*n)
  which_to_use <- which(runif(n=timesteps*n) <= exp(-kappa * latent_lambda))
  results[which_to_use] <- rexp(n = length(which_to_use),
                                rate = latent_lambda[which_to_use])
  return(matrix(results, ncol=n))
}
setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source('gaver_lewis_simulation.R')
library(evir)
lgp_sim <- rlgprocess(alpha = 3.00,
                      beta = 1.0,
                      kappa = 1.7,
                      rho = 0.4,
                      timesteps = 10000,
                      n=1)
plot(lgp_sim[,1], type = 'l')
lines(lgp_sim[,2], type = 'l', col="red")
1000000-length(which(lgp_sim == 0.0))
summary(lgp_sim[,1])
gpd(lgp_sim, nextremes = 10000-length(which(lgp_sim == 0.0)))
c(1/3,  (1+1.7)/3)

plot(density(lgp_sim[lgp_sim>0]))
lines(density(evir::rgpd(10000, xi=0.33, beta=(1+1.7)/3)), col = "red")

acf(lgp_sim[lgp_sim>0])

