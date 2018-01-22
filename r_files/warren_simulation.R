rwprocess <- function(alpha, beta, rho, kappa, timesteps, n){
  if(rho >= 1 | rho < 0) stop("Rho should be between 0 and 1.")
  
  latent_lambda <- rgamma(n = timesteps*n,
                          shape = alpha,
                          rate = beta)
  latent_lambda <- rpois(n = timesteps*n,
                         lambda = rho*latent_lambda*beta/(1-rho))

  latent_lambda <- rho*rgamma(n = timesteps*n,
                              shape = latent_lambda+alpha,
                              rate = (1-rho)/beta)
  results <- rep(0, timesteps*n)
  which_to_use <- which(runif(n=timesteps*n) <= exp(-kappa * latent_lambda))
  results[which_to_use] <- rexp(n = length(which_to_use),
                                rate = latent_lambda[which_to_use])
  return(matrix(results, ncol=n))
}

library('evir')
wp_sim <- rwprocess(alpha = 3.00,
                      beta = 1.0,
                      kappa = 1.7,
                      rho = 0.4,
                      timesteps = 100000,
                      n=1)
plot(wp_sim[,1], type = 'l')
lines(wp_sim[,2], type = 'l', col="red")
length(which(wp_sim == 0.0))
summary(lgp_sim[,1])
gpd(wp_sim, nextremes = 44)

plot(density(wp_sim[wp_sim>0]))
lines(density(evir::rgpd(1000, xi=0.33, beta=(1+1.7)/3)))
