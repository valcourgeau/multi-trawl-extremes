# Due to Walker (2000)

rlgprocess <- function(alpha, beta, rho, kappa, timesteps, n){
  # works only for n = 1
  if(rho >= 1 | rho <= 0) stop("Rho should be between 0 and 1.")
  latent <- matrix(0, timesteps, n)
  exceedances <- matrix(0, timesteps, n)
  
  latent_lambda <- rgamma(n = n, 
                          shape = alpha,
                          scale = beta)
  
  for(index in 1:timesteps){
    if(index > 1){
      latent_pt <- rgamma(n = n,
                          shape = alpha,
                          scale = 1.0)
      latent_pt <- rpois(n = n,
                             lambda = latent_pt*(1-rho)/rho)
      latent_pt <- rgamma(n = n, 
                              shape = latent_pt,
                              rate = beta/rho)
      latent_lambda <- rho*latent_lambda + latent_pt 
    }
    
    latent[index] <- latent_lambda
  }
  
  which_to_use <- which(runif(n=timesteps) <= exp(-kappa * latent))
  exceedances[which_to_use] <- rexp(n = length(which_to_use),
                                rate = latent[which_to_use])
  return(exceedances)
}
setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
library(evir)
lgp_sim <- rlgprocess(alpha = 3.00,
                      beta = 1.0,
                      kappa = 1.7,
                      rho = 0.4,
                      timesteps = 100000,
                      n=1)
length(which(lgp_sim>0.0)/length(lgp_sim))
plot(lgp_sim[,1], type = 'l')
lines(lgp_sim[,2], type = 'l', col="red")
summary(lgp_sim[,1])
gpd(lgp_sim, threshold = 1e-10)$par.est
c(1/3,  (1+1.7)/3)

plot(density(lgp_sim[lgp_sim>0]))
lines(density(evir::rgpd(10000, xi=0.33, beta=(1+1.7)/3)), col = "red")

acf(lgp_sim)

