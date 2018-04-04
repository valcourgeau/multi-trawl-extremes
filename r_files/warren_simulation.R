setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source('pairwise_latent_trawl.R')

rw_simulate <- function(alpha, beta, rho, kappa, timesteps, n){
  if(alpha <= 0) stop('alpha should be positive to simulate Warren process.')
  # TODO can only generate for n = 1
  if(rho >= 1 | rho < 0) stop("Rho should be between 0 and 1.")
  latent <- matrix(0, timesteps, 1)
  exceedances <- matrix(0, timesteps, 1)
  latent_lambda <- rgamma(n = n,
                          shape = alpha,
                          scale = beta)
  for(index in 1:timesteps){
    #print(latent_lambda)
    if(index > 1){
      latent_lambda <- rpois(n = n,
                             lambda = rho*latent_lambda*beta/(1-rho))
      latent_lambda <- rgamma(n = n,
                              shape=latent_lambda + alpha,
                              scale = (1-rho)/beta)
    }
    
    latent[index] <- latent_lambda
  }
  which_to_use <- which(runif(n=timesteps) <= exp(-kappa * latent))
  exceedances[which_to_use] <- rexp(n = length(which_to_use),
                                    rate = latent[which_to_use])
  return(exceedances)
}


rwprocess <- function(alpha, beta, rho, kappa, timesteps, n, transformation=F, n_moments=4){
  if(transformation){
    offset_shape <- n_moments + 1
    offset_scale <- trf_find_offset_scale(alpha = alpha, beta = beta, kappa = kappa, offset_shape = offset_shape)
    print(offset_shape)
    print(offset_scale)
    temp_sim <- rw_simulate(alpha = offset_shape, 
                            beta = offset_scale, 
                            rho = rho, 
                            kappa = 1.0, 
                            timesteps = timesteps, 
                            n = n)
    which_zero <- which(temp_sim < 1e-16)
    temp_sim[-which_zero] <- trf_g(x = temp_sim[-which_zero], 
                                   alpha = alpha, 
                                   beta = beta, 
                                   kappa = kappa, 
                                   offset_scale = offset_scale, 
                                   offset_shape = offset_shape)
    return(temp_sim)
  }else{
    return(rw_simulate(alpha = alpha, beta = beta, rho = rho, kappa = kappa, timesteps = timesteps, n = n))
  }
}

library('evir')

# without trf
alpha <- 3
beta <- 2
kappa <- 1.7
rho <- 0.4
wp_sim <- rwprocess(alpha = alpha,
                      beta = beta,
                      kappa = kappa,
                      rho = rho,
                      timesteps = 100000,
                      n=1)
(1+kappa/beta)^{-alpha}
length(which(wp_sim>0.0))/length(wp_sim)
fExtremes::gpdFit(wp_sim, u = 0.0)
acf(wp_sim, lag.max = 10)

plot(density(wp_sim[wp_sim > 0]), xlim = c(-0.5, 50))
hist(wp_sim[wp_sim > 0], breaks = 50, probability = T, xlim = c(0,50))
lines(density(fExtremes::rgpd(n = 10000, xi = 1/alpha, beta = (beta+kappa)/alpha)))

plot(density(wp_sim[wp_sim>0.0]))
lines(density(evir::rgpd(100000, xi=0.33, beta=(1+1.7)/3)))

# with trf
alpha <- -1
beta <- 2
kappa <- 1.7
rho <- 0.4
wp_sim_trf <- rwprocess(alpha = alpha,
                    beta = beta,
                    kappa = kappa,
                    rho = rho,
                    timesteps = 1000,
                    transformation = T,
                    n=1)
(1+kappa/beta)^{-5}
length(which(wp_sim>0.0))/length(wp_sim)
fExtremes::gpdFit(wp_sim, u = 0.0)
acf(wp_sim, lag.max = 10)

plot(density(wp_sim[wp_sim > 0]), xlim = c(-0.5, 50))
hist(wp_sim[wp_sim > 0], breaks = 50, probability = T, xlim = c(0,50))
lines(density(fExtremes::rgpd(n = 10000, xi = 1/alpha, beta = (beta+kappa)/alpha)))

plot(density(wp_sim[wp_sim>0.0]))
lines(density(evir::rgpd(100000, xi=0.33, beta=(1+1.7)/3)))
