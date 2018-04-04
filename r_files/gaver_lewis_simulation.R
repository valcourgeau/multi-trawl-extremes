# Due to Walker (2000)

rlg_simulate <- function(alpha, beta, rho, kappa, timesteps, n){
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


rlgprocess <- function(alpha, beta, rho, kappa, timesteps, n, transformation=F, n_moments=4){
  if(transformation){
    offset_shape <- n_moments + 1
    offset_scale <- trf_find_offset_scale(alpha = alpha, beta = beta, kappa = kappa, offset_shape = offset_shape)
    print(offset_shape)
    print(offset_scale)
    temp_sim <- rlg_simulate(alpha = offset_shape, 
                            beta = offset_scale, 
                            rho = rho, 
                            kappa = 1.0, 
                            timesteps = timesteps, 
                            n = n)
    which_non_zero <- which(temp_sim > 0.0)
    temp_sim[which_non_zero] <- trf_g(x = temp_sim[which_non_zero],
                                      alpha = alpha,
                                      beta = beta,
                                      kappa = kappa,
                                      offset_scale = offset_scale,
                                      offset_shape = offset_shape)
    return(temp_sim)
  }else{
    return(rlg_simulate(alpha = alpha, beta = beta, rho = rho, kappa = kappa, timesteps = timesteps, n = n))
  }
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


# with trf
alpha <- -5
beta <- 2
kappa <- 1.7
rho <- 0.4
wp_sim_trf <- rlgprocess(alpha = alpha,
                        beta = beta,
                        kappa = kappa,
                        rho = rho,
                        timesteps = 100000,
                        transformation = T,
                        n=1)
#(1+kappa/beta)^{-5}
(1+1/trf_find_offset_scale(alpha = alpha, beta = beta, kappa = kappa, offset_shape = 5))^{-5}
length(which(wp_sim_trf>0.0))/length(wp_sim_trf)
fExtremes::gpdFit(wp_sim_trf, u = 0.0)
acf(wp_sim_trf, lag.max = 10)

plot(density(wp_sim_trf[wp_sim_trf > 0]), xlim = c(-0.5, 50))
hist(wp_sim_trf[wp_sim_trf > 0], breaks = 50, probability = T, xlim = c(0,50))
lines(density(fExtremes::rgpd(n = 10000, xi = 1/alpha, beta = -(beta+kappa)/alpha)))

plot(density(wp_sim_trf[wp_sim_trf>0.0]))
lines(density(evir::rgpd(100000, xi=-0.2, beta=(2+1)/5)))
