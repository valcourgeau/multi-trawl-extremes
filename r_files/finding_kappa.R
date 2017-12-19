setwd("C:/Users/Valentin/Documents/MRes/Data")
source('pairwise_latent_trawl.R')

# exp used to have positive kappa yet global optimisation can be used
likelihood_find_kappa <- function(z, xi, sigma){
  return(function(params){
    x <- inv_g(x = z, xi = xi, sigma = sigma, kappa = exp(params))
    
    x <- dgpd(x = x, xi = 1, beta = 1+exp(params))
    return(-sum(log(x[which(x>0.0)]), na.rm = T))
  })
}
