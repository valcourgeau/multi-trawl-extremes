###
# Functions to infer values for latent processes given observations
###
setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source('pairwise_latent_trawl.R')

mode.non.zero <- function(x, alpha, beta, kappa) {
  if(sum(x <= 0) > 0) stop("x should be component-wise positive.")
  if(alpha <= 0) stop("alpha should be positive.")
  if(beta <= 0) stop("beta should be positive.")
  if(kappa <= 0) stop("kappa should be positive.")

  return(alpha/(x+beta+kappa))
}

mode.zero <- function(alpha, beta, kappa){
  if(alpha <= 1) stop("alpha should be greater than 1.")
  if(beta <= 0) stop("beta should be positive.")
  if(kappa <= 0) stop("kappa should be positive.")
  
  density.zero <- function(lambda){
    res <- 1-exp(-kappa*lambda)
    res <- beta^alpha*lambda^{alpha-1}*res
    res <- res*exp(-beta*lambda)/gamma(x = alpha)
    return(-res)
  }
  
  start.value <- (alpha-1)/beta
  opt.res <- optim(density.zero, par = start.value, method = "L-BFGS-B", lower=1e-2, upper=Inf)
  if(opt.res$convergence == 0)
    return(opt.res$par)
  else stop("Numerical optimisation failed. Mode cannot be computed.")
}

get.latent.values <- function(x, alpha, beta, kappa, randomise=F){
  if(alpha < 0){
    offset_shape <- 4
    offset_scale <- trf_find_offset_scale(alpha, beta, kappa, offset_shape)
    x[x>0] <- trf_inv_g(x[x>0], alpha = alpha, beta = beta, kappa = kappa, 
                        offset_scale = offset_scale, offset_shape = offset_shape)
    alpha <- 4
    beta <- 1
  }
  res <- x %x% rep(1, 1)
  if(randomise){
    res[x > 0] <- mode.non.zero(res[x>0], alpha = alpha, beta = beta, kappa = kappa)
    res[x == 0] <- mode.zero(alpha = alpha, beta = beta, kappa = kappa)
  }else{
    res[x > 0] <- mode.non.zero(x = x[x>0], alpha = alpha, beta = beta, kappa = kappa)
    res[x == 0] <- rgamma(n = length(which(x==0)), shape = alpha, rate = beta)
  }
  
  return(res)
}

get.latent.values.mat <- function(x, val_params, randomise=F){
  # val_params contains alpha beta kappa
  if(is.vector(x))
    return(get.latent.values(x, val_params))
  
  n_vars <- length(x[1,])
  res<- vapply(1:n_vars, function(i){return(get.latent.values(x = x[,i],
                                                        alpha = val_params[i,1],
                                                        beta = val_params[i,2],
                                                        kappa = val_params[i,3],
                                                        randomise=F))}, FUN.VALUE = as.double(1:length(x[,1])))
  return(res)
}


plgpd <- function(x, p.zero, alpha, beta, kappa){
  if(p.zero < 0 | p.zero > 1) stop("p.zero should be between 0 and 1.")
  if(x == 0)
    return(runif(n = 1, min = 0, max = p.zero))
  else{
    return(p.zero + (1-p.zero)*(1-max(0, (1+sign(alpha)*x/(beta+kappa))^{-alpha})))
  }
}

plgpd.row <- function(xs, p.zeroes, params.mat){
  # params.mat contains alpha beta rho kappa
  res <- rep(0, length(xs))
  for(i in 1:length(xs)){
    res[i] <- plgpd(x = xs[i],
                    p.zero = p.zeroes[i],
                    alpha = params.mat[i,1],
                    beta = params.mat[i,2],
                    kappa = params.mat[i,4])
  }
  
  return(res)
}


plgpd.mat <- function(xs, p.zeroes, params.mat){
  # params.mat contains alpha beta rho kappa
  res <- xs %x% rep(1, 1)
  n_vars <- length(xs[1,])
  return(vapply(1:n_vars, function(i){plgpd(xs[,i],
                                            p.zero = p.zeroes[i],
                                            alpha = params.mat[i,1],
                                            beta = params.mat[i,2],
                                            kappa = params.mat[i,4])},
  FUN.VALUE = as.double(1:length(xs[,1]))))
}

platent.mat <- function(latent.vals, params.mat){
  # params.mat contains alpha beta rho kappa
  res <- latent.vals %x% rep(1, 1)
  n_vars <- length(latent.vals[1,])
  return(vapply(1:n_vars, function(i){pgamma(latent.vals[,i], shape = params.mat[i,1],
                                      rate = params.mat[j,1])
                              },
         FUN.VALUE = as.double(1:length(latent.vals[,1]))))
}

