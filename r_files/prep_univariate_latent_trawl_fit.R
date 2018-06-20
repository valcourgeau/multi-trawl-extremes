
require(hypergeo)
zeta <- function(alpha, beta, kappa){
  res.zeta <- beta^alpha / ((beta+kappa)^{alpha-1})
  res.zeta <- res.zeta * hypergeo(A = alpha-1, B = alpha-1, C = alpha, z = - beta/(beta+kappa))
  return(Re(res.zeta))
}

get_t1 <- function(alpha, beta, kappa, rho){
  # TODO adapt to other trawl functions
  d_plus <- beta^2 / ((alpha-2)*(alpha-1))*(1+kappa/beta)^{2-alpha}
  d_plus <- d_plus * (log(1+2*kappa/beta)^{1-alpha} + (2*alpha-3)/((alpha-2)*(alpha-1)))
  d_times <- beta * log(1+kappa/beta)/(alpha*(alpha-1)-1)*(alpha*beta/(alpha-2)*log(1+kappa/beta)*(1+2*kappa/beta)^{2-alpha} + alpha/(alpha-2)*zeta(alpha-1,beta,kappa) + zeta(alpha+1,beta,kappa))
  return(-rho * (d_plus - d_times))
}

get_estimate_rho <- function(alpha, beta, kappa, index, data){
  if(alpha < 0){
    data[data>0] <- vapply(data[data>0], function(x){trf_inv_g(x1, alpha = alpha, 
                                                           beta = beta, kappa = kappa, 
                                                           offset_scale = 2, 
                                                           offset_shape = 1+kappa)}, FUN.VALUE = 1.0)
    alpha <- 4
    beta <- 1
  }
 
  
  d_plus <- beta^2 / ((alpha-2)*(alpha-1))*(1+2*kappa/beta)^{2-alpha}
  d_plus <- d_plus * (log(1+2*kappa/beta) + (2*alpha-3)/((alpha-2)*(alpha-1)))
  #d_times <- beta * log(1+kappa/beta)/(alpha*(alpha-1)-1)*(alpha*beta/(alpha-2)*log(1+kappa/beta)*(1+2*kappa/beta)^{2-alpha} + alpha/(alpha-2)*zeta(alpha-1,beta,kappa) + zeta(alpha+1,beta,kappa))
  
  d_times <- 2*beta/((alpha-2)*(alpha-1))*(beta*(1+2*kappa/beta)^{2-alpha}*log(1+kappa/beta)+zeta(alpha = alpha-1, beta = beta, kappa = kappa))
  return( abs(var(data) / (index * alpha * (d_plus-d_times))))
}

generate_parameters <- function(data, cluster.size){
  params_to_work_with <- rep(0, 4)
  fit_marginal <-  fExtremes::gpdFit(data[data > 0], u= 0)@fit$fit$par
  p_nz <- length(which(data > 0))/length(data)
  
  params_to_work_with <- rep(0, 4)
  params_to_work_with[1] <- 1/fit_marginal[1]
  params_to_work_with[2] <- abs(fit_marginal[2]*params_to_work_with[1])
  
  params_to_work_with[4] <- abs(params_to_work_with[2]*(p_nz^{-1/params_to_work_with[1]}-1) / (p_nz^{-1/params_to_work_with[1]}))
  params_to_work_with[2] <- params_to_work_with[2] - params_to_work_with[4]
  params_to_work_with[4] <- (params_to_work_with[4])
  
  # if(params_to_work_with[1] > 4.5){
  #   al <- params_to_work_with[1]
  #   ratio <- (params_to_work_with[2]+params_to_work_with[4])/params_to_work_with[1]
  #   params_to_work_with[1] <- 4 
  #   params_to_work_with[2] <- ratio * 4 - params_to_work_with[4]
  #   params_to_work_with[4] <- (ratio - params_to_work_with[2]/al)/4
  # }
  # beta_trf <- params_to_work_with[4]/(p_nz^{-1/4}-1)
  # params_to_work_with[3] <- get_estimate_rho(alpha = 4, beta = beta_trf, kappa = params_to_work_with[4], cluster.size, 
  #                                            trf_inv_g(z = latent_1[,i], alpha = params_to_work_with[1], beta = params_to_work_with[2],
  #                                                      kappa = params_to_work_with[4], offset_scale = beta_trf+params_to_work_with[4], offset_shape = 4))
  # params_to_work_with[3] <- min(params_to_work_with[3], 1.0)
  
  params_to_work_with[3] <- get_estimate_rho(alpha = params_to_work_with[1],
                                             beta =  params_to_work_with[2],
                                             kappa = params_to_work_with[4],
                                             cluster.size, 
                                             data)
  #params_to_work_with[3] <- (min(params_to_work_with[3], 1.0))
  
  if(is.na(params_to_work_with[3])){
    print("NA rho")
    params_to_work_with[3] <- runif(min = 0.1,0.5, n = 1)
  }
  return(params_to_work_with)
}