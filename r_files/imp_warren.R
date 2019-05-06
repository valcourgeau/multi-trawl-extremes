setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source('pairwise_latent_trawl.R')
source('warren_simulation.R')

pw_lp1 <- function(x, alpha, beta){
  return((beta / (beta + x))^alpha)
}

pw_lp1_dx <- function(x, alpha, beta){
  return(-alpha / (beta+x) * pw_lp1(x, alpha, beta))
}

pw_lp2_warren <- function(x1, x2, j, alpha, beta, rho){
  temp <- 1 + (x1+x2)/beta + (1-rho^{j})*x1*x2/beta^2
  return(temp^{-alpha})
}

pw_lp2_warren_dx1 <- function(x1, x2, j, alpha, beta, rho){
  temp <- -alpha * (1/beta + (1-rho^{j})*x2/beta^2)
  return(temp*pw_lp2_warren(x1 = x1, x2 = x2, j=j, alpha = alpha+1, beta = beta, rho = rho))
}

pw_lp2_warren_dx2 <- function(x1, x2, j, alpha, beta, rho){
  temp <- -alpha * (1/beta + (1-rho^{j})*x1/beta^2)
  return(temp*pw_lp2_warren(x1 = x1, x2 = x2, j=j, alpha = alpha+1, beta = beta, rho = rho))
}

pw_lp2_warren_dx1dx2 <- function(x1, x2, j, alpha, beta, rho){
  temp <- (1 - rho^{j}) / beta^2 * pw_lp2_warren(x1 = x1, x2 = x2,
                                                 j = j, alpha = alpha+1,
                                                 beta = beta, rho = rho)
  temp <- temp + (1.0/beta + (1-rho^j)*x1/beta^2)*pw_lp2_warren_dx2(x1 = x1, x2 = x2, j = j, alpha = alpha+1, beta = beta, rho = rho)
  return( alpha * temp)
}

# TODO custom transformed scale parameter
pairwise_warren <- function(x1, x2, i, j, alpha, beta, rho, kappa, tol=1e-16, transformation=F){
  offset_shape <- 2
  if(x1 < tol & x2 < tol){
    
    temp <- pw_lp2_warren(x1 = kappa,
                      x2 = kappa,
                      j = j-i,
                      alpha = alpha,
                      beta = beta,
                      rho = rho)
    return(1.0-2*pw_lp1(kappa, alpha = alpha, beta = beta)+temp)
  }else{
    if(x2 < tol){
      if(transformation){
        offset_scale <- trf_find_offset_scale(alpha = alpha, beta = beta, kappa = kappa, offset_shape = offset_shape)
        jac <- trf_jacobian(z = x1, alpha = alpha, beta = beta, kappa = kappa, offset_shape = 5, offset_scale = offset_scale)
        inv_x <- trf_inv_g(x1, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)
        return(jac*(-pw_lp1_dx(inv_x+kappa, alpha=offset_shape, beta=offset_scale)
                    +pw_lp2_warren_dx1(x1 = inv_x+kappa,
                                       x2 = kappa,
                                       j = j-i,
                                       alpha = offset_shape,
                                       beta = offset_scale,
                                       rho = rho)
                  )
               )
      }else{
        return(-pw_lp1_dx(x1+kappa, alpha, beta)
               +pw_lp2_warren_dx1(x1 = x1+kappa,
                                  x2 = kappa,
                                  j = j-i,
                                  alpha = alpha,
                                  beta = beta,
                                  rho = rho)
               )
      }
    }else{
      if(x1 < tol){
        if(transformation){
          offset_scale <- trf_find_offset_scale(alpha = alpha, beta = beta, kappa = kappa, offset_shape = offset_shape)
          jac <- trf_jacobian(z = x2, alpha = alpha, beta = beta, kappa = kappa, offset_shape = 5, offset_scale = offset_scale)
          inv_x <- trf_inv_g(x2, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)
          return(jac*(-pw_lp1_dx(inv_x+kappa, alpha=offset_shape, beta=offset_scale)
                      +pw_lp2_warren_dx1(x1 = kappa,
                                         x2 = inv_x+kappa,
                                         j = j-i,
                                         alpha = offset_shape,
                                         beta = offset_scale,
                                         rho = rho)
                     )
          )
        }else{
          return(-pw_lp1_dx(x2+kappa, alpha, beta)
                 +pw_lp2_warren_dx1(x1 = kappa,
                                    x2 = x2+kappa,
                                    j = j-i,
                                    alpha = alpha,
                                    beta = beta,
                                    rho = rho)
          )
        }
      }else{
        # both not zero
        if(transformation){
          offset_scale <- trf_find_offset_scale(alpha = alpha, beta = beta, kappa = kappa, offset_shape = offset_shape)
          jac <- trf_jacobian(z = x1, alpha = alpha, beta = beta, kappa = kappa, offset_shape = 5, offset_scale = offset_scale)
          jac <- jac * trf_jacobian(z = x2, alpha = alpha, beta = beta, kappa = kappa, offset_shape = 5, offset_scale = offset_scale)
          inv_x1 <- trf_inv_g(x1, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)
          inv_x2 <- trf_inv_g(x2, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)
          
          return(jac*(pw_lp2_warren_dx1dx2(x1 = inv_x1+kappa,
                                                  x2 = inv_x2+kappa,
                                                  j = j-i,
                                                  alpha = offset_shape,
                                                  beta = offset_scale,
                                                  rho = rho)
                      )
          )
        }else{
          return(pw_lp2_warren_dx1dx2(x1 = x1+kappa,
                                      x2 = x2+kappa,
                                      j = j-i,
                                      alpha = alpha,
                                      beta = beta,
                                      rho = rho))
        }
      }
    }
  }
}

pl_warren_full <- function(times, values, delta, params, logscale=T, transformation=F){
  pl <- 0.0
  n <- length(times)
  
  # Decode params
  alpha <- params[1]
  beta <- params[2]
  rho <- params[3]
  kappa <- params[4]
  i <- 0
  for(first_index in 1:(n-1)){
    for(second_index in (first_index+1):min(n, first_index+delta)){
      #cat(data[first_index], data[second_index], "\n")
      temp <- pairwise_warren(x1 = values[first_index],
                          x2 = values[second_index],
                          i = first_index,
                          j = second_index,
                          alpha = alpha,
                          beta = beta,
                          rho = rho,
                          kappa = kappa)
      #print(temp)
      if(is.na(temp) | is.nan(temp)){
        temp <- 0
      }
      if(temp > 0){
      # print(values[first_index])
      # print(values[second_index])
      # print(temp)
        pl <- pl + log(temp)
      }else{
        #i <- i + 1
        #cat(data[first_index], data[second_index], temp, "\n")
        #pl <- pl - 1000
      }
    }
  }
  #cat("Not accepted:", i/length(data), "\n")
  return(pl)
}

pl_univ_warren <- function(times, values, delta, fixed_names, fixed_params, params, model_vars_names, logscale=T, transformation=F){
  if(length(fixed_names) > length(model_vars_names)) stop('Too many fixed parameters compared to number of model params.')
  if(length(fixed_params) + length(params) != length(model_vars_names)) stop('Wrong number of params compared to model specs.')
  if(length(fixed_params) != length(fixed_names)) stop('fixed_params and fixed_names should have same length.')
  
  opti_params <- !(model_vars_names %in% fixed_names)
  opti_params_names <- model_vars_names[opti_params]
  params_all <- rep(0, length(model_vars_names))
  params_all[opti_params] <- params
  
  if(length(fixed_params) > 0){
    params_all[!opti_params] <- fixed_params
  }
  
  return(pl_warren_full(times = times,
                        values = values,
                        delta = delta,
                        params = params_all,
                        logscale = logscale,
                        transformation = transformation))
}

library(evir)
lw_sim <- rwprocess(alpha = 3.00,
                      beta = 1.0,
                      kappa = 1.7,
                      rho = 0.4,
                      timesteps = 5000,
                      n=1)
gpd(lw_sim, nextremes = length(which(lw_sim > 1e-8)))$par.est

# 
# plot(density(rgpd(n = 10000, xi = 1/169, beta = (2168+27)/169)))
# lines(density(lgp_sim[lgp_sim > 0.0]))
params_init <- c(3.00, 1.0, 0.4, 1.7) * 0.7
k <- pl_warren_full(lw_sim, params_init, 4)
-pl_warren_full(lw_sim, c(3.00, 1.0, 0.4, 1.7), 4)
lambda <- 100

fn_to_optim <- function(params){
  # epsilon_1 <- rep(0, length(params))
  # epsilon_2 <- rep(0, length(params))
  # H_n <- matrix(0, length(params), length(params))
  # loglik <- pl_warren_full(lw_sim, params, 4)
  # 
  # for(eps_ind1 in 1:length(params)){
  #   epsilon_1[eps_ind1] <- 1e-6
  #   for(eps_ind2 in 1:length(params)){
  #     epsilon_2[eps_ind2] <- 1e-6
  #     
  #     loglik_plus_plus <- pl_warren_full(lw_sim, params+epsilon_1+epsilon_2, 4)
  #     loglik_plus_minus <- pl_warren_full(lw_sim, params+epsilon_1-epsilon_2, 4)
  #     loglik_minus_plus <- pl_warren_full(lw_sim, params-epsilon_1+epsilon_2, 4)
  #     loglik_minus_minus <- pl_warren_full(lw_sim, params-epsilon_1-epsilon_2, 4)
  #     
  #     H_n[eps_ind1, eps_ind2] <- (loglik_plus_plus - loglik_plus_minus - loglik_minus_plus + loglik_minus_minus) / (4*1e-12)
  #     epsilon_2[eps_ind2] <- 0.0
  #   }
  #   epsilon[eps_index] <- 0.0
  # }
  # loglik <- pl_warren_full(lw_sim, params, 4)
  # loglik_plus <- pl_warren_full(lw_sim, params+epsilon, 4)
  # loglik_minus <- pl_warren_full(lw_sim, params-epsilon, 4)
  # H_n <- loglik_plus - 2*loglik + loglik_minus
  temp<- -pl_warren_full(lw_sim, params, 4) +lambda*sum(abs(1/params) * abs(sum((params-params_init)^2)))
  
  return(temp)
}
fn_to_optim(c(3.00, 1.0, 0.4, 1.7))
fn_to_optim(c(1.65, 0.96, 1.00, 1.51))
fn_to_optim(c(1.65, 0.1, 1.00, 1.51))

optim(fn = fn_to_optim, par=params_init, method="L-BFGS-B", control=list(trace=2, 
                                                                         pgtol = 1e-6),
      lower = rep(0.01, 4), upper = c(5, 5, 0.99, 5))
