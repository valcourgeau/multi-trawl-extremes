compute_A_exp <- function(rho){
  return(1 / rho)
}

compute_B3_exp <- function(rho, t1, t2){
  # Compute A_t_2 \ A_t_1
  if(t2 > t1){
    return((1 - exp(rho * (t1-t2)))/rho)
  }else{
    if(t2 == t1){
      return(0.0)
    }else{
      return((1 - exp(rho * (t2 - t1)))/rho)
    }
  }
}

compute_B1_exp <- function(rho, t1, t2){
  # Compute A_t_1 \ A_t_2
  return(compute_B3_exp(rho, t2, t1))
}

compute_B_inter_exp <- function(rho, t1, t2){
  if(t1 > t2){
    return(exp(rho * (t2-t1))/rho)
  }else{
    if(t1 == t2){
        return(1/rho)
    }else{
      return(exp(rho * (t1-t2))/rho)
    }
  }
}

inv_g <- function(x, xi, sigma, kappa){
  sigma <- abs(sigma)
  # proba from GPD(alpha, beta)
  temp <- evir::pgpd(x, xi = xi, beta = sigma)[1]
  
  # inverve of GPD(1, 1+kappa)
  temp <- evir::qgpd(temp, xi= 1, beta = (1+kappa))[1]
  return(temp) 
}

# Example:
rho <- 0.3
t1 <- 0.1
t2 <- 0.3
compute_A_exp(rho)
compute_B1_exp(rho, t1, t2) + compute_B_inter_exp(rho, t1, t2)
compute_B3_exp(rho, t1, t2) + compute_B_inter_exp(rho, t1, t2)


# Case 0-0 

pairwise_00_1 <- function(alpha, beta, kappa){
  return(-2 * (1 + kappa / beta)^{-alpha})
}

pairwise_00_2 <- function(t1, t2, alpha, beta, kappa, rho, B1, B2, B3){
  temp = (1 + 2 * kappa / beta)^{-alpha * rho * B2}
  return((1 + kappa / beta)^{-alpha * rho * (B1 + B3)} * temp)
}

pairwise_00_exp <- function(t1, t2, alpha, beta, kappa, rho){
  B1 <- compute_B1_exp(rho, t1, t2)
  B2 <- compute_B_inter_exp(rho, t1, t2)
  B3 <- compute_B3_exp(rho, t1, t2)
  
  temp <- pairwise_00_1(alpha, beta, kappa)
  temp <- temp + pairwise_00_2(t1, t2, alpha, beta, kappa, rho, B1, B2, B3)
  return(1 + temp)
}

# Example
t1 <- 0.0
t2 <- 1.0
alpha <- -1.0
beta <- 10
rho <- 1.0
kappa <- 1.0
B1 <- compute_B1_exp(rho, t1, t2)
B2 <- compute_B_inter_exp(rho, t1, t2)
B3 <- compute_B3_exp(rho, t1, t2)

pairwise_00_exp(t1 = t1, t2 = t2, 
                alpha = alpha, beta = beta, 
                kappa = kappa, rho = rho)
answer <- 1 - 2 * (1 + 0.1) + (1 + 0.1)^{B1 + B3}*(1 + 0.2)^{B2}

alpha <- 1.0
pairwise_00_exp(t1 = t1, t2 = t2, 
                alpha = alpha, beta = beta, 
                kappa = kappa, rho = rho)
answer <- 1 - 2 * (1 + 0.1)^{-1} + (1 + 0.1)^{-B1 - B3}*(1 + 0.2)^{-B2}
answer

# Case 1-0
pairwise_10_1 <- function(t1, x1, t2, alpha, beta, kappa, rho, trawlA){
  return(alpha * rho * trawlA / beta * (1 + (kappa + x1) / beta)^{-alpha * rho * trawlA - 1})
}

pairwise_10_2_1 <- function(t1, x1, t2, alpha, beta, kappa, rho, B1){
  return( - alpha * rho / beta * (1 + (kappa + x1) / beta)^{-alpha * rho * B1 - 1})
}

pairwise_10_2_2<- function(t1, x1, t2, alpha, beta, kappa, rho, B2){
  return((1 + (2*kappa + x1) / beta)^{-alpha * rho * B2 - 1})
}

pairwise_10_2_3 <- function(t1, x1, t2, alpha, beta, kappa, rho, B3){
  return((1 + kappa / beta)^{-alpha * rho * B3})
}

pairwise_10_2_4 <- function(t1, x1, t2, alpha, beta, kappa, rho, trawlA, B1){
  return(trawlA * (1 + (kappa + x1) / beta) + B1 * kappa / beta)
}

pairwise_10_2 <- function(t1, x1, t2, alpha, beta, kappa, rho, trawlA, B1, B2, B3, transA){
  temp <- pairwise_10_2_1(t1, x1, t2, alpha, beta, kappa, rho, B1)
  temp <- temp * pairwise_10_2_2(t1, x1, t2, alpha, beta, kappa, rho, B2)
  temp <- temp * pairwise_10_2_3(t1, x1, t2, alpha, beta, kappa, rho, B3)
  temp <- temp * pairwise_10_2_4(t1, x1, t2, alpha, beta, kappa, rho, trawlA, B1)
  
  return(temp)
}

pairwise_10_exp <- function(t1, x1, t2, alpha, beta, kappa, rho, transformation=F){
  # Marginal Transformation
  if(transformation){
    xi_mt <- 1/alpha
    sigma_mt <- abs(beta/alpha)
    inv_x <- inv_g(x1, xi = xi_mt, sigma = sigma_mt, kappa = kappa)
    jacobian <- evir::dgpd(x1, xi = xi_mt, beta = sigma_mt) / evir::dgpd(inv_x, xi = 1, beta = 1+kappa)
    jacobian <- 1.0/jacobian
    new_x <- inv_x
  }else{
    new_x <- x1
  }
 
  trawlA <- compute_A_exp(rho)
  B1 <- compute_B1_exp(rho, t1, t2)
  B2 <- compute_B_inter_exp(rho, t1, t2)
  B3 <- compute_B3_exp(rho, t1, t2)
  
  temp <- pairwise_10_1(t1, x1, t2, alpha, beta, kappa, rho, trawlA)
  temp <- temp + pairwise_10_2(t1, x1, t2, alpha, beta, kappa, rho, trawlA, B1, B2, B3)
  
  # Apply MT
  if(transformation){
    temp <- temp * jacobian
  }
  
  if(temp == 0.0 || is.na(temp) || is.nan(temp)){
    temp <- 1.0
  }
    
  return(temp)
}

# Example
t1 <- 0.0
t2 <- 1.0
x1 <- 1.0
alpha <- -1.0
beta <- 10
rho <- 1.0
kappa <- 1.0
B1 <- compute_B1_exp(rho, t1, t2)
B2 <- compute_B_inter_exp(rho, t1, t2)
B3 <- compute_B3_exp(rho, t1, t2)

pairwise_10_exp(t1 = t1, t2 = t2,
                x1 = x1,
                alpha = alpha, beta = beta,
                kappa = kappa, rho = rho)
answer <- -1/10 * (1+2/10)^{0} + 1/10*(1+2/10)^{B1-1}*(1+3/10)^{B2-1}*(1+0.1)^{B3}*((B1+B2)*(1+0.2)+B1/10)
answer

alpha <- 1.0
pairwise_10_exp(t1 = t1, t2 = t2,
                x1 = x1,
                alpha = alpha, beta = beta,
                kappa = kappa, rho = rho)
answer <- 1/10 * (1+2/10)^{-B1-B2-1} - 1/10*(1+2/10)^{-B1-1}*(1+3/10)^{-B2-1}*(1+0.1)^{-B3}*((B1+B2)*(1+0.2)+B1/10)
answer

# Case 1-1

pairwise_11_1_1 <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, B1){
  return(alpha^2 * rho^2 / beta^2 * (1+(kappa+x1)/beta)^{-alpha*rho*B1-1})
}

pairwise_11_1_2 <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, B2){
  return((1+(2*kappa+x1+x2)/beta)^{-alpha*rho*B2-1})
}

pairwise_11_1_3 <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, B3){
  return((1+(kappa+x2)/beta)^{-alpha*rho*B3-1})
}

pairwise_11_1 <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, B1, B2, B3){
  temp <- pairwise_11_1_1(t1, x1, t2, x2, alpha, beta, kappa, rho, B1)
  temp <- temp * pairwise_11_1_2(t1, x1, t2, x2, alpha, beta, kappa, rho, B2)
  temp <- temp * pairwise_11_1_3(t1, x1, t2, x2, alpha, beta, kappa, rho, B3)
  return(temp)
}

pairwise_11_2_1 <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, B1, B2){
  return(B1*B2*(1+(2*kappa+x1+x2)/beta)*(1+(kappa+x2)/beta))
}

pairwise_11_2_2 <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, B1, B3){
  return(B1*B3*(1+(2*kappa+x1+x2)/beta)^2)
}

pairwise_11_2_3 <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, B2){
  temp <- B2*(B2+1/(alpha*rho))
  temp <- temp*(1+(kappa+x1)/beta)*(1+(kappa+x2)/beta)
  return(temp)
}

pairwise_11_2_4 <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, B2, B3){
  return(B2*B3*(1+(kappa+x1)/beta)*(1+(2*kappa+x1+x2)/beta))
}

pairwise_11_2 <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, B1, B2, B3){
  temp <- pairwise_11_2_1(t1, x1, t2, x2, alpha, beta, kappa, rho, B1, B2)
  temp <- temp + pairwise_11_2_2(t1, x1, t2, x2, alpha, beta, kappa, rho, B1, B3)
  temp <- temp + pairwise_11_2_3(t1, x1, t2, x2, alpha, beta, kappa, rho, B2)
  temp <- temp + pairwise_11_2_4(t1, x1, t2, x2, alpha, beta, kappa, rho, B2, B3)
  return(temp)
}

pairwise_11_exp <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, transformation=F, epsilon=1e-8){
  # Marginal Transformation
  if(transformation){
    xi_mt <- 1/alpha
    sigma_mt <- abs(beta/alpha)
    inv_x1 <- inv_g(x1, xi = xi_mt, sigma = sigma_mt, kappa = kappa)
    inv_x2 <- inv_g(x2, xi = xi_mt, sigma = sigma_mt, kappa = kappa)
    new_x1 <- inv_x1
    new_x2 <- inv_x2
    jacobian1 <- evir::dgpd(x1, xi = xi_mt, beta = sigma_mt) / (evir::dgpd(inv_x1, xi = 1, beta = 1+kappa) + epsilon)
    jacobian2 <- evir::dgpd(x2, xi = xi_mt, beta = sigma_mt) / (evir::dgpd(inv_x2, xi = 1, beta = 1+kappa) + epsilon)
    temp <- jacobian1 * jacobian2
    #temp <- 1.0
    temp <- 1/temp
  }else{
    new_x1 <- x1
    new_x2 <- x2
    temp <- 1.0
  }
  
  if(temp == 0.0 || is.na(temp) || is.nan(temp)){
    temp <- 1e-16
  }
  
  trawlA <- compute_A_exp(rho)
  B1 <- compute_B1_exp(rho, t1, t2)
  B2 <- compute_B_inter_exp(rho, t1, t2)
  B3 <- compute_B3_exp(rho, t1, t2)
  
  temp <- temp * pairwise_11_1(t1, new_x1, t2, new_x2, alpha, beta, kappa, rho, B1, B2, B3)
  temp <- temp * pairwise_11_2(t1, new_x1, t2, new_x2, alpha, beta, kappa, rho, B1, B2, B3)

  return(temp)
}

# Example
t1 <- 0.0
t2 <- 1.0
x1 <- 1.0
x2 <- 2.0
alpha <- -1.0
beta <- 10
rho <- 1.0
kappa <- 1.0
B1 <- compute_B1_exp(rho, t1, t2)
B2 <- compute_B_inter_exp(rho, t1, t2)
B3 <- compute_B3_exp(rho, t1, t2)

pairwise_11_exp(t1 = t1, x1 = x1,
                t2 = t2, x2 = x2,
                alpha = alpha, beta = beta,
                rho = rho, kappa = kappa)

answer <- 1/100*(1+0.2)^{B1-1}*(1+0.5)^{B2-1}*(1+3/10)^{B3-1}
answer <- answer * (B1*B2*(1+0.5)*(1+0.3)+B1*B3*(1+0.5)^2
  +B2*(B2-1)*(1+0.2)*(1+0.3)+B2*B3*(1+0.2)*(1+0.5))
answer

pairwise_likelihood_single_pair <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, transformation=F){
  if(x1 < 1e-16){
    if(x2 < 1e-16){
      return(pairwise_00_exp(t1, t2, alpha, beta, kappa, rho))  
    }else{
      return(pairwise_10_exp(t2, x2, t1, alpha, beta, kappa, rho, transformation))
    }
  }else{
    if(x2 < 1e-16){
      return(pairwise_10_exp(t1, x1, t2, alpha, beta, kappa, rho, transformation))
    }else{
      return(pairwise_11_exp(t1, x1, t2, x2, alpha, beta, kappa, rho, transformation))
    }
  }
}

# Example
pairwise_likelihood_single_pair(0.1, 2.0, 0.3, 5.0, 2., 3., 30, 0.3, F)
pairwise_likelihood_single_pair(0.1, 2.0, 0.3, 1.0, -2., 3., 30, 0.3, T)

pairwise_likelihood_single_full <- function(times, values, alpha, beta, kappa, rho, delta, logscale=T, transformation=F){
  ok_ind <- which(!is.na(values))
  values <- values[ok_ind]
  times <- times[ok_ind]
  k <- length(values)
  
  temp <- 0.0
  upper <- pmin(1:(k-1)+delta, k)
  lower <- 2:k
  
  accepted <- 0
  total <- 0
  for(i in 1:(k-1)){
    ind <- (lower[i]):(upper[i])
    m <- 0
    total <- total + length(ind)
    for(j in ind){
      warnon <- pairwise_likelihood_single_pair(times[i], values[i], 
                                                times[j], values[j],
                                                alpha = alpha,
                                                beta = beta,
                                                kappa = kappa,
                                                rho = rho, 
                                                transformation=transformation)
      if(!is.na(warnon) & !is.nan(warnon)){
        if(warnon > 1e-12){
          # log the result
          accepted <- accepted + 1
          temp <- temp + log(warnon)
        }else{
          if(warnon >= 0.0){
            temp <- temp - 100
          }
        }
      }
    }
    
    if(temp > 1e14 | abs(temp) == Inf){
      temp <- 1e10
    }
  }
  
  cat("Accepted: ", accepted/total, "\n")
  if(logscale){
    return(temp)
  }else{
    return(exp(temp))
  }
}

pl_single_all_params <- function(times, values, delta, params, logscale=T, transformation=F){
  return(pairwise_likelihood_single_full(times, values, 
                                         alpha = params[1], 
                                         beta = params[2], 
                                         kappa = exp(params[4]), 
                                         rho = exp(params[3]), 
                                         delta = delta, 
                                         logscale = T, 
                                         transformation = transformation))
}

pl_single_all_params_with_kappa <- function(times, values, delta, kappa, params, logscale=T, transformation=F){
  return(pairwise_likelihood_single_full(times, values, 
                                         alpha=1/exp(params[1]), 
                                         beta=abs(exp(params[2])/exp(params[1])), 
                                         kappa=exp(kappa), 
                                         rho=exp(params[3]), 
                                         delta=delta, 
                                         logscale=T, 
                                         transformation=transformation))
}

pl_single_all_params_with_alpha_beta <- function(times, values, delta, alpha, beta, params, logscale=T, transformation=F){
  return(pairwise_likelihood_single_full(times, values, 
                                         alpha=alpha, 
                                         beta=beta, 
                                         kappa=exp(params[2]), 
                                         rho=exp(params[1]), 
                                         delta=delta, 
                                         logscale=T, 
                                         transformation=transformation))
}

pl_single_all_params_with_alpha_beta_kappa <- function(times, values, delta, alpha, beta, kappa, params, logscale=T, transformation=F){
  return(pairwise_likelihood_single_full(times, values, 
                                         alpha=alpha, 
                                         beta=beta, 
                                         kappa=exp(kappa), 
                                         rho=exp(params), 
                                         delta=delta, 
                                         logscale=T, 
                                         transformation=transformation))
}

pl_final <- function(times, values, deltas, params, logscale=T){
  d <- length(values[1,])
  mat_params <- t(matrix(params, nrow=4))
  print(mat_params)
  temp <- 0.0
  
  for(i in 1:d){
    if(mat_params[i,1] <= 0){
      trf <- T
    }else{
      trf <- F
    }
    temp <- temp + pl_single_all_params(times[,i], values[,i],
                                        deltas[i], 
                                        mat_params[i,],
                                        transformation=trf)
  }
  
  if(logscale){
    return(temp)
  }else{
    return(exp(temp))
  }
}

pl_final_univ <- function(times, values, delta, params, trf=T, logscale=T){
  temp <- pl_single_all_params(times = times, 
                               values = values,
                               delta = delta, 
                               params = params,
                               transformation=trf)

  if(logscale){
    return(temp)
  }else{
    return(exp(temp))
  }
}

pl_final_univ_with_kappa <- function(times, values, delta, kappa, params, trf=T, logscale=T){
  temp <- pl_single_all_params_with_kappa(times = times, 
                               values = values,
                               delta = delta,
                               kappa = kappa,
                               params = params,
                               transformation=trf)
  
  if(logscale){
    return(temp)
  }else{
    return(exp(temp))
  }
}

pl_final_univ_with_alpha_beta <- function(times, values, alpha, beta, delta, params, trf=T, logscale=T){
  temp <- pl_single_all_params_with_alpha_beta(times = times, 
                                                     values = values,
                                                     alpha = alpha,
                                                     beta = beta,
                                                     delta = delta,
                                                     params = params,
                                                     transformation=trf)
  
  if(logscale){
    return(temp)
  }else{
    return(exp(temp))
  }
}

pl_final_univ_with_alpha_beta_kappa <- function(times, values, alpha, beta, delta, kappa, params, trf=T, logscale=T){
  temp <- pl_single_all_params_with_alpha_beta_kappa(times = times, 
                                          values = values,
                                          alpha = alpha,
                                          beta = beta,
                                          delta = delta,
                                          kappa = kappa,
                                          params = params,
                                          transformation=trf)
  
  if(logscale){
    return(temp)
  }else{
    return(exp(temp))
  }
}



