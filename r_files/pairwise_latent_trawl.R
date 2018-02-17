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

# Example:
rho <- 0.3
t1 <- 0.1
t2 <- 0.3
compute_A_exp(rho)
compute_B1_exp(rho, t1, t2) + compute_B_inter_exp(rho, t1, t2)
compute_B3_exp(rho, t1, t2) + compute_B_inter_exp(rho, t1, t2)

is_vector_elem <- function(vec_to_check, var_name){
  return(var_name %in% vec_to_check)
}

# Example
test_vec <- c("elem")
is_vector_elem(test_vec, "elem") # Returns True
is_vector_elem(test_vec, "ele") # Returns False

trf_inv_g <- function(z, alpha, beta, kappa, offset=10){
  # From GPD(alpha, beta) to GPD(offset, offset+kappa)
  res <- (offset+kappa)*((1+z/beta)^{alpha/offset}-1)
  return(res)
}

trf_g <- function(x, alpha, beta, kappa, offset=10){
  # From GPD(offset, offset+kappa) to GPD(xi, sigma)
  res <- beta*((1+x/(offset+kappa))^{offset/alpha}-1)
  return(res)
}

# Example
kappa <- 1.7
alpha <- 3
beta <- 1.0
rho <- 0.4

## Verification trf_inv_g and trf_g are inverse functions of eachother
plot(1:10, trf_inv_g(trf_g(x = 1:10, alpha = alpha, beta = beta, kappa = kappa), alpha = alpha, beta = beta, kappa = kappa))
plot(1:10, trf_g(trf_inv_g(z = 1:10, alpha = alpha, beta = beta, kappa = kappa), alpha = alpha, beta = beta, kappa = kappa))

## Verification of output densities
n_samples <- 2^10
gpd_data_offset <- gPdtest::rgp(n = n_samples, shape = 1/10, scale = ((10+kappa)/10))
gpd_data_ab <- gPdtest::rgp(n_samples, shape = 1/alpha, scale=beta/alpha)

### from GPD(offset, offset+kappa) to GPD(alpha, beta)
hist(trf_g(x = gpd_data_offset, alpha = alpha, beta = beta, kappa = kappa), freq = F, breaks=80)

plot(density(trf_g(x = gpd_data_offset, alpha = alpha, beta = beta, kappa = kappa)))
lines(density(gpd_data_ab), type="l", col="red")
lines(seq(0,5,length.out = 100), alpha/beta*(1+seq(0,5,length.out = 100)/beta)^{-alpha-1}, lty = 2, col="green")

#### testing the marginal fits
gPdtest::gpd.fit(gpd_data_ab, method = "amle")
gPdtest::gpd.fit(trf_g(x = gpd_data_offset, alpha = alpha, beta = beta, kappa = kappa), method ="amle")

### from GPD(alpha, beta) to GPD(offset, offset+kappa)
plot(density(gpd_data_offset), col="red", probability = T)
lines(density(trf_inv_g(z = gpd_data_ab, alpha = alpha, beta = beta, kappa = kappa)))
lines(seq(0,5,length.out = 100), 10/(10+kappa)*(1+seq(0,5,length.out = 100)/(10+kappa))^{-10-1}, lty = 2, col="green", lwd=2)

#### testing the marginal fits
gPdtest::gpd.fit(gpd_data_offset, method = "amle")
gPdtest::gpd.fit(trf_inv_g(z = gpd_data_ab, alpha = alpha, beta = beta, kappa = kappa), method ="amle")

dlgpd <- function(x, alpha, beta){
  return(alpha/beta*(1+x/beta)^{-alpha-1.0})
}

plgpd <- function(x, alpha, beta, lower.tail=F){
  res <- 1-(1+x/beta)^{-alpha}
  if(lower.tail){
    res <- 1-res
  }
  return(res)
}

trf_jacobian <- function(z, alpha, beta, kappa, offset = 10){
  # TODO check whether it is numerically stable by division of pdfs
  inv_g_z <- trf_inv_g(z = z, alpha = alpha, beta = beta, kappa = kappa)
  res <- dlgamma(x = z, alpha = alpha, beta = beta) / dlgamma(x = inv_g_z, alpha = offset, beta = offset+kappa)
  return(res)
}

# Example
n_sample <- 1000
alpha <- 6
beta <- 12
kappa <- 4
proba_trf_k <- 1/(1+kappa)
proba_trf_a_b <- (1+1/beta)^{-alpha}

proba_no_trf <- (1+kappa/beta)^{-alpha}
cat("Proba trf:", proba_trf_k)
cat("Proba no trf:", proba_no_trf)
library(fExtremes)
library(evir)

# TODO rewrite example

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
    inv_x <- trf_inv_g(x1, alpha = alpha, beta = beta, kappa = kappa)
    jacobian <- trf_jacobian(z = x1, alpha = alpha, beta = beta, kappa = kappa)
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
    inv_x1 <- trf_inv_g(x1, alpha = alpha, beta = beta, kappa = kappa)
    inv_x2 <- trf_inv_g(x2, alpha = alpha, beta = beta, kappa = kappa)
    new_x1 <- inv_x1
    new_x2 <- inv_x2
    jacobian1 <- trf_jacobian(z = x1, alpha = alpha, beta = beta, kappa = kappa)
    jacobian2 <- trf_jacobian(z = x2, alpha = alpha, beta = beta, kappa = kappa)
    temp <- jacobian1 * jacobian2
    #temp <- 1.0
    #temp <- 1/temp
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
            temp <- temp - 1000
          }
        }
      }
    }
    
    if(temp > 1e14 | abs(temp) == Inf){
      temp <- 1e10
    }
  }
  
  #cat("Accepted: ", accepted/total, "\n")
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

# pl_single_all_params_with_kappa <- function(times, values, delta, kappa, params, logscale=T, transformation=F){
#   return(pairwise_likelihood_single_full(times, values, 
#                                          alpha=1/exp(params[1]), 
#                                          beta=abs(exp(params[2])/exp(params[1])), 
#                                          kappa=exp(kappa), 
#                                          rho=exp(params[3]), 
#                                          delta=delta, 
#                                          logscale=T, 
#                                          transformation=transformation))
# }
# 
# pl_single_all_params_with_alpha_beta <- function(times, values, delta, alpha, beta, params, logscale=T, transformation=F){
#   return(pairwise_likelihood_single_full(times, values, 
#                                          alpha=alpha, 
#                                          beta=beta, 
#                                          kappa=exp(params[2]), 
#                                          rho=exp(params[1]), 
#                                          delta=delta, 
#                                          logscale=T, 
#                                          transformation=transformation))
# }
# 
# pl_single_all_params_with_alpha_beta_kappa <- function(times, values, delta, alpha, beta, kappa, params, logscale=T, transformation=F){
#   return(pairwise_likelihood_single_full(times, values, 
#                                          alpha=alpha, 
#                                          beta=beta, 
#                                          kappa=exp(kappa), 
#                                          rho=exp(params), 
#                                          delta=delta, 
#                                          logscale=T, 
#                                          transformation=transformation))
# }
# 
# pl_final <- function(times, values, deltas, params, logscale=T){
#   d <- length(values[1,])
#   mat_params <- t(matrix(params, nrow=4))
#   print(mat_params)
#   temp <- 0.0
#   
#   for(i in 1:d){
#     if(mat_params[i,1] <= 0){
#       trf <- T
#     }else{
#       trf <- F
#     }
#     temp <- temp + pl_single_all_params(times[,i], values[,i],
#                                         deltas[i], 
#                                         mat_params[i,],
#                                         transformation=trf)
#   }
#   
#   if(logscale){
#     return(temp)
#   }else{
#     return(exp(temp))
#   }
# }
# 
# pl_final_univ <- function(times, values, delta, params, trf=T, logscale=T){
#   temp <- pl_single_all_params(times = times, 
#                                values = values,
#                                delta = delta, 
#                                params = params,
#                                transformation=trf)
# 
#   if(logscale){
#     return(temp)
#   }else{
#     return(exp(temp))
#   }
# }
# 
# pl_final_univ_with_kappa <- function(times, values, delta, kappa, params, trf=T, logscale=T){
#   temp <- pl_single_all_params_with_kappa(times = times, 
#                                values = values,
#                                delta = delta,
#                                kappa = kappa,
#                                params = params,
#                                transformation=trf)
#   
#   if(logscale){
#     return(temp)
#   }else{
#     return(exp(temp))
#   }
# }
# 
# pl_final_univ_with_alpha_beta <- function(times, values, alpha, beta, delta, params, trf=T, logscale=T){
#   temp <- pl_single_all_params_with_alpha_beta(times = times, 
#                                                      values = values,
#                                                      alpha = alpha,
#                                                      beta = beta,
#                                                      delta = delta,
#                                                      params = params,
#                                                      transformation=trf)
#   
#   if(logscale){
#     return(temp)
#   }else{
#     return(exp(temp))
#   }
# }
# 
# pl_final_univ_with_alpha_beta_kappa <- function(times, values, alpha, beta, delta, kappa, params, trf=T, logscale=T){
#   temp <- pl_single_all_params_with_alpha_beta_kappa(times = times, 
#                                           values = values,
#                                           alpha = alpha,
#                                           beta = beta,
#                                           delta = delta,
#                                           kappa = kappa,
#                                           params = params,
#                                           transformation=trf)
#   
#   if(logscale){
#     return(temp)
#   }else{
#     return(exp(temp))
#   }
# }


pl_univ <- function(times, values, delta, fixed_names, fixed_params, params, model_vars_names, logscale=T, transformation=F){
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
  
  return(pl_single_all_params(times = times,
                              values = values,
                              delta = delta,
                              params = params_all,
                              logscale = logscale,
                              transformation = transformation))
}
