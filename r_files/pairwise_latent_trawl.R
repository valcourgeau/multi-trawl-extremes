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

trf_inv_g <- function(z, alpha, beta, kappa, offset_scale, offset_shape){
  # From GPD(alpha, beta+kappa) to GPD(offset, offset)
  res <- (offset_scale)*((1+sign(alpha)*z/(beta+kappa))^{alpha/offset_shape}-1)

  return(res)
}

trf_g <- function(x, alpha, beta, kappa, offset_scale, offset_shape){
  # From GPD(offset, offset) to GPD(alpha, beta+kappa)
  res <- sign(alpha)*(beta+kappa)*((1+x/(offset_scale))^{offset_shape/alpha}-1)
  return(res)
}

trf_find_offset_scale <- function(alpha, beta, kappa, offset_shape){
  #extract_inverse_shape <- (1+kappa/beta)^{alpha/offset_shape} - 1
  #return(1/extract_inverse_shape - 1.0)
  
  # We conserve entropy
  # val_scale <- offset_shape * abs((beta+kappa)/alpha) * exp(1/alpha - 1/offset_shape)
  # return(val_scale - 1)
  #return(kappa/((1+sign(alpha)*kappa/beta)^{alpha/offset_shape}-1))
  return(1+kappa)
  #return(1+kappa/beta)
}

# Example
kappa <- 1.2
alpha <- 3
beta <- 3
rho <- 0.4
n_moments <- 4

offset_shape <- n_moments + 1
offset_scale <- trf_find_offset_scale(alpha = alpha, beta = beta, kappa = kappa, offset_shape = offset_shape)

## Verification of output densities
n_samples <- 2^16
gpd_data_offset <- gPdtest::rgp(n = n_samples, shape = 1/offset_shape, scale = (offset_scale)/offset_shape)
gpd_data_ab <- gPdtest::rgp(n_samples, shape = 1/alpha, scale=(beta+kappa)/alpha)

## Verification trf_inv_g and trf_g are inverse functions of eachother
plot(1:5, 1:5-trf_inv_g(trf_g(x = 1:5, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape), 
                     alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape), ylab = "id")
plot(1:5, 1:5-trf_g(trf_inv_g(z = 1:5, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape), 
                 alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape), ylab = "id")

#### testing the marginal fits
gPdtest::gpd.fit(gpd_data_ab, method = "amle") # should be 1/3 (1.2+3)/3
gPdtest::gpd.fit(trf_g(x = gpd_data_offset, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape), method ="amle")


# Example 2: negative alpha
kappa <- 1.2
alpha <- -3
beta <- 0.5
rho <- 0.4
n_moments <- 1

offset_shape <- n_moments + 1
offset_scale <- trf_find_offset_scale(alpha = alpha, beta = beta, kappa = kappa, offset_shape = offset_shape)

## Verification trf_inv_g and trf_g are inverse functions of eachother
plot(1:5, 1:5-trf_inv_g(trf_g(x = 1:5, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape), 
                        alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape), ylab = "id")
plot(1:5, 1:5/5-trf_g(trf_inv_g(z = 1:5/5, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape), 
                    alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape), ylab = "id")


## Verification of output densities
n_samples <- 2^16
gpd_data_offset <- gPdtest::rgp(n = n_samples, shape = 1/offset_shape, scale = (offset_scale)/offset_shape)
gpd_data_ab <- gPdtest::rgp(n_samples, shape = 1/alpha, scale=-(beta+kappa)/alpha)

### from GPD(offset, offset+kappa) to GPD(alpha, beta)
hist(trf_g(x = gpd_data_offset, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape), freq = F, breaks=80)

plot(density(trf_g(x = gpd_data_offset, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)))
lines(density(gpd_data_ab), type="l", col="red")
lines(seq(0,15,length.out = 100), -alpha/(beta+kappa)*(1-seq(0,15,length.out = 100)/(beta+kappa))^{-alpha-1}, lty = 2, col="green")

#### testing the marginal fits
gPdtest::gpd.fit(gpd_data_ab, method = "combined") # -1/3 (0.5+1.2)/3
gPdtest::gpd.fit(trf_g(x = gpd_data_offset, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape), method ="combined")

### from GPD(alpha, beta) to GPD(offset, offset+kappa)
hist(gpd_data_offset[gpd_data_offset < 10], col="red", xlim=c(0,10), breaks=10, probability = T)
lines(density(trf_inv_g(z = gpd_data_ab[gpd_data_ab < 1], alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)))
lines(seq(0,5,length.out = 100), offset_shape/(offset_scale)*(1+seq(0,5,length.out = 100)/(offset_scale))^{-offset_shape-1}, lty = 2, col="green", lwd=2)

#### testing the marginal fits
gPdtest::gpd.fit(gpd_data_offset, method = "amle") # 1/2 
gPdtest::gpd.fit(trf_inv_g(z = gpd_data_ab, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape), method ="amle")

dlgpd <- function(x, alpha, beta){
  return(abs(alpha)/beta*max(0,(1+sign(alpha)*x/beta))^{-alpha-1.0})
}

plgpd <- function(x, alpha, beta, lower.tail=F){
  res <- 1-(1+x/beta)^{-alpha}
  if(lower.tail){
    res <- 1-res
  }
  return(res)
}

trf_jacobian <- function(z, alpha, beta, kappa, offset_scale, offset_shape){
  # TODO check whether it is numerically stable by division of pdfs
  inv_g_z <- trf_inv_g(z = z, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)
  res <- dlgpd(x = z, alpha = alpha, beta = beta+kappa) / dlgpd(x = inv_g_z, alpha = offset_shape, beta = offset_scale)
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

pairwise_10_exp <- function(t1, x1, t2, alpha, beta, kappa, rho, transformation=F, n_moments=4){
  # Marginal Transformation
  if(transformation){
    offset_shape <- n_moments + 1
    offset_scale <- trf_find_offset_scale(alpha = alpha, beta = beta, kappa = kappa, offset_shape = offset_shape)
    inv_x <- trf_inv_g(x1, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)
    jacobian <- trf_jacobian(z = x1, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)
    new_x <- inv_x
  }else{
    new_x <- x1
  }
 
  trawlA <- compute_A_exp(rho)
  B1 <- compute_B1_exp(rho, t1, t2)
  B2 <- compute_B_inter_exp(rho, t1, t2)
  B3 <- compute_B3_exp(rho, t1, t2)
  
  
  if(transformation){
    temp <- pairwise_10_1(t1, new_x, t2, alpha=offset_shape, beta=1, kappa, rho, trawlA)
    temp <- temp + pairwise_10_2(t1, new_x, t2, alpha=offset_shape, beta=1, kappa, rho, trawlA, B1, B2, B3)
    temp <- temp * jacobian
  }else{
    temp <- pairwise_10_1(t1, new_x, t2, alpha, beta, kappa, rho, trawlA)
    temp <- temp + pairwise_10_2(t1, new_x, t2, alpha, beta, kappa, rho, trawlA, B1, B2, B3)
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

pairwise_11_exp <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, transformation=F, epsilon=1e-8, n_moments=4){
  # Marginal Transformation
  if(transformation){
    offset_shape <- n_moments + 1
    offset_scale <- trf_find_offset_scale(alpha = alpha, beta = beta, kappa = kappa, offset_shape = offset_shape)
    inv_x1 <- trf_inv_g(x1, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)
    inv_x2 <- trf_inv_g(x2, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)
    new_x1 <- inv_x1
    new_x2 <- inv_x2
    jacobian1 <- trf_jacobian(z = x1, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)
    jacobian2 <- trf_jacobian(z = x2, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)
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
  
  if(transformation){
    temp <- temp * pairwise_11_1(t1, new_x1, t2, new_x2, alpha = offset_shape, beta = 1, kappa, rho, B1, B2, B3)
    temp <- temp * pairwise_11_2(t1, new_x1, t2, new_x2, alpha = offset_shape, beta = 1, kappa, rho, B1, B2, B3)
  }else{
    temp <- temp * pairwise_11_1(t1, new_x1, t2, new_x2, alpha, beta, kappa, rho, B1, B2, B3)
    temp <- temp * pairwise_11_2(t1, new_x1, t2, new_x2, alpha, beta, kappa, rho, B1, B2, B3)
  }
  

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
  # TODO check whether t1 should be <= t2 or not
  #print(x1)
  if(x1 < 1e-16){
    if(x2 < 1e-16){
      return(pairwise_00_exp(t1, t2, alpha = 4, beta = 1, kappa, rho))  
    }else{
      return(pairwise_10_exp(t2, x2, t1, alpha, beta, kappa, rho, transformation, n_moments = 3))
    }
  }else{
    if(x2 < 1e-16){
      return(pairwise_10_exp(t1, x1, t2, alpha, beta, kappa, rho, transformation, n_moments = 3))
    }else{
      return(pairwise_11_exp(t1, x1, t2, x2, alpha, beta, kappa, rho, transformation, n_moments = 3))
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
          #cat("x1 ",values[i],"x2",values[j],"\n")
          #print(warnon)
          temp <- temp + log(warnon)
        }else{
          if(warnon >= 0.0){
            temp <- temp + log(warnon)
          }
        }
      }
    }
    
    # if(temp > 1e14 | abs(temp) == Inf){
    #   temp <- 1e10
    # }
  }
  
  #cat("Accepted: ", accepted/total, "\n")
  if(logscale){
    return(temp)
  }else{
    return(exp(temp))
  }
}

pl_single_all_params <- function(times, values, delta, params, logscale=T, transformation=F){
  # TODO add general model parameter names
  #print(params)
  return(pairwise_likelihood_single_full(times, values, 
                                         alpha = (params["alpha"][[1]]), 
                                         beta = (params["beta"][[1]]), 
                                         kappa = params["kappa"][[1]], 
                                         rho = params["rho"][[1]], 
                                         delta = delta, 
                                         logscale = T, 
                                         transformation = transformation))
}

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
  
  params_list <- list()
  params_list[fixed_names] <- fixed_params
  #params_list[opti_params_names] <- params
  for(i in 1:length(fixed_names)){
    params_list[fixed_names[i]] <- fixed_params[i]
  }
  
  params_list[opti_params_names] <- params
  
  return(pl_single_all_params(times = times,
                              values = values,
                              delta = delta,
                              params = params_list,
                              logscale = logscale,
                              transformation = transformation))
}

marginal_gpd_likelihood <- function(values, fixed_names, fixed_params, params, model_vars_names, logscale=T, transformation=F, n_moments=4){
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
  
  if(transformation){
    if(length(params_all) < 3) stop('Marginal GPD with transformation requires 3 parameters: alpha, beta and kappa.')
    print(params_all)
    lik <- vapply(values, 
                  function(x){
                    temp_alpha <- params_all[1]
                    temp_beta <- params_all[2]
                    temp_kappa <- params_all[3]
                    offset_shape <- n_moments + 1
                    offset_scale <- trf_find_offset_scale(alpha = temp_alpha, beta = temp_beta, 
                                                          kappa = temp_kappa, offset_shape = offset_shape)
                    inv_x1 <- trf_inv_g(x1, alpha = temp_alpha, beta = temp_beta, kappa = temp_kappa, 
                                        offset_scale = offset_scale, offset_shape = offset_shape)
                    jacobian1 <- trf_jacobian(z = x1, alpha = temp_alpha, beta = temp_beta, kappa = temp_kappa, 
                                              offset_scale = offset_scale, offset_shape = offset_shape)
                    temp <- jacobian1
                      return(dlgpd(x = inv_x1, alpha = offset_shape, beta = offset_scale+temp_kappa)*jacobian1)
                    },
                  1.0)
  }else{
    lik <- vapply(values, function(x){return(dlgpd(x = x, alpha = params_all[1], beta = (params_all[2]+params_all[3])))}, 1.0)
  }
  
  if(logscale){
    return(sum(log(lik[lik > 0.0])))
  }else{
    return(prod(lik))
  }
}

marginal_simple_lik <- function(values, params){
  alpha <- params[1]
  beta_kappa <- params[2]
  lik <- alpha / beta_kappa * (1+sign(alpha) * values / beta_kappa)^{-alpha-1}
  return(sum(log(lik)))
}

gpd_fit <- function(values, initial_guess){
  fn_to_optim <- function(x){return(-marginal_simple_lik(values = values, params = x))}
  ff <- optim(fn_to_optim, par = initial_guess, method = "L-BFGS-B", lower = c(0.1,0.1), upper = c(20,20))
  return(ff$par)
}

mom_gpd <- function(values_array){
  # workds under the assumption that alpha > 2
  # values_array contains the time series with first axis as time and second as # of time series
  n_dim <- length(values_array[1,])
  n_values <- length(values_array[,1])
  alphas_mom <- rep(0, n_dim)
  betas_mom <- rep(0, n_dim)
  kappas_mom <- rep(0, n_dim)
  
  for(index in 1:n_dim){
    var_mom <- var(values_array[,index][values_array[,index]>0])
    mean_mom <- mean(values_array[,index][values_array[,index]>0])
    p_mom <- length(values_array[,index][values_array[,index]>0])/n_values
    
    alphas_mom[index] <- 2*var_mom/(var_mom-mean_mom^2)
    betas_mom[index] <- mean_mom*(alphas_mom[index]-1)
    
    kappas_mom[index] <- betas_mom[index] * (1.0 - p_mom^{1/alphas_mom[index]})
    betas_mom[index] <- betas_mom[index] - kappas_mom[index]
  }
  
  return(list(alpha=alphas_mom, beta=betas_mom, kappa=kappas_mom,
              mean_alpha=mean(alphas_mom), mean_beta=mean(betas_mom), mean_kappa=mean(kappas_mom),
              sd_alpha=sd(alphas_mom), sd_beta=sd(betas_mom), sd_kappa=sd(kappas_mom)))
}

