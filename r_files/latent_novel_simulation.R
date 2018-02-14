setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
source('pairwise_latent_trawl.R')

compute_A_exp <- function(rho){
  return(1 / rho)
}

gamma_box <- function(alpha, beta, dx, dy, n){
  alpha <- abs(alpha)
  beta <- abs(beta)
  return(rgamma(shape=alpha*dx*dy, rate=beta, n=n))
}

gamma_sqbox <- function(alpha, beta, dt, n){
  return(gamma_box(alpha = alpha, beta = beta, dx = dt, dy = dt, n = n))
}

trawl_exp <- function(t, rho, max_value=1, min_value=1e-6){
  time_eta <- -log(min_value)/rho
  return(function(s){
    if(is.na(s[1])){
      return(list(rho=rho, 
                  time_eta=time_eta, 
                  max_value=max_value, 
                  trawl_time=t,
                  A=compute_A_exp(rho)))
    }
    
    s_to_use <- which(s <= t & s >= t-time_eta)
    results <- rep(0, length(s))
    results[s_to_use] <- exp(-rho * (t-s[s_to_use]))
    return(results)
  })
}

trawl_exp_primitive <- function(t, rho, zero_at=-Inf){
  # primitive of exponential trawl function which is zero at zero_at
  return(function(s){
    if(is.na(s[1])){
      return(list(trawl_time=t))
    }
    
    s_to_use <- which(s < t)
    results <- rep(1/rho - exp(-rho * (t-zero_at)) / rho, length(s))
    results[s_to_use] <- exp(-rho * (t-s)) / rho - exp(-rho * (t-zero_at)) / rho
    return(results)
  })
}

trawl_deterministic <- function(trawl_f, alpha, beta, dt, n = n){
  trawl_info <- trawl_f(NA)
 
  time_grid <- 0:(trawl_info$time_eta/dt + 1) * dt + trawl_info$trawl_time - trawl_info$time_eta
  space_grid <- 0:(trawl_info$max_value/dt + 1) * dt
  values <- trawl_f(time_grid)
  space_ind <- rep(0, length(values))

  for(index in 1:length(values)){
    space_ind[index] <- which(space_grid >= values[index])[1]
  }

  result <- rep(0, n)
  for(index in 1:n){
    accum <- 0
    for(nb_boxes in space_ind){
      if(accum > 2^20){
        result[index] <- result[index] + sum(gamma_sqbox(alpha = alpha / trawl_info$A, beta = beta, dt = dt, n = accum))
        accum <- 0
      }
      accum <- accum + nb_boxes
    }
    if(accum > 0){
      result[index] <- result[index] + sum(gamma_sqbox(alpha = alpha, beta = beta, dt = dt, n = accum))
    }
    result[index] <- result[index]
  }
  return(result)
}

slice_area <- function(i, j, times, trawl_f_prim){
  prim_info <- trawl_f_prim(NA)
  origin_time <- prim_info$trawl_time
  times <- sort(times)
  if(i != 1 & j != length(times)){
    temp <- trawl_f_prim(times[i] - times[j] + origin_time)
    temp <- temp - trawl_f_prim(times[i] - times[j+1] + origin_time)
    temp <- temp - trawl_f_prim(times[i-1] - times[j] + origin_time)
    temp <- temp + trawl_f_prim(times[i-1] - times[j+1] + origin_time)
  }else{
      temp <- trawl_f_prim(times[i] - times[j] + origin_time)
      if(j != length(times)){
        temp <- temp - trawl_f_prim(times[i] - times[j+1] + origin_time)
      }else{
        if(i == length(times)){
          temp <- trawl_f_prim(origin_time) - trawl_f_prim(origin_time-times[i]+times[i-1])
        }
      }
  }
  
  return(temp)
}

trawl_slice_sets_not_optim <- function(alpha, beta, times, n, trawl_f, trawl_f_prim){
  times <- sort(times)
  A <- trawl_f(NA)$A
  slice_mat <- matrix(0, nrow = length(times), ncol = length(times))
  gamma_sim <- matrix(0, length(times) * length(times), ncol= n)
  for(i in 1:length(times)){
    for(j in i:length(times)){
      slice_mat[i, j] <- slice_area(i, j, times, trawl_f_prim)
      gamma_sim[(i-1) * length(times) + j,] <- rgamma(shape = alpha * slice_mat[i,j] / A,
                                                      rate = beta,
                                                      n = n)
    }
  }
  
  n_trawl_forward <- trawl_f(NA)$time_eta
  results <- matrix(0, nrow = length(times), ncol = n)
  for(current_trawl in 1:length(times)){
    for(back_trawl in max(1, current_trawl-8+1):current_trawl){
      for(slice_piece in (current_trawl):min(floor(current_trawl+n_trawl_forward+1), length(times))){
        if(back_trawl != current_trawl | slice_piece > current_trawl){
          results[current_trawl,] <- results[current_trawl,] + gamma_sim[(back_trawl-1) * length(times) + slice_piece,]
        }
      }
    }
  }
  
  for(i in 1:(length(times))){
    results[i,] <- results[i,] + gamma_sim[(i-1) * length(times) + i,]
  }
  
  #results[length(times),] <- results[length(times), ] + gamma_sim[(length(times)-1) * length(times) + length(times),]
  
  return(results)
}

trawl_slice_sets <- function(alpha, beta, times, trawl_f, trawl_f_prim, n){
  times <- sort(times)
  trawl_info <- trawl_f(NA)
  A <- trawl_info$A
  results <- matrix(0, nrow = length(times), ncol = n)
  
  for(i in 1:length(times)){
    max_index <- suppressWarnings(min(length(times), min(which(times - times[i] > trawl_info$time_eta))))
    for(j in (i+1):length(times)){
      if(j <= length(times)){
        results[i,] <- results[i,] + rgamma(shape = alpha * slice_area(i, j, times, trawl_f_prim) / A,
                                            rate = beta,
                                            n = n)
      }
    }
    if(i < length(times)){
      results[i+1,] <- results[i+1,] + results[i,]
      results[i,] <- results[i,] + rgamma(shape = alpha * slice_area(i, i, times, trawl_f_prim) / A,
                                          rate = beta,
                                          n = n)
    }
  }
  
  results[length(times),] <- results[length(times), ] + rgamma(shape = alpha * slice_area(length(times), length(times), times, trawl_f_prim) / A,
                                                               rate = beta,
                                                               n = n)
  return(results)
}

rltrawl <- function(alpha, beta, times, trawl_f = trawl_exp, trawl_f_prim=trawl_exp_primitive, n, kappa = 0, transformation=T){
  if(!transformation){
    results <- trawl_slice_sets_not_optim(alpha = alpha,
                                          beta = beta+kappa,
                                          times = times,
                                          trawl_f = trawl_f,
                                          trawl_f_prim = trawl_f_prim,
                                          n = n)
  }else{
    results <- trawl_slice_sets_not_optim(alpha = 1.0,
                                          beta = 1.0+kappa,
                                          times = times,
                                          trawl_f = trawl_f,
                                          trawl_f_prim = trawl_f_prim,
                                          n = n)
  }
  # 
  # special_g <- function(x){
  #   return(trf_g(x, xi = 1/alpha, sigma = abs(beta/alpha), kappa = kappa))
  # }
  # 
  # if(!transformation){
  #   return(results)
  # }else{
  #   return(vapply(X = results, FUN = special_g, FUN.VALUE=1))
  # }
  return(results)
}

rlexceed <- function(alpha, beta, kappa, times, trawl_f, trawl_f_prim, n, transformation){
  # Generate Gamma
  gen_trawl <- rltrawl(alpha = alpha,
                       beta = beta,
                       kappa = kappa,
                       times = times,
                       trawl_f = trawl_f,
                       trawl_f_prim = trawl_f_prim,
                       n = n,
                       transformation = transformation)

  
  unif_samples <- runif(n=length(times)*n)
  if(n == 1){
    gen_exceedances <- rep(0, length(times))
  }else{
    gen_exceedances <- matrix(0, nrow = length(times), ncol = n)
  }
  
  #print(gen_trawl)
  prob_zero <- 1-exp(-kappa * gen_trawl)
  which_zero <- which(prob_zero < unif_samples)
  print(which_zero)
  #plot(rexp(n = 1, rate = trawl_not_optim[!which_zero]))
  if(transformation){
    gen_exceedances[-which_zero] <-  vapply(rexp(n = length(gen_trawl)-length(which_zero), rate = gen_trawl[-which_zero]),
                                          FUN.VALUE = 1.0,
                                          FUN = function(x){return(trf_g(x, xi = 1/alpha, sigma = beta/alpha, kappa = kappa))})
  }else{
    gen_exceedances[-which_zero] <-  rexp(n = length(gen_trawl)-length(which_zero), rate = gen_trawl[-which_zero])
  }
  return(gen_exceedances)
}


# Example
n_sims <- 50
times <- 1:1000
kappa <- 1.7
alpha <- 3
beta <- 1
rho <- 0.4


## Trawl process simulation
trawl_1 <- trawl_exp(t, rho)
trawl_1_prim <- trawl_exp_primitive(t, rho)

### no transformation
gen_trawl <- rltrawl(alpha = alpha,
                      beta = beta,
                      times = times,
                      n = 1,
                      trawl_f = trawl_1,
                      trawl_f_prim = trawl_1_prim,
                      kappa = 2,
                      transformation = F)
acf(gen_trawl, type = "covariance")[0]
(alpha)/(beta+kappa)^2

### no transformation
(1+kappa/beta)^{-alpha}
1/(1+kappa)
par_ests_sims_no_trf <- matrix(0, ncol = 2, nrow = n_sims)
for(i in 1:n_sims){
  gen_exc <- rlexceed(alpha = alpha,
                      beta = beta,
                      kappa = kappa,
                      trawl_f = trawl_exp(length(times), rho),
                      trawl_f_prim = trawl_exp_primitive(length(times), rho),
                      times = times,
                      n = 1,
                      transformation = F)
  par_ests_sims_no_trf[i,] <- fExtremes::gpdFit(gen_exc, u=1e-6)@fit$par.ests
}

#### xi
1/alpha
boxplot(par_ests_sims_no_trf[,1])
abline(h=1/alpha, col = "red")
mean(par_ests_sims_no_trf[,1])
sd(par_ests_sims_no_trf[,1])

#### sigma
(beta+kappa)/alpha
boxplot(par_ests_sims_no_trf[,2])
abline(h=(beta+kappa)/alpha, col = "red")
mean(par_ests_sims_no_trf[,2])
sd(par_ests_sims_no_trf[,2])

gen_exc_trf <- rlexceed(alpha = alpha,
                    beta = beta,
                    kappa = kappa,
                    trawl_f = trawl_exp(length(times), rho),
                    trawl_f_prim = trawl_exp_primitive(length(times), rho),
                    times = times,
                    n = 1,
                    transformation = T)
proba_exc <- (1+kappa/beta)^{-alpha}
hist(gen_exc_trf[gen_exc_trf>0.0&gen_exc_trf<50], breaks=50, probability = T, xlim = c(0, 50))
test_data_trf <- evir::rgpd(n = 1000, xi = 1/alpha, beta =(beta)/alpha)
hist(test_data[test_data >0 & test_data < 50], breaks=50, probability = T, xlim = c(0, 50))
fExtremes::gpdFit(gen_exc_trf, type = "mle", u = quantile(gen_exc_trf, proba_exc))
