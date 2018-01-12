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
    
    s_to_use <- which(s < t & s > t-time_eta)
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

trawl_slice_sets <- function(alpha, beta, times, trawl_f, trawl_f_prim, n){
  times <- sort(times)
  A <- trawl_f(NA)$A
  slice_mat <- matrix(0, nrow = length(times), ncol = length(times))
  for(i in 1:length(times)){
    for(j in i:length(times)){
      slice_mat[i,j] <- slice_area(i, j, times, trawl_f_prim)
    }
  }
  
  results <- matrix(0, nrow = length(times), ncol = n)
  for(i in 1:length(times)){
    for(j in (i+1):length(times)){
      if(j <= length(times)){
        results[i,] <- results[i,] + rgamma(shape = alpha * slice_mat[i,j] / A,
                                            rate = beta,
                                            n = n)
      }
    }
  }
  
  for(i in 1:(length(times)-1)){
    results[i+1,] <- results[i+1,] + results[i,]
    results[i,] <- results[i,] + rgamma(shape = alpha * slice_mat[i,i] / A,
                                        rate = beta,
                                        n = n)
  }
  
  results[length(times),] <- results[length(times), ] + rgamma(shape = alpha * slice_mat[length(times),length(times)] / A,
                                                               rate = beta,
                                                               n = n)
  
  return(results)
}

rltrawl <- function(alpha, beta, times, trawl_f, trawl_f_prim, n, kappa = 0, transformation=T){
  results <- trawl_slice_sets(alpha = alpha,
                              beta = beta,
                              times = times,
                              trawl_f = trawl_f,
                              trawl_f_prim = trawl_f_prim,
                              n = n)
  if(!transformation){
    return(results)
  }else{
    res_trf <- replicate(1, results)
    for(i in 1:length(times)){
      for(j in 1:n){
        res[i, j] <- inv_g(x = res[i, j], xi = 1/alpha, sigma = abs(beta/alpha), kappa = kappa)
      }
    }
  }
}

# Example
alpha <- 0.7
beta <- 8
t <- 2
rho <- 3
dtime <- 0.001
n <- 10
trawl_1 <- trawl_exp(t, rho)
trawl_1_prim <- trawl_exp_primitive(t, rho)
plot(-100:100/10, trawl_1(-100:100/10))
l <- trawl_deterministic(trawl_1, alpha = alpha, beta = beta, dt = dtime, n=1)
hist(l)
hist(rgamma(shape = alpha, rate = beta, n=50))


times <- c(1,2,3,4,5,6)
slsl <- trawl_slice_sets(alpha = alpha,
                 beta = beta,
                 times = times,
                 trawl_1,
                 trawl_1_prim,
                 1000)
hist(slsl[1,])
hist(rgamma(shape = alpha, rate = beta, 1000))
