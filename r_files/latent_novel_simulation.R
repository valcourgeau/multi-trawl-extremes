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

trawl_exp <- function(t, rho, max_value=1, min_value=1e-2){
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

collection_trawl <- function(times, params, type, prim=F){
  if(!is.list(params)) stop('params should be a list.')
  # TODO Add more than exp
  if(type=="exp"){
    if(! "rho" %in% names(params)) stop('rho should in list of parameters params.')
    params_rho <- params$rho
    results <- 
    for(time_index in 1:length(times)){
      if(prim){
        return(lapply(times, function(t){trawl_exp_primitive(t, params_rho)}))
      }else{
        return(lapply(times, function(t){trawl_exp(t, params_rho)}))
      }
    }
  }
  stop('Fatal error: no trawl functions list created')
}

funlist <- collection_trawl(times = 1:5, list("rho"=0.4), type="exp")
funlist[[2]](NA)

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

trawl_slice_sets_not_optim <- function(alpha, beta, times, n, trawl_fs, trawl_fs_prim){
  # TODO sort the trawl_fs and trawl_fs_prim as the times
  # TODO current_trawl and going further back? instead of 8
  
  if(!is.list(trawl_fs)) stop('Wrong type: trawl set should be a list.')
  
  times <- sort(times)
  A <- trawl_fs[[1]](NA)$A
  slice_mat <- matrix(0, nrow = length(times), ncol = length(times))
  gamma_sim <- matrix(0, length(times) * length(times), ncol= n)
  
  # Creating the matrix of gamma realisations
  for(main_index in 1:length(times)){
    for(second_index in main_index:length(times)){
      slice_mat[main_index, second_index] <- slice_area(main_index, second_index, times, trawl_fs_prim[[main_index]])
      gamma_sim[(main_index-1) * length(times) + second_index,] <- rgamma(shape = alpha * slice_mat[main_index, second_index] / A,
                                                                          rate = beta,
                                                                          n = n)
    }
  }
    
  # Going back in time to use the dep structure
  n_trawl_forward <- trawl_fs[[1]](NA)$time_eta
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
  
  # Using independent scattering of Levy basis to add time dependence to trawls
  for(main_index in 1:(length(times))){
    results[main_index,] <- results[main_index,] + gamma_sim[(main_index-1) * length(times)+main_index,]
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


rltrawl <- function(alpha, beta, times, trawl_fs, trawl_fs_prim, n, kappa = 0, transformation=T, offset_shape=NULL, offset_scale=NULL){
  if(length(trawl_fs) != length(times)){
    stop('Wrong number of trawl functions compared to timestamps.')
  }
  if(length(trawl_fs_prim) != length(times)){
    stop('Wrong number of trawl primitives compared to timestamps.')
  }
  if(transformation & (is.null(offset_scale) | is.null(offset_shape))){
    stop('When using marginal trf, indicate shape and scale offsets.')
  }
  
  
  if(!transformation){
    results <- trawl_slice_sets_not_optim(alpha = alpha,
                                          beta = beta+kappa,
                                          times = times,
                                          trawl_fs = trawl_fs,
                                          trawl_fs_prim = trawl_fs_prim,
                                          n = n)
  }else{
    results <- trawl_slice_sets_not_optim(alpha = offset_shape,
                                          beta = offset_scale,
                                          times = times,
                                          trawl_fs = trawl_fs,
                                          trawl_fs_prim = trawl_fs_prim,
                                          n = n)
  }
  
  return(results)
}

rlexceed <- function(alpha, beta, kappa, times, trawl_fs, trawl_fs_prim, n, transformation, n_moments = 4){
  offset_shape <- n_moments+1
  offset_scale <- trf_find_offset_scale(alpha = alpha,
                                       beta = beta,
                                       kappa = kappa,
                                       offset_shape = n_moments+1)
  # Generate Gamma
  gen_trawl <- rltrawl(alpha = alpha,
                       beta = beta,
                       kappa = 0.0,
                       times = times,
                       trawl_fs = trawl_fs,
                       trawl_fs_prim = trawl_fs_prim,
                       n = n,
                       transformation = transformation,
                       offset_shape = offset_shape,
                       offset_scale = offset_scale)

  # Uniform threshold
  unif_samples <- matrix(runif(n=length(times)*n), nrow = length(times))
  if(n == 1){
    gen_exceedances <- rep(0, length(times))
  }else{
    gen_exceedances <- matrix(0, nrow = length(times), ncol = n)
  }
  
  #print(gen_trawl)
  prob_zero <- 1.0-exp(-kappa * gen_trawl)
  which_zero <- which(prob_zero >= unif_samples)
  mean(gen_exceedances)
  #plot(rexp(n = 1, rate = trawl_not_optim[!which_zero]))
  if(transformation){
    gen_exceedances[-which_zero] <-  vapply(rexp(n = length(gen_trawl)-length(which_zero), rate = gen_trawl[-which_zero]),
                                          FUN.VALUE = 1.0,
                                          FUN = function(x){return(trf_g(x, alpha = alpha, 
                                                                         beta = beta,
                                                                         kappa = kappa, 
                                                                         offset_scale = offset_scale,
                                                                         offset_shape = offset_shape))})
  }else{
    gen_exceedances[-which_zero] <-  rexp(n = length(gen_trawl[-which_zero]), rate = gen_trawl[-which_zero])
  }
  #mean(gen_exceedances)
  return(gen_exceedances)
}

# Various tests
# {
#   # Example
#   n_sims <- 50
#   times <- 1:400
#   kappa <- 0.3
#   alpha <- 3
#   beta <- 6
#   rho <- 0.3
#   n_moments <- 4
#   
#   ## Find offset scale
#   offset_shape <- n_moments + 1
#   kappa / ((1+kappa/beta)^{alpha/offset_shape} - 1)
#   trf_find_offset_scale(alpha = alpha, beta = beta, kappa = kappa, offset_shape = offset_shape)
#   offset_scale  <- trf_find_offset_scale(alpha = alpha, beta = beta, kappa = kappa, offset_shape = offset_shape)
#   
#   cat("Prob non zero for non-trf",(1+kappa/beta)^{-alpha}, "\n")
#   cat("Prob non zero for trf",(1+kappa/offset_scale)^(-offset_shape), "\n")
#   
#   ## Trawl process simulation
#   library(gPdtest)
#   
#   ### Generating the functions
#   trawl_1 <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = F)
#   trawl_1_prim <- collection_trawl(times = times, params = list(rho=rho), type = "exp", prim = T)
#   
#   trl_slice <- trawl_slice_sets_not_optim(alpha = alpha, beta = beta, times = times, trawl_fs = trawl_1, trawl_fs_prim = trawl_1_prim, n = 1)
#   
#   ### no transformation
#   gen_trawl <- rltrawl(alpha = alpha,
#                         beta = beta,
#                         times = times,
#                         n = 1,
#                         trawl_fs = trawl_1,
#                         trawl_fs_prim = trawl_1_prim,
#                         kappa = 0.0,
#                         transformation = F)
#   acf(gen_trawl, type = "covariance")
#   (alpha)/(beta)^2
#   hist(gen_trawl, probability = T)
#   lines(seq(0.01, 8, length.out = 200),dgamma(seq(0.01, 8, length.out = 200), shape = alpha, scale = 1/(beta)), col = "red")
#   
#   par(mfrow=c(1,2))
#   #### ACF
#   acf(gen_trawl, main = paste("ACF trawl with rho =", rho))
#   lines(0:20, exp(-rho*0:20), col = "red")
#   
#   #### distribution
#   plot(density(gen_trawl), "Marginal density of trawls")
#   lines(density(rgamma(n = 1000, shape = alpha, rate = beta)), col="red")
#   
#   ### no transformation
#   (1+kappa/beta)^{-alpha}
#   1/(1+kappa)
#   par_ests_sims_no_trf <- matrix(0, ncol = 2, nrow = n_sims)
#   for(i in 1:n_sims){
#     gen_exc <- rlexceed(alpha = alpha,
#                         beta = beta,
#                         kappa = kappa,
#                         trawl_fs = trawl_1,
#                         trawl_fs_prim = trawl_1_prim,
#                         times = times,
#                         n = 1,
#                         transformation = F)
#     print(mean(gen_exc))
#     par_ests_sims_no_trf[i,] <- fExtremes::gpdFit(gen_exc, u =1e-6)@fit$par.ests
#   }
#   cat("mean:", (1+kappa/beta)^{-alpha}*(beta+kappa)/(alpha-1), "\n")
#   
#   #### xi
#   1/alpha
#   boxplot(par_ests_sims_no_trf[,1])
#   abline(h=1/alpha, col = "red")
#   mean(par_ests_sims_no_trf[,1])
#   sd(par_ests_sims_no_trf[,1])
#   
#   #### sigma
#   (beta+kappa)/alpha
#   boxplot(par_ests_sims_no_trf[,2])
#   abline(h=(beta+kappa)/alpha, col = "red")
#   mean(par_ests_sims_no_trf[,2])
#   sd(par_ests_sims_no_trf[,2])
#   # OBS: depending on fitting procedure used: over or under estimation happening gPdtest::gpd.fit and fExtremes::gpdFit
#   
#   
#   ### transformation
#   par_ests_sims_trf <- matrix(0, ncol = 2, nrow = n_sims)
#   for(i in 1:n_sims){
#     gen_exc_trf <- rlexceed(alpha = alpha,
#                             beta = beta,
#                             kappa = kappa,
#                             trawl_fs = trawl_1,
#                             trawl_fs_prim = trawl_1_prim,
#                             times = times,
#                             n = 1,
#                             transformation = T)
#     par_ests_sims_trf[i,] <- fExtremes::gpdFit(gen_exc_trf, u =1e-6)@fit$par.ests
#   }
#   
#   #### xi
#   10
#   1/alpha
#   boxplot(par_ests_sims_trf[,1])
#   abline(h=1/alpha, col = "red")
#   mean(par_ests_sims_trf[,1])
#   sd(par_ests_sims_trf[,1])
#   
#   #### sigma
#   (beta+kappa)/alpha
#   (beta)/alpha
#   boxplot(par_ests_sims_trf[,2])
#   abline(h=(beta)/alpha, col = "red")
#   mean(par_ests_sims_trf[,2])
#   sd(par_ests_sims_trf[,2])
# }
