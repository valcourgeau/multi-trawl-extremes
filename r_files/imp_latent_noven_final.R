  setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data")
  col_names <- c("Date", "Time",
                 "O3", "Status O3",
                 "NO", "Status NO",
                 "NO2", "Status NO2",
                 "NO2x", "Status NO2x",
                 "SO2", "Status SO2",
                 "SO2x", "Status SO2x",
                 "CO", "Status CO",
                 "PM10", "Status PM10",
                 "PM2.5", "Status PM2.5")
  col_classes <- c("Date", "Date", 
                   "numeric", "character", 
                   "numeric", "character", 
                   "numeric", "character", 
                   "numeric", "character", 
                   "numeric", "character", 
                   "numeric", "character", 
                   "numeric", "character",
                   "numeric", "character")
  
  bl <- read.csv("bloomsbury_1994_1998.csv", sep = ",",
                 as.is = T)
  bl <- bl[1:2000,]
  univ_model_vars_names <- c("alpha", "beta", "rho", "kappa")
  setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
  # list of packages
  library(spatstat)
  library(fExtremes)
  library(gPdtest)
  library(DEoptim)
  library(evir)
  library('optimr')
  library(matlib)
  
  # source files
  source('pairwise_latent_trawl.R')
  source('differentiate.R')
  
  # Function
  {
    initial_guess_trawl <- function(values_array){
        mom_estimates <- mom_gpd(values_array)
        d <- length(mom_estimates$alpha)
        rho <- rep(0, d)
        for(index in 1:d){
          acf_values <- acf(values_array[,index], type="covariance")$acf
          acf_values <- acf_values[2:(which.min(acf_values > 0.0)-1)]
          rho[index] <- -line(log(acf_values))$coefficients[2]
        }
        mom_estimates$rho <- rho
        mom_estimates$mean_rho <- mean(rho)
        mom_estimates$sd_rho <- sd(rho)
        return(mom_estimates)
    }
    
    create_ig_trawl <- function(dataset, quantiles, rho, kappa, trf=F)
    {
      d <- length(dataset[1,]) - 1
      params <- rep(0,4*d)
      
      for(i in 1:d){
        x_to_work_on <- dataset[,i]
        fit <- evir::gpd(x_to_work_on, method = "ml", threshold = quantiles[i])$par.ests
        fit[1] <- 1/fit[1]
        fit[2] <- fit[2]*fit[1]
        params[(4*(i-1)+1):(4*(i-1)+2)] <- fit
        
        p <- length(which(dataset[,i] > quantiles[i])) / length(dataset[,i])
        if(trf){
          temp_kappa <- (1-p)/p
        }else{
          temp_kappa <- fit[2]*(p^{-1.0/fit[1]} - 1.0)
          params[4*(i-1)+2] <- params[4*(i-1)+2] - temp_kappa
        }
        
        rho <- - line(x = 1:15, log(acf(dataset[,i], plot = F)[1:15][[1]]))$coefficients[2]
        #rho <- -log(fit[2]^2/fit[1]*as.numeric(acf(dataset[,i], type="covariance")[1][[1]]))
        params[(4*(i-1)+3):(4*(i-1)+4)] <- c(log(rho), log(temp_kappa))
        #print(rho)
      }
      
      return(params)
    }
    
    # loglik_pl <- function(times, values, deltas, logscale=T){
    #   d <- length(values[1,])
    #   return(function(params){
    #             return(pl_final(times = times,
    #                           values = values,
    #                           deltas = deltas,
    #                           params = params,
    #                           logscale = logscale))
    #           })
    # }
    # 
    # loglik_pl_univ <- function(times, values, delta, trf, logscale=T, lambda){
    #   return(function(params){
    #     temp <- pl_final_univ(times = times,
    #                           values = values,
    #                           delta = delta,
    #                           params = params,
    #                           trf = trf,
    #                           logscale = logscale)
    #     print(params)
    #     print(-temp+lambda*sum(params^2)+lambda/2*sum((1/params)^2))
    #     return(-temp+lambda*sum(params^2)+lambda*sum((1/params)^2))
    #   })
    # }
    
    loglikelihood_pl_univ_ic <- function(times, values, delta, model_vars_names, fixed_names, fixed_params, lambda, logscale=T, transformation=F){
      return(function(params){
        temp <- pl_univ(times = times,
                        values = values,
                        delta = delta,
                        fixed_names = fixed_names,
                        fixed_params = fixed_params,
                        params = params,
                        model_vars_names = model_vars_names,
                        logscale = logscale,
                        transformation = transformation)
        
        # J hat
        j_hat <- 0.0
        
        if(lambda != 0.0){
          # estimating the Information Criterion
          # H hat
          f_whole <- function(params){
            return(pl_univ(times = times,
                          values = values,
                          delta = delta,
                          fixed_names = fixed_names,
                          fixed_params = fixed_params,
                          params = params,
                          model_vars_names = model_vars_names,
                          transformation = transformation,
                          logscale = T))
          }
          
          h_hat <- -hessian.f(f = f_whole,
                              params = params)
     
         
          
          n_max <- length(values)
          #n_j_estimation <- n_max-delta+1
          n_j_estimation <- sqrt(length(times))/2
          #cat("njesti:", n_j_estimation, "\n")
          for(start_block in 1:(n_j_estimation)){
            f_block <- function(params){
              return(pl_univ(times = times[floor((start_block-1)*n_j_estimation + 1):floor((start_block)*n_j_estimation + 1)],
                             values = values[floor((start_block-1)*n_j_estimation + 1):floor((start_block)*n_j_estimation + 1)],
                             delta = delta,
                             fixed_names = fixed_names,
                             fixed_params = fixed_params,
                             params = params,
                             model_vars_names = model_vars_names,
                             transformation = transformation,
                             logscale = T))
            }
            tp <- grad.f(f = f_block, params)
            j_hat <- j_hat + tp %o% tp
  
            #cat("params:", params, "\n")
          }
         
            if(length(j_hat) > 1){
              j_hat <- lambda*sum(diag(j_hat %*% matlib::Inverse(h_hat, tol=1e-2)))
            }else{
              j_hat <- lambda * j_hat / h_hat
            }
        }
          #cat("j_hat:", (j_hat), "\n")
          #print(-temp + j_hat)
          #print(j_hat)
          #print(j_hat)
          return(-temp + j_hat)
        
      })
    }
    
    # loglik_pl_univ_ic <- function(times, values, delta, trf, logscale=T, lambda){
    #   return(function(params){
    #     temp <- pl_final_univ(times = times,
    #                           values = values,
    #                           delta = delta,
    #                           params = params,
    #                           trf = trf,
    #                           logscale = logscale)/1.0
    #     #print(params)
    #     
    #     # estimating the Information Criterion
    #     # H hat
    #     f_whole <- function(params){
    #       pl_final_univ(times = times,
    #                     values = values,
    #                     delta = delta,
    #                     params = params,
    #                     trf = trf,
    #                     logscale = logscale)/1.0
    #     }
    #     
    #     h_hat <- -hessian.f(f = f_whole,
    #                         params = params)
    #     
    #     # J hat
    #     j_hat <- 0.0
    #    
    #     n_max <- length(values)
    #     #n_j_estimation <- n_max-delta+1
    #     n_j_estimation <- sqrt(length(times))
    #       #cat("njesti:", n_j_estimation, "\n")
    #     for(start_block in 1:(n_j_estimation)){
    #       f_block <- function(params){
    #         return(pl_final_univ(times = times[floor((start_block-1)*n_j_estimation + 1):floor((start_block)*n_j_estimation + 1)],
    #                              values = values[floor((start_block-1)*n_j_estimation + 1):floor((start_block)*n_j_estimation + 1)],
    #                              delta = delta,
    #                              params = params,
    #                              trf = trf,
    #                              logscale = logscale)/1.0)
    #       }
    #       tp <- grad.f(f = f_block, params)
    #       j_hat <- j_hat + tp %o% tp
    #       #j_hat <- j_hat
    #     }
    #       #cat("params:", params, "\n")
    #       
    #       j_hat <- lambda*sum(diag(j_hat %*% matlib::inv(h_hat)))
    #       cat("j_hat:", (j_hat), "\n")
    #       #print(-temp + j_hat)
    #       return(-temp + j_hat)
    #     
    #   })
    # }
    
    mom_trawl <- function(data){
      # Hosking Wallis (1987) notation
      data_wo_zeros <- data
      k <- 0.5 * (mean(data)^2/var(data) - 1)
      alpha_original <- 0.5 * mean(data) * (mean(data)^2/var(data) + 1)
      
      beta <- -alpha_original/k
      alpha <- -1/k
      
      rho <- -log(as.numeric(acf(data, type = "covariance")[1][[1]]) * beta^2/alpha)

      p <- length(which(data > 0.0)) / length(data)
      kappa <- beta*(p^{-1.0/alpha} - 1.0)
      return(c(alpha, beta, rho, kappa))
    }
  }
  
  # Prepare data
  {
    colnames(bl) <- col_names
    bl <- subset(bl, select=c("Date", "O3", "NO2x", "SO2x", "NO", "CO"))
    invalid_data_indices <- which(is.na(bl)) %% length(bl$O3)
    
    bl_times <- matrix(rep(seq(1, length.out = length(bl[,1])), length(bl[1,])),
                       ncol=length(bl[1,]))
    bl_times <- bl_times[-invalid_data_indices,]
    bl_times <- bl_times[-39413,]
    bl <- bl[-invalid_data_indices,]
    bl <- bl[-39413,]
    bl_values <- bl[,-1]
    SAME_QUANTILE <- 0.7
    SAME_QUANTILE_G <- 0.7
    
    bl_q <-  c(
      quantile(bl$O3, SAME_QUANTILE, na.rm = T)[[1]],
      quantile(bl$NO2x, SAME_QUANTILE, na.rm = T)[[1]],
      quantile(bl$SO2x, SAME_QUANTILE, na.rm = T)[[1]],
      quantile(bl$NO, SAME_QUANTILE, na.rm = T)[[1]],
      quantile(bl$CO, SAME_QUANTILE, na.rm = T)[[1]]
    )
    
    bl_pdf <- c(
      approxfun(density(bl$O3, na.rm = T)),
      approxfun(density(bl$NO2x, na.rm = T)),
      approxfun(density(bl$SO2x, na.rm = T)),
      approxfun(density(bl$NO, na.rm = T)),
      approxfun(density(bl$CO, na.rm = T))
    )
    
    bl_cdf <- c(
      CDF(density(bl$O3, na.rm = T)),
      CDF(density(bl$NO2x, na.rm = T)),
      CDF(density(bl$SO2x, na.rm = T)),
      CDF(density(bl$NO, na.rm = T)),
      CDF(density(bl$CO, na.rm = T))
    )
    
    bl_tail <- list(
      O3=which(bl$O3 > bl_q[1]),
      NO2=which(bl$NO2x > bl_q[2]),
      SO2=which(bl$SO2x > bl_q[3]),
      NO=which(bl$NO > bl_q[4]),
      CO=which(bl$CO > bl_q[5])
    )
    
    bl_thres <- cbind(
      pmax(bl$O3 - bl_q[1], 0.0),
      pmax(bl$NO2x - bl_q[1], 0.0),
      pmax(bl$SO2x - bl_q[1], 0.0),
      pmax(bl$NO - bl_q[1], 0.0),
      pmax(bl$CO - bl_q[1], 0.0)
                      )

    # Fitting marginals
    bl_rho <- 0.3
    bl_kappa <- 10
    bl_deltas <- rep(4, 5)
    params_init <- create_ig_trawl(dataset = bl_values, 
                                   rho = bl_rho, 
                                   kappa = bl_kappa,
                                   quantiles = bl_q,
                                   trf = F)
  }

  fn_to_optim <- loglik_pl(times = bl_times,
                           values = bl_values,
                           deltas = bl_deltas)
  #fn_to_optim(params_init)
  #optim(par = params_init, fn = fn_to_optim)

  # Univariate case
  ## O3 without transfo
  params_to_work_with <- params_init[1:4]
  fn_to_optim <- loglikelihood_pl_univ_ic(times = bl_times[,1],
                                          values = bl_thres[,1],
                                          delta = bl_deltas[1],
                                          lambda = 1.0,
                                          model_vars_names = univ_model_vars_names,
                                          fixed_names = c(),
                                          fixed_params = c(),
                                          logscale = T,
                                          transformation = T)

  fn_to_optim(params_to_work_with)

  fn_to_optim <- loglikelihood_pl_univ_ic(times = bl_times[,1],
                                          values = bl_thres[,1],
                                          delta = bl_deltas[1],
                                          lambda = 1.0,
                                          model_vars_names = univ_model_vars_names,
                                          fixed_names = c("alpha", "beta", "kappa"),
                                          fixed_params = params_to_work_with[c(1,2,4)],
                                          logscale = T,
                                          transformation = T)

  fn_to_optim(params_to_work_with[3])

  # fn_to_optim <- loglik_pl_univ_ic(times = bl_times[,1],
  #                                  values = bl_thres[,1],
  #                                  delta = bl_deltas[1],
  #                                  lambda = 1.0,
  #
  #                                  logscale = T,
  #                                  trf = T)
  # fn_to_optim(params_to_work_with)
  #

  lower_limit <- c(
    1,
    0.001,
    -5,
    -5
  )
  lower_limit

  upper_limit <- c(
    200,
    10000.0,
    -0.01,
    5.0
  )
  upper_limit
  o3_univ <- optim(par=params_to_work_with[3],
                  fn = fn_to_optim,
                  control = list(trace=5, maxit=10, pgtol=1e-3, parscale=rep(1.0, 1)),
                  method = "L-BFGS-B",
                  lower = -10,
                  upper = 10)
  o3_univ$par
  # 
  # te <- seq(-5,5,length.out = 200)
  # tv <- rep(0, length(te))
  # for(index in 1:length(te)){
  #   tv[index] <- fn_to_optim(te[index])
  # }  
  # 
  # 
  # 