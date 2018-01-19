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
  bl <- bl[1:5000,]
  
  setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
  # list of packages
  library(spatstat)
  library(fExtremes)
  library(gPdtest)
  library(DEoptim)
  library(evir)
  library('optimr')
  
  
  # source files
  source('pairwise_latent_trawl.R')
  
  # Function
  {
    create_ig_trawl <- function(dataset, quantiles, rho, kappa)
    {
      d <- length(dataset[1,]) - 1
      params <- rep(0,4*d)
      
      for(i in 1:d){
        x_to_work_on <- dataset[,i]
        fit <- evir::gpd(x_to_work_on, method = "ml", threshold = quantiles[i])
        #print(fit@fit$par.ests)
        params[(4*(i-1)+1):(4*(i-1)+2)] <- log(fit$par.ests)
        params[(4*(i-1)+2)] <- log(params[(4*(i-1)+2)])
        alpha <- params[4*(i-1)+1]
        if(alpha >= 0){
          p <- length(which(dataset[,i] > quantiles[i])) / length(dataset[,i])
          cat(i, "P > 0 = ",p, "\n")
          temp_kappa <- (1-p)/p
        }else{
          temp_kappa <- kappa
        }
        
        p <- length(which(dataset[,i] > quantiles[i])) / length(dataset[,i])
        cat(i, "P > 0 = ",p, "\n")
        temp_kappa <- (1-p)/p
        
        params[(4*(i-1)+3):(4*(i-1)+4)] <- c(log(rho), log(temp_kappa))
      }
      
      return(params)
    }
    
    loglik_pl <- function(times, values, deltas, logscale=T){
      d <- length(values[1,])
      return(function(params){
                return(pl_final(times = times,
                              values = values,
                              deltas = deltas,
                              params = params,
                              logscale = logscale))
              })
    }
    
    loglik_pl_univ <- function(times, values, delta, trf, logscale=T){
      return(function(params){
        temp <- pl_final_univ(times = times,
                              values = values,
                              delta = delta,
                              params = params,
                              trf = trf,
                              logscale = logscale)
        print(params)
        return(-temp)
      })
    }
    
    loglik_pl_univ_with_kappa <- function(times, values, delta, kappa, trf, logscale=T){
      return(function(params){
        temp <- pl_final_univ_with_kappa(times = times,
                              values = values,
                              delta = delta,
                              kappa = kappa,
                              params = params,
                              trf = trf,
                              logscale = logscale)
        print(params)
        print(-temp)
        return(-temp)
      })
    }
  }
  
  # Prepare data
  {
    colnames(bl) <- col_names
    bl <- subset(bl, select=c("Date", "O3", "NO2x", "SO2x", "NO", "CO"))
    invalid_data_indices <- which(is.na(bl)) %% length(bl$O3)
    
    bl_times <- matrix(rep(seq(1, length.out = length(bl[,1])), length(bl_values[1,])),
                       ncol=length(bl_values[1,]))
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
                                   quantiles = bl_q)
  }
  
  fn_to_optim <- loglik_pl(times = bl_times,
                           values = bl_values,
                           deltas = bl_deltas)
  #fn_to_optim(params_init)
  #optim(par = params_init, fn = fn_to_optim)
  
  # Univariate case
  
  
  
  # O3
  params_to_work_with <- params_init[1:4]
  fn_to_optim <- loglik_pl_univ(times = bl_times[,1],
                                values = bl_thres[,1],
                                delta = bl_deltas[1],
                                trf = T)
  lower_limit <- c(
    -10,
    -1,
    -5,
    -5
  )
  lower_limit
  
  upper_limit <- c(
    -0.001,
    4.0,
    0.0,
    5.0
  )
  upper_limit
  o3_univ <- optim(par=params_to_work_with, 
        fn = fn_to_optim, 
        control = list(trace=1, maxit=30), 
        method = "L-BFGS-B",
        lower = lower_limit,
        upper = upper_limit)
  o3_univ$par
  # o3_univ <- DEoptim(fn = fn_to_optim, 
  #                  DEoptim.control(trace=TRUE, NP=20*10),
  #                  lower = lower_limit,
  #                  upper = upper_limit)
  # o3_univ <- optimr(par = params_init[1:4],
  #                   fn = fn_to_optim,
  #                   method = "BFGS")
  
  
  hist(bl_thres[,1], breaks=50, probability = T)
  lines(0:5000/10, evir::dgpd(0:5000/10, xi = 1/params_to_work_with[1], beta = abs(exp(params_to_work_with[2])/params_to_work_with[1])) 
        / evir::dgpd(inv_g(0:500/10, xi = 1/params_to_work_with[1], sigma = abs(exp(params_to_work_with[2])/params_to_work_with[1]), kappa = params_to_work_with[4]), xi = 1, beta = 1+params_to_work_with[4]))  
  lines(0:5000/10, evir::dgpd(0:5000/10, xi = 1/o3_univ$par[1], beta = abs(exp(o3_univ$par[2])/o3_univ$par[1])) / evir::dgpd(inv_g(0:500/10, xi = 1/o3_univ$par[1], sigma = abs(exp(o3_univ$par[2])/o3_univ$par[1]), kappa = o3_univ$par[4]), xi = 1, beta = 1+o3_univ$par[4]))  
  
  fn_to_optim(params_to_work_with)
  fn_to_optim(o3_univ$par)
  acf(bl_thres[,1])
  
  # O3 with Kappa
  
  params_to_work_with <- params_init[1:3]
  fn_to_optim <- loglik_pl_univ_with_kappa(times = bl_times[,1],
                                values = bl_thres[,1],
                                delta = bl_deltas[1],
                                kappa = params_init[4],
                                trf = F)
  lower_limit <- c(
    -1,
    -2.0,
    -2
  )
  lower_limit
  
  upper_limit <- c(
    20.0,
    5.0,
    1.0
  )
  upper_limit
  o3_univ <- optim(par=params_to_work_with, 
                   fn = fn_to_optim, 
                   control = list(trace=1, maxit=30), 
                   method = "L-BFGS-B",
                   lower = lower_limit,
                   upper = upper_limit)
  
  # NO2x
  params_to_work_with <- params_init[5:8]
  fn_to_optim <- loglik_pl_univ(times = bl_times[,1],
                                values = bl_thres[,2],
                                delta = bl_deltas[2],
                                trf = T)
  lower_limit <- c(
    0.0001,
    -7.0,
    -5,
    0.0001
  )
  lower_limit
  
  upper_limit <- c(
    2.0,
    5.0,
    1.0,
    5
  )
  upper_limit
  no2_univ <- optim(par=params_to_work_with, 
                   fn = fn_to_optim, 
                   control = list(trace=1, maxit=30), 
                   method = "L-BFGS-B",
                   lower = lower_limit,
                   upper = upper_limit)
  
  
  hist(bl_thres[,2], breaks=50, probability = T)
  lines(0:5000/10, evir::dgpd(0:5000/10, xi = 1/params_to_work_with[1], beta = abs(exp(params_to_work_with[2])/params_to_work_with[1])) 
        / evir::dgpd(inv_g(0:500/10, xi = 1/params_to_work_with[1], sigma = abs(exp(params_to_work_with[2])/params_to_work_with[1]), kappa = params_to_work_with[4]), xi = 1, beta = 1+params_to_work_with[4]))  
  lines(0:5000/10, evir::dgpd(0:5000/10, xi = 1/no2_univ$par[1], beta = abs(exp(no2_univ$par[2])/no2_univ$par[1])) / evir::dgpd(inv_g(0:500/10, xi = 1/no2_univ$par[1], sigma = abs(exp(no2_univ$par[2])/no2_univ$par[1]), kappa = no2_univ$par[4]), xi = 1, beta = 1+no2_univ$par[4]))  
  
  fn_to_optim(params_to_work_with)
  fn_to_optim(o3_univ$par)
  acf(bl_thres[,1])
  
  hist(bl_thres[,1], breaks=10, probability = T)
  lines(0:500/10, evir::dgpd(0:500/10, xi = -1/-0.11, beta = abs(exp(0.69)/-0.11)) / evir::dgpd(inv_g(0:500/10, xi = -1/-0.11, sigma = abs(exp(0.69)/-0.11), kappa=exp(2.3)), xi = 1, beta = 1+exp(2.3)))  
