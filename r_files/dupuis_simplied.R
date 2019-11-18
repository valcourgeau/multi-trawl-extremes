# packages
require(hypergeo)
require(ev.trawl)
require(lubridate)
require(magrittr)
require(rlist)

# setwd("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/r_files/")
# source("utils.R")
# source("multi_ev.R")

setwd('~/GitHub/multi-trawl-extremes/r_files/')
source('multi_ev.R')
source('draw_acf.R')

load("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data/merged-datasets/clean_east_light_data_v2.Rda")



DupuisSimplified <- function(data_u, mult_fac=c(0.3, 3), cl=NULL){
  params <- rep(0, 6)
  params[1:2] <- CustomMarginalMLE(data_u)
  p <- length(which(data_u > 0))/length(data_u)
  params[4] <- (1 - p^{params[1]}) * params[2]/abs(params[1])
  params[5] <- 1.0/params[1]
  params[6] <- params[2]/ abs(params[1]) - params[4]
  
  depth <- 10
  kk <- acf(data_u, lag.max = depth-1, plot=F)
  
  # print(params[1])
  alpha_tmp <- params[5]
  beta_tmp <- params[6]
  kappa_tmp <- params[4]
  # if(params[1] > 0){
  #   alpha_tmp <- params[5]
  #   beta_tmp <- params[6]
  #   kappa_tmp <- params[4]
  # }else{
  #   alpha_tmp <- 2
  #   beta_tmp <- 1
  #   kappa_tmp <- (p^{-1/alpha_tmp} - 1) * beta_tmp
  # }
  
  # cat('alpha', alpha_tmp, '\n')
  # cat('beta', beta_tmp, '\n')
  # cat('kappa', kappa_tmp, '\n')
  # 
  
  n_trials <- 10
  mae_tab <- rep(0, n_trials)
  mse_tab <- rep(0, n_trials)
  # rho_tab <- seq(0.01, (3), length.out = n_trials)
  index <- 1
  
  lin_rho <- line(x = c(0, 1:(depth-2)), log(kk$acf[1:(depth-1)] %>% abs))
  rho_tmp <- abs(lin_rho$coefficients[2])
  rho_tab <- seq(log(rho_tmp*mult_fac[1]), log(rho_tmp*mult_fac[2]), length.out = n_trials) %>% exp
  
  # plot(c(0.05, 1:(depth-1)), kk$acf, ylim=c(0,1))
  # print(kk$acf)
  
  for(rho_iter in rho_tab){
    # print(index)
    
    if(!is.null(cl)){
      parallel::clusterExport(cl, c('acf_trawl', 'alpha_tmp', 'beta_tmp', 'kappa_tmp', 'rho_iter'))
      acf_vals <- parallel::parLapply(cl, X = c(0.05, 1:(depth-1)), fun = function(h){
        return(acf_trawl(h, alpha = alpha_tmp, beta = beta_tmp, kappa = kappa_tmp, 
                         rho = rho_iter, delta = 0.1, end_seq = 50))})
      acf_vals <- unlist(acf_vals)
    }else{
      acf_vals <- vapply(c(0.05, 1:(depth-1)), function(h){
        acf_trawl(h, alpha = alpha_tmp, beta = beta_tmp, kappa = kappa_tmp, 
                  rho = rho_iter, delta = 0.5, end_seq = 50)}, 1)
    }
    
    # plot(x = c(0.05, 1:(depth-1)), log(acf_vals))
    # lines(lin_rho$fitted.values)
    mae_tab[index] <- sum(abs(kk$acf - acf_vals)) 
    mse_tab[index] <- sum(((kk$acf - acf_vals)^2)) #+ 

    # cat('rho', rho_iter, '\n')
    # cat('MSE', (mse_tab[index]), '\n')
    # cat('MAE', (mae_tab[index]), '\n')
    # 
    # if(index %% 2 == 0){
    #   acf_vals %>% (function(x){lines(c(0.05, 1:(depth-1)), x, col = 'red')})
    # }
    index <- index + 1
  }
  
  params[3] <- rho_tab[which.min(mse_tab)]
  names(params) <- c('xi', 'sigma', 'rho', 'kappa', 'alpha', 'beta')
  return(params)
}

# pp <- DupuisSimplified(exc_test_jittered[,12])
# acf(exc_test_jittered[,12])
# acf_vals <- vapply(c(0.05, 1:(depth-1)), function(h){
#   acf_trawl(h, alpha = pp[5], beta = pp[6], kappa = pp[4], 
#             rho = pp[3], delta = 0.1, end_seq = 30)}, 1) %>% lines

EstimateBootstrap <- function(data, n, method='subsample'){
  n_sample <- nrow(data)
  n_vars <- ncol(data)
  n_sub <- 10*round(n_sample^{2/3})
  # print(n_sub)
  start_pts <- sample(x = 1:(n_sample - n_sub), size = n, replace = T)
  cores <- parallel::detectCores(logical = FALSE)
  cl <- parallel::makeCluster(cores*2-1)
 
  # cores <- parallel::detectCores(logical = FALSE)
  # cl <- parallel::makeCluster(cores*2-1)
  # kappa_tmp <- 10
  # beta_tmp <- 5
  # alpha_tmp <- 1
  # rho_iter <-1 
  # parallel::clusterExport(cl, c('acf_trawl', 'alpha_tmp', 'beta_tmp', 'kappa_tmp', 'rho_iter'))
  # system.time(
  #   res1.p <- parallel::parLapply(cl, X = 1:(10-1), fun = function(h){
  #     print(alpha_tmp)
  #     return(acf_trawl(h, alpha = alpha_tmp, beta = beta_tmp, kappa = kappa_tmp, 
                       # rho = rho_iter, delta = 0.1, end_seq = 50))}))
  
  sol <- array(0, dim = c(n, n_vars, 6))
  for(i in 1:n){
    cat('Sample', i, '\n')
    sol[i,,] <- t(apply(data[start_pts[i]:(start_pts[i]+n_sub),], MARGIN = 2, FUN = function(x){DupuisSimplified(data_u = x, cl = NULL)})) 
    # print(sol[i,,])
  }
  
  parallel::stopCluster(cl)
  
  # res_mean <- matrix(0, nrow = n_vars, ncol=6)
  # for(i in 1:n_vars){
  #   res_mean[i,] <- apply(sol[,i,], MARGIN = 2, FUN = mean)
  # }
  # return(res_mean)
  
  res <- list()
  res[['mean']] <- apply(sol, MARGIN = c(2,3), mean)
  res[['sd']] <- apply(sol, MARGIN = c(2,3), sd)
  return(res)
}

boot_east <- EstimateBootstrap(exc_test_jittered, n = 500) # 13:54
# rlist::list.save(boot_east, 'bootstrap_32_east_n_500.RData')
boot_east[['sd']]
params_east <- apply(exc_test_jittered, MARGIN = 2, FUN = function(x){DupuisSimplified(data_u = x)})
# rlist::list.save(list(mle=params_east), 'final_params_32_east.RData')

# par(mfrow=c(5,5))
# for(i in 1:25){
#   init_guess <- eva::gpdFit(exc_test_jittered[,i], 0.0)$par.sum$Estimate[2:1]
#   hist(exc_test_jittered[,i][exc_test_jittered[,i] > 0], breaks = 100, probability = T)
#   lines(0:1500/100, eva::dgpd(0:1500/100, scale = params_east[2,i], shape = params_east[1,i]), col = 'red')
#   lines(0:1500/100, eva::dgpd(0:1500/100, scale = init_guess[2], shape = init_guess[1]), col = 'blue')
# }


tail_coefficients <- function(data, h){
  data_bin <- (data > 0.0)
  data_ar <- array(0, dim = c(2, nrow(data)[1]-h, ncol(data)))
  
  
  data_ar[1,,] <- data_bin[1:(nrow(data)-h),]
  data_ar[2,,] <- data_bin[(h+1):nrow(data_bin),]
  
  sub_sum <- apply(data_ar, MARGIN = 2, FUN = function(sub_arr){
    return(sub_arr[1,] %*% t(sub_arr[2,]))
  })
  
  dim_sum <- dim(sub_sum)
  sub_sum_aug <- array(0, dim=c(sqrt(dim_sum[1]), sqrt(dim_sum[1]), dim_sum[2]))

  for(i in 1:dim_sum[2]){
    sub_sum_aug[,,i] <- matrix(sub_sum[,i], byrow = T, nrow = dim(data_ar)[3])
  }
  
  rm(sub_sum)
  
  # mean_mat <- matrix(0, nrow=dim(sub_sum_aug)[1], ncol=dim(sub_sum_aug)[2])
  # for(i in 1:dim_sum[2]){
  #   mean_mat <- 
  # }

  sub_sum_aug <- apply(sub_sum_aug, MARGIN = c(1,2), FUN = mean)
  sub_sum_aug <- sub_sum_aug / apply(data_bin, MARGIN = 2, FUN = mean)
  colnames(sub_sum_aug) <- colnames(data)
  rownames(sub_sum_aug) <- colnames(data)
  
  res <- list()
  return(sub_sum_aug)
}

tail_coefficients_stack <- function(data, h.max=1){
  if(length(h.max) == 1){
    return(vapply(1:h.max,FUN = function(h){tail_coefficients(data = data, h = h)}, matrix(0, nrow = ncol(data), ncol = ncol(data))))  
  }else{
    return(vapply(h.max,FUN = function(h){tail_coefficients(data = data, h = h)}, matrix(0, nrow = ncol(data), ncol = ncol(data))))
  }
  
}

res_tc <- tail_coefficients(exc_test_jittered, 5)
tc_stack <- tail_coefficients_stack(exc_test_jittered, 10)

horizons_tron_east <- c(1,2,3,6,12,24,48)
tc_stack <- tail_coefficients_stack(exc_test_jittered, horizons_tron_east)


tron_east_loaded <- rlist::list.load('~/GitHub/multi-trawl-extremes/r_files/tron_east_ligh_1_to_72_4th_tron.RData')
AEP_cond_prob <- vapply(horizons_tron_east, function(h){tron_east_loaded[[h]]$mean[19,]}, rep(0, length(tron_east_loaded[[1]]$mean[19,]))) %>% t
AEP_tc_stack <- tc_stack[,31,] %>% t

par(mfrow=c(3,5))
for(i in 1:15){
  plot(horizons_tron_east, AEP_cond_prob[,i], type='b', lwd=3,
       ylim=c(0,1.2*max(c(AEP_cond_prob[,i], AEP_tc_stack[,i]))),
       col = 'red', ylab=colnames(AEP_tc_stack)[i], 
       cex.lab=1.5, cex=1.3, cex.axis=1.5,
       xlab='Horizon')
  lines(horizons_tron_east, AEP_tc_stack[,i], type='b', lwd=3, cex=1.3)
  abline(h=0.04, lty=3, col='darkgrey', lwd=2)
}
