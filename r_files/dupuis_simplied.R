setwd('~/GitHub/multi-trawl-extremes/r_files/')
source('multi_ev.R')
source('draw_acf.R')

DupuisSimplified <- function(data_u){
  params <- rep(0, 6)
  params[1:2] <- CustomMarginalMLE(data_u)
  p <- length(which(data_u > 0))/length(data_u)
  params[4] <- (1 - p^{params[1]}) * params[2]/abs(params[1])
  params[5] <- 1.0/params[1]
  params[6] <- params[2]/ abs(params[1]) - params[4]
  
  depth <- 15
  kk <- acf(data_u, lag.max = depth-1, plot=F)
  
  alpha_tmp <- params[5]
  beta_tmp <- params[6]
  kappa_tmp <- params[4]
  
  n_trials <- 25
  mae_tab <- rep(0, n_trials)
  mse_tab <- rep(0, n_trials)
  rho_tab <- seq(log(1e-2), log(2), length.out = n_trials) %>% exp
  # rho_tab <- seq(0.01, (3), length.out = n_trials)
  index <- 1
  
  plot(c(0.05, 1:(depth-1)), kk$acf)
  # print(kk$acf)
  for(rho_iter in rho_tab){
    print(rho_iter)
    acf_vals <- vapply(c(0.05, 1:(depth-1)), function(h){
      acf_trawl(h, alpha = alpha_tmp, beta = beta_tmp, kappa = kappa_tmp, 
                rho = rho_iter, delta = 0.1, end_seq = 30)}, 1)
    mae_tab[index] <- sum(abs(kk$acf - acf_vals)) 
    mse_tab[index] <- sum(((kk$acf - acf_vals)^2)) #+ 

    if(index %% 2 == 0){
      acf_vals %>% (function(x){lines(c(0.05, 1:(depth-1)), x, col = 'red')})
    }
    index <- index + 1
  }
  
  cat('MSE', (mse_tab), '\n')
  cat('MSE', (mae_tab), '\n')
  params[3] <- rho_tab[which.min(mse_tab)]
  return(params)
}

DupuisSimplified(exc_test_jittered[,1])

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
  return(vapply(1:h.max,FUN = function(h){tail_coefficients(data = data, h = h)}, matrix(0, nrow = ncol(data), ncol = ncol(data))))
}

res_tc <- tail_coefficients(exc_test_jittered, 5)
tc_stack <- tail_coefficients_stack(exc_test_jittered, 10)


res_tc %>% dim
