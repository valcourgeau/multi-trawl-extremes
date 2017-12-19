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
bl <- bl[1:10000,]
# list of packages
library(spatstat)
library(evir)
library(gPdtest)
library(DEoptim)
library(fExtremes)

# Functions
{
  hat_pdf <- function(x, xi, beta, q, pdf_f, ecdf_f){
    x <- as.numeric(x)
    if(x > q){
      if(xi < 0 & x >= -beta/xi){
        return(0.0)
      }
      return((1.0 - ecdf_f(q)) * (evir::dgpd(x-q, xi = xi, beta = beta)[1]))
    } else{
      return(pdf_f(x))
    }
  }
  
  hat_cdf <- function(x, xi, beta, q, ecdf_f){
    x <- as.numeric(x)
    if(x > q){
      return(1 - (1.0 - ecdf_f(q)) * (1-evir::pgpd(x-q, xi = xi, beta = beta)[1]))
    } else{
      return(ecdf_f(x))
    }
  }
  
  vectorised_hat_pdf <- function(x_vec, xi, beta, q, pdf_f, ecdf_f, without_zero=F){
    res <- rep(0, length(x_vec))
    
    for(i in 1:length(res)){
      if(is.na(x_vec[i])){
        res[i] <- NA
      } else {
        res[i] <- hat_pdf(x_vec[i], xi, beta, q, pdf_f, ecdf_f)
      }
    }
    
    if(without_zero){
      return(res[res>0])  
    } else {
      return(res)
    }
  }
  
  vectorised_hat_cdf <- function(x_vec, xi, beta, q, ecdf_f, without_zero=F){
    res <- rep(0, length(x_vec))
    
    for(i in 1:length(res)){
      if(is.na(x_vec[i])){
        res[i] <- NA
      } else {
        res[i] <- hat_cdf(x_vec[i], xi, beta, q, ecdf_f)
      }
    }
    
    if(without_zero){
      return(res[res>0])  
    } else {
      return(res)
    }
  }
  
  create_initial_guess_marginals <- function(dataset, tail_ind, qs)
  {
    d <- length(dataset[1,]) - 1
    params <- rep(0,2*d)
    
    for(i in 1:d){
      fit <- evir::gpd(dataset[,i+1][tail_ind[i][[1]]], method = "ml", qs[i])
      params[(2*(i-1)+1):(2*(i-1)+2)] <- fit$par.ests
    }
    
    return(params)
  }
  
  loglikelihood_marginals <- function(dataset, quantiles, pdfs, cdfs, tail_index, neg=T)
  {
    return(function(params){
              result <- 0
              temp <- 0
              d <- length(dataset[1,]) - 1
              for(i in 1:2)
              {
                temp <- vectorised_hat_pdf(x_vec = dataset[,i+1][bl_tail[i][[1]]], xi = (params[2*(i-1)+1]), 
                                           beta = (params[2*(i-1)+2]),
                                           q = quantiles[i], pdf_f = pdfs[i][[1]], 
                                           ecdf_f = cdfs[i][[1]], without_zero = F)
                temp[which(temp == 0.0)] <- 1e-16
                #temp <- temp[which(temp > 0)]
                temp <- temp[which(!is.nan(temp))]
                #print(length((temp)))
                result <- result + sum(log(temp)) 
              }
              
              if(neg){
                result <- -result
              }
              
              return(result)
          }
    )
  }
  
  create_gumbel <- function(dataset, quantiles, cdfs, tail_index, params){
    d <- length(dataset[1,]) - 1
    tranformed <- dataset
    for(i in 1:d){
      tranformed[,i+1] <- -log(-log(cdfs[i][[1]](dataset[,i+1])))
    }
    
    return(tranformed)
  }
  
  gumbel_single_loglik <- function(dataset, quantiles, params, tail_index, index){
    # Should be given Gumbel marginals
    # index: index on which we condition
    
    indices <- tail_index[index][[1]]
    d <- length(dataset[1,])-1
    data_to_use <- dataset[indices,2:(d+1)]
    n <- length(data_to_use[,1])
    
    mat_params <- t(matrix(params, nrow=d-1))
    
    result <- 0
    temp <- 0
    temp_mat_b <- 0
    temp_mat_mu <- 0
    temp_mat_sigma <- 0


    ai <-       mat_params[(index-1)*4+1,]
    bi <-       mat_params[(index-1)*4+2,]
    sigmai <-   mat_params[(index-1)*4+3,]
    mui <-      mat_params[(index-1)*4+4,]
    
    #over_threshold_indices <- which(data_to_use[,index] > quantiles[index])
    
    # b|i(y)
    temp_mat_y <- matrix(rep(data_to_use[,index], d-1), ncol = d-1) # create a matrix n x (d-1) with d-1 the same col
    temp_mat_b <- t(t(temp_mat_y)^rep(bi, each=n))
    
    # Computing \mu_{|i}(y) = a_{|i}*y|i + mu|i * y|i^b|i and sigma|i(y)
    temp_mat_mu <- t(t(temp_mat_y) * rep(ai, each=n) + t(temp_mat_b) * rep(mui, each=n))
    temp_mat_sigma <- t(t(temp_mat_b) * rep(sigmai, each=n))


    temp_mat_b <- as.vector(temp_mat_b)
    temp_mat_mu <- as.vector(temp_mat_mu)
    temp_mat_sigma <- as.vector(temp_mat_sigma)
    indices_pos_sigma <- which(temp_mat_sigma > 1e-6)
    
    data_vector <- as.vector(as.matrix(data_to_use[,-index]))

    temp <- sum(log(temp_mat_sigma[indices_pos_sigma]))
    temp <- temp + 10*(n-length(indices_pos_sigma))
    result <- temp + sum(((data_vector[indices_pos_sigma] - temp_mat_mu[indices_pos_sigma])/temp_mat_sigma[indices_pos_sigma])^2)

    return(result)
  }
  
  loglikelihood_gumbel <- function(dataset, quantiles, tail_index, neg=F){
    d <- length(dataset[1,]) - 1
    return(function(params){
      result <- 0.0
      for(i in 1:d){
        result <- result + gumbel_single_loglik(dataset, quantiles, params, tail_index, i)
      }
      
      if(neg){
        result <- -result
      }
      
      return(result)
    })
  }
}

# Prepare data
{
  colnames(bl) <- col_names
  bl <- subset(bl, select=c("Date", "O3", "NO2x", "SO2x", "NO", "CO"))
  valid_data_indices <- which(is.na(bl)) %% length(bl$O3)
  bl <- bl[-valid_data_indices,]
  bl <- bl[-39413,]
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
  
  # Fitting marginals
  params_init <- create_initial_guess_marginals(bl, bl_tail, bl_q)
  
  # Create Gumbel marginals
  gb <- create_gumbel(dataset = bl, quantiles = bl_q, cdfs = bl_cdf, params = params_init)
  
  gb_q <-  c(
    quantile(gb$O3, SAME_QUANTILE_G, na.rm = T)[[1]],
    quantile(gb$NO2x, SAME_QUANTILE_G, na.rm = T)[[1]],
    quantile(gb$SO2x, SAME_QUANTILE_G, na.rm = T)[[1]],
    quantile(gb$NO, SAME_QUANTILE_G, na.rm = T)[[1]],
    quantile(gb$CO, SAME_QUANTILE_G, na.rm = T)[[1]]
  )
  
  gb_tail <- list(
    O3=which(gb$O3 > gb_q[1]),
    NO2=which(gb$NO2x > gb_q[2]),
    SO2=which(gb$SO2x > gb_q[3]),
    NO=which(gb$NO > gb_q[4]),
    CO=which(gb$CO > gb_q[5])
  )
}

# Testing the functions
{
  # Testing with SO2x
  params_init <- create_initial_guess_marginals(bl, bl_tail, bl_q)
  fit <- gpd(bl$SO2x[bl_tail[3][[1]]], method="ml", threshold = bl_q[3])
  qqparetoPlot(bl$SO2x[bl_tail[3][[1]]] - bl_q[3])
  qqparetoPlot(bl$SO2x[bl_tail[3][[1]]], xi = fit$par.ests[1][[1]], threshold = bl_q[3])
  xi_SO2 <- fit$par.ests[1][[1]]
  beta_SO2 <- fit$par.ests[2][[1]]
  plot(sort(bl$SO2x), vectorised_hat_cdf(sort(bl$SO2x),
                                   xi = xi_SO2,
                                   beta = beta_SO2,
                                   q = bl_q[3],
                                   ecdf_f = bl_cdf[3][[1]]),
       type = 'l')
  plot(sort(bl$SO2x), vectorised_hat_pdf(sort(bl$SO2x),
                                         xi = xi_SO2,
                                         beta = beta_SO2,
                                         q = bl_q[3],
                                         pdf_f = bl_pdf[3][[1]],
                                         ecdf_f = bl_cdf[3][[1]]),
       type = 'l')
  plot(sort(bl$SO2x[bl_tail[3][[1]]]), vectorised_hat_pdf(sort(bl$SO2x[bl_tail[3][[1]]]),
                                         xi = xi_SO2,
                                         beta = beta_SO2,
                                         q = bl_q[3],
                                         pdf_f = bl_pdf[3][[1]],
                                         ecdf_f = bl_cdf[3][[1]]),
       type = 'l')
  
  gumbel_single_loglik(dataset = gb, 
                       quantiles = gb_q, 
                       params = seq(0,1,length.out=80),
                       tail_index = gb_tail,
                       index=1)
  
  fn_to_optim <- loglikelihood_gumbel(dataset = gb, 
                       quantiles = gb_q,
                       tail_index = gb_tail, neg=F)
  fn_to_optim(seq(0.1, 20., length.out = 80))
  
  res <- DEoptim(fn = fn_to_optim, lower=rep(c(rep(0.001,4), rep(-1,4), rep(0.001,4), rep(-5,4)), 5),
                 upper=rep(c(rep(1,4), rep(1,4), rep(2,4), rep(5,4)), 5), DEoptim.control(itermax = 20,
                                                          parallelType = 1,
                                                          NP = 15*80,
                                                          parVar = c("loglikelihood_marginals",
                                                                     "vectorised_hat_pdf",
                                                                     "bl_tail",
                                                                     "bl_q",
                                                                     "hat_pdf",
                                                                     "gumbel_single_loglik")))
  
  
  result_param <- t(matrix(res$optim$bestmem, nrow = 4))
  result_param
  
  # ai
  result_param[4*(0:4) + 1,]
  
  # bi
  result_param[4*(0:4) + 2,]
  
  # sigmai
  result_param[4*(0:4) + 3,]
  
  # mui
  result_param[4*(0:4) + 4,]
  
  #write.csv(result_param, "result_on_10000_pts.csv")
  result_param <- read.csv("C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/results/tawn/result_on_10000_pts.csv")
  result_param <- matrix(as.matrix(result_param[,-1]), ncol = 4)
  
  # page 519: bootstrap convex hull?
  res$counts
  
  a21 <- result_param[1,]
  b21 <- result_param[2,]
  
  a21y <- a21 * gb[gb_tail[1][[1]], 2] 
  b21y <- gb[gb_tail[1][[1]], 2]^b21
  
  # write.csv(a21y, file = "a21y.csv")
  # write.csv(b21y, file = "a21y.csv")
  
  a21y <- read.csv("a21y.csv")[,2]
  b21y <- read.csv("b21y.csv")[,2]
  # the higher y_1, the less volatile y_2
  
  cor(a21y, b21y)
  
  summary(a21y)
  summary(b21y)
  plot(a21y, b21y)
  
  
  # Given O3
  a_given_1 <- result_param[1,]
  b_given_1 <- result_param[2,]
  
  a_given_1y <- matrix(rep(gb[gb_tail[1][[1]], 2], 4), ncol=4)*rep(a_given_1, length(gb_tail[1][[1]]))
  b_given_1y <- t(t(matrix(rep(gb[gb_tail[1][[1]], 2], 4), ncol=4))^rep(b_given_1, length(gb_tail[1][[1]])))
  
  par(mfrow=c(1,4), mai=c(.9,0.8,0.3,0.1))
  for(i in 1:4)
    plot(a_given_1y[,i], b_given_1y[,i])
  par(mfrow=c(1,1), mai=c(.9,0.8,0.3,0.1))
  # write.csv(a_given_1y, file = "a_given_1y.csv")
  # write.csv(b_given_1y, file = "b_given_1y.csv")
  
  cor(a_given_1y, b_given_1y)
  
  # Given NO2
  a_given_2 <- result_param[5,]
  b_given_2 <- result_param[6,]
  
  a_given_2y <- matrix(rep(gb[gb_tail[2][[1]], 3], 4), ncol=4)*rep(a_given_2, length(gb_tail[2][[1]]))
  b_given_2y <- t(t(matrix(rep(gb[gb_tail[2][[1]], 3], 4), ncol=4))^rep(b_given_2, length(gb_tail[2][[1]])))
  
  par(mfrow=c(1,4), mai=c(.9,0.8,0.3,0.1))
  for(i in 1:4)
    plot(a_given_2y[,i], b_given_2y[,i])
  par(mfrow=c(1,1), mai=c(.9,0.8,0.3,0.1))
  
  # write.csv(a_given_2y, file = "a_given_2y.csv")
  # write.csv(b_given_2y, file = "b_given_2y.csv")
  
  cor(a_given_2y, b_given_2y)
  
  # Given SO2
  a_given_3 <- result_param[9,]
  b_given_3 <- result_param[10,]
  
  a_given_3y <- matrix(rep(gb[gb_tail[3][[1]], 4], 4), ncol=4)*rep(a_given_3, length(gb_tail[3][[1]]))
  b_given_3y <- t(t(matrix(rep(gb[gb_tail[3][[1]], 4], 4), ncol=4))^rep(b_given_3, length(gb_tail[3][[1]])))
  
  par(mfrow=c(1,4), mai=c(.9,0.8,0.3,0.1))
  for(i in 1:4)
    plot(a_given_3y[,i], b_given_3y[,i], ylim = c(0,max(b_given_3y[,i])))
  par(mfrow=c(1,1), mai=c(.9,0.8,0.3,0.1))
  
  # write.csv(a_given_3y, file = "a_given_3y.csv")
  # write.csv(b_given_3y, file = "b_given_3y.csv")

  cor(a_given_3y, b_given_3y)
}

# testing to fit marginals
{
  # # Fitting marginals
  # params_init <- create_initial_guess_marginals(bl, bl_tail, bl_q)
  # 
  # params_init[1:4]
  # # params_init[2*(0:4) + 2] <- log(params_init[2*(0:4) + 2])
  # params_init <- 0.2*(params_init)
  # params_init
  # fn_to_optim <- loglikelihood_marginals(dataset = bl,
  #                         quantiles = bl_q,
  #                         pdfs = bl_pdf,
  #                         cdfs = bl_cdf,
  #                         tail_index = bl_tail,
  #                         neg=T)
  # fn_to_optim(params_init)
  # 
  # res <- optim(par = params_init, fn = fn_to_optim,
  #       control=list(trace=2, maxit=200), method = "BFGS")
  # res$par
  # 
  # res <- DEoptim(fn = fn_to_optim, lower=rep(c(-1,0.001), 2),
  #                upper=rep(c(1,150), 2), DEoptim.control(itermax = 100,
  #                                                         parallelType = 1,
  #                                                         NP = 20*10,
  #                                                         parVar = c("loglikelihood_marginals",
  #                                                                    "vectorised_hat_pdf",
  #                                                                    "bl_tail",
  #                                                                    "bl_q",
  #                                                                    "hat_pdf")))
}

