setwd("C:/Users/Valentin/Documents/MRes/Data")
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

data_bloomsbury <- read.csv("bloomsbury_1994_1998.csv", sep = ",",
                            as.is = T)
colnames(data_bloomsbury) <- col_names
data_bloomsbury <- subset(data_bloomsbury, select=c("Date", "O3", "NO2x", "SO2x", "NO", "CO"))
valid_data_indices <- which(is.na(data_bloomsbury)) %% length(data_bloomsbury$O3)
SAME_QUANTILE <- 0.7
# Cleaning the dataset from NA values.
# data_bloomsbury <- data_bloomsbury[valid_data_indices,]

require('evir')
hat_pdf <- function(x, xi, beta, q, pdf_f, ecdf_f){
  x <- as.numeric(x)
  if(x > q){
    return((1.0 - ecdf_f(q)) * dgpd(x, xi = xi, beta = beta)[1])
  } else{
    return(pdf_f(x))
  }
}

vec_hat_pdf <- function(x_vec, xi, beta, q, pdf_f, ecdf_f, without_zero=F){
  data_to_use <- x_vec[which(x_vec > q)]
  res <- rep(0, length(data_to_use))
  for(i in 1:length(res)){
    res[i] <- hat_pdf(data_to_use[i], xi, beta, q, pdf_f, ecdf_f)
  }
  if(without_zero){
    return(res[res>0])  
  } else {
    return(res)
  }
}



vec_full_hat_pdf <- function(x_vec, xi, beta, q, pdf_f, ecdf_f, without_zero=F){
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


hat_cdf <- function(x, xi, beta, q, ecdf_f){
  x <- as.numeric(x)
  if(x > q ){
    print("OTHER USED")
    return((1.0 - ecdf_f(q)) * pgpd(x, xi = xi, beta = beta)[1])
  } else{
    return(ecdf_f(x))
  }
}

vec_full_hat_cdf <- function(x_vec, xi, beta, q, ecdf_g, without_zero=F){
  res <- rep(0, length(x_vec))
  for(i in 1:length(res)){
    if(is.na(x_vec[i])){
      res[i] <- NA
    } else {
      res[i] <- hat_cdf(x=x_vec[i], xi=xi, beta=beta, q=q, ecdf_f=ecdf_g)
    }
  }
  if(without_zero){
    return(res[res>0])  
  } else {
    return(res)
  }
}
library(lubridate)
# PDF, CDF, QUANTILEs
{
  month_data <- month(as.POSIXlt(data_bloomsbury$Date, format="%Y-%m-%d"))
  winter_points <- which(month_data >= 2 & month_data >= 11)
  early_summer_points <- which(month_data >= 4 & month_data <= 7)
  
  bloom_winter <- data_bloomsbury[winter_points,]
  bloom_summer <- data_bloomsbury[early_summer_points,]
  library('fExtremes')
  par(mfrow=c(1,1), mai=c(.9,0.8,0.3,0.1))
  mePlot(bloom_winter$NO2x)
  plot(bloom_winter$O3, bloom_winter$NO2x)
  plot(bloom_summer$O3, bloom_summer$NO2x)
  
  q_winter_O3 <- quantile(bloom_winter$O3, SAME_QUANTILE, na.rm = T)[[1]]
  q_winter_NO2 <- quantile(bloom_winter$NO2x, SAME_QUANTILE, na.rm = T)[[1]]
  q_winter_SO2 <- quantile(bloom_winter$SO2x, SAME_QUANTILE, na.rm = T)[[1]]
  q_winter_NO <- quantile(bloom_winter$NO, SAME_QUANTILE, na.rm = T)[[1]]
  q_winter_CO <- quantile(bloom_winter$CO, SAME_QUANTILE, na.rm = T)[[1]]
  
  q_winter_full <- c(q_winter_O3,
                     q_winter_NO2,
                     q_winter_SO2,
                     q_winter_NO,
                     q_winter_CO)
  
  q_summer_O3 <- quantile(bloom_summer$O3, SAME_QUANTILE, na.rm = T)[[1]]
  q_summer_NO2 <- quantile(bloom_summer$NO2x, SAME_QUANTILE, na.rm = T)[[1]]
  q_summer_SO2 <- quantile(bloom_summer$SO2x, SAME_QUANTILE, na.rm = T)[[1]]
  q_summer_NO <- quantile(bloom_summer$NO, SAME_QUANTILE, na.rm = T)[[1]]
  q_summer_CO <- quantile(bloom_summer$CO, SAME_QUANTILE, na.rm = T)[[1]]
  
  q_summer_full <- c(q_summer_O3,
                     q_summer_NO2,
                     q_summer_SO2,
                     q_summer_NO,
                     q_summer_CO)
  
  ecdf_winter_O3 <- CDF(density(bloom_winter$O3, na.rm = T))
  ecdf_winter_NO2 <- CDF(density(bloom_winter$NO2x, na.rm = T))
  ecdf_winter_SO2 <- CDF(density(bloom_winter$SO2x, na.rm = T))
  ecdf_winter_NO <- CDF(density(bloom_winter$NO, na.rm = T))
  ecdf_winter_CO <- CDF(density(bloom_winter$CO, na.rm = T))
  
  ecdf_winter_full <- list(ecdf_winter_O3, ecdf_winter_NO2, ecdf_winter_SO2, ecdf_winter_NO, ecdf_winter_CO)
  
  pdf_winter_O3 <- approxfun(density(bloom_winter$O3, na.rm = T))
  pdf_winter_NO2 <- approxfun(density(bloom_winter$NO2x, na.rm = T))
  pdf_winter_SO2 <- approxfun(density(bloom_winter$SO2x, na.rm = T))
  pdf_winter_NO <- approxfun(density(bloom_winter$NO, na.rm = T))
  pdf_winter_CO <- approxfun(density(bloom_winter$CO, na.rm = T))
  
  pdf_winter_full <- list(pdf_winter_O3, pdf_winter_NO2, pdf_winter_SO2, pdf_winter_NO, pdf_winter_CO)
  library(spatstat)
  # ecdf_summer_O3 <- ecdf(bloom_summer$O3)
  # ecdf_summer_NO2 <- ecdf(bloom_summer$NO2x)
  # ecdf_summer_SO2 <- ecdf(bloom_summer$SO2x)
  # ecdf_summer_NO <- ecdf(bloom_summer$NO)
  # ecdf_summer_CO <- ecdf(bloom_summer$CO)
  
  ecdf_summer_O3 <- CDF(density(bloom_summer$O3, na.rm = T))
  ecdf_summer_NO2 <- CDF(density(bloom_summer$NO2x, na.rm = T))
  ecdf_summer_SO2 <- CDF(density(bloom_summer$SO2x, na.rm = T))
  ecdf_summer_NO <- CDF(density(bloom_summer$NO, na.rm = T))
  ecdf_summer_CO <- CDF(density(bloom_summer$CO, na.rm = T))
  
  ecdf_summer_full <- list(ecdf_summer_O3, ecdf_summer_NO2, ecdf_summer_SO2, ecdf_summer_NO, ecdf_summer_CO)
  
  pdf_summer_O3 <- approxfun(density(bloom_summer$O3, na.rm = T))
  pdf_summer_NO2 <- approxfun(density(bloom_summer$NO2x, na.rm = T))
  pdf_summer_SO2 <- approxfun(density(bloom_summer$SO2x, na.rm = T))
  pdf_summer_NO <- approxfun(density(bloom_summer$NO, na.rm = T))
  pdf_summer_CO <- approxfun(density(bloom_summer$CO, na.rm = T))
  
  pdf_summer_full <- list(pdf_summer_O3, pdf_summer_NO2, pdf_summer_SO2, pdf_summer_NO, pdf_summer_CO)
  
  tail_w_indices_O3 <- which(bloom_winter$O3 > q_winter_O3 & !is.na(bloom_winter$O3))
  tail_w_indices_NO2 <- which(bloom_winter$NO2x > q_winter_NO2 & !is.na(bloom_winter$NO2x))
  tail_w_indices_SO2 <- which(bloom_winter$SO2x > q_winter_SO2 & !is.na(bloom_winter$SO2x))
  tail_w_indices_NO <- which(bloom_winter$NO > q_winter_NO & !is.na(bloom_winter$NO))
  tail_w_indices_CO <- which(bloom_winter$CO > q_winter_CO & !is.na(bloom_winter$CO))
  
  tail_s_indices_O3 <- which(bloom_summer$O3 > q_summer_O3 & !is.na(bloom_summer$O3))
  tail_s_indices_NO2 <- which(bloom_summer$NO2x > q_summer_NO2 & !is.na(bloom_summer$NO2x))
  tail_s_indices_SO2 <- which(bloom_summer$SO2x > q_summer_SO2 & !is.na(bloom_summer$SO2x))
  tail_s_indices_NO <- which(bloom_summer$NO > q_summer_NO & !is.na(bloom_summer$NO))
  tail_s_indices_CO <- which(bloom_summer$CO > q_summer_CO & !is.na(bloom_summer$CO))
  
  tail_w_full <- list(tail_w_indices_O3,
                   tail_w_indices_NO2,
                   tail_w_indices_SO2,
                   tail_w_indices_NO,
                   tail_w_indices_CO)
  
  tail_s_full <- list(tail_s_indices_O3,
                   tail_s_indices_NO2,
                   tail_s_indices_SO2,
                   tail_s_indices_NO,
                   tail_s_indices_CO)
}

library('gPdtest')
xi_beta_gpd <- function(data, quantile)
{
  gpd.fit(data[data > q])
}

params <- rep(0, 10)

### Initial guess

par(mfrow=c(1,5), mai=c(.9,0.8,0.3,0.1))
for(i in 1:5)
  mePlot(bloom_winter[,i+1][tail_w_full[i][[1]]])

# Creating the initial guess for marginal MLE in params in winter
create_initial_guess_marginals <- function(dataset, tail_ind, test_fit=F)
  {
  params <- rep(0,10)
  # O3
  fit <- gpd.fit(dataset$O3[tail_ind[1][[1]]], method = "combined")
  params[1:2] <- as.numeric(fit)
  
  # NO2
  fit <- gpd.fit(dataset$NO2x[tail_ind[2][[1]]], method = "amle")
  params[3:4] <- as.numeric(fit)
  
  # SO2
  fit <- gpd.fit(dataset$SO2x[tail_ind[3][[1]]], method = "amle")
  params[5:6] <- as.numeric(fit)
  
  # NO
  fit <- gpd.fit(dataset$NO[tail_ind[4][[1]]], method = "amle")
  params[7:8] <- as.numeric(fit)
  
  # CO
  fit <- gpd.fit(dataset$CO[tail_ind[5][[1]]], method = "combined")
  params[9:10] <- as.numeric(fit)
  
  if(test_fit)
  {
    gpd.test(dataset$O3[tail_ind[1][[1]]])
    gpd.test(dataset$O3[tail_ind[2][[1]]])
    gpd.test(dataset$O3[tail_ind[3][[1]]])
    gpd.test(dataset$O3[tail_ind[4][[1]]])
    gpd.test(dataset$NO[tail_ind[5][[1]]])
    
  }
  
  return(params)
}

create_y <- function(dataset, quantiles)
{
  dataset_y <- data.frame(dataset$Date)
  colnames(dataset_y) <- "Date"
  dataset_y$O3 <- pmax(dataset$O3 - quantiles[1], 0.0)
  dataset_y$NO2x <- pmax(dataset$NO2x - quantiles[2], 0.0)
  dataset_y$SO2x <- pmax(dataset$SO2x - quantiles[3], 0.0)
  dataset_y$NO <- pmax(dataset$NO - quantiles[4], 0.0)
  dataset_y$CO <- pmax(dataset$CO - quantiles[5], 0.0)
  
  return(dataset_y)
}

bloom_winter_y <- create_y(bloom_winter, q_winter_full)
bloom_summer_y <- create_y(bloom_summer, q_summer_full)

loglikelihood_marginals <- function(params, dataset, quantiles, pdfs, ecdfs, tail_index)
{
  result <- 0
  temp <- 0
  d <- length(dataset[1,]) - 1
  for(i in 2:2)
  {
    temp <- vec_full_hat_pdf(x_vec = dataset[,i+1][tail_index[i][[1]]], xi = params[2*(i-1)+1], beta = params[2*i],
                            q = quantiles[i], pdf_f = pfds[i][[1]], ecdf_f = ecdfs[i][[1]], without_zero = T)
    temp <- temp[temp > 0]
    result <- result + sum(log(temp)) 
  }
  
  return(result)
}

# Fitting with winter data:
params_w <- create_initial_guess_marginals(bloom_winter, tail_w_full)

encap_loglikelihood <- function(params){
  return(-loglikelihood_marginals(params, bloom_winter, q_winter_full, pdf_winter_full, ecdf_winter_full, tail_w_full))
}

optim_winter <- optim(par = params_w[3:4], fn = encap_loglikelihood, method = "L-BFGS-B",
                      lower = rep(c(0.1,0.1), 2), upper = rep(c(15,120), 2),
                      control = list(trace=1))
params_winter <- optim_winter$par
params_winter

library('DEoptim')
optim_winter <- DEoptim(encap_loglikelihood, 
                        lower = rep(c(0.001,0.001), 5),
                        upper = rep(c(50,10.0), 5),
                        DEoptim.control(NP = 10*20,
                                        itermax = 400, 
                                        F = 1.2, 
                                        CR = 0.7,
                                        trace=1,
                                        parallelType = 1,
                                        parVar=c("loglikelihood_marginals",
                                                 "bloom_winter",
                                                 "vec_full_hat_pdf",
                                                 "tail_w_full",
                                                 "hat_pdf",
                                                 "q_winter_full",
                                                 "ecdf_winter_full",
                                                 "dgpd")))
params_winter <- optim_winter$par


# Fitting with summer data:
params <- create_initial_guess_marginals(bloom_summer, tail_s_full)

encap_loglikelihood <- function(params){
  return(-loglikelihood_marginals(params, bloom_summer, q_summer_full, pdf_summer_full, ecdf_summer_full, tail_s_full))
}

#optim_summer <- optim(par = params, fn = encap_loglikelihood)
#params_summer <- optim_summer$par


create_z <- function(params, dataset, quantiles, ecdfs)
{
  result <- data.frame(dataset$Date)
  
  d <- length(dataset[1,]) - 1
  for(i in 1:d)
  {
    temp <- vec_full_hat_cdf(x_vec = dataset[,i+1], xi = params[2*(i-1)+1], beta = params[2*i],
                             q = quantiles[i], ecdf_g = ecdfs[i][[1]])
    result[,1+i] <- -log(-log(temp))
  }
  
  colnames(result) <- c("Date", "O3", "NO2x", "SO2x", "NO", "CO")
  
  return(result)
}



bloom_winter_z <- create_z(params=params_w, dataset=bloom_winter, quantiles=q_winter_full, ecdfs=ecdf_winter_full)
bloom_summer_z <- create_z(params_summer, bloom_summer, q_summer_full, ecdf_summer_full)


{
  q_gumbel_winter_O3 <- quantile(bloom_winter_z$O3, SAME_QUANTILE, na.rm = T)[[1]]
  q_gumbel_winter_NO2 <- quantile(bloom_winter_z$NO2x, SAME_QUANTILE, na.rm = T)[[1]]
  q_gumbel_winter_SO2 <- quantile(bloom_winter_z$SO2x, SAME_QUANTILE, na.rm = T)[[1]]
  q_gumbel_winter_NO <- quantile(bloom_winter_z$NO, SAME_QUANTILE, na.rm = T)[[1]]
  q_gumbel_winter_CO <- quantile(bloom_winter_z$CO, SAME_QUANTILE, na.rm = T)[[1]]
  
  q_gumbel_w_full <- c(q_gumbel_winter_O3,
                       q_gumbel_winter_NO2,
                       q_gumbel_winter_SO2,
                       q_gumbel_winter_NO,
                       q_gumbel_winter_CO)
  
  # ecdf_gumbel_winter_O3 <- ecdf(bloom_winter_z$O3)
  # ecdf_gumbel_winter_NO2 <- ecdf(bloom_winter_z$NO2x)
  # ecdf_gumbel_winter_SO2 <- ecdf(bloom_winter_z$SO2x)
  # ecdf_gumbel_winter_NO <- ecdf(bloom_winter_z$NO)
  # ecdf_gumbel_winter_CO <- ecdf(bloom_winter_z$CO)
  
  ecdf_gumbel_winter_O3 <- CDF(density(bloom_winter_z$O3, na.rm = T))
  ecdf_gumbel_winter_NO2 <- CDF(density(bloom_winter_z$NO2x, na.rm = T))
  ecdf_gumbel_winter_SO2 <- CDF(density(bloom_winter_z$SO2x, na.rm = T))
  ecdf_gumbel_winter_NO <- CDF(density(bloom_winter_z$NO, na.rm = T))
  ecdf_gumbel_winter_CO <- CDF(density(bloom_winter_z$CO, na.rm = T))
  
  ecdfs_winter_gumbel <- c(ecdf_gumbel_winter_O3,
                           ecdf_gumbel_winter_NO2,
                           ecdf_gumbel_winter_SO2,
                           ecdf_gumbel_winter_NO,
                           ecdf_gumbel_winter_CO)
  
  pdf_gumbel_winter_O3 <- approxfun(density(bloom_winter_z$O3, na.rm = T))
  pdf_gumbel_winter_NO2 <- approxfun(density(bloom_winter_z$NO2x, na.rm = T))
  pdf_gumbel_winter_SO2 <- approxfun(density(bloom_winter_z$SO2x, na.rm = T))
  pdf_gumbel_winter_NO <- approxfun(density(bloom_winter_z$NO, na.rm = T))
  pdf_gumbel_winter_CO <- approxfun(density(bloom_winter_z$CO, na.rm = T))
  
  pdfs_winter_gumbel <- c(pdf_gumbel_winter_O3,
                          pdf_gumbel_winter_NO2,
                          pdf_gumbel_winter_SO2,
                          pdf_gumbel_winter_NO,
                          pdf_gumbel_winter_CO)
}


ecdfs_gumbel_w_hat_cdf <- c(function(x){return(hat_cdf(x, xi = params_winter[1], beta = params_winter[2],
                                                       q = q_gumbel_w_full[1], ecdf_f = ecdf_winter_full[1][[1]]))})

for(i in 2:5){
  ecdfs_gumbel_w_hat_cdf <- c(ecdfs_gumbel_w_hat_cdf, 
                              function(x){return(hat_cdf(x, xi = params_winter[2*(i-1)+1], beta = params_winter[2*i],
                                                         q = q_gumbel_w_full[i], ecdf_f = ecdf_winter_full[i][[1]]))})
  
}


make_indices <- function(dataset, quantiles, ecdf_fs){
  C_indexes <- which(dataset$O3 > quantiles[1] 
                     | dataset$NO2x > quantiles[2]
                     | dataset$SO2x > quantiles[3]
                     | dataset$NO > quantiles[4]
                     | dataset$CO > quantiles[5])
  
  C1_indexes <- C_indexes[which(ecdf_fs[1][[1]](dataset$O3[C_indexes]) 
                                > ecdf_fs[2][[1]](dataset$NO2x[C_indexes]) 
                                & 
                                  ecdf_fs[1][[1]](dataset$O3[C_indexes])
                                > ecdf_fs[3][[1]](dataset$SO2x[C_indexes])
                                & 
                                  ecdf_fs[1][[1]](dataset$O3[C_indexes])
                                > ecdf_fs[4][[1]](dataset$NO[C_indexes])
                                &
                                  ecdf_fs[1][[1]](dataset$O3[C_indexes])
                                > ecdf_fs[5][[1]](dataset$CO[C_indexes]))]
  
  
  C2_indexes <- C_indexes[which(ecdf_fs[2][[1]](dataset$NO2x[C_indexes]) 
                                > ecdf_fs[1][[1]](dataset$O3[C_indexes]) 
                                & 
                                  ecdf_fs[2][[1]](dataset$NO2x[C_indexes]) 
                                > ecdf_fs[3][[1]](dataset$SO2x[C_indexes])
                                & 
                                  ecdf_fs[2][[1]](dataset$NO2x[C_indexes]) 
                                > ecdf_fs[4][[1]](dataset$NO[C_indexes])
                                &
                                  ecdf_fs[2][[1]](dataset$NO2x[C_indexes]) 
                                > ecdf_fs[5][[1]](dataset$CO[C_indexes]))]
  
  C3_indexes <- C_indexes[which(ecdf_fs[3][[1]](dataset$SO2x[C_indexes]) 
                                > ecdf_fs[1][[1]](dataset$O3[C_indexes]) 
                                & 
                                  ecdf_fs[3][[1]](dataset$SO2x[C_indexes]) 
                                > ecdf_fs[2][[1]](dataset$NO2x[C_indexes])
                                & 
                                  ecdf_fs[3][[1]](dataset$SO2x[C_indexes]) 
                                > ecdf_fs[4][[1]](dataset$NO[C_indexes])
                                &
                                  ecdf_fs[3][[1]](dataset$SO2x[C_indexes]) 
                                > ecdf_fs[5][[1]](dataset$CO[C_indexes]))]
  
  C4_indexes <- C_indexes[which(ecdf_fs[4][[1]](dataset$NO[C_indexes])
                                > ecdf_fs[1][[1]](dataset$O3[C_indexes]) 
                                & 
                                  ecdf_fs[4][[1]](dataset$NO[C_indexes])
                                > ecdf_fs[2][[1]](dataset$NO2x[C_indexes])
                                & 
                                  ecdf_fs[4][[1]](dataset$NO[C_indexes])
                                > ecdf_fs[3][[1]](dataset$SO2x[C_indexes]) 
                                &
                                  ecdf_fs[4][[1]](dataset$NO[C_indexes])
                                > ecdf_fs[5][[1]](dataset$CO[C_indexes]))]
  
  C5_indexes <- C_indexes[which(ecdf_fs[5][[1]](dataset$CO[C_indexes])
                                > ecdf_fs[1][[1]](dataset$O3[C_indexes]) 
                                & 
                                  ecdf_fs[5][[1]](dataset$CO[C_indexes])
                                > ecdf_fs[2][[1]](dataset$NO2x[C_indexes])
                                & 
                                  ecdf_fs[5][[1]](dataset$CO[C_indexes])
                                > ecdf_fs[3][[1]](dataset$SO2x[C_indexes]) 
                                &
                                  ecdf_fs[5][[1]](dataset$CO[C_indexes])
                                > ecdf_fs[4][[1]](dataset$NO[C_indexes]))]
  return(list("C"=C_indexes,"C1"=C1_indexes,"C2"=C2_indexes,"C3"=C3_indexes,"C4"=C4_indexes,"C5"=C5_indexes))
}



multi_ind_w_full <- make_indices(bloom_winter_z, q_gumbel_w_full, ecdfs_winter_gumbel)
multi_ind_w_full

single_cond_loglikelihood <- function(params, dataset, indices_list, quantiles, index){
  indices <- indices_list[index+1][[1]]
  # Should be given Gumbel marginals
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


  over_threshold_indices <- which(data_to_use[,index] > quantiles[index])
  # b|i(y)
  
  temp_mat_y <- matrix(rep(data_to_use[,index], d-1), ncol = d-1) # create a matrix n x (d-1) with d-1 the same col
  temp_mat_b <- temp_mat_y^rep(bi, each=n)

  # Computing \mu_{|i}(y) = a_{|i}*y|i + mu|i * y|i^b|i and sigma|i(y)
  temp_mat_mu <- temp_mat_y * rep(ai, each=n) + temp_mat_b * rep(mui, each=n)
  temp_mat_sigma <- temp_mat_b * rep(sigmai, each=n)
  
  
  temp_mat_b <- as.vector(temp_mat_b[over_threshold_indices,])
  temp_mat_mu <- as.vector(temp_mat_mu[over_threshold_indices,])
  temp_mat_sigma <- as.vector(temp_mat_sigma[over_threshold_indices,])
  indices_pos_sigma <- which(temp_mat_sigma > 0)
  
  data_vector <- as.vector(as.matrix(data_to_use[,-index]))
  
  
  temp <- sum(log(temp_mat_sigma[indices_pos_sigma]))
  result <- temp + sum(((data_vector[indices_pos_sigma] - temp_mat_mu[indices_pos_sigma])/temp_mat_sigma[indices_pos_sigma])^2)

  return(result)
}

full_cond_loglikelihood <- function(params, dataset, indices_list, quantiles)
{
  result <- 0
  d <- length(dataset[1,])-1
  
  for(i in 1:d)
  {
    result <- result + single_cond_loglikelihood(params=params,
                              dataset=bloom_winter_z,
                              indices_list=multi_ind_w_full,
                              quantiles=q_gumbel_w_full,
                              index=i)  
  }
  
  return(result)
}


par(mfrow=c(1,1), mai=c(.9,0.8,0.3,0.1))
ok <- full_cond_loglikelihood(params=seq(from = 0.1, to=0.51, length.out=80),
                                dataset=bloom_winter_z,
                                indices_list=multi_ind_w_full,
                                quantiles=q_gumbel_w_full) 

full_dependence_structure <- function(params){
  return(full_cond_loglikelihood(params=params,
                          dataset=bloom_winter_z,
                          indices_list=multi_ind_w_full,
                          quantiles=q_gumbel_w_full)) 
}

#test_params <- create_initial_guess_marginals(bloom_winter_z, multi_ind_w_full[1][[1]])

full_params <- seq(from = 0, to=2, length.out=80)
# one set of params:
one_set_params <- c(rep(0.5,4), rep(.5,4), rep(.5,4), rep(.5,4))
full_params <- rep(one_set_params, 5)


de_lower <- c(rep(0, 4), rep(-1, 4), rep(0, 4), rep(-10, 4))
de_lower <- rep(de_lower, 5)

de_upper <- c(rep(1, 4), rep(1, 4), rep(+10, 4), rep(+10, 4))
de_upper <- rep(de_upper, 5)

system.time(
optim_full <- optim(fn=full_dependence_structure, par=full_params, lower = de_lower, upper = de_upper,
                    method = "L-BFGS-B", control=list(maxit=60, trace=1)))

tic = Sys.time()
optim_full2 <- DEoptim(fn=full_dependence_structure, lower=de_lower, upper=de_upper,
                    control=DEoptim.control(NP=80*15, itermax=100, parallelType=1, parVar=c("full_cond_loglikelihood",
                                                                                 "bloom_winter_z",
                                                                                 "multi_ind_w_full",
                                                                                 "q_gumbel_w_full",
                                                                                 "single_cond_loglikelihood"), trace=T))
print(Sys.time() - tic) # 14.37 mins


result_param <- t(matrix(optim_full2$optim$bestmem, nrow = 4))
result_param

# ai
result_param[4*(0:4) + 1,]

# bi
result_param[4*(0:4) + 2,]

# sigmai
result_param[4*(0:4) + 3,]

# mui
result_param[4*(0:4) + 4,]


optim_full$counts

a21 <- result_param[1,1]
b21 <- result_param[2,1]

a21y <- a21 * bloom_winter_z[intersect(multi_ind_w_full[2][[1]], which(bloom_winter_z$O3 > q_gumbel_w_full[1]+.1)), 2] 
b21y <- bloom_winter_z[intersect(multi_ind_w_full[2][[1]], which(bloom_winter_z$O3 > q_gumbel_w_full[1]+.1)), 2]^b21
plot(a21y, b21y)
