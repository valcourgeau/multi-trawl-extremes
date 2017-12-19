setwd("C:/Users/Valentin/Documents/MRes/Data")
library("fExtremes")
library("evir")
library('gPdtest')
library(lubridate)
#source("implementation_new_heffernan_tawn.R")

# Colnames etc
{
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
}


data_bloomsbury <- read.csv("bloomsbury_1994_1998.csv", sep = ",",
                            as.is = T)

colnames(data_bloomsbury) <- col_names
data_bloomsbury <- subset(data_bloomsbury, select=c("Date", "O3", "NO2x", "SO2x", "NO", "CO"))
valid_data_indices <- which(is.na(data_bloomsbury)) %% length(data_bloomsbury$O3)

SAME_QUANTILE <- 0.70

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
  
  ecdf_winter_O3 <- ecdf(bloom_winter$O3)
  ecdf_winter_NO2 <- ecdf(bloom_winter$NO2x)
  ecdf_winter_SO2 <- ecdf(bloom_winter$SO2x)
  ecdf_winter_NO <- ecdf(bloom_winter$NO)
  ecdf_winter_CO <- ecdf(bloom_winter$CO)
  
  ecdf_winter_full <- list(ecdf_winter_O3, ecdf_winter_NO2, ecdf_winter_SO2, ecdf_winter_NO, ecdf_winter_CO)
  
  pdf_winter_O3 <- approxfun(density(bloom_winter$O3, na.rm = T))
  pdf_winter_NO2 <- approxfun(density(bloom_winter$NO2x, na.rm = T))
  pdf_winter_SO2 <- approxfun(density(bloom_winter$SO2x, na.rm = T))
  pdf_winter_NO <- approxfun(density(bloom_winter$NO, na.rm = T))
  pdf_winter_CO <- approxfun(density(bloom_winter$CO, na.rm = T))
  
  pdf_winter_full <- list(pdf_winter_O3, pdf_winter_NO2, pdf_winter_SO2, pdf_winter_NO, pdf_winter_CO)
  
  ecdf_summer_O3 <- ecdf(bloom_summer$O3)
  ecdf_summer_NO2 <- ecdf(bloom_summer$NO2x)
  ecdf_summer_SO2 <- ecdf(bloom_summer$SO2x)
  ecdf_summer_NO <- ecdf(bloom_summer$NO)
  ecdf_summer_CO <- ecdf(bloom_summer$CO)
  
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

create_y <- function(dataset, quantiles){
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

source('pairwise_latent_trawl.R')

d <- length(bloom_winter_y[1,]) - 1
bloom_winter_y_values <- bloom_winter_y[,2:(d+1)]
for(i in 1:d){
  bloom_winter_y_values <- bloom_winter_y_values[!is.na(bloom_winter_y_values[,i]),]
}
n <- length(bloom_winter_y_values[,1])
bloom_winter_y_times <- matrix(rep(seq(0,100,length.out = n), d), ncol=d)

default_rho <- 0.3
{
  indication_params <- rep(0,15)
  indication_params[1:2] <- gpd.fit(bloom_winter_y_values$O3[which(bloom_winter_y_values$O3 > 0.0)], method="combined")
  qqparetoPlot(bloom_winter_y_values$O3[which(bloom_winter_y_values$O3 > 0.0)], xi = indication_params[1])
  indication_params[3] <- default_rho
  
  indication_params[4:5] <- gpd.fit(bloom_winter_y_values$NO2x[which(bloom_winter_y_values$NO2x > 0.0)], method="combined")
  qqparetoPlot(bloom_winter_y_values$NO2x[which(bloom_winter_y_values$NO2x > 0.0)], xi = indication_params[4])
  indication_params[6] <- default_rho
  
  indication_params[7:8] <- gpd.fit(bloom_winter_y_values$SO2x[which(bloom_winter_y_values$SO2x > 0.0)], method="combined")
  qqparetoPlot(bloom_winter_y_values$SO2x[which(bloom_winter_y_values$SO2x > 0.0)], xi = indication_params[7])
  indication_params[9] <- default_rho
  
  indication_params[10:11] <- gpd.fit(bloom_winter_y_values$NO[which(bloom_winter_y_values$NO > 0.0)], method="combined")
  qqparetoPlot(bloom_winter_y_values$NO[which(bloom_winter_y_values$NO > 0.0)], xi = indication_params[10])
  indication_params[12] <- default_rho
  
  indication_params[13:14] <- gpd.fit(bloom_winter_y_values$CO[which(bloom_winter_y_values$CO > 0.0)], method="combined")
  qqparetoPlot(bloom_winter_y_values$CO[which(bloom_winter_y_values$CO > 0.0)], xi = indication_params[13])
  indication_params[15] <- default_rho
  
  # transform (xi, sigma) into (alpha, beta)
  cat("Original:\n")
  print(t(matrix(indication_params, nrow=3)))
  indication_params[0:4*3 + 1] <- 1/indication_params[0:4*3 + 1] # alpha = 1/xi
  indication_params[0:4*3 + 2] <- abs(indication_params[0:4*3 + 2] * indication_params[0:4*3 + 1]) # beta = sigma * alpha
  #indication_params[0:4*3 + 2][which(indication_params[0:4*3 + 2] < 0)] <- - indication_params[0:4*3 + 2][which(indication_params[0:4*3 + 2] < 0)] 
  cat("Transformed:\n")
  print(t(matrix(indication_params, nrow=3)))
}

source('finding_kappa.R')
for(i in 1:5){
  data_to_use <- bloom_winter_y_values[,i]
  to_optim <- likelihood_find_kappa(z = data_to_use, xi = indication_params[(i-1)*3+1], sigma = indication_params[(i-1)*3+2])
  to_optim
  to_optim(3)
  
  print(optim(to_optim, par = 0.5, method = "BFGS")$par)
}



initialise_marginal_transform <- function(dataset, kappas, methods, indication_params){
  par <- rep(0,15)
  d <- length(dataset[1,])

  for(i in 1:d){
    data_to_use <- dataset[,i][which(dataset[,i] > 0.0)]
    if(indication_params[(i-1)*3+1] < 0.0){
      # Construct F_(xi, simga)(data)
      fit_partly_transformed <- gpd.fit(data_to_use, method=methods[i])
      data_to_use_transformed <- pgpd(data_to_use, xi = fit_partly_transformed[1], beta = fit_partly_transformed[2])
  
      # Taking the inverse with 
      data_to_use_transformed <- qgpd(data_to_use_transformed, xi = 1, beta = 1/(1+kappas[i]))
      
      par[(i-1)*3 + 1] <- 1/fit_partly_transformed[1] # alpha = 1/xi
      par[(i-1)*3 + 2] <- abs(fit_partly_transformed[2] * par[(i-1)*3 + 1]) # beta = sigma * alpha
      par[(i-1)*3 + 3] <- default_rho
      dataset[,i][which(dataset[,i] > 0.0)] <- data_to_use_transformed
    }else{
      par[(i-1)*3 + 1] <- indication_params[(i-1)*3 + 1] # alpha = 1/xi
      par[(i-1)*3 + 2] <- indication_params[(i-1)*3 + 2]  # beta = sigma * alpha
      par[(i-1)*3 + 3] <- default_rho
    }
  }
  return(list("par" = par, "data" = dataset))
}


deltas <- rep(2, 5)
kappas <- rep(0, 5)
methods_gpd <- c("combined", "combined", "combined", "combined", "combined")
init_mt <- initialise_marginal_transform(bloom_winter_y_values, kappas, methods_gpd, indication_params)
t(matrix(init_mt$par, nrow=3))
hist(init_mt$data[,4][which(init_mt$data[,4] > 0)], breaks=50)

# Parameters init by component-wise fit
fn_to_optim <- function(params){
  temp <- -pairwise_full(bloom_winter_y_times,bloom_winter_y_values,
                         deltas, kappas, params)
  print(temp)
  return(temp)
}

fn_to_optim_univ <- function(params){
  temp <- -pairwise_likelihood_single_full_params(bloom_winter_y_times[,4][which(bloom_winter_y_values[,4] > 0.0)],
                                                  bloom_winter_y_values[,4][which(bloom_winter_y_values[,4] > 0.0)],
                         deltas[4], kappas[4], params)
  return(temp)
}

de_lower <- rep(c(0.001, 0.001, 0.001), 1)
de_upper <- rep(c(+5, 50, 1), 1)

library('DEoptim')

fn_to_optim_univ_all <- function(params){
  temp <- -pl_single_all_params(bloom_winter_y_times[,4][which(bloom_winter_y_values[,4] > 0.0)],
                                                  bloom_winter_y_values[,4][which(bloom_winter_y_values[,4] > 0.0)],
                                                  deltas[4], params)
  return(temp)
}

get_general_controls_univ <- function(itmax, n_NP)
{return(DEoptim.control(NP=n_NP, itermax=itmax, parallelType=1, parVar=c("pairwise_likelihood_single_full",
                                                                         "bloom_winter_y_values",
                                                                         "bloom_winter_y_times",
                                                                         "kappas",
                                                                         "transfos",
                                                                         "deltas",
                                                                         "pairwise_likelihood_single_pair",
                                                                         "pairwise_00_exp",
                                                                         "compute_B1_exp",
                                                                         "compute_B_inter_exp",
                                                                         "compute_B3_exp",
                                                                         "pairwise_00_1",
                                                                         "pairwise_00_2",
                                                                         "pairwise_10_exp",
                                                                         "inv_g",
                                                                         "pgpd",
                                                                         "qgpd",
                                                                         "compute_A_exp",
                                                                         "pairwise_10_1",
                                                                         "pairwise_10_2",
                                                                         "pairwise_10_2_1",
                                                                         "pairwise_10_2_2",
                                                                         "pairwise_10_2_3",
                                                                         "pairwise_10_2_4",
                                                                         "dgpd",
                                                                         "pairwise_11_exp",
                                                                         "pairwise_11_1",
                                                                         "pairwise_11_1_1",
                                                                         "pairwise_11_1_2",
                                                                         "pairwise_11_1_3",
                                                                         "pairwise_11_2",
                                                                         "pairwise_11_2_1",
                                                                         "pairwise_11_2_2",
                                                                         "pairwise_11_2_3",
                                                                         "pairwise_11_2_4",
                                                                         "pl_single_all_params"/,
                                                                         "pl_single_kappa_rho",
                                                                         "pairwise_likelihood_single_full_params_with_kappa"), trace=T))
}

run_deoptim_univ <- function(dataset, delta){
  fn_to_optim_univ_all <- function(params){
    temp <- -pl_single_all_params(bloom_winter_y_times[,1][which(bloom_winter_y_values[,1] > 0.0)],
                                  bloom_winter_y_values[,1][which(bloom_winter_y_values[,1] > 0.0)],
                                  deltas[4], params, F)
    return(temp)
  }
  de_lower <- rep(c(-10,    0.0, 0.001, 0.0001), 1)
  de_upper <- rep(c(0.001, 50,  10,    6), 1)
  system.time(
    optim_1 <- DEoptim(fn=fn_to_optim_univ_all, lower=de_lower, upper=de_upper,
                       control=get_general_controls_univ(20, 4*20))
    
  )
  
  return(optim_1)
}

run_deoptim_univ(bloom_winter_y_values, 4)
pl_single_kappa_rho

run_deoptim_univ_kappa_rho <- function(dataset, delta, alpha, beta){
  fn_to_optim_univ_all <- function(params){
    temp <- -pl_single_kappa_rho(bloom_winter_y_times[,4][which(bloom_winter_y_values[,4] > 0.0)],
                                  bloom_winter_y_values[,4][which(bloom_winter_y_values[,4] > 0.0)],
                                  deltas[4], params)
    return(temp)
  }
  de_lower <- rep(c(0.001,    0.0, 0.001, 0.0001), 1)
  de_upper <- rep(c(10, 50,  10,    6), 1)
  system.time(
    optim_1 <- DEoptim(fn=fn_to_optim_univ_all, lower=de_lower, upper=de_upper,
                       control=get_general_controls_univ(20, 4*20))
    
  )
  
  return(optim_1)
}

run_deoptim_univ(bloom_winter_y_values, 4)

### UNIVARIATE DEOPTIM
source('pairwise_latent_trawl.R')
# Define controls
get_general_controls <- function(itmax, n_NP)
  {return(DEoptim.control(NP=n_NP, itermax=itmax, parallelType=1, parVar=c("pairwise_likelihood_single_full",
                                                                                  "bloom_winter_y_values",
                                                                                  "bloom_winter_y_times",
                                                                                  "kappas",
                                                                                  "transfos",
                                                                                  "deltas",
                                                                                  "pairwise_likelihood_single_pair",
                                                                                  "pairwise_00_exp",
                                                                                  "compute_B1_exp",
                                                                                  "compute_B_inter_exp",
                                                                                  "compute_B3_exp",
                                                                                  "pairwise_00_1",
                                                                                  "pairwise_00_2",
                                                                                  "pairwise_10_exp",
                                                                                  "inv_g",
                                                                                  "pgpd",
                                                                                  "qgpd",
                                                                                  "compute_A_exp",
                                                                                  "pairwise_10_1",
                                                                                  "pairwise_10_2",
                                                                                  "pairwise_10_2_1",
                                                                                  "pairwise_10_2_2",
                                                                                  "pairwise_10_2_3",
                                                                                  "pairwise_10_2_4",
                                                                                  "dgpd",
                                                                                  "pairwise_11_exp",
                                                                                  "pairwise_11_1",
                                                                                  "pairwise_11_1_1",
                                                                                  "pairwise_11_1_2",
                                                                                  "pairwise_11_1_3",
                                                                                  "pairwise_11_2",
                                                                                  "pairwise_11_2_1",
                                                                                  "pairwise_11_2_2",
                                                                                  "pairwise_11_2_3",
                                                                                  "pairwise_11_2_4",
                                                                                  "pairwise_likelihood_single_full_params_with_kappa"), trace=T))
}

run_deoptim_univ <- function(indices, deltas, kappas){
  
  for(i in indices){
    fn_to_optim_univi <- function(params){
      temp <- -pairwise_likelihood_single_full_params(bloom_winter_y_times[,i][which(bloom_winter_y_values[,i] > 0.0)],
                                                     bloom_winter_y_values[,i][which(bloom_winter_y_values[,i] > 0.0)],
                                                     deltas[i], kappas[i], params, transformation=F)
      return(temp)
    }
    
    if(i == 1 || i == 3){
      de_lower <- rep(c(-10, 0.0, 0.001), 1)
      de_upper <- rep(c(-0.001, 50, 1), 1)
    }else{
      de_lower <- rep(c(0.01, 0.0, 0.001), 1)
      de_upper <- rep(c(10, 50, 1), 1)
    }

    
    system.time(
      optim_1 <- DEoptim(fn=fn_to_optim_univi, lower=de_lower, upper=de_upper,
                         control=get_general_controls(20, 3*20))
      
    )
  }
  
}


deltas <- rep(4, 5)
kappas <- rep(30, 5)

run_deoptim_univ(2, deltas, kappas)
run_deoptim_univ(4, deltas, kappas)
run_deoptim_univ(5, deltas, kappas)


run_deoptim_multivar <- function(indices, deltas, transfo){
  
    fn_to_optim_mult <- function(params){
      temp <- 0.0
      for(i in indices){
        temp <- temp-pairwise_likelihood_single_full_params_with_kappa(times=bloom_winter_y_times[,i][which(bloom_winter_y_values[,i] > 0.0)],
                                                        values=bloom_winter_y_values[,i][which(bloom_winter_y_values[,i] > 0.0)],
                                                        delta=deltas[i], params=params, transformation=transfo[i])
      }
      return(temp)
    }
    
    de_lower <- rep(c(0.01, 0.0, 0.001, 1.0), length(indices))
    de_upper <- rep(c(10, 50, 1, 80.0), length(indices))
    
    system.time(
      optim_1 <- DEoptim(fn=fn_to_optim_mult, lower=de_lower, upper=de_upper,
                         control=get_general_controls(50, 12*20))
      
    )
    
    return(optim_1)
}

transfos <- c(T,F,T,F,F)
res5 <- run_deoptim_multivar(c(2,4,5), deltas, transfos)
param_res5 <- t(matrix(res5$optim$bestmem, ncol=3))

par(mfrow=c(2,3), mai=c(.9,0.8,0.3,0.1))
hist(bloom_winter_y_values[,2][which(bloom_winter_y_values[,2] > 0.0)], breaks=15, probability = T)
lines(1:15000/10, dgpd(1:15000/10, xi = param_res5[1,1], beta = (param_res5[1,2]+param_res5[1,4])))

hist(bloom_winter_y_values[,4][which(bloom_winter_y_values[,4] > 0.0)], breaks=15, probability = T)
lines(1:25000/10, dgpd(1:25000/10, xi = param_res5[2,1], beta = (param_res5[2,2]+param_res5[2,4])))

hist(bloom_winter_y_values[,5][which(bloom_winter_y_values[,5] > 0.0)], breaks=15, probability = T)
lines(1:100/10, dgpd(1:100/10, xi = param_res5[3,1], beta = (param_res5[3,2] + param_res5[3,4])))

qqparetoPlot(bloom_winter_y_values[,2][which(bloom_winter_y_values[,2] > 0.0)], xi = param_res5[1,1])
qqparetoPlot(bloom_winter_y_values[,4][which(bloom_winter_y_values[,4] > 0.0)], xi = param_res5[2,1])
qqparetoPlot(bloom_winter_y_values[,5][which(bloom_winter_y_values[,5] > 0.0)], xi = param_res5[3,1])

res <- c(.390081,0.048149,    0.522284,    0.210688 ,  22.903407   , 0.548827   , 9.596196 ,  21.612637 ,   0.436637)
res


run_deoptim_multivar(1, deltas, transfos)
